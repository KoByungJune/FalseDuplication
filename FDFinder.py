######################################
#  False duplication identification  #
######################################


import numpy as np
import time
import pickle
import os
import sys
import parmap
import subprocess as sp
import random
from cigar import Cigar
from multiprocessing import Manager


def Main0_ConfigLoad(ConfigPath):
    PD = {}
    with open(ConfigPath, 'r') as cr: Params = cr.read().split('\n')
    while Params.count('') >0: Params.remove('')
    for i in Params:
        DT = i.split('\t')
        Key = DT[0]
        if len(DT) == 2: Val = DT[1]
        else: Val = DT[1:]
        PD[Key] = Val
    print(PD)
    InputMaf = PD['InputMaf'] 
    f0 = open(InputMaf, 'r')
    Maf = f0.read().split('\na\n')
    Maf = Maf[1:]
    print("#MafBlocks: " + str(len(Maf)))
    f0.close()
    return PD, Maf

### [Simple Tools]
def BinSaving(JobPath, PythonData, FileName):#c     
    with open(JobPath + FileName+'.Pybin', 'wb') as fb: pickle.dump(PythonData, fb)

def BinLoading(JobPath, FileName):#c      
    with open(JobPath + FileName+'.Pybin', 'rb') as rb: return pickle.load(rb)

def Overlapping(QS, QE, TS, TE):#c
    OvLen = 0
    QS, QE, TS, TE = int(QS), int(QE), int(TS), int(TE)
    if QS <= TE and TS <= QE:
        OvList = [QS, QE, TS, TE]
        OvList.remove(max(OvList))
        OvList.remove(min(OvList))
        OvLen = max(OvList) - min(OvList) + 1
    else: pass
    return OvLen

def CEDCalculation(ContigName, FDLoci):#c
    ContigPos, FDPos = Str2Pos(ContigName), Str2Pos(FDLoci)
    ContigChr, FDChr = ContigPos[0], FDPos[0]
    if ContigChr != FDChr: return 'Diff Chr'
    ContigS, ContigE = int(ContigPos[1]), int(ContigPos[2])
    Start, End = int(FDPos[1]), int(FDPos[2])
    CED = min(Start - ContigS, ContigE - End)
    return CED

def CEDCalculation_Sol(ContigName, Pos):#c
    ContigPos = Str2Pos(ContigName)
    ContigS, ContigE = int(ContigPos[1]), int(ContigPos[2])
    CED = min(int(Pos) - ContigS, ContigE - int(Pos))
    return CED

def PosLenSum(DoubleList):#c
    Len = 0
    for i in DoubleList:
        Start, End = int(i[0]), int(i[1])
        Len += End - Start + 1
    return Len


### [LociStr Processing]
def Pos2Str(Chr, Start, End):#c
    Loci = Chr + ':' + str(Start) + '-' + str(End)
    return Loci
def Str2Pos(Loci):#c
    Chr, Start, End = Loci.split(':')[0], Loci.split(':')[1].split('-')[0], Loci.split('-')[1]
    return [Chr, Start, End]

### [Read Inf Parsing]
def DepthCalling(Loci, BamPath):#c
    LL = Str2Pos(Loci)
    Chr, Start, End = LL[0], LL[1], LL[2]
    Len = int(End) - int(Start) + 1
    TotalDepth = sp.getoutput("samtools depth -r {0} {1} | awk 'BEGIN{{Depth = 0}} {{Depth += $3}} END{{print Depth}}'".format(Loci, BamPath))
    MeanDepth = int(TotalDepth) / Len
    return MeanDepth

def ReadParser(Loci, Bam):#c
    Reads = sp.getoutput("samtools view {0} {1} | awk '{{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7\"\t\"$8\"\t\"$9}}'".format(Bam, Loci)).split('\n')
    return Reads
def PMQPFSParser(RN, Bam, Loci):#c
    Reads = sp.getoutput("samtools view {0} {1} | awk '$1 == \"{2}\" {{print $1\"\t\"$2\"\t\"$5\"\t\"$6}}'".format(Bam, Loci, RN)).split('\n')
    for i in Reads:
        DT = i.split('\t')
        PRN, PFS, PMQ, CS = DT[0], DT[1], DT[2], DT[3]
        if PRN == RN: return [PFS, PMQ, CS]

### [Read Inf Processing]
def Cigar4Len(CigarString):#c
    CigarList = list(Cigar(CigarString).items())
    SoftClip = 0
    for i in CigarList:
        if i[1] == 'S': SoftClip = i[0]
    ReadLen = len(Cigar(CigarString)) - SoftClip
    return ReadLen

def FSDecimal(FS):#c
    DecimalDic = {} # {Decinal : binary}
    FS, Decimal = int(FS), 2048
    while Decimal > 0.9:
        if FS >= Decimal:
            DecimalDic[round(Decimal, 0)] = 1
            FS -= Decimal
            Decimal = Decimal / 2
        else:
            DecimalDic[round(Decimal, 0)] = 0
            Decimal = Decimal / 2
    return DecimalDic

def FFRRFinder(DecimalDic1, DecimalDic2):#c # RR = 32 and 16; LL = non of 32, 16
    DD1, DD2 = DecimalDic1, DecimalDic2
    if [DD1[32], DD1[16], DD2[32], DD2[16]] == [1, 1, 1, 1]: Type = "RR"
    elif [DD1[32], DD1[16], DD2[32], DD2[16]] == [0, 0, 0, 0]: Type = "FF"
    else: Type = "Diff"
    return Type


### [AGDG preocessing]
def GCA2AG(NCBIGapPath, ChrConvDic, JobPath, AsmNo):
    os.system("awk -F \"\t\" 'index($1,\"#\") == 0 {split($1,arr,\".\"); print arr[1]\"\t\"$2\"\t\"$3}' " + NCBIGapPath + " > {0}{1}AsmGap.txt".format(JobPath, AsmNo))

def AllDepth4DepthGap(BamPath, JobPath, AsmNo):
    os.system("samtools depth -a {0} > {1}{2}AllDepth.txt".format(BamPath, JobPath, AsmNo))
    os.system("awk -F \"\t\" '$3 == 0 {print $1\"\t\"$2-1\"\t\"$2}' " + "{0}{1}AllDepth.txt > {0}{1}ZeroDepth.bed".format(JobPath, AsmNo))
    os.system("bedtools merge -i {0}{1}ZeroDepth.bed | awk -F \"\t\" '{{print $1\"\t\"$2 + 1\"\t\"$3}}' > {0}{1}DepGap.txt".format(JobPath, AsmNo))

### [Kmer Dic Building]
def Merqury2Loci(AssemblyPath, AssemblyKmerDBPath, ReadKmerDBPath, JobPath, AsmNo):
    os.system("meryl-lookup -dump -sequence {0} -mers {1} > {2}{3}AssemblyKmer.bed".format(AssemblyPath, AssemblyKmerDBPath, JobPath, AsmNo))
    #print("meryl-lookup -dump -sequence {0} -mers {1} > {2}{3}AssemblyKmer.bed".format(AssemblyPath, AssemblyKmerDBPath, JobPath, AsmNo))
    os.system("meryl-lookup -dump -sequence {0} -mers {1} > {2}{3}ReadKmer.bed".format(AssemblyPath, ReadKmerDBPath, JobPath, AsmNo))
    os.system("paste {0}{1}AssemblyKmer.bed {0}{1}ReadKmer.bed > {0}{1}_kmer.txt".format(JobPath, AsmNo))

def KmerDicBulding(JobPath, AsmNo):
    KmerBedPath, KmerDicPath = '{0}{1}_kmer.txt'.format(JobPath, AsmNo), '{0}{1}_kmerDic/'.format(JobPath, AsmNo)
    os.system('mkdir {0}'.format(KmerDicPath))
    f0 = open(KmerBedPath,'r')
    KmerDic, Step = {'null':0}, 0
    while True:
        Step += 1
        #if Step % 1000000 == 0:  print(Step)
        line = f0.readline()
        if not line:
            with open(KmerDicPath + '{0}.dic'.format(PreChr), 'wb') as fw: pickle.dump(KmerDic, fw)
            break
        DT = line.split('\t')
        Chr, Start, FC, RC, FD, RD, Pal  = DT[0], DT[2], int(DT[5]), int(DT[7]), int(DT[13]), int(DT[15]), 'N'
        if FC == RC: # Palindromic
            CN, Depth, Pal = FC, FD, 'Y'
        else: CN, Depth = FC + RC, FD + RD

        if Chr in KmerDic:
            KmerDic[Chr][int(Start)] = [CN, Depth, Pal]
        else:
            if KmerDic == {'null':0}: pass
            else:
                with open(KmerDicPath + '{0}.dic'.format(PreChr), 'wb') as fw:
                    pickle.dump(KmerDic, fw)
            KmerDic = {}
            KmerDic[Chr] = {int(Start):[CN, Depth, Pal]}
        PreChr = Chr
    

### [Pre-Processing]
def Fai2Bed(PreFaiPath, VGPFaiPath, JobPath):#c # Preprocessing; contig formation; both assemblies
    os.system("awk '{{print $1\"\t\"0\"\t\"$2}}' {0} > {1}1Faibed.bed".format(PreFaiPath, JobPath)) # Fai2bed
    os.system("awk '{{print $1\"\t\"0\"\t\"$2}}' {0} > {1}2Faibed.bed".format(VGPFaiPath, JobPath))


def Fai2Dic(PreFaiPath, VGPFaiPath):
    FaiDic = {} # {Chr :[Start-End]..}
    PreFai = sp.getoutput("awk '{{print $1\"\t\"0\"\t\"$2}}' {0}".format(PreFaiPath)).split('\n')
    VGPFai = sp.getoutput("awk '{{print $1\"\t\"0\"\t\"$2}}' {0}".format(VGPFaiPath)).split('\n')
    if PreFai.count('') > 0: PreFai.remove('')
    if VGPFai.count('') > 0: VGPFai.remove('')
    for i in PreFai + VGPFai:
        DT = i.split('\t')
        Chr, Start, End = DT[0], DT[1], DT[2]
        FaiDic[Chr] = [Start, End]
    return FaiDic


def ChrConversion(ChrPairPath):#c # Preprocessing; sca. name conversion; both assemblies
    f3 = open(ChrPairPath, 'r')
    ChrPair = f3.read().split('\n')
    f3.close()
    if ChrPair.count('') >0: ChrPair.remove('')
    ChrPairDic = {}
    for i in ChrPair:
        DT = i.split('\t')
        CM, NC = DT[0], DT[1]
        ChrPairDic[CM] = NC
    return ChrPairDic


def AGDGmerging(JobPath, ChrPairDic):#c # Preprocessing; Gap load; Both assemblies 
    FinalGap = {} # {'Chr:Start-End':Type}
    f_PreAG, f_PreDG, = open(JobPath + '1AsmGap.txt' ), open(JobPath + '1DepGap.txt' ) # KDic = {'Chr:Start' : Depth}
    f_VGPAG, f_VGPDG, = open(JobPath + '2AsmGap.txt' ), open(JobPath + '2DepGap.txt' )
    PreAG, PreDG, VGPAG, VGPDG = f_PreAG.read().split('\n'), f_PreDG.read().split('\n'), f_VGPAG.read().split('\n'), f_VGPDG.read().split('\n')
    if PreAG.count('') > 0: PreAG.remove('')
    if PreDG.count('') > 0: PreDG.remove('')
    if VGPAG.count('') > 0: VGPAG.remove('')
    if VGPDG.count('') > 0: VGPDG.remove('')
    file_list_Pre, file_list_VGP = os.listdir(JobPath + '1_kmerDic/'), os.listdir(JobPath + '2_kmerDic/')
    KGapDic = {}
    def ZeroKmerParsing(file_list, DicPath, JobPath):
        Step = 0
        for i in file_list:
            Step += 1
            #print('zeroK parsing: ' + str(len(file_list) - Step))
            fk = open(JobPath + DicPath + i,'rb')
            Dic = pickle.load(fk)
            for j in Dic:
                Chr, DD = j, Dic[j]
                for k in DD:
                    Pos, Copy, Multi = k, DD[k][0], DD[k][1]
                    if Multi == 0:
                        if Chr in KGapDic: KGapDic[Chr].append([str(int(Pos) + 1), str(int(Pos)+20)])
                        else: KGapDic[Chr] = [[str(int(Pos) + 1), str(int(Pos)+20)]] # K-mer = bed format (Pos+1 - Pos+k)
    
    ZeroKmerParsing(file_list_Pre, '1_kmerDic/', JobPath) 
    ZeroKmerParsing(file_list_VGP, '2_kmerDic/', JobPath) #Zero K-mer save
    ##BinSaving(KGapDic, '04KGapDic') 
    ##KGapDic = BinLoading('04KGapDic')
    AGDic = {}
    def AGSaving(GapList):
        for i in GapList:
            DT = i.split('\t')
            Chr, Start, End = ChrPairDic[DT[0]], DT[1], DT[2]
            Loci = Pos2Str(Chr, Start, End) #none-bed start
            FinalGap[Loci] = 'AG'
            if Chr in AGDic: AGDic[Chr].append([Start, End])
            else: AGDic[Chr] = [[Start, End]]
    
    AGSaving(PreAG)
    AGSaving(VGPAG)
    #print("AGSaving")
    def DGSaving(GapList):
        Step = 0
        for i in GapList:
            #print(len(GapList) - Step)
            Step += 1
            DT = i.split('\t')
            Chr, Start, End = DT[0], DT[1], DT[2]
            Loci = Pos2Str(Chr, Start, End) #none-bed start
            if Chr in AGDic: AGList = AGDic[Chr]
            else: AGList = [[-1, -1]]
            if Chr in KGapDic: KGapList = KGapDic[Chr]
            else: KGapList = [[-1, -1]]
            KGap = 0
            for j in KGapList:
                KGapStart, KGapEnd = j[0], j[1]
                if int(Start) <= int(KGapEnd) and int(KGapStart) <= int(End): 
                    KGap = 1 # KGap Support DepthGap
                    break
            if KGap == 0: continue # Repeat avoid (DG by repeat) 
            AGOVLP = 0
            for j in AGList:
                AGStart, AGEnd = j[0], j[1]
                if int(Start) <= int(AGEnd) and int(AGStart) <= int(End): 
                    AGOVLP = 1
                    break
            if AGOVLP == 1: continue # AG in DG
            else: FinalGap[Loci] = 'DG'
    
    DGSaving(PreDG)
    #print("PreDGSaving")
    DGSaving(VGPDG)
    #print("VGPDGSaving")
    with open(JobPath + 'BothGapbed.bed', 'w') as fG:
        for i in FinalGap: 
            Chr, Start, End, Type = i.split(':')[0], i.split(':')[1].split('-')[0], i.split('-')[1], FinalGap[i]
            fG.write('\t'.join([Chr, str(int(Start)-1), End, Type]) + '\n') # Gap2bed


def ContigbyAGDG(AssemblyType, JobPath):#c # Preprocessing; contig segmenting; factors for assembly.
    ContigBed = sp.getoutput("bedtools subtract -a {0} -b {1}".format(JobPath + AssemblyType + 'Faibed.bed', JobPath + 'BothGapbed.bed')).split('\n')
    if ContigBed.count('') >0: ContigBed.remove('') 
    ContigDic = {}
    for i in ContigBed: # bed2back
        DT = i.split('\t')
        Chr, Start, End = DT[0], str(int(DT[1]) + 1), DT[2]
        Loci = Chr + ':' + Start + '-' + End
        if Chr in ContigDic: ContigDic[Chr][Loci] = {'FDPos':[],'FDBlock':[]}
        else: ContigDic[Chr] = {Loci:{'FDPos':[],'FDBlock':[]}} # {Chr:{Loci1:{'FDPos':[Start, End], 'FDBlock':[BLockNo1, BLockNo2]},Loci2:{}}
    
    return ContigDic


def PurgeDupPafLoad(PafPath):#c
    with open(PafPath, 'r') as fr:
        Paf = fr.read().split('\n')
        if Paf.count('') > 0: Paf.remove('')
        PafDic = {} # Region: other inf
        for i in Paf:
            DT = i.split('\t')
            Query, QueryDL, Target = DT[0], DT, DT[5]
            if Query in PafDic: 
                if Target in PafDic[Query]: PafDic[Query][Target].append(QueryDL)
                else: PafDic[Query][Target] = [QueryDL]
            else: 
                PafDic[Query] = {Target:[QueryDL]}
    return PafDic


def PurgeDupLociSpecify(PurgeList, PurgePafDic):
    PADic = {}
    for i in PurgeList:
        if len(i) < 5: continue
        PurgeNo, FChr, FStart, FEnd, TChr, TStart, TEnd, Type = i[0], i[1], int(i[2]), int(i[3]), i[6], int(i[4]), int(i[5]), i[7] # non-bed
        FLoci, TLoci = Pos2Str(FChr, FStart, FEnd), Pos2Str(TChr, TStart, TEnd)
        #print(PurgeNo)
        if Type == "HAPLOTIG":
            ADL = PurgePafDic[FLoci][TLoci]
            for j in ADL: # non-bed loci info
                #Alignment = sp.getoutput("awk '$1 == {0} && $6 == {1} {{print $4-$3\"\t\"$9-$8\"\t\"$0}}' {2}".format(FLoci, TLoci, PafPath)).split('\n')
                FAStart, FAEnd, TAStart, TAEnd, Direction = int(j[2]), int(j[3]), int(j[7]), int(j[8]), j[4] # bed
                FALen, TALen = FAEnd - FAStart, TAEnd - TAStart
                ReFStart, ReFEnd, ReTStart, ReTEnd = FStart + FAStart, FStart + FAEnd - 1, TStart + TAStart, TStart + TAEnd - 1
                if PurgeNo in PADic:
                    if (FALen + TALen) / 2 > PADic[PurgeNo][3]: PADic[PurgeNo] = [Pos2Str(FChr, ReFStart, ReFEnd), Pos2Str(TChr, ReTStart, ReTEnd), Direction, (FALen + TALen) / 2] # bed conversion considered
                    else: pass
                else: PADic[PurgeNo] = [Pos2Str(FChr, ReFStart, ReFEnd), Pos2Str(TChr, ReTStart, ReTEnd), Direction, (FALen + TALen) / 2]
        if Type == "OVLP": PADic[PurgeNo] = [Pos2Str(FChr, FStart, FEnd), Pos2Str(TChr, TStart, TEnd), 'n', ((FEnd - FStart + 1) + (TEnd - TStart + 1)) / 2] # bed conversion considered

    return PADic


#[Main Start]
def Main1_CactusCalling(ProcessedMaf): # Main Script; Cactus Loci Parsing; Both assemblies
    MafDic, Step, DirectionDic = {}, 0, {} #{Loci:relative direction from reference}
    for i in ProcessedMaf:
        Step += 1
        #if (len(ProcessedMaf) - Step) % 100 == 0: print('Main1\t' + str(len(ProcessedMaf) - Step))
        Breaker = 0
        DL = i.split('\n')
        while DL.count('') > 0: DL.remove('')

        LenList, RefList, TarList, AltList = [], [], [], []
        RefPos, TarPos, AltPos = [], [], []
        for j in DL:
            DT = j.split('\t')
            Node, Chr, BedStart, Len, Direction, ChrTotal, Seq = DT[1].split(".")[0], DT[1].split(".")[1] , DT[2], DT[3], DT[4], DT[5], DT[6].upper()
            if Direction == '+':
                Start = int(BedStart) + 1
                End = int(BedStart) + int(Len) #Start + int(Len) : This is cause of +1 End
                DirectionDic[Pos2Str(Chr, Start, End)] = '+'
            elif Direction == '-':
                Start = int(ChrTotal) - int(BedStart) - int(Len) + 1
                End = int(ChrTotal) - int(BedStart)
                DirectionDic[Pos2Str(Chr, Start, End)] = '-'
            NCal = Seq.count('N') # N = Total - 
            NLen = int(Len) - NCal
            if Node == Ref: 
                RefList.append(Chr)  
                LenList.append(int(NLen)) 
                RefPos.append([Start, End])

            elif Node == Tar:
                TarList.append(Chr)
                LenList.append(int(NLen))
                TarPos.append([Start, End])

            elif Node == Alt:
                AltList.append(Chr)
                LenList.append(int(NLen))
                AltPos.append([Start, End])
        
        if max(LenList) == 0: continue
        if min(LenList)/max(LenList) >= LenCutoff: pass
        else: continue
        if min(LenList) >= Lenlimit: pass
        else: continue
        if len(RefList) == 0 or len(TarList) == 0: continue 

        MafType = str(len(RefList)) + '_' + str(len(AltList)) + '_' + str(len(TarList))        
        MafNo = str(Step)
        MafDic[MafNo] = {} #{MafNo:{'NoSeq':[],'Dat':{'R1':{'Chr':,'Start':,'End':,'CS':},..}}
        MafDic[MafNo]['NoSeq'] = [len(RefList), len(AltList), len(TarList)]
        MafDic[MafNo]['Dat'] = {}
                
        Idx = 0
        for k in RefList:
            Chr, Start, End, Name = k, RefPos[Idx][0], RefPos[Idx][1], 'R' + str(Idx + 1) #ChrPairDic[k] collapsed by conversed inputMaf
            MafDic[MafNo]['Dat'][Name] = {'Chr':Chr,'Start':str(Start),'End':str(End)}
            Key = '_'.join([Chr, str(Start), str(End)])
            Idx += 1
          
        Idx = 0
        for k in TarList:
            Chr, Start, End, Name = k, TarPos[Idx][0], TarPos[Idx][1], 'T' + str(Idx + 1) #ChrPairDic[k] collapsed by conversed inputMaf
            MafDic[MafNo]['Dat'][Name] = {'Chr':Chr,'Start':str(Start),'End':str(End)}
            Key = '_'.join([Chr, str(Start), str(End)])
            Idx += 1

    return MafDic, DirectionDic


def Main2_ContigCalling(InputMafDic): # Main Script; Contig analyzing; Both assemblies
    FDDic = {} # {MafNo:{'NoSEq':[0,0,0], 'Dat':{'R1':{'Chr':,'Start':,'End':,'CS':,'Depth':,},'T1'..}}}
    DepthInput = []
    Step, LenLen = 0, len(InputMafDic)
    ContigPreListDic, ContigVGPListDic = {}, {}
    for j in ContigPreDic:
        for k in ContigPreDic[j].keys():
            Chr, Start, End = k.split(':')[0], int(k.split(':')[1].split('-')[0]), int(k.split('-')[1])
            if Chr in ContigPreListDic: ContigPreListDic[Chr].append([Start, End])
            else: ContigPreListDic[Chr] = [[Start, End]]
    for j in ContigVGPDic:
        for k in ContigVGPDic[j].keys():
            Chr, Start, End = k.split(':')[0], int(k.split(':')[1].split('-')[0]), int(k.split('-')[1])
            if Chr in ContigVGPListDic: ContigVGPListDic[Chr].append([Start, End])
            else: ContigVGPListDic[Chr] = [[Start, End]]
    for j in ContigPreListDic: ContigPreListDic[j].sort()
    for j in ContigVGPListDic: ContigVGPListDic[j].sort()
    for i in InputMafDic:
        Step += 1
        #if (LenLen - Step) % 100 == 0: print('Main2\t' + str(LenLen - Step))
        #if (LenLen - Step) % 1000 == 0: print(str(LenLen - Step) + '\t' + "AddInf")
        MafNo, DataDic, PosDic = i, InputMafDic[i], InputMafDic[i]['Dat']
        Reflen, Tarlen = DataDic['NoSeq'][0], DataDic['NoSeq'][2]
        if Reflen == 1 and Tarlen > 1:
            FDDic[MafNo] = {'NoSeq':DataDic['NoSeq'],'Dat':{}}
            for j in PosDic:
                if j[0] == 'T':
                    SubDic = PosDic[j]
                    Chr, Start, End, Name = SubDic['Chr'], SubDic['Start'], SubDic['End'], j
                    FDDic[MafNo]['Dat'][j] = PosDic[j]
                    
                    CC = {} # Candidate Contigs : {Overlap:[CS,CCD, CED], Overlap2:CS,..}
                    if Chr in ContigPreListDic: # 20210201 (Full DG Scaffold in CalannPre)
                        for k in ContigPreListDic[Chr]: # Contig Searching; Must be sorted
                            ContigS, ContigE = int(k[0]), int(k[1])
                            CN, LociName = Pos2Str(Chr, ContigS, ContigE), Pos2Str(Chr, Start, End)
                            if Overlapping(int(Start), int(End), ContigS, ContigE) > 0:
                                CS = ContigE - ContigS + 1
                                CCD = abs(((ContigE + ContigS)/ 2) - ((int(Start) + int(End))/2)) # Distance of Block Center from Contig Center
                                CED = CEDCalculation(CN, LociName) # FD Distance from gap (Contig tip) # i think min(int(Start) - contigS, ContigE - int(End)) will be better. (over area => minus)
                                if ContigE < int(End):  
                                    CC[Overlapping(int(Start), int(End), ContigS, ContigE)] = [CS, CCD, CED, CN]
                                else:  
                                    CC[Overlapping(int(Start), int(End), ContigS, ContigE)] = [CS, CCD, CED, CN]
                                    break
                            else:  pass
                        
                        if len(CC) == 0: # FD on DG operation 
                            LociName == Pos2Str(Chr, Start, End)
                            CN = ('DG' + LociName)
                            CS = int(End) - int(Start) + 1
                            CCD = 0
                            CED = 0
                            CC[0] = [CS, CCD, CED, CN]
                            '''
                            Dist = 1000000000
                            for k in ContigPreListDic[Chr]:
                                ContigS, ContigE = int(k[0]), int(k[1])
                                CN, LociName = Pos2Str(Chr, ContigS, ContigE), Pos2Str(Chr, Start, End)
                                CS = ContigE - ContigS + 1
                                CCD = abs(((ContigE + ContigS)/ 2) - ((int(Start) + int(End))/2)) # Distance of Block Center from Contig Center
                                CED = CEDCalculation(CN, LociName)
                                
                                ##### Don't know what this mean
                                if CED > Dist: break
                                CC[int(End) - int(Start) + 1] = [CS, CCD, CED, CN]
                                Dist = CED
                            '''
                    else: # Full DG Scaffold 20210201
                        LociName = Pos2Str(Chr, Start, End)
                        CN = ('DG' + LociName)
                        CS = int(End) - int(Start) + 1
                        CCD = 0
                        CED = 0
                        CC[0] = [CS, CCD, CED, CN] # actually no overlap with real 'contig' (because whole DG area on contig)

                    CS, CCD, CED, CN = CC[max(CC)][0], CC[max(CC)][1], CC[max(CC)][2], CC[max(CC)][3] # Largest Overllaped Contig
                    FDDic[MafNo]['Dat'][j]['CS'], FDDic[MafNo]['Dat'][j]['CCD'], FDDic[MafNo]['Dat'][j]['CED'], FDDic[MafNo]['Dat'][j]['CN'] = CS, CCD, CED, CN #NewWrite
                    DepthInput.append((Chr, Start, End, PreBamPath, '1_xTar'))

                elif j == 'R1':
                    SubDic = PosDic[j]
                    Chr, Start, End, Name = SubDic['Chr'], SubDic['Start'], SubDic['End'], j
                    FDDic[MafNo]['Dat'][j] = PosDic[j]
                    for k in ContigVGPListDic[Chr]:
                        ContigS, ContigE = int(k[0]), int(k[1])
                        CN = Pos2Str(Chr, ContigS, ContigE)
                        if Overlapping(int(Start), int(End), ContigS, ContigE) > 0:
                            CS = ContigE - ContigS + 1
                            break
                        else: pass

                    FDDic[MafNo]['Dat'][j]['CS'] = CS
                    DepthInput.append((Chr, Start, End, VGPBamPath, '1_xRef'))

        elif Reflen > 1 and Tarlen == 1: # VGP duplication
            FDDic[MafNo] = {'NoSeq':DataDic['NoSeq'],'Dat':{}} # Load and Write
            for j in PosDic:
                if j[0] == 'R':
                    SubDic = PosDic[j]
                    Chr, Start, End, Name = SubDic['Chr'], SubDic['Start'], SubDic['End'], j
                    FDDic[MafNo]['Dat'][j] = PosDic[j]

                    CC = {}
                    if Chr in ContigVGPListDic:
                        for k in ContigVGPListDic[Chr]:
                            ContigS, ContigE = int(k[0]), int(k[1])
                            CN, LociName = Pos2Str(Chr, ContigS, ContigE), Pos2Str(Chr, Start, End)
                            if Overlapping(int(Start), int(End), ContigS, ContigE) > 0:
                                CS = ContigE - ContigS + 1
                                CCD = abs(((ContigE + ContigS)/ 2) - ((int(Start) + int(End))/2))
                                CED = CEDCalculation(CN, LociName)
                                if ContigE < int(End):  
                                    CC[Overlapping(int(Start), int(End), ContigS, ContigE)] = [CS, CCD, CED, CN]
                                else:  
                                    CC[Overlapping(int(Start), int(End), ContigS, ContigE)] = [CS, CCD, CED, CN]
                                    break
                            else: pass
                        
                        if len(CC) == 0: # FD on DG operation 
                            LociName == Pos2Str(Chr, Start, End)
                            CN = ('DG' + LociName)
                            CS = int(End) - int(Start) + 1
                            CCD = 0
                            CED = 0
                            CC[0] = [CS, CCD, CED, CN]
                            '''
                            Dist = 1000000000
                            for k in ContigVGPListDic[Chr]:
                                ContigS, ContigE = int(k[0]), int(k[1])
                                CN, LociName = Pos2Str(Chr, ContigS, ContigE), Pos2Str(Chr, Start, End)
                                CS = ContigE - ContigS + 1
                                CCD = abs(((ContigE + ContigS)/ 2) - ((int(Start) + int(End))/2)) # Distance of Block Center from Contig Center
                                CED = CEDCalculation(CN, LociName)
                                #print([CN, CS, CED])
                                if CED > Dist: break
                                CC[int(End) - int(Start) + 1] = [CS, CCD, CED, CN]
                                Dist = CED
                            '''
                    else:
                        LociName = Pos2Str(Chr, Start, End)
                        CN = ('DG' + LociName)
                        CS = int(End) - int(Start) + 1
                        CCD = 0
                        CED = 0
                        CC[0] = [CS, CCD, CED, CN] # actually no overlap with real 'contig' (because whole DG area on contig)

                    CS, CCD, CED, CN = CC[max(CC)][0], CC[max(CC)][1], CC[max(CC)][2], CC[max(CC)][3]
                    FDDic[MafNo]['Dat'][j]['CS'], FDDic[MafNo]['Dat'][j]['CCD'], FDDic[MafNo]['Dat'][j]['CED'], FDDic[MafNo]['Dat'][j]['CN'] = CS, CCD, CED, CN
                    DepthInput.append((Chr, Start, End, VGPBamPath, 'x_1Ref'))

                elif j == 'T1':
                    SubDic = PosDic[j]
                    Chr, Start, End, Name = SubDic['Chr'], SubDic['Start'], SubDic['End'], j
                    FDDic[MafNo]['Dat'][j] = PosDic[j]
                    for k in ContigPreListDic[Chr]:
                        ContigS, ContigE = int(k[0]), int(k[1])
                        CN = Pos2Str(Chr, ContigS, ContigE)
                        if Overlapping(int(Start), int(End), ContigS, ContigE) > 0:
                            CS = ContigE - ContigS + 1
                            break
                        else: pass
                    
                    FDDic[MafNo]['Dat'][j]['CS'] = CS
                    DepthInput.append((Chr, Start, End, PreBamPath, 'x_1Tar'))

        elif Reflen == 1 and Tarlen == 1: # singleton
            FDDic[MafNo] = {'NoSeq':DataDic['NoSeq'],'Dat':DataDic['Dat']}
            for j in PosDic:
                if j == 'R1':
                    SubDic = PosDic[j]
                    Chr, Start, End, Name = SubDic['Chr'], SubDic['Start'], SubDic['End'], j
                    FDDic[MafNo]['Dat'][j]['CS'] = 'ND'
                    DepthInput.append((Chr, Start, End, VGPBamPath, '1_1Ref'))
                elif j == 'T1':
                    SubDic = PosDic[j]
                    Chr, Start, End, Name = SubDic['Chr'], SubDic['Start'], SubDic['End'], j
                    FDDic[MafNo]['Dat'][j]['CS'] = 'ND'
                    DepthInput.append((Chr, Start, End, PreBamPath, '1_1Tar'))
    
    return FDDic, DepthInput, ContigPreListDic, ContigVGPListDic


### [Main3 SubTool: Depth Data load]###
'''
def M3_Saved_DepthLoad(): # Preprocessing; Depth Load; Both assemblies; from previous output
    f_PreD, f_VGPD = open(DepthPrePath), open(DepthVGPPath)
    PD, VD = f_PreD.read().split('\n'), f_VGPD.read().split('\n')
    if PD.count('') > 0: PD.remove('')
    if VD.count('') > 0: VD.remove('')
    Load_DepthDic = {}
    for i in PD + VD:
        DT = i.split('\t')
        Chr1, Start1, End1, Chr2, Start2, End2, Depth1, Depth2, RChr, RStart, REnd, RDepth = DT[1], DT[2], DT[3], DT[6], DT[4], DT[5], DT[11], DT[12], DT[7], DT[8], DT[9], DT[13]
        Loci1, Loci2, RLoci = Pos2Str(Chr1, Start1, End1), Pos2Str(Chr2, Start2, End2), Pos2Str(RChr, RStart, REnd)
        Load_DepthDic[Loci1] = Depth1
        Load_DepthDic[Loci2] = Depth2
        Load_DepthDic[RLoci] = RDepth

    if os.listdir(LogPath).count('New_DepthDic.Pybin') == 1:
        with open(LogPath + 'New_DepthDic.Pybin', 'rb') as fr:
            Saved_NewDepthDic = pickle.load(fr)
            return {**Load_DepthDic, **Saved_NewDepthDic}
    return Load_DepthDic # 'Chr:Start-End':Depth
'''

def M3_M_DepthMeanOperator(DataList, Cores):
    input_data = np.array_split(DataList, Cores)
    CalledDepthDic = {}
    Result = parmap.map(M3_M_DepthMean, input_data, pm_pbar = True)
    for i in Result:
        CalledDepthDic = {**CalledDepthDic, **i}
    return CalledDepthDic


def M3_M_DepthMean(DataList): # Methods; Depth calculation; factors for assembly
    Step, ResultDic = 0, {}
    for i in DataList:
        Step += 1
        Loci, BamPath = i[0], i[1]
        #if (len(DataList) - Step) % 100 == 0: print(len(DataList) - Step)
        LL = Str2Pos(Loci)
        Chr, Start, End = LL[0], LL[1], LL[2]
        Len = int(End) - int(Start) + 1
        TotalDepth = sp.getoutput("samtools depth -r {0} {1} | awk 'BEGIN{{Depth = 0}} {{Depth += $3}} END{{print Depth}}'".format(Loci, BamPath))
        MeanDepth = int(TotalDepth) / Len
        ResultDic[Loci] = str(round(float(MeanDepth),1))
    return ResultDic

#### END Main3-Subtool ####

def Main3_DetphAdd2FDDic(FDDic, Loaded_DepthDic, Cores, PreBamPath, VGPBamPath):
    InFDDic = FDDic
    Step = 0
    NewCall_DepthDic = {}
    Need2Call = []
    for i in InFDDic:
        NoSeq, PosData = InFDDic[i]['NoSeq'], InFDDic[i]['Dat']
        for j in PosData:
            SubDic = PosData[j]
            Chr, Start, End = SubDic['Chr'], SubDic['Start'], SubDic['End']
            Region = Pos2Str(Chr, Start, End)
            if NoSeq[0] > 1 and NoSeq[2] > 1: continue
            if NoSeq[0] == 1 and NoSeq[2] == 1: pass
            else:
                if Region in Loaded_DepthDic: pass
                else:
                    if j[0] == 'R': Need2Call.append([Region, VGPBamPath])
                    elif j[0] == 'T': Need2Call.append([Region, PreBamPath])
    
    CalledDepthDic = M3_M_DepthMeanOperator(Need2Call, Cores)
    New_DepthDic = {**CalledDepthDic, **Loaded_DepthDic}
    for i in InFDDic:
        Step += 1
        #if (len(InFDDic) - Step) % 100 == 0: print('Main3\t' + str(len(InFDDic) - Step))
        MafNo, NoSeq = i, InFDDic[i]['NoSeq']
        PosData = InFDDic[MafNo]['Dat']
        for j in PosData:
            SubDic = PosData[j]
            Chr, Start, End = SubDic['Chr'], SubDic['Start'], SubDic['End']
            Region = Pos2Str(Chr, Start, End)
            if NoSeq[0] > 1 and NoSeq[2] > 1: continue
            if NoSeq[0] == 1 and NoSeq[2] == 1: pass
            else: 
                InFDDic[MafNo]['Dat'][j]['Depth'] = New_DepthDic[Region]
    return InFDDic, New_DepthDic


def Main4_CactusFalseTrueAllocating(FDDic): # Analyzing; False True Decision; Both assembly
    InFDDic = FDDic
    Step = 0
    for i in InFDDic:
        Step += 1
        #if (len(InFDDic) - Step) % 100 == 0: print('Main4\t' + str(len(InFDDic) - Step))
        MafNo, NoSeqList = i, InFDDic[i]['NoSeq']
        Reflen, Altlen, Tarlen = NoSeqList[0], NoSeqList[1],NoSeqList[2]
        PosData = InFDDic[i]['Dat']
        if Reflen == 1 and Tarlen > 1:
            TarCSDic, TarSEDic = {}, {}
            for j in PosData: 
                if j[0] == 'T':
                    SubDic = PosData[j]
                    Chr, Start, End, CS, Depth, CCD = SubDic['Chr'], SubDic['Start'], SubDic['End'], SubDic['CS'], SubDic['Depth'], SubDic['CCD']  
                    if float(Depth) > PreSECutoff: 
                        TarCSDic[(int(CS), - int(CCD))] = [j, int(CCD), float(Depth)] # contig 중간 artificial duplication 피하기 
                        InFDDic[MafNo]['Dat'][j]['Bool'] = 'F' # evryblocks F allocating
                    else:
                        TarSEDic[(int(CS), - int(CCD))] = [j, int(CCD), float(Depth)]
                        InFDDic[MafNo]['Dat'][j]['Bool'] = 'F'
            if len(TarCSDic) > 0: 
                InFDDic[MafNo]['Dat'][TarCSDic[max(TarCSDic)][0]]['Bool'] = 'T' # non-SE maxed CS is True in one to many
            else: 
                InFDDic[MafNo]['Dat'][TarSEDic[max(TarSEDic)][0]]['Bool'] = 'T' # SE maxed CS is True in one to many of all SE


        elif Tarlen == 1 and Reflen > 1:
            RefCSDic, RefSEDic = {}, {}
            for j in PosData:
                if j[0] == 'R':
                    SubDic = PosData[j]
                    Chr, Start, End, CS, Depth, CCD = SubDic['Chr'], SubDic['Start'], SubDic['End'], SubDic['CS'], SubDic['Depth'], SubDic['CCD']
                    if float(Depth) > VGPSECutoff:
                        RefCSDic[(int(CS), - int(CCD))] = [j, int(CCD), float(Depth)]
                        InFDDic[MafNo]['Dat'][j]['Bool'] = 'F'
                    else:
                        RefSEDic[(int(CS), - int(CCD))] = [j, int(CCD), float(Depth)]
                        InFDDic[MafNo]['Dat'][j]['Bool'] = 'F'
            if len(RefCSDic) > 0:
                InFDDic[MafNo]['Dat'][RefCSDic[max(RefCSDic)][0]]['Bool'] = 'T'
            else:
                InFDDic[MafNo]['Dat'][RefSEDic[max(RefSEDic)][0]]['Bool'] = 'T'
    return InFDDic

    
def Main56_DepthAndGapFiltering(FDDic):
    InFDDic = FDDic
    PreFDList, VGPFDList = [], []
    TFDNo, RFDNo = 0, 0
    Step = 0
    for i in InFDDic:
        Step += 1
        #if (len(InFDDic) - Step) % 100 == 0: print('Main56\t' + str(len(InFDDic) - Step))
        #print(str(i) + '\t' + str(InFDDic[i]))
        MafNo, NoSeqList = i, FDDic[i]['NoSeq']
        Reflen, Altlen, Tarlen = NoSeqList[0], NoSeqList[1],NoSeqList[2]
        Type = str(Reflen) + '_' + str(Altlen) + '_' + str(Tarlen)
        PosData = FDDic[MafNo]['Dat']
        if Reflen == 1 and Tarlen > 1:
            TTarDat, FTarDat = [], []
            for j in PosData:
                Node = j
                if Node[0] == 'R':
                    DD = PosData[Node]
                    Chr, Start, End, CS, Depth = DD['Chr'], DD['Start'], DD['End'], DD['CS'], DD['Depth']
                    RefDat = [Chr, Start, End, CS, Depth]

                if Node[0] == 'T':
                    DD = PosData[Node]
                    Chr, Start, End, CS, Depth, FDBool, CED, CN = DD['Chr'], DD['Start'], DD['End'], DD['CS'], DD['Depth'], DD['Bool'], DD['CED'], DD['CN']
                    if FDBool == 'T':  TTarDat.append([Chr, Start, End, CS, Depth, FDBool, CED, CN])
                    else:              FTarDat.append([Chr, Start, End, CS, Depth, FDBool, CED, CN])
            
            for j in FTarDat: 
                DL = j
                TChr, TStart, TEnd, TCS, TDepth, TCED, TCN = TTarDat[0][0], TTarDat[0][1], TTarDat[0][2], TTarDat[0][3], TTarDat[0][4], TTarDat[0][6], TTarDat[0][7]
                FChr, FStart, FEnd, FCS, FDepth, FCED, FCN = DL[0], DL[1], DL[2], DL[3], DL[4], DL[6], DL[7]
                if HapChr.count(FChr) == 0 and float(FDepth) < PreCutoff and TCN != FCN: pass
                elif HapChr.count(FChr) == 1 and float(FDepth) < PreHapCutoff and TCN != FCN: pass
                else: continue
                RChr, RStart, REnd, RCS, RDepth = RefDat[0], RefDat[1], RefDat[2], RefDat[3], RefDat[4]
                TFDNo += 1
                PreFilterPass = ['Cactus '+ str(TFDNo), FChr, str(FStart), str(FEnd), str(TStart), str(TEnd), TChr, RChr, str(RStart), str(REnd), str(int(REnd) - int(RStart) + 1), str(FDepth), str(TDepth), str(RDepth), MafNo, str(FCED), str(FCN), str(TCED), str(TCN)]
                PreFDList.append(PreFilterPass)

        elif Tarlen == 1 and Reflen > 1:
            TRefDat, FRefDat = [], []
            for j in PosData:
                Node = j
                if Node[0] == 'T':
                    DD = PosData[Node]
                    Chr, Start, End, CS, Depth = DD['Chr'], DD['Start'], DD['End'], DD['CS'], DD['Depth']
                    TarDat = [Chr, Start, End, CS, Depth]

                if Node[0] == 'R':
                    DD = PosData[Node]
                    Chr, Start, End, CS, Depth, FDBool, CED, CN = DD['Chr'], DD['Start'], DD['End'], DD['CS'], DD['Depth'], DD['Bool'], DD['CED'], DD['CN']
                    if FDBool == 'T':  TRefDat.append([Chr, Start, End, CS, Depth, FDBool, CED, CN])
                    else:              FRefDat.append([Chr, Start, End, CS, Depth, FDBool, CED, CN])

            for j in FRefDat:
                DL = j
                TChr, TStart, TEnd, TCS, TDepth, TCED, TCN = TRefDat[0][0], TRefDat[0][1], TRefDat[0][2], TRefDat[0][3], TRefDat[0][4], TRefDat[0][6], TRefDat[0][7]
                FChr, FStart, FEnd, FCS, FDepth, FCED, FCN = DL[0], DL[1], DL[2], DL[3], DL[4], DL[6], DL[7]
                if HapChr.count(FChr) == 0 and float(FDepth) < VGPCutoff and TCN != FCN: pass
                elif HapChr.count(FChr) == 1 and float(FDepth) < VGPHapCutoff and TCN != FCN: pass
                else: continue
                RChr, RStart, REnd, RCS, RDepth = TarDat[0], TarDat[1], TarDat[2], TarDat[3], TarDat[4]
                RFDNo += 1
                VGPFilterPass = ['Cactus ' + str(RFDNo), FChr, str(FStart), str(FEnd), str(TStart), str(TEnd), TChr, RChr, str(RStart), str(REnd), str(int(REnd) - int(RStart) + 1), str(FDepth), str(TDepth), str(RDepth), MafNo, str(FCED), str(FCN), str(TCED), str(TCN)]
                VGPFDList.append(VGPFilterPass)
    return PreFDList, VGPFDList


def M7_OverlapMerge(List):
    SortList = []
    for i in List:
        Start, End = int(i[0]), int(i[1])
        SortList.append([Start, End])
    SortList.sort()
    MergeList = []
    for i in SortList:
        if len(MergeList) == 0: MergeList.append(i)
        else:
            Start, End, PrevStart, PrevEnd = i[0], i[1], MergeList[-1][0], MergeList[-1][1]
            if Start <= PrevEnd and PrevStart <= End:
                if Start <= PrevStart and PrevEnd <= End: MergeList[-1] = [Start, End]
                elif PrevStart <= Start and PrevEnd <= End: MergeList[-1] = [PrevStart, End]
                elif Start <= PrevStart and End <= PrevEnd: MergeList[-1] = [Start, PrevEnd]
                elif PrevStart <= Start and End <= PrevEnd: pass
            else:
                MergeList.append(i)
    return MergeList


def Extra_NonBedNonOverlap_ListLen(List):
    Len = 0
    for i in List: Len += int(i[1]) - int(i[0]) + 1
    return Len


def Main7_CactusXPurgeDup(CactusBlockList, PurgeDupPath, ContigDic, ContigListDic):
    f_Purge = open(PurgeDupPath, 'r')
    PD = f_Purge.read().split('\n')
    InContigDic = ContigDic
    PurgeBlockList = []
    if PD.count('') > 0: PD.remove('')
    
    PurgeNo = 1
    for i in PD:
        DT = i.split('\t')
        if len(DT) == 4: continue
        PurgeFDNo, FChr, FStart, FEnd, Type, TChr ,TStart, TEnd = 'PurgeDup ' + str(PurgeNo), DT[0], str(int(DT[1]) + 1), DT[2], DT[3], DT[4], str(int(DT[5]) + 1), DT[6]
        if Type == "HAPLOTIG" or Type == "OVLP": pass
        else: continue
        if Type == "OVLP":
            Saver = ''
            FCN, TCN = '', '' # 20210201
            for k in ContigListDic[FChr]:
                ContigS, ContigE = int(k[0]), int(k[1])
                ConLen = ContigE - ContigS + 1
                FDLen  = Overlapping(int(FStart), int(FEnd), ContigS, ContigE)
                if FDLen != ConLen and FDLen > 0: 
                    FCN = Pos2Str(FChr, ContigS, ContigE)
                    break
                if FDLen == ConLen: Saver = 'pass'
                if FDLen != ConLen and Saver == 'pass': 
                    FCN = Pos2Str(FChr, ContigS, ContigE) 
                    break

            Saver = ''
            for k in ContigListDic[TChr]:
                ContigS, ContigE = int(k[0]), int(k[1])
                ConLen = ContigE - ContigS + 1
                TDLen  = Overlapping(int(TStart), int(TEnd), ContigS, ContigE)
                if TDLen != ConLen and TDLen > 0: 
                    TCN = Pos2Str(TChr, ContigS, ContigE)
                    break
                if TDLen == ConLen: Saver = 'pass'
                if TDLen != ConLen and Saver == 'pass':
                    TCN = Pos2Str(TChr, ContigS, ContigE)
                    break

            PurgeBlockList.append([PurgeFDNo, FChr, FStart, FEnd, TStart, TEnd, TChr, Type, FCN, TCN])
            PurgeNo += 1
        else: 
            PurgeBlockList.append([PurgeFDNo, FChr, FStart, FEnd, TStart, TEnd, TChr, Type])
            PurgeNo += 1

    BlockDic = {} # {BlockName : DL}
    Step = 0
    DGBlockDic = {}
    for i in CactusBlockList + PurgeBlockList:
        Step += 1
        #if (len(CactusBlockList + PurgeBlockList) - Step) % 100 == 0: print('Main7\t' + str(len(CactusBlockList + PurgeBlockList) - Step))
        DL = i
        BlockNo, FChr, FStart, FEnd = DL[0], DL[1], int(DL[2]), int(DL[3])
        BlockDic[BlockNo] = DL[1:]
        if FChr in InContigDic: 
            SubContigDic = InContigDic[FChr] # {Chr:{Loci1:{'FDPos':[Start, End], 'FDBlock':[BLockNo1, BLockNo2]},Loci2:{}} 
            DGSave = 0
            for j in SubContigDic:
                Loci, LL = j, Str2Pos(j)
                CStart, CEnd = int(LL[1]), int(LL[2])
                if CStart <= FEnd and FStart <= CEnd: #1. Narrow: break, Over: continue
                    DGSave = 1
                    InContigDic[FChr][Loci]['FDBlock'].append(BlockNo)
                    if FStart <= CStart and CEnd <= FEnd:
                        InContigDic[FChr][Loci]['FDPos'].append([CStart, CEnd])
                    elif CStart <= FStart and FEnd <= CEnd: # Narrow
                        InContigDic[FChr][Loci]['FDPos'].append([FStart, FEnd])
                        break
                    elif FStart <= CStart and FEnd <= CEnd:
                        InContigDic[FChr][Loci]['FDPos'].append([CStart, FEnd])
                        break
                    elif CStart <= FStart and CEnd <= FEnd:
                        InContigDic[FChr][Loci]['FDPos'].append([FStart, CEnd])
            if DGSave == 0: DGBlockDic[BlockNo] = DL[1:] + [{'Pass':['FDonDGCon'], 'Filtered':[]}] # No ContigLoci in ContigDic by FD between contig (DG)
        else: DGBlockDic[BlockNo] = DL[1:] + [{'Pass':['FDonDGSca'], 'Filtered':[]}] # No Sca in ContigDic by FD on DG
    
    for i in InContigDic:
        Chr, LociDic = i, InContigDic[i]
        for j in LociDic:
            Loci, DD = j, LociDic[j]
            FDPosList, FDBlockList = DD['FDPos'], DD['FDBlock']
            MergedFDPos = M7_OverlapMerge(FDPosList)
            InContigDic[Chr][Loci]['FDPos'] = MergedFDPos
    return InContigDic, BlockDic, PurgeBlockList, DGBlockDic


def M8_ContigFDPercentAndCED(Loci, List):
    CStart, CEnd = int(Loci.split(':')[1].split('-')[0]), int(Loci.split('-')[1])
    CLen = CEnd - CStart + 1
    CFDLen, CEDList = 0, []
    for i in List:
        Start, End = int(i[0]), int(i[1])
        RegionLen = End - Start + 1
        CFDLen += RegionLen
        #CED = 1 - (min(Start - CStart, CEnd - End) / (CLen / 2))
        CED = min(Start - CStart, CEnd - End)
        CEDList.append(CED)
    CFDPercent = CFDLen / CLen
    return CFDPercent, CEDList


def Main8_ContigFDTypeAllocation(ContigDic, InsertSizeCutoff):  # {Chr:{Loci1:{'FDPos':[Start, End], 'FDBlock':[BLockNo1, BLockNo2], 'CFDType':['Whole' or 'Edge','Central'..]'},Loci2:{}}
    InContigDic = ContigDic
    for i in InContigDic:
        Chr, LociDic = i, InContigDic[i]
        for j in LociDic:
            Loci, DD = j, LociDic[j]
            FDPosList = DD['FDPos']
            CFDCED = M8_ContigFDPercentAndCED(Loci, FDPosList)
            FDPercent, CEDList = CFDCED[0], CFDCED[1]
            InContigDic[Chr][Loci]['CFDType'] = []
            if FDPercent >= 0.5: InContigDic[Chr][Loci]['CFDType'].append('Whole')
            elif FDPercent < 0.5:
                for k in CEDList:
                    CED = k
                    if CED <= InsertSizeCutoff: InContigDic[Chr][Loci]['CFDType'].append('Edge')
                    if CED > InsertSizeCutoff: InContigDic[Chr][Loci]['CFDType'].append('Central')
    return InContigDic            


def Main8_2_PurgeTLociTrans(BlockDic, PADic):
    InBlockDic = BlockDic
    for i in InBlockDic:
        BlockNo, DL = i, InBlockDic[i]
        Tool, FChr, FStart, FEnd, TStart, TEnd, TChr, Left = BlockNo.split(' ')[0], DL[0], DL[1], DL[2], DL[3], DL[4], DL[5], DL[6:]
        if Tool == 'Cactus': pass
        else:
            PADL = PADic[BlockNo]
            PATPos = Str2Pos(PADL[1])
            TStart, TEnd = PATPos[1], PATPos[2]
        InBlockDic[i] = [FChr, FStart, FEnd, TStart, TEnd, TChr] + Left  + [{'Filtered':[], 'Pass':[]}]#Conversion to True Alignment
    return InBlockDic


###################################################################################################
#####################################################################################################


### [M9-Sub tool] Pair Connection Counting
def M9_Sub_ReadLogger(FLoci, TLoci, ParsingRegion, BlockNo, RN, Status, JobPath):
    with open(JobPath + 'ReadLogger.txt', 'a') as fwl:
        NowTime = time.ctime(time.time())
        fwl.write('\t'.join([NowTime, FLoci, TLoci, ParsingRegion, BlockNo, RN, Status]) + '\n')


#### USED
def M9_Tool_SuppleReadFilter(ReadList):
    PassReads = []
    if ReadList.count('') > 0: ReadList.remove('')
    for i in ReadList:
        DT = i.split('\t')
        if len(DT) > 1: pass#20210202
        else:               #
            #print(DT)       #
            #print(ReadList) #
            continue        #
        RN, FS = DT[0], DT[1]
        if FSDecimal(FS)[256] == 1: continue
        if FSDecimal(FS)[1024] == 1: continue
        PassReads.append(i)
    return PassReads

def M9_ColChr(PairChr, Chr):
    if PairChr == '=': return Chr
    else: return PairChr

#### USED 
def M9_OverlapExtract(List1, List2):#c #Only Overlap Region Extractor
    SortList1, SortList2 = [], []
    for i in List1:
        Start, End = int(i[0]), int(i[1])
        SortList1.append([Start, End])
    for i in List2:
        Start, End = int(i[0]), int(i[1])
        SortList2.append([Start, End])
    SortList1.sort()
    SortList2.sort()
    OverlapList = []
    for i in SortList1:
        Start1, End1 = i[0], i[1]
        for j in SortList2:
            Start2, End2 = j[0], j[1]
            if Start1 <= End2 and Start2 <= End1:
                OverlapPos = [Start1, Start2, End1, End2]
                OverlapPos.remove(max(OverlapPos))
                OverlapPos.remove(min(OverlapPos))
                OverlapList.append(OverlapPos)
    return OverlapList

def M9_PosClustering(DoubleList, ClusteringSize): # Clustering by distance
    List = []
    for i in DoubleList:
        Start, End = int(i[0]), int(i[1])
        List.append([Start, End])
    List.sort()
    Clustering = []
    for i in List:
        if len(Clustering) == 0: Clustering.append(i)
        else:
            Start, End, PrevStart, PrevEnd = i[0], i[1], MergeList[-1][0], MergeList[-1][1]
            if Start - ClusteringSize <= PrevEnd: Clustering[-1] = [PrevStart, End]
            else: MergeList.append(i)
    return Clustering

#### Used
def M9_BlockClustering_Stepwise(BlockList, BlockDic, ClusteringSize):#c # 1. FD sort; 2. Clustering: both same Chr?; both distance < Flanking? + Both TDirection
    FPosDic = {} #
    for i in BlockList:
        BlockNo, DL = i, BlockDic[i]
        FStart, FEnd, TStart, TEnd, TChr = int(DL[1]), int(DL[2]), int(DL[3]), int(DL[4]), DL[5] # Chr will be same
        FPosDic[(FStart, FEnd)] = {'BlockNo': BlockNo, 'TPos':[TChr, TStart, TEnd]}
    Keys = list(FPosDic.keys())
    Keys.sort()
    NewKeys = [] # {NewKey: Included BlockNOs}
    for i in Keys:
        FStart, FEnd = i[0], i[1]
        Inner = 0
        for j in Keys:
            Start, End = j[0], j[1]
            if Start <= FStart and FEnd <= End and (FStart, FEnd) != (Start, End):
                Inner = 1
                del FPosDic[i] # del redundancy
                continue
        if Inner == 0: NewKeys.append(i) # only add no redundant for i

    FCluster, TCluster, ClusterDic, BlockList  = [], [], {}, [] # TRegion, FRegion, BlockNo
    for i in NewKeys: # Inner blocks are not analysed.
        if len(FCluster) == 0: 
            FStart, FEnd, TChr, TStart, TEnd, BlockNo = i[0], i[1], FPosDic[i]['TPos'][0], FPosDic[i]['TPos'][1], FPosDic[i]['TPos'][2], FPosDic[i]['BlockNo']
            FCluster = [FStart, FEnd]
            TCluster = [TChr, TStart, TEnd]
            BlockList.append(BlockNo)
        else:
            Start, End, TChr, TStart, TEnd, BlockNo = i[0], i[1], FPosDic[i]['TPos'][0], FPosDic[i]['TPos'][1], FPosDic[i]['TPos'][2], FPosDic[i]['BlockNo']
            PrevStart, PrevEnd = FCluster[0], FCluster[1]
            PrevTChr, PrevTStart, PrevTEnd = TCluster[0], TCluster[1], TCluster[2]
            if TChr == PrevTChr and Start - ClusteringSize <= PrevEnd: # FClustering
                if TStart - ClusteringSize <= PrevTEnd or PrevTStart - ClusteringSize <= TEnd: # TClustering for two direction
                    FCluster = [PrevStart, End]
                    TCluster = [TChr, min([TStart, TEnd, PrevTStart, PrevTEnd]), max([TStart, TEnd, PrevTStart, PrevTEnd])]
                    BlockList.append(BlockNo)
            else: #End Clustering
                ClusterDic[(TCluster[0], TCluster[1], TCluster[2])] = {'FDPos':FCluster, 'BlockNo':BlockList} # Cluster Saving
                FCluster, TCluster, BlockList = [], [], [] # reset saved data
                FCluster, TCluster = [Start, End], [TChr, TStart, TEnd]
                BlockList.append(BlockNo)
    ClusterDic[(TCluster[0], TCluster[1], TCluster[2])] = {'FDPos':FCluster, 'BlockNo':BlockList} # Last cluster saving
    #M9_Sub_ReadLogger(str(TCluster), str(FCluster), str(BlockList), '-', '-', 'Clustered')
    return ClusterDic # {(ClusterChr, ClsuteredTStart, ClusteredTEnd): {'FDPos':[ClusteredFDStart, ClusteredFDEnd], 'BlockNo':[Block1, Block2..]}


#### USED  
def PurgeContigLociFinding(ContigSubDic, PurgeHaplotigLoci):
    ContigList = list(ContigSubDic.keys())
    Start, End = int(PurgeHaplotigLoci[0]), int(PurgeHaplotigLoci[1])
    Nearest = ['init', 100000000000]
    for i in ContigList:
        Pos = Str2Pos(i)
        CChr, CStart, CEnd = Pos[0], int(Pos[1]), int(Pos[2])
        if CStart <= End and Start <= CEnd: return i
        else:
            if Nearest[1] > abs(CStart - Start): Nearest = [i, abs(CStart - Start)]
    return Nearest[0]

#### USED 
def M9_EdgePosbyCED(FDBlockList, BlockDic, ContigLoci, FlankingSize, InsertSizeCutoff, BamPath, JobPath, ContigDic): # Clustering embedded
    TContigList, FDPosCED, TPosCED = [], {}, {}
    for k in FDBlockList: # TLoci Deciding by Largest FD 
        BlockNo, DL, Tool = k, BlockDic[k], k.split(' ')[0]
        FChr, FStart, FEnd, TStart, TEnd, TChr = DL[0], int(DL[1]), int(DL[2]), int(DL[3]), int(DL[4]), DL[5]
        FLociStr, TLociStr = Pos2Str(FChr, FStart, FEnd), Pos2Str(TChr, TStart, TEnd)
        if Tool == 'Cactus': TLociContig = DL[17]
        else: 
            if len(DL) < 9: TLociContig = PurgeContigLociFinding(ContigDic[TChr], [TStart, TEnd])  #20210202
            else: TLociContig = DL[8] 
            #if len(DL) <9: # mis-contiging by non-included assembly gap (e.g. Ns < 10 in VGP himmingbird: CM012116:112728208-113254958 (112737392 will be end of contig: actually not 'Edge' but 'WCD') 
        FCED, TCED = CEDCalculation(ContigLoci, FLociStr), CEDCalculation(TLociContig, TLociStr)
        TContigList.append(TLociContig)
        FDPosCED[FCED], TPosCED[TCED] = [FChr, FStart, FEnd], [TChr, TStart, TEnd]
   
    if len(set(TContigList)) == 1: # Colinear: Clustered
        FInnerFD, TInner, FOutter, TOutter = FDPosCED[max(FDPosCED.keys())], TPosCED[max(TPosCED.keys())], FDPosCED[min(FDPosCED.keys())], TPosCED[min(TPosCED.keys())]
        FChr, FStart, FEnd, TChr, TStart, TEnd = FInnerFD[0], FInnerFD[1], FInnerFD[2], TInner[0], TInner[1], TInner[2]
        FOStart, FOEnd, TOStart, TOEnd = FOutter[1], FOutter[2], TOutter[1], TOutter[2]
        CEDFLociStr, CEDTLociStr = Pos2Str(FChr, min(FStart, FOStart), max(FEnd, FOEnd)), Pos2Str(TChr, min(TStart, TOStart), max(TEnd, TOEnd)) 
        if CEDCalculation_Sol(ContigLoci, FStart) > CEDCalculation_Sol(ContigLoci, FEnd): FFlanking = Pos2Str(FChr, FStart - FlankingSize, FStart - 1) # F Flanking Direction
        else: FFlanking = Pos2Str(FChr, FEnd + 1, FEnd + FlankingSize)
        if CEDCalculation_Sol(TContigList[0], TStart) > CEDCalculation_Sol(TContigList[0], TEnd): TFlanking = Pos2Str(TChr, TStart - FlankingSize, TStart - 1) # T Flanking Direction
        else: TFlanking = Pos2Str(TChr, TEnd + 1, TEnd + FlankingSize)
        #print(CEDFLociStr + '\t' + CEDTLociStr + '\t' + FFlanking + '\t' + ContigLoci)
        #print(CEDTLociStr + '\t' + CEDFLociStr + '\t' + TFlanking + '\t' + str(TContigList))
        FDiscDic = M9_Sub_PairConnectionOVLP4NormalFinal(CEDFLociStr, CEDTLociStr, InsertSizeCutoff, FlankingSize, BamPath, 'Edge', FFlanking, JobPath)
        TDiscDic = M9_Sub_PairConnectionOVLP4NormalFinal(CEDTLociStr, CEDFLociStr, InsertSizeCutoff, FlankingSize, BamPath, 'Edge', TFlanking, JobPath)
        if FDiscDic['OverlapSize'] + TDiscDic['OverlapSize'] > 0: return FDBlockList, [CEDFLociStr, CEDTLociStr, CEDTLociStr, CEDFLociStr], [FDiscDic, TDiscDic] # return FDBlock, Evidence
        else: return 'nonFD', [CEDFLociStr, CEDTLociStr, CEDTLociStr, CEDFLociStr], [FDiscDic, TDiscDic]
    else: # non-Colinnear
        Index = 0
        DiscBlockList, DiscBlockDic = [], {}
        for l in FDBlockList:
            BlockNo, DL, Tool = l, BlockDic[l], l.split(' ')[0]
            FChr, FStart, FEnd, TStart, TEnd, TChr = DL[0], int(DL[1]), int(DL[2]), int(DL[3]), int(DL[4]), DL[5]
            FLoci, TLoci = Pos2Str(FChr, FStart, FEnd), Pos2Str(TChr, TStart, TEnd)
            if Tool == 'Cactus': TLociContig = DL[17]
            else: 
                if len(DL) < 9: TLociContig = PurgeContigLociFinding(ContigDic[TChr], [TStart, TEnd])  #20210202
                else: TLociContig = DL[8]
            if CEDCalculation_Sol(ContigLoci, FStart) > CEDCalculation_Sol(ContigLoci, FEnd): FFlanking = Pos2Str(FChr, FStart - FlankingSize, FStart - 1) # F Flanking Direction
            else: FFlanking = Pos2Str(FChr, FEnd + 1, FEnd + FlankingSize)
            if CEDCalculation_Sol(TContigList[Index], TStart) > CEDCalculation_Sol(TContigList[Index], TEnd): TFlanking = Pos2Str(TChr, TStart - FlankingSize, TStart - 1) # T Falnking Direction
            else: TFlanking = Pos2Str(TChr, TEnd + 1, TEnd + FlankingSize)
            FDiscDic = M9_Sub_PairConnectionOVLP4NormalFinal(FLoci, TLoci, InsertSizeCutoff, FlankingSize, BamPath, 'Edge', FFlanking, JobPath)
            TDiscDic = M9_Sub_PairConnectionOVLP4NormalFinal(TLoci, FLoci, InsertSizeCutoff, FlankingSize, BamPath, 'Edge', TFlanking, JobPath)
            if FDiscDic['OverlapSize'] + TDiscDic['OverlapSize'] > 0: 
                DiscBlockList.append(BlockNo)
                DiscBlockDic[BlockNo] = [FDiscDic, TDiscDic]
            Index += 1
        if len(DiscBlockList) == 0: return 'nonFD', ['Seperately'], []
        else: return DiscBlockList, ['Seperately'], DiscBlockDic
                            
####USED 
def M9_Sub_PairConnectionOVLP4NormalFinal(TLoci, FLoci, InsertSizeCutoff, FlankingSize, BamPath, FlankingStat, O1ParsingLoci, JobPath): #c
    TPos, FPos, DiscDic, DiscReadRegion, NormalReadRegion = Str2Pos(TLoci), Str2Pos(FLoci), {'DisRead':[], 'NormRead':[], 'OverlapSize':0}, [], [] # 1. Read Parsing on Flanking of TLoci; 2. Region Merging 4 Disc Read to FLoci and Normal Read to TLoci; 3. Overlap between 2 region
    TChr, TStart, TEnd, FChr, FStart, FEnd = TPos[0], int(TPos[1]), int(TPos[2]), FPos[0], int(FPos[1]), int(FPos[2])
    if FlankingStat == "Both": # 4 WCD and Central 
        ParsingLoci1, ParsingLoci2 = Pos2Str(TChr, TEnd + 1, TEnd + FlankingSize), Pos2Str(TChr, TStart - FlankingSize, TStart - 1) # Both Direction of Falnking of TLoci
        if TStart - 1 == 0: ParsingLoci2 = Pos2Str(TChr,10000000000,10000000000) # 20210204
        ParsingRead, TLociRead, FLociRead = ReadParser(ParsingLoci1, BamPath) + ReadParser(ParsingLoci2, BamPath), ReadParser(TLoci, BamPath), ReadParser(FLoci, BamPath)      
    elif FlankingStat == "Edge": # 4 WED
        ParsingRead, TLociRead, FLociRead = ReadParser(O1ParsingLoci, BamPath), ReadParser(TLoci, BamPath), ReadParser(FLoci, BamPath) # One Direction of Falnking of TLoci
    ParsingRead, TLociRead, FLociRead = M9_Tool_SuppleReadFilter(ParsingRead), M9_Tool_SuppleReadFilter(TLociRead), M9_Tool_SuppleReadFilter(FLociRead)
    TLociReadDic, FLociReadDic = {}, {}
    for i in TLociRead:
        DT = i.split('\t')
        ReadName = DT[0]
        TLociReadDic[ReadName] = DT
    for i in FLociRead:
        DT = i.split('\t')
        ReadName = DT[0]
        FLociReadDic[ReadName] = DT
    for i in ParsingRead:
        DT = i.split('\t')
        ReadName, FlagStat, ReadChr, ReadStart, MappingQuality, CigarString, PairReadChr, PairReadPos, InsertSize = DT[0], DT[1], DT[2], int(DT[3]), DT[4], DT[5], M9_ColChr(DT[6], DT[2]), int(DT[7]), abs(int(DT[8]))
        if CigarString == '*': continue
        if ReadName in TLociReadDic and int(TLociReadDic[ReadName][3]) == PairReadPos: # Find Paired Read on normal status
            if InsertSize <= InsertSizeCutoff and PairReadChr == TChr:
                DiscDic['NormRead'].append(ReadName)
                ReadEnd = ReadStart + Cigar4Len(CigarString) - 1
                NormalReadRegion.append([ReadStart, ReadEnd])
                #if FlankingStat == "Both": M9_Sub_ReadLogger(TLoci, FLoci, ParsingLoci1 + '|' + ParsingLoci2, '-', ReadName, 'Normal', LogPath)
                #elif FlankingStat == "Edge": M9_Sub_ReadLogger(TLoci, FLoci, O1ParsingLoci, '-', ReadName, 'Normal', LogPath)
        elif ReadName in FLociReadDic and int(FLociReadDic[ReadName][3]) == PairReadPos: # Find paired Read on Discordant Status
            if (InsertSize > InsertSizeCutoff and PairReadChr == FChr) or (InsertSize == 0 and PairReadChr == FChr and FChr != TChr):
                DiscDic['DisRead'].append(ReadName)
                ReadEnd = ReadStart + Cigar4Len(CigarString) - 1
                DiscReadRegion.append([ReadStart, ReadEnd])
                #if FlankingStat == "Both": M9_Sub_ReadLogger(TLoci, FLoci, ParsingLoci1 + '|' + ParsingLoci2, '-', ReadName, 'Discord', LogPath)
                #elif FlankingStat == "Edge": M9_Sub_ReadLogger(TLoci, FLoci, O1ParsingLoci, '-', ReadName, 'Discord', LogPath)
        else: continue
    NormalMergeRegion, DiscMergeRegion = M7_OverlapMerge(NormalReadRegion), M7_OverlapMerge(DiscReadRegion) # Inner Region merge
    OverlapRegion = M9_OverlapExtract(NormalMergeRegion, DiscMergeRegion) # Find Overlap between Normal and Discordant mapped read region
    DiscDic['OverlapSize'] = PosLenSum(OverlapRegion)
    #M9_Sub_ReadLogger(TLoci, FLoci, '-', str(NormalMergeRegion), str(DiscMergeRegion), str(DiscDic['OverlapSize']), LogPath)
    return DiscDic


def M9_Sub_PairConnection(ParsingLoci, CheckLoci, BamPath): #1. ReadParsing; 2. Discordance Check
    ParsingReads, CheckReads = ReadParser(ParsingLoci, BamPath), ReadParser(CheckLoci, BamPath)
    ParsingReads, CheckReads = M9_Tool_SuppleReadFilter(ParsingReads), M9_Tool_SuppleReadFilter(CheckReads)
    DiscRead, CheckDic = 0, {}
    for i in CheckReads:
        DT = i.split('\t')
        RN = DT[0]
        CheckDic[RN] = i
    for i in ParsingReads:
        DT = i.split('\t')
        ReadName, FlagStat, ReadChr, ReadStart, MappingQuality, CigarString, PairReadChr, PairReadPos, InsertSize = DT[0], DT[1], DT[2], int(DT[3]), DT[4], DT[5], M9_ColChr(DT[6], DT[2]), int(DT[7]), abs(int(DT[8]))
        if InsertSize > InsertSizeCutoff and PairReadChr == TChr:    
            if ReadName in CheckDic:
                PairReadStart = int(CheckDic[ReadName].split('\t')[3])
                if PairReadPos == PairReadStart: DiscRead += 1
    return DiscRead

        
def M9_Misassigned_Central_Pacbio(Loci1, FlankingSize, FlankingAvoid): # 3 spot pacbio read connnection checking tool ( mis-assignment test )
    FLoci, Pos = Loci1, Str2Pos(Loci1)
    FChr, Start, End = Pos[0], Pos[1], Pos[2]
    ForStart, ForEnd = int(Start) - FlankingSize, int(Start) - FlankingSize + FlankingAvoid - 1 
    RevStart, RevEnd = int(End) + FlankingSize, int(End) + FlankingSize - FlankingAvoid + 1
    ForSpotLoci, RevSpotLoci, FDSpotLoci = Pos2Str(FChr, ForStart, ForEnd), Pos2Str(FChr, RevStart, RevEnd), Pos2Str(FChr, Start, End)
    ReadDic = {}
    ForSpotReads = sp.getoutput("samtools view {0} {1} | awk '{{print $1}}' ".format(PacBam, ForSpotLoci)).split('\n')
    RevSpotReads = sp.getoutput("samtools view {0} {1} | awk '{{print $1}}' ".format(PacBam, RevSpotLoci)).split('\n')
    FDSpotReads = sp.getoutput("samtools view {0} {1} | awk '{{print $1}}' ".format(PacBam, FDSpotLoci)).split('\n')
    for i in ForSpotReads: ReadDic[i] = ['For']
    for i in RevSpotReads:
        if i in ReadDic: ReadDic[i].append('Rev')
        else: ReadDic[i] = ['Rev']
    for i in FDSpotReads:
        if i in ReadDic: ReadDic[i].append('FD')
        else: ReadDic[i] = ['FD']

    Connection = 0
    for i in ReadDic:
        if len(set(ReadDic[i])) == 3:
            Connection += 1
            M9_Sub_ReadLogger(FLoci, FLoci, i)
    if Connection == 0: return "FD", Connection
    else: return "nonFD", Connection


###############################################################################################################################################
###########                                           [OUTTER: FDTD Reads Analysis Tool]                                            ###########

def DepthCalling(Loci, BamPath):
    LL = Str2Pos(Loci)
    Chr, Start, End = LL[0], LL[1], LL[2]
    Len = int(End) - int(Start) + 1
    TotalDepth = sp.getoutput("samtools depth -r {0} {1} | awk 'BEGIN{{Depth = 0}} {{Depth += $3}} END{{print Depth}}'".format(Loci, BamPath))
    MeanDepth = int(TotalDepth) / Len
    return MeanDepth

### [OUTTER-Main]
def M_TDFDLociParsingbyDepth(M_Keys): 
    InFDDic = MainDic1
    TDVGP, TDPre, FDPre, FDVGP, M_Log = [], [], [], [], [] #[Loci1, Loci2..]
    Step = 0
    for i in M_Keys:
        Step += 1
        #if (len(M_Keys) - Step) % 10 == 0: print(len(M_Keys) - Step)

        MafNo, DD = i, InFDDic[i]
        NoSeqList, DD2 = DD['NoSeq'], DD['Dat']
        if NoSeqList[0] == 2 and NoSeqList[2] == 2 and len(TDVGP) < (MQSampling / Cores) and len(TDPre) < (MQSampling / Cores):
            VGP1Chr, VGP1Start, VGP1End, VGP2Chr, VGP2Start, VGP2End = DD2['R1']['Chr'], DD2['R1']['Start'], DD2['R1']['End'], DD2['R2']['Chr'], DD2['R2']['Start'], DD2['R2']['End']
            Pre1Chr, Pre1Start, Pre1End, Pre2Chr, Pre2Start, Pre2End = DD2['T1']['Chr'], DD2['T1']['Start'], DD2['T1']['End'], DD2['T2']['Chr'], DD2['T2']['Start'], DD2['T2']['End']
            VGP1Loci, VGP2Loci, VGP1Len, VGP2Len = VGP1Chr + ':' + VGP1Start + '-' + VGP1End, VGP2Chr + ':' + VGP2Start + '-' + VGP2End, int(VGP1End)-int(VGP1Start) + 1, int(VGP2End)-int(VGP2Start) + 1
            Pre1Loci, Pre2Loci, Pre1Len, Pre2Len = Pre1Chr + ':' + Pre1Start + '-' + Pre1End, Pre2Chr + ':' + Pre2Start + '-' + Pre2End, int(Pre1End)-int(Pre1Start) + 1, int(Pre2End)-int(Pre2Start) + 1
            VGP1Depth, VGP2Depth = DepthCalling(VGP1Loci, VGPBamPath), DepthCalling(VGP2Loci, VGPBamPath)
            Pre1Depth, Pre2Depth = DepthCalling(Pre1Loci, PreBamPath), DepthCalling(Pre2Loci, PreBamPath)
            fLog.write(str([MafNo, VGP1Loci, VGP2Loci, Pre1Loci, Pre2Loci, VGP1Depth, VGP2Depth, Pre1Depth, Pre2Depth]) + '\n')
            if VGPCutoff < VGP1Depth and VGPCutoff < VGP2Depth and PreCutoff < Pre1Depth and PreCutoff < Pre2Depth:
                TDVGP.append([VGP1Loci, VGP2Loci])
                TDPre.append([Pre1Loci, Pre2Loci])
        elif NoSeqList[0] == 2 and NoSeqList[2] == 1 and len(FDVGP) < (MQSampling / Cores):
            VGP1Chr, VGP1Start, VGP1End, VGP2Chr, VGP2Start, VGP2End = DD2['R1']['Chr'], DD2['R1']['Start'], DD2['R1']['End'], DD2['R2']['Chr'], DD2['R2']['Start'], DD2['R2']['End']
            VGP1Loci, VGP2Loci = VGP1Chr + ':' + VGP1Start + '-' + VGP1End, VGP2Chr + ':' + VGP2Start + '-' + VGP2End
            VGP1Depth, VGP2Depth = DepthCalling(VGP1Loci, VGPBamPath), DepthCalling(VGP2Loci, VGPBamPath)
            fLog.write(str([MafNo, VGP1Loci, VGP2Loci, VGP1Depth, VGP2Depth]) + '\n')
            if VGPCutoff >= VGP1Depth and VGPCutoff >= VGP2Depth: FDVGP.append([VGP1Loci, VGP2Loci])

        elif NoSeqList[0] == 1 and NoSeqList[2] == 2 and len(FDPre) < (MQSampling / Cores):
            Pre1Chr, Pre1Start, Pre1End, Pre2Chr, Pre2Start, Pre2End = DD2['T1']['Chr'], DD2['T1']['Start'], DD2['T1']['End'], DD2['T2']['Chr'], DD2['T2']['Start'], DD2['T2']['End']
            Pre1Loci, Pre2Loci = Pre1Chr + ':' + Pre1Start + '-' + Pre1End, Pre2Chr + ':' + Pre2Start + '-' + Pre2End
            Pre1Depth, Pre2Depth = DepthCalling(Pre1Loci, PreBamPath), DepthCalling(Pre2Loci, PreBamPath)
            fLog.write(str([MafNo, Pre1Loci, Pre2Loci, Pre1Depth, Pre2Depth]) + '\n')
            if PreCutoff >= Pre1Depth and PreCutoff >= Pre2Depth: FDPre.append([Pre1Loci, Pre2Loci])
    PID = str(os.getpid())
    with open(JobPath + 'Duplication/' + PID + 'Duplication.list', 'wb') as fw:
        pickle.dump([TDVGP, TDPre, FDPre, FDVGP], fw)


### [OUTTER-Pipeline]
def M_TDFDLociParsingbyDepth_Operator(TopFDDic, Cores):
    M_FDDic = MainDic1
    os.system("mkdir {0}Duplication/".format(LogPath))
    KeyList = list(M_FDDic.keys())
    random.shuffle(KeyList)
    KeyList = KeyList[0:100000]  #####Subset
    input_keys = np.array_split(KeyList, Cores)
    parmap.map(M_TDFDLociParsingbyDepth, input_keys, pm_pbar = True)
    TDVGP, TDPre, FDPre, FDVGP = [], [], [], []
    file_list = os.listdir(JobPath + 'Duplication/')
    for i in file_list:
        with open(JobPath + 'Duplication/' + i, 'rb') as lr:
            LoadList = pickle.load(lr)
            TDVGP += LoadList[0]
            TDPre += LoadList[1]
            FDPre += LoadList[2]
            FDVGP += LoadList[3]
    
    with open(JobPath + 'Merged_Duplication.list','wb') as fw: pickle.dump([TDVGP, TDPre, FDPre, FDVGP], fw)
    return [TDVGP, TDPre, FDPre, FDVGP]

### [OUTTER-Loader]
def M_TDFDLociParsingbyDepth_Loader():
    with open(JobPath + 'Merged_Duplication.list','rb') as fr: 
        List = pickle.load(fr)
        return List 

### [OUTTER-Saver]
def MQList2File(List, Name):
    with open(JobPath + 'MQResult.txt', 'a') as fw:
        for i in List: # RegionTotal, Both60 Discord, Both60 normal, Both0 Discord, Both0 Normal
            #print(i)
            fw.write(Name + '\t' + '\t'.join(map(str, i)) + '\n')

### [OUTTER-Main2]
def ReadMQParsing(LociList, BamPath):
    ReadInfList = [] # RegionTotal, Both60 Discord, Both60 normal, Both0 Discord, Both0 Normal; <20 == 0, >=20 == 60
    Step = 0
    for i in LociList[0:100]:
        #print(len(LociList) - Step)
        Step += 1
        #print(Step)
        Loci1, Loci2 = i[0], i[1]
        Chr1, Start1, End1, Chr2, Start2, End2 = Loci1.split(':')[0], Loci1.split(':')[1].split('-')[0], Loci1.split('-')[1], Loci2.split(':')[0], Loci2.split(':')[1].split('-')[0], Loci2.split('-')[1]
        if Chr1 == Chr2: 
            if Overlapping(int(Start1) - 550, int(End1) + 550, int(Start2) - 550, int(End2) + 550) > 0: continue # Overlap Controller
        FLoci1, FLoci2 = Chr1 +':' + str(int(Start1) - 550) + '-' + str(int(End1) + 550), Chr2 + ':' + str(int(Start2) - 550) + '-' + str(int(End2) + 550)
        Reads1 = sp.getoutput("samtools view {1} {0} | awk '{{print $1\"\t\"$3\"\t\"$5\"\t\"$7\"\t\"$8\"\t\"$9\"\t\"$4\"\t\"$2}}'".format(Loci1, BamPath)).split('\n')
        Reads2 = sp.getoutput("samtools view {1} {0} | awk '{{print $1\"\t\"$3\"\t\"$5\"\t\"$7\"\t\"$8\"\t\"$9\"\t\"$4\"\t\"$2}}'".format(Loci2, BamPath)).split('\n')
        ReadsF = sp.getoutput("samtools view {2} {0} {1}| awk '{{print $1\"\t\"$3\"\t\"$5\"\t\"$7\"\t\"$8\"\t\"$9\"\t\"$4\"\t\"$2}}'".format(FLoci1, FLoci2, BamPath)).split('\n')
        MQDic = {} # {ReadsName:{Chr-StartLoci:PMQ}}
        FSDic = {} # {Readname:{Chr-StartLoci:FlagStat}}
        for j in ReadsF:
            DT = j.split('\t')
            ReadName, Chr, Loc, MQ, FS = DT[0], DT[1], DT[6], DT[2], DT[7]
            Loc = Chr + ':' + Loc
            if ReadName in MQDic: pass
            else: MQDic[ReadName], FSDic[ReadName] = {}, {}
            if Loc in MQDic[ReadName]: 
                del MQDic[ReadName][Loc]
                del FSDic[ReadName][Loc]
            else: MQDic[ReadName][Loc], FSDic[ReadName][Loc] = MQ, FS
        #print("ReadParsing") 
        MaxDepth1 = sp.getoutput("samtools depth -r {0} {1} | awk 'BEGIN{{Max = 0}} $3 > Max {{Max = $3}} END{{print Max}}'".format(Loci1, BamPath))
        MaxDepth2 = sp.getoutput("samtools depth -r {0} {1} | awk 'BEGIN{{Max = 0}} $3 > Max {{Max = $3}} END{{print Max}}'".format(Loci2, BamPath))
        #print('DepthCal')
        if int(MaxDepth1) > 200 or int(MaxDepth2) > 200: continue
        ReadType = [0,0,0,0,0, Loci1, Loci2]
        if Reads1.count('') > 0: Reads1.remove('')
        if Reads2.count('') > 0: Reads2.remove('')
        for j in Reads1:
            DT = j.split('\t')
            RN, Chr, MQ, PairChr, PairLoc, InsertSize, Loc, FS = DT[0], DT[1], DT[2], DT[3], DT[4], DT[5], DT[6], DT[7]
            if PairChr == '=' and Loc == PairLoc: continue
            if PairChr == '=': 
                PairChr = Chr
                PairLoc2 = PairChr + ':' + PairLoc
            else: 
                if PairChr == '*': continue
                else: PairLoc2, InsertSize = PairChr + ':' + PairLoc, 9999
            if PairLoc2 in MQDic[RN]: PMQ, PFS = MQDic[RN][PairLoc2], FSDic[RN][PairLoc2]
            else: continue
            MappingDirection = FFRRFinder(FSDecimal(FS), FSDecimal(PFS))
            if MappingDirection == 'Diff': pass
            else: continue
            ReadType[0] += 1
            if MQCutoff <= int(MQ) and MQCutoff <= int(PMQ): #60, 60
                fLog.write(str([RN, Chr, Start1, Loc, End1, Chr2, PairChr, Start2, PairLoc, End2, InsertSize, MQ, PMQ, FS, PFS, MappingDirection]) + '\n')
                if InsertSizeCutoff < abs(int(InsertSize)) and PairChr == Chr2 and (int(Start2) - 550) <= int(PairLoc) <= (int(End2) + 550):  
                    ReadType[1] += 1
                    fLog.write(str([RN, Chr, Start1, Loc, End1, Chr2, PairChr, Start2, PairLoc, End2, InsertSize, MQ, PMQ, FS, PFS, MappingDirection]) + '\tReadType1' + '\n')
                else:  ReadType[2] += 1
            elif MultiMQCutoff >= int(MQ) and MultiMQCutoff >= int(MQ): #0, 0
                fLog.write(str([RN, Chr, Start1, Loc, End1, Chr2, PairChr, Start2, PairLoc, End2, InsertSize, MQ, PMQ, FS, PFS, MappingDirection]) + '\n')
                if InsertSizeCutoff < abs(int(InsertSize)) and PairChr == Chr2 and (int(Start2) - 550) <= int(PairLoc) <= (int(End2) + 550):  ReadType[3] += 1
                else:  ReadType[4] += 1
        for j in Reads2:
            DT = j.split('\t')
            RN, Chr, MQ, PairChr, PairLoc, InsertSize, Loc, FS = DT[0], DT[1], DT[2], DT[3], DT[4], DT[5], DT[6], DT[7]
            if PairChr == '=' and Loc == PairLoc: continue
            if PairChr == '=': 
                PairChr = Chr
                PairLoc2 = PairChr + ':' + PairLoc
            else: 
                if PairChr == '*': continue
                else: PairLoc2, InsertSize = PairChr + ':' + PairLoc, 9999
            if PairLoc2 in MQDic[RN]: PMQ, PFS = MQDic[RN][PairLoc2], FSDic[RN][PairLoc2]
            else: continue
            MappingDirection = FFRRFinder(FSDecimal(FS), FSDecimal(PFS))
            if MappingDirection == 'Diff': pass
            else: continue
            ReadType[0] += 1
            if MQCutoff <= int(MQ) and MQCutoff <= int(PMQ): #60, 60
                fLog.write(str([RN, Chr, Start2, Loc, End2, Chr1, PairChr, Start1, PairLoc, End1, InsertSize, MQ, PMQ, MappingDirection]) + '\n')
                if InsertSizeCutoff < abs(int(InsertSize)) and PairChr == Chr1 and int(Start1) - 550 <= int(PairLoc) <= int(End1) + 550:  
                    ReadType[1] += 1
                    fLog.write(str([RN, Chr, Start2, Loc, End2, Chr1, PairChr, Start1, PairLoc, End1, InsertSize, MQ, PMQ, MappingDirection]) + '\tReadType1' + '\n')
                else:  ReadType[2] += 1
            elif MultiMQCutoff >= int(MQ) and MultiMQCutoff >= int(MQ): #0, 0
                fLog.write(str([RN, Chr, Start2, Loc, End2, Chr1, PairChr, Start1, PairLoc, End1, InsertSize, MQ, PMQ, MappingDirection]) + '\n')
                if InsertSizeCutoff < abs(int(InsertSize)) and PairChr == Chr1 and int(Start1) - 550 <= int(PairLoc) <= int(End1) + 550:  ReadType[3] += 1
                else:  ReadType[4] += 1
        ReadInfList.append(ReadType)
    return ReadInfList

###########                                                                                                                                 #############
#########################################################################################################################################################
def M9_MultiProc(ContigDic, BlockDic, InsertSizeCutoff, FlankingSize, ClusteringSize, BamPath, SECutoff, Cores, JobPath, New_DepthDic):
    Keys, DataLen = [], 0 
    for i in ContigDic: 
        Keys += list(ContigDic[i].keys())
        for j in ContigDic[i]:
            DD = ContigDic[i][j]
            if len(DD['FDPos']) > 0: DataLen += 1

    random.shuffle(Keys)
    Calling_Keys = np.array_split(Keys, Cores)
    Result = parmap.map(Main9_ContigFDDiscordance, Calling_Keys, ContigDic, BlockDic, InsertSizeCutoff, FlankingSize, ClusteringSize, BamPath, SECutoff, 'Multi', JobPath, New_DepthDic, pm_pbar = True)
    MBlockDic, MContigDic, ResultAmount = {}, ContigDic, 0
    for i in Result: # zero to sum
        BD, CD = i[0], i[1]
        for j in BD:
            BlockNo, DL = j, BD[j]
            FilteredContig, PassedContig = len(DL[-1]['Filtered']), len(DL[-1]['Pass'])
            if FilteredContig + PassedContig != 0: # Data is.
                if BlockNo in MBlockDic:
                    MBlockDic[BlockNo][-1]['Filtered'] += DL[-1]['Filtered']
                    MBlockDic[BlockNo][-1]['Pass'] += DL[-1]['Pass']
                else: MBlockDic[BlockNo] = DL
            else: pass # not analyzed in multi-process

        for j in CD: # Changing
            Chr, LociDic = j, CD[j]
            for k in LociDic:
                Loci, DD = k, LociDic[k]
                if 'Pass' in DD: 
                    MContigDic[Chr][Loci] = DD # Analyzed in Multi-processing
                    ResultAmount += 1
                else: pass # not analyzed in multi-process

    return MBlockDic, MContigDic, DataLen, ResultAmount 


def Main9_ContigFDDiscordance(Keys, ContigDic, BlockDic, InsertSizeCutoff, FlankingSize, ClusteringSize, BamPath, SECutoff, MultiProc, JobPath, New_DepthDic): #(ContigDic, BlockDic, BamPath, SECutoff, PADic)
    InContigDic, InBlockDic = ContigDic, BlockDic
    Step = 0
    if MultiProc == 'Multi': Keys = Keys                        ###
    elif MultiProc == 'Single': Keys = list(InContigDic.keys()) ### dosen't work now
    for i in Keys:                                              ###
        #if (len(Keys) - Step) % 10 == 0: print(len(Keys) - Step)
        Step += 1                                             
        ContigLoci, Chr = i, Str2Pos(i)[0]
        DD = InContigDic[Chr][ContigLoci]
        if len(DD['FDBlock']) == 0: continue # more fast than below?
        InContigDic[Chr][ContigLoci]['Pass'] = []
        InContigDic[Chr][ContigLoci]['Filtered'] = []
        InContigDic[Chr][ContigLoci]['Cluster'] = {}
        FDPosList, FDBlockList, CFDType = DD['FDPos'], DD['FDBlock'], list(set(DD['CFDType'])) # FDPos = merged FD Pos
        #if len(FDPosList) == 0: continue
        FDSumLen, FDDepths = 0, 0
        for k in FDPosList: # SE Check
            FDStart, FDEnd = k[0], k[1]
            FDLen = FDEnd - FDStart + 1
            FDLoci = Pos2Str(Chr, FDStart, FDEnd)
            if FDLoci in New_DepthDic: FDDepth = float(New_DepthDic[FDLoci])
            else: FDDepth = DepthCalling(FDLoci, BamPath)
            FDSumLen += FDLen
            FDDepths += FDDepth * FDLen 
        if FDDepths / FDSumLen > SECutoff: # non-SE
            if len(CFDType) == 1:
                if CFDType[0] == 'Whole' or CFDType[0] == 'Central':
                    ClusterDic = M9_BlockClustering_Stepwise(FDBlockList, BlockDic, ClusteringSize)
                    for k in ClusterDic:
                        TPosKey = k
                        TLoci, FLoci = Pos2Str(TPosKey[0], TPosKey[1], TPosKey[2]), Pos2Str(Chr, ClusterDic[TPosKey]['FDPos'][0], ClusterDic[TPosKey]['FDPos'][1])
                        ClusterBlockList = ClusterDic[TPosKey]['BlockNo']
                        DisDic = M9_Sub_PairConnectionOVLP4NormalFinal(TLoci, FLoci, InsertSizeCutoff, FlankingSize, BamPath, "Both", '-', JobPath)
                        if DisDic['OverlapSize'] > 0:
                            InContigDic[Chr][ContigLoci]['Pass'] += ClusterBlockList
                            InContigDic[Chr][ContigLoci]['Cluster'][str(ClusterBlockList)] = 'Pass'
                            for BlockNo in ClusterBlockList:    InBlockDic[BlockNo][-1]['Pass'].append(ContigLoci)
                        else: 
                            InContigDic[Chr][ContigLoci]['Filtered'] += ClusterBlockList
                            InContigDic[Chr][ContigLoci]['Cluster'][str(ClusterBlockList)] = 'Filtered'
                            for BlockNo in ClusterBlockList:    InBlockDic[BlockNo][-1]['Filtered'].append(ContigLoci)

                elif CFDType[0] == 'Edge':
                    DiscResult = M9_EdgePosbyCED(FDBlockList, BlockDic, ContigLoci, FlankingSize, InsertSizeCutoff, BamPath, JobPath, ContigDic) # Clustering and falnking finding by CED and resulting
                    if DiscResult[0] == 'nonFD': 
                        InContigDic[Chr][ContigLoci]['Filtered'] += FDBlockList
                        for BlockNo in FDBlockList:     InBlockDic[BlockNo][-1]['Filtered'].append(ContigLoci)
                    else: # return Whole BlockList or Partial BlockList
                        for k in DiscResult[0]:
                            BlockNo = k
                            InContigDic[Chr][ContigLoci]['Pass'] += [BlockNo]
                            InBlockDic[BlockNo][-1]['Pass'].append(ContigLoci)
            
            elif len(CFDType) == 2: # Hybrid type: ex> OVLP found by cactus..
                ClusterDic = M9_BlockClustering_Stepwise(FDBlockList, BlockDic, ClusteringSize)
                for k in ClusterDic:
                    TPosKey = k
                    TLoci, FLoci = Pos2Str(TPosKey[0], TPosKey[1], TPosKey[2]), Pos2Str(Chr, ClusterDic[TPosKey]['FDPos'][0], ClusterDic[TPosKey]['FDPos'][1])
                    ClusterBlockList = ClusterDic[TPosKey]['BlockNo']
                    CED =  CEDCalculation(ContigLoci, FLoci)
                    if CED <= 550: # Edge
                        DiscResult = M9_EdgePosbyCED(ClusterBlockList, BlockDic, ContigLoci, FlankingSize, InsertSizeCutoff, BamPath, JobPath, ContigDic)
                        if DiscResult[0] == 'nonFD': 
                            InContigDic[Chr][ContigLoci]['Filtered'] += ClusterBlockList
                            InContigDic[Chr][ContigLoci]['Cluster'][str(ClusterBlockList)] = 'Filtered'
                            for BlockNo in ClusterBlockList:    InBlockDic[BlockNo][-1]['Filtered'].append(ContigLoci)
                        else:
                            for l in DiscResult[0]:
                                BlockNo = l
                                InContigDic[Chr][ContigLoci]['Pass'] += [BlockNo]
                                InContigDic[Chr][ContigLoci]['Cluster'][str(BlockNo)] = 'Pass'
                                InBlockDic[BlockNo][-1]['Pass'].append(ContigLoci)
                            

                    if CED > 550: # Central
                        DisDic = M9_Sub_PairConnectionOVLP4NormalFinal(TLoci, FLoci, InsertSizeCutoff, FlankingSize, BamPath, "Both", '-', JobPath) #Central
                        if DisDic['OverlapSize'] > 0:
                            InContigDic[Chr][ContigLoci]['Pass'] += ClusterBlockList
                            InContigDic[Chr][ContigLoci]['Cluster'][str(ClusterBlockList)] = 'Pass'
                            for BlockNo in ClusterBlockList:    InBlockDic[BlockNo][-1]['Pass'].append(ContigLoci)
                        else:
                            InContigDic[Chr][ContigLoci]['Filtered'] += ClusterBlockList
                            InContigDic[Chr][ContigLoci]['Cluster'][str(ClusterBlockList)] = 'Filtered'
                            for BlockNo in ClusterBlockList:    InBlockDic[BlockNo][-1]['Filtered'].append(ContigLoci)
        else: 
            InContigDic[Chr][ContigLoci]['Pass'] += FDBlockList
            for BlockNo in FDBlockList: InBlockDic[BlockNo][-1]['Pass'].append(ContigLoci)
        #print(ContigLoci + '\t' + str(InContigDic[Chr][ContigLoci]))
    
    return InBlockDic, InContigDic


def Main10_BlockDic2txt(BlockDic, Name):
    with open(JobPath + '{0}CactusXPurgeDup_Filtered.txt'.format(Name),'w') as fw:
        for i in BlockDic:
            BlockNo, DD, Tool = i, BlockDic[i], i.split(' ')[0]
            ##print(str(len(DD)) + '\t' + str(DD))
            if Tool == 'Cactus': 
                FChr, FStart, FEnd, TStart, TEnd, TChr, RChr, RStart, REnd, Filtered = DD[0], DD[1], DD[2], DD[3], DD[4], DD[5], DD[6], DD[7], DD[8], DD[-1]
                DL = [FChr, FStart, FEnd, TStart, TEnd, TChr, RChr, RStart, REnd, Filtered]
            else: 
                FChr, FStart, FEnd, TStart, TEnd, TChr, PurgeType, Filtered = DD[0], DD[1], DD[2], DD[3], DD[4], DD[5], DD[6], DD[-1]
                DL = [FChr, FStart, FEnd, TStart, TEnd, TChr, PurgeType, '-', '-', Filtered]
            if len(Filtered['Pass']) == 0: continue
            fw.write(BlockNo + '\t' + '\t'.join(list(map(str, DL))) + '\n')


def Main11_FDLoci2Bed(FDLociPath):
    os.system("awk -F \"\t\" '{print $2\"\t\"$3-1\"\t\"$4\"\t\"$7\":\"$5-1\"-\"$6}' " + FDLociPath + " | sort -n -k1,1 -k2,2 > {0}".format(FDLociPath.split("Cactus")[0] + "FD.bed" ))
    os.system("bedtools merge -d 1 -i {0} > {1}".format(FDLociPath.split("Cactus")[0] + "FD.bed", FDLociPath.split("Cactus")[0] + "FD_merged.bed"))


if __name__ == "__main__":
    Parameter = sys.argv
    ConfigPath= Parameter[1]
    FirstImport = Main0_ConfigLoad(ConfigPath)
    PD, Maf = FirstImport[0], FirstImport[1]
    # Parameter setting
    JobPath = PD['JobPath']
    Lenlimit, LenCutoff, Cores, HapChr, InsertSizeCutoff, Ref, Tar, Alt = 0, 0, int(PD['Cores']), PD['HapChr'], int(PD['InsertSizeCutoff']), PD['Ref'], PD['Tar'], '-'
    Asm1DepthCutoff, Asm1HapDepthCutoff, Asm2DepthCutoff, Asm2HapDepthCutoff, Asm1SECutoff, Asm2SECutoff = float(PD['Asm1DepthCutoff']), float(PD['Asm1HapDepthCutoff']), float(PD['Asm2DepthCutoff']), float(PD['Asm2HapDepthCutoff']), float(PD['Asm1SECutoff']), float(PD['Asm2SECutoff'])
    Asm1FaiPath, Asm2FaiPath, InputMaf, Asm1BamPath, Asm2BamPath, Asm1GapPath, Asm2GapPath =  PD['Asm1FaiPath'], PD['Asm2FaiPath'], PD['InputMaf'], PD['Asm1BamPath'], PD['Asm2BamPath'], PD['Asm1GapPath'], PD['Asm2GapPath']
    Asm1Path, Asm2Path, Asm1AsmKmerDirPath, Asm2AsmKmerDirPath, Asm1ReadKmerDirPath, Asm2ReadKmerDirPath =  PD['Asm1Path'], PD['Asm2Path'], PD['Asm1AsmKmerDirPath'], PD['Asm2AsmKmerDirPath'], PD['Asm1ReadKmerDirPath'], PD['Asm2ReadKmerDirPath']
    Asm1PurgeDupPath, Asm2PurgeDupPath, Asm1PurgePafPath, Asm2PurgePafPath =  PD['Asm1PurgeDupPath'], PD['Asm2PurgeDupPath'], PD['Asm1PurgePafPath'], PD['Asm2PurgePafPath']
    ChrNameFile = PD['ChrPairPath']
    PreFaiPath, VGPFaiPath, PreBamPath, VGPBamPath, PreGapPath, VGPGapPath = Asm1FaiPath, Asm2FaiPath, Asm1BamPath, Asm2BamPath, Asm1GapPath, Asm2GapPath
    ClusteringSize, FlankingSize = InsertSizeCutoff, InsertSizeCutoff
    PreCutoff, PreHapCutoff, VGPCutoff, VGPHapCutoff, PreSECutoff, VGPSECutoff = Asm1DepthCutoff, Asm1HapDepthCutoff, Asm2DepthCutoff, Asm2HapDepthCutoff, Asm1SECutoff, Asm2SECutoff
    PrePurgeDupPath, VGPPurgeDupPath, PrePurgePafPath, VGPPurgePafPath = Asm1PurgeDupPath, Asm2PurgeDupPath, Asm1PurgePafPath, Asm2PurgePafPath

    #-preprocessing-
    ChrPairDic = ChrConversion(ChrNameFile)
    Fai2Bed(Asm1FaiPath, Asm2FaiPath, JobPath)
    GCA2AG(Asm1GapPath, ChrPairDic, JobPath, '1')
    GCA2AG(Asm2GapPath, ChrPairDic, JobPath, '2')
    AllDepth4DepthGap(PreBamPath, JobPath, '1')
    AllDepth4DepthGap(VGPBamPath, JobPath, '2')
    Merqury2Loci(Asm1Path, Asm1AsmKmerDirPath, Asm1ReadKmerDirPath, JobPath, '1')
    Merqury2Loci(Asm2Path, Asm2AsmKmerDirPath, Asm2ReadKmerDirPath, JobPath, '2')
    KmerDicBulding(JobPath, '1')
    KmerDicBulding(JobPath, '2')
    AGDGmerging(JobPath, ChrPairDic) # to bed file
    FaiDic = Fai2Dic(Asm1FaiPath, Asm2FaiPath)
    ContigPreDic = ContigbyAGDG('1', JobPath)
    ContigVGPDic = ContigbyAGDG('2', JobPath)
    #BinSaving(ChrPairDic, 'ChrPairDic')
    #BinSaving(FaiDic, 'FaiDic')
    print("PreProcessing Done")
    print("MainStep 1")
    #-MainStep-
     
    Main1 = Main1_CactusCalling(Maf) # return MafDic
    MainDic, DirectionDic = Main1[0], Main1[1]
    Maf, Main1 = 0, 0
    #BinSaving(JobPath, MainDic, 'MainDic')
    #BinSaving(JobPath, MainDic, 'DirectionDic')
    
    print("MainStep 2")
    Main2 = Main2_ContigCalling(MainDic)
    MainDic, DepthInput, ContigPreListDic, ContigVGPListDic = Main2[0], Main2[1], Main2[2], Main2[3]
    #BinSaving(JobPath, MainDic, 'MainDic2')
    #BinSaving(JobPath, ContigPreListDic, 'Contig1ListDic')
    #BinSaving(JobPath, ContigVGPListDic, 'Contig2ListDic')
    
    #MainDic = BinLoading(JobPath, 'MainDic2')
    #ContigPreListDic = BinLoading(JobPath, 'Contig1ListDic')
    #ContigVGPListDic = BinLoading(JobPath, 'Contig2ListDic')
    Loaded_DepthDic = {}
    print("MainStep 3")
    Main3 = Main3_DetphAdd2FDDic(MainDic, Loaded_DepthDic, Cores, PreBamPath, VGPBamPath)
    MainDic, New_DepthDic = Main3[0], Main3[1]
    #BinSaving(JobPath, MainDic, 'MainDic_DepthLoad')
    #BinSaving(JobPath, New_DepthDic, 'New_DepthDic')
    print("MainStep 4")
    MainDic = Main4_CactusFalseTrueAllocating(MainDic)
    #BinSaving(JobPath, MainDic, 'MainDic_Bool')
    print("MainStep 5,6")
    #MainDic = BinLoading(JobPath, 'MainDic_Bool')
    MainFDList = Main56_DepthAndGapFiltering(MainDic)
    MainPreFDList, MainVGPFDList = MainFDList[0], MainFDList[1]
    #BinSaving(JobPath, MainPreFDList, 'MainAsm1FD')
    #BinSaving(JobPath, MainVGPFDList, 'MainAsm2FD')
    print("MainStep 7")
    #MainPreFDList, MainVGPFDList, ContigPreListDic, ContigVGPListDic = BinLoading(JobPath, 'MainAsm1FD'), BinLoading(JobPath, 'MainAsm2FD'), BinLoading(JobPath, 'Contig1ListDic'), BinLoading(JobPath, 'Contig2ListDic')
    Main7Pre = Main7_CactusXPurgeDup(MainPreFDList, PrePurgeDupPath, ContigPreDic, ContigPreListDic)
    Main7VGP = Main7_CactusXPurgeDup(MainVGPFDList, VGPPurgeDupPath, ContigVGPDic, ContigVGPListDic) 
        
    MainPreContigDic, MainVGPContigDic = Main7Pre[0], Main7VGP[0]
    MainPreBlockDic, MainVGPBlockDic = Main7Pre[1], Main7VGP[1]
    PurgeBlockListPre, PurgeBlockListVGP = Main7Pre[2], Main7VGP[2]
    DGPreBlockDic, DGVGPBlockDic = Main7Pre[3], Main7VGP[3]
    #BinSaving(JobPath, MainPreContigDic, 'MainAsm1CD')
    #BinSaving(JobPath, MainVGPContigDic, 'MainAsm2CD')
    #BinSaving(JobPath, MainPreBlockDic, 'MainAsm1BD')
    #BinSaving(JobPath, MainVGPBlockDic, 'MainAsm2BD')
    #BinSaving(JobPath, PurgeBlockListPre, 'PurgeAsm1')
    #BinSaving(JobPath, PurgeBlockListVGP, 'PurgeAsm2')
    #BinSaving(JobPath, DGPreBlockDic, 'DGAsm1Block')
    #BinSaving(JobPath, DGVGPBlockDic, 'DGAsm2Block')
    
    #PurgeBlockListPre = BinLoading(JobPath, 'PurgeAsm1')
    #PurgeBlockListVGP = BinLoading(JobPath, 'PurgeAsm2')
    #MainPreContigDic = BinLoading(JobPath, 'MainAsm1CD')
    #MainVGPContigDic = BinLoading(JobPath, 'MainAsm2CD')
    #MainPreBlockDic = BinLoading(JobPath, 'MainAsm1BD')
    #MainVGPBlockDic = BinLoading(JobPath, 'MainAsm2BD')
    
    PurgePafPreDic = PurgeDupPafLoad(PrePurgePafPath)
    PurgePafVGPDic = PurgeDupPafLoad(VGPPurgePafPath)

    PAPreDic, PAVGPDic = PurgeDupLociSpecify(PurgeBlockListPre, PurgePafPreDic), PurgeDupLociSpecify(PurgeBlockListVGP, PurgePafVGPDic)
    #BinSaving(JobPath, PAPreDic, 'PAAsm1Dic')
    #BinSaving(JobPath, PAVGPDic, 'PAAsm2Dic')
    print("MainStep 8")
    MainPreContigDic = Main8_ContigFDTypeAllocation(MainPreContigDic, InsertSizeCutoff)
    MainVGPContigDic = Main8_ContigFDTypeAllocation(MainVGPContigDic, InsertSizeCutoff)
    #BinSaving(JobPath, MainPreContigDic, 'MainAsm1CD+FDType')
    #BinSaving(JobPath, MainVGPContigDic, 'MainAsm2CD+FDType')
      
    #MainPreContigDic = BinLoading(JobPath, 'MainAsm1CD+FDType')
    #MainVGPContigDic = BinLoading(JobPath, 'MainAsm2CD+FDType')
    #MainPreBlockDic = BinLoading(JobPath, 'MainAsm1BD')
    #MainVGPBlockDic = BinLoading(JobPath, 'MainAsm2BD')
    #PAPreDic = BinLoading(JobPath, 'PAAsm1Dic')
    #PAVGPDic = BinLoading(JobPath, 'PAAsm2Dic')
    MainPreBlockDic = Main8_2_PurgeTLociTrans(MainPreBlockDic, PAPreDic)
    MainVGPBlockDic = Main8_2_PurgeTLociTrans(MainVGPBlockDic, PAVGPDic)
    #BinSaving(JobPath, MainPreBlockDic, 'Main8_2Asm1BD')
    #BinSaving(JobPath, MainVGPBlockDic, 'Main8_2Asm2BD')
     
    #MainPreContigDic = BinLoading(JobPath, 'MainAsm1CD+FDType')
    #MainVGPContigDic = BinLoading(JobPath, 'MainAsm2CD+FDType')
    #MainPreBlockDic = BinLoading(JobPath, 'Main8_2Asm1BD')
    #MainVGPBlockDic = BinLoading(JobPath, 'Main8_2Asm2BD')

    f0, Maf, ContigPreDic, ContigVGPDic, PurgePafPreDic, PurgePafVGPDic = 0, 0, 0, 0, 0, 0
    Main1, MainDic, DirectionDic, Main2, DepthInput, ContigPreListDic, ContigVGPListDic = 0, 0, 0, 0, 0, 0, 0
    Main3, MainFDList, MainPreFDList, MainPreVGPFDList, Main7Pre, Main7VGP, PAPreDic, PAVGPDic = 0, 0, 0, 0, 0, 0, 0, 0
    print("MainStep 9")
    #New_DepthDic = BinLoading(JobPath, 'New_DepthDic')
    Main9Pre = M9_MultiProc(MainPreContigDic, MainPreBlockDic, InsertSizeCutoff, FlankingSize, ClusteringSize, PreBamPath, PreSECutoff, Cores, JobPath, New_DepthDic)
    MainPreContigDic, MainPreBlockDic = Main9Pre[1], Main9Pre[0]
    #BinSaving(JobPath, MainPreContigDic, 'Main9Asm1ContigDic')
    #BinSaving(JobPath, MainPreBlockDic, 'Main9Asm1BlockDic') 
     
    Main9VGP = M9_MultiProc(MainVGPContigDic, MainVGPBlockDic, InsertSizeCutoff, FlankingSize, ClusteringSize, VGPBamPath, VGPSECutoff, Cores, JobPath, New_DepthDic)
    MainVGPContigDic, MainVGPBlockDic = Main9VGP[1], Main9VGP[0]
    #BinSaving(JobPath, MainVGPContigDic, 'Main9Asm2ContigDic')
    #BinSaving(JobPath, MainVGPBlockDic, 'Main9Asm2BlockDic')

    #print("InputLoci" + '\t' + str(Main9Pre[2]))
    #print("MergedLoci" + '\t' + str(Main9Pre[3]))
    #print("InputLoci" + '\t' + str(Main9VGP[2]))
    #print("MergedLoci" + '\t' + str(Main9VGP[3]))
    
    print("Saving to Bed")
    #DGPreBlockDic = BinLoading(JobPath, 'DGAsm1Block')
    #DGVGPBlockDic = BinLoading(JobPath, 'DGAsm2Block')
    #PreBD = {**BinLoading(JobPath, 'Main9Asm1BlockDic'), **DGPreBlockDic} 
    #VGPBD = {**BinLoading(JobPath, 'Main9Asm2BlockDic'), **DGVGPBlockDic}
    PreBD = {**MainPreBlockDic, **DGPreBlockDic} 
    VGPBD = {**MainVGPBlockDic, **DGVGPBlockDic}

    Main10_BlockDic2txt(PreBD, 'Asm1') #Prefix name
    Main10_BlockDic2txt(VGPBD, 'Asm2')
    Main11_FDLoci2Bed(JobPath + 'Asm1CactusXPurgeDup_Filtered.txt')
    Main11_FDLoci2Bed(JobPath + 'Asm2CactusXPurgeDup_Filtered.txt')


