# 1. Block Type Defining
# 2. Duplication Calling : One to many, Many to Many + 1
# 3. Dictionary Save : {BlockNo: {T1:{Chr:,Start,End},T2:..
import numpy as np
import time
import pickle
import os
import parmap
import subprocess as sp
import random
from cigar import Cigar
from multiprocessing import Manager


if __name__ == "__main__": 
    #Fixed Parameter
    Lenlimit = 0#20 # MAF Length cutoff (Reference)
    LenCutoff = 0#0.8 #0.9 # Old Length % cut off 
    Cores = 40
    HapChr = ['NC_009114','NC_009115','NC_009116','NC_009118'] + ['CM014221', 'CM014222', 'CM014223', 'CM014224', 'CM014225', 'RZJT01000027', 'RZJT01000028', 'RZJT01000029', 'RZJT01000030', 'RZJT01000031'] # only VGPOrnana + PreOrnanaSex( VGP has half of sex chromsome reads for pre. assembly)
    MQCutoff = 20
    MultiMQCutoff = 3
    InsertSizeCutoff = 550
    ClusteringSize = 550
    FlankingSize = 550
    MQSampling = 1000
    ###################

    #Operation parameter : Zebra finch     
    Path, KPath, PreFaiPath, VGPFaiPath = '/home/ko/cactus/02/', '/disk3/ko_d3/', '/home/ko/rawdata/pre_taegut.fasta.fai', '/home/ko/rawdata/vgp_taegut.fasta.fai'
    Input = 'Zebrafinch2_2'
    Ref, Tar, Alt = 'Zebrafinch1', 'Zebrafinch0', 'Zebrafinch2'
    PreBamPath, VGPBamPath = '/disk1/ko_d1/merged_L001_to_L008.sorted.bam', '/disk1/ko_d1/merged_L001_to_L008.toVGP.sorted.bam'
    PreAG, VGPAG, PreDG, VGPDG, PreKmerDic, VGPKmerDic = '02AssemblyGap_parsing.txt', '02VGPAssemblyGap_parsing.txt', '02PreDepthGap_3TNew.txt', '02VGPDepthGap_3TNew.txt', '02PreKmerDic/', '02VGPKmerDic/'
    PreCutoff = 45 # between minimum 4 bimodal distribution
    PreHapCutoff = '-'
    VGPCutoff = 67 * 0.75 #
    PreSECutoff, VGPSECutoff = 5, 2
    ChrNameFile = 'ChrName'
    PrePurgeDupPath, VGPPurgeDupPath = 'Predups.bed', 'VGPdups.bed'
    PrePurgePafPath, VGPPurgePafPath = 'pre_taegut.fasta.split.self.paf', 'vgp_taegut.fasta.split.self.paf'
    
    #Operation parameter : Hummingbird
    '''
    Path, KPath, PreFaiPath, VGPFaiPath = '/home/ko/cactus/04/', '/disk3/ko_d3/', '/home/ko/rawdata/Calann1_1_SH.fna.fai', '/home/ko/rawdata/Calann0_1_SH.fna.fai'
    Input = 'Calann'
    Ref, Tar, Alt = 'Calann0', 'Calann2', 'Calann1'
    PreBamPath, VGPBamPath = '/disk3/ko_d3/0410x/04_10x_pre/Calann10x2Pre.bam', '/disk3/ko_d3/0410x/04_10x_VGP/Calann10x2VGP.bam'
    PreAG, VGPAG, PreDG, VGPDG, PreKmerDic, VGPKmerDic = '04AssemblyGap_parsing.txt', '04VGPAssemblyGap_parsing.txt', '04PreDepthGap_3TNew.txt', '04VGPDepthGap_3TNew.txt', '04PreKmerDic/', '04VGPKmerDic/'
    PreCutoff = 37 * 0.75 # between minimum 4 bimodal distribution
    PreHapCutoff = '-'
    VGPCutoff = 37 * 0.75 #
    PreSECutoff, VGPSECutoff = 8, 2
    ChrNameFile = 'ChrName_OldOld'
    #ChrNameFile = 'ChrName_OldNew' #onlyused for FDResulting.py ChrPairDic producing
    PrePurgeDupPath, VGPPurgeDupPath = 'Predups.bed', 'VGPdups.bed'
    PrePurgePafPath, VGPPurgePafPath = 'Calann1_1_SH.fna.split.self.paf', 'Calann0_1_SH.fna.split.self.paf'
    '''

    #Operation parameter : platypus
    '''
    Path, KPath, PreFaiPath, VGPFaiPath = "/home/ko/cactus/09_2/", '/disk3/ko_d3/', '/home/ko/rawdata/Ornana1_1DMHD.fna.fai', '/home/ko/rawdata/Ornana0_1HD.fna.fai'
    Input = "Ornana_2"
    Ref, Tar, Alt = 'Ornana0_1', 'Ornana1_1DM', 'Ornana2'
    PreBamPath, VGPBamPath = '/disk3/ko_d3/0910x/09_10x_pre/Ornana10X2Old.bam', '/disk3/ko_d3/0910x/09_10x_VGP/Ornana10X2VGP.bam'
    PreAG, VGPAG, PreDG, VGPDG, PreKmerDic, VGPKmerDic = '09AssemblyGap_parsing.txt', '09VGPAssemblyGap_parsing.txt', '09PreDepthGap_3TNew.txt', '09VGPDepthGap_3TNew.txt', '09PreKmerDic/', '09VGPKmerDic/'
    PreCutoff =  65 # lowest delta depth coverage #68 # between minimum 4 bimodal distribution
    PreHapCutoff = 65 * 0.75
    VGPHapCutoff = 62 * 0.75 # haploid maximum
    VGPCutoff = 68 # lowest depth between bimodal distribution #119 * 0.75 #
    PreSECutoff, VGPSECutoff = 22, 9
    ChrNameFile = 'ChrName_OldOld' #used for Main
    #ChrNameFile = 'ChrName_OldNew' #only used for FDResulting.py ChrPairDic Producing 
    PrePurgeDupPath, VGPPurgeDupPath = 'Predups.bed', 'VGPdups.bed'
    PrePurgePafPath, VGPPurgePafPath = 'Ornana1_1DMHD.fna.split.self.paf', 'Ornana0_1HDDS.fna.split.self.paf'
    '''

    ############################
    DepthPrePath = Path + 'PreFDLocibyDic.txt'
    DepthVGPPath = Path + 'VGPFDLocibyDic.txt'
    PrePurgeDupPath, VGPPurgeDupPath = 'Predups_test.bed', 'VGPdups_test.bed'
    LogPath = Path + 'Log/'                            ###
    
    f0 = open(Path+Input+".maf", 'r')                  ###
    f3 = open(Path + ChrNameFile + '.txt', 'r')
    
    print(f0.readline())
    print(f0.readline())
    print(f0.readline())
    print(f0.readline())
     
    Maf = f0.read().split('\na\n')
    f0.close()
    if Maf.count('') > 0: Maf.remove('')


### [Normal Tool]
def BinSaving(PythonData, FileName):#c     
    with open(LogPath + FileName+'.Pybin', 'wb') as fb: pickle.dump(PythonData, fb)

def BinLoading(FileName):#c      
    with open(LogPath + FileName+'.Pybin', 'rb') as rb: return pickle.load(rb)

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


### [Pre-Processing]
def Fai2Bed():#c # Preprocessing; contig formation; both assemblies
    os.system("awk '{{print $1\"\t\"0\"\t\"$2}}' {0} > {1}PreFaibed.bed".format(PreFaiPath, LogPath)) # Fai2bed
    os.system("awk '{{print $1\"\t\"0\"\t\"$2}}' {0} > {1}VGPFaibed.bed".format(VGPFaiPath, LogPath))


def Fai2Dic():
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


def ChrConversion():#c # Preprocessing; sca. name conversion; both assemblies
    ChrPair = f3.read().split('\n')
    if ChrPair.count('') >0: ChrPair.remove('')
    ChrPairDic = {}
    for i in ChrPair:
        DT = i.split('\t')
        CM, NC = DT[0], DT[1]
        ChrPairDic[CM] = NC
    return ChrPairDic


def AGDGmerging(PreAG, PreDG, VGPAG, VGPDG, PreKmerDic, VGPKmerDic):#c # Preprocessing; Gap load; Both assemblies 
    FinalGap = {} # {'Chr:Start-End':Type}
    f_PreAG, f_PreDG, = open(Path + '../../Bam/' + PreAG ), open(Path + '../../Bam/' + PreDG ) # KDic = {'Chr:Start' : Depth}
    f_VGPAG, f_VGPDG, = open(Path + '../../Bam/' + VGPAG ), open(Path + '../../Bam/' + VGPDG )
    PreAG, PreDG, VGPAG, VGPDG = f_PreAG.read().split('\n'), f_PreDG.read().split('\n'), f_VGPAG.read().split('\n'), f_VGPDG.read().split('\n')
    if PreAG.count('') > 0: PreAG.remove('')
    if PreDG.count('') > 0: PreDG.remove('')
    if VGPAG.count('') > 0: VGPAG.remove('')
    if VGPDG.count('') > 0: VGPDG.remove('')
    file_list_Pre, file_list_VGP = os.listdir(KPath + PreKmerDic), os.listdir(KPath + VGPKmerDic)
    KGapDic = {}
    def ZeroKmerParsing(file_list, DicPath):
        Step = 0
        for i in file_list:
            Step += 1
            print('zeroK parsing: ' + str(len(file_list) - Step))
            fk = open(KPath + DicPath + i,'rb')
            Dic = pickle.load(fk)
            for j in Dic:
                Chr, DD = j, Dic[j]
                for k in DD:
                    Pos, Copy, Multi = k, DD[k][0], DD[k][1]
                    if Multi == 0:
                        if Chr in KGapDic: KGapDic[Chr].append([str(int(Pos) + 1), str(int(Pos)+20)])
                        else: KGapDic[Chr] = [[str(int(Pos) + 1), str(int(Pos)+20)]] # K-mer = bed format (Pos+1 - Pos+k)
    
    ZeroKmerParsing(file_list_Pre, PreKmerDic) 
    ZeroKmerParsing(file_list_VGP, VGPKmerDic) #Zero K-mer save
    ##BinSaving(KGapDic, '04KGapDic') 
    ##KGapDic = BinLoading('04KGapDic')
    AGDic = {}
    def AGSaving(GapList):
        for i in GapList:
            DT = i.split('\t')
            Chr, Start, End = DT[0], DT[1], DT[2]
            Loci = Pos2Str(Chr, Start, End) #none-bed start
            FinalGap[Loci] = 'AG'
            if Chr in AGDic: AGDic[Chr].append([Start, End])
            else: AGDic[Chr] = [[Start, End]]
    
    AGSaving(PreAG)
    AGSaving(VGPAG)
    print("AGSaving")
    def DGSaving(GapList):
        Step = 0
        for i in GapList:
            print(len(GapList) - Step)
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
    print("PreDGSaving")
    DGSaving(VGPDG)
    print("VGPDGSaving")
    with open(LogPath + 'BothGapbed.bed', 'w') as fG:
        for i in FinalGap: 
            Chr, Start, End, Type = i.split(':')[0], i.split(':')[1].split('-')[0], i.split('-')[1], FinalGap[i]
            fG.write('\t'.join([Chr, str(int(Start)-1), End, Type]) + '\n') # Gap2bed


def ContigbyAGDG(AssemblyType):#c # Preprocessing; contig segmenting; factors for assembly.
    ContigBed = sp.getoutput("bedtools subtract -a {0} -b {1}".format(LogPath + AssemblyType + 'Faibed.bed', LogPath + 'BothGapbed.bed')).split('\n')
    if ContigBed.count('') >0: ContigBed.remove('') 
    ContigDic = {}
    for i in ContigBed: # bed2back
        DT = i.split('\t')
        Chr, Start, End = DT[0], str(int(DT[1]) + 1), DT[2]
        Loci = Chr + ':' + Start + '-' + End
        if Chr in ContigDic: ContigDic[Chr][Loci] = {'FDPos':[],'FDBlock':[]}
        else: ContigDic[Chr] = {Loci:{'FDPos':[],'FDBlock':[]}} # {Chr:{Loci1:{'FDPos':[Start, End], 'FDBlock':[BLockNo1, BLockNo2]},Loci2:{}}
    
    return ContigDic


def PurgeDupPafLoad(Path, PafPath):#c
    with open(Path + PafPath, 'r') as fr:
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
        print(PurgeNo)
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
        if (len(ProcessedMaf) - Step) % 100 == 0: print('Main1\t' + str(len(ProcessedMaf) - Step))
        Breaker = 0
        DL = i.split('\n')
        if DL.count('') > 0: DL.remove('')

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
            Chr, Start, End, Name = ChrPairDic[k], RefPos[Idx][0], RefPos[Idx][1], 'R' + str(Idx + 1)
            MafDic[MafNo]['Dat'][Name] = {'Chr':Chr,'Start':str(Start),'End':str(End)}
            Key = '_'.join([Chr, str(Start), str(End)])
            Idx += 1
          
        Idx = 0
        for k in TarList:
            Chr, Start, End, Name = ChrPairDic[k], TarPos[Idx][0], TarPos[Idx][1], 'T' + str(Idx + 1)
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
        if (LenLen - Step) % 100 == 0: print('Main2\t' + str(LenLen - Step))
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
        if (len(DataList) - Step) % 100 == 0: print(len(DataList) - Step)
        LL = Str2Pos(Loci)
        Chr, Start, End = LL[0], LL[1], LL[2]
        Len = int(End) - int(Start) + 1
        TotalDepth = sp.getoutput("samtools depth -r {0} {1} | awk 'BEGIN{{Depth = 0}} {{Depth += $3}} END{{print Depth}}'".format(Loci, BamPath))
        MeanDepth = int(TotalDepth) / Len
        ResultDic[Loci] = str(round(float(MeanDepth),1))
    return ResultDic

#### END Main3-Subtool ####

def Main3_DetphAdd2FDDic(FDDic, Loaded_DepthDic, Cores):
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
        if (len(InFDDic) - Step) % 100 == 0: print('Main3\t' + str(len(InFDDic) - Step))
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

#input_data = np.array_split(DepthLack, Cores)
#parmap.map(M_DepthMean, input_data, pm_pbar = True)

def Main4_CactusFalseTrueAllocating(FDDic): # Analyzing; False True Decision; Both assembly
    InFDDic = FDDic
    Step = 0
    for i in InFDDic:
        Step += 1
        if (len(InFDDic) - Step) % 100 == 0: print('Main4\t' + str(len(InFDDic) - Step))
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
        if (len(InFDDic) - Step) % 100 == 0: print('Main56\t' + str(len(InFDDic) - Step))
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
    f_Purge = open(Path + PurgeDupPath, 'r')
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
        if (len(CactusBlockList + PurgeBlockList) - Step) % 100 == 0: print('Main7\t' + str(len(CactusBlockList + PurgeBlockList) - Step))
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


def Main8_ContigFDTypeAllocation(ContigDic):  # {Chr:{Loci1:{'FDPos':[Start, End], 'FDBlock':[BLockNo1, BLockNo2], 'CFDType':['Whole' or 'Edge','Central'..]'},Loci2:{}}
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
                    if CED <= 550: InContigDic[Chr][Loci]['CFDType'].append('Edge')
                    if CED > 550: InContigDic[Chr][Loci]['CFDType'].append('Central')
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


### [M9-Sub tool] Pair Connection Counting
def M9_Sub_ReadLogger(FLoci, TLoci, ParsingRegion, BlockNo, RN, Status, LogPath):
    with open(LogPath + 'ReadLogger.txt', 'a') as fwl:
        NowTime = time.ctime(time.time())
        fwl.write('\t'.join([NowTime, FLoci, TLoci, ParsingRegion, BlockNo, RN, Status]) + '\n')

def M9_Tool_SuppleReadFilter(ReadList):
    PassReads = []
    if ReadList.count('') > 0: ReadList.remove('')
    for i in ReadList:
        DT = i.split('\t')
        if len(DT) > 1: pass#20210202
        else:               #
            print(DT)       #
            print(ReadList) #
            continue        #
        RN, FS = DT[0], DT[1]
        if FSDecimal(FS)[256] == 1: continue
        if FSDecimal(FS)[1024] == 1: continue
        PassReads.append(i)
    return PassReads

def M9_ColChr(PairChr, Chr):
    if PairChr == '=': return Chr
    else: return PairChr

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


def M9_EdgePosbyCED(FDBlockList, BlockDic, ContigLoci, FlankingSize, InsertSizeCutoff, BamPath, LogPath, ContigDic): # Clustering embedded
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
        FDiscDic = M9_Sub_PairConnectionOVLP4NormalFinal(CEDFLociStr, CEDTLociStr, InsertSizeCutoff, FlankingSize, BamPath, 'Edge', FFlanking, LogPath)
        TDiscDic = M9_Sub_PairConnectionOVLP4NormalFinal(CEDTLociStr, CEDFLociStr, InsertSizeCutoff, FlankingSize, BamPath, 'Edge', TFlanking, LogPath)
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
            FDiscDic = M9_Sub_PairConnectionOVLP4NormalFinal(FLoci, TLoci, InsertSizeCutoff, FlankingSize, BamPath, 'Edge', FFlanking, LogPath)
            TDiscDic = M9_Sub_PairConnectionOVLP4NormalFinal(TLoci, FLoci, InsertSizeCutoff, FlankingSize, BamPath, 'Edge', TFlanking, LogPath)
            if FDiscDic['OverlapSize'] + TDiscDic['OverlapSize'] > 0: 
                DiscBlockList.append(BlockNo)
                DiscBlockDic[BlockNo] = [FDiscDic, TDiscDic]
            Index += 1
        if len(DiscBlockList) == 0: return 'nonFD', ['Seperately'], []
        else: return DiscBlockList, ['Seperately'], DiscBlockDic
                            

def M9_Sub_PairConnectionOVLP4NormalFinal(TLoci, FLoci, InsertSizeCutoff, FlankingSize, BamPath, FlankingStat, O1ParsingLoci, LogPath): #c
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


def M9_Sub_PairConnectionOVLP4Normal_BackTrace4WCD(ParsingLoci, TLociContig, InsertSizeCutoff, FlankingSize, BamPath): # 1. DiscReadParsing on Edge of False Contig; 2. Pair Region Merging; 3. Normal Read Parsing in Merged region == FD evidence
    ParsingPos, TPos, DiscDic, DiscReadRegion, TLociDic = Str2Pos(ParsingLoci), Str2Pos(TLociContig), {'DisRead':[], 'NormRead':[]}, [], {}
    TChr, TStart, TEnd = TPos[0], TPos[1], TPos[2]
    ParsingLoci1, ParsingLoci2 = Pos2Str(ParsingPos[0], ParsingPos[1], int(ParsingPos[1] + FlankingSize - 1)), Pos2Str(ParsingPos[0], int(ParsingPos[2]) - 550 + 1, ParsingPos[2])
    ParsingReads, TLociReads = ReadParser(ParsingLoci1, BamPath) + ReadParser(ParsingLoci2, BamPath), ReadParser(TLociContig, BamPath)
    ParsingReads, TLociReads = M9_Tool_SuppleReadFilter(ParsingReads), M9_Tool_SuppleReadFilter(TLociReads)
    for i in TLociReads: #DiscPairedRead + NormalRead
        DT = i.split('\t')
        ReadName, FlagStat, ReadChr, ReadStart, MappingQuality, CigarString, PairReadChr, PairReadPos, InsertSize = DT[0], DT[1], DT[2], int(DT[3]), DT[4], DT[5], M9_ColChr(DT[6], DT[2]), int(DT[7]), abs(int(DT[8]))
        TLociDic[ReadName] = i

    for i in ParsingReads:
        DT = i.split('\t')
        ReadName, FlagStat, ReadChr, ReadStart, MappingQuality, CigarString, PairReadChr, PairReadPos, InsertSize = DT[0], DT[1], DT[2], int(DT[3]), DT[4], DT[5], M9_ColChr(DT[6], DT[2]), int(DT[7]), abs(int(DT[8]))
        if InsertSize > InsertSizeCutoff and PairReadChr == TChr and TStart <= PairReadPos <= TEnd:
            PairInfo = TLociDic[ReadName]
            PairCigar, PairStart = PairInfo[5], PairInfo[3]
            if PairReadPos != PairStart: continue
            PairEnd = PairStart + Cigar4Len(PairCigar) - 1
            DiscReadRegion.append([PairStart, PairEnd])
            DiscDic['DisRead'].append(ReadName)

    DisMergeRegion = M7_OverlapMerge(DiscReadRegion)
    DisLociStrList = []
    for i in DisMergeRegion:
        DisChr, DisStart, DisEnd = TChr, i[0], i[1]
        DisLociStrList.append(Pos2Str(DisChr, DisStart, DisEnd))
    DisMergeRegionReads = ReadParser(' '.join(DisLociStrList), BamPath) # Merging parsing
    DisMergeRegionReads = M9_Tool_SuppleReadFilter(DisMergeRegionReads)
    for i in DisMergeRegionReads:
        ReadName, FlagStat, ReadChr, ReadStart, MappingQuality, CigarString, PairReadChr, PairReadPos, InsertSize = DT[0], DT[1], DT[2], int(DT[3]), DT[4], DT[5], M9_ColChr(DT[6], DT[2]), int(DT[7]), abs(int(DT[8]))
        if InsertSize <= InsertSizeCutoff and PairReadChr == TChr:    DiscDic['NormRead'].append(ReadName)
    return DiscDic


def M9_Sub_PairConnectionOVLP4Normal_BackTrace(ParsingLoci, TLociContig, InsertSizeCutoff, FlankingMulti, BamPath): # 1. Read Parsing; 2. FD Paired Disc. Read region Merging; 3. Merge region read parsing; 4. non-disc. read pair of 3 in TLoci? 
    ParsingPos, TPos, DiscDic, DiscReadRegion, TLociDic = Str2Pos(ParsingLoci), Str2Pos(TLociContig), {'DisRead':[], 'NormRead':[]}, [], {}
    TChr, TStart, TEnd = TPos[0], TPos[1], TPos[2]
    ParsingReads, TLociReads = ReadParser(ParsingLoci, BamPath), ReadParser(TLociContig, BamPath)
    ParsingReads, TLociReads = M9_Tool_SuppleReadFilter(ParsingReads), M9_Tool_SuppleReadFilter(TLociReads)
    for i in TLociReads:
        DT = i.split('\t')
        ReadName, FlagStat, ReadChr, ReadStart, MappingQuality, CigarString, PairReadChr, PairReadPos, InsertSize = DT[0], DT[1], DT[2], int(DT[3]), DT[4], DT[5], M9_ColChr(DT[6], DT[2]), int(DT[7]), abs(int(DT[8]))
        if InsertSize <= InsertSizeCutoff and PairReadChr == TChr:    TLociDic[ReadName] = i

    for i in ParsingReads:
        DT = i.split('\t')
        ReadName, FlagStat, ReadChr, ReadStart, MappingQuality, CigarString, PairReadChr, PairReadPos, InsertSize = DT[0], DT[1], DT[2], int(DT[3]), DT[4], DT[5], M9_ColChr(DT[6], DT[2]), int(DT[7]), abs(int(DT[8]))
        if InsertSize > InsertSizeCutoff and PairReadChr == TChr: 
            PairInfo = PMQPFSParser(ReadName, BamPath, Str2Pos(ParReadChr, PairReadPos, PairReadPos))
            PairCigar = PairInfo[2]
            PairStart, PairEnd = PairReadPos, PairReadPos + Cigar4Len(PairCigar) - 1
            if PairStart <= TEnd and TStart <= PairEnd: continue
            DiscReadRegion.append([PairStart, PairEnd])

    DisMergeRegion = M7_OverlapMerge(DiscReadRegion) # Pair Read Flanking of TLoci Merging
    DisLociStrList = []
    for i in DisMergeRegion: 
        DisChr, DisStart, DisEnd = TChr, i[0], i[1]
        DisLociStrList.append(Pos2Str(DisChr, DisStart, DisEnd))
    DisMergeRegionReads = ReadParser(' '.join(DisLociStrList), BamPath) # Merging parsing
    DisMergeRegionReads = M9_Tool_SuppleReadFilter(DisMergeRegionReads)
    for i in DisMergeRegionReads:
        ReadName, FlagStat, ReadChr, ReadStart, MappingQuality, CigarString, PairReadChr, PairReadPos, InsertSize = DT[0], DT[1], DT[2], int(DT[3]), DT[4], DT[5], M9_ColChr(DT[6], DT[2]), int(DT[7]), abs(int(DT[8]))
        DiscDic['DisRead'].append(ReadName)
        if ReadName in TLociDic:
            TRDT = TLociDic[ReadName].split('\t')
            TReadStart = int(TRDT[3])
            if PairReadPos == TRDTStart: DiscDic['NormRead'].append(ReadName)
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


def M9_Sub_PairConnectionOVLP4Normal(ParsingLoci, ReadParsingList, CheckLoci, NormalLoci):
    DL, DL2 = Str2Pos(CheckLoci), Str2Pos(NormalLoci)
    CheckChr, CheckStart, CheckEnd = DL[0], int(DL[1]), int(DL[2])
    NormalChr, NormalStart, NormalEnd = DT2[0], int(DL2[1]), int(DL2[2])
    CheckReads, CheckReadDic, NormalReads, NormalReadDic, DisRegion, DisReadDic = ReadParser(CheckLoci, Bam), {}, ReadParser(NormalLoci, Bam), {}, [], {'DisRead':0, 'NomReadOn':0}  # 
    for i in CheckReads: # BUG? CHECK? i = list or Str?
        RN, DL, FS = i[0], i, i[1]
        if FSDecimal(FS)[256] == 1: continue
        CheckReadDic[RN] = DL
    for i in NormalReads:
        RN, DL, FS = i[0], i, i[1]
        if FSDecimal(FS)[256] == 1: continue
        NormalReadDic[RN] = DL
    for i in ReadParsingList:
        DT = i.split('\t')
        ReadName, FlagStat, ReadChr, ReadStart, MappingQuality, CigarString, PairReadChr, PairReadPos = DT[0], DT[1], DT[2], int(DT[3]), DT[4], DT[5], DT[6], int(DT[7])
        if FSDecimal(FS)[256] == 1: continue
        if ReadName in CheckReadDic: PairInfo = CheckReadDic[ReadName]
        else: continue
        PairCigar = PairInfo[4]
        PairReadStart, PairReadEnd, ReadEnd = PairReadPos, PairReadPos + Cigar4Len(PairCigar) - 1, ReadStart + Cigar4Len(PairCigar) - 1
        if PairReadChr == CheckChr and CheckStart <= PairReadEnd and PairReadStart <= CheckEnd: 
            DisRegion.append([int(ReadStart), int(ReadEnd)])
            DisReadDic['DisRead'] += 1
            #M9_Sub_ReadLogger(ParsingLoci, CheckLoci, ReadName)
    DisMergeRegion = M7_OverlapMerge(DisRegion)
    DisLociStrList = []
    for i in DisMergeRegion: # Pos2Str
        DisChr, DisStart, DisEnd = Str2Pos(ParsingLoci)[0], i[0], i[1]
        DisLociStrList.append(Pos2Str(DisChr, DisStart, DisEnd))
    DisLociRead = ReadParser(' '.join(DisLociStrList), Bam)
    for i in DisLociRead:
        ReadName, FlagStat, ReadChr, ReadStart, MappingQuality, CigarString, PairReadChr, PairReadPos = DT[0], DT[1], DT[2], int(DT[3]), DT[4], DT[5], DT[6], int(DT[7])
        if FSDecimal(FS)[256] == 1: continue
        if ReadName in NormalReadDic: PairInfo = NormalReadDic[ReadName]
        else: continue
        PairCigar, PairPairPos = PairInfo[4], int(PairInfo[7])
        if PairReadPos == PairPairPos: continue # Same Read
        PairStart, PairEnd = PairReadPos, PairReadPos + Cigar4Len(PairCigar) - 1, 
        if PairStart <= NormalEnd and NormalStart <= PairEnd: 
            DisReadDic['NomReadOn'] += 1
            #M9_Sub_ReadLogger(FLoci, CheckLoci, ReadName)
    return DisReadDic


def M9_Sub_PairConnectionOVLP(FLoci, ReadParsingList, CheckLoci, Direction, DirectionMode): #Paired read in CheckLoci?
    DL = Str2Pos(CheckLoci)
    CheckChr, CheckStart, CheckEnd = DL[0], int(DL[1]), int(DL[2]) # Check Target
    CheckReads, CheckReadDic = ReadParser(CheckLoci, Bam), {}
    for i in CheckReads:
        RN, DL, FS = i[0], i, i[1]
        FSDecimalDic = FSDecimal(FS)
        if FSDecimalDic[256] == 1: continue
        CheckReadDic[RN] = DL #only primary mapping read
    Connection = 0
    for i in ReadParsingList:
        DT = i.split('\t')
        ReadName, FlagStat, MappingQuality, CigarString, PairReadChr, PairReadPos = DT[0], DT[1], DT[4], DT[5], DT[6], int(DT[7])
        # BUG : PairReadEnd Calculation; PairReadStart, PairReadEnd = PairReadPos, PairReadPPos + Cigar4Len(CigarString) - 1
        if PairReadChr == CheckChr and CheckStart <= PairReadEnd and PairReadStart <= CheckEnd: 
            PairInfo = CheckReadDic[ReadName]
            PairFlagStat, PairMappingQuality = PairInfo[1], PairInfo[4]
            if MQCutoff <= int(MappingQuality) and MQCutoff <= int(PairMappingQuality):
                if DirectionMode == "On":
                    if Direction == 'Same':
                        if FFRRFinder(FSDecimal(FlagStat), FSDecimal(PairFlagStat)) == 'Diff': 
                            Connection += 1
                            M9_Sub_ReadLogger(FLoci, CheckLoci, RN)
                    else: 
                        if FFRRFinder(FSDecimal(FlagStat), FSDecimal(PairFlagStat)) != 'Diff': 
                            Connection += 1
                            M9_Sub_ReadLogger(FLoci, CheckLoci, RN)
                else: 
                    Connection += 1
                    M9_Sub_ReadLogger(FLoci, CheckLoci, RN)
    return Connection
    

def M9_Misassigned_WholeContig_suppressed(FLoci, TLoci, Bam, FlankingSize): #Misassignment: ReadParsing on neighbor contigs and check connection; HE: Discordant read check
    FLociPos, TLociPos = Str2Pos(FLoci), Str2Pos(TLoci)
    FChr, FStart, FEnd = FLociPos[0], FLociPos[1], FLociPos[2]
    TChr, TStart, TEnd = TLociPos[0], TLociPos[1], TLociPos[2]
    ContigList = ContigDic[FChr].keys()
    Index = ContigList.index(Loci)
    LContig, RContig = ContigList[Index - 1], ContigList[Index + 1]
    LFlanking, RFlanking = int(LContig.split('-')[1]) - FlankingSize, int(RContig.split(':')[1].split('-')[0] + FlankingSize)
    LFFStart, LFFEnd, RFFStart, RFFEnd = LFlanking, LFlanking + FlankingSize, RFlanking - FlankingSize, RFlanking
    LFlanking, RFlanking = Pos2Str(FChr, LFFStart, LFFEnd), Pos2Str(FChr, RFFStart, RFFEnd)
    TFalnking = Pos2Str(TChr, TStart - FlankingSize, TEnd + FlankingSize)

    if Overlapping(LFFStart, LFFEnd, TStart, TEnd) > 0: Type = 'LOVLP'
    elif Overlapping(RFFStart, RFFEnd, TStart, TEnd) > 0: Type = 'ROVLP'
    else: Type = 'Away'
    if DirectionDic[FLoci] == DirectionDic[TLoci]: FTDirection = 'Same'
    else: FTDirection = 'Diff'

    if Type == 'LOVLP':
        RFlankingReads = ReadParser(RFlanking, Bam)
        Connection = M9_Sub_PairConnectionOVLPNormal(RFlanking, RFlankingReads, FLoci, FTDirection, "Off")
        if Connection == 0: return "FD", Connection
    elif Type == 'ROVLP':  
        LFlankingReads = ReadParser(LFlanking, Bam) 
        Connection = M9_Sub_PairConnectionOVLPNormal(LFlanking, LFlankingReads, FLoci, FTDirection, "Off")
        if Connection == 0: return "FD", Connection
    elif Type == 'Away':
        LFlankingReads, RFlankingReads = ReadParser(LFlanking, Bam), ReadParser(RFlanking, Bam)
        Connection = M9_Sub_PairConnectionOVLPNormal(RFlanking +'|'+ LFlanking, LFlankingReads + RFlankingReads, FLoci, FTDirection, "Off")
        if Connection == 0: return "FD", Connection

    FLociReads = ReadParser(FLoci, Bam) # Mis-assignment first, 
    Connection = M9_Sub_PairConnectionOVLPNormal(FLociReads, TLoci, FTDirection, "On")
    if Connection >0: return "FD", Connection
    return "nonFD", Connection


def M9_Misassigned_Edge(FLoci, TLoci, FDirection, TDirection, Type, FlankingSize, Bam): # SE: FFLoci to FT~TLoci connection check; HE: FLoci to FT~TLoci Connection Check
    FLociPos, TLociPos = Str2Pos(FLoci), Str2Pos(TLoci)
    FChr, FStart, FEnd, FSeqDirection = FLociPos[0], FLociPos[1], FLociPos[2], DirectionDic[FLoci] 
    TChr, TStart, TEnd, TSeqDirection = TLociPos[0], TLociPos[1], TLociPos[2], DirectionDic[TLoci]
    # 3 type of Edge: 1. Low depth (SE) 2. High depth (HE)    
    if FDirection == 'Left': FFStart, FFEnd = FStart - FlankingSize - 1, FStart - 1
    elif FDirection == 'Right': FFStart, FFEnd = FStart + 1, FEnd + FlankingSize + 1
    if TDirection == 'Left': FTStart, FTEnd = TStart - FlankingSize - 1, TStart - 1
    elif TDirection == 'Right': FTStart, FTEnd = TStart + 1, TEnd + FlankingSize + 1
    
    if DirectionDic[FLoci] == DirectionDic[TLoci]: FTDirection = 'Same'
    else: FTDirection = 'Diff'
    if Type == 'SE': 
        FFLoci = Pos2Str(FChr, FFStart, FFEnd)
        FFLociReads = ReadParser(FFLoci, Bam)
        Connection = M9_Sub_PairConnectionOVLPNormal(FFLoci, FFLociReads, Pos2Str(TChr, min(TStart, TEnd, FTStart, FTEnd), max(TStart, TEnd, FTStart, FTEnd)), FTDirection, "On")
    elif Type == 'HE': 
        FLociReads = ReadParser(FLoci, Bam)
        Connection = M9_Sub_PairConnectionOVLPNormal(FLoci, FLociReads, Pos2Str(TChr, min(TStart, TEnd, FTStart, FTEnd), max(TStart, TEnd, FTStart, FTEnd)), FTDirection, "On")
    if Connection == 0: return "nonFD", Connection
    else: return "FD", Connection

        
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
        if (len(M_Keys) - Step) % 10 == 0: print(len(M_Keys) - Step)

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
    with open(LogPath + 'Duplication/' + PID + 'Duplication.list', 'wb') as fw:
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
    file_list = os.listdir(LogPath + 'Duplication/')
    for i in file_list:
        with open(LogPath + 'Duplication/' + i, 'rb') as lr:
            LoadList = pickle.load(lr)
            TDVGP += LoadList[0]
            TDPre += LoadList[1]
            FDPre += LoadList[2]
            FDVGP += LoadList[3]
    
    with open(LogPath + 'Merged_Duplication.list','wb') as fw: pickle.dump([TDVGP, TDPre, FDPre, FDVGP], fw)
    return [TDVGP, TDPre, FDPre, FDVGP]

### [OUTTER-Loader]
def M_TDFDLociParsingbyDepth_Loader():
    with open(LogPath + 'Merged_Duplication.list','rb') as fr: 
        List = pickle.load(fr)
        return List 

### [OUTTER-Saver]
def MQList2File(List, Name):
    with open(LogPath + 'MQResult.txt', 'a') as fw:
        for i in List: # RegionTotal, Both60 Discord, Both60 normal, Both0 Discord, Both0 Normal
            print(i)
            fw.write(Name + '\t' + '\t'.join(map(str, i)) + '\n')

### [OUTTER-Main2]
def ReadMQParsing(LociList, BamPath):
    ReadInfList = [] # RegionTotal, Both60 Discord, Both60 normal, Both0 Discord, Both0 Normal; <20 == 0, >=20 == 60
    Step = 0
    for i in LociList[0:100]:
        #print(len(LociList) - Step)
        Step += 1
        print(Step)
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
        print("ReadParsing") 
        MaxDepth1 = sp.getoutput("samtools depth -r {0} {1} | awk 'BEGIN{{Max = 0}} $3 > Max {{Max = $3}} END{{print Max}}'".format(Loci1, BamPath))
        MaxDepth2 = sp.getoutput("samtools depth -r {0} {1} | awk 'BEGIN{{Max = 0}} $3 > Max {{Max = $3}} END{{print Max}}'".format(Loci2, BamPath))
        print('DepthCal')
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
def M9_MultiProc(ContigDic, BlockDic, InsertSizeCutoff, FlankingSize, ClusteringSize, BamPath, SECutoff, Cores, LogPath, New_DepthDic):
    Keys, DataLen = [], 0 
    for i in ContigDic: 
        Keys += list(ContigDic[i].keys())
        for j in ContigDic[i]:
            DD = ContigDic[i][j]
            if len(DD['FDPos']) > 0: DataLen += 1

    random.shuffle(Keys)
    Calling_Keys = np.array_split(Keys, Cores)
    Result = parmap.map(Main9_ContigFDDiscordance, Calling_Keys, ContigDic, BlockDic, InsertSizeCutoff, FlankingSize, ClusteringSize, BamPath, SECutoff, 'Multi', LogPath, New_DepthDic, pm_pbar = True)
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


def Main9_ContigFDDiscordance(Keys, ContigDic, BlockDic, InsertSizeCutoff, FlankingSize, ClusteringSize, BamPath, SECutoff, MultiProc, LogPath, New_DepthDic): #(ContigDic, BlockDic, BamPath, SECutoff, PADic)
    InContigDic, InBlockDic = ContigDic, BlockDic
    Step = 0
    if MultiProc == 'Multi': Keys = Keys                        ###
    elif MultiProc == 'Single': Keys = list(InContigDic.keys()) ### dosen't work now
    for i in Keys:                                              ###
        if (len(Keys) - Step) % 10 == 0: print(len(Keys) - Step)
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
                        DisDic = M9_Sub_PairConnectionOVLP4NormalFinal(TLoci, FLoci, InsertSizeCutoff, FlankingSize, BamPath, "Both", '-', LogPath)
                        if DisDic['OverlapSize'] > 0:
                            InContigDic[Chr][ContigLoci]['Pass'] += ClusterBlockList
                            InContigDic[Chr][ContigLoci]['Cluster'][str(ClusterBlockList)] = 'Pass'
                            for BlockNo in ClusterBlockList:    InBlockDic[BlockNo][-1]['Pass'].append(ContigLoci)
                        else: 
                            InContigDic[Chr][ContigLoci]['Filtered'] += ClusterBlockList
                            InContigDic[Chr][ContigLoci]['Cluster'][str(ClusterBlockList)] = 'Filtered'
                            for BlockNo in ClusterBlockList:    InBlockDic[BlockNo][-1]['Filtered'].append(ContigLoci)

                elif CFDType[0] == 'Edge':
                    DiscResult = M9_EdgePosbyCED(FDBlockList, BlockDic, ContigLoci, FlankingSize, InsertSizeCutoff, BamPath, LogPath, ContigDic) # Clustering and falnking finding by CED and resulting
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
                        DiscResult = M9_EdgePosbyCED(ClusterBlockList, BlockDic, ContigLoci, FlankingSize, InsertSizeCutoff, BamPath, LogPath, ContigDic)
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
                        DisDic = M9_Sub_PairConnectionOVLP4NormalFinal(TLoci, FLoci, InsertSizeCutoff, FlankingSize, BamPath, "Both", '-', LogPath) #Central
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
    with open(LogPath + '{0}CactusXPurgeDup_Filtered.txt'.format(Name),'w') as fw:
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


if __name__ == "__main__":
     
    #-preprocessing-
    ChrPairDic = ChrConversion()
    Fai2Bed()
    AGDGmerging(PreAG, PreDG, VGPAG, VGPDG, PreKmerDic, VGPKmerDic) # to bed file
    FaiDic = Fai2Dic()
    ContigPreDic = ContigbyAGDG('Pre')
    ContigVGPDic = ContigbyAGDG('VGP')
    BinSaving(ChrPairDic, 'ChrPairDic')
    BinSaving(FaiDic, 'FaiDic')
     
    #-MainStep-
    
    Main1 = Main1_CactusCalling(Maf) # return MafDic
    MainDic, DirectionDic = Main1[0], Main1[1]
    Maf, Main1 = 0, 0
    BinSaving(MainDic, 'MainDic')
    BinSaving(MainDic, 'DirectionDic')
     
    Main2 = Main2_ContigCalling(MainDic)
    MainDic, DepthInput, ContigPreListDic, ContigVGPListDic = Main2[0], Main2[1], Main2[2], Main2[3]
    BinSaving(MainDic, 'MainDic2')
    BinSaving(ContigPreListDic, 'ContigPreListDic')
    BinSaving(ContigVGPListDic, 'ContigVGPListDic')
    
    MainDic = BinLoading('MainDic2')
    ContigPreListDic = BinLoading('ContigPreListDic')
    ContigVGPListDic = BinLoading('ContigVGPListDic')
    Loaded_DepthDic = M3_Saved_DepthLoad()
    BinSaving(Loaded_DepthDic, 'Loaded_Depth')
    Main3 = Main3_DetphAdd2FDDic(MainDic, Loaded_DepthDic, Cores)
    MainDic, New_DepthDic = Main3[0], Main3[1]
    BinSaving(MainDic, 'MainDic_DepthLoad')
    BinSaving(New_DepthDic, 'New_DepthDic')
    
    MainDic = Main4_CactusFalseTrueAllocating(MainDic)
    BinSaving(MainDic, 'MainDic_Bool')
     
    MainDic = BinLoading('MainDic_Bool')
    MainFDList = Main56_DepthAndGapFiltering(MainDic)
    MainPreFDList, MainVGPFDList = MainFDList[0], MainFDList[1]
    BinSaving(MainPreFDList, 'MainPreFD')
    BinSaving(MainVGPFDList, 'MainVGPFD')
     
    MainPreFDList, MainVGPFDList, ContigPreListDic, ContigVGPListDic = BinLoading('MainPreFD'), BinLoading('MainVGPFD'), BinLoading('ContigPreListDic'), BinLoading('ContigVGPListDic')
    Main7Pre = Main7_CactusXPurgeDup(MainPreFDList, PrePurgeDupPath, ContigPreDic, ContigPreListDic)
    Main7VGP = Main7_CactusXPurgeDup(MainVGPFDList, VGPPurgeDupPath, ContigVGPDic, ContigVGPListDic) 
        
    MainPreContigDic, MainVGPContigDic = Main7Pre[0], Main7VGP[0]
    MainPreBlockDic, MainVGPBlockDic = Main7Pre[1], Main7VGP[1]
    PurgeBlockListPre, PurgeBlockListVGP = Main7Pre[2], Main7VGP[2]
    DGPreBlockDic, DGVGPBlockDic = Main7Pre[3], Main7VGP[3]
    BinSaving(MainPreContigDic, 'MainPreCD')
    BinSaving(MainVGPContigDic, 'MainVGPCD')
    BinSaving(MainPreBlockDic, 'MainPreBD')
    BinSaving(MainVGPBlockDic, 'MainVGPBD')
    BinSaving(PurgeBlockListPre, 'PurgePre')
    BinSaving(PurgeBlockListVGP, 'PurgeVGP')
    BinSaving(DGPreBlockDic, 'DGPreBlock')
    BinSaving(DGVGPBlockDic, 'DGVGPBlock')
    
    PurgeBlockListPre = BinLoading('PurgePre')
    PurgeBlockListVGP = BinLoading('PurgeVGP')
    MainPreContigDic = BinLoading('MainPreCD')
    MainVGPContigDic = BinLoading('MainVGPCD')
    MainPreBlockDic = BinLoading('MainPreBD')
    MainVGPBlockDic = BinLoading('MainVGPBD')
    
    PurgePafPreDic = PurgeDupPafLoad(Path, PrePurgePafPath)
    PurgePafVGPDic = PurgeDupPafLoad(Path, VGPPurgePafPath)

    PAPreDic, PAVGPDic = PurgeDupLociSpecify(PurgeBlockListPre, PurgePafPreDic), PurgeDupLociSpecify(PurgeBlockListVGP, PurgePafVGPDic)
    BinSaving(PAPreDic, 'PAPreDic')
    BinSaving(PAVGPDic, 'PAVGPDic')

    MainPreContigDic = Main8_ContigFDTypeAllocation(MainPreContigDic)
    MainVGPContigDic = Main8_ContigFDTypeAllocation(MainVGPContigDic)
    BinSaving(MainPreContigDic, 'MainPreCD+FDType')
    BinSaving(MainVGPContigDic, 'MainVGPCD+FDType')
      
    MainPreContigDic = BinLoading('MainPreCD+FDType')
    MainVGPContigDic = BinLoading('MainVGPCD+FDType')
    MainPreBlockDic = BinLoading('MainPreBD')
    MainVGPBlockDic = BinLoading('MainVGPBD')
    PAPreDic = BinLoading('PAPreDic')
    PAVGPDic = BinLoading('PAVGPDic')
    MainPreBlockDic = Main8_2_PurgeTLociTrans(MainPreBlockDic, PAPreDic)
    MainVGPBlockDic = Main8_2_PurgeTLociTrans(MainVGPBlockDic, PAVGPDic)
    BinSaving(MainPreBlockDic, 'Main8_2PreBD')
    BinSaving(MainVGPBlockDic, 'Main8_2VGPBD')
     
    MainPreContigDic = BinLoading('MainPreCD+FDType')
    MainVGPContigDic = BinLoading('MainVGPCD+FDType')
    MainPreBlockDic = BinLoading('Main8_2PreBD')
    MainVGPBlockDic = BinLoading('Main8_2VGPBD')

    f0, Maf, ContigPreDic, ContigVGPDic, PurgePafPreDic, PurgePafVGPDic = 0, 0, 0, 0, 0, 0
    Main1, MainDic, DirectionDic, Main2, DepthInput, ContigPreListDic, ContigVGPListDic = 0, 0, 0, 0, 0, 0, 0
    Main3, MainFDList, MainPreFDList, MainPreVGPFDList, Main7Pre, Main7VGP, PAPreDic, PAVGPDic = 0, 0, 0, 0, 0, 0, 0, 0
    
    New_DepthDic = BinLoading('New_DepthDic')
    Main9Pre = M9_MultiProc(MainPreContigDic, MainPreBlockDic, InsertSizeCutoff, FlankingSize, ClusteringSize, PreBamPath, PreSECutoff, Cores, LogPath, New_DepthDic)
    MainPreContigDic, MainPreBlockDic = Main9Pre[1], Main9Pre[0]
    BinSaving(MainPreContigDic, 'Main9PreContigDic')
    BinSaving(MainPreBlockDic, 'Main9PreBlockDic') 
     
    Main9VGP = M9_MultiProc(MainVGPContigDic, MainVGPBlockDic, InsertSizeCutoff, FlankingSize, ClusteringSize, VGPBamPath, VGPSECutoff, Cores, LogPath, New_DepthDic)
    MainVGPContigDic, MainVGPBlockDic = Main9VGP[1], Main9VGP[0]
    BinSaving(MainVGPContigDic, 'Main9VGPContigDic')
    BinSaving(MainVGPBlockDic, 'Main9VGPBlockDic')

    print("InputLoci" + '\t' + str(Main9Pre[2]))
    print("MergedLoci" + '\t' + str(Main9Pre[3]))
    print("InputLoci" + '\t' + str(Main9VGP[2]))
    print("MergedLoci" + '\t' + str(Main9VGP[3]))
     
    DGPreBlockDic = BinLoading('DGPreBlock')
    DGVGPBlockDic = BinLoading('DGVGPBlock')
    PreBD = {**BinLoading('Main9PreBlockDic'), **DGPreBlockDic} # being not merged area possibly?
    VGPBD = {**BinLoading('Main9VGPBlockDic'), **DGVGPBlockDic}
    Main10_BlockDic2txt(PreBD, 'Pre') #Prefix name
    Main10_BlockDic2txt(VGPBD, 'VGP')
    


    #Main1 = Main1_CactusCalling(Maf)
    #MainDic = Main1[0]
    #BinSaving(MainDic, 'Full_MainDic')
