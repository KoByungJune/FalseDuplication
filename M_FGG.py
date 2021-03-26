import os
import parmap
import numpy
import time
import pickle
import subprocess as sp
from M_FDFinder_Dic_HapDep_GapFilter_TDFDReadDistTest import * 
import re

if __name__=="__main__":
    ##### Operation Parameter #####
    '''         
    Path = '/home/ko/cactus/02/'
    AnnFile = "GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomiclargest_noEx"
    AnnFile2 = "GCF_003957565.1_bTaeGut1_v1.p_genomiclargest_noEx"
    HalName = 'Zebrafinch2' #except expansion
    Maf = 'Zebrafinch2_2'
    ChrNamePath = '/home/ko/cactus/02/ChrName.txt' #FullPath
    RefSeq = 'Zebrafinch1'
    TarSeq = 'Zebrafinch0'
    AT = 'Pre'
    ChrMafPath = '/disk2/home/ko/rawdata/Zebrafinch0.fna'
    PurgeChrMode = 'N2C'
    '''
    '''     
    Path = '/home/ko/cactus/04/'
    AnnFile = "GCF_000699085.1_ASM69908v1_genomiclargest_noEx"
    AnnFile2 = "GCF_003957555.1_bCalAnn1_v1.p_genomiclargest_noEx"
    HalName = 'Calann' #except expansion
    Maf = 'Calann'
    ChrNamePath = '/home/ko/cactus/04/ChrName_OldNew.txt' #FullPath
    RefSeq = 'Calann0'
    TarSeq = 'Calann2'
    AT = 'Pre'
    ChrMafPath = '/home/ko/rawdata/Calann1.fna'
    PurgeChrMode = 'C2N'
    '''
    '''    
    Path = '/home/ko/cactus/09_2/'
    AnnFile = "GCF_000002275.2_Ornithorhynchus_anatinus_5.0.1_genomiclargest_noEx"
    AnnFile2 = "GCF_004115215.1_mOrnAna1.p.v1_genomiclargest_noEx"
    HalName = 'Ornana_2' #except expansion
    Maf = 'Ornana_2'
    ChrNamePath = '/home/ko/cactus/09_2/ChrName_OldNew.txt' #FullPath
    RefSeq = 'Ornana0_1'
    TarSeq = 'Ornana1_1DM'
    AT = 'Pre'
    ChrMafPath = '/home/ko/rawdata/Ornana1_1DM.fna'
    PurgeChrMode = 'C2N'
    '''
    '''      
    Path = '/home/ko/cactus/02/'
    AnnFile2 = "GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomiclargest_noEx"
    AnnFile = "GCF_003957565.1_bTaeGut1_v1.p_genomiclargest_noEx"
    HalName = 'Zebrafinch2' #except expansion
    Maf = 'Zebrafinch2_2'
    ChrNamePath = '/home/ko/cactus/02/ChrName.txt' #FullPath
    TarSeq = 'Zebrafinch1'
    RefSeq = 'Zebrafinch0'
    AT = 'VGP'
    ChrMafPath = '/disk2/home/ko/rawdata/Zebrafinch1.fna'
    PurgeChrMode = 'N2C'
    '''
    '''    
    Path = '/home/ko/cactus/04/'
    AnnFile2 = "GCF_000699085.1_ASM69908v1_genomiclargest_noEx"
    AnnFile = "GCF_003957555.1_bCalAnn1_v1.p_genomiclargest_noEx"
    HalName = 'Calann' #except expansion
    Maf = 'Calann'
    ChrNamePath = '/home/ko/cactus/04/ChrName_OldNew.txt' #FullPath
    TarSeq = 'Calann0'
    RefSeq = 'Calann2'
    AT = 'VGP'
    ChrMafPath = '/home/ko/rawdata/Calann0_1.fna'
    PurgeChrMode = 'C2N'
    '''
    ''' 
    Path = '/home/ko/cactus/09_2/'
    AnnFile2 = "GCF_000002275.2_Ornithorhynchus_anatinus_5.0.1_genomiclargest_noEx"
    AnnFile = "GCF_004115215.1_mOrnAna1.p.v1_genomiclargest_noEx"
    HalName = 'Ornana_2' #except expansion
    Maf = 'Ornana_2'
    ChrNamePath = '/home/ko/cactus/09_2/ChrName_OldNew.txt' #FullPath
    TarSeq = 'Ornana0_1'
    RefSeq = 'Ornana1_1DM'
    AT = 'VGP'
    ChrMafPath = '/home/ko/rawdata/Ornana0_1.fna'
    PurgeChrMode = 'C2N'
    

    LogPath = Path + 'Log/'
    #f0 = open(Path + AnnFile + '.txt','r') #annotation 정보
    #f0_1 = open(Path + AnnFile2 + '.txt','r')
    f1 = open(LogPath + AT + "CactusXPurgeDup_Filtered_Type.txt", 'r')
    Cores = 40
    FGGCutoff = 0.5
    FDBlockList = f1.read().split('\n')
    if FDBlockList.count('') > 0: FDBlockList.remove('')
    f1.close()
    '''
#####
def BinSaving(PythonData, FileName):#c     
    with open(LogPath + FileName+'.Pybin', 'wb') as fb: pickle.dump(PythonData, fb)
         
def BinLoading(FileName):#c      
    with open(LogPath + FileName+'.Pybin', 'rb') as rb: return pickle.load(rb)

def List2File(List, FileName):
    with open(LogPath + FileName +'.txt', 'w') as fw:
        List = list(set(List))
        List.sort()
        for i in List: fw.write(str(i) + '\n')

def Dic2File(Dic, FileName):
    with open(LogPath + FileName + '.txt', 'w') as fw:
        for i in Dic: fw.write(str(i) + '\t' + str(Dic[i]) + '\n')

####Pre Processing ####
def ChrConversion(ChrPath):
    with open(ChrPath,'r') as fr:
        ChrPair = fr.read().split('\n')
        if ChrPair.count('') >0: ChrPair.remove('')
        ChrPairDic, ReChrPairDic = {}, {}
        for i in ChrPair:
            DT = i.split('\t')
            CM, NC = DT[0], DT[1]
            ChrPairDic[CM] = NC
            ChrPairDic[NC] = CM
        return ChrPairDic


def MafChrDic(Path, MafName, FnaPath):
    MafChrName = sp.getoutput("awk '$1 == \"s\" {{print $2}}' {0}{1}.maf | sort | uniq".format(Path, MafName))
    MafChrNameList = MafChrName.split('\n')
    if MafChrNameList.count('') > 0: MafChrNameList.remove('')
    MafChrNameList = MafChrNameList[2:] # what is mean of this.
    MafChrDic = {}
    MafChrNameList = list(set(MafChrNameList))
    print(len(MafChrNameList))
    for i in MafChrNameList:
        DT = i.split('.')
        if len(DT) == 3: CN, MCN = DT[1], DT[1] + '.' + DT[2]
        elif len(DT) >3: CN, MCN = DT[1], '.'.join(DT[1:])     #DT[1] + '.' + DT[2] + '.' + DT[3]
        MafChrDic[CN] = MCN
    with open(FnaPath, 'r') as fr: # data without maf
        FnaLine = fr.readlines()
        for i in FnaLine:
            if i[0] == '>':
                ChrName = i[1:].split('.')[0]
                if ChrName in MafChrDic: pass
                else:
                    MafChrDic[ChrName] = i[1:].split('\n')[0]
    return MafChrDic


def M_PurgeDupCacRef_Operator(FDBlockList, HalName, Maf, RefSeq, TarSeq, ChrDic, ChrMafDic, LogPath, Cores, PurgeChrMode):
    FDBlockList = FDBlockList
    random.shuffle(FDBlockList)
    input_data = np.array_split(FDBlockList, Cores)
    Result = parmap.map(PurgeDupCacRef, input_data, HalName, Maf, RefSeq, TarSeq, ChrDic, ChrMafDic, LogPath, PurgeChrMode, pm_pbar = True)
    MFDBlockList = []
    for i in Result: MFDBlockList += i
    return MFDBlockList


def PurgeDupCacRef(FDBlockList, HalName, Maf, RefSeq, TarSeq, ChrDic, ChrMafDic, LogPath, PurgeChrMode):
    InBlockList, Step, PID = [], 0, str(os.getpid())
    for i in FDBlockList:
        Step += 1
        if (len(FDBlockList) - Step) % 100 == 0: print(len(FDBlockList) - Step)
        DT = i.split('\t')
        #print(DT[1])
        #print(ChrDic[DT[1]])
        if PurgeChrMode == 'N2C': BlockNo, FChr, FStart, FEnd, TStart, TEnd, TChr, Tool, Left = DT[0], ChrMafDic[ChrDic[DT[1]]], int(DT[2]), int(DT[3]), int(DT[4]), int(DT[5]), DT[6], DT[0].split(' ')[0], DT[7:]
        if PurgeChrMode == 'C2N': BlockNo, FChr, FStart, FEnd, TStart, TEnd, TChr, Tool, Left = DT[0], ChrMafDic[DT[1]], int(DT[2]), int(DT[3]), int(DT[4]), int(DT[5]), DT[6], DT[0].split(' ')[0], DT[7:]
        if Tool == "Cactus": 
            InBlockList.append(DT)
            continue
        with open(Path + PID + 'Purge2Cac.bed', 'w') as Pur: 
            Pur.write('\t'.join([FChr, str(FStart - 1), str(FEnd)]))
            #print([FChr, str(FStart - 1), str(FEnd)])
        os.system("/home/vgp/hal/bin/hal2maf --noAncestors --maxBlockLen 10000000 --refGenome {2} --refTargets {0}{3}Purge2Cac.bed {0}{1}.hal {0}{3}Purge2Cac.maf".format(Path, HalName, TarSeq, PID))
        RefData = sp.getoutput("awk '{{if ($2 ~ {1}) {{ print $0 }} }}' {0}{2}Purge2Cac.maf".format(Path, '"' + RefSeq + '"', PID))
        RefData = RefData.split('\n')
        if RefData.count('') > 0: RefData.remove('')
        RefLoci = []
        for j in RefData:
            DT2 = j.split('\t')
            Chr, BedStart, Len, Direction, Total = ChrDic[DT2[1].split('.')[1]], int(DT2[2]), int(DT2[3]), DT2[4], DT2[5]
            if Direction == '+':
                Start = int(BedStart) + 1
                End = int(BedStart) + int(Len) #Start + int(Len) : This is cause of +1 End
            elif Direction == '-':
                Start = int(Total) - int(BedStart) - int(Len) + 1
                End = int(Total) - int(BedStart)
            RefLoci.append([Chr, Start, End])
        DT[8] = RefLoci
        InBlockList.append(DT)
        #print(DT)
    return InBlockList


#### TOOL ####
def Len_Cal(J, E):
    if E == 'B1':
        List = [int(J[6]), int(J[7]), int(J[8]), int(J[9])]
        if List[2] <= List[0] and List[1] <= List[3]:
            Length = List[1] - List[0] + 1
        else:
            List.remove(max(List))
            List.remove(min(List))
            Length = max(List) - min(List) + 1
    elif E == 'B2':
        List = [int(J[6]), int(J[7]), int(J[10]), int(J[11])]
        if List[2] <= List[0] and List[1] <= List[3]:
            Length = List[1] - List[0] + 1
        else:
            List.remove(max(List))
            List.remove(min(List))
            Length = max(List) - min(List) + 1
    elif E == 'B1B2': # CDS 1개 partial duplication -> 사실상 별로 의미 없음.
        List = [int(J[6]), int(J[7]), int(J[8]), int(J[9]), int(J[10]), int(J[11])]
        Length = 0
        if List[0] <= List[2] and List[5] <= List[1]:
            Length = List[3] - List[2] + List[5] - List[4] + 2
        else:
            Length = List[3] - List[0] + List[1] - List[4] + 2
    elif E == 'R1':
        List = [int(J[6]), int(J[7]), int(J[8]), int(J[9])]
        if List[2] <= List[0] and List[1] <= List[3]:
            Length = List[1] - List[0] + 1
        else:
            List.remove(max(List))
            List.remove(min(List))
            Length = max(List) - min(List) + 1
    return Length


def Annt2Dic(AnnList):
    AnnDic = {} # {Chr : data..}
    for i in AnnList:
        DT = i.split('\t')
        Chr = DT[0].split('.')[0]
        if Chr in AnnDic: AnnDic[Chr].append(i)
        else: AnnDic[Chr] = [i]
    return AnnDic

#### Main Start ####


def Main1_M_CDSOverlapping(FDBlockList, PreAnnFilePath, VGPAnnFilePath, Cores):
    input_data = FDBlockList
    random.shuffle(input_data)
    input_data = np.array_split(FDBlockList, Cores)
    input_data = [i.tolist() for i in input_data]
    PreAnntResult = parmap.map(M1_ABBCBBannotation, input_data, PreAnnFilePath, pm_pbar = True)
    RefAnntResult = parmap.map(M1_ABBCAannotation, input_data, VGPAnnFilePath, pm_pbar = True)
    MPreAnnt, MRefAnnt = [], []
    for i in PreAnntResult: MPreAnnt += i 
    for i in RefAnntResult: MRefAnnt += i
    return MPreAnnt + MRefAnnt

def M1_ABBCAannotation(FDBlock, AnnFilePath):
    with open(Path + AnnFilePath + '.txt','r') as fr:
        AnnVGP = fr.read().split('\n')
        if AnnVGP.count('') > 0: AnnVGP.remove('')
    AnnDic = Annt2Dic(AnnVGP)
    RefAnnt, Step  = [], 0
    for i in FDBlock:
        Step += 1
        if (len(FDBlock) - Step) % 100 == 0: print(len(FDBlock) - Step)
        DT = i
        BlockNo = DT[0]
        if DT[0].find('Cactus') > -1: RefList = [[DT[7], DT[8], DT[9]]]
        else: RefList = DT[8]
        for j in RefList:
            RefChr, RefStart, RefEnd = j[0], int(j[1]), int(j[2])
            if RefChr in AnnDic: AnnList = AnnDic[RefChr]
            else: AnnList = []
            for l in AnnList:
                R1, DT = 0, l.split('\t')
                CDSChr, CDSStart, CDSEnd, DetailTab = DT[0].split('.')[0], int(DT[1]), int(DT[2]), DT[5].split(';')
                if RefChr == CDSChr:
                    for Detail in DetailTab:
                        if Detail.find("GeneID:") > -1:  Gene = Detail.split(':')[1].split(',')[0]
                        elif Detail.split('=')[0] == 'product':  Product = Detail.split('=')[1]
                        elif Detail.split('=')[0] == 'gene':  Sym = Detail.split('=')[1]
                    if RefStart <= CDSEnd and CDSStart <= RefEnd:  R1 = 1
                    if R1 == 1: RefAnnt.append([BlockNo, RefChr, RefEnd - RefStart + 1, '-', 'R1', 'CDS', CDSStart, CDSEnd, RefStart, RefEnd, '-', '-', Gene, Product, Sym])
    return RefAnnt

def M1_ABBCBBannotation(FDBlock, AnnFilePath):
    with open(Path + AnnFilePath + '.txt','r') as fr:
        AnnPre = fr.read().split('\n')
        if AnnPre.count('') > 0: AnnPre.remove('')
    AnnDic = Annt2Dic(AnnPre)
    PreAnnt, Step = [], 0
    for i in FDBlock:
        Step += 1
        if (len(FDBlock) - Step) % 100 == 0: print(len(FDBlock) - Step)
        DT = i
        BlockNo, FD1Chr, FD2Chr, FD1Start, FD1End, FD2Start, FD2End = DT[0], DT[1], DT[6], int(DT[2]), int(DT[3]), int(DT[4]), int(DT[5])
        if FD1Chr == FD2Chr: 
            if FD1Chr in AnnDic: AnnList = AnnDic[FD1Chr]
            else: AnnList = []
        else: 
            if FD1Chr in AnnDic and FD2Chr in AnnDic: AnnList = AnnDic[FD1Chr] + AnnDic[FD2Chr]
            elif FD1Chr in AnnDic: AnnList = AnnDic[FD1Chr] # elif로 갔어야.. 20210208 >>>B1 죄다 씹힘. in hummingbird
            elif FD2Chr in AnnDic: AnnList = AnnDic[FD2Chr]
            else: AnnList = []
        for l in AnnList:
            B1, B2 = 0, 0
            DT = l.split('\t')
            CDSChr, CDSStart, CDSEnd, DetailTab = DT[0].split('.')[0], int(DT[1]), int(DT[2]), DT[5].split(';')
            if FD1Chr == CDSChr or FD2Chr == CDSChr:
                for Detail in DetailTab:
                    if Detail.find("GeneID:") > -1:  Gene = Detail.split(':')[1].split(',')[0]
                    elif Detail.split('=')[0] == 'product':  Product = Detail.split('=')[1]
                    elif Detail.split('=')[0] == 'gene':  Sym = Detail.split('=')[1]
                if FD1Chr == CDSChr:
                    if FD1Start <= CDSEnd and CDSStart <= FD1End:  B1 = 1
                if FD2Chr == CDSChr:
                    if FD2Start <= CDSEnd and CDSStart <= FD2End:  B2 = 1
                if B1 == 1 and B2 == 1:
                    PreAnnt.append([BlockNo, FD1Chr, FD1End-FD1Start + 1, FD2End-FD2Start + 1, 'B1B2', 'CDS', CDSStart, CDSEnd, FD1Start, FD1End, FD2Start, FD2End, Gene, Product, Sym, FD2Chr])
                elif B1 == 1:
                    PreAnnt.append([BlockNo, FD1Chr, FD1End-FD1Start + 1, FD2End-FD2Start + 1, 'B1', 'CDS', CDSStart, CDSEnd, FD1Start, FD1End, FD2Start, FD2End, Gene, Product, Sym, FD2Chr])
                elif B2 == 1:
                    PreAnnt.append([BlockNo, FD1Chr, FD1End-FD1Start + 1, FD2End-FD2Start + 1, 'B2', 'CDS', CDSStart, CDSEnd, FD1Start, FD1End, FD2Start, FD2End, Gene, Product, Sym, FD2Chr])
    return PreAnnt


def Extra_TEannotation(LogPath, FDBed, TEFilePath, Suffix, Cores):
    with open(LogPath + FDBed + '.bed','r') as fr:
        FDBedMerged = fr.read().split('\n')
        if FDBedMerged.count('') > 0: FDBedMerged.remove('')
    with open(LogPath + '../' + TEFilePath + '.txt','r') as fr:
        TEAnn = fr.read().split('\n')
        if TEAnn.count('') > 0: TEAnn.remove('')
    TEDic = Annt2Dic(TEAnn)
    random.shuffle(FDBedMerged)
    input_data = np.array_split(FDBedMerged, Cores)
    Result = parmap.map(Extra_TEFinding, input_data, TEDic, pm_pbar = True)
    TEFDDic = {}
    for i in Result:
        ML = i
        for j in ML:
            TENo, TEStart, TEEnd, OverlapLen, FDLoci, Name, Feature = j[0], int(j[1]), int(j[2]), j[3], j[4], j[5], j[6]
            TELen = TEEnd - TEStart + 1
            TEChr = FDLoci.split(':')[0]
            TEPos = TEChr + ':' + str(TEStart) + '-' + str(TEEnd)
            if TENo in TEFDDic: 
                TEFDDic[TENo]['FDOVLP'] += OverlapLen
                TEFDDic[TENo]['FDPos'].append(FDLoci)
            else: TEFDDic[TENo] = {'Len':TELen, 'FDOVLP':OverlapLen, 'FDPos':[FDLoci], 'TEPos':TEPos, 'Inf':Name + '\t' + Feature}
    with open(LogPath + Suffix + '.txt', 'w') as fw:
        for i in TEFDDic: fw.write(i + '\t' + str(TEFDDic[i]['FDOVLP'] / TEFDDic[i]['Len']) + '\t' + str(TEFDDic[i]['Len']) + '\t' + str(TEFDDic[i]['FDOVLP']) + '\t' + str(TEFDDic[i]['FDPos']) + '\t' + str(TEFDDic[i]['TEPos']) + '\t' + TEFDDic[i]['Inf'] + '\n')


def Extra_TEFinding(FDList, TEDic):
    TEFDList, Step = [], 0
    for i in FDList:
        Step += 1
        if (len(FDList) - Step) % 100 == 0: print(len(FDList) - Step)
        DT = i.split('\t')
        FChr, FStart, FEnd = DT[0], int(DT[1]) + 1, int(DT[2])
        FDLoci = Pos2Str(FChr, FStart, FEnd)
        if FChr in TEDic: TEList = TEDic[FChr]
        else: continue
        for l in TEList:
            DT = l.split('\t')
            TEChr, TEStart, TEEnd, Name, Feature, TENo = DT[0].split('.')[0], int(DT[1]), int(DT[2]), DT[3], DT[4], DT[5]
            if FStart <= TEEnd and TEStart <= FEnd:
                TELen = TEEnd - TEStart + 1
                OverlapLen = Overlapping(TEStart, TEEnd, FStart, FEnd)
                TEFDList.append([TENo, TEStart, TEEnd, OverlapLen, FDLoci, Name, Feature])
    return TEFDList


def SegFDCalling(FDBedMerged, SegDup):
    ReturnOVLP = [] #
    Step = 0
    for i in FDBedMerged:
        Step += 1
        if (len(FDBedMerged) - Step) % 100 == 0: print(len(FDBedMerged) - Step)
        DT = i.split('\t')
        FChr, FStart, FEnd = DT[0], int(DT[1]) + 1, int(DT[2])
        FDLoci = Pos2Str(FChr, FStart, FEnd)
        for j in SegDup:
            DT = j.split('\t')
            Chr, Start, End, Similarity, SegNo = DT[0].split('.')[0], int(DT[2]), int(DT[3]), round(float(DT[5]), 4), DT[6]
            SegLoci = Pos2Str(Chr, Start, End)
            if FChr == Chr and FStart <= End and Start <= FEnd:
                OVLPLen = Overlapping(FStart, FEnd, Start, End)
                ReturnOVLP.append([SegNo, SegLoci, FDLoci, OVLPLen])
    return ReturnOVLP


def Extra_SegmentalDuplication(LogPath, FDBed, SegDupFilePath, Cores): #
    with open(LogPath + FDBed + '.bed','r') as fr:
        FDBedMerged = fr.read().split('\n')
        if FDBedMerged.count('') > 0: FDBedMerged.remove('')
    with open(LogPath + '../' + SegDupFilePath + '.txt','r') as fr:
        SegDup = fr.read().split('\n')
        if SegDup.count('') > 0: SegDup.remove('')
    SegDic = {} # {SDNo:{SD1Pos:'', SD1FDOVLP:0, SD2Pos:'', SD2FDOVLP:0, Similarity:''}
    for i in SegDup:
        DT = i.split('\t')
        Chr, Start, End, Similarity, SegNo = DT[0].split('.')[0], int(DT[2]), int(DT[3]), round(float(DT[5]), 4), DT[6]
        if SegNo in SegDic: 
            SegDic[SegNo]['SD2Pos'] = Pos2Str(Chr, Start, End)
            SegDic[SegNo]['SD2OVLP'] = 0
        else: SegDic[SegNo] = {'SD1Pos':Pos2Str(Chr, Start, End), 'SD1OVLP':0, 'Similarity':Similarity, 'FDLoci':[]}
    input_data = numpy.array_split(FDBedMerged, Cores)
    Result = parmap.map(SegFDCalling, input_data, SegDup, pm_pbar = True)
    for i in Result:
        for j in i:
            SegNo, SegLoci, FDLoci, OVLPLen = j[0], j[1], j[2], j[3]
            SegDic[SegNo]['FDLoci'].append(FDLoci)
            if SegDic[SegNo]['SD1Pos'] == SegLoci: SegDic[SegNo]['SD1OVLP'] += OVLPLen
            elif SegDic[SegNo]['SD2Pos'] == SegLoci: SegDic[SegNo]['SD2OVLP'] += OVLPLen

    with open(LogPath + 'SegmentlaDuplicationFDOVLP.txt','w') as fw:
        for i in SegDic:
            print(i)
            SegNo, DD = i, SegDic[i]
            SD1Loci, SD2Loci, SD1OVLPLen, SD2OVLPLen, FDLoci = DD['SD1Pos'], DD['SD2Pos'], DD['SD1OVLP'], DD['SD2OVLP'], DD['FDLoci']
            SD1Pos, SD2Pos = Str2Pos(SD1Loci), Str2Pos(SD2Loci)
            SD1Len, SD2Len = int(SD1Pos[2]) - int(SD1Pos[1]) + 1, int(SD2Pos[2]) - int(SD2Pos[1]) + 1
            fw.write('\t'.join(list(map(str, [SegNo, SD1Pos, SD2Pos, SD1OVLPLen / SD1Len, SD2OVLPLen / SD2Len]))) + '\n')


def Extra_V1RFGG(LogPath, FGD, GeneList, SymbolDic, AnntList, GeneLenDic):
    with open(LogPath + FGD, 'rb') as fr: FGD = pickle.load(fr)
    with open(LogPath + SymbolDic, 'rb') as fr: SymbolDic = pickle.load(fr)
    with open(LogPath + AnntList, 'rb') as fr: Annt = pickle.load(fr)
    with open(LogPath + GeneList, 'r') as fr: GeneList = fr.read().split('\n')
    with open(LogPath + GeneLenDic, 'rb') as fr: GeneLenDic = pickle.load(fr)
    if GeneList.count('') > 0: GeneList.remove('')
    V1ROVLP = []
    for i in  FGD:
        Gene1, Gene2, DD = i[0], i[1], FGD[i]
        if Gene1 in SymbolDic: Symbol1 = SymbolDic[Gene1]
        else: Symbol1 = Gene1
        if Gene2 in SymbolDic: Symbol2 = SymbolDic[Gene2]
        else: Symbol2 = Gene2
        G1OVLPR, G2OVLPR, Type = DD['G1CDSOLenR'], DD['G2CDSOLenR'], DD['Type']
        for j in GeneList:
            if j == Gene1 and G1OVLPR > 0: V1ROVLP.append([j, i, (Symbol1, Symbol2), G1OVLPR, G2OVLPR, Type, 'G1OVLP'])
            elif j == Gene2 and G2OVLPR > 0: V1ROVLP.append([j, i, (Symbol1, Symbol2), G1OVLPR, G2OVLPR, Type, 'G2OVLP'])
    OVLPDic = {}
    for i in Annt:
        OVLPGene, BlockNo, CStart, CEnd, FDType = i[12], i[0], int(i[6]), int(i[7]), i[4]
        if FDType == 'B1': FStart, FEnd = int(i[8]), int(i[9])
        else: continue
        if GeneList.count(OVLPGene) > 0:
            Loci = [CStart, CEnd, FStart, FEnd]
            Loci.remove(max(Loci))
            Loci.remove(min(Loci))
            Symbol = SymbolDic[OVLPGene]
            if OVLPGene in OVLPDic: 
                OVLPDic[OVLPGene]['OVLPLoci'].append(Loci)
                OVLPDic[OVLPGene]['FDBlockNo'].append(BlockNo)
            else: OVLPDic[OVLPGene] = {'OVLPLoci':[Loci], 'FDBlockNo':[BlockNo]}
    with open(LogPath + 'V1RinPreFGD.txt', 'w') as fw:
        for i in V1ROVLP: fw.write('\t'.join(list(map(str, i))) + '\n')
        for i in OVLPDic: fw.write('\t'.join(list(map(str, [i, SymbolDic[i], str(Extra_NonBedNonOverlap_ListLen(M7_OverlapMerge(OVLPDic[i]['OVLPLoci'])))+ " / " + str(GeneLenDic[i]), Extra_NonBedNonOverlap_ListLen(M7_OverlapMerge(OVLPDic[i]['OVLPLoci'])) /  GeneLenDic[i],  str(OVLPDic[i]['FDBlockNo'])]))) + '\n')


def M1_C2N(FDBlockList, PurgeChrMode, ChrDic):
    FDBlockList_Conv = []
    if PurgeChrMode == 'C2N':
        for i in FDBlockList:
            DL = i
            Tool = i[0].split(' ')[0]
            if Tool == 'Cactus': DL[1], DL[6], DL[7] = ChrDic[DL[1]], ChrDic[DL[6]], ChrDic[DL[7]]
            else:  DL[1], DL[6] = ChrDic[DL[1]], ChrDic[DL[6]]
            FDBlockList_Conv.append(DL)
    else: FDBlockList_Conv = FDBlockList
    return FDBlockList_Conv


def Main2_BlockDicConstruction(FDBlock, AnntOverlapping):
    BlockDic, PreProductDic, PreSymbolDic, VGPProductDic, VGPSymbolDic = {}, {}, {}, {}, {}
    AffectedGene, AffectedCDS, Step = [], [], 0
    FDBlockDic = {}
    for i in FDBlock: FDBlockDic[i[0]] = i[1:]
    for i in AnntOverlapping:
        DL = i
        Step += 1
        print(len(AnntOverlapping) - Step)
        BlockNo, Chr, Type, Gene, Product, Sym = DL[0], DL[1], DL[4], DL[12], DL[13], DL[14]
        if len(i) == 16: Chr2 = DL[15] # F2
        else: pass # R1
        Overlap = Len_Cal(i, Type)
        if Type[0] == "R": VGPProductDic[Gene], VGPSymbolDic[Gene] = Product, Sym
        else: PreProductDic[Gene], PreSymbolDic[Gene] = Product, Sym
        CDSStart, CDSEnd, B1Start, B1End, B2Start, B2End = DL[6], DL[7], DL[8], DL[9], DL[10], DL[11]
        B1Pos, B2Pos, R1Pos = [CDSStart, CDSEnd, B1Start, B1End], [CDSStart, CDSEnd, B2Start, B2End], [CDSStart, CDSEnd, B1Start, B1End, Chr]
        if BlockNo in BlockDic:
            if Type == 'B1':
                AffectedGene.append(Gene) #B2는...? Bug..? 실수..? -> BUG 아님!! B1가 False duplication이니까.! -> 근데 else에 안넣음..
                AffectedCDS.append((Chr, CDSStart, CDSEnd))
                BlockDic[BlockNo]['B1CDS'].append(B1Pos)
                if Gene in BlockDic[BlockNo]['B1']:  BlockDic[BlockNo]['B1'][Gene] += Overlap
                else:  BlockDic[BlockNo]['B1'][Gene] = Overlap
            elif Type == 'B2':
                BlockDic[BlockNo]['B2CDS'].append(B2Pos)
                if Gene in BlockDic[BlockNo]['B2']:  BlockDic[BlockNo]['B2'][Gene] += Overlap
                else:  BlockDic[BlockNo]['B2'][Gene] = Overlap
            elif Type == 'B1B2': # One CDS Two Duplication Region
                if Gene in BlockDic[BlockNo]['B1B2']:  BlockDic[BlockNo]['B1B2'][Gene] += Overlap
                else:  BlockDic[BlockNo]['B1B2'][Gene] = Overlap
            elif Type == 'R1':
                BlockDic[BlockNo]['R1CDS'].append(R1Pos)
                if Gene in BlockDic[BlockNo]['R1']:  BlockDic[BlockNo]['R1'][Gene] += Overlap
                else:  BlockDic[BlockNo]['R1'][Gene] = Overlap
        else:
            if Type == 'B1':
                BlockDic[BlockNo] = {'B1':{Gene:Overlap},'B2':{},'B1B2':{},'R1':{}, 'B1CDS':[B1Pos], 'B2CDS':[], 'R1CDS':[]}
                AffectedGene.append(Gene) # 20201106
                AffectedCDS.append((Chr, CDSStart, CDSEnd))
            elif Type == 'B2':  BlockDic[BlockNo] = {'B1':{},'B2':{Gene:Overlap},'B1B2':{},'R1':{}, 'B1CDS':[], 'B2CDS':[B2Pos], 'R1CDS':[]}
            elif Type == 'B1B2':  BlockDic[BlockNo] = {'B1':{},'B2':{},'B1B2':{Gene:Overlap},'R1':{}, 'B1CDS':[], 'B2CDS':[], 'R1CDS':[]}
            elif Type == 'R1':  BlockDic[BlockNo] = {'B1':{},'B2':{},'B1B2':{},'R1':{Gene:Overlap}, 'B1CDS':[], 'B2CDS':[], 'R1CDS':[R1Pos]}
        DD = FDBlockDic[BlockNo]
        B1Chr, B1Start, B2Start, B2Chr, B2Start, B2End, RChr, RStart, REnd = DD[0], DD[1], DD[2], DD[5], DD[3], DD[4], DD[6], DD[7], DD[8]
        BlockDic[BlockNo]['B1LociPos'] = [B1Chr, B1Start, B1End] 
        BlockDic[BlockNo]['B2LociPos'] = [B2Chr, B2Start, B2End]
        BlockDic[BlockNo]['RLociPos'] = [RChr, RStart, REnd]
        if Type == 'B1':  BlockDic[BlockNo]['B1Loci'], BlockDic[BlockNo]['B2Loci'] = Chr, Chr2
    
    return BlockDic, PreProductDic, PreSymbolDic, VGPProductDic, VGPSymbolDic, AffectedGene, AffectedCDS


def M3_CDSDic(AnnFilePath):
    with open(Path + AnnFilePath + '.txt','r') as fr:
        Ann = fr.read().split('\n')
    if Ann.count('') > 0: Ann.remove('')
    CDSDic, ProductDic, SymDic, GeneLenDic = {}, {}, {}, {}
    for i in Ann:
        DT = i.split('\t')
        CDSChr, CDSStart, CDSEnd = DT[0].split('.')[0], int(DT[1]), int(DT[2])
        CDSLen = CDSEnd - CDSStart + 1
        DetailTab = DT[5].split(';')
        for Detail in DetailTab:
            if Detail.find("GeneID:") > -1:  Gene = Detail.split(':')[1].split(',')[0]
            elif Detail.split('=')[0] == 'product':  Product = Detail.split('=')[1]
            elif Detail.split('=')[0] == 'gene':  Sym = Detail.split('=')[1]
        if Gene in CDSDic: 
            CDSDic[Gene].append([CDSStart, CDSEnd])
            GeneLenDic[Gene] += CDSLen
        else: 
            CDSDic[Gene] = [[CDSStart, CDSEnd]]
            GeneLenDic[Gene] = CDSLen
        ProductDic[Gene], SymDic[Gene] = Product, Sym
    return CDSDic, ProductDic, SymDic, GeneLenDic


def M3_GeneDic(CDSDic):
    GeneDic = {}
    for i in CDSDic:
        GeneID, DD = i, CDSDic[i]
        GStart, GEnd = DD[0][0], DD[-1][1]
        GeneDic[GeneID] = [GStart, GEnd]
    return GeneDic


def M3_NOLCal(FGDDD, BlockNoDL, Target): #FGDDD = FGDDic[G1G2Pair]
    InDD = FGDDD
    for i in BlockNoDL:
         BCStart, BCEnd, FStart, FEnd = int(i[0]), int(i[1]), int(i[2]), int(i[3])
         for j in InDD[Target]:
             CDSNo = j
             DD = InDD[Target][CDSNo]
             CDSStart, CDSEnd = int(DD['Pos'][0]), int(DD['Pos'][1])
             if BCStart == CDSStart: 
                 if InDD[Target][CDSNo]['NOL'] == []: continue
                 Index = 0
                 for r in InDD[Target][CDSNo]['NOL']:
                     CRStart, CREnd = r[0], r[1]
                     if FStart <= CRStart and CREnd <= FEnd:
                         #print(InDD[Target][CDSNo]['NOL'])
                         #print(Index)
                         InDD[Target][CDSNo]['NOL'][Index] = ''
                     elif CRStart < FStart and FEnd < CREnd:
                         IndexPos = [[CRStart, FStart - 1], [FEnd + 1, CREnd]]
                         InDD[Target][CDSNo]['NOL'][Index] = ''
                         InDD[Target][CDSNo]['NOL'] += IndexPos
                         break
                     elif FStart <= CREnd and CRStart <= FEnd: #NOL segmenting
                         IndexPos = [0,0]
                         if FStart <= CRStart:  IndexPos[0] = FEnd + 1
                         else: IndexPos[0] = CRStart
                         if CREnd <= FEnd: IndexPos[1] = FStart - 1
                         else: IndexPos[1] = CREnd
                         InDD[Target][CDSNo]['NOL'][Index] = IndexPos
                     Index += 1
                 while InDD[Target][CDSNo]['NOL'].count('') != 0:  InDD[Target][CDSNo]['NOL'].remove('')
    return InDD


def Main3_BlockDic2FGPairDic(BlockDic, PreCDSDic, VGPCDSDic, PreProductDic, VGPProductDic, PreSymDic, VGPSymDic, PreGeneLenDic, VGPGeneLenDic, PreGeneDic, VGPGeneDic):
    FGDDic = {} # {(Gene1, Gene2)}: {..}
    for k in BlockDic:
        BlockNo, Tool = k, k.split(' ')[0]
        B1GeneList, B2GeneList, R1GeneList = list(BlockDic[BlockNo]['B1'].keys()), list(BlockDic[BlockNo]['B2'].keys()), list(BlockDic[BlockNo]['R1'].keys())
        DD = BlockDic[BlockNo]
        Change = 0
        #print(BlockNo + '\t' + str(BlockDic[k]))
        if len(B1GeneList) == 1 and len(B2GeneList) == 1:
            B1Gene, B2Gene, R1Gene = B1GeneList[0], B2GeneList[0], R1GeneList
            G1G2Pair = [int(B1Gene), int(B2Gene)]
            if G1G2Pair[0] <= G1G2Pair[1]: pass
            else:
                Change = 1
                G1G2Pair.sort()
            G1Gene, G2Gene = str(G1G2Pair[0]), str(G1G2Pair[1])
            G1G2Pair = (G1Gene, G2Gene)
            # Dictionary format initialization
            G1Product, G2Product = PreProductDic[B1Gene], PreProductDic[B2Gene]
            if G1G2Pair in FGDDic: pass
            else:
                FGDDic[G1G2Pair] = {'G1Prod':G1Product,'G2Prod':G2Product,'G1Total':PreGeneLenDic[G1Gene],'G2Total':PreGeneLenDic[G2Gene],'R1List':R1Gene, 'G2List':G2Gene, 'G1CDS':{}, 'G2CDS':{}, 'R1CDS':{}, 'BlockList':[], 'G1Chr':'', 'G2Chr':''}
                CDSName = 1
                for i in PreCDSDic[G1Gene]:
                    FGDDic[G1G2Pair]['G1CDS'][str(CDSName)] = {'Pos':i,'NOL':[[i[0], i[1]]]}
                    CDSName += 1
                CDSName = 1
                for i in PreCDSDic[G2Gene]:
                    FGDDic[G1G2Pair]['G2CDS'][str(CDSName)] = {'Pos':i,'NOL':[[i[0], i[1]]]}
                    CDSName += 1
                if len(R1GeneList) == 0:
                    FGDDic[G1G2Pair]['R1CDS'] = 'ND'
                    continue
                else: # Ref is used in FCG, FCG = (B1 != B2, B1 != Ref, B2 == Ref; therefore G1allocating first.)
                    if len(R1GeneList) == 1: R1Gene = R1GeneList[0]
                    elif G1Gene == G2Gene:
                        if G1Gene in R1Gene: R1Gene = G1Gene
                        else: 
                            FGDDic[G1G2Pair]['R1CDS'] = 'ND' # non-homologs
                            continue
                    elif G1Gene != G2Gene:
                        if G1Gene in R1GeneList and G2Gene in R1GeneList: # Purge 2 Cactus multi-homologs
                            if DD['R1'][G1Gene] > DD['R1'][G2Gene]: R1Gene = G1Gene # R1Gene specifying
                            else: R1Gene = G2Gene
                        elif G1Gene in R1GeneList: R1Gene = G1Gene
                        elif G2Gene in R1GeneList: R1Gene = G2Gene
                        else: 
                            FGDDic[G1G2Pair]['R1CDS'] = 'ND' # non-homologs
                            continue
                    
                    CDSName = 1
                    for i in VGPCDSDic[R1Gene]:
                        FGDDic[G1G2Pair]['R1CDS'][str(CDSName)] = {'Pos':i,'NOL':[[i[0], i[1]]]}
                        CDSName += 1

            # OVLP Data to Gene Dictionary
            if G1G2Pair in FGDDic:
                FGDDic[G1G2Pair]['BlockList'].append(BlockNo)
                if Change == 0:
                    FGDDic[G1G2Pair]['G1Chr'], FGDDic[G1G2Pair]['G2Chr'] = BlockDic[BlockNo]['B1Loci'], BlockDic[BlockNo]['B2Loci']
                    DL = BlockDic[BlockNo]['B1CDS']
                    FGDDic[G1G2Pair] = M3_NOLCal(FGDDic[G1G2Pair], DL, 'G1CDS')
                elif Change == 1: # Data positioning
                    FGDDic[G1G2Pair]['G1Chr'], FGDDic[G1G2Pair]['G2Chr'] = BlockDic[BlockNo]['B2Loci'], BlockDic[BlockNo]['B1Loci']
                    DL = BlockDic[BlockNo]['B1CDS'] #right == B1 = false , B2 = true(except to save)
                    FGDDic[G1G2Pair] = M3_NOLCal(FGDDic[G1G2Pair], DL, 'G2CDS')

                DL = BlockDic[BlockNo]['R1CDS'] #R1 Data inputting
                if DL == 'ND' or FGDDic[G1G2Pair]['R1CDS'] == 'ND': pass
                else: FGDDic[G1G2Pair] = M3_NOLCal(FGDDic[G1G2Pair], DL, 'R1CDS')
        
        elif (len(B1GeneList) > 1 or len(B2GeneList) > 1) and (len(B1GeneList) != 0 and len(B2GeneList) != 0): # one to many ~ many to many => only for False position
            for j in B1GeneList:
                B1Gene = j
                GPos = PreGeneDic[B1Gene]
                GStart, GEnd = int(GPos[0]), int(GPos[1])
                FStart, FEnd = int(BlockDic[BlockNo]['B1LociPos'][1]), int(BlockDic[BlockNo]['B2LociPos'][2])
                if FStart <= GStart and GEnd <= FEnd:   G1G2Pair = (B1Gene, str('Whole')) # >>> FGG
                else:                                   
                    G1G2Pair = (B1Gene, str('Partial')) # >>> if B1Gene in G2GeneList: FED; elif No B1Gene in R1Gene List: FCG
                G1Product = PreProductDic[B1Gene]
                FGDDic[G1G2Pair] = {'G1Prod':G1Product,'G2Prod':'Many','G1Total':PreGeneLenDic[B1Gene],'G2Total':'Many', 'R1List':R1GeneList, 'G2List':B2GeneList, 'G1CDS':{}, 'G2CDS':'Many', 'R1CDS':'Many', 'BlockList':[], 'G1Chr':'', 'G2Chr':''}
                CDSName = 1
                for i in PreCDSDic[B1Gene]:
                    FGDDic[G1G2Pair]['G1CDS'][str(CDSName)] = {'Pos':i,'NOL':[[i[0], i[1]]]}
                    CDSName += 1

                FGDDic[G1G2Pair]['BlockList'].append(BlockNo)
                FGDDic[G1G2Pair]['G1Chr'], FGDDic[G1G2Pair]['G2Chr'] = BlockDic[BlockNo]['B1Loci'], BlockDic[BlockNo]['B2Loci']
                DL = BlockDic[BlockNo]['B1CDS']
                FGDDic[G1G2Pair] = M3_NOLCal(FGDDic[G1G2Pair], DL, 'G1CDS') # G2 position = whole or partial overlap of a gene
    return FGDDic
        

def M4_OVLPAdd(SubDD, NewDD, VN): #
    InDD = NewDD
    for j in SubDD:
        CDSNo, Start, End, RemainingList = j, SubDD[j]['Pos'][0], SubDD[j]['Pos'][1], SubDD[j]['NOL']
        CDSLen = int(End) - int(Start) + 1
        InDD[VN + 'CDSLen'] += CDSLen
        if len(RemainingList) == 0: # FD fully overlap
            InDD[VN + 'WFEN'] += 1
            InDD[VN + 'CDSOL'] += CDSLen                
            InDD[VN + 'EPosbyO'].append(j) # 20201012 4 FD exon pos
            if len(SubDD) == 1: pass
            else:
                InDD[VN + 'EPos'].append((int(j) - 1) / (len(SubDD) - 1))
        else:
            RemainingLen = 0
            for k in RemainingList:
                RStart, REnd = k[0], k[1]
                RemainingLen += int(REnd) - int(RStart) + 1
            InDD[VN + 'CDSOL'] += CDSLen - RemainingLen
    return InDD


def M4_TextSimilarity4FCG(Prod1, Prod2):
    ZeroInf = ['LOW', 'QUALITY', 'PROTEIN:','like',' ']
    Words1 = re.split(' |-', Prod1)
    Words2 = re.split(' |-', Prod2)
    for i in ZeroInf:
        while Words1.count(i) >0: Words1.remove(i)
        while Words2.count(i) >0: Words2.remove(i)

    Intersection = 0
    for i in Words1:
        if Words2.count(i) >0: Intersection += 1
    for i in Words2:
        if Words1.count(i) >0: Intersection += 1
    Similarity = Intersection / (len(Words1) + len(Words2))
    return Similarity


def M4_TextSimilarity4FCG_print(Prod1, Prod2):
    ZeroInf = ['LOW', 'QUALITY', 'PROTEIN:','like',' ']
    Words1 = re.split(' |-', Prod1)
    Words2 = re.split(' |-', Prod2)
    for i in ZeroInf:
        while Words1.count(i) >0: Words1.remove(i)
        while Words2.count(i) >0: Words2.remove(i)

    Intersection = 0
    for i in Words1:
        if Words2.count(i) >0: Intersection += 1
    for i in Words2:
        if Words1.count(i) >0: Intersection += 1
    Similarity = Intersection / (len(Words1) + len(Words2))
    return Similarity, Words1, Words2, Intersection




def Main4_FGDTyper(FGDDic, FGGCutoff):
    FGD, Textmining, FCGProduct = {}, [], [] # FGDDic2
    for i in FGDDic:
        # Dictionary formatting, Data inputing
        Pair, DD = i, FGDDic[i]
        G1, G2, R1List, G2List, BlockList = i[0], i[1], DD['R1List'], DD['G2List'], DD['BlockList']
        G1Product, G2Product = DD['G1Prod'], DD['G2Prod']   
        G1CDSDic, G2CDSDic, R1CDSDic = DD['G1CDS'], DD['G2CDS'], DD['R1CDS']
        G1CDSNo, G2CDSNo, R1CDSNo = len(G1CDSDic), len(G2CDSDic), len(R1CDSDic)
        FGD[i] = {'G1EN':G1CDSNo, 'G2EN':G2CDSNo, 'R1EN':R1CDSNo, 'G1CDSLen':0, 'G2CDSLen':0, 'R1CDSLen':0, 'G1CDSOL':0, 'G2CDSOL':0, 'R1CDSOL':0, 'G1CDSOLenR':0, 'G2CDSOLenR':0, 'R1CDSOLenR':0, 'G1WFEN':0, 'G2WFEN':0, 'R1WFEN':0,'Type':'', 'G1FER':0, 'G2FER':0, 'R1FER':0,'VGP':'', 'Evenness':0, 'G1EPos':[], 'G2EPos':[], 'R1EPos':[], 'G1EPosbyO':[], 'G2EPosbyO':[], 'R1EPosbyO':[]}
        FGD[i] = M4_OVLPAdd(G1CDSDic, FGD[i], 'G1')
        if G2 == 'Whole' or G2 == "Partial":
            FGD[i]['G2EN'], FGD[i]['G2CDSLen'], FGD[i]['G2CDSOL'], FGD[i]['G2WFEN'], FGD[i]['G2CDSOLenR'], FGD[i]['G2FER'], FGD[i]['G2EPos'], FGD[i]['G2EPosbyO'] = 'ND', 'ND', 'ND', 'ND', 'ND', 'ND', 'ND', 'ND'
        else: FGD[i] = M4_OVLPAdd(G2CDSDic, FGD[i], 'G2')
        if R1CDSDic == 'ND':
            FGD[i]['R1CDSLen'], FGD[i]['R1WFEN'], FGD[i]['R1CDSOL'], FGD[i]['R1EPos'] = 'ND', 'ND', 'ND', 'ND'
        elif G2 == 'Whole' or G2 == 'Partial':
            FGD[i]['R1CDSLen'], FGD[i]['R1WFEN'], FGD[i]['R1CDSOL'], FGD[i]['R1EPos'] = 'Many', 'Many', 'Many', 'Many'
        else:
            FGD[i] = M4_OVLPAdd(R1CDSDic, FGD[i], 'R1')
       
        # Data Calculation #
        FGD[i]['G1CDSOLenR'] = FGD[i]['G1CDSOL'] / FGD[i]['G1CDSLen']
        if G2 != 'Whole' and G2 != 'Partial': FGD[i]['G2CDSOLenR'] = FGD[i]['G2CDSOL'] / FGD[i]['G2CDSLen']
        if FGD[i]['R1CDSLen'] == 'ND': FGD[i]['R1CDSOLenR'] = 'ND'
        elif FGD[i]['R1CDSLen'] == 'Many': FGD[i]['R1CDSOLenR'] = 'Many'
        else: FGD[i]['R1CDSOLenR'] = FGD[i]['R1CDSOL'] / FGD[i]['R1CDSLen']
        FGD[i]['G1FER'] = FGD[i]['G1WFEN'] / FGD[i]['G1EN']
        if G2 != 'Whole' and G2 != 'Partial': FGD[i]['G2FER'] = FGD[i]['G2WFEN'] / FGD[i]['G2EN']
        if FGD[i]['R1CDSLen'] == 'ND': FGD[i]['R1FER'] = 'ND'
        elif FGD[i]['R1CDSLen'] == 'Many': FGD[i]['R1FER'] = 'Many'
        else: FGD[i]['R1FER'] = FGD[i]['R1WFEN'] / FGD[i]['R1EN']
        if G2 != 'Whole' and G2 != 'Partial': FGD[i]['Evenness'] = 1 - abs(FGD[i]['G1FER'] - FGD[i]['G2FER'])
        #### FGD classification
        # 0. Anyoverlap                                                >>> Affected Gene    Clear.
        # 1. [Whole exon conver >0] 
        # 2. G1 == G2                                                  >>> FED              Clear.
        # 3. G1P == G2P; G > 50% of CDSLen                             >>> FGG              Cutoff..?; Same Function definition..?
        # 4. (G1 or G2) == 100%                                        >>> FGG              Clear.
        # 5. G1P != G2P; (G1P or G2P) != R1P and (G1P or G2P == R1P)   >>> FCG              Best.
        # 6. G2 == Whole                                               >>> FGG
        # 7. G2 == Partial and G1 in G2                                >>> FED
        # 8. G2 == Partial and G1 !in R1                               >>> FCG
       
        Textmining.append([(G1Product, G2Product), M4_TextSimilarity4FCG(G1Product, G2Product)])
        ProductSimilarity = M4_TextSimilarity4FCG(G1Product, G2Product)
        G1ProductProc, G2ProductProc = G1Product.split('-like')[0], G2Product.split('-like')[0]
        if G1ProductProc[-1] == ' ': G1ProductProc[:-1]
        if G2ProductProc[-1] == ' ': G2ProductProc[:-1]
        if G2 != 'Whole' and G2 != 'Partial':
            if G1 == G2:
                if FGD[i]['G1WFEN'] > 0 or FGD[i]['G2WFEN'] > 0: FGD[i]['Type'] = 'FED' # >>> FED
            elif ProductSimilarity > 0.5:
                if FGD[i]['G1WFEN'] > 0 and FGD[i]['G2WFEN'] > 0: 
                    if (FGD[i]['G1CDSOL'] + FGD[i]['G2CDSOL']) / min(FGD[i]['G1CDSLen'], FGD[i]['G2CDSLen']) > FGGCutoff: FGD[i]['Type'] = 'MixedFGG' # >>> G1G2FGG >>> 20210208 method changed
                elif FGD[i]['G1WFEN'] > 0 and FGD[i]['G1CDSOLenR'] > FGGCutoff: FGD[i]['Type'] = 'G1FGG' # >>> G1FGG
                elif FGD[i]['G2WFEN'] > 0 and FGD[i]['G2CDSOLenR'] > FGGCutoff: FGD[i]['Type'] = 'G2FGG' # >>> G2FGG
            else:
                if FGD[i]['G1WFEN'] > 0 or FGD[i]['G2WFEN'] > 0:
                    if FGD[i]['G1FER'] == 1: FGD[i]['Type'] = 'G1FGG' # 100% overlap >>> G1FGG
                    elif FGD[i]['G2FER'] == 1: FGD[i]['Type'] = 'G2FGG' # 100% overlap >>> G2FGG
                    else:
                        FCGProduct.append((i, G1Product,G2Product, ProductSimilarity))
                        if FGD[i]['G1WFEN'] > 0 and R1List.count(G1) == 0: FGD[i]['Type'] = 'G1FCG' # >>> G2FCG
                        elif FGD[i]['G2WFEN'] > 0 and R1List.count(G2) == 0: FGD[i]['Type'] = 'G2FCG' # >>> G1FCG
        else:
            if G2 == 'Whole': FGD[i]['Type'] = 'G1FGG' # >>> G1FGG
            elif FGD[i]['G1WFEN'] > 0 and G1 in G2List: FGD[i]['Type'] = 'FED' # >>> FED
            elif FGD[i]['G1WFEN'] > 0 and R1List.count(G1) == 0: 
                FunctionSame = 0
                for j in G2List:
                    G2Product = PreProductDic[j] # VGPProductDic -> PreProductDic 20210205
                    ProductSimilarity = M4_TextSimilarity4FCG(G1Product, G2Product)
                    if ProductSimilarity > 0.5:
                        FunctionSame = 1
                        break
                if FunctionSame == 0: 
                    FGD[i]['Type'] = 'G1FCG' # >>> G1FCG
                    FCGProduct.append((i, G1Product, G2Product, ProductSimilarity))
            elif FGD[i]['G1WFEN'] > 0 and FGD[i]['G1CDSOLenR'] > FGGCutoff: FGD[i]['Type'] = 'G1FGG' # >>> G1FGG
    return FGD, Textmining, FCGProduct


def Main5_FGDResulting(FGD, AT):
    FGDTypeDic = {'FGG':[], 'FED':[], 'FCG':[], 'Combination':{}, 'PhasingError':{}}
    for i in FGD:
        G1, G2 = i[0], i[1]
        if FGD[i]['Type'] == 'MixedFGG':
            FGDTypeDic['PhasingError'][(G1, G2)] = [FGD[i]['G1EPosbyO'], FGD[i]['G2EPosbyO']]
            FGDTypeDic['FGG'].append(G1)
            FGDTypeDic['FGG'].append(G2)
        if FGD[i]['Type'] == 'G1FGG': FGDTypeDic['FGG'].append(G1)
        if FGD[i]['Type'] == 'G2FGG': FGDTypeDic['FGG'].append(G2)
        if FGD[i]['Type'] == 'G1FCG': FGDTypeDic['FCG'].append(G1)
        if FGD[i]['Type'] == 'G2FCG': FGDTypeDic['FCG'].append(G2)
        if FGD[i]['Type'] == 'FED': FGDTypeDic['FED'].append(G1)
    
    TotalGeneSet = list(set(FGDTypeDic['FGG'] + FGDTypeDic['FCG'] + FGDTypeDic['FED']))
    for i in TotalGeneSet:
        Gene, Type = i, []
        if Gene in FGDTypeDic['FGG']: Type.append('FGG')
        if Gene in FGDTypeDic['FCG']: Type.append('FCG')
        if Gene in FGDTypeDic['FED']: Type.append('FED')
        if len(Type) > 1: FGDTypeDic['Combination'][i] = Type
    with open(LogPath + AT + 'FGDStat.txt','w') as fw:
        #Row = i <= ?????
        for i in FGDTypeDic: 
            if i == 'Combination' or i == 'PhasingError': fw.write(i + '\t' + str(FGDTypeDic[i]) + '\n')
            else: fw.write(i + '\t' + str(len(list(set(FGDTypeDic[i])))) + '\t' + '\t'.join(list(set(FGDTypeDic[i]))) + '\n')
    return FGDTypeDic


def Extra_CDSOverlap(CDSDic, Annt):
    GCNOL = {}
    for i in CDSDic:
        Gene, Pos  = i, CDSDic[i]
        GCNOL[Gene] = {}
        Index = 1
        for j in Pos:
            Start, End = int(j[0]), int(j[1])
            GCNOL[Gene][str(Index)] = {'Pos':[Start, End], 'NOL':[[Start, End]]}
            Index += 1
    for i in Annt:
        Gene, Type, FStart, FEnd = i[12], i[4], int(i[8]), int(i[9])
        if Type != 'B1': continue
        DDL = GCNOL[Gene]
        for j in DDL: # j =
            CDSNo, DD = j, DDL[j]
            if len(DD['NOL']) == 0: continue
            else:
                Index = 0
                for k in DD['NOL']:
                    NStart, NEnd = k[0], k[1]
                    if FStart <= NStart and NEnd <= FEnd:
                        GCNOL[Gene][CDSNo]['NOL'][Index] = ''
                    elif NStart < FStart and FEnd < NEnd:
                        IndexPos = [[NStart, FStart - 1], [FEnd + 1, NEnd]]
                        GCNOL[Gene][CDSNo]['NOL'][Index] = ''
                        GCNOL[Gene][CDSNo]['NOL'] += IndexPos
                        break
                    elif FStart <= NEnd and NStart <= FEnd:
                        IndexPos = [0, 0]
                        if FStart <= NStart: IndexPos[0] = FEnd + 1
                        else: IndexPos[0] = NStart
                        if NEnd <= FEnd: IndexPos[1] = FStart - 1
                        else: IndexPos[1] = NEnd
                        GCNOL[Gene][CDSNo]['NOL'][Index] = IndexPos
                    Index ++ 1
                while GCNOL[Gene][CDSNo]['NOL'].count('') != 0: GCNOL[Gene][CDSNo]['NOL'].remove('')
    WholeCDSFDList, WholeCDSFDRatioDic = [], {}
    for i in GCNOL:
        Gene, DD = i , GCNOL[i]
        GeneCDSNo = len(DD)
        WholeCDSFDNo = 0
        for j in DD:
            CDSNo, Pos, NOL = j, DD[j]['Pos'], DD[j]['NOL']
            CDSLoci = Pos2Str(Gene, Pos[0], Pos[1])
            if len(NOL) >0: continue 
            WholeCDSFDList.append(CDSLoci)
            WholeCDSFDNo += 1
        WholeCDSFDRatio = WholeCDSFDNo / len(DD)
        WholeCDSFDRatioDic[Gene] = WholeCDSFDRatio
    return GCNOL, WholeCDSFDList, WholeCDSFDRatioDic


#####################################
## Preprocessing
if __name__=="__main__":
    '''
    ChrDic = ChrConversion(ChrNamePath)
    ChrMafDic = MafChrDic(Path, Maf, ChrMafPath)
    BinSaving(ChrMafDic, 'ChrMafDic')
    FDBlockList = M_PurgeDupCacRef_Operator(FDBlockList, HalName, Maf, RefSeq, TarSeq, ChrDic, ChrMafDic, LogPath, Cores, PurgeChrMode)
    FDBlockList = M1_C2N(FDBlockList, PurgeChrMode, ChrDic) # Chr Conv
    BinSaving(FDBlockList, AT + 'PurgeTLoci_FDBlockList')

    ## Main Start
    FDBlockList = BinLoading(AT + 'PurgeTLoci_FDBlockList')
    AnntDat = Main1_M_CDSOverlapping(FDBlockList, AnnFile, AnnFile2, Cores)
    BinSaving(AnntDat, AT + 'AnntDat')
    
    FDBlockList = BinLoading(AT + 'PurgeTLoci_FDBlockList')
    AnntDat = BinLoading(AT + 'AnntDat')
    Main2 = Main2_BlockDicConstruction(FDBlockList, AnntDat)
    BlockDic, PreProductDic, PreSymbolDic, VGPProductDic, VGPSymbolDic, AffectedGene, AffectedCDS = Main2[0], Main2[1], Main2[2], Main2[3], Main2[4], Main2[5] , Main2[6]
    BinSaving(BlockDic, AT + 'FGD_Main2_BlockDic')
    BinSaving(AffectedGene, AT + 'FGD_Main2_AffectedGeneList')
    BinSaving(AffectedCDS, AT + 'FGD_Main2_AffectedCDSList')
    List2File(AffectedGene, AT + 'FGD_Main2_AffectedGeneList')
    List2File(AffectedCDS, AT + 'FGD_Main2_AffectedCDSList')

    BlockDic = BinLoading(AT + 'FGD_Main2_BlockDic')

    AnnPreDic, AnnVGPDic = M3_CDSDic(AnnFile), M3_CDSDic(AnnFile2)
    PreCDSDic, PreProductDic, PreSymDic, PreGeneLenDic = AnnPreDic[0], AnnPreDic[1], AnnPreDic[2], AnnPreDic[3]
    VGPCDSDic, VGPProductDic, VGPSymDic, VGPGeneLenDic = AnnVGPDic[0], AnnVGPDic[1], AnnVGPDic[2], AnnVGPDic[3]
    BinSaving(PreCDSDic, AT +  'FGD_PreCDSDic')
    BinSaving(PreProductDic, AT + 'FGD_PreProductDic')
    BinSaving(PreSymDic, AT + 'FGD_PreSymDic')
    BinSaving(PreGeneLenDic, AT + 'FGD_PreGeneLenDic')
    BinSaving(VGPCDSDic, AT + 'FGD_VGPCDSDic')
    BinSaving(VGPProductDic, AT + 'FGD_VGPProductDic')
    BinSaving(VGPSymDic, AT + 'FGD_VGPSymDic')
    BinSaving(VGPGeneLenDic, AT + 'FGD_VGPGeneLenDic')

    PreGeneDic, VGPGeneDic = M3_GeneDic(PreCDSDic), M3_GeneDic(VGPCDSDic)
    FGDDic = Main3_BlockDic2FGPairDic(BlockDic, PreCDSDic, VGPCDSDic, PreProductDic, VGPProductDic, PreSymDic, VGPSymDic, PreGeneLenDic, VGPGeneLenDic, PreGeneDic, VGPGeneDic)
    BinSaving(FGDDic, AT + 'FGD_FGDDic')
    Main4 = Main4_FGDTyper(FGDDic, FGGCutoff)
    FGD, ProductSimilarity, FCGProductPair = Main4[0], Main4[1], Main4[2]
    BinSaving(FGD, AT + 'FGD_FGD')
    BinSaving(ProductSimilarity, AT + 'FGD_ProductSimilarity')
    BinSaving(FCGProductPair, AT + 'FGD_FCGProductPair')

    FGD = BinLoading(AT + 'FGD_FGD')
    FGDTypeDic = Main5_FGDResulting(FGD, AT)

    CDSDic = BinLoading(AT + 'FGD_PreCDSDic')
    Annt = BinLoading(AT + 'AnntDat')
    CDSOVLP = Extra_CDSOverlap(CDSDic, Annt)
    GCNOL, WholeCDSFDList, WholeCDSFDRatioGene = CDSOVLP[0], CDSOVLP[1], CDSOVLP[2]

    List2File(WholeCDSFDList, AT + 'FGD_WholeCDSOverlap')
    Dic2File(WholeCDSFDRatioGene, AT + 'FGD_WholeCDSOverlapGeneRatio')
    
    Extra_TEannotation('./02/Log/', 'Pre02All_merged', '02PreLTR_Final', '02PreLTROVLP', 30)
    Extra_TEannotation('./09_2/Log/', 'Pre09All_merged', '09PreSINE_final', '09PreSINEOVLP', 30)
    
    Extra_SegmentalDuplication('./09_2/Log/', 'Pre09All_merged', '09_2SegmentalDuplication', 30) #2465 deleted by mtChr
    Extra_V1RFGG('./09_2/Log/', 'PreFGD_FGD.Pybin', '../V1RgeneList.txt', 'PreFGD_PreSymDic.Pybin', 'PreAnntDat.Pybin', 'PreFGD_PreGeneLenDic.Pybin')
    '''

    Extra_TEannotation('./04/Log/', 'Pre04All_merged_Chr', '04PreTEs', '04PreTEsOVLP', 50)
    Extra_TEannotation('./02/Log/', 'Pre02All_merged', '02PreTEs', '02PreTEsOVLP', 50)
    Extra_TEannotation('./09_2/Log/', 'Pre09All_merged', '09PreTEs', '09PreTEsOVLP', 50)

