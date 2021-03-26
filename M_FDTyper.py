# Classifying by K-mer

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
     
    HapChr = ['NC_009114','NC_009115','NC_009116','NC_009118'] + ['CM014221', 'CM014222', 'CM014223', 'CM014224', 'CM014225', 'RZJT01000027', 'RZJT01000028', 'RZJT01000029', 'RZJT01000030', 'RZJT01000031'] 
    ''' 
    Path, KPath, KmerDic = '/home/ko/cactus/02/', '/disk3/ko_d3/', '02PreKmerDic/'
    SECutoff =  5
    SEMCutoff = 5
    LogPath = Path + 'Log/'
    HighCutoff = 59 * 1.5 #  by cn-hist and fd.sh in merqury
    EHighCutoff = 59 * 3
    Cores  = 30
    FDFileName = 'PreCactusXPurgeDup_Filtered.txt' 
    BamPath = '/disk1/ko_d1/merged_L001_to_L008.sorted.bam'
    D_SECutoff, K_HECutoff = SECutoff, 0
    AT = 'Pre'
    '''
    ''' 
    Path, KPath, KmerDic = '/home/ko/cactus/02/', '/disk3/ko_d3/', '02VGPKmerDic/'
    SECutoff =  2
    SEMCutoff = 2
    LogPath = Path + 'Log/'
    HighCutoff = 59 * 1.5 #  by cn-hist and fd.sh in merqury
    EHighCutoff = 59* 3
    Cores  = 30
    FDFileName = 'VGPCactusXPurgeDup_Filtered.txt' 
    BamPath = '/disk1/ko_d1/merged_L001_to_L008.toVGP.sorted.bam'
    D_SECutoff, K_HECutoff = SECutoff, 0
    AT = 'VGP' 
    '''
    ''' 
    Path, KPath, KmerDic = '/home/ko/cactus/04/', '/disk3/ko_d3/', '04PreKmerDic/'
    SECutoff =  8
    SEMCutoff = 2
    LogPath = Path + 'Log/'
    HighCutoff = 32 * 1.5 #  by cn-hist and fd.sh in merqury
    EHighCutoff = 32 * 3
    Cores  = 30
    FDFileName = 'PreCactusXPurgeDup_Filtered.txt' 
    BamPath = '/disk3/ko_d3/0410x/04_10x_pre/Calann10x2Pre.bam'
    D_SECutoff, K_HECutoff = SECutoff, 0
    AT = 'Pre'
    '''
    '''  
    Path, KPath, KmerDic = '/home/ko/cactus/04/', '/disk3/ko_d3/', '04VGPKmerDic/'
    SECutoff =  2
    SEMCutoff = 2
    LogPath = Path + 'Log/'
    HighCutoff = 32 * 1.5 #  by cn-hist and fd.sh in merqury
    EHighCutoff = 32* 3
    Cores  = 30
    FDFileName = 'VGPCactusXPurgeDup_Filtered.txt' 
    BamPath = '/disk3/ko_d3/0410x/04_10x_VGP/Calann10x2VGP.bam'
    D_SECutoff, K_HECutoff = SECutoff, 0
    AT = 'VGP' 
    '''
    '''     
    Path, KPath, KmerDic = '/home/ko/cactus/09_2/', '/disk3/ko_d3/', '09PreKmerDic/'
    SECutoff =  22
    SEMCutoff = 17
    LogPath = Path + 'Log/'
    HighCutoff = 107 * 1.5 #  by cn-hist and fd.sh in merqury
    EHighCutoff = 107 * 3
    Cores  = 30
    FDFileName = 'PreCactusXPurgeDup_Filtered.txt' 
    BamPath = '/disk3/ko_d3/0910x/09_10x_pre/Ornana10X2Old.bam'
    D_SECutoff, K_HECutoff = SECutoff, 0
    AT = 'Pre'
    HapChrCutoff = 65 * 0.75
    '''
    
    Path, KPath, KmerDic = '/home/ko/cactus/09_2/', '/disk3/ko_d3/', '09VGPKmerDic/'
    SECutoff =  9
    SEMCutoff = 9
    LogPath = Path + 'Log/'
    HighCutoff = 108 * 1.5 #  by cn-hist and fd.sh in merqury
    EHighCutoff = 108 * 3
    Cores  = 10
    FDFileName = 'VGPCactusXPurgeDup_Filtered.txt' 
    BamPath = '/disk3/ko_d3/0910x/09_10x_VGP/Ornana10X2VGP.bam'
    D_SECutoff, K_HECutoff = SECutoff, 0
    AT = 'VGP' 
    HapChrCutoff = 62 * 0.75
    

### Operator ###
def BinSaving(PythonData, FileName):#c     
    with open(LogPath + FileName+'.Pybin', 'wb') as fb: pickle.dump(PythonData, fb)

def BinLoading(FileName):#c      
    with open(LogPath + FileName+'.Pybin', 'rb') as rb: return pickle.load(rb)


### Tools ###
def Pos2Str(Chr, Start, End):#c
    Loci = Chr + ':' + str(Start) + '-' + str(End)
    return Loci

def Str2Pos(Loci):#c
    Chr, Start, End = Loci.split(':')[0], Loci.split(':')[1].split('-')[0], Loci.split('-')[1]
    return [Chr, Start, End]
                                                    

### Main ###
def Main1_FDLoad_Sort(FileName):
    with open(LogPath + FileName, 'r') as f0:
        Union = f0.read().split('\n')
        if Union.count('') > 0: Union.remove('')
        RegionDic = {}  # {Chr:[[Start,End,DataNo],[Start,End,DataNo]]}
        for i in Union:
            DT = i.split('\t')
            if len(DT) > 5:
                BlockNo, FChr, FStart, FEnd, TStart, TEnd, TChr = DT[0], DT[1], DT[2], DT[3], DT[4], DT[5], DT[6] # Data to List
            else: # Single region parsing
                FChr, FStart, FEnd, BlockNo = DT[0], DT[1], DT[2], DT[3]
            if FChr in RegionDic: RegionDic[FChr].append((int(FStart), int(FEnd)))#, Tool + str(DataNo), 'F'])
            else: RegionDic[FChr] = [(int(FStart), int(FEnd))] #Tool + str(DataNo), 'F']] 
        RegionDic2 = RegionDic
        RegionDic = {}
        for i in RegionDic2:
            Chr, FLoci =  i, RegionDic2[i]
            DelDouble = list(set(FLoci))
            DelDouble.sort()
            RegionDic[Chr] = DelDouble
        return RegionDic, Union


def Main2_MultiProc(MainDic, KDicPath, Cores):
    InMainDic = MainDic
    ChrKeys = list(InMainDic.keys())
    random.shuffle(ChrKeys)
    Calling_Keys = ChrKeys
    Calling_Keys = np.array_split(ChrKeys, Cores)
    Calling_Keys = [i.tolist() for i in Calling_Keys]
    Result = parmap.map(Main2_KmerAllocating, Calling_Keys, InMainDic, KDicPath, pm_pbar = True, pm_processes = Cores )
    MKRDic = {}
    for i in Result: MKRDic = {**MKRDic, **i}
    return MKRDic


def Main2_KmerAllocating(KeyList, MainDic, KDicPath):
    InRegionDic = MainDic
    KRDic = {} # {LociStr: {Len:, CN1:, CN2:, Shared:, CN2<:, Zero:, SE:, High:, EHigh:, CN1NSE:}
    KeyLen = len(KeyList)
    CompleteList = []
    for i in KeyList:
        KeyLen -= 1
        Chr, Pos = i, InRegionDic[i]
        if Chr[:2] == 'NC' or Chr[:2] == 'CM': CompleteList.append(Chr)
        print(str(KeyLen) + '\t' + str(CompleteList))
        if HapChr.count(Chr) == 1: continue
        with open(KDicPath + Chr + '.dic', 'rb') as fd:
            InnerLen = len(Pos)
            KmerDic = pickle.load(fd)
            KmerDic = KmerDic[Chr]
            KResultDic = {}
            for j in Pos:
                InnerLen -= 1
                print(str(KeyLen) + '\t' + str(InnerLen))
                #if InnerLen % 100 == 0: print(str(KeyLen) + '\t' + str(InnerLen))
                Start, End = j[0], j[1]
                Len = End - Start + 1
                KmerLen = End - Start + 1 - 19
                Index = Start
                Region = {'KmerLen':KmerLen, 'Len':Len, 'CN1':0, 'CN2':0, 'Shared':0,  'CN2<':0, 'FSRS':[], 'Zero':0, 'SE':0, 'High':0, 'EHigh':0, 'CN1NSE':0} #CN1, [FS+RS], zero, high
                KmerMDepth, CN1MDepth, CN2MDepth, CN3MDepth, CN1NSEMDepth = [], [], [], [], []
                RegionInf = Chr + ':' + str(Start) + '-' + str(End)
                while Index -1 < End - 19:
                    if Index - 1 in KmerDic:
                        KDat = KmerDic[Index - 1]
                        CN, MDepth, Pal = int(KDat[0]), int(KDat[1]), KDat[2]
                        if CN == 1:
                            Region['CN1'] += 1
                            CN1MDepth.append(MDepth)
                            if MDepth > SEMCutoff: # 20200901
                                Region['CN1NSE'] += 1
                                CN1NSEMDepth.append(MDepth)
                        elif CN == 2: Region['CN2'] += 1
                        elif CN >2:
                            Region['CN2<'] += 1
                            CN3MDepth.append(MDepth)

                        if MDepth <= SEMCutoff: Region['SE'] += 1
                        elif MDepth > HighCutoff: Region['High'] += 1
                        if MDepth > EHighCutoff: Region['EHigh'] += 1
                        Index += 1
                        KmerMDepth.append(MDepth)
                    else:
                        Index += 1
                        Region['SE'] += 1
                        Region['Zero'] += 1
                
                Region['MeanDepth'] = sum(KmerMDepth) / KmerLen
                if len(CN1MDepth) > 0: Region['CN1Depth'] = sum(CN1MDepth) / len(CN1MDepth)
                else: Region['CN1Depth'] = 0
                if len(CN2MDepth) > 0: Region['CN2Depth'] = sum(CN2MDepth) / len(CN2MDepth)
                else: Region['CN2Depth'] = 0
                if len(CN3MDepth) > 0: Region['CN3Depth'] = sum(CN3MDepth) / len(CN3MDepth)
                else: Region['CN3Depth'] = 0
                if len(CN1NSEMDepth) > 0: Region['CN1NSEDepth'], Region['CN1NSEDepthSTDEV'] = sum(CN1NSEMDepth) / len(CN1NSEMDepth), np.std(CN1NSEMDepth)
                else: Region['CN1NSEDepth'], Region['CN1NSEDepthSTDEV'] = 0, 0
                KRDic[RegionInf] = Region
                #print(str(KeyLen) + '\t' + str(Start) + '\t' + str(End) + '\t' + str(Region['CN1']) + '\t' + str(Region['CN2']) + '\t' + str(Region['Zero']) + '\t' + str(Region['SE']) + '\t' + str(Region['High']) + '\t' + str(Region['EHigh']) + '\t' + str(Region['MeanDepth']))
    
    return KRDic


def M3Pre_FDKmerDistribution(MainDic):
    InMainDic = MainDic
    CN1NSERatio = []
    CN2Ratio = []
    CN1NSEvsCN2 = []
    CN1SERatio = []
    High = []
    EHigh = []
    for i in InMainDic: # {'Len':Len, 'CN1':0, 'CN2':0, 'Shared':0,  'CN2<':0, 'FSRS':[], 'Zero':0, 'SE':0, 'High':0, 'EHigh':0, 'CN1NSE':0, 'MeanDepth':, 'CN1Depth'CN2Depth,CN3Depth,Cn1NSEDepth}
        DD = InMainDic[i]
        CN1NSERatio.append(DD['CN1NSE'] / DD['KmerLen'])
        CN2Ratio.append(DD['CN2'] / DD['KmerLen'])
        if DD['CN2'] != 0: CN1NSEvsCN2.append(DD['CN1NSE'] / DD['CN2'])
        CN1SERatio.append((DD['CN1'] - DD['CN1NSE']) / DD['KmerLen'])
        High.append(DD['High'] / DD['KmerLen'])
        EHigh.append(DD['EHigh'] / DD['KmerLen'])
    return CN1NSERatio, CN2Ratio, CN1NSEvsCN2, CN1SERatio, High, EHigh


def M3Pre_Hist(List, RStart, REnd, Bins):
    bins = np.arange(RStart, REnd, Bins)
    hist, bins = np.histogram(List, bins)
    Index = 1
    HistDic = {}
    for i in list(hist):
        Freq, Bin = i, list(bins)[Index]
        HistDic[Bin] = Freq
        Index += 1
    return HistDic


def M3_Multi_DepthCalling(DataList, BamPath, TargetTool, Cores):
    InDataList = []
    for i in DataList:
        DT = i.split('\t')
        Tool = DT[0].split(' ')[0]
        if Tool == TargetTool: InDataList.append(i)
    random.shuffle(InDataList)
    InDataList = np.array_split(InDataList, Cores)
    Result = parmap.map(M3_DepthMean, InDataList, BamPath, pm_pbar = True)
    MResultDic = {} # {Loci1:Depth1, ..}
    for i in Result: MResultDic = {**MResultDic, **i}
    return MResultDic


def M3_DepthMean(DataList, BamPath): # Methods; Depth calculation; factors for assembly
    Step, Len, ResultDic = 0, len(DataList), {}
    for i in DataList:
        if (Len - Step) % 10 == 0: print(str(Len - Step))
        Step += 1
        DT = i.split('\t')
        BlockNo, FChr, FStart, FEnd = DT[0], DT[1], DT[2], DT[3]
        FLoci = Pos2Str(FChr, FStart, FEnd)
        SeqLen = int(FEnd) - int(FStart) + 1
        MeanDepth = sp.getoutput("samtools depth -r {0} {1} | awk 'BEGIN {{Total = 0}} {{Total += $3}} END {{print Total / {2}}}'".format(FLoci, BamPath, SeqLen))
        ResultDic[FLoci] = round(float(MeanDepth), 2)
    return ResultDic


def Main3_FDTypeClassifier(MainDic, D_SECutoff, K_HECutoff, DepthDic): # Kmer
    InMainDic = MainDic
    for i in InMainDic:
        Loci, DD = i, InMainDic[i]
        CN1NSERatio = DD['CN1NSE'] / DD['KmerLen']
        CN1SERatio = (DD['CN1'] - DD['CN1NSE']) / DD['KmerLen']
        MeanDepth = float(DepthDic[Loci])
        Chr = Str2Pos(Loci)[0]
        if MeanDepth > D_SECutoff:
            if CN1NSERatio > K_HECutoff: Type = 'HE'
            else: Type = "HO" 
        else: Type = "SE"
        InMainDic[Loci]['Type'] = Type
    return InMainDic


def Main4_TypeAdding(BlockList, FDKmerDic, DepthDic):
    FDBlockList = BlockList
    FDBlockDic = {}
    with open(LogPath + FDFileName.split('.')[0] + '_Type.txt', 'w') as fw:
        for i in FDBlockList:
            DT = i.split('\t')
            BlockNo, FChr, FStart, FEnd = DT[0], DT[1], DT[2], DT[3]
            Loci = Pos2Str(FChr, FStart, FEnd)
            if HapChr.count(FChr) == 1:
                if float(DepthDic[Loci]) > HapChrCutoff: continue # Depth Filtering for PurgeDup
                else: 
                    if HapChrCutoff > D_SECutoff: FDType = 'HO'
                    else: FDType = 'SE'
            else: FDType = FDKmerDic[Loci]['Type']
            Kmer, KmerLen = FDKmerDic[Loci], FDKmerDic[Loci]['KmerLen']
            CN1NSERatio, SERatio, CN2Ratio, EHighRatio = round(Kmer['CN1NSE'] / KmerLen, 3) , round((Kmer['CN1'] - Kmer['CN1NSE']) / KmerLen, 3), round(Kmer['CN2'] / KmerLen, 3), round(Kmer['EHigh'] / KmerLen, 3)
            NewLine = DT + [FDType, str(CN1NSERatio), str(SERatio), str(CN2Ratio), str(EHighRatio)]
            fw.write('\t'.join(NewLine) + '\n')
            FDBlockDic[BlockNo] = NewLine
        return FDBlockDic



#### Main Program Operation #####

Main1 = Main1_FDLoad_Sort(FDFileName)
MainDic, BlockList = Main1[0], Main1[1]
BinSaving(Main1, AT+'FD_TyperMain1')
print('Main1Done')
MKRDic = Main2_MultiProc(MainDic, KPath + KmerDic, Cores)
BinSaving(MKRDic, AT+'MKRDic')
'''
BlockList = BinLoading(AT+'FD_TyperMain1')[1]
New_DepthDic = BinLoading('New_DepthDic')
DepthDic = M3_Multi_DepthCalling(BlockList, BamPath, 'PurgeDup', Cores)
New_DepthDic = {**New_DepthDic, **DepthDic}
BinSaving(New_DepthDic, AT + 'FD_TyperDepthDic')
'''
MKRDic = BinLoading(AT + 'MKRDic')
New_DepthDic = BinLoading(AT + 'FD_TyperDepthDic')
KmerRatioList = M3Pre_FDKmerDistribution(MKRDic) # return CN1NSERatio, CN2Ratio, CN1NSEvsCN2, CN1SERatio
CN1NSERatio, CN2Ratio, CN1NSEvsCN2, CN1SERatio, HighRatio, EHighRatio = KmerRatioList[0], KmerRatioList[1], KmerRatioList[2], KmerRatioList[3], KmerRatioList[4], KmerRatioList[5]
CN1NSEHist = M3Pre_Hist(CN1NSERatio, -0.1, 1.1, 0.01)
CN2Hist = M3Pre_Hist(CN2Ratio, -0.1, 1.1, 0.01)
CN1NSEvsCn2Hist = M3Pre_Hist(CN1NSEvsCN2, -0.1, 1.1, 0.01)
CN1SEHist = M3Pre_Hist(CN1SERatio, -0.1, 1.1, 0.01)
High = M3Pre_Hist(HighRatio, -0.1, 1.1, 0.01)
EHigh = M3Pre_Hist(EHighRatio, -0.1, 1.1, 0.01)


BinSaving(CN1NSEHist, AT+'CN1NSEHist')
BinSaving(CN2Hist, AT+'CN2Hist')
BinSaving(CN1SEHist, AT+'CN1SEHist')
BinSaving(High, AT+'HighHist')
BinSaving(EHigh, AT+'EHighHist')

BlockList = BinLoading(AT+'FD_TyperMain1')[1]
New_DepthDic = BinLoading(AT+'FD_TyperDepthDic')
MKRDic = BinLoading(AT + 'MKRDic')
MKRDic_3 = Main3_FDTypeClassifier(MKRDic, D_SECutoff, K_HECutoff, New_DepthDic)
BlockDic = Main4_TypeAdding(BlockList, MKRDic_3, New_DepthDic) 
BinSaving(BlockDic, AT+"FD_TyperBlockDic")


#### Single Region KRDic Operation #### 4 SE cal
'''
Cores = 50
LogPath, KPath, KmerDic, SinFileName = '/home/ko/cactus/02/Log/', '/disk3/ko_d3/', '02PreKmerDic/', '../../../Bam/02Single_Pre.txt'
SECutoff, HighCutoff, EHighCutoff = 5, 59 * 1.5, 59 * 3
MainDic = Main1_FDLoad_Sort(SinFileName)[0]
SinMKRDic = Main2_MultiProc(MainDic, KPath + KmerDic, Cores)
BinSaving(SinMKRDic, '02PreSinMKRDic')

LogPath, KPath, KmerDic, SinFileName = '/home/ko/cactus/04/Log/', '/disk3/ko_d3/', '04PreKmerDic/', '../../../Bam/04Single_Pre.txt'
SECutoff, HighCutoff, EHighCutoff = 8, 32 * 1.5, 32 * 3
MainDic = Main1_FDLoad_Sort(SinFileName)[0]
SinMKRDic = Main2_MultiProc(MainDic, KPath + KmerDic, Cores)
BinSaving(SinMKRDic, '04PreSinMKRDic')


LogPath, KPath, KmerDic, SinFileName = '/home/ko/cactus/09_2/Log/', '/disk3/ko_d3/', '09PreKmerDic/', '../../../Bam/09Single_Pre.txt'
SECutoff, HighCutoff, EHighCutoff = 22, 107 * 1.5, 107 * 3
MainDic = Main1_FDLoad_Sort(SinFileName)[0]
SinMKRDic = Main2_MultiProc(MainDic, KPath + KmerDic, Cores)
BinSaving(SinMKRDic, '09PreSinMKRDic')


LogPath, KPath, KmerDic, SinFileName = '/home/ko/cactus/02/Log/', '/disk3/ko_d3/', '02VGPKmerDic/', '../../../Bam/02Single.txt'
SECutoff, HighCutoff, EHighCutoff = 2, 59 * 1.5, 59 * 3
MainDic = Main1_FDLoad_Sort(SinFileName)[0]
SinMKRDic = Main2_MultiProc(MainDic, KPath + KmerDic, Cores)
BinSaving(SinMKRDic, '02VGPSinMKRDic')


LogPath, KPath, KmerDic, SinFileName = '/home/ko/cactus/04/Log/', '/disk3/ko_d3/', '04VGPKmerDic/', '../../../Bam/04Single.txt'
SECutoff, HighCutoff, EHighCutoff = 2, 32 * 1.5, 32 * 3
MainDic = Main1_FDLoad_Sort(SinFileName)[0]
SinMKRDic = Main2_MultiProc(MainDic, KPath + KmerDic, Cores)
BinSaving(SinMKRDic, '04VGPSinMKRDic')


LogPath, KPath, KmerDic, SinFileName = '/home/ko/cactus/09_2/Log/', '/disk3/ko_d3/', '09VGPKmerDic/', '../../../Bam/09Single.txt'
SECutoff, HighCutoff, EHighCutoff = 9, 108 * 1.5, 108 * 3 
MainDic = Main1_FDLoad_Sort(SinFileName)[0]
SinMKRDic = Main2_MultiProc(MainDic, KPath + KmerDic, Cores)
BinSaving(SinMKRDic, '09VGPSinMKRDic')
'''

