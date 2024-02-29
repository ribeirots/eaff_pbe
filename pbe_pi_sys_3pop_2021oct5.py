# calculate PBE and pi per window defined by # of SNPs (imported from PoPoolation2 - SNP by SNP analysis and processed with my own python script)
# Author: Tiago Ribeiro (tiaaagosr@gmail.com)

import re, math, statistics, sys

# Argument 1 = windowsBoundaries.csv
# Argument 2 = bre or mie data.sync (sync file with three pops, focal pop should be the third. Adding a 4th population or more will break the script.
# Chrm/Scaffold/Contig need to be numbers only. Also accepts 'Scaffold###' in *.sync file.
# Sync row example1 (\t delimited): 1        186        G        0:0:0:0:0:0        0:0:0:0:0:0        0:0:0:3:0:0 
# Sync row example2 (\t delimited): 1        186        G        0:0:0:0                0:0:0:0                0:0:0:3 
# Both files should be ordered by Chrm/Scaffold/Contig which themselves are ordered by SNP.

filename = str(sys.argv[2])[:-5]

MaxFST = 0.95
FreqThresh = 0.05 # Minimum frequency of the Minor Allele for FST calculation

sitesInTheWindow = 0
currentSite = []
snpDATA = []
windowDATA = []

with open(sys.argv[1]) as inputData: # read the file into an object (a list), each line as a single str
    windBound = inputData.readlines()


for i in range(0, len(windBound)): #editing windBound data
    windBound[i] = re.split(r'\t',windBound[i])
    windBound[i][-1] = windBound[i][-1].rstrip()
    windBound[i] = list(map(int,windBound[i]))

Pi1Sum = 0
Pi2Sum = 0
Pi3Sum = 0
Num12Sum = 0
Den12Sum = 0
Num13Sum = 0
Den13Sum = 0
Num23Sum = 0
Den23Sum = 0
MaxPBS = -999
MaxFST12 = -999
MaxFST13 = -999
MaxFST23 = -999
i = 0
lastSNP = [-999, -999]
snp_skip = 0
scaff_check = -999
    
with open(sys.argv[2]) as inputData:
    for fresh in inputData:
        fresh = re.split(r'\t',fresh[:-1])
        fresh[0] = re.sub('Scaffold','',fresh[0])
        fresh[0] = int(fresh[0])
        fresh[1] = int(fresh[1])
        if fresh[0] != scaff_check:
            print('Scaffold: '+str(fresh[0]))
            scaff_check = fresh[0]

        # Check if window stats should be calculated
        if i == len(windBound):
            print('All windows analyzed. Waiting, reading the final SNPs - SKIP WINDOW.')  
        elif lastSNP[0] != -999:
            if ((fresh[0] != windBound[i][0]) or (fresh[1] > windBound[i][2])): # different scaffold or SNP outside the current window.
                if lastSNP[0] == windBound[i][0] and lastSNP[1] <= windBound[i][2] and lastSNP[1] > windBound[i][1]: # calc window stats when SNP is out of the current window and last snp was part of it. increment window by 1
                    if sitesInTheWindow != 0:
                        Pi1 = Pi1Sum/sitesInTheWindow
                        Pi2 = Pi2Sum/sitesInTheWindow
                        Pi3 = Pi3Sum/sitesInTheWindow
                        if Den12Sum != 0:
                            FST12_window = Num12Sum/Den12Sum
                        else:
                            FST12_window = 0
                        if Den13Sum != 0:
                            FST13_window = Num13Sum/Den13Sum
                        else:
                            FST13_window = 0
                        if Den23Sum != 0:
                            FST23_window = Num23Sum/Den23Sum
                        else:
                            FST23_window = 0
                        winPBS = ((-1 * math.log(1 - FST13_window)) + (-1 * math.log(1 - FST23_window)) - (-1 * math.log(1 - FST12_window))) / 2
                        # windowData columns: Scaffold, Starting Position, Ending Position, Window_Length, pi_MAE, pi_VIE, pi_Fresh, T_MAExVIE, T_MAExFresh, T_VIExFresh, FST_MAExVIE, MaxFST12, Position of MaxFST12,  Max_PBS, position from Max_PBS_SNP, window PBS
                        windowDATA.append([windBound[i][0]] + [windBound[i][1]] + [windBound[i][2]] + [windBound[i][2] - windBound[i][1]] + [Pi1] + [Pi2] + [Pi3] + [(-1 * math.log(1 - FST12_window))] + [(-1 * math.log(1 - FST13_window))] + [(-1 * math.log(1 - FST23_window))] + [FST12_window] + [MaxFST12] + [MaxFST12_info] + [FST13_window] + [MaxFST13] + [MaxFST13_info] + [FST23_window] + [MaxFST23] + [MaxFST23_info] + [MaxPBS] + [MaxPBS_info] + [winPBS])
                    else:
                        print('Window on scaffold '+str(windBound[i][0])+', ranging from '+str(windBound[i][1])+' to '+str(windBound[i][2])+' was not analyzed. Check whether all the SNPs in the regions didnt differ among the populations.')

                    i += 1

                    Pi1Sum = 0
                    Pi2Sum = 0
                    Pi3Sum = 0
                    Num12Sum = 0
                    Den12Sum = 0
                    Num13Sum = 0
                    Den13Sum = 0
                    Num23Sum = 0
                    Den23Sum = 0
                    MaxPBS = -999
                    MaxFST12 = -999
                    MaxFST13 = -999
                    MaxFST23 = -999
                    sitesInTheWindow = 0
                    snp_skip = 0

                elif snp_skip == 1:
                    print('Scaffold skipped. Likely a scaffold without any SNP being skip.')
                    print(windBound[i])
                    print('Last SNP:')
                    print(lastSNP)
                    print('Current SNP')
                    print(fresh)
                    
                else:
                    print('Error. Switching windows without analyzing the previous window. The last SNP should have been in the last window but it was not.')
                    print('Last SNP: ')
                    print(lastSNP)
                    print('Current window: ')
                    print(windBound[i])
                    print('Current SNP')
                    print(fresh)
                    exit()

        # Check if SNP is inside the current window
        if i == len(windBound):
            print('All windows analyzed. Waiting reading the final SNPs - SKIP SNP.')
        elif fresh[0] > windBound[i][0]:
            while fresh[0] != windBound[i][0] or fresh[0] > windBound[i][0]: # If the script breaks out of bound here it might also be an ordering issue.
                i += 1
            if fresh[0] == windBound[i][0]:
                snp_skip = 0
            else:
                snp_skip = 1
        elif fresh[0] < windBound[i][0]:
            snp_skip = 1
        else: # scaff_snp == window
            snp_skip = 0
            


        if snp_skip == 0:
            if i == len(windBound):
                print('All windows analyzed. Waiting reading the final SNPs - SKIP SNP.')                
            else:
                if fresh[1] > windBound[i][2]: 
                    while fresh[1] > windBound[i][2] and fresh[0] == windBound[i][0]:
                        i += 1
                    if fresh[0] != windBound[i][0]: # Error, the SNP was outside all the available windows.
                        print('Error, the SNP '+str(fresh[1])+' on Scaffold '+str(fresh[0])+' was outside all the available windows.')
                        exit()
            if fresh[1] <= windBound[i][2]: 
                nucleotides = []
                currentSite = fresh[:2]+fresh[3:6]
                
                nucleotides.append(list(map(int,re.split(r':',currentSite[2])))) # this stores the info from the current site.
                nucleotides.append(list(map(int,re.split(r':',currentSite[3]))))
                nucleotides.append(list(map(int,re.split(r':',currentSite[4]))))

                nucleotides[0] = nucleotides[0][:4]
                nucleotides[1] = nucleotides[1][:4]
                nucleotides[2] = nucleotides[2][:4]

                SampleSize1 = sum(nucleotides[0])
                SampleSize2 = sum(nucleotides[1])
                SampleSize3 = sum(nucleotides[2])
                
                if ((SampleSize1 != 0) and (SampleSize2 != 0) and (SampleSize3 != 0)):                   
                    #Define major and minor alleles (mask third/fourth alleles)
                    AFreqSum = (nucleotides[0][0]/SampleSize1) + (nucleotides[1][0]/SampleSize2) + (nucleotides[2][0]/SampleSize3)
                    TFreqSum = (nucleotides[0][1]/SampleSize1) + (nucleotides[1][1]/SampleSize2) + (nucleotides[2][1]/SampleSize3)
                    CFreqSum = (nucleotides[0][2]/SampleSize1) + (nucleotides[1][2]/SampleSize2) + (nucleotides[2][2]/SampleSize3)
                    GFreqSum = (nucleotides[0][3]/SampleSize1) + (nucleotides[1][3]/SampleSize2) + (nucleotides[2][3]/SampleSize3)

                    FreqSums = [AFreqSum, TFreqSum, CFreqSum, GFreqSum]
                    MajorAllele = sorted(range(len(FreqSums)), key=lambda k: FreqSums[k])[-1] # gets index of major freq allele
                    MinorAllele = sorted(range(len(FreqSums)), key=lambda k: FreqSums[k])[-2] # gets index of minor freq allele

                    SampleSize1 = nucleotides[0][MajorAllele] + nucleotides[0][MinorAllele]
                    SampleSize2 = nucleotides[1][MajorAllele] + nucleotides[1][MinorAllele]
                    SampleSize3 = nucleotides[2][MajorAllele] + nucleotides[2][MinorAllele]

                if ((SampleSize1 != 0) and (SampleSize2 != 0) and (SampleSize3 != 0)): # I am including this check point  because there was 1 SNP in which VIE didn't have any of the Major or Minor freq allele. SNP info: [85, 332364, '0:11:29:1', '19:0:0:14', '0:24:40:2']
                    # Pi calculations (adding site data to window sum)
                    sitesInTheWindow += 1
                    MajorFreq1 = nucleotides[0][MajorAllele] / SampleSize1
                    MajorFreq2 = nucleotides[1][MajorAllele] / SampleSize2
                    MajorFreq3 = nucleotides[2][MajorAllele] / SampleSize3
                    MinorFreq1 = nucleotides[0][MinorAllele] / SampleSize1
                    MinorFreq2 = nucleotides[1][MinorAllele] / SampleSize2
                    MinorFreq3 = nucleotides[2][MinorAllele] / SampleSize3

                    Pi1Sum += 2 * MajorFreq1 * MinorFreq1
                    Pi2Sum += 2 * MajorFreq2 * MinorFreq2
                    Pi3Sum += 2 * MajorFreq3 * MinorFreq3

                    # If conditions matched to call a SNP, calculate FST
                    if (MinorFreq1 + MinorFreq2 + MinorFreq3 )/3 > FreqThresh:
                        
                        # Reynolds FST calculations (pops 1 & 2)
                        denom = MinorFreq1 + MinorFreq2 - 2 * MinorFreq1 * MinorFreq2
                        SharedNum = SampleSize1 * (2 * MinorFreq1 - 2 * MinorFreq1 ** 2) + SampleSize2 * (2 * MinorFreq2 - 2 * MinorFreq2 ** 2)
                        NumA = (MinorFreq1 - MinorFreq2) ** 2
                        FracNum = (SampleSize1 + SampleSize2) * SharedNum
                        FracDen = 4 * SampleSize1 * SampleSize2 * (SampleSize1 + SampleSize2 - 1)
                        frac = FracNum / FracDen
                        WholeNum = NumA - frac
                        DenFracNum = (4 * SampleSize1 * SampleSize2 - SampleSize1 - SampleSize2) * SharedNum
                        DenFrac = DenFracNum / FracDen
                        WholeDen = NumA + DenFrac
                        Num12Sum += WholeNum
                        Den12Sum += WholeDen
                        if WholeDen != 0:
                            FST12 = WholeNum / WholeDen
                        else:
                            FST12 = 0
                        if FST12 == 1:
                            FST12 = MaxFST
                               
                        # Reynolds FST calculations (pops 1 & 3)
                        denom = MinorFreq1 + MinorFreq3 - 2 * MinorFreq1 * MinorFreq3
                        SharedNum = SampleSize1 * (2 * MinorFreq1 - 2 * MinorFreq1 ** 2) + SampleSize3 * (2 * MinorFreq3 - 2 * MinorFreq3 ** 2)
                        NumA = (MinorFreq1 - MinorFreq3) ** 2
                        FracNum = (SampleSize1 + SampleSize3) * SharedNum
                        FracDen = 4 * SampleSize1 * SampleSize3 * (SampleSize1 + SampleSize3 - 1)
                        frac = FracNum / FracDen
                        WholeNum = NumA - frac
                        DenFracNum = (4 * SampleSize1 * SampleSize3 - SampleSize1 - SampleSize3) * SharedNum
                        DenFrac = DenFracNum / FracDen
                        WholeDen = NumA + DenFrac
                        Num13Sum += WholeNum
                        Den13Sum += WholeDen
                        if WholeDen != 0:
                            FST13 = WholeNum / WholeDen
                        else:
                            FST13 = 0
                        if FST13 == 1:
                            FST13 = MaxFST

                        # Reynolds FST calculations (pops 2 & 3)
                        denom = MinorFreq2 + MinorFreq3 - 2 * MinorFreq2 * MinorFreq3
                        SharedNum = SampleSize2 * (2 * MinorFreq2 - 2 * MinorFreq2 ** 2) + SampleSize3 * (2 * MinorFreq3 - 2 * MinorFreq3 ** 2)
                        NumA = (MinorFreq2 - MinorFreq3) ** 2
                        FracNum = (SampleSize2 + SampleSize3) * SharedNum
                        FracDen = 4 * SampleSize2 * SampleSize3 * (SampleSize2 + SampleSize3 - 1)
                        frac = FracNum / FracDen
                        WholeNum = NumA - frac
                        DenFracNum = (4 * SampleSize2 * SampleSize3 - SampleSize2 - SampleSize3) * SharedNum
                        DenFrac = DenFracNum / FracDen
                        WholeDen = NumA + DenFrac
                        Num23Sum += WholeNum
                        Den23Sum += WholeDen
                        if WholeDen != 0:
                            FST23 = WholeNum / WholeDen
                        else:
                            FST23 = 0
                        if FST23 == 1:
                                FST23 = MaxFST

                        # Calculating PBS
                        SitePBS = ((-1 * math.log(1 - FST13)) + (-1 * math.log(1 - FST23)) - (-1 * math.log(1 - FST12))) / 2
                        if SitePBS > MaxPBS:
                            MaxPBS = SitePBS
                            MaxPBS_info = currentSite[1]
                                 
                        if FST12 > MaxFST12:
                            MaxFST12 = FST12
                            MaxFST12_info = currentSite[1]

                        if FST13 > MaxFST13:
                            MaxFST13 = FST13
                            MaxFST13_info = currentSite[1]

                        if FST23 > MaxFST23:
                            MaxFST23 = FST23
                            MaxFST23_info = currentSite[1]

                        # Check if site pass threshold to be printed out
                        # snpDATA colums: Scaffold, snp position, alleles in MAE, VIE and Fresh, T12, T13, T23, FST12, PBS
                        snpDATA.append(currentSite + [(-1 * math.log(1 - FST12)), (-1 * math.log(1 - FST13)), (-1 * math.log(1 - FST23))] + [FST12] + [FST13] + [FST23] + [SitePBS])
            lastSNP = [fresh[0], fresh[1]]

    print('calc win pbe')
    # Calculating window PBE    
    medianPBS = statistics.median([pbs[-1] for pbs in windowDATA]) # Used for both Window and SNP PBE calculation
    medianT_MAExVIE = statistics.median([pbs[7] for pbs in windowDATA])
    for i in range(0, len(windowDATA)):
        windowDATA[i] = windowDATA[i] + [windowDATA[i][-1] - (windowDATA[i][7] * (medianPBS/medianT_MAExVIE))] # add window PBE to the end of windowDATA

    # Calculating SNP PBE
    print('calc snp pbe')
    for i in range(0, len(snpDATA)):
        snpDATA[i] = snpDATA[i] + [snpDATA[i][-1] - (snpDATA[i][5] * (medianPBS/medianT_MAExVIE))] # add snp PBE to the end of snpDATA




snpPrint = open(filename + 'data_persnp_PBE.txt', 'w')
header_s = 'Scaffold\tsnp position\talleles_1\talleles_2\talleles_3\tT1.2\tT1.3\tT2.3\tFST1.2\tFST1.3\tFST2.3\tPBS\tPBE\n'
snpPrint.write(header_s)
for obs in snpDATA:
  snpPrint.write('\t'.join(str(summStats) for summStats in obs) + '\n')
snpPrint.close()

print('calc win max pbe')
# Calculating window MaxPBE
j = 0
for i in range(0,len(windowDATA)):
    MaxPBE = -999
    while ((j < len(snpDATA)) and (snpDATA[j][0] == windowDATA[i][0]) and (snpDATA[j][1] <= windowDATA[i][2])):
        if snpDATA[j][-1] > MaxPBE:
            MaxPBE = snpDATA[j][-1]
            MaxPBE_info = snpDATA[j][1]
        j += 1
    if ((snpDATA[j-1][0] == windowDATA[i][0]) and (snpDATA[j-1][1] <= windowDATA[i][2]) and (snpDATA[j-1][1] > windowDATA[i][1])):
        windowDATA[i] = windowDATA[i] + [MaxPBE] + [MaxPBE_info]
    else:
        windowDATA[i] = windowDATA[i] + [-999] + [-999]

for wdata in windowDATA:
    if wdata[11] == -999: # FST 12
        wdata[11] = 'NA'
        wdata[10] = 'NA'
    if wdata[14] == -999: # FST 13
        wdata[14] = 'NA'
        wdata[13] = 'NA'
    if wdata[17] == -999: # FST 23
        wdata[17] = 'NA'
        wdata[16] = 'NA'
    if wdata[19] == -999: # PBS
        wdata[19] = 'NA'
        wdata[21] = 'NA'
    if wdata[23] == -999: # PBE
        wdata[23] = 'NA'
        wdata[22] = 'NA'


windowPrint = open(filename +'_window_summary_PBE', 'w')
# Scaff, Start, End, Length, Pi1, Pi2, Pi3, t12, t13, t23, win_fst_12, max_fst_12, pos_max_fst_12, maxpbs, maxpbs_pos, win_pbs
header_w = 'Scaffold\tStart\tEnd\tWindowLength\tpi_1\tpi_2\tpi_3\tT_1.2\tT_1.3\tT_2.3\tWinFST_1.2\tMaxSNP_FST_1.2\tPositionMaxSNP_FST_1.2\tWinFST_1.3\tMaxSNP_FST_1.3\tPositionMaxSNP_FST_1.3\tWinFST_2.3\tMaxSNP_FST_2.3\tPosition_MaxSNP_FST_2.3\tMaxSNP_PBS\tPositionMaxSNP_PBS\tWindowPBS\tWindowPBE\tMaxSNP_PBE\tPosition_MaxSNP_PBE\n'
windowPrint.write(header_w)
for obs in windowDATA:
  windowPrint.write('\t'.join(str(summStats) for summStats in obs) + '\n')
windowPrint.close()
