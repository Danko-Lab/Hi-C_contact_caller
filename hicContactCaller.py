import os
import timing

import sys
import statsmodels.stats.multitest as mt
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector as ivect
import rpy2.robjects as robjects
import rpy2

mass = importr('MASS')
stats = importr('stats')
base = importr('base')

# Set testing flags
testMode = False
testingNumberOfRows = 10000000

peakFile = sys.argv[1] # File of "bait" peaks to analyze (format: Chromosome <\t> peak-center-position)
tssFile = sys.argv[2] # File of "prey" TSSs to analyze (format: Chromosome <\t> TSS-center-position)
conFile = sys.argv[3] # Contact file in juicer format 
outputFile = sys.argv[4] # Output file
chromToProcess = sys.argv[5] # Chromosome to process

# Initialize
output_folder  = ''
outputFile     = output_folder + outputFile + '_' + chromToProcess + '.txt'

out = open(outputFile, "w")

# Display arguments
print ''
print 'peakFile is:', peakFile
print 'tssFile is:', tssFile
print 'conFile is:', conFile
print 'outputFile is:', outputFile
print 'chromToProcess is:', chromToProcess
print 'Test mode is:', str(testMode)
print ''

DIST = 300000 # Search window length
CAP = 2500 # Capture window around locus

#f = robjects.r['expFunction'] # apparently not used

# Get contact data from file

loci, contacts, peaks, TSSs, conList, contactProbabilities, pvalues = [],[],[],[],[],[],[]

print "Getting loci for contact analysis"

rowNumber = 0
for locus in open(peakFile): # locus line format: <chr>\t<start of peak>\t<end of peak>

    if testMode:
        rowNumber += 1
        if rowNumber >= testingNumberOfRows: break
        
    locus = locus.strip()
    locus = locus.split()
    midPoint = int(locus[1]) + ((int(locus[2]) - int(locus[1])) / 2)
    if midPoint >= DIST:
        start = midPoint - DIST
    else:
        start = 0
    stop = midPoint + DIST
    lociChromo = locus[0]
    if lociChromo[:3] == 'chr': lociChromo = lociChromo[3:]

    loci.append((lociChromo,start,stop,midPoint))
    contacts.append([])
    TSSs.append([])
    
# Get TSSs from file

print "Sorting TSSs"

rowNumber = 0
for tss in open(tssFile): # TSS line format: <chr>\t<start of TSS peak>\t<end of TSS peak>\t<gene name>

    if testMode:
        rowNumber += 1
        if rowNumber >= testingNumberOfRows: break
        
    split = tss.strip()
    split = tss.split()
    chrom, pos, geneName = split[0], int(split[1]), split[3]
    if chrom[:3] == 'chr': chrom = chrom[3:]

    if chromToProcess <> 'all' and chromToProcess <> chrom: continue

    for i in range(len(loci)):
        if chrom == loci[i][0]:
            if loci[i][1] < pos < loci[i][2]:
                distance = pos - loci[i][3]
                TSSs[i].append((distance,chrom,pos,geneName))
                print 'TSS:', distance, chrom, pos, geneName
            else:
                continue
        else:
            continue
            

print "Scanning",len(loci),"loci for interactions!"

rowNumber = 0
for line in open(conFile):

    rowNumber += 1
    if testMode:
        if rowNumber >= testingNumberOfRows: break
        
    if rowNumber%100000 == 0: 
        sys.stdout.write('.')
        sys.stdout.flush()

    split = line.split()

    chromo1 = split[1]
    if chromo1[:3] == 'chr': chromo1 = chromo1[3:]
    chromo2 = split[5]
    if chromo2[:3] == 'chr': chromo2 = chromo2[3:]
    if chromo1 == chromo2:
        if int(split[2]) <= int(split[6]):
            conList.append([chromo1,split[2],chromo2,split[6]])
        else:
            conList.append([chromo2,split[6],chromo1,split[2]])
    else:
        continue

# only keep contacts around loci

for entry in conList: 
    
    chr1,pos1,chr2,pos2 = entry[0],int(entry[1]),entry[2],int(entry[3])
    for i in range(len(loci)):
        if chr1 == loci[i][0]:      
            if int(loci[i][1]) <= pos1 <= int(loci[i][2]):
                contacts[i].append([chr1,pos1,pos2])
            else:
                continue
        else:
            continue
            
    for i in range(len(loci)):
        if chr1 == loci[i][0]:
            if int(loci[i][1]) <= pos2 <= int(loci[i][2]):
                contacts[i].append([chr1,pos1,pos2])
            else:
                continue 
        else:
            continue

print len(contacts[0])

# Calculate observed and expected distributions for each locus 


counter = 0
for locus in open(peakFile):
    expect, plus, minus = [],[],[]
    locus = locus.strip().split()
    chrom, position = locus[0], int(locus[1])
    if chrom[:3] == 'chr': chrom = chrom[3:]
    
    # Calculate parameters for expected distribution
    for contact in contacts[counter]:
        distance = int(contact[2]) - int(contact[1])
	expect.append(distance)


    if len(expect) < 1:
        expect.append(1)

    empDist = ivect(expect)
    robjects.r.assign("empDist",empDist)
    
    P = robjects.r("ecdf(empDist)")
    robjects.r.assign("P",P)

    
    # Get empirical distribution for locus
    
    baitStart, baitStop = (position-CAP), (position+CAP)
    
    for contact in contacts[counter]: 
        
        if baitStart <= int(contact[1]) <= baitStop:
            distance = contact[2]-contact[1]
            plus.append(distance)
        else:
            continue

    for contact in contacts[counter]: 
            
        if baitStart <= int(contact[2]) <= baitStop:
            distance = contact[2]-contact[1]
            minus.append(distance)
        else:
            continue
            
    #//---print len(plus), len(minus)
    
        # test plus strand contacts 
    
    for tss in TSSs[counter]:
        if tss[0] > 0:
            print "Testing plus!"
            tssChrom, tssPosition, tssDist, tssGene = tss[1],tss[2], tss[0], tss[3]
            preyStart, preyStop = (tssDist - CAP), (tssDist + CAP)
#            print tssPosition, preyStart, preyStop
            
            
            if (tssDist - (CAP*1.5)) > 0:
                expStart, expStop = (tssDist - (CAP*1.5)), (tssDist + (CAP*1.5))
            else:
                expStart, expStop = 0, (tssDist + (CAP*1.5))
            
            robjects.r.assign("expStart",expStart)
            robjects.r.assign("expStop", expStop)
            z = robjects.r("seq(expStart,expStop,by=1)")
            robjects.r.assign("z",z)
            p = robjects.r("P(z)")
            exp = p[-1]-p[0]
        
            observed = 0
            for distance in plus:
                if preyStart <= distance <= preyStop:
                    observed += 1
                else:
                    continue            

            
            expected = exp*len(plus)
#            print "Expected:",str(exp),str(expected)

                
#            print observed, expected
            v = robjects.FloatVector([observed, expected, (len(plus)-observed), (len(plus)-expected)])
            robjects.r.assign("v", v)
            test = robjects.r("matrix(v, nrow = 2)")
            robjects.r.assign("matrix", test)
 #           print robjects.r("matrix")
            try:
                robjects.r("exact <- fisher.test(matrix, alternative = 'greater')")
                result = robjects.r("exact")
                p_val = result[0][0]
            except:
                p_val = 999.99
#                print "Exception!"
#            print p_val
            contactProbabilities.append([chrom, position, tssChrom, tssPosition, tssDist, p_val, observed, expected, tssGene])
            pvalues.append(p_val)
        
    # test minus strand contacts
    
    for tss in TSSs[counter]:
        if tss[0] <= 0:
            print "Testing minus!"
            tssChrom, tssPosition, tssDist, tssGene = tss[1],tss[2], tss[0], tss[3]
            preyStart, preyStop = (tssDist - CAP), (tssDist + CAP)
#            print tssPosition, preyStart, preyStop

            observed = 0
            for distance in minus:
                if preyStart <= (-1*distance) <= preyStop:
                    observed += 1
                else:
                    continue
            
            preyStart, preyStop = abs(preyStop), abs(preyStart)            
            
            robjects.r.assign("preyStart",preyStart)
            robjects.r.assign("preyStop", preyStop)
            z = robjects.r("seq(preyStart,preyStop,by=1)")
            robjects.r.assign("z",z)
            p = robjects.r("P(z)")
            exp = p[-1]-p[0]
            
            expected = exp*len(plus)
#           print "Expected:",str(exp),str(expected)

#            print observed, expected
            v = robjects.FloatVector([observed, expected, (len(minus)-observed), (len(minus)-expected)])
            robjects.r.assign("v", v)
            test = robjects.r("matrix(v, nrow = 2)")
            robjects.r.assign("matrix", test)
#            print robjects.r("matrix")
            try:
                robjects.r("exact <- fisher.test(matrix, alternative = 'greater')")
                result = robjects.r("exact")
                p_val = result[0][0]
            except:
                p_val = 999.99
#                print "Exception!"

#            print p_val
            contactProbabilities.append([chrom, position, tssChrom, tssPosition, tssDist, p_val, observed, expected, tssGene])
            pvalues.append(p_val)
    
    counter += 1
    if counter%1000 == 0: print 'Counter:', str(counter)



corrected = mt.multipletests(pvalues, method='fdr_bh')

p_count = 0
for prob in contactProbabilities:
    outStr = str(prob[0]) +"\t"+ str(prob[1]) +"\t"+ str(prob[2]) +"\t"+ str(prob[3]) +"\t"+ str(prob[4]) +"\t"+ str(prob[5]) +"\t"+ str(prob[6]) +"\t"+ str(prob[7]) +"\t"+ str(corrected[1][p_count]) +"\t"+ str(prob[8]) +"\n"
    out.write(outStr)
    p_count += 1

#print p_count
