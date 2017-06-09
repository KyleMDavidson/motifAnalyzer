import scipy.stats as sps
import altschulEriksonDinuclShuffle
# import tkFileDialog
# from Tkinter import Tk
from sys import argv
import itertools
import numpy as np
import matplotlib
import random
from sklearn import mixture
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib.text import Text as plttext
from collections import Counter 
from sys import argv
from scipy.stats import ks_2samp


normVals = np.random.normal 
#module has its own alphabet dependencies
dinuclShuffle = altschulEriksonDinuclShuffle.dinuclShuffle

fpmers = '/Users/kyledavidson/Documents/Syllabi/1640/merdata/'

fp1640 = '/Users/kyledavidson/Documents/Syllabi/1640/'
file = 'statFile'
statFile = fp1640 + file
uv = 'U4v2.fasta'
sourceFile = fp1640 + uv
posFile = fp1640 + 'temp'
alphabet = 'u', 'g', 'c', 'a'

# Tk().withdraw()

# Filename = tkFileDialog.askopenfilename()
# inputFile = open(uiFilename, 'r')


# uiScanFile = tkFileDialog.askopenfilename()
# sourceFile = open(uiScanFile, 'r')

#statFile = tkFileDialog.askopenfile()
# statfile = filepath


def catalogue():
    """makes exhaustive list of 4mer strings"""
    catalogue = []
    for i in (alphabet):
        for j in (alphabet):
            for k in (alphabet):
                for h in (alphabet):
                    mer = i + j + k + h
                    catalogue.append(mer)
    return catalogue


def catalog(length):
    catalogue = itertools.product(alphabet, repeat=length)
    catalog = []
    for element in catalogue:
        catalog.append(element)
    return catalog


def poswrite(line, mer, start, outFile):
    """does not provide ordered output - deprecated by poswrite_iter"""
    pos = line.find(mer, start)
    if pos != -1:
        outFile.write(str(pos) + ',')
        poswrite(line, mer, pos + 1, outFile)


def poswrite_iter(line, mer, charCount, outFile):
    """iterates through lines making list of distances to each instance of a given string k-mer(i.e. relative position from start of string)"""
    pos = line.find(mer)
    count = 0
    while pos != -1:
        outFile.write(str(pos + charCount) + ',')
        pos = line.find(mer, pos + 1)
        count += 1
    return count


def pos_calc(mer, inputFile):
    """takes in a FASTA and searches for mers, given as param"""
    outFile = open('temp', 'a+')
    inputFile = open(inputFile, 'r')
    globalMerCount = 0
    for line in inputFile:
        if '>' in line:
            if 'merCount' in locals():
                outFile.write(str(merCount))
                outFile.write('\n')
            csvLine = line.replace(',', '')
            outFile.write(csvLine.strip() + ',')
            charCount = 0
            merCount = 0
        else:
            lineMerCount = poswrite_iter(line, mer, charCount, outFile)
            charCount += len(line)
            merCount += lineMerCount
            globalMerCount += lineMerCount
    outFile.write(str(globalMerCount) + '\n')
    outFile.write('<' + mer + ',' + str(globalMerCount) + '\n')


def merCount(posFile, separator):
    """DEPRECATED by distanceCount. takes positions from temp file and produces distances statistics, ugly"""
    statFile = open('statFile', 'w')
    temp = open(posFile, 'r')
    next(temp)
    for line in temp:
        length = len(line)
        if line[length - i - 1] != separator:
            i += 1
        count = int(line[(length - i):length])
        i = line.find(',')
        j = line.find(',', i + 1)
        statFile.write(line[0:i] + ',')
        lowerMer = int(line[i + 1:j])
        for i in (1, count):
            k = line.find(',', j + 1)
            upperMer = int(line[j + 1:k])
            distance = upperMer - lowerMer
            statFile.write(str(distance) + ',')
            j = k


def distanceCount(posFile, statFile, separator):
    """takes positions from temp file and produces distances statistics"""
    statFile = open(statFile, 'a+')
    temp = open(posFile, 'r')
    for line in temp:
        positionArray = line.split(",")
        statFile.write(positionArray[0] + ',')

        length = len(positionArray)
        if length > 1:
            for distance in range(2, length - 1):
                if length == 2:
                    statFile.write(MAX_LEN + ',')

                else:
                    d = (
                        str(int(positionArray[distance]) - int(positionArray[distance - 1])) + ',')
                    statFile.write(d)
        statFile.write('\n')
    statFile.close()


def recurseDistances(statFile, array):
    statFile = open(statFile, 'r')
    try:
        line = next(statFile)
        data = line.split(",")
        thisArray = map(int, data[1:len(data) - 1])
        qs(thisArray)
        array = recurseDistances(statFile, array)
        bigArray = mergeArrays(thisArray, array)
        return bigArray
    except:
        return array
# pos_calc('ug',input)
# distanceCount(sourceFile, 'statFile.txt', ',')
# def statCalc(sortedFile):


def percentiles(statFile, writeFile, mer):
    """default behavior is to write for mer a list of distances that define  5th, 10th, 25th, 50th, 75th, 90th, and 95th
    percentiles"""
    writeFile = open(writeFile)
    for line in statFile:
        if '<' in line:
            test = line.split(',')
            if mer in test:
                test[1] = count
    array = []
    percentiles = recurseDistances(statFile, array)
    elements = len(percentiles)
    writeFile.write(mer + ',')
    writeFile.write(count + ',')
    writeFile.write(percentiles[elements * .5] + ',')
    writeFile.write(percentiles[elements * .10] + ',')
    writeFile.write(percentiles[elements * .25] + ',')
    writeFile.write(percentiles[elements * .50] + ',')
    writeFile.write(percentiles[elements * .75] + ',')
    writeFile.write(percentiles[elements * .90] + ',')
    writeFile.write(percentiles[elements * .95] + ',')


def getdistance(sequence, mer):
    """. Takes sequence
    as a string, and mer also as a string to be identified within the string."""
    ABSENT = -2
    SINGLETON = -1
    pos1 = sequence.find(mer)
    if pos1 != -1:

        pos2 = sequence.find(mer, pos1 + len(mer))
        if pos2 != -1:
            return (pos2 - pos1, pos1)
        else:
            return (SINGLETON, pos1)
    else:
        return (ABSENT,)


def recursedistances2(sequence, mer):
    """wrapper for getdistances. Sets things up for recursion."""
    dList = []
    pos1 = sequence.find(mer)
    if pos1 == -1:
        return dList
    else:
        getdistances(sequence, mer, dList, pos1)
    return dList


def getdistances(sequence, mer, dList, priorPos):
    """. Takes sequence
    as a string, and mer also as a string to be identified within the string."""
    SINGLETON = 1000000
    ABSENT = 2000000
    pos = sequence.find(mer, priorPos + len(mer))
    if pos == -1:
        return
    else:
        d = pos - priorPos
        dList.append(d)
        getdistances(sequence, mer, dList, pos)


def simulateddistances(sequence, mer):
    """"makes a distanceRank array to store randomized values of simulated distances"""
    distanceRank = []
    for i in range(1000):
        distanceRank.append(getdistance(dinuclShuffle(sequence), mer))
        distanceRank.sort()
    return distanceRank


def drawscatter(x, y):
    """x is optimally distance and y is frequency of distance."""
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.ticker import NullFormatter

    plt.figure(1, figsize=(8, 8))
    nullfmt = NullFormatter()

#defines sizing and axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    #makes rectangle
    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x, y)

    axScatter.set_xlim((-lim, lim))
    axScatter.set_ylim((-lim, lim))
    plt.draw()
    plt.show()

def simulated_sequence_array(sequence, shuffles):
    simArray = []
    for i in range(shuffles):
        simArray.append(dinuclShuffle(sequence))
    return simArray


def get_simulated_distances(array, mer):
    distanceArray = []
    for sequence in array:
        distanceArray.append(getdistance(sequence, mer))
    distanceArray.sort()
    return distanceArray


def getsimulateddistances2(array, mer):
    distanceArray = []
    distanceFromSingleton = 0
    singleton = False
    for sequence in array:
        d = getdistance(sequence, mer)
        if singleton:
            if d[0] != -2:
                distanceArray.append(distanceFromSingleton + d[1])
                singleton = False
            else:
                distanceFromSingleton = distanceFromSingleton + len(sequence)
                try:
                    continue
                except StopIteration:
                    break

        if d[0] == -2:
            continue
        if d[0] == -1:
            distanceFromSingleton = len(sequence) - d[1]
            singleton = True
        else:
                distanceArray.append(d[0])
        
    distanceArray.sort()
    return distanceArray


def getmeans(array):
    meanList = []
    for distances in array:
        mean = sum(distances)/float(len(distances))
        meanList.append(mean)
    return meanList





def rankobs(sequence, mer):
    """works by way of indices, and one is added after calculation of index of rank."""
    d = getdistance(sequence, mer)
    distanceRank = simulateddistances(sequence, mer)
    try:
        rank = distanceRank.index(d)
        return len(distanceRank) - rank

    except ValueError:
        i = 0
        while i < len(distanceRank):
            if distanceRank[i] <= d:
                return len(distanceRank) - i
            i += 1
            return 1


def rankobslist(sequence, mer, dList, distanceRank):
    pvals = []
    for d in dList:
        try:
            rank = distanceRank.index(d)
            pvals.append(len(distanceRank) - rank)
        except ValueError:
            i = 0
            done = False
            while i < len(distanceRank) and not done:
                if distanceRank[i] <= d:
                    pvals.append(len(distanceRank) - i)
                    done = True
                i += 1
            if done:
                pvals.append(1)

    return pvals


def sequence_intake(filename):
    seqFile = open(filename, 'r')
    seqs = {}
    seqStr = ''
    for line in reversed(seqFile.readlines()):
        if '>' in line:
            seqStr = seqStr.upper()
            seqStr = transcription(seqStr)
            seqs[line] = seqStr.replace('\n', '')

            seqStr = ''
        else:
            seqStr += line
    return seqs


def boxplot(array):
    plt.boxplot(array)


def boxplot2(array):

    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.canvas.set_window_title('boxplot')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

    bp = plt.boxplot(array, notch=0, sym='+', vert=1, whis=1.5)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], color='red', marker='+')
    plt.show()


def transcription(sequence):
    sequence = sequence.replace('T', 'U')
    return sequence


def catalogue_mers(sequence):
    mers = catalogue()
    pvals = {}
    for motif in mers:
        p = rankobs(sequence, motif)
        pvals[motif] = p
    return pvals


def saveList(filename, aList):
    np.save(filename, aList)


def loadlist(filename):
    alist = np.load(filename + '.npy').tolist()
    return alist




def meanskew_sims(mers, sequenceList,  filename, shuffles = 100):
    """simulates and calculates distro factors of the given mer's spacing from a list of sequences"""
    dList = []
    dDict = {}
    sDict = {}
    for mer in mers:
        sDict[mer] = []
        dDict[mer] = {}
    simString = ""
    for seq in sequenceList:
        simString = simString + seq
    simArray = simulated_sequence_array(simString, shuffles)
    for mer in mers:
        dList = []
        for seq in sequenceList:
            simDistances = getsimulateddistances2(simArray, mer)
            dList.extend(simDistances)
            dDict[mer][seq] = simDistances
        if len(dList) != 0:
            mean = sum(dList)/float(len(dList))
            skewness = sps.skew(np.array(dList))
            variance = np.var(np.array(dList))
            # return (mean, skewness, variance)
            sDict[mer] = [mean, skewness, variance]
        else:
            print "no distances found in simulated sequences for given mer(s)"
    
    saveList(filename + '.d.sim', dList)
    saveList(filename + '.d.seq.sim', dDict)
    saveList(filename + '.stats.sim', sDict)
    return dList
    

def meanskew_obs(mers, sequenceList, filename):
    dList = []
    dDict = {}
    for mer in mers:
        dDict[mer] = {}
    for mer in mers:
        for seq in sequenceList:
            distances = recursedistances2(seq, mer)
            dList.extend(distances)
            dDict[mer][seq] = distances
    if len(dList) != 0:
        mean = sum(dList)/float(len(dList))
        skewness = sps.skew(np.array(dList))
        variance = np.var(np.array(dList))
        saveList(filename + '.d', dList)
        saveList(filename + '.d.seq', dDict)
        saveList(filename + '.stats', [mean, skewness, variance])
        return dList
    else:
        return

def distancefrequencies(distances1, distances2):
    freq = Counter(distances1)
    freq2 = Counter(distances2)
    # valFreqs = sorted(freq.items())
    drawmultiscatter(freq.keys(), freq2.keys(), freq.values(), freq2.values())



def distancefrequencies2(distances):
    freq = []
    for d in distances:
        c = Counter(d)
        freq.append([c.keys(), c.values()])
        
    drawmultiscatter2(freq)


    
def drawmultiscatter2(freq, names, stats, ks=None):
    #fo many from dictionary
    cmp=plt.get_cmap('coolwarm')

    # biggest = max(itertools.chain.from_iterable(freq))
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    j = 0
    for x, z in zip(freq, stats):
        c = []
        for i in range(4):
            y = random.randrange(100)
            c.append(float(y)/100)
        type (x)
        s = [str(round(i, 2)) for i in z]
        if ks is not None:
            plt.scatter(x[0], x[1], color=c, label=names[j] + ' mean:'+s[0] + ' var:'+s[1] + ' skew:'+s[2] + ' KSpval:' + str(round(ks[1], 2)),  marker='o', edgecolors=[0,0,0,1])
        else:
            plt.scatter(x[0], x[1], color=c, label=names[j] + ' mean:'+s[0] + ' var:'+s[1] + ' skew:'+s[2], marker='o', edgecolors=[0,0,0,1])
        # places the stats data at the respective maxima of the observed dataset along both axis
        # plt.annotate(names[j] + " mean: " + str(y), xy=(max(x[0]), max(x[1])))
        j += 1
    plt.text(.9, .9, str(y))
    plt.ylabel('frequency')
    plt.xlabel('distance')
    plt.legend()
    plt.show()
    
    # plt.scatter(x1, y1, s=100,  cmap = cmp, vmin=0, vmax=(biggest))
    # plt.scatter(x2, y2, s=100, cmap = cmp, vmin=0, vmax=(biggest))

    # plt.show()

    
def drawscatter_barren(x,y):
    cmp=plt.get_cmap('coolwarm')
    plt.scatter(x,y)



def binnednrmlzdfrequencies(distances, ks=None, view=None):
    freq = []
    names = []
    s = {}
    stats = []
    for d in distances: 
        binned = np.histogram(d[0], bins=10)
        binnorms = binned[0]/float(len(d[0]))
        binnedx = binned[1][0:len(binned[1])-1]
        freq.append([binnedx, binnorms])
        mean = sum(d[0])/len(d[0])
        skew = sps.skew(d[0])
        var = np.var(d[0])
        s[d[1]] = [mean, var, skew]
        stats.append([mean, var, skew])

    if ks is not None:
        ks = ks_2samp(distances[0][0], distances[1][0])
    for i in distances:
        names.append(i[1])
    if view is not None:
        drawmultiscatter2(freq, names, stats, ks)
        
    print 'this is stats:'
    print stats
    if ks is not None:
        return s, ks
    else:
        return s



def binnedfrequencies2(distances):
    #fo many
    freq = []
    names = []
    s = {}

    for d in distances: 
        binned = np.histogram(d[0], bins = 10)
        binnedx = binned[1][0:len(binned[1])-1]
        freq.append([binnedx, binned[0]])
        mean = sum(d[0])/len(d[0])
        skew = sps.skew(d[0])
        var = np.var(d[0])
        s[d[1]] = [mean, var, skew]
            
    for i in distances:
        names.append(i[1])
    drawmultiscatter2(freq, names)


def drawmultiscatter(x1, x2, y1, y2):
    cmp=plt.get_cmap('coolwarm')

    biggest = max(y1 + y2)
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)
    plt.scatter(x1,y1,color='red', marker = 's')
    plt.scatter(x2,y2,color='blue', marker = 'o')
    plt.show()
    
    # plt.scatter(x1, y1, s=100,  cmap = cmp, vmin=0, vmax=(biggest))
    # plt.scatter(x2, y2, s=100, cmap = cmp, vmin=0, vmax=(biggest))

    # plt.show()





def binnedfrequencies(distances1, distances2):
    binned = np.histogram(distances1, bins = 10)
    binnedx = binned[1][0:len(binned[1])-1]

    binned2 = np.histogram(distances2, bins = 10)
    binnedx2 = binned[1][0:len(binned2[1])-1]
    # drawscatter(binnedx, binned[0])
    drawmultiscatter(binnedx, binnedx2, binned[0], binned2[0])



def reducelists(listOfLists):
    bigList = list(itertools.chain.from_iterable(listOfLists))
    # bigList = [x for x in bigList]
    return bigList

    
    


def mer_ranks(mers, sequence, filename, shuffles):
    mersRankDistro = {}
    simArray = simulated_sequence_array(sequence, shuffles)
    distanceArray = {}
    for mer in mers:
        simDistances = get_simulated_distances(simArray, mer)
        dList = recursedistances2(sequence, mer)
        rankDistro = rankobslist(sequence, mer, dList, simDistances)
        mersRankDistro[mer] = rankDistro
        distanceArray[mer] = dList

    saveList(filename + 'distance', distanceArray)
    saveList(filename, mersRankDistro)

def synth_distances(distances):
    mean = np.mean(distances)
    var = np.var(distances)
    std = np.sqrt(var/4)
    length = len(distances)
    synthVals = normVals((mean-std), std, length)
    synthVals = [0+abs(x) for x in synthVals]
    return synthVals


def insert_string(insert, string, pos):
    head = string[:pos]
    tail = string[pos:]
    outString = head + insert + tail 
    return outString


def delete_substring(target, string):
    string = string.replace(target, '')
    return string


    
def testsequence(distances, sequence, mer):
    d = synth_distances(distances)
    pos = 0
    # delete_substring(mer, sequence)
    for d in distances:
        d = int(d)
        try:
            sequence = insert_string(mer, sequence, pos + d)
            pos += d
        except IndexError:
            return sequence
    return sequence

def comparemachine(sequence, insert, mers=None):
    # for the case of comparative mers
    d = []
    s = {}
    if mers is not None:
        for mer in mers:
            d.extend(recursedistances2(sequence, mer))
        d = []
        for mer in mers:
            spacing = recursedistances2(sequence, mer)
            a = [spacing, mer]
            mean = sum(spacing)/float(len(spacing))
            skew = sps.skew(spacing)
            var = np.var(spacing)
            s[mer] = [mean, var, skew]
            d.append(a)
        b = (recursedistances2(sequence, insert), insert)
        mean = sum(b)/float(len(b))
        skew = sps.skew(spacing)
        var = np.var(spacing)
        s[insert] = [mean, var, skew]
        d.append(b)


def testmachine(sequence, insert, mers=None, distances=None, ks=1, view=None):
    """the case of no comparative mers. Inserts a mer, then simulates, and compares observered distances to shuffled distances."""
    s = {}
    dList = []
    mean = sum(distances)/len(distances)
    skew = sps.skew(distances)
    var = np.var(distances)
    s[insert] = [mean, var, skew]
    sequence = testsequence(distances, sequence, insert)
    distances = recursedistances2(sequence, insert)
    print 'obs #'
    print len(distances)
    seq = simulated_sequence_array(sequence, 1000)
    # create simulated shuffles by taken the distances specified, planting the insert there,
    # shuffling it at those istances, and then finding the insert again. Not the same numbers
    # as exiting for the insert(if written) d.sim or d.seq.sim files.
    dList.append([getsimulateddistances2(seq, insert), 's' + insert])
    dList.append([distances, insert])
    
    mean = sum(dList[0][0])/len(dList[0])
    skew = sps.skew(distances)
    var = np.var(distances)
   
    print 'sim d #:'
    print len(dList[0][0])
    # mean = sum(dList[0][0])/len(dList[0][0])
    # skew = sps.skew(dList[0][0])
    # var = np.var(dList[0][0])
    # s['s' + insert] = [mean, var, skew]
    # print 'this is mean, skew, var:'
    # print s
    if view is not None:
        stats = binnednrmlzdfrequencies(dList, ks, view=1)
        return stats

    else:
         statsAndKs = binnednrmlzdfrequencies(dList, ks)
         return statsAndKs



#for the case of distances given to be interspersed.
    # else: 
    #     d = []
    #     p = 0
    #     for x in distances:
    #         try:
    #             string = insert_string(insert, sequence, p+x)
    #             p += x
    #         except IndexError:
    #             for mer in mers:
    #                 a = [recursedistances2(string, mer), mer]
    # #                 d.append(a)
    #             b = (recursedistances2(string, insert), insert)
    #             d.append(b)
    #             binnednrmlzdfrequencies
    #     for mer in mers:
    #         a = [recursedistances2(string, mer), mer]
    #         d.append(a)
    #     b = (recursedistances2(string, insert), insert)
    #     d.append(b)
    #     binnednrmlzdfrequencies
    # binnednrmlzdfrequencies(d)
    # print s


seqDict = sequence_intake(fp1640 + 'nanos_3utr.fasta.txt')
tseq = seqDict[seqDict.keys()[0]]
tseq = transcription(tseq)

seqDict = sequence_intake(fp1640 + 'sdad1.fasta.txt')
lseq = seqDict[seqDict.keys()[0]]
lseq = transcription(lseq)

seqDict = sequence_intake(fp1640 + 'cdk1.fasta.txt')
cseq = seqDict[seqDict.keys()[0]]
cseq = transcription(lseq)


seqDict = sequence_intake(fp1640 + 'ptma.fasta.txt')
pseq = seqDict[seqDict.keys()[0]]
pseq = transcription(lseq)


seqDict = sequence_intake(fp1640 + 'hcluster3h2a.fasta.txt')
mseq = seqDict[seqDict.keys()[0]]
mseq = transcription(lseq)


mer4 = catalog(4)

# print mer4[1]
# mers4 = []
# for i in mer4:
#     mers4.append(''.join(i))

# mers4 = [x.upper() for x in mers4]
# # mer_ranks(mers4, tseq, 'testseq', 10)
# mer6 = loadlist(fp6+ '6mers')
# mers6 = []
# for i in mer6:
#     mers6.append("".join(i))
# mers6 = [x.upper() for x in mers6]
#current state of mer files.proper list of capitalized strings.

# # saveList(fp1640 + '6mers', cata)
# # print '6mers saved to file 6mers.'

# nanoM = np.load(filepath + '6mersh3h2a.npy')
# h3h2aM = nanoM[()]

# edenbp = np.load(filepath + '4merscdk1.npy')
# edenbp = edenbp[()]

# pum = np.load(filepath + '4mersdad.npy')
# pum = pum[()]

# smaug = np.load(filepath + '4mernanos10000.npy')
# smaug = smaug[()]

# smaug6 = np.load(filepath + '6mernanos.npy')
# smaug6 = smaug6[()]


simArray = simulated_sequence_array(cseq, 100)
simDistances = getsimulateddistances2(simArray, 'CUGG')

# simfile = open('simDistances', 'w')
# for line in simDistances:
#     simfile.write(str(line) + '\n')

# with open(filepath + '4merscdk1distance')as g:
   # print g.next()

# meanskew_obs(mers4, [tseq, mseq, pseq, cseq, lseq], 'mers4onControls')
# x =  meanskew_sims(['CUGG'],[tseq, lseq, cseq, pseq, mseq], 'CUGG')
# y =  meanskew_sims(['AAAA'],[tseq, lseq, cseq, pseq, mseq], 'CUGG')

# meanskew_sims(mers6, [tseq, lseq, cseq, pseq, mseq], 10000)
# meanskew_obs(mers4,  [tseq, lseq, cseq, pseq, mseq], '4mersonTests')
# meanskew_obs(mers6,  [tseq, lseq, cseq, pseq, mseq], '6mersonTests')

# x = dDict['UUUU'].values()
# y = dDict['UAUA'].values()
# x = reducelists(x)
# y = reducelists(y)

# test = reducelists(eq.keys())

# sauuu = [retrieve_distances('AUUU', sDict), 'SAUUU']

# cugg = retrieve_distances('CUGG', dDict)

# x =  binnednrmlzdfrequencies([[scugg, 'SCUGG'],[cugg, 'CUGG']])
# binnednrmlzdfrequencies([auuu, sauuu])

dDict = loadlist(fpmers + '4mersonTests.d.seq')
# sDict = loadlist(fpmers + '4mersonTests.d.seq.sim')


def retrieve_distances(mer, Dict):
    x = Dict[mer].values()
    d = reducelists(x)
    return d


def retrieve_sequences(mer, Dict):
    x = Dict[mer].keys()
    seq = reduce((lambda x, y: x + y), x)
    return seq


testSeq = dDict['GCAG'].keys()
tests = ""
testSeq = "".join(testSeq)


distances = retrieve_distances('GCAG', dDict)
print 'obs dist:'
print distances
d = synth_distances(distances)
d = d[0:2]

kslist = []

for i in range(1,2):
    x = (testmachine(testSeq, 'GCAG', distances=d, view = 1))
    x2 = x[1]
# gmm = mixture.GaussianMixture(n_components=2).fit(d)
# gmmy = gmm.predict(d)
# drawscatter_barren(gmm,gmmy)

