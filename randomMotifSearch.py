#!/usr/bin/env python
# coding: utf-8

# **CLASS FASTAREADER**
#
# This class is for reading fasta files and creating a manipulatable file.
#
# Functions:
#     **\__init__** - class constructor
#

# In[1]:


########################################################################
# Class: fastaReader
#
# Purpose:class to read .fa files
#
#
# Author: Chris Nguyen
# History:      cdn 10/10/21 Created
#
########################################################################


class fastaReader:

    def __init__(self, inFile):
        '''
        fastaReader constructor. Takes a file as input, creates a fastaReader object.
        '''
        self.file = inFile
        self.fastaFile = []
        self.headerFile = []
        header = ''
        sequence = ''

        with open(inFile) as fileR:
            for line in fileR:
                if line[0] == ">":

                    line = line.split()
                    header = line[0][1: len(line[0])]
                    self.headerFile.append(header)
                elif line[0] != ">":

                    sequence = line.rstrip("\n")
                    self.fastaFile.append(sequence)


# **CLASS RANDOMIZEDSEARCHMOTIFER**
# This class performs the actual searching.
#
# Functions:
#
#     __init__(self, DNA, kmer, iterate, pseudo,  gibbs = False) : class constructor
#     makeNull(self) : makes the null model
#     countAndProfile(self, motifMatrix) : creates count and profile matrix
#     replaceMotif(self, profileMatrix, motifMatrix, index) : for gibbs sampling, replaces a random motif
#     fullCalc(self, profileMatrix) : replaces all motifs with better ones
#     scoreMatrix(self, profileMatrix) : scores a profile matrix
#     getConsensus(self, motifMatrix) : gets consensus sequence
#     performSearch(self) : does the full search

# In[2]:


########################################################################
# Class: randomizedSearchMotifer
#
# Purpose:class to do the randomized searching for motifs
#
#
# Author: Chris Nguyen
# History:      cdn 10/10/21 Created
#
########################################################################
class randomizedSearchMotifer:

    def __init__(self, DNA, kmer, iterate, pseudo, gibbs=False):
        '''
        Constructor for the randomizedSearchMotifer class.
        Takes in parameters and stores them in member variables.

        Parameters:
            DNA (list) - set of data to search from
            kmer (int) - length of kmer to use
            iterate (int) - how many times to do the algorithm
            pseudo (float) - which value to use for pseudocount
            gibbs (bool) - boolean to see if gibbs sampling should be done
        '''
        self.d = DNA
        self.g = gibbs
        self.k = kmer
        self.i = iterate
        self.p = pseudo

        self.null = self.makeNull()
        self.nullScore = 0

        # creates null model to compare values

    def makeNull(self):
        '''
        Class function that makes the null model.

        Returns the null model as a list.

        '''
        # initializes local variables
        numA = numC = numT = numG = 0
        # loop that counts how many instances of a base there are
        for x in self.d:
            numA += x.count('a')
            numC += x.count('c')
            numT += x.count('t')
            numG += x.count('g')

        # after which, a total is yielded from the counts
        total = numA + numC + numG + numT
        pA = ((numA + self.p) / (total + self.p * len(self.d)))
        pC = ((numC + self.p) / (total + self.p * len(self.d)))
        pT = ((numT + self.p) / (total + self.p * len(self.d)))
        pG = ((numG + self.p) / (total + self.p * len(self.d)))

        return [pA, pC, pT, pG]

    def countAndProfile(self, motifMatrix):
        '''
        Class function that creates count matrix and subsequently creates profile matrix.

        Parameters:
            motifMatrix (list) - motifs to count and profile

        Return:
            profileMatrix (list) - profile matrix
        '''
        # initialization of local variables
        profileMatrix = []
        numA = 0
        numC = 0
        numT = 0
        numG = 0

        # loop that counts instances of nucleotide
        for j in range(0, self.k):
            for i in motifMatrix:
                if i[j] == 'a':
                    numA += 1
                elif i[j] == 'c':
                    numC += 1
                elif i[j] == 't':
                    numT += 1
                elif i[j] == 'g':
                    numG += 1

            # gets profile probability for a row, and then appends it

            numA = (numA + self.p) / ((4 * self.p) + len(motifMatrix))
            numC = (numC + self.p) / ((4 * self.p) + len(motifMatrix))
            numT = (numT + self.p) / ((4 * self.p) + len(motifMatrix))
            numG = (numG + self.p) / ((4 * self.p) + len(motifMatrix))
            profileMatrix.append([numA, numC, numT, numG])

            # zeros out the variables for re-use
            numA = numC = numT = numG = 0

        return profileMatrix

    def replaceMotif(self, profileMatrix, motifMatrix, index):
        '''
        Class function that replaces a random motif in the motif matrix with the most profile-supported one.
        Only used for Gibbs Sampling.

        Parameters:
            profileMatrix (list) - matrix of profile probabilities
            motifMatrix (list) - matrix of motifs being used
            index (int) - location of motif being replaced in motifMatrix

        Return:
            motifMatrix (list) - motif matrix with one motif mutated

        '''
        # initialization of local variables
        bestKmer = ''
        probNew = 0

        # loop chooses kmers for consideration
        for i in range(0, (len(self.d[index]) - self.k + 1)):
            kmer = self.d[index][i: (i + self.k)]
            prob = 1

            # probability created by comparing it to profile matrix
            for x in range(0, (len(kmer))):
                if kmer[x].upper() == "A":

                    prob = prob * profileMatrix[x][0]
                elif kmer[x].upper() == "C":

                    prob = prob * profileMatrix[x][1]
                elif kmer[x].upper() == "T":

                    prob = prob * profileMatrix[x][2]
                elif kmer[x].upper() == "G":

                    prob = prob * profileMatrix[x][3]

            # if it's better, it replaces the best one so far
            if prob > probNew:
                probNew = prob

                bestKmer = kmer

        # replaces motif in matrix once overall best one is found
        motifMatrix[index] = bestKmer
        return motifMatrix

    # calculates k-mer probabilities for all k-mers in the matrix
    # input takes profileMatrix and motifMatrix and returns a better motifMatrix
    def fullCalc(self, profileMatrix):

        '''
        Class function that replaces all motifs in the motif matrix with the most profile-supported one.


        Parameters:
            profileMatrix (list) - matrix of profile probabilities


        Return:
            motifMatrix (list) - new set of motifs

        '''
        # initialization of local variables
        bestKmer = ''
        probBest = 0
        motifMatrix = []

        # loop for choosing a string from the whole data set
        for string in self.d:
            # loop for picking all kmers from a string
            for i in range(0, (len(string) - self.k + 1)):

                kmer = string[i: (i + self.k)]

                prob = 1

                # sees what its probability is by comparing it to profile matrix
                for x in range(0, (len(kmer) - 1)):

                    if kmer[x].upper() == "A":

                        prob = prob * profileMatrix[x][0]

                    elif kmer[x].upper() == "C":

                        prob = prob * profileMatrix[x][1]
                    elif kmer[x].upper() == "T":

                        prob = prob * profileMatrix[x][2]
                    elif kmer[x].upper() == "G":

                        prob = prob * profileMatrix[x][3]

                # if better, replaces the best one so far
                if prob > probBest:
                    probBest = prob

                    bestKmer = kmer

            # best one gets put into the matrix and probability is zeroed out for next calculation
            motifMatrix.append(bestKmer)
            probBest = 0

        # replaces motif in matrix

        return motifMatrix

    def scoreMatrix(self, profileMatrix):
        '''
        Class function that scores a profile matrix using relative entropy.

        Parameters:
            profileMatrix (list) - matrix of profile probabilities

        Return:
            relEnt (float) - score of matrix

        '''
        # initialization of local variable
        relEnt = 0
        # loop that calculates relative entropy by visiting all elements in matrix
        for x in range(0, len(profileMatrix)):
            for y in range(0, len(profileMatrix[x])):
                relEnt += profileMatrix[x][y] * math.log2((profileMatrix[x][y]) / (self.null[y]))

        return relEnt

    def getConsensus(self, motifMatrix):
        '''
        Class function that produces a consensus sequence for a given motif matrix.

        Parameters:
            motifMatrix (list) - matrix of motifs being used

        Return:
            conString (str) - string containing consensus sequence

        '''
        # initialization of local variables
        numA = numC = numT = numG = 0
        conString = ''

        # counts nucleotide occurence for motif matrix
        for j in range(0, self.k):
            for i in motifMatrix:

                if i[j] == 'a':
                    numA += 1
                elif i[j] == 'c':
                    numC += 1
                elif i[j] == 't':
                    numT += 1
                elif i[j] == 'g':
                    numG += 1

            consCount = max(numA, numC, numT, numG)

            # if the largest expression number matches the count for a given base, it's written to the consensus
            if consCount == numA:

                conString += 'A'
            elif consCount == numC:

                conString += 'C'
            elif consCount == numT:

                conString += 'T'
            elif consCount == numG:

                conString += 'G'

            # zeroes out variables for next round
            numA = numC = numT = numG = 0

        return conString

    def performSearch(self):
        '''
        Class function that does the entire search, using Gibbs if specified


        Return:
            finalScore (float) - final true best score for the best motifs chosen from the data
            finalConsensus (str) - final true best consensus for the file
            finalMotifs (list) - final true best set of motifs across all iterations

        '''

        # initialization of local variables
        finalScore = 0
        finalConsensus = ''
        finalMotifs = []
        # for a given number of iterations, performs the search
        for i in range(0, self.i):
            # iteration-local variables
            motifs = []
            newMat = []

            # loops through each sequence in the DNA
            for seq in self.d:
                # makes indices for k-mers

                start = random.randint(0, len(seq) - self.k)
                end = start + self.k + 1  # this makes it inclusive
                kmer = seq[start:end]

                motifs.append(kmer)

            # then, puts them all into a set of motifs called "bestMotifs"
            bestMotifs = motifs

            # now, we have everything we need to do either method.
            # if self.g is True, do Gibbs
            if self.g:

                # for gibbs, the iterations are used again for the inner calculation
                for i in range(0, self.i):

                    # randomly choosing a motif to discard

                    bString = random.randint(0, len(self.d) - 1)

                    # make new matrix without the string to pass to function
                    for h in range(0, len(self.d)):

                        # if index is not equal to the random int, put into new matrix
                        if h != bString:
                            newMat.append(bestMotifs[h])

                    # after, we pass into the function to make a profile
                    profile = self.countAndProfile(newMat)

                    # use the resulting profile to replace a motif
                    motifs = self.replaceMotif(profile, motifs, bString)

                    # if the score of the new matrix is better than the current best, it replaces the current best
                    if self.scoreMatrix(self.countAndProfile(motifs)) > self.scoreMatrix(
                            self.countAndProfile(bestMotifs)):
                        bestMotifs = motifs

                bestProfile = self.countAndProfile(bestMotifs)

                # at the end, the overall best ones get chosen and returned
                finalConsensus = self.getConsensus(bestMotifs)
                finalScore = self.scoreMatrix(bestProfile)
                finalMotifs = bestMotifs
                return finalScore, finalConsensus, finalMotifs


            # if self.g is False, do the search without Gibbs
            else:

                bestProfile = self.countAndProfile(bestMotifs)

                while True:

                    # use profile to choose another set of motifs
                    newMotifs = self.fullCalc(bestProfile)
                    newProfile = self.countAndProfile(newMotifs)

                    # if the new matrix score is better, it replaces the previous best
                    if self.scoreMatrix(newProfile) > self.scoreMatrix(bestProfile):

                        bestMotifs = newMotifs
                        bestProfile = newProfile

                    # otherwise, the best has been found
                    else:

                        bestMotifs = newMotifs
                        bestProfile = newProfile
                        bestConsensus = self.getConsensus(bestMotifs)
                        bestScore = self.scoreMatrix(bestProfile)
                        break

                # if the best score from an iteration is better than the overall best, it gets replaced
                if finalScore == 0 or finalScore < bestScore:
                    finalMotifs = bestMotifs
                    finalScore = bestScore
                    finalConsensus = bestConsensus

        return finalScore, finalConsensus, finalMotifs


# **CLASS COMMANDLINE**
#
# This class specifies a commandline to use for the program.
#
# Functions:
#
#     __init__(self, inOpts=None) : class constructor

# In[3]:


########################################################################
# Class: CommandLine
#
# Purpose:class to represent a command line
#
#
# Author: Chris Nguyen
# History:      cdn 10/10/21 Created
#
########################################################################
import argparse


class CommandLine():

    def __init__(self, inOpts=None):
        '''
        Class constructor that creates the command line.

        '''

        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        # these arguments can be accessed or seen using "args.(argument name)" so args.list would hold the list of values passed to the option
        self.parser.add_argument('-i', '--iterations', type=int, action='store', required=True,
                                 help='iterations to perform')
        self.parser.add_argument('-k', '--motifLength', type=int, choices=range(5, 48), action='store', required=True,
                                 help='length')
        self.parser.add_argument('-p', '--pseudocount', type=float, action='store', required=True,
                                 help='designates pseudocount value')
        self.parser.add_argument('-g', '--gibbs', action='store_true', help='performs gibbs sampling')
        self.parser.add_argument('-m', '--printMotif', action='store_true', help='print motifs of consensus')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


# **MAIN**
#
# This contains main, which makes object instances and runs the search.

# In[4]:


########################################################################
# File: Main
#
# Purpose:main function for randomizedMotifSearch program
#
#
# Author: Chris Nguyen
# History:      cdn 10/10/21 Created
#
########################################################################
import random
import math


def main(inFile, myCommandLine=None):
    import sys

    cLine = CommandLine(myCommandLine)
    kmer = cLine.args.motifLength
    pseudo = cLine.args.pseudocount
    iterate = cLine.args.iterations
    gibbs = cLine.args.gibbs
    printOut = cLine.args.printMotif

    reader = fastaReader(inFile)
    DNA = reader.fastaFile
    headers = reader.headerFile

    # need to read values in from command line to affect constructor
    motifer = randomizedSearchMotifer(DNA, kmer, iterate, pseudo, gibbs)
    score, consensus, motifs = motifer.performSearch()
    print("Best score: " + str(score))
    print("Consensus: " + consensus)

    if printOut:
        for x in range(0, len(DNA)):
            print("Name of organism: " + headers[x])
            print("DNA sequence: " + DNA[x])
            print("Motif: " + motifs[x])


if __name__ == "__main__":
    main('allFilesCrisprs.fa', ['-i 1000', '-p 1', '-k 13'])

