# This file is part of Allele-Wrangler.
#
# Allele-Wrangler is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Allele-Wrangler is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Allele-Wrangler. If not, see <http://www.gnu.org/licenses/>.
from matplotlib.sphinxext.tests.test_tinypages import HERE

# Version 1.0 

SoftwareVersion = "Allele-Wrangler Version 1.0"

import sys
import pysam
import os
from os.path import split, join, isdir
from os import mkdir, makedirs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from shutil import copyfile
#from shutil import copy

class AlleleWrangler():   
    
    def __init__(self, readsFile, outputDir, referenceFile, numIters, numThreads):
         
        print ('Setting up the Allele Wrangler...')
        self.readInput          = readsFile
        self.inputDirectory = False # Is the input a directory? If not, it is a single file.
        self.outputRootDirectory = outputDir
        self.referenceSequenceFileName     = referenceFile
        self.totalIterations      = numIters
        self.numberThreads         = numThreads
        
        # Determine if input is a file, or a directory
        if (os.path.isfile(self.readInput)):
            #print ('Read input is a file that exists.')
            pass
            
            # Determine Input File Type
            if (".fasta" == self.readInput[-6:] or ".fa" == self.readInput[-3:]):
                self.readInputFormat = "fasta"
            elif (".fastq"== self.readInput[-6:] or ".fq" == self.readInput[-3:]):
                self.readInputFormat = "fastq"
            else:
                print('I expect a .fasta or .fastq format for read input. Alternatively, specify a directory containing read inputs. Please check your input.')
                raise Exception('Bad Read Input Format')
                
        elif (os.path.isdir(self.readInput)):
            #print ('Read input is a directory that exists.')
            self.inputDirectory = True
            
        else :
            print('I expect a .fasta or .fastq format for read input. Alternatively, specify a directory containing read inputs. Please check your input.')
            raise Exception('Bad Read Input Format')
        
    # TODO Copy constructor
    # Could come in handy for recursing...


    def wrangle(self):
        
        print ('Wrangling.')
        # Create a LogFile
        self.wranglerLog = createOutputFile(join(self.outputRootDirectory,'wrangling_summary.txt'))
        
        self.wranglerLog.write('Read Input:' + self.readInput + '\n')
        self.wranglerLog.write('Result directory:' + str(self.outputRootDirectory) + '\n')
        self.wranglerLog.write('Number Threads:' + str(self.numberThreads) + '\n')
        
        if(self.inputDirectory):
            self.wranglerLog.write('Read Input is a Directory:' + str(self.inputDirectory) + '\n')
            for currentInputReadFile in os.listdir(self.readInput):
                self.wranglerLog.write(currentInputReadFile + '\n')
            # for each file name, print the file
            self.wranglerLog.write('\t:' + str(self.inputDirectory) + '\n')
        else:
            self.wranglerLog.write('Read Input Format:' + self.readInputFormat + '\n')
            

      
            #lklkjvsldkfj
            #need to make new reference for each subgroup in the loop.
        
        if(self.inputDirectory):
        # If it's a directory, then for each file
            for currentInputReadFile in os.listdir(self.readInput):
                self.wranglerLog.write(currentInputReadFile + '\n')
                print('I will start a wrangler for these reads:' + str(currentInputReadFile))
                #consensusSequenceFileName = self.wrangleRecurser(None, 1)            
                #self.summarizeWrangling(consensusSequenceFileName)
                
                
                
                
                
            # Speficy output directory
            
            # Specify file format
            # Do the normal stuff.
        
        else:
            
            #if(self.referenceSequenceFileName is None):
            #    print('No Reference was provided.  No problem, I will use the first read as a reference.')
            #    self.createFirstGuessReferenceFromReads()            
            #else:
                #print('A Reference sequence was provided.')
            
            # If read file input, do all the normal stuff
            #self.currentIteration = 1
            consensusSequenceFileName = self.wrangleRecurser(None, 1)            
            self.summarizeWrangling(consensusSequenceFileName)

        
        
        
        
        # Was a consensus sequence provided?
        #if(self.referenceSequenceFileName is None):
        #    print('No Reference was provided.  No problem, I will use the first read as a reference.')
        #    self.createFirstGuessReferenceFromReads()            
        #else:
        #    print('A Reference sequence was provided.')
            #self.prepareReferenceInput() 
            #Copy the Reference to the input folder
            #consensusSequence = self.openConsensusSequence()
            
        #self.currentIteration = 1
        #consensusSequenceFileName = self.wrangleRecurser(self.referenceSequenceFileName)
        
        #self.summarizeWrangling(consensusSequenceFileName)
        
        self.wranglerLog.close()
        
        
        #TODO: After iteration, do A "FINAL" alignment with results.
        #I can report heterozygous positions and all that in there.
     
    def summarizeWrangling(self, finalConsensusSequence):  
        print('\n\nSummarizing Wrangling Results.') 
        try:
            summaryDirectory = join(self.outputRootDirectory, 'FinalConsensus')
                        
            self.alignReads(finalConsensusSequence,self.readInput,summaryDirectory)
            
            copyfile(finalConsensusSequence, join(self.outputRootDirectory, 'ConsensusSequence.fasta'))
            
            # copy alignment reference and call it "final" 
            
        except Exception:
            print ('Exception performing the summarize alignment')                  
            raise 
    
    # Input = location of Reference Sequence, iteration
    # Output = Location of the Consensus Sequence
    def wrangleRecurser(self, currentReferenceSequence, currentIteration):
        print('\n\nWRANGLING SOME ALLELES, ITERATION ' + str(currentIteration))
        try:
            
            # Was a consensus sequence provided?
            if(currentReferenceSequence is None):
                print('No Reference was provided.  No problem, I will use the first read as a reference.')
                currentReferenceSequence = self.createFirstGuessReferenceFromReads()            
            else:
                print('A Reference sequence was provided.')
                #self.prepareReferenceInput() 
                #Copy the Reference to the input folder
                #consensusSequence = self.openConsensusSequence()            
                       
            alignmentSubdir = join(self.outputRootDirectory,'alignments')
            if not isdir(alignmentSubdir):
                mkdir(alignmentSubdir)
            currentIterationSubdirectory = join(alignmentSubdir,'iter_'+ str(currentIteration))
            if not isdir(currentIterationSubdirectory):
                mkdir(currentIterationSubdirectory)
            
            # Write a bit about this iteration to the Log.

            self.wranglerLog.write('\nIteration # (' + str(currentIteration) + '/' + str(self.totalIterations) + ')\n')
            self.wranglerLog.write('Iteration Alignment Directory:' + currentIterationSubdirectory + '\n')
            self.wranglerLog.write('Reference Filename:' + str(currentReferenceSequence) + '\n')
            
    
            
                
            self.alignReads(currentReferenceSequence,self.readInput,currentIterationSubdirectory)
            self.analyzeAlignment(currentIterationSubdirectory)
            
            # If we want more iterations, I should Recurse and try again.
            if (int(self.totalIterations) > int(currentIteration)):
                # On the next iteration, we want to use the new consensus sequence 
                # as the reference. 
                print('I am on iteration (' + str(currentIteration) + '/' + self.totalIterations + ') I will continue...')
                newReferenceSequenceFileName = join(currentIterationSubdirectory, 'Consensus.fasta')
                #currentIteration += 1
                
                # Return the consensus from one layer deeper.
                return self.wrangleRecurser(newReferenceSequenceFileName, currentIteration + 1)
                
                
            else:
                print('That was the last iteration (' + str(currentIteration) + '/' + self.totalIterations + '), I am done now.')
                
                # Return the consensus
                return join(currentIterationSubdirectory, 'Consensus.fasta')
                
            
            
            
        except Exception:
            print ('Exception encountered in wrangle()')                  
            print sys.exc_info()[0]
            print sys.exc_info()[1]
            print sys.exc_info()[2]        
            raise 
        
        
    # Perform BW Alignment.  Align all reads against the Reference.
    def alignReads(self, referenceLocation, readFileLocation, alignmentOutputDirectory):
        print('\nStep 1.) Aligning reads against the reference.')
        
        if not isdir(alignmentOutputDirectory):
            mkdir(alignmentOutputDirectory)
        
        # Part 1 Index the Reference        
        try:
            # Copy the reference sequence to the alignment directory. This is a complicated way to do it.
            newReferenceLocation = join(alignmentOutputDirectory,'AlignmentReference.fasta')
            refSequence = list(SeqIO.parse(referenceLocation, 'fasta'))[0]
            refSequence.id = 'AlignmentReference'
            sequenceWriter = createOutputFile(newReferenceLocation)
            SeqIO.write([refSequence], sequenceWriter, 'fasta')
            sequenceWriter.close()
                        
            # Index The Reference
            cmd = ("bwa index " + newReferenceLocation)
            os.system(cmd)
            
        except Exception:
            print ('Exception indexing alignment reference. Is bwa installed? folder writing permission issue?')                  
            raise 
        
        # Part 2 Align
        try:
            # align | sam->bam | sort
            tempAlignmentName = join(alignmentOutputDirectory,'alignment')
            bwaMemArgs = "-t " + str(self.numberThreads) + " -x ont2d"
            cmd = ("bwa mem " + 
                bwaMemArgs + " " +  
                newReferenceLocation + " " +
                readFileLocation + 
                " | samtools view  -Sb - | samtools sort - "
                + tempAlignmentName)
            #print ('alignment command:\n' + cmd)
            os.system(cmd)
            alignmentOutputName = tempAlignmentName + '.bam'
            
        except Exception:
            print ('Exception aligning reads against reference. Are bwa and samtools installed?')                  
            raise 
        
        # Part 3 Index Alignment
        try:
            cmd = ("samtools index " + alignmentOutputName)
            #print ('alignment index command:\n' + cmd)
            os.system(cmd)
            #print ('index command:\n' + cmd)
        except Exception:
            print ('Exception indexing alignment reference. Is bwa installed?')                  
            raise 
  
    

    def analyzeAlignment(self, alignmentOutputDirectory):
        print ('\nStep 2.) Parse the alignment and create a new consensus sequence.')
        
        # Load up the Alignment Reference file, we'll need it.
        alignmentReferenceFileName = join(alignmentOutputDirectory,'AlignmentReference.fasta')
        alignmentRef = list(SeqIO.parse(alignmentReferenceFileName, 'fasta'))[0]
        
        # Count the reads in the input file
        totalReadCount = len(list(SeqIO.parse(self.readInput, self.readInputFormat)))
        #self.readInputFormat
        #self.readInput
                
        # We generate a new consensus sequence from the alignment results.
        newConsensusSequence = ""
        
        # Open the bam file
        bamfile = pysam.AlignmentFile(join(alignmentOutputDirectory,'alignment.bam'), 'rb')  
        
        # Open alignment analysis text file
        alignmentSummaryFile = createOutputFile(join(alignmentOutputDirectory,'AlignmentSummary.csv')) 
        alignmentSummaryFile.write('Ref_Position,Ref_Base,Reference_Adjustment,Aligned_Count,Unaligned_Count,Match_Count,Mismatch_Count,In_Count,Del_Count,A_Count,G_Count,C_Count,T_Count\n')
        
        # A smaller log. I will provide human-readable descriptions of the
        # bases that were adjusted in the new consensus sequence.
        # TODO: Provide surrounding sequence as well, maybe it's a repeat region....
        # Acutally NAH, I want to just put it in the wrangler log. 
        #adjustedBasesSummaryFile = createOutputFile(join(alignmentOutputDirectory,'AdjustedBases.txt')) 
        
        # Keep a running total of adjustments made to the reference.
        # If this total is 0, then theoretically the consensus matches the alignment reference, and we're done.
        totalSequenceAdjustments = 0
        
        # Iterate the reference sequence column by column.
        pileupIterator = bamfile.pileup(alignmentRef.id)
        # TODO: Check where this pileup iterator starts. Are there reads mapped before or after the reference begins/ends?
        for pileupColumn in pileupIterator:
            
            #
            referencePosition = 0
            referenceBase = ''
            referenceAdjustment = '?'
            alignedCount = 0
            unalignedCount = 0
            matchCount = 0
            mismatchCount = 0
            inCount = 0
            delCount = 0
            aCount = 0
            gCount = 0
            cCount = 0
            tCount = 0
            
            referencePosition = pileupColumn.reference_pos
            referenceBase = alignmentRef[pileupColumn.reference_pos].upper()
            alignedCount = pileupColumn.nsegments
            unalignedCount = totalReadCount - alignedCount
            
            # Iterate the Reads at this position           
            for pileupRead in pileupColumn.pileups:
                
                # If this read is a deletion
                if(pileupRead.is_del == 1):
                    delCount += 1
                # else if this read is an insertion
                elif(pileupRead.indel > 0):
                    #print ('INSERTION DETECTED, INDEL=' + str(pileupRead.indel))  
                    inCount += 1                   
                # Else if it is a refskip (TODO What does this mean? no read aligned? Count these?)
                elif(pileupRead.is_refskip):
                    print('This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)
                    raise Exception('This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)
                # else this means we have a base aligned at this position for this read.
                else:    
                    currentBase = pileupRead.alignment.query_sequence[pileupRead.query_position].upper()                    
                    #print('Reference,Current:' + referenceBase + ',' + currentBase)
                    #print('Curr')
                    if(currentBase == referenceBase):
                        matchCount += 1
                    else:
                        mismatchCount += 1
                   
                # Count the nucleotide 
                if (currentBase == 'A'):
                    aCount += 1
                elif (currentBase == 'G'):
                    gCount += 1
                elif (currentBase == 'C'):
                    cCount += 1
                elif (currentBase == 'T'):
                    tCount += 1
                else:
                    print('Unknown Base found in Alignment at position ' + str(referencePosition) + ':' + currentBase)
                    raise Exception('Unknown Base in Alignment')
                
                
                # TODO: What if the query insertion sequence is longer than one base?
                # Maybe I can only adjust one base per iteration, is that okay? Probably for the Best, actually..
                # Don't worry bout it for now.
            
            # Calculate highest frequency base
            # I hope this algorithm makes sense, probably there is a smarter way to do it.
            if(aCount >= gCount and aCount >= cCount and aCount >= tCount):
                mostFrequentBase = 'A'
                mostFrequentBaseCount = aCount
            elif(gCount >= cCount and gCount >= tCount):
                mostFrequentBase = 'G'
                mostFrequentBaseCount = gCount
            elif(cCount >= tCount):
                mostFrequentBase = 'C'
                mostFrequentBaseCount = cCount
            else:
                mostFrequentBase = 'T'
                mostFrequentBaseCount = tCount
            
            #TODO: Detect heterozygosity here
            # Do the base frequencies look "normal"?
            # High proportion of Inserts or Deletions?
            # Maybe I don't care, because I want to build a read-clusterer tool.
            
            
            # Add the next base to the new consensus sequence            
            if (matchCount >= mismatchCount and matchCount >= inCount and matchCount >= delCount):
                # Aligned bases match the reference, add reference base to the consensus.
                referenceAdjustment='-'
                newConsensusSequence += referenceBase
                
            elif (inCount >= mismatchCount and inCount >= delCount):
                # Aligned bases show an insertion.
                # Add the Reference Base and the Insertion Base to the consensus.  
                totalSequenceAdjustments += 1 
                referenceAdjustment='I'  
                newConsensusSequence += referenceBase + mostFrequentBase         
                
                self.wranglerLog.write(str(referencePosition) + ':Insertion' +
                    '\n(' + str(inCount) + '/' + str(alignedCount) + ') = ' + str((100.0 * inCount) / alignedCount) + '% of aligned reads'
                    '\n(' + referenceBase + ' > ' + referenceBase + mostFrequentBase + ')' +
                    '\n')
                
                #TODO: I need to insert multiple bases, if that is waht the alignment suggests.

            elif (delCount >= mismatchCount):
                # Reads show a deletion.
                # Don't add anything to the consensus.
                totalSequenceAdjustments += 1
                referenceAdjustment='D'
                
                self.wranglerLog.write(str(referencePosition) + ':Deletion' +
                    '\n(' + str(delCount) + '/' + str(alignedCount) + ') = ' + str((100.0 * delCount) / alignedCount) + '% of aligned reads'
                    '\n(' + referenceBase + ' > _)' +
                    '\n')
                
            else:
                # Mismatch base.
                # Add the highest read count base to the reference.
                # It might actually be the same base as the reference,
                # Because this just means there are more mismatches than matches.
                # Problematic base, at least we'll notice here.
                # TODO: What to do with highly heterozygous Positions?
                # I should report those that look particularly heterozygous, somewhere.
                newConsensusSequence += mostFrequentBase 
                totalSequenceAdjustments += 1     
                referenceAdjustment='M'   
                
                self.wranglerLog.write(str(referencePosition) + ':Mismatch' +
                    '\n(' + str(mostFrequentBaseCount) + '/' + str(alignedCount) + ') = ' + str((100.0 * mostFrequentBaseCount) / alignedCount) + '% of aligned reads'
                    '\n(' + referenceBase + ' > ' + mostFrequentBase + ')' +
                    '\n')
              

            # Write a line to the alignment Summary 
            alignmentSummaryFile.write(str(referencePosition) + 
                ',' + str(referenceBase) +
                ',' + str(referenceAdjustment) + 
                ',' + str(alignedCount) + 
                ',' + str(unalignedCount) + 
                ',' + str(matchCount) + 
                ',' + str(mismatchCount) + 
                ',' + str(inCount) + 
                ',' + str(delCount) + 
                ',' + str(aCount) + 
                ',' + str(gCount) + 
                ',' + str(cCount) + 
                ',' + str(tCount) +
                '\n')
            
        print('\nTotal Sequence Adjustments:' + str(totalSequenceAdjustments) + ' (How many bases the consensus differs from the reference.)\n')    
        
        # Write the newly constructed consensus sequence.
        currentConsensusSequenceFileName = join(alignmentOutputDirectory, 'Consensus.fasta')        
        consensusWriter = createOutputFile(currentConsensusSequenceFileName)          
           
        SeqIO.write([SeqRecord(Seq(newConsensusSequence,
            IUPAC.unambiguous_dna),
            id="GeneratedConsensusSequence|Coverage=GarbageInformation", description="") ], consensusWriter, 'fasta')
        consensusWriter.close()
            
        self.wranglerLog.write('Total Sequence Adjustments:' + str(totalSequenceAdjustments) + '\n')
            
        # Close Summary Files
        alignmentSummaryFile.close()
        #adjustedBasesSummaryFile.close()
        
        #return totalSequenceAdjustments
   
    def createFirstGuessReferenceFromReads(self):        
        # TODO: I should choose better than just the first sequence.  The longest?  Average size? Consensus of a handful of reads?
        try:

            # Load Reads from File
            parsedReads = list(SeqIO.parse(self.readInput, self.readInputFormat))
            firstRead = parsedReads[0]
            
            #print('First Read:' + str(firstRead))
            
            # Write the first read, for use as an alignment reference.
            referenceDirectory = join(self.outputRootDirectory,'Initial_Reference')
            if not isdir(referenceDirectory):
                mkdir(referenceDirectory)
                
            referenceFileName = join(referenceDirectory, 'FirstGuessReference.fasta')
            
            firstGuessRefFileWriter = createOutputFile(referenceFileName)        
            SeqIO.write([firstRead], firstGuessRefFileWriter, 'fasta')
            firstGuessRefFileWriter.close()
            
            return referenceFileName
       
        except Exception:
            print ('Exception encountered in createFirstGuessReferenceFromReads()') 
            raise  




# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput