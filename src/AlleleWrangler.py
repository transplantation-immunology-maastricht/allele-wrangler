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
import random
from os.path import split, join, isdir
from os import mkdir, makedirs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
#from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Sequencing.Applications import BwaIndexCommandline
from Bio.Align import AlignInfo
from Bio import AlignIO

from subprocess import Popen, PIPE, STDOUT


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


    def analyzeReads(self):
        
        print ('Beginning analysis of reads.')
        # Create a LogFile
        self.wranglerLog = createOutputFile(join(self.outputRootDirectory,'analysis_summary.txt'))
        
        self.wranglerLog.write('Read Input:' + self.readInput + '\n')
        self.wranglerLog.write('Result directory:' + str(self.outputRootDirectory) + '\n')
        self.wranglerLog.write('Number Threads:' + str(self.numberThreads) + '\n')
        

        # Was a consensus sequence provided?
        if(self.referenceSequenceFileName is None):
            print('No Reference was provided.  No problem, I can generate my own reference.')
            self.createFirstGuessReferenceFromReads()            
>>>>>>> referencemethod
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
            

        self.currentIteration = 1
        consensusSequenceFileName = self.analysisRecurser(self.referenceSequenceFileName)
        
        # TODO: splitHeterozygotes should be commandline parameter
        self.summarizeAnalysis(consensusSequenceFileName, True)
>>>>>>> referencemethod
        
        self.wranglerLog.close()
        
        
        #TODO: After iteration, do A "FINAL" alignment with results.
        #I can report heterozygous positions and all that in there.
     
    def summarizeAnalysis(self, finalConsensusSequence, splitHeterozygotes):  
        print('\n\nSummarizing Analysis Results.') 
        try:
            summaryDirectory = join(self.outputRootDirectory, 'FinalConsensus')
                        
            self.alignReads(finalConsensusSequence,self.readFileName,summaryDirectory, False)
>>>>>>> referencemethod
            
            if(splitHeterozygotes):
                # Find heterozy
                self.splitHeterozygousPositions()
                
            
            else:
                # write consensus to file. 
                # copy alignment reference and call it "final" 
                
                pass
                        
        except Exception:
            print ('Exception performing the summarize alignment')                  
            raise 
    

    
    
    def splitHeterozygousPositions(self):
        print('Splitting reads by heterozygous positions')
        # Heterozygous base list
        heterozygousBasesSummaryFile = createOutputFile(join(self.outputRootDirectory,'HeterozygousBases.txt')) 
        
        heterozygousBasesSummaryFile.close()
        # Find positions
        
        # Open the bam file
        #bamfile = pysam.AlignmentFile(join(alignmentOutputDirectory,'alignment.bam'), 'rb')  
        
        # Iterate the reference sequence column by column.
        #pileupIterator = bamfile.pileup(alignmentRef.id)
        # TODO: Check where this pileup iterator starts. Are there reads mapped before or after the reference begins/ends?
        #for pileupColumn in pileupIterator:
        
        """:
        #TODO: Detect heterozygosity here
            # Do the base frequencies look "normal"?
            # High proportion of Inserts or Deletions?
            # 
            
            # TODO: Params?            
            matchCutoff = .90
            heterozygousCutoff = .30
            
            isHeterozygousBase = False
            
            if (1.0 * matchCount / alignedCount) > matchCutoff:
                # Looks fine.
                pass
            elif (1.0 * mismatchCount / alignedCount) > heterozygousCutoff:
                isHeterozygousBase = True
            elif (1.0 * inCount / alignedCount) > heterozygousCutoff:
                isHeterozygousBase = True
            elif (1.0 * delCount / alignedCount) > heterozygousCutoff:
                isHeterozygousBase = True
            
            if isHeterozygousBase:
                heterozygousBasesSummaryFile.write(
                    'Position:' + str(referencePosition)
                    + '\nReferenceBase:' + str(referenceBase)
                    + '\nAlignedCount:' + str(alignedCount)
                    + '\nMatchCount:' + str(matchCount)
                    + '\nMismatchCount:' + str(mismatchCount)
                    + '\nInsertion Count:' + str(inCount) 
                    + '\nDeletion Count:' + str(delCount)
                    + '\nA Count:' + str(aCount)
                    + '\nG Count:' + str(gCount)
                    + '\nC Count:' + str(cCount)
                    + '\nT Count:' + str(tCount)
                    + '\n\n'
                    
                    
                )
        """
        # Iterate through reads
        # this read has an array of "polymorphisms" the lenght of our heterozyoug positions
        # #store a 1 for  polymorphic position 1, -1 for the other one.            
            # Iterate through positions

                # assign 1, 0, -1 based on 
                # match, notsure, mismatch/insert/del
                
        # iterate reads
            # if there are no clusters
                # this read starts a cluster
            
            # if there is only one cluster
                # Am I kinda like that cluster?
                    # add to cluster
                # else
                    # new cluster
                    
            # if there is 2 (or more?) clusters
                #iterate clusters
                    #Am I kinda like that cluster? and not liek others?
                        # add to cluster
                
        
    
    
    # Input = location of Reference Sequence
    # Output = Location of the Consensus Sequence
    def analysisRecurser(self, currentReferenceSequence):
        print('\n\nAttempting a read alignment and consensus polish, Iteration ' + str(self.currentIteration))
>>>>>>> referencemethod
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


            self.alignReads(currentReferenceSequence,self.readFileName,currentIterationSubdirectory, True)
            self.analyzeAlignment(currentIterationSubdirectory)
            
            # If we want more iterations, I should Recurse and try again.
            if (int(self.totalIterations) > int(currentIteration)):
                # On the next iteration, we want to use the new consensus sequence 
                # as the reference. 
                print('I am on iteration (' + str(currentIteration) + '/' + self.totalIterations + ') I will continue...')
                newReferenceSequenceFileName = join(currentIterationSubdirectory, 'Consensus.fasta')
                #currentIteration += 1
                
                # Return the consensus from one layer deeper.

                return self.analysisRecurser(newReferenceSequenceFileName)
                
                
            else:
                print('That was the last iteration (' + str(currentIteration) + '/' + self.totalIterations + '), I am done now.')
                
                # Return the consensus
                return join(currentIterationSubdirectory, 'Consensus.fasta')
                
            
            
            
        except Exception:
            print ('Exception encountered in analyzeReads()')                  
            print sys.exc_info()[0]
            print sys.exc_info()[1]
            print sys.exc_info()[2]        
            raise 
        
        
    # Perform BW Alignment.  Align all reads against the Reference.
    def alignReads(self, referenceLocation, readFileLocation, alignmentOutputDirectory, useReadSubset):
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
            indexCmd = BwaIndexCommandline(infile=newReferenceLocation, algorithm="bwtsw")
            indexCmd()

            
        except Exception:
            print ('Exception indexing alignment reference. Is bwa installed? folder writing permission issue?')                  
            raise
        
        # TODO: Make this a commandline parameter.  Lower = faster. Higher = more accurate consensus correction
        alignmentReadSubsetCount = 200
        try:
            if useReadSubset:
                # load Reads
                parsedReads = list(SeqIO.parse(readFileLocation, self.readInputFormat)) 
                
                # choose random subset
                randomIndexes = range(0, len(parsedReads))
                random.shuffle(randomIndexes)                
                sampledReads = []
                for i in range(0,alignmentReadSubsetCount):
                    sampledReads.append(parsedReads[randomIndexes[i]])
                
                # write random reads to alignment directory
                # Reassign the reads we'll use downstream
                readFileLocation = join(alignmentOutputDirectory, 'ReadSample.fasta')        
                readSampleWriter = createOutputFile(readFileLocation)          
                   
                SeqIO.write(sampledReads, readSampleWriter, 'fasta')
                readSampleWriter.close()

            else:
                # We'll use the whole read file.
                pass
                        
        except Exception:
            print ('Exception selecting a read subset.')                  
            raise
        
        # Part 2 Align
        try:
            # TODO: How can i put this into biopython?  Pipelines are hard.
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
        
        # Todo: I should keep a more structured array of info for these alignments.
        # Store this info into an object
        #class columnStats():
        
        # Keep a running total of adjustments made to the reference.
        # If this total is 0, then theoretically the consensus matches the alignment reference, and we're done.
        totalSequenceAdjustments = 0
        
        # Iterate the reference sequence column by column.
        pileupIterator = bamfile.pileup(alignmentRef.id)
        # TODO: Check where this pileup iterator starts. Are there reads mapped before or after the reference begins/ends?
        for pileupColumn in pileupIterator:
            
            #columnResults = None
           # columnResults.name='ll'
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
           
        # TODO This header is nonsense
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
        #TODO: I should make this a commandline parameter. More = MSA takes longer. Less = worse reference
        msaReadCount = 3
        
        print ('I choose ' + str(msaReadCount) + ' random reads.'
            + '\nThese are aligned to form a rough initial consensus sequence. Here:'
            + '\n' + join(self.outputRootDirectory,'Initial_Reference')
            + '\nPerforming ClustalO Multiple Sequence Alignment Now...')
        try:            
            # Load Reads from File

            parsedReads = list(SeqIO.parse(self.readFileName, self.readInputFormat))            
            referenceSequence = None

            
            # Reference Directory
            referenceDirectory = join(self.outputRootDirectory,'Initial_Reference')
            if not isdir(referenceDirectory):
                mkdir(referenceDirectory)
                        
            if (len(parsedReads) > msaReadCount):
                

                # Select a subset of reads for Multiple SequneceAlignment. Randomly, i guess.
                randomIndexes = range(0, len(parsedReads))
                random.shuffle(randomIndexes)                
                rawClustalReads = []
                for i in range(0,msaReadCount):
                    rawClustalReads.append(parsedReads[randomIndexes[i]])
              
                rawClustalReadsFilename = join(referenceDirectory, 'MSARaw.fasta')                
                rawClustalReadsFileWriter = createOutputFile(rawClustalReadsFilename)        
                SeqIO.write(rawClustalReads, rawClustalReadsFileWriter, 'fasta')
                rawClustalReadsFileWriter.close()
            
                #Perform Clustal MSA
                clustalOAlignmentOutputFileName = join(referenceDirectory, 'clustalOAlignment.fasta')
                clustalOCommandLine = ClustalOmegaCommandline(infile=rawClustalReadsFilename, outfile=clustalOAlignmentOutputFileName, verbose=True, auto=True, force=True, threads=int(self.numberThreads))
                clustalOCommandLine()                
        
                # Calculate consensus 
                # A dumb consensus has lots of ambiguous nucleotides.  We'll polish those out later.
                alignmentType = 'fasta'    
                alignmentObject = AlignIO.read(clustalOAlignmentOutputFileName, alignmentType)           
                alignmentSummaryInfo = AlignInfo.SummaryInfo(alignmentObject)                
                dumbConsensus = alignmentSummaryInfo.dumb_consensus(threshold=.5)
                
                referenceSequence = SeqRecord(Seq(str(dumbConsensus) , IUPAC.IUPACUnambiguousDNA),
                    id='Initial_Consensus',
                    description='Initial_Consensus')

                
            # Else
            else:
                # Select the first read, use it as the reference. It's something.
                referenceSequence = parsedReads[0]        
             
            #Write reference to file
            self.referenceSequenceFileName = join(referenceDirectory, 'FirstGuessReference.fasta')            
            firstGuessRefFileWriter = createOutputFile(self.referenceSequenceFileName)        
            SeqIO.write([referenceSequence], firstGuessRefFileWriter, 'fasta')

            firstGuessRefFileWriter.close()
            
            return referenceFileName
       
       
            print ('Done making initial consensus sequence.')
      
                                    
                                     
        except Exception:
            print ('Exception encountered in createFirstGuessReferenceFromReads()') 
            print sys.exc_info()[0]
            print sys.exc_info()[1]
            print sys.exc_info()[2] 
            raise    
            


# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput
