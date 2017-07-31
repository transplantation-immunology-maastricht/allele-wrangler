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

# Version 1.0 

SoftwareVersion = "Allele-Wrangler Version 1.0"

import sys
import pysam
import os
#import random
#from os.path import split, join, isdir
#from os import mkdir
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
#from Bio.Align.Applications import ClustalwCommandline
#from Bio.Align.Applications import ClustalOmegaCommandline
#from Bio.Sequencing.Applications import BwaIndexCommandline
#from Bio.Align import AlignInfo
#from Bio import AlignIO

#from subprocess import Popen, PIPE, STDOUT
#from shutil import copy

class AlignmentInfo():  
    
    def __init__(self):
        self.sequence = ''
        self.alignmentInfo = []



class AlignmentColumn():   
    
    def __init__(self):

        self.referencePosition = 0
        self.referenceBase = ''
        self.referenceAdjustment = '?'
        self.alignedCount = 0
        self.unalignedCount = 0
        self.matchCount = 0
        self.mismatchCount = 0
        self.inCount = 0
        self.delCount = 0
        self.aCount = 0
        self.gCount = 0
        self.cCount = 0
        self.tCount = 0
            
       