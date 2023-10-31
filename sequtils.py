"""
Collection of simple tools for dealing with biological data 

"""
## STDLIB
import os
import sys
import gzip
import string
import re
import random
import tempfile
import types
import shelve
import subprocess
import string
import glob
import csv
import urllib2
import StringIO
import itertools

from collections import namedtuple
from xml.etree import ElementTree

## CONTRIB
import lumpy
import pandas
from matplotlib import mlab
from scipy import stats
import scipy.optimize
import pysam
import cruzdb, cruzdb.models
import utils


NT_HASH = {'A':'00', 
           'a':'00',
           'C':'01',
           'c':'01',
           'T':'10',
           't':'10',
           'G':'11',
           'g':'11'}
    
NT_REV_HASH = {'00':'A',  
               '01':'C',
               '10':'T',
               '11':'G'}

NT_HASH_TRANSLATE_TABLE = string.maketrans('ACTGactg', '01230123')
NT_REV_HASH_TRANSLATE_TABLE = string.maketrans('0123', 'ACTG')


CHR2ACC_MAP={}


#
# Fasta I/O
####################################

def readFasta(file):
 """
 Read a fasta file into a sequence (or list of sequences) for consumtion.
 returns a list of tuples (annotation, seq) for each sequence in the fasta
 file.

 >>> seqs = readFasta('tests/data/sample.fa')
 >>> seqs
 """
 if type(file) == types.StringType:
   if os.path.splitext(file)[-1] == '.gz':
     IF = gzip.GzipFile(file)
   else:
     IF = open(file)
 else:
   IF = file
 seqs = []
 for seq in [l.split('\n') for l in  IF.read().split('>')]:
   if len(seq) > 1:
     seqs.append((seq[0].strip(), ''.join([s.strip() for s in seq[1:]])))
 return(seqs)

def fastaIter(file, seqSeperator=None):
  """
  iterate over the sequences in a fasta file returning yielding annotation,
  sequence pairs.  

  seqSeperator allows sequences to be split at the seperator- useful for paired read data
  """
  if os.path.splitext(file)[-1] == '.gz':
    IF = gzip.GzipFile(file)
  else:
    IF = open(file)
  seq =[]
  for line in IF:
    if line.strip() == '':
      continue
    elif line[0]=='>':
      if len(seq)>0 and header:
        if seqSeperator is not None:
          seq = ''.join(seq)
          seqs = seq.split(seqSeperator)
          seqs.insert(0, header)
          yield(tuple(seqs))
        else:
	  yield(header, ''.join(seq))
      header = line[1:].strip()
      seq = []
    else:
      seq.append(line.strip())
  if len(seq)>0 and header:
    if seqSeperator is not None:
      seq = ''.join(seq)
      seqs = seq.split(seqSeperator)
      seqs.insert(0, header)
      yield(tuple(seqs))
    else:
      yield(header, ''.join(seq))

def writeFasta(filehandle, seq, header, newLines=False):
   """
   write a fasta file
   """
   filehandle.write('>%s\n'%(header)) 
   if newLines:
     filehandle.write('\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)]))
   else:
     filehandle.write('%s\n'%(seq))

def saveToFasta(seqs, filename):
  """
  quick and direty script to save sequences to a file in fasta format
  """
  OF = open(filename, 'w')
  for i,s in enumerate(seqs):
    if type(s)==type(tuple()):
      name, seq = s
    else:
      name = 'seq-%i'%(i)
      seq = s
    writeFasta(OF, seq, name)
  OF.close()

def seqAsTmpFasta(s):
    """
    return a tmp fasta file containing s
    """
    tf =tempfile.mktemp('from-bioutils')+'.fa'
    saveToFasta(s, tf)
    return(tf)

#
# Simple Sequence Tools
####################################

def complement(seq, forceType=True):
  """returns the complement of a DNA or RNA sequence (IUPAC nomenclature
  supported)"""
  letterTlower = string.find(seq, 't')
  letterTupper = string.find(seq, 'T')
  letterUlower = string.find(seq, 'u')
  letterUupper = string.find(seq, 'U')
  if forceType:
    if (letterTlower >= 0 or letterTupper >= 0) and \
	   (letterUlower >= 0 or letterUupper >= 0):
      msg = 'Sequence should contain Tt\'s or Uu\'s, but not both.'
      raise ValueError, msg
  #mRNA mode
  if letterUlower >= 0 or letterUupper >= 0:
    t = string.maketrans('ACBDGHKMNSRUWVYXacbdghkmnsruwvyx',  # full IUPAC data
                         'UGVHCDMKNSYAWBRXugvhcdmknsyawbrx')
  #DNA mode
  else:
    t = string.maketrans('ACBDGHKMNSRTWVYXacbdghkmnsrtwvyx',  # full IUPAC data
                         'TGVHCDMKNSYAWBRXtgvhcdmknsyawbrx')
  return(seq.translate(t))

def reverse(seq):
  """reverse a sequence"""
  return seq[::-1]

def reverseComplement(seq): 
  """reverse complement a sequence"""
  return(reverse(complement(seq)))

def countOccurance(pattern, ref, edits=0):
  """
  count the number times pattern occurs in the reference allowing for N edits (including gaps)
  """
  counts = ref.count(pattern)
  for i in range(edits+1):
    for ev in variantIter(pattern, i):
      counts+=len(re.findall(ev, ref))
  return(counts)

def ntCount(seq, nts):
  return(numpy.sum([seq.count(nt) for nt in nts]))

def variantIter(seq, edits, nucs=['A','C','T','G'], mask=None, ins=True, dels=True, subs=True, previous=None):
  """
  iterate through all unique varients of a sequences that are up to N edits from initial sequence.

  mask is a boolean array or list of length len(seq) which positions are allowed to vary
  ins:  controls if insertional variants will be included
  dels: controls if deletion variants will be included

  previous -- any variants included in this set will be excluded from the output
  """
  if previous is None:
      previous = set()
  if seq not in previous:
    yield seq
    previous.add(seq)
  if edits>0:
      if mask is None:
          mask = [True]*len(seq)
      else:
          assert len(seq)==len(mask), "ERROR mask!=len(seq) // len([%s])!=len(%s)"%(str(mask), seq)
      if subs:
      ## render all posible subsitutions
          for pos,pm in enumerate(mask):
              if pm:
                  for nt in nucs:
                      varSeq = '%s%s%s'%(seq[:pos],nt,seq[pos+1:])
                      for v in variantIter(varSeq, edits-1, nucs, mask, ins, dels, subs, previous):
                        yield v
                      previous.add(varSeq)
      ## render all posible insertions
      if ins:
          for nt in nucs:
              for pos,pm in enumerate(mask):
                  if pm:
                      ## insertions pre-position 
                      varSeq = '%s%s%s'%(seq[:pos],nt,seq[pos:])
                      newMask = mask[:pos]+[True]+mask[pos:]
                      for v in variantIter(varSeq, edits-1, nucs, newMask, ins, dels, subs, previous):
                        yield v
                      previous.add(varSeq)
                      ## insertions post-position  
                      varSeq = '%s%s%s'%(seq[:pos+1],nt,seq[pos+1:])
                      newMask = mask[:pos]+[True]+mask[pos:]
                      for v in variantIter(varSeq, edits-1, nucs, newMask, ins, dels, subs, previous):
                        yield v
                      previous.add(varSeq)
      ## render all posible deletions
      if dels:
          for pos,pm in enumerate(mask):
              if pm:
                  varSeq = "%s%s"%(seq[:pos],seq[pos+1:])
                  newMask = mask[:pos]+mask[pos+1:]
                  for v in variantIter(varSeq, edits-1, nucs, newMask, ins, dels, subs, previous):
                    yield v
                  previous.add(varSeq)

def randomReadIter(ref, numReads, meanLength):
  """
  randomly generate sequences with lengths derived from a poisson distribution defined by meanLength
  """
  N = len(ref)
  for i in xrange(numReads):
    readStart = random.randint(0, N)
    length = numpy.random.poisson(25, 1)[0]
    yield(ref[readStart:readStart+length])



def summarizeGeneHits(db, chrom, start, end=None, strand=None):
    """ given a ucsc db name (e.g. hg19, mm10) a chromosome and a start and end summarize what's there 
    """
    if type(db)==types.StringType or type(db)==types.UnicodeType:
        ucscGenome = cruzdb.Genome(db, host=LOCALUCSC.get('host','genome-mysql.cse.ucsc.edu'), 
                                       password=LOCALUCSC.get('password',''), 
                                       user=LOCALUCSC.get('user', 'genome'))
    else:
        ucscGenome = db

    results = {'transcripts':[],
               'polyA':[],
               }
    if end is None:
        end = start+1
    ii = cruzdb.models.Interval(start, end, chrom=chrom)
    ## First RefGenes
    refGenes = ucscGenome.knearest('refGene', ii) 
    for gene in refGenes:
        if strand is not None and  gene.strand != strand:
            continue
        if not utils.positionsOverlap((gene.start, gene.end), (start, end)):
            continue
        features = set([s.upper() for s in gene.features(start, end)])
        overlappingExons = [i for i, (exonStart, exonEnd) in enumerate(gene.exons) if utils.positionsOverlap((start,end), (exonStart,exonEnd))]
        if strand == '-':
            totalExonCount = len(gene.exons)
            overlappingExons = [totalExonCount-i for i in overlappingExons]

        results['transcripts'].append( {
                                        'symbol': gene.name2, 
                                        'accession': gene.name,
                                        '3UTR': 'UTR3' in features,
                                        '5UTR': 'UTR5' in features,
                                        'exonic': ('EXON' in features) or ('CDS' in features),
                                        'coding': 'CDS' in features,
                                        'intronic':  'INTRON' in features,
                                        'hit_exons': overlappingExons,
                                        }
                                     )
    if ucscGenome.db=='hg19':
        polyA = ucscGenome.knearest('polyaDb', ii)
        for pa in polyA:
            if strand is not None and  pa.strand != strand:
                continue
            if not utils.positionsOverlap((pa.start, pa.end), (start, end)):
                continue
            results['poly A'].append((pa.name, pa.score))


    return(results)

def summarizeSNPs(db, chrom, start, end=None, strand=None):
    """ summaryize the snp information from UCSC """
    raise NotImplementedError
    if type(db)==types.StringType or type(db)==types.UnicodeType:
        ucscGenome = cruzdb.Genome(db, host=LOCALUCSC.get('host','genome-mysql.cse.ucsc.edu'), 
                                       password=LOCALUCSC.get('password',''), 
                                       user=LOCALUCSC.get('user', 'genome'))
    else:
        ucscGenome = db
    if ucscGenome.db =='hg19' or ucscGenome.db == 'mm10':
        snps = ucscGenome.knearest('snp137Common', ii)
        for snp in snps:
            if strand is not None and  snp.strand != strand:
                continue
            if not utils.positionsOverlap((snp.start, snp.end), (start, end)):
                continue
            results['snps'].append(
                                   {'name': snp.name,
                                    'validation': snp.valid,
                                    'class': snp.__getattribute__('class'),  ## this is a hack
                                    'moltype': snp.molType,
                                    'alleles': zip(snp.alleles, map(float, snp.alleleFreqs.split(',')[:-1]))
                                    }
                                   )



# Other Biology-Related I/O
####################################

class BedTableIter:
  """
  return an iterator that iterates over records from a bed file
  """
  def __init__(self, file):
    self.IF = open(file)
  def next(self):
    line = self.IF.readline()
    tokens = self.IF.readline()[:-1].split('\t')
    if len(tokens)<12:
      raise StopIteration
    else:
      r = utils.Bunch(chrom = tokens[0],
                chromStart = int(tokens[1]),
                chromEnd = int(tokens[2]),
                name = tokens[3],
                score = tokens[4],
                strand = tokens[5],
                thickStart = map(int, tokens[6].split(',')),
                thickEnd = map(int, tokens[7].split(',')),
                itemRGB = tokens[8],
                blockCount = tokens[9],
                blockSizes = map(int, filter(None, tokens[10].split(','))),
                blockStarts = map(int,filter(None,  tokens[11].split(','))))
      return(r)
  def __iter__(self):
    return(self)

def wiggleToNumpy(file):
  """
  load a wiggle file into a numpy array
  """
  dt = [('chr', 'S1'), ('start', numpy.int64), ('end', numpy.int64), ('score', numpy.float)]
  IF = open(file)
  IF.readline() # skip header
  data = numpy.loadtxt(IF, dtype=dt)
  return(data)

def wiggleWritter(file, data, name=None, description=None, visibility='full', color=None, 
                  itemRGB=None, useScore=1, group=None, priority=None, db=None, 
                  offset=None, url=None, htmlUrl=None, writeHeader=True, sort=True):
  """
  add data to a wiggle format, data is expected to be iterable of 4-tuples containing
  chrom, chromStart, chromEnd, dataValue
  
  See for more information:
  http://genome.ucsc.edu/goldenPath/help/wiggle.html
  """
  if type(file)==types.StringType:
     fileHandle = open(file, 'w')
  else: ##FIXME: add check for file type
    fileHandle = file
  definition ={'name':name,
               'description':description,
               'visibility':visibility ,
               'color':color,
               'itemRGB':itemRGB,
               'useScore':useScore,
               'group':group,
               'priority':priority,
               'db':db,
               'offset': offset,
               'url':url,
               'htmlUrl': htmlUrl}
  if writeHeader:
    fileHandle.write("track type=wiggle_0 ")
    fileHandle.write(' '.join(['%s=%s'%(key,value) for key,value in definition.items() if value is not None])+'\n')
  if sort:
    data.sort(key=lambda x: (int(x[0].replace('chr','')), x[1]))
  for chrom, chromStart, chromEnd, dataValue in data:
    fileHandle.write('%s\t%i\t%i\t%f\n'%(chrom, chromStart,chromEnd, dataValue))

def pslToDF(input, isAString=False):
    """
    parse PSL records
    """
    columns = ['matches',    # - Number of bases that match that aren't repeats
              'misMatches', # - Number of bases that don't match
              'repMatches', # - Number of bases that match but are part of repeats
              'nCount',     # - Number of 'N' bases
              'qNumInsert', # - Number of inserts in query
              'qBaseInsert',# - Number of bases inserted in query
              'tNumInsert', # - Number of inserts in target
              'tBaseInsert',# - Number of bases inserted in target
              'strand',     # - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
              'qName',      # - Query sequence name
              'qSize',      # - Query sequence size
              'qStart',     # - Alignment start position in query
              'qEnd',       # - Alignment end position in query
              'tName',      # - Target sequence name
              'tSize',      # - Target sequence size
              'tStart',     # - Alignment start position in target
              'tEnd',       # - Alignment end position in target
              'blockCount', # - Number of blocks in the alignment (a block contains no gaps)
              'blockSizes', # - Comma-separated list of sizes of each block
              'qStarts',    # - Comma-separated list of starting positions of each block in query
              'tStarts']    # - Comma-separated list of starting positions of each block in target 

    if isAString:
        input = StringIO.StringIO(input)
    try:
        df = pandas.read_table(input, names=columns, header=None)
    except StopIteration:
        df = pandas.DataFrame(None, columns=columns)
    return(df)

def embossSimpleAlignmentParser(file):
    """
    parse the simple aformat for emboss alignments.  Return a list of dictionaries with the following keys:
    {
     seq1: 
     seq2: 
     seq1Name:
     seq2Name:
     seq1_start:
     seq1_end:
     seq2_start:
     seq2_end:
     score:
     length:
     identity:
     similarity:
     gaps:
     matchstr:
    }
    """
    alignments = []
    inblock=False
    ## scan through the header
    for line in file:
        if line[0]=='#':
            continue
        else:
            break
    ## Now we read the blocks
    inheader=False
    alignmentRecord = {}
    while file:
        line = file.next()
        if line.strip()=='':
            continue
        if line[:5]=="#====" and not inheader:
            if len(alignmentRecord)!=0:
                alignments.append(alignmentRecord)
            inheader=True
            alignmentRecord = {}
            continue
        elif line[:5]=="#====" and inheader:
            inheader=False
            continue
        if inheader:
            if line[1:].strip() =='':
                continue
            key, value = map(string.strip, line[1:].split(':'))
            if key in ['Identity', 'Similarity', 'Gaps']:
                percentValue = value.split('(')[1].split('%')[0]
                alignmentRecord[key]=utils.convertStringToNumber(percentValue)
            else:
                try:
                    alignmentRecord[key]=utils.convertStringToNumber(value)
                except ValueError:
                    alignmentRecord[key]=value

        elif line[:5]!='#----':
            seqname, start, alnStr, end = line.split()
            alignmentRecord.update({'seq1name':seqname, 'seq1_start':int(start), 'seq1_end':int(end), 'seq1':alnStr})
            line = file.next()
            alignmentRecord['matchstr'] = line.strip()
            line = file.next()
            seqname, start, alnStr, end = line.split()
            alignmentRecord.update({'seq2name':seqname, 'seq2_start':int(start), 'seq2_end':int(end), 'seq2':alnStr})
            alignments.append(alignmentRecord)
            alignmentRecord={}
        else:
            return(alignments)

def embossAlign(aseq, bseq, cmd='water', parseoutput=True, aformat='fasta', **kwargs):
  """
  use an EMBOSS aligner (tested with water and needle) to align the sequences
  passed in as a tuple of sequences
  """
  tf =tempfile.mktemp('emboss-from-bioutils')
  if type(bseq) != type([]):
      bseq = [bseq]
  saveToFasta([aseq], tf+'_aseq.fa')
  saveToFasta(bseq, tf+'_bseq.fa')
  cmdTemplate = '%s -auto -asequence %s -bsequence %s -aformat %s -outfile %s  %s'
  cmd = cmdTemplate%(cmd,
                     tf+'_aseq.fa',
                     tf+'_bseq.fa',
                     aformat,
                     tf+'.out',
                     ' '.join(['-%s %s'%(k,v) for k,v in kwargs.items()]))
  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
  ret = os.waitpid(p.pid,0)
  if parseoutput and aformat=='fasta':
      alignedSeqs = readFasta(open(tf+'.out'))
      groupedSeqs = []
      for i in range(0, len(alignedSeqs),2):
          aseq1 = alignedSeqs[i]
          aseq2 = alignedSeqs[i+1]
          groupedSeqs.append((aseq1, aseq2))
      results=groupedSeqs
  elif parseoutput and aformat=='simple':
      results = embossSimpleAlignmentParser(open(tf+'.out'))
  else:
      results = open(tf+'.out').read()
  os.remove(tf+'_aseq.fa')
  os.remove(tf+'_bseq.fa')
  os.remove(tf+'.out')
  return(results)

def sigmaAlign(seqs):
  """
  use sigma align to align two simple DNA or protein sequences
  """
  tf =tempfile.mktemp('signma-from-bioutils')
  saveToFasta(seqs, tf)
  p = subprocess.Popen('sigma -F '+tf, stdout=subprocess.PIPE, shell=True)
  ret = os.waitpid(p.pid,0)
  os.remove(tf)
  alignedSeqs = readFasta(p.stdout)
  return(alignedSeqs)

def loadSeqsInSeaview(seqs):
  """
  seqs is a tuple of (name, value) pairs, value is a python string
  """
  tmpname = '/tmp/tmp.fa'
  OF = open(tmpname, 'w')
  OF = open(tmpname, 'w')
  for name,seq in seqs:
    writeFasta(OF, seq, name)
  OF.close()
  os.system('seaview %s &'%(tmpname))

class htWrapper:
    """
    quick and dirty wrapper around a numpy kmer array to have a .get method
    """
    def __init__(self, a,k):
        self.a = a
        self.k =k
    def get(self, seq):
        return(self.a[ntHash(seq)])

def quickAssemble(ht, seed,k, maxContig=1000, ignoreMaxedEdges=True, minCount=0, maxCount=200, stopAtMaxedEdges=False):
    """
    quick assemble hack
    """
    reason=None
    path =[(s, numpy.zeros((4,), dtype=[('nt', 'S1'), ('depth', int)])) for s in seed]
    nts = [(0,'A'),(1,'C'),(2,'T'),(3,'G')] ## cached list of nts
    done = False
    curNode = seed[1:]
    i = 0
    prevNodes = set()
    while not done and len(path)<maxContig:
        i+=1
        putativePaths = numpy.zeros((4,), dtype=[('nt', 'S1'), ('depth', int)])
        for i,nt in nts:
            nextEdge = '%s%s'%(curNode, nt)
            putativePaths[i]=(nt, ht.get(nextEdge))
        maxedPaths = (putativePaths['depth']>maxCount) 
        ignoreSel = maxedPaths | (putativePaths['depth']<minCount) 
        if numpy.any(maxedPaths) and stopAtMaxedEdges:
            done=True
            reason="Path-maxed"
            break
        putativePaths['depth'] = numpy.where(ignoreSel, 0, putativePaths['depth'])
        maxWeight = numpy.argmax(putativePaths['depth'])
        if putativePaths['depth'][maxWeight]==0:
            done = True
            reason = "END"
            break
        if sum(putativePaths['depth']==putativePaths['depth'][maxWeight])>1:
            done = True
            reason = "BRANCH"
            break
        prevNodes.add(curNode)
        curNode = '%s%s'%(curNode[1:], putativePaths[maxWeight]['nt'])
        if curNode in prevNodes:
            done=True
            reason = "CYCLE"
            break
        path.append((putativePaths[maxWeight]['nt'],putativePaths))
    contig = ''.join([p[0] for p in path])
    weights =  numpy.array([p[1] for p in path], dtype=[('nt', 'S1'), ('depth', int)])
    return(contig, weights, reason)

def removeSubsequences(seqs):
    """ remove seqs that are a complete subsequence of another """
    #FIXME: this is a slow and dumb algorithm for this
    done =False
    compressedSeqs = set(seqs)
    for qs in seqs:
        for ts in compressedSeqs:
            if ts.find(qs)>-1 and ts!=qs:
                compressedSeqs.remove(qs)
                break
    return(list(compressedSeqs))

def ntHash(seq, safe=True):
    """
    render an integer hash of a seq based on a simple 2-bit encoding
    """
    try:
        h =int(string.translate(seq[::-1], NT_HASH_TRANSLATE_TABLE),4)
    except ValueError:
        if safe:
            return(None)
        else:
            raise ValueError('%s can not be hashed')
    return(h)

def ntRevHash(hash, seqLength):
    """
    return the sequence encoded for my the ntHash
    """ 
    ##FIXME: this is still kinda slow
    s = bin(hash)[2:].rjust(seqLength*2, '0')
    seq = ''.join([NT_REV_HASH[s[i-2:i]] for i in range(len(s),0,-2)])
    #seq = string.translate(numpy.base_repr(hash, 4)[::-1].ljust(seqLength,'0'), NT_REV_HASH_TRANSLATE_TABLE)
    return(seq)

def cigar2aln(cigar, baseseq):
    """ return the base sequence with variations introduced into it as specified in the cigar string """
    raise NotImplementedError 
    if type(cigar)==types.StringType:
        raise NotImplementedError 
        #need code to convert cigar string to bam cigar list like in pysam
        cigarOps={'M', 0,  #Match to ref.
                  'I', 1,  #Insert from ref.
                  'D', 2,  #Delete from ref.
                  'N', 3,  #skipped region in ref
                  'S', 4,  #soft clipping (clipped region in seq)
                  'H', 5,  #hard clipping (clipped region not in seq)
                  'P', 6,  #padding (silent deletion from padded ref.)
                  '=', 7,  #sequence match
                  'X', 8,  # sequence mismatch
                  }
    newSeq= list(baseseq)
    idx=0
    for op, N in cigar:
        if op==0:
            idx+=N
        elif op==1:
            newSeq.insert(idx, '-')
            idx+=N
        elif op==2:
            for i in range(N):
                newSeq.pop(idx)
        elif op==3:
            for i in range(N):
                newSeq[idx+i]=''
    return(''.join(newSeq))

def simpleBTSamParser(sf):
    """ simplified samfile parser for bowtie2 outputs"""
    columns = ['read', 'flagsum', 'ref', 'refpos', 'mapquality',
               'cigar', 'matesref', 'matepos', 'fragsize', 'readseq',
               'quals','optional']
    samrecord = namedtuple('samrecord', columns)
    reader=csv.reader(open(sf), delimiter='\t')
    dataTypeDict={'A':str, 'i': int, 'f':float, 'Z': str,
                  # FIXME: implement these
                  #'H': <some function to decode hex arrays>,
                  #'B': <some function to decode integer arrays>,
                  }
    for data in reader:
        record = data[:len(columns)-1]
        optionalFieldsStr=data[len(columns)-1:]
        optionalFields ={}
        for tag,type,value in [block.split(':') for block in optionalFieldsStr]:
            if type in dataTypeDict:
                optionalFields[tag]=dataTypeDict[type](value)
            else:
                optionalFields[tag]='%s:%s'%(type,value)
        record.append(optionalFields)
        yield(samrecord(*record))

def blatPSL2rec(file):
    """
    load blat PSL file into a record array
    """
    if type(file) == types.StringType:
        ## this reads the column headers that are split across two lines
        IF = open(filename)
    else:
        IF = file
    head = [IF.next() for i in range(4)]
    IF.close()
    names = [l1.strip().replace('-','')+l2.strip().replace('-','') for l1,l2 in zip(head[2].split('\t'), head[3].split('\t'))]
    ## this reads data into the rec array
    blatHits = mlab.csv2rec(filename, skiprows=5, delimiter='\t', names=names, checkrows=100)
    return(blatHits)

def gtfIter(filename, groupFields=None, groupFieldDelimiter=';', groupKeyValueSeperator=' '):
    """ iterate through records in a gtf file.  If group fields are provided the group section will be also parsed"""
    recIter = csv.reader(open(filename), delimiter='\t')

    fields = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'group']
    if groupFields is not None:
        fields =  fields[:-1]+groupFields
    gtfRecord = namedtuple('gtf_record', fields)
    if groupFields is None:
        for rec in recIter:
            rec[3]=int(rec[3])
            rec[4]=int(rec[3])
            if rec[5]!='.':
                rec[5]=float(rec[5])
            else:
                rec[5]=float('nan')
            gtfRec = gtfRecord(*rec)
            yield gtfRec
    else:
        for rec in recIter:
            rec[3]=int(rec[3])
            rec[4]=int(rec[3])
            if rec[5]!='.':
                rec[5]=float(rec[5])
            else:
                rec[5]=float('nan')
            groupDataIter = (groupData.strip().split(groupKeyValueSeperator) for groupData in rec[8].split(groupFieldDelimiter) if groupData != '')
            groupDict = dict([(name,value.replace('"','')) for name,value in groupDataIter])
            gtfData = rec[:-1]+[groupDict.get(k, None) for k in groupFields]
            gtfRec = gtfRecord(*gtfData)
            yield gtfRec

def bowtie2(seqs, options, usepysam=False):
  """
  return alignments for the sequences stored in the list seqs.   Alignments are returned as dictionarys

  All bowtie arguments are required except for input/output parameters.

  Notice:  this is intended for "quick" searches of a few sequences.  
  """
  tfOut =tempfile.mktemp('_pysis.bioutils-bowtie.sam')
  tfIn = tempfile.mktemp('_pysis.bioutils-bowtie.fa')
  if type(seqs) != type([]):
      seqs = [seqs]
  saveToFasta(seqs, tfIn)
  cmd =['bowtie2']
  cmd.extend(['-S', tfOut])
  cmd.extend(['-f', '-U', tfIn])
  if not usepysam:
    cmd.extend(['--sam-no-hd', '--sam-no-sq'])
  cmd.extend(options)
  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  ret = os.waitpid(p.pid,0)
  sys.stdout.write('----- BOWTIE STDERR -----\n')
  sys.stdout.write(p.stderr.read())
  sys.stdout.write('-------------------------\n')
  sys.stdout.flush()
  if usepysam:
    hits = [h for h in pysam.Samfile(tfOut)]
  else:
    hits = list(simpleBTSamParser(tfOut))
  os.remove(tfIn)
  os.remove(tfOut)
  return(hits)
 
def getBlastHitsFromOutputFormat4(filehandle, refOrder):
    """
    return a numpy array with the hits for a single query blast result formated
    with option 4.  RefOrder is a list of references as ordered in the referenceDB
    used for blast
    They dtype is:
    [('chr', 'S5'), ('start', int8), ('end', int8), 
     ('matches', int8), ('mismatches', int8), 
     ('gaps', int8), ('ins', int8), ('matchstring', bool*len(aso))] 
    """
    headerInfo = ''
    for line in filehandle:
        if line[:3] == '1_0':
            break
        else:
            headerInfo+=line
    match = re.search('(?P<id>\S+)\s+(?P<start>\S+)\s+(?P<aso>\S+)\s+(?P<end>\S+)',line)
    if match.group('start')!='1':
        raise ValueError('No matches found to whole ASO')
    gappedAso = match.group('aso')
    asoStartCol, asoEndCol =match.span('aso')
    startStartCol = match.span('start')[0]
    endStartCol = match.span('end')[0]
    ## Count the number of lines remaining in the file
    filehandle, temp = itertools.tee(filehandle)
    numberOfAlignments = 0
    for line in temp: 
        if re.match('^[0-9]{1,2}', line):
            numberOfAlignments+=1
    ## Define the numpy datastructure
    dt = numpy.dtype([('chr', 'S5'), ('start', int), ('end', int), 
                      ('matches', int8), ('mismatches', int8), 
                      ('gap', int8), ('ins', int8), ('matchstring', 'S1', len(gappedAso))])
    def lineParser(line):
        chrIdx = int(line[:startStartCol])
        start =int(line[startStartCol:asoStartCol])
        matchString = line[asoStartCol:asoEndCol]
        end = int(line[asoEndCol:])
        chr = refOrder[int(chrIdx)]
        start = int(start)
        end = int(end)
        matches = 0
        mismatches = 0
        ins = 0
        gap = 0 
        for pos,i in enumerate(matchString):
            if i=='.':
                matches+=1
            else:
                aso = gappedAso[pos]
                if i=='-' and aso == '-':
                    pass
                elif i=='-':
                    ins+=1
                elif aso=='-':
                    gap+=1
                elif i!=aso:
                    mismatches+=1
                else:
                    print "warning"
        return(chr, start, end, matches, mismatches, gap, ins, fromstring(matchString, dtype='S1'))
    blastResults = zeros(numberOfAlignments, dtype=dt)
    for i,line in enumerate(filehandle):
        if re.match('^[0-9]{1,2}', line):
            blastResults[i] = lineParser(line)
    return(blastResults, gappedAso)

def pprintAlignedSeqs(s1,s2, stream=sys.stdout):
    """
    for two sequences s1 and s2 that are aligned (e.g. len(s1)==len(s2)) 
    """
    if type(s1)!=type('') and len(s1)==2:
        a1, s1=s1
    else:
        a1 ='s1'
    if type(s2)!=type('') and len(s2)==2:
        a2, s2=s2
    else:
        a2 = 's2'
    fmtStr = "{:>%i}{} {:<}\n"%(max(len(a1), len(a2), 12)+3)
    stream.write(fmtStr.format(a1,':',s1))
    matchString = fmtStr.format('',' ', ''.join(['|' if x==y else " " for x,y in zip(s1.upper(), s2.upper())]))
    matchCount = matchString.count('|')
    stream.write(matchString)
    stream.write(fmtStr.format(a2,':', s2))
    stream.write(fmtStr.format("%i/%i matches"%(matchCount, len(s1)), '',''))


class GEO:
    """ parse a GEO soft file into pandas dataframe"""
    def __init__(self, query):
        """ """
        pass
    def _proc_new_entity(self,line):
        """ processor for entity lines"""
        label, value = line[1:].split('=')
        if label in self.__dict__:
            self.__dict__[label][value] = self._proc_entity_attributes()
        else:
            self.__dict__[label] = {value : self._proc_entity_attributes()}
    def _proc_entity_attributes(self):
        """ iterate through enity attributes yield name, value pairs"""
        entityDict = {}
        last_pos = self.IF.tell()
        line = self.IF.readline()
        while line[0]=='!':
            name, value = line[1:].split('=')
            if name not in entityDict:
                entityDict[name]=value
            elif type(entityDict[name])==types.ListType:
                entityDict[name].append(value)   
            else:
                entityDict[name]=[entityDict[name], value]  
            last_pos = self.IF.tell()
            line = self.IF.readline()
        self.IF.seek(last_pos)
        self._proc_new_datatable()
    def _proc_new_datatable(self):
        """ processor for attribute lines"""
        ## First process the header
        last_pos = self.IF.tell()
        line = self.IF.readline()
        while line[0]=='#':
            name, value = line[1:].split('=')
            if name not in entityDict:
                entityDict[name]=value
            elif type(entityDict[name])==types.ListType:
                entityDict[name].append(value)   
            else:
                entityDict[name]=[entityDict[name], value]  
            last_pos = self.IF.tell()
            line = self.IF.readline()
        self.IF.seek(last_pos)
        self._proc_new_datatable()
    def _load(self):
        self.IF.seek(0)
        line = self.IF.readline()
        while line:
            if line[0]=='^':
                self._proc_new_entity(line)
            elif line[0]=='#':
                self._proc_new_datatable(line)
            else:
                raise ValueError, 'unexpected format'
            line = self.IF.readline()
