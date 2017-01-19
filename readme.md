
#Bam to Fastq Docker

## Goal: converting bam files to processable fastq 

### Linh Lam January 10, 2017

## Overview
This docker image generates fastq file(s) from bam input. It uses bamUtil (https://github.com/statgen/bamUtil) to identify correctly formatted bam then run samtools (http://samtools.sourceforge.net/) to convert the file into fastqs. The Samtools bam to fastq con- version will remove reads with unknown paired ids. Once the conversion is done, the docker also performs additional steps to remove duplicate sequence ids from fastq files. The repository is https://github.com/UCSC-Treehouse/b2f_docker

Docker pull command : docker pull linhvoyo/btfv2

Usage: docker run -v /path/to/work:/data -e input=filename.bam linhvoyo/btfv2

### Implementation

docker run

– Command to activate the docker image.

-v /path/to/work

– Place the input.bam into a new working directory. The docker command -v option will mount the folder into the docker container. User must give full path to current working directory.

-e input=filename.bam

– Provide the name of input.bam. Do not state full path of the input file ie. /path/to/input.bam.

#### Example run command:

input: /pod/home/linhvoyo/work/gerald_C2DBEACXX_3.bam

docker run command: docker run -v /pod/home/linhvoyo/work:/data \

-e input=gerald_C2DBEACXX_3.bam linhvoyo/btfv2

#### Expected Output:

Output files will be in the mounted directory containing input.bam. The docker image will create a log, results.txt from bamUtil analysis, and one or two fq.gz files from single-end or pair-end reads, respectively.

#### Example output from gerald_C2DBEACXX_3.bam

##### gerald_C2DBEACXX_3.R1.fq.gz

##### gerald_C2DBEACXX_3.R2.fq.gz

##### gerald_C2DBEACXX_3.results.txt

```
Wed Dec 14 07:33:45 UTC 2016 Bam validate - gerald_C2DBEACXX_3.bam

Number of records read = 340099754

 Number of valid records = 340099754

TotalReads(e6)  340.10

MappedReads(e6) 250.30

PairedReads(e6) 340.10

ProperPair(e6)  208.33

DuplicateReads(e6)  0.00

QCFailureReads(e6)  0.00

MappingRate(%)  73.60

PairedReads(%)  100.00

ProperPair(%)   61.25

DupRate(%)  0.00

QCFailRate(%)   0.00

TotalBases(e6)  34009.98

BasesInMappedReads(e6)  25030.33

Returning: 0 (SUCCESS)
```

Note the "MappingRate" column. If input bam has <100% "Mapping Rate", it would suggest that file contains both mapped and unmapped reads. This information would be useful for other downstream analyses e.g fusion detection.

##### gerald_C2DBEACXX_3.log

```

Wed Dec 14 08:08:43 UTC 2016 Samtools fastq conversion

gerald_C2DBEACXX_3.bam

Wed Dec 14 08:52:28 UTC 2016 Removing all duplicate read ids -

   gerald_C2DBEACXX_3.bam

Wed Dec 14 09:33:51 UTC 2016 Combining paired end reads -

   gerald_C2DBEACXX_3.bam

gerald_C2DBEACXX_3.bam: is paired end reads

Wed Dec 14 09:36:19 UTC 2016 Compressing file ./c.ns.bam_output/c.ns.R1.fq.

   perl_pairs_R1.fastq

Wed Dec 14 09:38:48 UTC 2016 Compressing file ./c.ns.bam_output/c.ns.R2.fq.

   perl_pairs_R2.fastq

Wed Dec 14 09:38:48 UTC 2016 gerald_C2DBEACXX_3.bam - Conversion done
```

The log file will show "input.bam - Conversion done" once sample has finished processing.

# Script contents were included in this document for historical purposes; we should remove them so they do not get out of sync with updated code

## Scripts

Dockerfile

```
build command:

docker build . -t linhvoyo/btfv2

push command:

docker push linhvoyo/btv2

FROM ubuntu:14.04

MAINTAINER linh lam

    
RUN apt-get update && apt-get install -y --force-yes --no-install-

   recommends \

    python \

    python-pip \

    python-dev \

    build-essential \

    zlib1g-dev \

    libcurl4-gnutls-dev \

    libssl-dev \

pigz

WORKDIR /root

ADD https://github.com/samtools/samtools/releases/download/1.3.1/samtools

   -1.3.1.tar.bz2 /root

RUN tar xvf /root/samtools-1.3.1.tar.bz2

RUN make -C /root/samtools-1.3.1/htslib-1.3.1/

WORKDIR /root/samtools-1.3.1

RUN ./configure --without-curses

RUN make

RUN make install

ADD ./bamUtil /root/bamUtil

ADD ./libStatGen /root/libStatGen

ADD ./scripts /root/scripts

WORKDIR /root/libStatGen

RUN make

RUN make install

WORKDIR /root/bamUtil

RUN make

RUN make install

ADD ./run.sh /root

WORKDIR /data

ENTRYPOINT ["/root/run.sh"]
```
#### splitlines.pl

http://bridgecrest.blogspot.com/2011/08/sort-fastq-file-by-sequence.html

```
#!/usr/bin/perl -w

#splitlines.pl

use strict;

my $count = 0;

while(my $line = <STDIN>){

  chomp($line);

  my @vals = split(/\t/, $line);

  print $vals[0]."\n";

  print $vals[1]."\n";

  print $vals[2]."\n";

  print $vals[3]."\n";

}
```
#### mergelines.pl

   http://bridgecrest.blogspot.com/2011/08/sort-fastq-file-by-sequence.html

```
 #!/usr/bin/perl -w

#mergelines.pl

use strict;

my $count = 0;

while(my $line = <STDIN>){

  chomp($line);

  print $line;

  $count = ($count + 1) % 4;

  if($count == 0){

    print "\n";

  }else{

    print "\t";

} }
```

#### fastqCombinePairedEnd.py

https://github.com/enormandeau/Scripts/blob/master/fastqCombinePairedEnd.py

```
#!/usr/bin/env python

"""Resynchronize 2 fastq or fastq.gz files (R1 and R2) after they have been

trimmed and cleaned

WARNING! This program assumes that the fastq file uses EXACTLY four lines

per sequence

Three output files are generated. The first two files contain the reads of

   the

pairs that match and the third contains the solitary reads.

Usage:

python fastqCombinePairedEnd.py input1 input2 separator

input1 = LEFT  fastq or fastq.gz file (R1)

input2 = RIGHT fastq or fastq.gz file (R2)

separator = character that separates the name of the read from the part

   that

describes if it goes on the left or right, usually with characters ’1’ or

’2’.  The separator is often a space, but could be another character. A

space is used by default.

"""

# Importing modules

import gzip

import sys

# Parsing user input

try:

  in1 = sys.argv[1]

in2 = sys.argv[2]

except:

  print(__doc__)

sys.exit(1)

try:

  separator = sys.argv[3]

  
except:

  separator = " "

# Defining classes

class Fastq(object):

  """Fastq object with name and sequence

"""

def __init__(self, name, seq, name2, qual):

  self.name = name

self.seq = seq

self.name2 = name2

self.qual = qual

def getShortname(self, separator):

  self.temp = self.name.split(separator)

del(self.temp[-1])

return separator.join(self.temp)

def write_to_file(self, handle):

  handle.write(self.name + "\n")

handle.write(self.seq + "\n")

handle.write(self.name2 + "\n")

handle.write(self.qual + "\n")

# Defining functions

def myopen(infile, mode="r"):

  if infile.endswith(".gz"):

  return gzip.open(infile, mode=mode)

else:

  return open(infile, mode=mode)

def fastq_parser(infile):

  """Takes a fastq file infile and returns a fastq object iterator

"""

with myopen(infile) as f:

  while True:

  name = f.readline().strip()

if not name:

  break

seq = f.readline().strip()

name2 = f.readline().strip()

qual = f.readline().strip()

yield Fastq(name, seq, name2, qual)

# Main

if __name__ == "__main__":

  seq1_dict = {}

seq2_dict = {}

seq1 = fastq_parser(in1)

seq2 = fastq_parser(in2)

s1_finished = False

s2_finished = False

if in1.endswith(’.gz’):

  outSuffix=’.fastq.gz’

else:

  outSuffix=’.fastq’

with myopen(in1 + "_pairs_R1" + outSuffix, "w") as out1:

  with myopen(in2 + "_pairs_R2" + outSuffix, "w") as out2:

  with myopen(in1 + "_singles" + outSuffix, "w") as out3:

  while not (s1_finished and s2_finished):

try:

  s1 = seq1.next()

except:

  s1_finished = True

try:

  s2 = seq2.next()

except:

  s2_finished = True

# Add new sequences to hashes

if not s1_finished:

  seq1_dict[s1.getShortname(separator)] = s1

if not s2_finished:

  seq2_dict[s2.getShortname(separator)] = s2

if not s1_finished and s1.getShortname(separator) in seq2_dict:

  seq1_dict[s1.getShortname(separator)].write_to_file(out1)

seq1_dict.pop(s1.getShortname(separator))

seq2_dict[s1.getShortname(separator)].write_to_file(out2)

seq2_dict.pop(s1.getShortname(separator))

if not s2_finished and s2.getShortname(separator) in seq1_dict:

  seq2_dict[s2.getShortname(separator)].write_to_file(out2)

seq2_dict.pop(s2.getShortname(separator))

seq1_dict[s2.getShortname(separator)].write_to_file(out1)

seq1_dict.pop(s2.getShortname(separator))

# Treat all unpaired reads

for r in seq1_dict.values():

  r.write_to_file(out3)

for r in seq2_dict.values():

  r.write_to_file(out3)
```
run.sh

```
#!/bin/bash

bamInput=$input

#making temp dir.

mkdir $bamInput’_’output

#bam validate

  
echo ‘date‘ "Bam validate" >> $bamInput.results.txt

bam validate --in $bamInput 2>>$bamInput.results.txt

printf "\n" >>$bamInput.results.txt

printf "%*s" $COLUMNS | tr " " "=" >>$bamInput.results.txt

printf "\n" >>$bamInput.results.txt

set -eu -o pipefail

#set -eu here because the bam validate part will return a non zero value

then crahs

#samtools

echo ‘date‘ "Samtools fastq conversion">>$bamInput.log

echo -e " \t\t\t\t $bamInput " >>$bamInput.log

samtools fastq -1 ./$bamInput’_’output/${bamInput/.bam}.R1.fq -2 ./

   $bamInput’_’output/${bamInput/.bam}.R2.fq $bamInput

#dedup

echo ‘date‘ "Removing duplicate read ids - $bamInput">>$bamInput.log

for j in {1..2}; do

cat ./$bamInput’_’output/${bamInput/.bam}.R${j}.fq | perl /root/scripts/

   mergelines.pl | sort  -k1,1 -t " "  --stable --parallel=10 -T ./ -S 10G

   | uniq | perl /root/scripts/splitlines.pl > ./$bamInput’_’output/${

   bamInput/.bam}.R${j}.fq.perl

done

#fastqCombinePairedEnd

a=./$bamInput’_’output/${bamInput/.bam}.R1.fq.perl

b=./$bamInput’_’output/${bamInput/.bam}.R2.fq.perl

size=$(wc -c < $b)

if [ $size -ge 10 ]; then

echo ‘date‘ "Combining paired end reads" >>$bamInput.log

echo -e "\t\t\t\t $bamInput: is paired end reads">>$bamInput.log

python /root/scripts/fastqCombinePairedEnd.py $a $b

else

echo -e "\t\t\t\t $bamInput: is single end reads">>$bamInput.log

mv ./$bamInput’_’output/${bamInput/.bam}.R1.fq.perl ./$bamInput’_’output/${

   bamInput/.bam}.R1.fq.perl_pairs_R1.fastq

fi

#pigz

for uncompressedFq in ./$bamInput’_’output/*[0-9].fastq;do echo ‘date‘ "

   Compressing file $uncompressedFq">>$bamInput.log; pigz $uncompressedFq;

   done

#rename

for compressedFq in ./$bamInput’_’output/*fastq.gz;do    mv $compressedFq $

   {compressedFq/.perl_pairs_R[0-9].fastq.gz}.gz; done

echo ‘date‘ "$bamInput - Conversion done"  >>$bamInput.log

#moving gz. files to work directory

for gz in ./$bamInput’_’output/*.gz; do mv $gz ./; done


#

rm -r $bamInput’_’output

#chown output files

finish() {

    # Fix ownership of output files

    uid=$(stat -c ’%u:%g’ /data)

    chown $uid /data/*${bamInput/.bam}.R[0-9].fq.gz

    chown $uid /data/$bamInput.log

    chown $uid /data/$bamInput.results.txt

}

trap finish EXIT
```
Example command using GNU Parallel (https://www.gnu.org/software/parallel):

parallel -j 2 ’docker run -v /data/btfdockerv2/work:/data -e input={} linhvoyo/btfv2’ ::: e.ns.bam f.ns.bam

 
