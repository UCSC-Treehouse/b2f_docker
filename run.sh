#!/bin/bash
bamInput=$input

#making temp dir.
mkdir $bamInput'_'output

#bam validate
echo `date` "Bam validate" >> $bamInput.results.txt
bam validate --in $bamInput 2>>$bamInput.results.txt
printf "\n" >>$bamInput.results.txt
printf "%*s" $COLUMNS | tr " " "=" >>$bamInput.results.txt
printf "\n" >>$bamInput.results.txt

set -eu -o pipefail
#set -eu here because the bam validate part will return a non zero value then crahs

#samtools
echo `date` "Samtools fastq conversion">>$bamInput.log
echo -e " \t\t\t\t $bamInput " >>$bamInput.log
samtools fastq -1 ./$bamInput'_'output/${bamInput/.bam}.R1.fq -2 ./$bamInput'_'output/${bamInput/.bam}.R2.fq $bamInput


#dedup
echo `date` "Removing duplicate read ids - $bamInput">>$bamInput.log
for j in {1..2}; do
cat ./$bamInput'_'output/${bamInput/.bam}.R${j}.fq | perl /root/scripts/mergelines.pl | sort  -k1,1 -t " "  --stable --parallel=10 -T ./ -S 10G | uniq | perl /root/scripts/splitlines.pl > ./$bamInput'_'output/${bamInput/.bam}.R${j}.fq.perl
done


#fastqCombinePairedEnd
a=./$bamInput'_'output/${bamInput/.bam}.R1.fq.perl
b=./$bamInput'_'output/${bamInput/.bam}.R2.fq.perl
size=$(wc -c < $b)
if [ $size -ge 10 ]; then
echo `date` "Combining paired end reads" >>$bamInput.log
echo -e "\t\t\t\t $bamInput: is paired end reads">>$bamInput.log
python /root/scripts/fastqCombinePairedEnd.py $a $b
else
echo -e "\t\t\t\t $bamInput: is single end reads">>$bamInput.log
mv ./$bamInput'_'output/${bamInput/.bam}.R1.fq.perl ./$bamInput'_'output/${bamInput/.bam}.R1.fq.perl_pairs_R1.fastq
fi

#pigz
for uncompressedFq in ./$bamInput'_'output/*[0-9].fastq;do echo `date` "Compressing file $uncompressedFq">>$bamInput.log; pigz $uncompressedFq; done

#rename
for compressedFq in ./$bamInput'_'output/*fastq.gz;do    mv $compressedFq ${compressedFq/.perl_pairs_R[0-9].fastq.gz}.gz; done

echo `date` "$bamInput - Conversion done"  >>$bamInput.log

#moving gz. files to work directory
for gz in ./$bamInput'_'output/*.gz; do mv $gz ./; done

#
rm -r $bamInput'_'output

#chown output files
finish() {
    # Fix ownership of output files
    uid=$(stat -c '%u:%g' /data)
    chown $uid /data/*${bamInput/.bam}.R[0-9].fq.gz
    chown $uid /data/$bamInput.log
    chown $uid /data/$bamInput.results.txt
}
trap finish EXIT
