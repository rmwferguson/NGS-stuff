#!/bin/bash
#methods paper miseq analysis. 
# chmod +x /path/to/yourscript.sh 


# copy seq files into current directory 
cp /home/rob/methodsPaper/seqs/**/*.gz . 

# unzip
gunzip *.gz

# rename something use full to remove rubbish
#this includes removing indexes

rename 's/AAGAGGCA-//' *.fastq # forward indexes
rename 's/GTAGAGGA-//' *.fastq
rename 's/GCTCATGA-//' *.fastq
rename 's/AGGCAGAA-//' *.fastq
rename 's/TCCTGAGC-//' *.fastq
rename 's/GGACTCCT-//' *.fastq
rename 's/ATCTCAGG-//' *.fastq
rename 's/TAAGGCGA-//' *.fastq
rename 's/CGTACTAG-//' *.fastq
rename 's/TAGGCATG-//' *.fastq
rename 's/CTCTCTAC-//' *.fastq
rename 's/CGAGGCTG-//' *.fastq
rename 's/_CGTCTAAT_L002//' *.fastq #rev indexes
rename 's/_TCTCTCCG_L002//' *.fastq

# quality trimming with Sickle (add sickle to folder first)

mkdir qualTrimmed

#./sickle pe -f M_C_C2_R1.fastq -r M_C_C2_R2.fastq -o qualTrimmed/M_C_C2_R1_trimmed.fastq -p qualTrimmed/M_C_C2_R2_trimmed.fastq -s qualTrimmed/M_C_C2.fastq -t sanger -n

# for loop to automate

for file in *_R1.fastq; do ./sickle pe -f $file -r ${file%_R1.fastq}_R2.fastq -t sanger -o qualTrimmed/${file%.fastq}_trimmed.fastq -p qualTrimmed/${file%R1.fastq}R2_trimmed.fastq -s qualTrimmed/${file%_R1.fastq}_singles.fastq; done	

#error correction with SPAdes and BayesHammer

#mkdir qualTrimmed/errorCorrected

#~rob/spades/SPAdes-3.7.1-Linux/bin/spades.py -o qualTrimmed/errorCorrected --only-error-correction -1 qualTrimmed/M_C_C1_R1_trimmed.fastq -2 qualTrimmed/M_C_C1_R2_trimmed.fastq -t 8 -m 32 --disable-gzip-output

cd qualTrimmed

for file in *_R1_trimmed.fastq; do ~rob/spades/SPAdes-3.7.1-Linux/bin/spades.py -o ${file%_R1_trimmed.fastq}_errorCorrected --only-error-correction -1 $file -2 ${file%_R1_trimmed.fastq}_R2_trimmed.fastq -t 8 -m 32 --disable-gzip-output; done 

#house keeping, move files into a central folder and modify seq headers so they work downstream

mkdir allErrorCorrected

cd ..

cp qualTrimmed/**/corrected/*R[0-9]*.fastq qualTrimmed/allErrorCorrected

cd qualTrimmed/allErrorCorrected

for f in *R1*.fastq; do sed -i "/$(grep -m 1 "^@" $f | egrep -o "^[^:]+")/ s/$/ 1:N:0:/" $f; done 


for f in *R2*.fastq; do sed -i "/$(grep -m 1 "^@" $f | egrep -o "^[^:]+")/ s/$/ 2:N:0:/" $f; done

# Join reeds with PANDAseq (this removes primers which you might need so check this)

mkdir pairedSamples

#single set of seqs
#pandaseq -f sample1_R1_trimmed.00.0_0.cor.fastq -r sample1_R2_trimmed.00.0_0.cor.fastq -A pear -B -p CCTACGGGNGGCWGCAG -q GACTACHVGGGTATCTAATCC -w pairedSamples/sample1_paired.fasta -g pairedSamples/sample1_log.txt

#with loop

for f in *R1*.fastq; do pandaseq -f $f -r ${f%R1_trimmed.00.0_0.cor.fastq}R2_trimmed.00.0_0.cor.fastq -A pear -B -p CCTACGGGNGGCWGCAG -q GACTACHVGGGTATCTAATCC -w pairedSamples/${f%_R1_trimmed.00.0_0.cor.fastq}_paired.fasta -g pairedSamples/${f%_R1_trimmed.00.0_0.cor.fastq}_log.txt; done 

cd pairedSamples 

# adjust seq headers and then combine into one file

for f in *.fasta;
do sample=">$(echo $f | awk -F "." '{print $1}')_" 
sed -i "s/>/$sample/" $f
done

cat *.fasta > seqs.fna

cp seqs.fna /home/rob/methodsPaper/seqs/
cd /home/rob/methodsPaper/seqs/

now ready for Qimme analysis. 






