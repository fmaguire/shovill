CPU=4
MEM=4
MINQ=3
R1='R1.fastq.gz'
R2='R2.fastq.gz'
TRIM='ILLUMINACLIP:adapters.fasta:1:30:11 LEADING:$(MINQ) TRAILING:$(MINQ) MINLEN:30 TOPHRED33'

all: contigs.fa

# pilon correct 
contigs.fa: contigs.fasta contigs.bam
	pilon --genome contigs.fasta --frags contigs.bam --output pilon --threads $(CPU) --changes --mindepth 0.5 2>&1 | tee -a 80-pilon.log

# make bam (only uses 2 CPUs hard coded) samtools bit may be wrong /tmp/samtools.$$ in original
contigs.bam: contigs.fasta R1.trim.fq.gz R2.trim.fq.gz
	bwa mem -v 3 -x intractg -t $(CPU) contigs.fasta R1.trim.fq.gz R2.fq.trim.gz | samtools sort --threads 2 -m $(MEM)G --reference contigs.fasta -T /tmp/samtools -o contigs.bam 2>&1 | tee -a 70-bwa.log 

# generate assembly need to configure the kmers settings
contigs.fasta: flash.notCombined_1.fastq.gz flash.notCombined_2.fastq.gz flashextendedFrags.fastq.gz kmers.txt
	spades.py -1 flash.notCombined_1.fastq.gz -2 flash.notCombined_2.fastq.gz -s flash.extendedFrags.fastq.gz --only-assembler --threads $(CPU) --memory $(MEM) -o . --tmp-dir /tmp -k $(cat kmers.txt)

# merge reads
flash.notCombined_1.fastq.gz flash.notCombined_2.fastq.gz flashextendedFrags.fastq.gz:
	R1.cor.fq.gz R2.cor.fq.gz
	flash -d . -o flash -z -M 300 -t $(CPU) R1.cor.fq.gz R2.cor.fq.gz 2>&1 | tee -a 50-flash.log

# correct gsize to do
R1.cor.fq.gz R2.cor.fq.gz: R1.trim.fq.gz R2.trim.fq.gz genome_size.txt
	lighter -od . -r R1.trim.fq.gz -r R2.trim.fq.gz -K 32 $(cat genome_size.txt) -t $(CPU) -maxcor 1 2>&1 | tee -a 40-lighter.log

# trim and discard singletons
#  need to do singletons
R1.trim.fq.gz R2.trim.fq.gz: $(R1) $(R2)
	trimmomatic PE -threads $(CPU) -phred33 $(R1) $(R2) R1.trim.fq.gz /dev/null R2.trim.fq.gz /dev/null $(TRIM) 2>&1 | tee -a 30-trimmomatic.log

genome_size.txt kmers.txt: $(R1)
	perl read_stats.pl 



