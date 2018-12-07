# Summarize partitions of a genome assembly

Download most current assembly

```
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.24.fasta.gz
```

### Calculate the following for all sequences less than or equal to 100kb and all sequences greater than 100 kb

Calculating total number of nucleotides, total number of Ns, and total number of sequences.

**Sequences greater than 100000**
```
module load jje/jjeutils
module load jje/kent

bioawk -c fastx 'length($seq) > 100000{ print ">"$name; print $seq }'  dmel-all-chromosome-r6.24.fasta.gz > bigger.FASTA

faSize bigger.FASTA

137547960 bases (490385 N's 137057575 real 137057575 upper 0 lower) in 7 sequences in 1 files
Total size: mean 19649708.6 sd 12099037.5 min 1348131 (4) max 32079331 (3R) median 23542271
N count: mean 70055.0 sd 92459.2
U count: mean 19579653.6 sd 12138278.9
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real

```
Total number of nucleotides: 137547960, total number of N's: 490385, and total number of sequences: 7.


**Sequences less than or equal to 100000**

```
bioawk -c fastx 'length($seq) <= 100000{ print ">"$name; print $seq }'  dmel-all-chromosome-r6.24.fasta.gz > small.fasta

faSize small.fasta

6178042 bases (662593 N's 5515449 real 5515449 upper 0 lower) in 1863 sequences in 1 files
Total size: mean 3316.2 sd 7116.2 min 544 (211000022279089) max 88768 (Unmapped_Scaffold_8_D1580_D1567) median 1567
N count: mean 355.7 sd 1700.6
U count: mean 2960.5 sd 6351.5
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real
```
Total number of nucleotides: 6178042, total number of N's: 662593, and total number of sequences: 1863.


### Plots of the following for the whole genome, for all sequences less than or equal to 100kb and all sequences greater than 100 kb

Plot the sequence length distribution, the sequence GC% distribution, and the cumulative genome size sorted from largest to smallest sequences

**Whole Genome**

Sequence length distribution for whole genome
```
bioawk -c fastx '{ print $name, length($seq) }' dmel-all-chromosome-r6.24.fasta.gz > Length_wholegenome.txt
```

To make the plot in R
```
Length <-read.table("Length_wholegenome.txt", header=TRUE)

Length$Cut <-cut(x=Length$Length, breaks=10000)

Lengthplot <-ggplot(data=Length)

Lengthplot + geom_bar(mapping = aes(x=Cut)) + labs(title="Sequence Length Distribution", x="Length", y="Count (Number of Contigs)") + theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))
```

![seq_whole](https://github.com/tiffanybatarseh/EE282_HW4/blob/master/seqlength_whole_10000.png)

Sequence GC distribution for whole genome
```
bioawk -c fastx '{ print $name, gc($seq) }' dmel-all-chromosome-r6.24.fasta.gz > GC_wholegenome.txt
```

To make the plot in R
```
GC2 <-read.table("GC_wholegenome.txt", header=TRUE)

GC2$GC_Percentcut <-cut(x=GC2$GC_Percent, breaks = 11)

library(ggplot2)

GC <-ggplot(data=GC2)

GC + geom_bar(mapping = aes(x=GC2$GC_Percentcut)) + labs(title="Sequence GC Distribution", x="GC Percentage", y="Count (Number of Contigs)") + theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))

```

![GC_wholegenome](https://github.com/tiffanybatarseh/EE282_HW4/blob/master/GC_percent.png)

Cumulative genome size from largest to smallest sequences for whole genome

```
bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.fasta.gz \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nkblength_Ctg\t0" } { print "kblength_Ctg\t" $1 } ' \
>  seq_dmel_all.lengths

plotCDF2 seq_dmel_all.lengths seq_all.png
```

![seq_whole](https://github.com/tiffanybatarseh/EE282_HW4/blob/master/seq_all.png)

**Sequences > 100kb**

Sequence length distribution for sequences > 100kb

```
bioawk -c fastx '{ print $name, length($seq) }' bigger.FASTA > Length_bigger.txt
```

To make the plot in R
```
big <-read.table("Length_bigger.txt", header=TRUE)

big$cut <-cut(x=big$Length, breaks=50)

bigplot <-ggplot(data=big)

bigplot + geom_bar(mapping = aes(x=cut)) + labs(title="Sequence Length Distribution", x="Length", y="Count (Number of Contigs)") + theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))
```

![lengths_bigger](https://github.com/tiffanybatarseh/EE282_HW4/blob/master/Bigseqlength.png)

Sequence GC distribution for sequences > 100kb

```
bioawk -c fastx '{ print $name, gc($seq) }' bigger.FASTA > GC_bigger.txt
```

To make the plot in R
```
GCbig <-read.table("GC_bigger.txt", header=TRUE)

GCbig$PercentCut <-cut(x=GCbig$Percent, breaks = 5)

library(ggplot2)

GCbigplot <-ggplot(data=GCbig)

GCbigplot + geom_bar(mapping = aes(x=GCbig$PercentCut)) + labs(title="Sequence GC Distribution for Sequences >100kb", x="GC Percentage", y="Count (Number of Contigs)") + theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))
```

![GC_bigger](https://github.com/tiffanybatarseh/EE282_HW4/blob/master/GC_percent_bigger.png)

Cumulative genome size from largest to smallest sequences for sequences > 100kb

```
bioawk -c fastx ' { print length($seq) } ' bigger.FASTA \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nkblength_big_Ctg\t0" } { print "kblength_big_Ctg\t" $1 } ' \
>  seq_bigger2.lengths

plotCDF2 seq_bigger2.lengths seq_bigger2.png
```

![seq_bigger](https://github.com/tiffanybatarseh/EE282_HW4/blob/master/seq_bigger2.png)

**Sequences < or equal to 100kb**

Sequence length distribution for sequences < or equal to 100kb

```
bioawk -c fastx '{ print $name, length($seq) }' small.fasta > Length_small.txt
```

To make the plot in R

```
small <-read.table("Length_small.txt", header=TRUE)

small$cut <-cut(x=small$Length, breaks=100)

smallplot <-ggplot(data=small)

smallplot + geom_bar(mapping = aes(x=cut)) + labs(title="Sequence Length Distribution", x="Length", y="Count (Number of Contigs)") + theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))
```

![lengths_smaller](https://github.com/tiffanybatarseh/EE282_HW4/blob/master/Small_seqlength.png)

Sequence GC distribution for sequences < or equal to 100kb

```
bioawk -c fastx '{ print $name, gc($seq) }' small.fasta > GC_small.txt
```

To make the plot in R

```
GCsmall <-read.table("GC_small.txt", header=TRUE)

GCsmall$PercentCut <-cut(x=GCsmall$Percent, breaks = 11)

library(ggplot2)

GCsmallplot <-ggplot(data=GCsmall)

GCsmallplot + geom_bar(mapping = aes(x=GCsmall$PercentCut)) + labs(title="Sequence GC Distribution for Sequences <=100kb", x="GC Percentage", y="Count (Number of Contigs)") + theme_bw()+ theme(axis.text.x = element_text(angle = 60, hjust = 1))

```
![GC_wholegenome](https://github.com/tiffanybatarseh/EE282_HW4/blob/master/GC_percent_Small.png)

Cumulative genome size from largest to smallest sequences for sequences < or equal to 100kb

```
bioawk -c fastx ' { print length($seq) } ' small.fasta \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nkblength_small_Ctg\t0" } { print "kblength_small_Ctg\t" $1 } ' \
>  seq_small2.lengths

plotCDF2 seq_small2.lengths seq_small2.png
```

![seq_small](https://github.com/tiffanybatarseh/EE282_HW4/blob/master/seq_small2.png)

Sequence length distribution comparison

```
plotCDF2 seq_dmel_all.lengths seq_bigger2.lengths seq_small2.lengths seqlengths.png
```

![seq_together](https://github.com/tiffanybatarseh/EE282_HW4/blob/master/seqlengths.png)

# Genome Assembly

**Assemble a genome from MinION reads**

Submitted as a job.

```
#!/bin/bash
#$ -N ONP_hw4
#$ -q abio128,abio,bsg,bsg2
#$ -pe openmp 32
#$ -m beas

module load jje/jjeutils
minimap=$(which minimap)
miniasm=$(which miniasm)
basedir=/pub/jje/ee282/$USER
projname=nanopore_assembly
basedir=$basedir/$projname
raw=$basedir/$projname/data/raw
processed=$basedir/$projname/data/processed
figures=$basedir/$projname/output/figures
reports=$basedir/$projname/output/reports

createProject $projname $basedir
ln -sf /bio/share/solarese/hw4/rawdata/iso1_onp_a2_1kb.fastq $raw/reads.fq

$minimap -t 32 -Sw5 -L100 -m0 $raw/reads.fq{,} \
| gzip -1 \
> $processed/onp.paf.gz

$miniasm -f $raw/reads.fq $processed/onp.paf.gz \
> $processed/reads.gfa

awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' $processed/reads.gfa \
| fold -w 60 \
> $processed/unitigs.fa
```

**Calculate N50**

```
bioawk -c fastx ' {li=length($seq); l=li+l; print li; } END {print l;} ' unitigs.fa \
| sort -rn \
| gawk ' NR == 1 {l = $1;} NR >1 {li=$1; lc=li+lc; if(lc/l >= 0.5) {print li; exit;} }' \
| less -S
```
The N50 = 4494246

Compare the N50 to the Flybase's contig N50

Flybase assembly of ISO1 N50 = 21,485,538

ONT assembly of ISO1 N50 = 4,492,246

**Compare assembly to contig assembly from FlyBase using dotplot with MUMmer**

Using faSplitByN to change the assembly from the scaffold to contig level. Checking the N50 to make sure it split into contigs correctly.

```
module load jje/jjeutils perl
faSplitByN dmel-all-chromosome-r6.24.fasta.gz dmel-all-chromosome-r6.24.ctg.fa 10

bioawk -c fastx ' {li=length($seq); l=li+l; print li; } END {print l;} ' dmel-all-chromosome-r6.24.ctg.fa \
| sort -rn \
| gawk ' NR == 1 {l = $1;} NR >1 {li=$1; lc=li+lc; if(lc/l >= 0.5) {print li; exit;} }' \
| less -S
```
The contig N50 was 21,485,538 which matches the contig N50 on NCBI.

Submitted the MUMmer analysis as a job through the HPC

```
#!/bin/bash
#$ -N MUMmer_hw4
#$ -q abio128,abio,bsg,bsg2
#$ -pe openmp 32-128
#$ -R Y
#$ -m beas

###Loading of binaries via module load or PATH reassignment
source /pub/jje/ee282/bin/.qmbashrc
module load gnuplot

###Query and Reference Assignment. State my prefix for output filenames
REF="/pub/jje/ee282/tbatarse/hw4/dmel-all-chromosome-r6.24.ctg.fa"
PREFIX="flybase"
SGE_TASK_ID=1
QRY=$(ls /pub/jje/ee282/tbatarse/nanopore_assembly/nanopore_assembly/data/processed/unitigs.fa | head -n $SGE_TASK_ID | tail -n 1)
PREFIX=${PREFIX}_$(basename ${QRY} .fa)

nucmer -l 100 -c 125 -d 10 -banded -D 5 -prefix ${PREFIX} ${REF} ${QRY}
mummerplot --fat --layout --filter -p ${PREFIX} ${PREFIX}.delta \
  -R ${REF} -Q ${QRY} --postscript

```

**Compare assembly to both the contig assembly and scaffold assembly from FlyBase using contiguity plot**

```
module load rstudio/0.99.9.9
module load perl
module load jje/jjeutils/0.1a
module load jje/kent

bioawk -c fastx ' { print length($seq) } ' /pub/jje/ee282/tbatarse/nanopore_assembly/nanopore_assembly/data/processed/unitigs.fa \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nUnitigs_Ctg\t0" } { print "Unitigs_Ctg\t" $1 } ' \
> unitigs.sorted.text

bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.fasta.gz \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\ndmel_scaffold\t0" } { print "dmel_scaffold\t" $1 } ' \
> dmel_scaffold.sorted.text

bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.ctg.fa \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\ndmel_ctg\t0" } { print "dmel_ctg\t" $1 } ' \
> dmel_ctg.sorted.text

plotCDF2 *.sorted.text /dev/stdout \
| tee hw4_plotcdf2.png \
| display
```

**Calculate BUSCO scores of both assemblies and compare them**

To run BUSCO on the nanopore assembly of ISO1, I submitted it as a job on the HPC.

```
#!/bin/bash
#$ -N BUSCO_hw4
#$ -q abio128,abio,bsg,bsg2
#$ -pe openmp 32-128
#$ -R Y
#$ -m beas

module load augustus/3.2.1
module load blast/2.2.31 hmmer/3.1b2 boost/1.54.0
source /pub/jje/ee282/bin/.buscorc

INPUTTYPE="geno"
MYLIBDIR="/data/users/tbatarse/bin/busco/lineages/"
MYLIB="diptera_odb9"
OPTIONS="-l ${MYLIBDIR}${MYLIB}"
QRY="/pub/jje/ee282/tbatarse/nanopore_assembly/nanopore_assembly/data/processed/unitigs.fa"
###Please change this based on your qry file. I.e. .fasta or .fa or .gfa
MYEXT=".fa"

BUSCO.py -c 128 -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}

```
This outputs a folder with the results and a short summary text file that gives the overall results of the analysis:

```
C:0.5%[S:0.5%,D:0.0%],F:1.1%,M:98.4%,n:2799

       13      Complete BUSCOs (C)
       13      Complete and single-copy BUSCOs (S)
       0       Complete and duplicated BUSCOs (D)
       32      Fragmented BUSCOs (F)
       2754    Missing BUSCOs (M)
       2799    Total BUSCO groups searched
```

To run BUSCO on the Flybase assembly of ISO1, I submitted it as a job on the HPC.

```
#!/bin/bash
#$ -N BUSCO_flybase
#$ -q abio128,abio,bsg,bsg2
#$ -pe openmp 32-128
#$ -R Y
#$ -m beas

module load augustus/3.2.1
module load blast/2.2.31 hmmer/3.1b2 boost/1.54.0
source /pub/jje/ee282/bin/.buscorc

INPUTTYPE="geno"
MYLIBDIR="/data/users/tbatarse/bin/busco/lineages/"
MYLIB="diptera_odb9"
OPTIONS="-l ${MYLIBDIR}${MYLIB}"
QRY="/pub/jje/ee282/tbatarse/hw4/dmel-all-chromosome-r6.24.ctg.fa"
###Please change this based on your qry file. I.e. .fasta or .fa or .gfa
MYEXT=".fa"

BUSCO.py -c 128 -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}

```
This outputs a folder with the results and a short summary text file that gives the overall results of the analysis:

```
C:98.7%[S:98.2%,D:0.5%],F:0.8%,M:0.5%,n:2799

        2763    Complete BUSCOs (C)
        2749    Complete and single-copy BUSCOs (S)
        14      Complete and duplicated BUSCOs (D)
        21      Fragmented BUSCOs (F)
        15      Missing BUSCOs (M)
        2799    Total BUSCO groups searched
```

Overall, the Flybase assembly appears more complete do to the higher score for complete BUSCOs at 98.7% whereas the assembly with ONP data was at 0.5%.
