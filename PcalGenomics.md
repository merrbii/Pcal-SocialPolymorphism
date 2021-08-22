# PcalGenomics
# Bioinformatic pipeline used to analyze genomic data of _Pogonomyrmex californicus_
###### By Mohammed Errbii "Simo"

This markdown document gives a detailed description of the bioinformatic pipeline used to analyze WGS of queens and RADseq of males of _Pogonomyrmex californicus_. Check out our [manuscript](https://www.biorxiv.org/content/10.1101/2021.03.21.436260v1.abstract) for further details on our research hypotheses and major findings.

### I- Mapping, variant calling and filtering
##### 1) Filter raw reads for paired reads with a BQ > 20 and a minimum length > 40bp using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

```bash

java -jar trimmomatic-0.38.jar PE \
reads_R1.fastq.gz \
reads_R2.fastq.gz \
reads_R1_paired.fastq.gz \
reads_R1_unpaired.fastq.gz \
reads_R2_paired.fastq.gz \
reads_R2_unpaired.fastq.gz \
SLIDINGWINDOW:4:20 \
MINLEN:40
```

##### 2) map paired reads with [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)

```bash

bwa mem \
-t 16 \
Pcal2.reference.fa \
reads_R1_paired.fastq.gz \
reads_R2_paired.fastq.gz > alignment.sam
```

##### 3) pre-processing of the alignments using [Picard](https://broadinstitute.github.io/picard/) and [SAMtools](http://samtools.sourceforge.net/)

```bash
# a- clean sam files
java -jar picard.jar CleanSam I=alignment.sam O=alignment.C.sam

# b- covert to bam format
samtools view -S -b alignment.C.sam > alignment.C.bam

# c- sort bam file (sort by default assume sorting by coordinates):
samtools sort alignment.C.bam > alignment.CS.bam

# d- fix mate information if needed
java -jar picard.jar FixMateInformation \
I=alignment.CS.bam \
O=alignment.CSF.bam \
ADD_MATE_CIGAR=true ASSUME_SORTED=true

# e- fix read group information
java -jar picard.jar AddOrReplaceReadGroups \
I=alignment.CSF.bam \
O=alignment.CSFR.bam \
RGID=id \
RGLB=library \
RGPL=illumina \
RGPU=lane \
RGSM=sample

# f- locates and tags duplicate reads in bam file
java -jar picard.jar MarkDuplicates \
I=alignment.CSFR.bam \
O=alignment.CSFRD.bam \
M=alignment.CSFRD.txt \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
ASSUME_SORTED=true

# g- index bam file
samtools index alignment.CSFRD.bam
```


##### 4) Variant calling using [GATK](https://gatk.broadinstitute.org/hc/en-us)

```bash

# a- generate gVCFs
gatk HaplotypeCaller \
-R Pcal2.reference.fa \
-I alignment.CSFRD.bam \
-O alignment.CSFRD.g.vcf.gz \
-ERC GVCF -ploidy 2

# b- combine gVCFs
gatk CombineGVCFs \
-R Pcal2.reference.fa \
-V input.list \ # this file contains the list of all g.vcf.gz files (one per sample)
-O all.g.vcf.gz

# c- genotype the combined gVCF
gatk GenotypeGVCFs \
-R Pcal2.reference.fa \
-V all.g.vcf.gz \
-O all.vcf.gz
```

##### 5) Apply hard filters using [GATK](https://gatk.broadinstitute.org/hc/en-us)

```bash
# a- it is better to visualize the variants attributes and decide for thresholds.
bcftools query -f '%INFO/DP\t%INFO/QD\t%INFO/FS\t%INFO/SOR\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\n' all.vcf.gz > info_density.tsv
# visualize each attribute in R (i.e. make density plots or histograms)

# b- get SNPs only
gatk SelectVariants -V all.vcf.gz --select-type-to-include SNP -O all_snps.vcf.gz

# c- apply hard filters:
gatk VariantFiltration \
-V all_snps.vcf.gz  \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -5.0" --filter-name "MQRankSum-5" \
-filter "ReadPosRankSum < -5.0" --filter-name "ReadPosRankSum-5" \
-O all_snps_filtered.vcf.gz

# d- extract variants that have passed the applied filters
gatk SelectVariants -V all_snps_filtered.vcf.gz -select 'vc.isNotFiltered()' -O all_snpsPASS.vcf.gz
```
##### 6) More filtering using [VCFtools](https://vcftools.github.io/index.html)

Similarly, generate and visualize annotations associated with each variant to decide for cutoffs. Use this [script](./scripts/get.QCmetrics.from.vcf.sh) to generate a set of metrics that you can inspect using [R](https://www.r-project.org/). Now use VCFtools to apply filters.

```bash
vcftools \
--gzvcf all_snpsPASS.vcf.gz \
--recode --recode-INFO-all \
--max-missing X \
--minDP Y \
--maxDP YY \
--mac Z \
--min-alleles 2 \
--max-alleles 2 \
--out all_snpsPASS_filtered

# replace X, Y, YY and Z with appropriate values based on the determined cutoffs. In our study: X=0.85, Y=3, YY=24 and Z=3.
```

### II- Population structure analyses

##### 1) Convert vcf to plink format using [VCFtools](https://vcftools.github.io/index.html) and [PLINK](https://www.cog-genomics.org/plink/)

```bash
# generate chromosomes map:
grep -v "^#" all_snpsPASS_filtered.vcf|cut -f1 | uniq| sort -V|awk '{print $0"\t"$0}' > chrom.map

# vcf >> plink
vcftools \
--vcf all_snpsPASS_filtered.vcf \
--plink \
--chrom-map chrom.map\
--out all_snpsPASS_filtered

# generate bed file
plink \
--file all_snpsPASS_filtered \
--make-bed \
--aec \
--out all_snpsPASS_filtered
```
##### 2) Linkage Disequilibrium (LD) pruning using [PLINK](https://www.cog-genomics.org/plink/)

```bash
# first generate a list of position to remove (or keep)
plink \
--bfile all_snpsPASS \
--aec \
--indep-pairwise 1 kb 1 0.2 \
--out all_snpsPASS.prune
# outputs two files: all_snpsPASS.prune.in (markers to keep) and all_snpsPASS.prune.out (markers to remove)

# then
plink \
--bfile all_snpsPASS \
--exclude all_snpsPASS.prune.out
--aec \
--out all_snpsPASS.LDpruned \
--make-bed
```

##### 3) Phasing using [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)

```bash
# extract positions with no missing data. Imputation of missing data is also an option. However in the absence of a reference panel, it will introduce more bias!
vcftools --gzvcf all_snpsPASS_filtered.vcf.gz \
--max-missing 1 \
--recode --recode-INFO-all \
--out all_snpsPASS_nomissing.vcf

# compress and index
bgzip all_snpsPASS_nomissing.vcf

tabix -p vcf all_snpsPASS_nomissing.vcf.gz

# get call in each scaffold individually
for i in {1..25}; do bcftools view -r scaffold_$i all_snpsPASS_nomissing.vcf.gz -Oz -o scaffold_$i.vcf.gz; done

# perform phasing (3 steps)
for i in {1..25}; do shapeit -check --input-vcf scaffold_$i.vcf.gz --output-log scaffold_$i.alignments;done # check the vcf files for missingness, polymorphism and so on.

for i in {1..25}; do shapeit -phase -V scaffold_$i.vcf.gz -O scaffold_$i.phased.haps.gz scaffold_$i.phased.samples --output-log scaffold_$i.main -T 10;done

for i in {1..25}; do shapeit -convert --input-haps scaffold_$i.phased.haps.gz scaffold_$i.phased.samples --output-vcf scaffold_$i.Onlyphased.vcf --output-log scaffold_$i.convert -T 10;done
```

##### 4) Population structure analyses

######a) PCA using [PLINK](https://www.cog-genomics.org/plink/)

```bash
# PCA
plink \
--bfile all_snpsPASS.LDpruned \
--aec \
--out all_snpsPASS.LDpruned \
--pca 4 # to limit the analysis to the first 4 eigenvectors
```

######b) admixture using [ADMIXTURE](http://dalexander.github.io/admixture/)

```bash
# ADMIXTURE
$ for K in `seq 2 4`; do admixture --cv all_snpsPASS.LDpruned.bed $K | tee log${K}.out ; done
# to inspect CV error values: grep -h CV log*.out
```

###### c) Relationships using [fineSTRUCTURE](https://people.maths.bris.ac.uk/~madjl/finestructure/)

A detailed [manual](https://people.maths.bris.ac.uk/~madjl/finestructure/manual.html) is available for more specific needs and further description.

```bash

#You will need an id file with three columns: <NAME> <POPULATION> <INCLUSION>
# It should look like: inclusion: 1 include 0 exclude a particular individual.

P102R Pleometrotic 1
H104W Haplometrotic 1
H108G Haplometrotic 1
H109G Haplometrotic 1
... ... ...

# Unzip the phased haplotype files: (note that we limited the analysis to the 20 largest scaffolds here)
for i in {1..20}; do zcat scaffold_$i.phased.haps.gz > ../../fs/scaffold_$i.phased.haps; done

# Generate phased genetic data in PHASE format:
for i in {1..20}; do perl -I /home/m/merrbii/programs/lib/perl5/share/perl/5.26.1 ~/programs/fs_4.1.1/impute2chromopainter.pl scaffold_$i.phased.haps scaffold_$i.phase; done

# Generate recombination rate: Genetic distances between each SNP were calculated
# In our case, we assumed a uniform recombination rate, based on the genome‐wide estimate of 14 cM/mb (0.14 cM/bp) obtained by Sirviö et al. (2011) from AFLP markers in Pogonomyrmex rugosus.
# !!! edit the perl script makeuniformrecfile.pl !!!


for i in {1..20}; do perl -I /home/m/merrbii/programs/lib/perl5/share/perl/5.26.1 ~/programs/fs_4.1.1/makeuniformrecfile.pl scaffold_$i.phase scaffold_$i.rec ; done


# Then run fs

fs HaploPleo.20scf.cp \
-idfile idfile.txt \
-phasefiles scaffold_1.phase \
-phasefiles scaffold_2.phase \
-phasefiles scaffold_3.phase \
-phasefiles scaffold_4.phase \
-phasefiles scaffold_5.phase \
-phasefiles scaffold_6.phase \
-phasefiles scaffold_7.phase \
-phasefiles scaffold_8.phase \
-phasefiles scaffold_9.phase \
-phasefiles scaffold_10.phase \
-phasefiles scaffold_11.phase \
-phasefiles scaffold_12.phase \
-phasefiles scaffold_13.phase \
-phasefiles scaffold_14.phase \
-phasefiles scaffold_15.phase \
-phasefiles scaffold_16.phase \
-phasefiles scaffold_17.phase \
-phasefiles scaffold_18.phase \
-phasefiles scaffold_19.phase \
-phasefiles scaffold_20.phase \
-recombfiles scaffold_1.rec \
-recombfiles scaffold_2.rec \
-recombfiles scaffold_3.rec \
-recombfiles scaffold_4.rec \
-recombfiles scaffold_5.rec \
-recombfiles scaffold_6.rec \
-recombfiles scaffold_7.rec \
-recombfiles scaffold_8.rec \
-recombfiles scaffold_9.rec \
-recombfiles scaffold_10.rec \
-recombfiles scaffold_11.rec \
-recombfiles scaffold_12.rec \
-recombfiles scaffold_13.rec \
-recombfiles scaffold_14.rec \
-recombfiles scaffold_15.rec \
-recombfiles scaffold_16.rec \
-recombfiles scaffold_17.rec \
-recombfiles scaffold_18.rec \
-recombfiles scaffold_19.rec \
-recombfiles scaffold_20.rec \
-go

# The analysis should not take too long time. After 1-2 min you should see
# FineStructure completed successfully.
# Finestructure complete!
# Get started by running the GUI with ...
```
Finally, the fineSTRUCTURE output can then be processed and visualized in R using the following [finestructureR scripts](https://people.maths.bris.ac.uk/~madjl/finestructure/finestructureR.html) available from the tool's developpers.
### III- Demographic population history analysis

##### 2) Demographic population history analysis [MSMC2](https://github.com/stschiff/msmc2)
To run MSMC2, check ([Schiffels & Wang, 2020](https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_7)) for a very nice and detailed description of how to estimate effective population size (_N_<sub>e</sub>) and the relative cross-coalescence rate (rCCR). The [GitHub](https://github.com/stschiff/msmc2) repository has also some nice resources.

### IV- Population genomic analyses

##### 1) Calculating nucleotide diversity (&#960;), Tajima's _D_, heterozygosity (and inbreeding) and linkage disequilibrium (LD) from WGS of queens using [VCFtools](https://vcftools.github.io/index.html)

For the first four metrics, use this [script](./scripts/) to calculate these metrics. You will need the filtered vcf file generated after filtering the raw variant calls (i.e. `all_snpsPASS_filtered.vcf`), a text file listing individuals belonging to each population, and finally to specify the window/step size. Then use R (or your favorite software) to make plots and do the stats.

For LD, minor allele frequency can affect LD estimates (see for example [Otyama et al. 2019](https://doi.org/10.1186/s12864-019-5824-9)). To minimize the effects of rare variants, only SNPs with a minor allele frequency (MAF) of 0.2 are considered for LD calculation here.

```bash
# get SNPs with a maf >= 0.2
vcftools --gzvcf all_snpsPASS_filtered.vcf.gz --recode --recode-INFO-all --maf 0.2 --out all_snpsPASS_filtered.maf0.2

# then calculate pairwise LD
#for population1: Haplometrotic
vcftools --gzvcf all_snpsPASS_filtered.maf0.2.vcf.gz --geno-r2 --ld-window-bp 200000 --min-r2 0.001 --out all_snpsPASS_filtered.maf0.2.LD.haplo --keep vcfs/haplo.txt
# and population2: Pleometrotic
vcftools --gzvcf all_snpsPASS_filtered.maf0.2.vcf.gz --geno-r2 --ld-window-bp 200000 --min-r2 0.001 --out all_snpsPASS_filtered.maf0.2.LD.pleo --keep vcfs/pleo.txt

# check https://vcftools.github.io/man_latest.html for information on each of the parameters used above!
```

The output files would still need to be processed in R before plotting. The goal is to obtain LD measures averaged for markers within 100 bp (or more) in each population in order to make plotting easier.

```R
# load libraries
library(tidyverse)
library(readr)

# read pairwise ld for each population
h <- read_delim("all_snpsPASS_filtered.maf0.2.LD.haplo.geno.ld", delim = "\t")
p <- read_delim("all_snpsPASS_filtered.maf0.2.LD.pleo.geno.ld", delim = "\t")

# get average r2 per 100 bp
h$dist <- ceiling((h$POS2 - h$POS1)/100)*100 #round up to the nearest 100 bp
p$dist <- ceiling((p$POS2 - p$POS1)/100)*100 #round up to the nearest 100 bp

h2 <- group_by(h,CHR, dist) %>%
+   summarise(meanR2 = mean(`R^2`))
p2 <- group_by(p,CHR, dist) %>%
+   summarise(meanR2 = mean(`R^2`))

# store data
dt <- bind_rows(h2,p2) #Combine data for two populations
dt$population <- rep(c("Haplometrotic","Pleometrotic"),c(nrow(h2),nrow(p2))) #Adds a population column
write_delim(dt, "maf0.2.windowed.ld.forR.tsv", delim = "\t") #save as a tsv file


# produce plots

# read data
ld0.2 <- read_delim("maf0.2.windowed.ld.forR.tsv", delim = "\t")
ld0.2$CHROM<- as.numeric(gsub("\\D", "", ld0.2$CHR)) # keep only scaffold number
sub0.2 <- ld0.2%>%filter(ld0.2$CHROM <= 25) # only 25 scaffolds here!

# Plot LD between markers up to 10 kb, averaged across the genome for each population

decay <- ggplot(sub0.2 %>% filter(sub0.2$dist <= 10000),aes(x = dist/1000, y = meanR2, col = pop)) +
  geom_smooth(aes(color=pop, fill = pop), method = "loess", level=0.95, span= .1, se=T, size =.6) +
  scale_fill_manual(values = c("#0075DC", "#2BCE48"),name = "Population")+
  scale_color_manual(values = c("#0075DC", "#2BCE48"), name = "Population")+
  theme_bw()+ #xlim(c(0,20))+
  xlab("Distance (Kb)")+
  ylab(expression(r^2))+
  theme(
    legend.position = c(0.9, 0.6),
    legend.title = element_blank(),
    plot.title = element_text(size=10, face = "bold"),
    legend.text = element_text(size = 8),
    axis.text.y=element_text(size=8),
    axis.text.x=element_text(size=7),
    axis.title.y=element_text(size=8, angle = 0, vjust = 0.5),
    axis.title.x=element_text(size=8),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
  ) + ggtitle("")


decay


# Plot decay for each scaffold independently

p <- ggplot(sub0.2%>%filter(sub0.2$dist <= 10000),aes(x = dist/1000, y = meanR2, col = pop)) +
  geom_smooth(aes(color=pop, fill = pop), method = "loess", level=0.95, span= .1, se=T, size =.6) +
  scale_fill_manual(values = c("#0075DC", "#2BCE48"),name = "Population")+
  scale_color_manual(values = c("#0075DC", "#2BCE48"), name = "Population")+
  theme_bw()+ #xlim(c(0,20))+
  xlab("Distance (Kb)")+
  ylab(expression(r^2))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
  )+
  facet_wrap(~CHROM, nrow = 5, scales = "free_x")+
  theme(panel.spacing.x=unit(.5,"lines"),
        panel.spacing.y=unit(.5,"lines"),
        panel.background = element_blank(),
        axis.text.x = element_text(size=7, angle=45),
        axis.title.y=element_text(size=8, angle = 0, vjust = 0.5),
        #axis.ticks.x = element_blank(),
        panel.spacing = unit(.1, "cm"),
        #panel.background=element_rect(fill='white', colour='black'),
        strip.background=element_rect(fill='white', colour='black')
  )

p


```


##### 2) Computing genetic differentiation (_F_<sub>ST</sub>) using [VCFtools](https://vcftools.github.io/index.html)

Here, since highly related individuals can lead to inaccurate _F_<sub>ST</sub> estimates between popiulations, it is highly advised to exclude individuals from pairs of high relatedness before calculating _F_<sub>ST</sub>. You can do so using [PLINK](https://www.cog-genomics.org/plink/) as follow:

```bash
# first generate a relationship matrix
plink --bfile all_snpsPASS_filtered --make-rel square --aec --out all_snpsPASS_filtered.rel
# open/visualize the matrix to decide for a threshold then


# exclude one individual per pair of samples with relatedness greater than the given cutoff (here 0.4)
plink --bfile all_snpsPASS_filtered --rel-cutoff 0.4 --aec --out relatedness.0.4.cutoff

# outputs a list of individuals with relatedness lower than the given cutoff
# To be used when performing genome scans for selection.
```
Now that we have generated a list of unrelated individuals, we can proceed with the _F_<sub>ST</sub> calculation
```bash
vcftools --gzvcf all_snpsPASS_filtered.vcf.gz --fst-window-size 100000 	--fst-window-step 100000 --weir-fst-pop unrelated.haplo.txt --weir-fst-pop unrelated.pleo.txt --out all_snpsPASS_filtered.w100ks100k.unrel.fst
```
##### 3) Extra filtering

Depending on individual project and question needs, one can check the distribution of each of the calculated metrics (and others like the number of variant per window, coverage and mapping quality) to see if there are any extremes that need to be considered. For example regions with high repeat content (high coverage and/or low mapping quality) can give rise to windows with extremely high number of variants, which could lead to some bias. Therefore, it might be wise to check for this and eventually to filter any extremes. To do so, first estimate average mapping quality, coverage and SNP density per window. You will need [BEDtools](https://bedtools.readthedocs.io/en/latest/), [SAMtools](http://samtools.sourceforge.net/) and [VCFtools](https://vcftools.github.io/index.html).

```bash
#######
# MQ #
#######

## get chromosome size
awk -v OFS='\t' {'print $1,$2'} Pcal2.reference.fa.fai > Pcal2.genome.txt

## make windows
bedtools makewindows -g Pcal2.genome.txt -w 100000 -s 100000 > Pcal2.100kbwindows.bed

## calculate mean MQ per window from each alignment file. Here I use parallel to speed up the process
ls *.bam|parallel --jobs 3 'bedtools map -a Cobs2.1.100kbwindows.bed -b <(bedtools bamtobed -i {}) -c 5 -o mean > {.}.mapping' # {} will be replaced by an alignment file to be processed. See (https://www.gnu.org/software/parallel/parallel_tutorial.html) for more information on parallel.

## prepare for R
for i in *.mapping; do sed -i "s/scaffold//g" $i ;done # keep only scaffold numbers
for i in *.mapping; do base=$(echo $i|cut -f1 -d.); awk 'BEGIN {OFS=FS="\t"} $4 == "." {$4="NA"}1' $i > $base.forR;done #replace missing (.) with NA


#############
# Coverage #
#############

## calculate mean coverage per window for each alignment file
ls *.bam|parallel --jobs 30 'samtools bedcov Cobs2.1.100kbwindows.bed {} > {.}.bedcov' # same, {} will be replaced by a file to be processed from ls *.bam

## prepare for R
for i in *.bedcov; do base=$(echo $i|cut -f1 -d.); awk 'OFS=FS="\t" {print $1,$2,$3,$4,$4/100000}' $i > $base.forR;done


################
# SNP density #
################

# normally you don't need this as VCFtools would calculate this automatically while computing FST for example.

vcftools --gzvcf all_snpsPASS_filtered.vcf.gz --SNPdensity 100000 --out all_snpsPASS_filtered.SNPden100k

```

Now in R, combine all these metrics (i.e. &#960;, Tajima's _D_, _F_<sub>ST</sub>, MQ and mean coverage) and then visualize the distribution of each parameter to decide whether to exclude some windows or not. This will output a file that gives for each genomic window, its &#960;, Tajima's _D_, _F_<sub>ST</sub>, MQ and mean coverage estimates.
```R

# load libraries
library(tidyverse)
library(Rmisc)

# Load the original files

# Pi (for each population; i.e. haplo- and pleometrotic)
pi_haplo <- read_delim("haplo.pi.windowed.pi",delim = "\t")
pi_haplo$CHROM <- as.numeric(gsub("\\D", "", pi_haplo$CHROM))

pi_pleo <- read_delim("pleo.pi.windowed.pi",delim = "\t")
pi_pleo$CHROM <- as.numeric(gsub("\\D", "", pi_pleo$CHROM))


# tajimaD
taj_haplo <- read_delim("haplo.tajD.Tajima.D",delim = "\t")
taj_haplo$CHROM <- as.numeric(gsub("\\D", "", taj_haplo$CHROM))

taj_pleo  <- read_delim("pleo.tajD.Tajima.D",delim = "\t")
taj_pleo$CHROM <- as.numeric(gsub("\\D", "", taj_pleo$CHROM))


# FST
fst <- read_delim("unrel.fst.windowed.weir.fst", delim = "\t")
colnames(fst) <- c("CHROM", "BIN_START","BIN_END","N_VARIANTS", "WEIGHTED_FST", "MEAN_FST")
fst$CHROM <- as.numeric(gsub("\\D", "", fst$CHROM))


# Mapping quality
# add them to the list

path = "~/path/to/folder/containing/mappingQ.files"
qual<-list()
file.names <- dir(path, pattern ="^.*.mappingQ.forR$")

for(i in file.names){
  qual[[i]]<-read_delim(paste(path,i,sep=""), col_names=F, delim = "\t")
}

ci.qual <- as.data.frame(aaply(laply(qual, as.matrix), c(2, 3), CI))

colnames(ci.qual)<- c("X1.upper", "X2.upper", "X3.upper", "X4.upper", "CHROM", "BIN_START", "BIN_END",  "meanQual", "X1.lower", "X2.lower", "X3.lower", "X4.lower")


# Coverage

path = "~/path/to/folder/containing/coverage.files/"
cov<-list()
file.names <- dir(path, pattern ="^.*.cov.forR$")

for(i in file.names){
  cov[[i]]<-read_delim(paste(path,i,sep=""), col_names=F, delim = "\t")
}

ci.cov <- as.data.frame(aaply(laply(cov, as.matrix), c(2, 3), CI))

colnames(ci.cov)<- c("X1.upper", "X2.upper", "X3.upper", "X4.upper", "CHROM", "BIN_START", "BIN_END",  "meanCov", "X1.lower", "X2.lower", "X3.lower", "X4.lower")



################################################
# combine everything                          #
# from now on, this one file produced here can #
# be used to make plots and to do the stats    #
################################################

d <- left_join(dplyr::select(fst,c(CHROM, BIN_START, BIN_END, N_VARIANTS, weightedFst = WEIGHTED_FST,  meanFst= MEAN_FST)),
               dplyr::select(pi_haplo,c(CHROM, BIN_START, BIN_END,  piHap = PI)))

d1 <- left_join(dplyr::select(d,c(CHROM, BIN_START, BIN_END, N_VARIANTS, weightedFst,  meanFst, piHap)),
                dplyr::select(pi_pleo,c(CHROM, BIN_START, BIN_END,  piPle = PI)))


d1$BIN_START <- d1$BIN_START-1


d2 <- left_join(dplyr::select(d1,c(CHROM, BIN_START, BIN_END, N_VARIANTS, weightedFst,  meanFst, piHap, piPle)),
                dplyr::select(taj_haplo,c(CHROM, BIN_START,  tajDHap = TajimaD)))

d3 <- left_join(dplyr::select(d2,c(CHROM, BIN_START, BIN_END, N_VARIANTS, weightedFst,  meanFst, piHap, piPle, tajDHap)),
                dplyr::select(taj_pleo,c(CHROM, BIN_START,  tajDPle = TajimaD)))


d4 <- left_join(dplyr::select(d3,c(CHROM, BIN_START, BIN_END, N_VARIANTS, weightedFst,  meanFst, piHap, piPle, tajDHap, tajDPle)),
                 dplyr::select(ci.qual,c(CHROM, BIN_START, BIN_END,  AverageMQ = meanQual)))

d5 <- left_join(dplyr::select(d4,c(CHROM, BIN_START, BIN_END, N_VARIANTS, weightedFst,  meanFst, piHap, piPle, tajDHap, tajDPle, AverageMQ)),
                dplyr::select(ci.cov,c(CHROM, BIN_START, BIN_END,  AverageCov = meanCov)))


colnames(d5) <- c("chr", "start", "end", "variants", "Fst", "meanFst", "piHap", "piPle", "tajDHap", "tajDPle", "meanMQ", "meanCoverage")

write_delim(d5, "combined.fst.pi.tajD.MQ.Cov.forR.tsv", delim = "\t")

# readable file in excel
write.table(d5, "combined.fst.pi.tajD.MQ.Cov.forexcel.tsv", sep = "\t", dec = ",", col.names = T, row.names = F)


# after visualizing the distribution of these metrics in our dataset
# we decided to drop out windows with less than a SNP/kb, coverage > 20, MQ < 40.
data <- df[df$variants >=100 &
             df$variants < quantile(df$variants, 0.98) & #keep windows with the number of variants bellow the top 2%
             df$meanMQ >=40 &
             df$meanCoverage <= 20,]
```


### IV- Genome scan for selective sweeps

To detect signatures of selective sweeps, we  applied an _F_<sub>ST</sub> outlier approach (more in [Kelley et al. 2006](https://genome.cshlp.org/content/16/8/980) and [Akey et al. 2010](https://doi.org/10.1073/pnas.0909918107)) and used the cross-population extended haplotype homozygosity (_xpEHH_) and the site-specific extended haplotype homozygosity (_EHHS_) tests ([Sabeti et al. 2007](https://doi.org/10.1038/nature06250) and [Tang et al. 2007](
https://doi.org/10.1371/journal.pbio.0050171)).

We only used unrelated individuals from each population to perform the scans (see **2) Computing genetic differentiation (_F_<sub>ST</sub>) using VCFtools**)

##### 1) Getting _F_<sub>ST</sub> outlier windows

We calculated _F_<sub>ST</sub> above and have generated a single file giving not only _F_<sub>ST</sub> estimates, but also others (e.g. π, Tajima's D and mean MQ) in 100 kb non-overlapping windows across the genome. Now in R:

```R
# read data
df <-read_delim("combined.fst.pi.tajD.MQ.Cov.forR.tsv", delim = "\t", col_names = T)
df<- na.omit(df)

# apply filters based on variant number, MQ and coverage
data <- df[df$variants >=100 &
             df$variants < quantile(df$variants, 0.98) &
             df$meanMQ >=40 &
             df$meanCoverage <= 20,]

# get Fst estimates from the combined data file
fst <- data[,c(1:6)]
colnames(fst) <- c("CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST")
fst$Population<- c("Haplometrotic-Pleometrotic")

cutoff <-  quantile(abs(fst$WEIGHTED_FST), 0.95)# get cutoff value

# plot distribution
fhist <- ggplot(fst,aes(x=WEIGHTED_FST, fill=Population, color =Population,after_stat(count))) +
  theme_bw()+ ylab("Count")+
  geom_histogram(aes(fill=WEIGHTED_FST >cutoff, colour =WEIGHTED_FST >cutoff), binwidth = .01, alpha =.9, size=.1)+
  geom_density(aes(y=.01 * ..count..), alpha=0, color='black')+
  geom_vline(xintercept = cutoff , color = "blue", linetype = 2)+
  scale_colour_manual(name = '', values = c('grey','grey','red'))+
  scale_fill_manual(name = '', values = c('grey','grey','red'))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,200))+
  scale_x_continuous(name = expression(paste (italic("F")["st"])),breaks=c(0,.2,.4,.6,.8), limits = c(0,.8))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(size=10,face = "bold"),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.line.x = element_blank(),
        legend.position = "none")+
  ggtitle("")

fhist


# get candidate regions now: Fst estimates in the top 5%
outliersd <- data[data$Fst >quantile(abs(data$Fst), 0.95),]

# in which population there is selection? based on pi and Tajima's D estimates
# logic here is that selection will lead to reduced diversity (tajima's D) in the population where it's acting
outliersd <- mutate(outliersd, selection= ifelse(outliersd$tajDPle < outliersd$tajDHap & outliersd$piPle < outliersd$piHap, "pleometrotic", ifelse(outliersd$tajDPle > outliersd$tajDHap & outliersd$piPle > outliersd$piHap, "haplometrotic", "ambiguous")))

# store this file!!

write_delim(outliersd, "~/candidate.region.q95.forR.tsv", delim = "\t")

# ambiguous windows where it's hard to decide were droped out
```
Adjacent windows can then be merged into one region subject of selection using [BEDtools' merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html) function.

##### 2) Detecting selection by applying extended haplotype homozygosity statistics

The _xpEHH_ analyses were performed using the phased data from unrelated queens with the R package [rehh](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html). Using the phased data produced above using SHAPEIT, we will need to generate a vcf file for each scaffold and in each population independently.

```bash
# get phased data for unrelated queens in the pleometrotic population
vcftools --gzvcf Onlyphased.merged.vcf.gz --keep unrelated.pleo.txt --recode --out Onlyphased.merged.unrel.pleo

# get phased data for unrelated queens in the haplometrotic population
vcftools --gzvcf Onlyphased.merged.vcf.gz --keep unrelated.haplo.txt --recode --out Onlyphased.merged.unrel.haplo

# get data for each scaffold individually
for i in {1..27}; do bcftools view -r scaffold_$i Onlyphased.merged.unrel.pleo.vcf.gz -Ov -o ~/rehh.run/pleo.scaffold_$i.vcf; done
for i in {1..27}; do bcftools view -r scaffold_$i Onlyphased.merged.unrel.haplo.vcf.gz -Ov -o ~/rehh.run/haplo.scaffold_$i.vcf; done
```
* ###### Running the _xpEHH_ analysis using rehh in R

```R
# load libraries
library(tidyverse)
library(rehh)
library(dplyr)

## for the haplometrotic population using the largest 25 scaffolds (loop)

for(i in 1:25) {

  hap_file = paste("~/rehh.run/haplo.scaffold_", i, ".vcf", sep = "")

  haplo_hh <- data2haplohh(hap_file = hap_file,
                     polarize_vcf = FALSE)
  haplo_hhsub <- subset(haplo_hh, min_maf = 0.05)#only SNPs with maf of 0.05

  haplo_scan <- scan_hh(haplo_hhsub,
                  discard_integration_at_border = F,
                  polarized = FALSE, scalegap = 20000, maxgap = 200000) # to avoid false positives due to large gaps between markers

  if (i == 1) {
    haplo_wgscan <- haplo_scan
  } else {
    haplo_wgscan <- rbind(haplo_wgscan, haplo_scan)
  }
}



## for the haplometrotic population using the largest 25 scaffolds (loop)

for(i in 1:25) {

  hap_file = paste("~/rehh.run/pleo.scaffold_", i, ".vcf", sep = "")

  pleo_hh <- data2haplohh(hap_file = hap_file,
                     polarize_vcf = FALSE)
  pleo_hhsub <- subset(pleo_hh, min_maf = 0.05)

  pleo_scan <- scan_hh(pleo_hhsub,
                  discard_integration_at_border = F,
                  polarized = FALSE, scalegap = 20000, maxgap = 200000) # to avoid false positives due to large gaps between markers


  if (i == 1) {
    pleo_wgscan <- pleo_scan
  } else {
    pleo_wgscan <- rbind(pleo_wgscan, pleo_scan)
  }
}




# Now perform an xp-ehh analysis: The xpEHH compares the length of haplotypes between populations
# to detect selective sweeps in which the selected allele has approached or reached fixation in
# at least one population

haplo_pleo <- ies2xpehh(pleo_wgscan, haplo_wgscan,
                       popname1 = "pleo", popname2 = "haplo",
                       include_freq = T,p.adjust.method = "BH") # adjust for FDR

# highly positive values indicate selection in population 1 (pleo),
# whereas negative values suggest selection in population 2 (halpo in this case)


# Reformat data: good for plotting
wgscan.haplo.pleo.xpehh.qqman <- data.frame(
  CHR = as.integer(factor(haplo_pleo$CHR,
                          levels = unique(haplo_pleo$CHR))),

  BP = haplo_pleo$POSITION,         # base pairs
  P = 10**(-haplo_pleo$LOGPVALUE),  # transform back to p-values
  PP = haplo_pleo$XPEHH_pleo_haplo,
  SNP = row.names(haplo_pleo)       # SNP names
)

  # write out haplo pleo xpEHH
tib <- as_tibble(wgscan.haplo.pleo.xpehh.qqman)

# save the final results for further processing
write_tsv(tib, "~/rehh.main/haplo.pleo.25scf.xpEHH.tsv")


# now extract outlier SNPs (p < 0.05) and then get candidate regions

# identify and cluster outliers:
myData <- tib
myData <- filter(myData, P < 0.05)
myData <- arrange(myData, CHR, BP)

# cluster outlier SNPs into regions (if they are less than 100 kb away)
threshold <- 100000
clust <- unlist(sapply(unique(myData$CHR), function(w){
y <- filter(myData, CHR == w) %>% .$BP
y <- cumsum(c(1, diff(y) > threshold))
paste0(w, "_", y)
}))
names(clust) <- NULL
myData$cluster <- factor(clust)

# get clusters and their span
myCluster <- myData %>% dplyr::group_by(cluster) %>% dplyr::summarise(start = min(BP), stop = max(BP), size = stop - start)


# identify highest peak in each cluster
myData2 <- myData %>% group_by(cluster) %>% filter(P == min(P)) # get most significant SNP(s) in each cluster
myData2 <- dplyr::mutate(as_tibble(myData2), selection= ifelse(myData2$PP >0, "pleometrotic","haplometrotic")) # positive XPEHH scores indicate selection in the pleometrotic population

# the myData2 contains the position and associated stats for each top significant SNP in each cluster.

# save it!!
write_delim(myData2, "top.sig.snp.in.each.cluster.25scf.forR.tsv", delim = "\t")

```

The resulting file can then be used to extract candidates for positive selective sweeps (e.g. regions that are within -/+ 100 kb around each significant SNP of each cluster).

* ###### Computing _EHHS_ for two regions one on scaffold_15 and the other on scaffold_16
**Scaffold_15**

```R
library(rehh)

#scaffold 15: id174340 (most sgnificant SNP)

# read in data for each population
# haplo
haplo_hh <- data2haplohh(hap_file = "~/rehh.main/haplo.scaffold_15.vcf",
                       polarize_vcf = FALSE)
# pleo
pleo_hh <- data2haplohh(hap_file = "~/rehh.main/pleo.scaffold_15.vcf",
                      polarize_vcf = FALSE)


# filter on MAF - here 0.05
haplo_hh_f <- subset(haplo_hh, min_maf = 0.05)
pleo_hh_f <- subset(pleo_hh, min_maf = 0.05)


# haplo
res <- calc_ehhs(haplo_hh,
               mrk = "id174340",
               include_nhaplo = TRUE,
               discard_integration_at_border = FALSE)

# reformat data
nehhs_haplo <- as_tibble(res$ehhs)


# pleo
res1 <- calc_ehhs(pleo_hh,
                mrk = "id174340",
                include_nhaplo = TRUE,
                discard_integration_at_border = FALSE)

# reformat data
nehhs_pleo <- as_tibble(res1$ehhs)

colnames(nehhs_haplo) <- c("POSITION",  "hEHHS",  "hNEHHS", "hNHAPLO")
colnames(nehhs_pleo) <- c("POSITION",  "pEHHS",  "pNEHHS", "pNHAPLO")



# make plot
scf15 <- ggplot()+
geom_line(data = nehhs_pleo, aes(x=POSITION/1e6, y=pNEHHS), colour="#2BCE48")+
geom_line(data = nehhs_haplo, aes(x=POSITION/1e6, y=hNEHHS), colour="#0075DC")+
labs(y = expression(italic("EHHS"))) + xlab(NULL)+
theme_bw() +
theme(legend.position = "none",
      panel.border = element_blank(),
      plot.title = element_text(size=10,face = "bold"),
      plot.subtitle = element_text(size = 8, hjust = .55),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"))+
xlim(c(.1,.3))+ ggtitle(label = "C", subtitle = expression(italic("scaffold 15:211998")))

scf15

```

**Scaffold_16**
```R
# Scaffold 16 : id184640

# read in data for each population
# haplo
haplo_hh <- data2haplohh(hap_file = "~/rehh.main/haplo.scaffold_16.vcf",
                       polarize_vcf = FALSE)
# pleo
pleo_hh <- data2haplohh(hap_file = "~/rehh.main/pleo.scaffold_16.vcf",
                      polarize_vcf = FALSE)


# filter on MAF - here 0.05
haplo_hh_f <- subset(haplo_hh, min_maf = 0.05)
pleo_hh_f <- subset(pleo_hh, min_maf = 0.05)


res <- calc_ehhs(haplo_hh,
               mrk = "id184640",
               include_nhaplo = TRUE,
               discard_integration_at_border = FALSE)
plot(res, nehhs = TRUE)



nehhs_haplo <- as_tibble(res$ehhs)


res1 <- calc_ehhs(pleo_hh,
                mrk = "id184640",
                include_nhaplo = TRUE,
                discard_integration_at_border = FALSE)


nehhs_pleo <- as_tibble(res1$ehhs)

colnames(nehhs_haplo) <- c("POSITION",  "hEHHS",  "hNEHHS", "hNHAPLO")
colnames(nehhs_pleo) <- c("POSITION",  "pEHHS",  "pNEHHS", "pNHAPLO")


# make plot
scf16 <- ggplot()+
geom_line(data = nehhs_pleo, aes(x=POSITION/1e6, y=pNEHHS), colour="#2BCE48")+
geom_line(data = nehhs_haplo, aes(x=POSITION/1e6, y=hNEHHS), colour="#0075DC")+
labs(x = "Position (Mb)",
     y = expression(italic("EHHS"))) +
theme_bw() +
theme(legend.position = "none",
      plot.title = element_text(size = 8, hjust = .55),
      panel.border = element_blank(),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"))+ #xlim(c(2,2.3))+
ggtitle(label = expression(italic("scaffold 16:2154248")))

scf16
```


##### 3) Characterization of candidate genes and SNPs
* **_Get candidate regions by overlapping results from both selection analyses_**


```bash
# get candidate regions from the xpEHH analysis
# (e.g. regions that are within -/+ 100 kb around each significant SNP of each cluster)
##pleo
cat top.sig.snp.in.each.cluster.25scf.forR.tsv |grep -v "^C"| grep "pleometrotic"|awk 'OFS=OF="\t" {print "scaffold_"$1, $2-100000, $2+100000, $2}' | bedtools merge -i - > candidate.regions.within.100kb.25scf.pleo.merged.bed
##haplo
cat top.sig.snp.in.each.cluster.25scf.forR.tsv |grep -v "^C"| grep "haplometrotic"|awk 'OFS=OF="\t" {print "scaffold_"$1, $2-100000, $2+100000, $2}' | bedtools merge -i - > candidate.regions.within.100kb.25scf.pleo.merged.bed


# get candidate regions from Fst analysis:
## Pleo
cat candidate.region.q95.forR.tsv |grep -v "^chr"| grep "pleometrotic"|awk 'OFS=OF="\t" {print "scaffold_"$1, $2, $3}'| bedtools merge -i - > candidate.regions.in.pleo.merged.bed

##Haplo
cat candidate.region.q95.forR.tsv |grep -v "^chr"| grep "haplometrotic"|awk 'OFS=OF="\t" {print "scaffold_"$1, $2, $3}'| bedtools merge -i - > candidate.regions.in.haplo.merged.bed


# overlap candidates regions from both tests
bedtools intersect -a candidate.regions.in.pleo.merged.NEW.bed -b ../../xpehh/regions/candidate.regions.within.100kb.25scf.pleo.merged.bed  > 17.candidates.regions.from.both.pleo.analyses.bed

# Only regions identified as candidates for selective sweeps by the two scan approaches were kept.
# This includes 17 regions in the pleometrotic population only.

```
* **_SNP annotation using [SnpEff](http://pcingola.github.io/SnpEff/)_**
```bash
# exctract SNPs found in candidates regions
vcftools --gzvcf all_snpsPASS_filtered.vcf.gz --bed 17.candidates.regions.from.both.pleo.analyses.bed --recode --recode-INFO-all --out candidate.snps.merged

# Then run the annotation
java -jar ~/softwares/snpEff/snpEff.jar -v -stats candidate.snps.merged.html Pcal.v2 candidate.snps.merged.vcf.gz >snpeffout/candidate.snps.merged.ann.vcf
```
* **_Putative function of candidate genes_**
```bash
# extract candidate genes
bedtools  getfasta -fi Pcal.genome.clean.sort.masked.fa -bed 17.candidates.regions.from.both.pleo.analyses.bed >  candidate.genes.flagged.fa -name

# then run blast
blastx -query candidate.genes.flagged.fa -db ncbi-nr/2020_01_08/nr -max_target_seqs 1 -max_hsps 2 -out candidate.genes.flagged.blastx.results -num_threads 25 -outfmt '6 qseqid sseqid evalue bitscore sgi sacc qaccver saccver staxids scomnames stitle sskingdom'
```

* **_Drosophila melanogaster orthologs to the candidate genes using [OrthoFinder](https://github.com/davidemms/OrthoFinder)_**

```bash
# Dowload Dmel transcripts from flybase (r6.37): http://ftp.flybase.net/releases/FB2020_06/dmel_r6.37/fasta/

# extract longest isoform Dmel:
seqkit fx2tab -l dmel-all-translation-r6.37.header.edited.fasta |sort -k2,2 -k4,4nr |sort -k2,2 -u -s |awk '{print ">"$1";length="$3"\n"$2}'> Dmel.longestIsoform.fa
# extract longest isoform Pcal:
seqkit fx2tab -l geneAnnotation/Pogonomyrmex_californicus.proteins.fa |sort -k2,2 -k4,4nr |sort -k2,2 -u -s |awk '{print ">"$1";symbol="$2";length="$4"\n"$3}'> Pcal.longestIsoform.fa

# put the two newly produced fa files in a new folder: OrthoInputs and then
mkdir OrthoInputs
mv *longestIsoform.fa OrthoInputs

# run orthofinder
orthofinder -f OrthoInputs/ -t 20
```


* **_Search for candidate genes differentially expressed in [Helmkampf et al. (2016)](https://doi.org/10.1111/mec.13700)_**

```bash
# Download the 7890 assembled and filtered transcripts (Appendix S1) from https://onlinelibrary.wiley.com/doi/10.1111/mec.13700
# make blast db (edited header sed -e 's/\ .*//' mec13700-sup-0002-appendixs1.fas)

makeblastdb -in mec13700-sup-0002-appendixs1.fas -dbtype nucl -title "db" -out "transcripts"

# extract exons from candidate regions
bedtools intersect -wa -a Pogonomyrmex_californicus.gff3 -b 17.candidates.regions.from.both.pleo.analyses.bed|awk '{if ($3=="exon") print $0}'|gff2bed|sort|uniq|sort -k1,1V -k2,2n > candidate.exons.pleo.bed

# and Then
bedtools  getfasta -fi Pcal.genome.clean.sort.masked.fa -bed candidate.exons.pleo.bed -fo candidate.exons.pleo.fa -name

# blast
blastn -query candidate.exons.pleo.fa -db db/transcripts.fa -max_target_seqs 1 -max_hsps 5 -out candidate.exons.pleo.tsv -outfmt 6
```

* **_Genotype phenotype association_**
Prepare input files for R.

```bash
# Use SnpEff's output and the list of gene (145 in total) in candidate regions
# Extract SNP annotations of these genes using bcftools
bcftools view -R 145genes.bed -Ov -o shared.snps.ann.in.145genes.vcf snpeffout/candidate.snps.merged.ann.vcf

# reformat annoations
bcftools query -f '%ANN\n' shared.snps.ann.in.145genes.haplo.vcf |tr -s "|"|tr "," "\n"|tr "|" "\t" > annotationsForR.tsv


# get genotypes for haplo queens
view -S haplo.txt -Ov -o shared.snps.ann.in.145genes.haplo.vcf shared.snps.ann.in.145genes.vcf # haplo.txt lists queen names
# next
bcftools query -f '%CHROM\t%POS\t%ALT[\t%GT]\n' shared.snps.ann.in.145genes.haplo.vcf > genotypes.haplo.forR.tsv

# get genotypes for pleo queens
view -S pleo.txt -Ov -o shared.snps.ann.in.145genes.pleo.vcf shared.snps.ann.in.145genes.vcf # pleo.txt lists queen names
# next
bcftools query -f '%CHROM\t%POS\t%ALT[\t%GT]\n' shared.snps.ann.in.145genes.pleo.vcf > genotypes.pleo.forR.tsv

# get list of genes with effects
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\n' shared.snps.ann.in.145genes.haplo.vcf| tr -s "|"|tr "," "\n" > listOfgenesWithEffect.tmp.tsv
# had to manually edit the listOfgenesWithEffect.tmp.tsv file
# some entries will be missing the first four columns ()<scaffold><position><REF><ALT>)
# this is because some SNPs affect more than one gene.
# after doing so
cat listOfgenesWithEffect.tmp.tsv |tr "|" "\t"|cut -f1-8 > listOfgenesWithEffect.tmp.tsv
```
Now in R combine all relevant information into one file and get genotype count.
```r
# get for each gene the number of effects (e.g. nonsense, synonymous, intron variants...)
## load data
ann <- read.table("annotationsForR.tsv", header = F, sep = "\t", fill = T)

colnames(ann) <- c("alt", "ann", "impact", "genename", "geneid", "feature",
                   "featureid","prot", "score", "direction", "AAchange",
                   "score2", "score3", "score4", "warnings")
head (ann)

list <- ann %>%
  dplyr::count(geneid, ann,sort = TRUE)

mdata <- tidyr::pivot_wider(list, names_from = ann, values_from = n, values_fill=0)

# save it!
write.table(mdata, "geneListWithEffectCounts.pivoted.tsv", quote = F, col.names = T, row.names = F, sep = "\t")


# get genotypes/allele counts (used to get allele frequencies at variants affecting candidate genes)
## load data
h <- read.table("genotypes.haplo.forR.tsv", header = F, sep = "\t", fill = T)

h$Href <- apply(h, 1, function(x) length(which(x=="0/0"|x=="0/1")))
h$Halt <- apply(h, 1, function(x) length(which(x=="1/1"|x=="0/1")))
h$Hmissing <- apply(h, 1, function(x) length(which(x=="./.")))


p <- read.table("genotypes.pleo.forR.tsv", header = F, sep = "\t", fill = T)

p$Pref <- apply(p, 1, function(x) length(which(x=="0/0"|x=="0/1")))
p$Palt <- apply(p, 1, function(x) length(which(x=="1/1"|x=="0/1")))
p$Pmissing <- apply(p, 1, function(x) length(which(x=="./.")))


eff <- read.table("listOfgenesWithEffect.tsv", header = F, sep = "\t", fill = T)

com <- full_join(dplyr::select(eff,c(V1:V8)),
                 dplyr::select(h,c(V1:V3,Href,Halt,Hmissing)))

com2 <- full_join(dplyr::select(com,c(V1:V8, Href,Halt,Hmissing)),
                  dplyr::select(p,c(V1:V3,Pref,Palt,Pmissing)))


# save file containing for each candidate gene: the position, ref and alt alleles, type of effect, gene ID and the number of queens carrying the ref/alt allele
write.table(com2, "genotypephenotype.withmissingCount.tsv", quote = F,col.names = T, row.names = F, sep = "\t")

# This files was used to filter the candidate genes and to narrow down the list to only 38 genes (with six very interesing ones)
```

### V- RAD sequencing and analysis

Raw RADseq reads were first filtered from adapter sequences using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) `ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:100` and then demultiplexed using [process_radtags](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php). Mapping was produced using BWA-MEM with default parameters (see above). After generating and sorting the BAM files, use  [STACKS](https://catchenlab.life.illinois.edu/stacks/manual/) to call variants:

```bash
# after generating bamfiles and sorting them:
# inside the folder with all alignment files:
perl ~/programs/STACKS/bin/ref_map.pl -T 30 --popmap ./popmap2.tsv -o ./stacks.2/ --samples ./
#popmap2.tsv is a two-columns files column 1 is the bam alignment file name for each sample
#the second column is population (for DNA_1.bam it will be <DNA_1><pop1>)

# Once the ref_map.pl is done
mkdir vcf
# then 	
~/programs/STACKS/bin/populations -P ./stacks.2/ -O ./vcf/ vcf
```
Next, to obtain a set of reliable markers from the raw calls, use VCFtools (to keep biallelic markers without missing data and a MAF ≥ 0.45)

```bash
vcftools --vcf populations.snps.1pop.vcf --recode-INFO-all --maf 0.45 --max-alleles 2 --min-alleles 2 --max-missing 1 --out populations.snps.1pop.biall.no-miss.maf.45 --recode
```
Finally, we excluded heterozygous calls, likely to be sequencing errors (males are haploid and therefore should not present heterozygous calls). To properly handle the vcf file using bcftools, first add missing contig information to the header (there are probably more efficient ways to do it).
```bash
# get scaffold information
awk '{printf("##contig=<ID=%s,length=%d>\\n",$1,$2);}' Pcal2.reference.fa.fai # can be obtained by running samtools faidx Pcal2.reference.fa

awk '/^#CHROM/ { printf("REPLACE WITH OUTPUT FROM PREVIOUS COMMAND");} {print;}' populations.snps.1pop.biall.no-miss.maf.45.vcf > populations.snps.1pop.biall.no-miss.maf.45.header.vcf

# get only homozygous calls
grep -v -w "0/1" populations.snps.1pop.biall.no-miss.maf.45.header.vcf > populations.snps.1pop.biall.no-miss.maf.45.header.Hom.vcf

# prepare marker for Multipoint
bcftools query -f '%CHROM\t%POS\t[\t%GT]\n' populations.snps.1pop.biall.no-miss.maf.45.header.Hom.vcf > populations.snps.1pop.biall.no-miss.maf.45.header.Hom.tsv

# diploid to haploid markers
sed -i '' "s@1/1@1@g" populations.snps.1pop.biall.no-miss.maf.45.header.Hom.tsv
sed -i '' "s@0/0@0@g" populations.snps.1pop.biall.no-miss.maf.45.header.Hom.tsv
```
### VI- Characterization of the supergene

**_Visualizing genotypes of males and queens at LG14 using [VIVA](https://compbiocore.github.io/VariantVisualization.jl/latest/)_**

```bash
# extract SNPs shared between males (linkage map) and queens (pop gen)
bcftools isec -p shared.snps haplo.linkageMap/populations.snps.1pop.biall.no-miss.maf.45.header.Hom.vcf.gz ../vcfs/Pcal.BQSR.snpsPass.biall.mac3.miss.85.3DP24.vcf.gz # requires vcf files to be bgzip/tabix
# the files we need are the two vcfs that hove records from both males and queens vcfs. See the README.txt inside shared.snps

# merge SNPs from males (also found in queens) and queens (also found in males)
bcftools merge -Oz -o merged.1064.snps.vcf.gz 0002.vcf.gz 0003.vcf.gz # 0002.vcf.gz and 0003.vcf.gz can be found in shared.snps

# extract markers on LG14
bcftools view -r scaffold_14,scaffold_22,scaffold_23,scaffold_38,scaffold_42,scaffold_48,scaffold_6 merged.1064.snps.vcf.gz -Ov -o Supergene.LM.PopGen.snps.vcf

```
Viva requires some pre formatting of the vcf file e.g. doesn't accept characters in contig names. Use something like `sed -i '' 's/scaffold_//g' Supergene.LM.PopGen.snps.vcf` to remove characters and keep only contig numbers and then re order if needed. You also need a popmap file (pop.meta.csv) that lists all indidvidual names and the population they belong to. More on file format (and more) can be found at [VIVA's GitHub](https://github.com/compbiocore/VariantVisualization.jl) repository. The tool can be ran as follow:
```bash
viva -f Supergene.LM.PopGen.snps.vcf --heatmap genotype -g pop.meta.csv LM,PopGen -x --save_remotely -o genotype_supergene.LM.PopGen.snps
```
To visually inspect the genotypes of queens at the supergene using VIVA:

```bash
# extract and thin SNPs on LG14.Thinning is needed to reduce the number of SNPs for easy paining using viva. one can also exclude rare alleles i.e. --maf 0.2
vcftools --gzvcf Pcal.BQSR.snpsPass.biall.mac3.miss.85.3DP24.vcf.gz --chr scaffold_14 --chr scaffold_22 --chr scaffold_23 --chr scaffold_38 --chr scaffold_42 --chr scaffold_48 --chr scaffold_6 --thin 1000 --recode Pcal.BQSR.snpsPass.biall.mac3.miss.85.3DP24.LGsupergene.thinned

# Now run Viva after formatting the vcf file i.e. remove characters and keep only contig numbers
viva -f Pcal.BQSR.snpsPass.biall.mac3.miss.85.3DP24.LGsupergene.thinned.vcf --heatmap genotype -g pop.meta.csv haplometrotic,pleometrotic -x --save_remotely -o genotype_supergene.queens.snps

```
**_PCA on SNPs found in the non-recombining supergene_**

```bash
# convert vcf to plink format
vcftools --gzvcf merged.1064.snps.vcf.gz --chrom-map chrom-map.txt --out merged.snps --plink # chrom-map.txt lists contigs see above how to get such file
# convert to bed format
plink --file merged.snps --aec --make-bed --out merged.snps.b

# run PCA on SNPs found on the supergene area
plink --bfile ../merged.1064.snps.b --extract range range.txt --pca 4 --aec --out merged.1064.snps.supergene.pca
# ranges.txt is tab delimited files with four columns giving the supergene coordinates for each scaffold (rows)
# <Chrom><start><end><id(random)>
```
**_Pairwise LD at LG14 using [LDheatmap](https://sfustatgen.github.io/LDheatmap/)_**
Use [flo]() to generate a chain file. initially flo is designed to lift over genome annotations, but if no genome annotation is provided, the program will stop after generating a chain file. Check the [flo]() repository for a detailed description of how to proceed. Next use [Picard's LiftoverVcf](https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard-) to lift over SNPs from scaffolds 22, 23, 38, 42, 48, and parts of 6 and 14 to one linkage group (LG14).
```bash
java -jar picard.jar LiftoverVcf I=biall.mac3.miss.85.3DP24.recode.vcf O=LG14.biall.mac3.miss.85.3DP24.liftedover.vcf CHAIN=liftover.chn REJECT=rejected_variants.vcf R=LG14.fa

# biall.mac3.miss.85.3DP24.recode.vcf: SNP file with markers on scaffolds 22, 23, 38, 42, 48, and parts of 6 and 14
# liftover.chn: chain file produced using flo
# rejected_variants.vcf: conatains SNPs that cannot be lifted over
# LG14.fa: genome file as one sequence (in this order scaffolds 14, 22, 23, 38, 42, 48 and 6) with one header line >LG14. i.e. this is the sequence of LG14
# LG14.biall.mac3.miss.85.3DP24.liftedover.vcf: markers with now coordinates on LG14
```
Now in R:

```R
# load libraries
library(vcfR)
library(LDheatmap)

# choose colors:
mycols <- c("sienna3", "sienna2", "bisque1", "gray85", "lightsteelblue", "skyblue1")
mypalette <- colorRampPalette(mycols)(18)


#read the VCF file (with Sh queens and a maf of 0.2)

# the vcfs are also thinned to make it possible to compute pairwise r2
vcf <- read.vcfR("lifted/LG14.biall.mac3.miss.85.3DP24.liftedover.MAF0.2.Sh.thinned.recode.vcf")

#extract genetic distances, subject IDs and genotypes for Sh queens from the vcfR object
sh <- vcfR2SnpMatrix(vcf)

#draw the heatmap
LDheatmap(sh$data, sh$genetic.distance, title='Sh queens', add.map=T, flip = T, color = mypalette)


#Same with Sp queens
vcf <- read.vcfR("lifted/LG14.biall.mac3.miss.85.3DP24.liftedover.MAF0.2.Sp.thinned.recode.vcf")
sp <- vcfR2SnpMatrix(vcf)
LDheatmap(sp$data, sp$genetic.distance, title='Sp queens', add.map=T, flip = T, color = mypalette)


#Same with pooled Sh and Sp queens
vcf <- read.vcfR("lifted/LG14.biall.mac3.miss.85.3DP24.liftedover.MAF0.2.ShSp.thinned.recode.vcf")
shsp <- vcfR2SnpMatrix(vcf)

LDheatmap(shsp$data, shsp$genetic.distance, title='pooled Sh & Sp queens', add.map=T, flip = T, color = mypalette)
```
**_TE and gene content_**

```bash

# get scaffold sizes
cat Pcal2.reference.fa.fai |awk '{print $1"\t"$2}'|sort -V -k1,1 > chrom-size.txt
# make bed window file (100 kb sliding windows)
bedtools makewindows -g chrom-size.txt -w 100000 |bedtools sort > 100kbwindows.bed
```
###### For transposable elements (TEs)
```bash
# sort TE annotation file
cat PcalTE.gff3|bedtools sort -i - > PcalTE.sorted.gff3

# calculate content in 100 kb windows
sortBed -i PcalTE.sorted.gff3 | gff2bed  |bedmap --echo --bases-uniq-f 100kbwindows.bed - |tr "|" "\t" > 100kbwindows.TEcontent.bed
```

###### For genes

```bash
# generate a gff file with genes
cat Pogonomyrmex_californicus.gff3 | awk '{if ($3=="gene") print $0}'|bedtools sort -i - > GENE.sorted.gff

#then calculate coverage
$ sortBed -i GENE.sorted.gff | gff2bed  |bedmap --echo --bases-uniq-f 100kbwindows.bed - |tr "|" "\t" > 100kbwindows.GENEcontent.bed
```
### VII- Dating the supergene
###### computing the rCCR bteween Sh and Sp
MSMC2 was used to compute rCCR between Sh and Sp using four Sh queens and four Sp queens. To do so we extracted SNPs found on scaffolds 22, 23 and parts of 6 (parts of the supergene that are > 1 Mb).

```bash
# Sh queens (H*) and Sp queens (P*)
samples=(H104W H38W H73W H98W P114W P50R P84P P97R)

# Get SNPs on each scaffold and for each queen
for i in 6 22 23; do for j in ${samples[*]}; do  bcftools view -s $j scaffold_$i.Onlyphased.header.vcf -Oz -o supergene.Ne/$j.scaffold_$i.phased.vcf.gz; done; done
# scaffold_$i.Onlyphased.header.vcf: SHAPEIT output with headers

# get only SNPs in the supergene area on scaffold 6
for j in ${samples[*]}; do bcftools view -r scaffold_6:0-2020227 -Oz -o supergene.Ne/$j.scaffold_6.phased.superegene.vcf.gz supergene.Ne/$j.scaffold_6.phased.vcf.gz; done
```
Then run MSMC2 ([Schiffels & Wang, 2020](https://link.springer.com/protocol/10.1007/978-1-0716-0199-0_7)) to calculate the relative cross-coalescence rate (rCCR) as indicated by the developers.

###### phylogenetic analysis of Sh and Sp

First, call SNPs using ``g.vcf.gz`` file generated for Pcal (see variant calling above) and the ``g.vcf.gz`` file of _Pogonomyrmex subnitidus_ generated as described above for _P. californicus_.

##### 4) Variant calling using [GATK](https://gatk.broadinstitute.org/hc/en-us)

```bash
# combine gVCFs
gatk CombineGVCFs \
-R Pcal2.reference.fa \
-V input.list \ # this file contains two lines: 1) all.g.vcf.gz and 2) Psub.g.vcf.gz
-O Pcal.Psub.g.vcf.gz

# genotype the combined gVCF
gatk GenotypeGVCFs \
-R Pcal2.reference.fa \
-V Pcal.Psub.g.vcf.gz \
-O Pcal.Psub.vcf.gz


# Apply hard filters using
# get SNPs only
gatk SelectVariants -V Pcal.Psub.vcf.gz --select-type-to-include SNP -O Pcal.Psub.snps.vcf.gz

# c- apply hard filters:
gatk VariantFiltration \
-V Pcal.Psub.snps.vcf.gz  \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -5.0" --filter-name "MQRankSum-5" \
-filter "ReadPosRankSum < -5.0" --filter-name "ReadPosRankSum-5" \
-O Pcal.Psub.snps.filtered.vcf.gz

# extract variants that have passed the applied filters
gatk SelectVariants -V Pcal.Psub.snps.filtered.vcf.gz -select 'vc.isNotFiltered()' -O Pcal.Psub.snpsPASS.vcf.gz

# get SNPs on the supergene
bcftools view -R supergene.bed -Oz -o Pcal.Psub.BQSR.snps.supergeneRegion.vcf.gz Pcal.Psub.snpsPASS.vcf.gz
# supergene.bed: coordinates of the superegene on each scaffold (all of scaffolds 22, 23, 38, 42, 48, and parts of 6 and 14)

# get bi-allelic an polymorphic SNPs with a missingness < 0.2 for six Sh and six Sp, and Psup (outgroup)
bcftools view -s H104W,H33G,H34G,H38W,H73W,H98W,P22P,P102R,P114W,P109R,P34R,P84P,Psub -e 'AC==0 // AC==AN // F_MISSING > 0.2' -m2 -M2 -Ov -o filtered.dataset.vcf Pcal.Psub.BQSR.snps.supergeneRegion.vcf.gz
```
Convert vcf to nex format using [convert_vcf_to_nexus.rb](https://github.com/ForBioPhylogenomics/tutorials/blob/main/week2_src/convert_vcf_to_nexus.rb) by [Michael Matschiner](https://github.com/mmatschiner) and then manually reformat to match the [phy](https://cme.h-its.org/exelixis/resource/download/hands-on/dna.phy) file format. Then run [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/hands_on.html) with default parameters.
```bash
# first
raxmlHPC -m GTRGAMMA -p 12345 -# 20 -s Pcal.Psub.BQSR.snps.supergeneRegion.edited.fasta -n tree1

# then
raxmlHPC -m GTRGAMMA -p 12345 -b 12345 -# 100 -s Pcal.Psub.BQSR.snps.supergeneRegion.edited.fasta -n tree2

# finally
raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.tree1 -z RAxML_bootstrap.tree2 -n tree3
```
The tree is then visualized using [iTOL](https://itol.embl.de/)
