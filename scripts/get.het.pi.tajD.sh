#!/bin/bash
# This script will compute heterozygosity, Nucl. diversity, Tajima's D and Fst using VCFtools.
# Errbii M 2021

#get all variables you need..

echo "Type the name of the vcf file you are working on now, followed by [ENTER]:"

read vcf

echo "Type the window size here, followed by [ENTER]:"
read windsize

echo "Type the window step here, followed by [ENTER]:"
read windstep

echo "Type the name of the file containing indv. of population 1, followed by [ENTER]:"
read pop1
pop1base=$(basename $pop1 .txt)
echo "$pop1base"

echo "Type the name of the file containing indv. of population 2, followed by [ENTER]:"
read pop2
pop2base=$(basename $pop2 .txt)
echo "$pop2base"

base=$(basename $vcf .vcf.gz)

echo "$base"

mkdir $(echo "$base"."$windsize"bp"$windstep")

out=$(echo "$base"."$windsize"bp"$windstep")

#no_snps=$(bcftools view $vcf|grep -v "#"|wc -l)

#echo "This file contains $no_snps"

echo "calculating heterozygosity and IF..."

#het
$(vcftools --gzvcf $vcf --het --out ./$out/$base.het)



echo "calculating nucleotide diversity in pop1..."
#PI
$(vcftools --gzvcf $vcf --window-pi $windsize --window-pi-step $windstep --keep $pop1 --out ./$out/$base.$pop1base.pi)

echo "calculating nucleotide diversity in pop2..."
$(vcftools --gzvcf $vcf --window-pi $windsize --window-pi-step $windstep --keep $pop2 --out ./$out/$base.$pop2base.pi)




echo "calculating Tajima's D in pop1..."
#Tajima's D
$(vcftools --gzvcf $vcf --TajimaD $windsize --keep $pop1 --out ./$out/$base.$pop1base.tajD)

echo "calculating Tajima's D in pop2..."
$(vcftools --gzvcf $vcf --TajimaD $windsize --keep $pop2 --out ./$out/$base.$pop2base.tajD)



echo "calculating genetic differentiation.... COMING SOON"

echo "done with calculating 5 metrics! Now got to R and visualize the data"
