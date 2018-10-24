# LEAP_PARIS_IMPUTATION
Please note, this is work in progress

By and large, this follows the same format as here: https://github.com/autism-research-centre/SSC_liftover_imputation

But, I am writing the scripts down for ease of use.

## LEAP

QC

```bash
./plink --bfile LEAP_forVarun_20181005  --geno 0.1 --mind 0.05 --hwe 0.000001 --me 0.05 0.1  --make-bed --out QC1output
./plink --bfile QC1output --het --out QC1checkhet
./plink --bfile QC1output --check-sex --out QC1checksex
```

```R
sex = read.delim("QC1checksex.sexcheck", sep = "")
head(sex)
sexproblem = subset(sex, STATUS == "PROBLEM")


het = read.delim("QC1checkhet.het", sep = "")
het$HET = (het$N.NM. - het$O.HOM.)/het$N.NM. #create heterozygosity stats

mean = mean(het$HET)
sd = sd(het$HET)
het$Z = (het$HET - mean)/sd #create Z scores of heterozygosity

hetoutlier = subset(het, abs(Z) > 3)

het2 = hetoutlier[,c(1:2)]
sex2 = sexproblem[,c(1:2)]
failedsample = rbind(het2, sex2)

write.table(failedsample, file = "failedsample.txt", row.names = F, col.names = T, quote = F)
```

```bash
./plink --bfile QC1output --remove failedsample.txt --make-bed --out QC2output
```

Unlike the SSC file, we impute all the samples and then later, post-imputation, remove ancestry outliers

```bash
./plink --bfile QC2output --freq --out QC2output

perl ~/SFARI/liftOverPlink/HRC-1000G-check-bim.pl -b ~/LEAP_PARIS/QC2output.bim -f ~/LEAP_PARIS/QC2output.frq -r ~/SFARI/liftOverPlink/1000GP_Phase3_combined.legend -g -p EUR
```

Create files for imputation

```bash
./plink --bfile QC2output --exclude Exclude-QC2output-1000G.txt --make-bed --out TEMP1
./plink --bfile TEMP1 --update-map Chromosome-QC2output-1000G.txt --update-chr --make-bed --out TEMP2
./plink --bfile TEMP2 --update-map Position-QC2output-1000G.txt --make-bed --out TEMP3
./plink --bfile TEMP3 --flip Strand-Flip-QC2output-1000G.txt --make-bed --out TEMP4
./plink --bfile TEMP4 --reference-allele Force-Allele1-QC2output-1000G.txt --make-bed --out QC2output-updated
rm TEMP*
  
for i in {1..22}; do ./plink --bfile QC2output-updated --reference-allele Force-Allele1-QC2output-1000G.txt --chr ${i} --recode-vcf --out LEAP_file_chr${i}; done
for i in {1..22}; do vcf-sort LEAP_file_chr${i}.vcf | bgzip -c > LEAP_file_chr${i}.vcf.gz; done
for i in {1..22}; do rm LEAP_file_chr${i}.vcf | rm LEAP_file_chr${i}.log; done
```

Before continuing, delete files as necessary, and move needed files to seperate folder (LEAP_files)


## PARIS OmniExpress


QC

```bash
./plink --bfile Omniexpres_559ind_ForVarun  --geno 0.1 --mind 0.05 --hwe 0.000001 --me 0.05 0.1  --make-bed --out QC1output
./plink --bfile QC1output --het --out QC1checkhet
./plink --bfile QC1output --check-sex --out QC1checksex
```

```R
sex = read.delim("QC1checksex.sexcheck", sep = "")
head(sex)
sexproblem = subset(sex, STATUS == "PROBLEM")


het = read.delim("QC1checkhet.het", sep = "")
het$HET = (het$N.NM. - het$O.HOM.)/het$N.NM. #create heterozygosity stats

mean = mean(het$HET)
sd = sd(het$HET)
het$Z = (het$HET - mean)/sd #create Z scores of heterozygosity

hetoutlier = subset(het, abs(Z) > 3)

het2 = hetoutlier[,c(1:2)]
sex2 = sexproblem[,c(1:2)]
failedsample = rbind(het2, sex2)

write.table(failedsample, file = "failedsample.txt", row.names = F, col.names = T, quote = F)
```

```bash
./plink --bfile QC1output --remove failedsample.txt --make-bed --out QC2output
```

Unlike the SSC file, we impute all the samples and then later, post-imputation, remove ancestry outliers

```bash
./plink --bfile QC2output --freq --out QC2output

perl ~/SFARI/liftOverPlink/HRC-1000G-check-bim.pl -b ~/LEAP_PARIS/QC2output.bim -f ~/LEAP_PARIS/QC2output.frq -r ~/SFARI/liftOverPlink/1000GP_Phase3_combined.legend -g -p EUR
```

Create files for imputation

```bash
./plink --bfile QC2output --exclude Exclude-QC2output-1000G.txt --make-bed --out TEMP1
./plink --bfile TEMP1 --update-map Chromosome-QC2output-1000G.txt --update-chr --make-bed --out TEMP2
./plink --bfile TEMP2 --update-map Position-QC2output-1000G.txt --make-bed --out TEMP3
./plink --bfile TEMP3 --flip Strand-Flip-QC2output-1000G.txt --make-bed --out TEMP4
./plink --bfile TEMP4 --reference-allele Force-Allele1-QC2output-1000G.txt --make-bed --out QC2output-updated
rm TEMP*
  
for i in {1..22}; do ./plink --bfile QC2output-updated --reference-allele Force-Allele1-QC2output-1000G.txt --chr ${i} --recode-vcf --out LEAP_file_chr${i}; done
for i in {1..22}; do vcf-sort LEAP_file_chr${i}.vcf | bgzip -c > LEAP_file_chr${i}.vcf.gz; done
for i in {1..22}; do rm LEAP_file_chr${i}.vcf | rm LEAP_file_chr${i}.log; done
```

## PARIS

First, seperate the files into the two files based on their arrays.
```R
info = fread("HOE-HO5_PARIS_For-Varun_Info (1)")
keep1 = subset(info, Techno_version == "HumanOmni5Exome-4v1-1_A")
keep2 = subset(info, Techno_version == "humanomni5exome-4v1_a")

keep1 = keep1[,2:3]
keep2 = keep2[,2:3]

write.table(keep1, file = "pariskeepexome1.txt", row.names = F, col.names = F, quote = F)
write.table(keep2, file = "pariskeepexome2.txt", row.names = F, col.names = F, quote = F)


```

Next, seperate the two files in bash

```bash
./plink --bfile HOE-HO5_PARIS_For-Varun --keep pariskeepexome1.txt --make-bed --out paris_exome1
./plink --bfile HOE-HO5_PARIS_For-Varun --keep pariskeepexome2.txt --make-bed --out paris_exome2

```

Now QC and all of that

```bash
./plink --bfile paris_exome1  --geno 0.1 --mind 0.1 --hwe 0.000001 --me 0.05 0.1  --make-bed --out QC1outputparisexome1
./plink --bfile QC1outputparisexome1 --het --out QC1checkhetparisexome1
./plink --bfile QC1outputparisexome1 --check-sex --out QC1checksexparisexome1




```

```R
sex = read.delim("QC1checksexparisexome1.sexcheck", sep = "")
head(sex)
sexproblem = subset(sex, STATUS == "PROBLEM")


het = read.delim("QC1checkhetparisexome1.het", sep = "")
het$HET = (het$N.NM. - het$O.HOM.)/het$N.NM. #create heterozygosity stats

mean = mean(het$HET)
sd = sd(het$HET)
het$Z = (het$HET - mean)/sd #create Z scores of heterozygosity

hetoutlier = subset(het, abs(Z) > 3)

het2 = hetoutlier[,c(1:2)]
sex2 = sexproblem[,c(1:2)]
failedsample = rbind(het2, sex2)

write.table(failedsample, file = "failedsampleparisexome1.txt", row.names = F, col.names = T, quote = F)
```

```bash
./plink --bfile QC1outputparisexome1 --remove failedsampleparisexome1.txt --make-bed --out QC2outputparisexome1

./plink --bfile QC1outputparisexome2 --remove failedsampleparisexome2.txt --make-bed --out QC2outputparisexome2
```
