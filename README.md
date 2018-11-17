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


## PARIS OmniExpress and PARIS Exome

First, seperate the files and create the files for updating the parents

```R
info = fread("HOE-HO5_PARIS_For-Varun_Info (1)")
keep_exome = subset(info, Techno == "HumanOmni5Exome")
keep_exome = keep_exome[,2:3]
write.table(keep_exome, file = "pariskeepexome.txt", row.names = F, col.names = F, quote = F)


express = subset(info, Techno == "humanOmniExpress")
express_updateID = express[,c("Barcodes", "Barcodes", "FID", "IID")]
express_updatesex = express[,c("FID", "IID", "sex")]
express_updateparents = express[,c("FID", "IID", "FatherID", "MotherID")]

write.table(express_updateID, file = "parisexpressupdateID.txt", row.names = F, col.names = F, quote = F)
write.table(express_updatesex, file = "parisexpressupdatesex.txt", row.names = F, col.names = F, quote = F)
write.table(express_updateparents, file = "parisexpressupdateparents.txt", row.names = F, col.names = F, quote = F)
```

In bash, create the exome file

```bash
./plink --bfile HOE-HO5_PARIS_For-Varun --keep pariskeepexome.txt --make-bed --out paris_exome
```

Now follow the QC pipeline for LEAP, and create the imputation files. Remember to rename the files at the last stage. This creates the files for the Exome. Remove it to a seperate file and proceed to the OmniExpress

```bash
./plink --bfile  Omniexpres_559ind_ForVarun --update-ids parisexpressupdateID.txt --make-bed --out parisexpress
./plink --bfile  parisexpress --update-sex parisexpressupdatesex.txt --update-parents parisexpressupdateparents.txt --make-bed --out parisexpress
```

Now follow the QC pipeline for LEAP, and create the imputation files. Remember to rename the files at the last stage. This creates the files for the Express. Remove it to a seperate file.

## Post imputation

Follow some of the commands from SPARK, and for remaining commands, follow this. 

```bash

for i in {1..22};do ./plink --bfile ./PARISexome_files/Imputed/PARISexome_chr${i}  --make-bed --biallelic-only --maf 0.05 --geno 0.05 --hwe 0.000001 --out ./PARISexome_files/Imputed/PARISexomeround1_chr${i}; done
for i in {1..22};do ./plink --bfile ./PARISexpress_files/Imputed/PARISexome_chr${i}  --make-bed --biallelic-only --maf 0.05 --geno 0.05 --hwe 0.000001 --out ./PARISexpress_files/Imputed/PARISexpressround1_chr${i}; done
for i in {1..22};do ./plink --bfile ./LEAP_files/Imputed/LEAP_chr${i}  --make-bed --biallelic-only --maf 0.05 --geno 0.05 --hwe 0.000001 --out ./LEAP_files/Imputed/LEAPround1_chr${i}; done



./plink --bfile ./PARISexome_files/Imputed/PARISexomeround1_chr22 --merge-list ./PARISexome_files/Imputed/mergelist.txt --make-bed  --out ./PARISexome_files/Imputed/Parisexome_merged
./plink --bfile ./PARISexpress_files/Imputed/PARISexpressround1_chr22 --merge-list ./PARISexpress_files/Imputed/mergelist.txt --make-bed  --out ./PARISexpress_files/Imputed/Parisexpress_merged
./plink --bfile ./LEAP_files/Imputed/LEAPround1_chr22 --merge-list ./LEAP_files/Imputed/mergelist.txt --make-bed  --out ./LEAP_files/Imputed/LEAP_merged


for i in {1..22}; do ./plink --bfile ./PARISexome_files/Imputed/PARISexomeround1_chr${i} --exclude ./PARISexome_files/Imputed/Parisexome_merged-merge.missnp --make-bed --out ./PARISexome_files/Imputed/PARISexomeround1_chr${i}; done
for i in {1..22}; do ./plink --bfile ./PARISexpress_files/Imputed/PARISexpressround1_chr${i} --exclude ./PARISexpress_files/Imputed/Parisexpress_merged-merge.missnp --make-bed --out ./PARISexpress_files/Imputed/PARISexpressround1_chr${i}; done
for i in {1..22}; do ./plink --bfile ./LEAP_files/Imputed/LEAPround1_chr${i} --exclude ./LEAP_files/Imputed/LEAP_merged-merge.missnp --make-bed --out ./LEAP_files/Imputed/LEAPround1_chr${i}; done


./plink --bfile ./PARISexome_files/Imputed/PARISexomeround1_chr22 --merge-list ./PARISexome_files/Imputed/mergelist.txt --make-bed  --out ./PARISexome_files/Imputed/Parisexome_merged2
./plink --bfile ./PARISexpress_files/Imputed/PARISexpressround1_chr22 --merge-list ./PARISexpress_files/Imputed/mergelist.txt --make-bed  --out ./PARISexpress_files/Imputed/Parisexpress_merged2
./plink --bfile ./LEAP_files/Imputed/LEAPround1_chr22 --merge-list ./LEAP_files/Imputed/mergelist.txt --make-bed  --out ./LEAP_files/Imputed/LEAP_merged2


./plink --bfile ./PARISexome_files/Imputed/Parisexome_merged2 --maf 0.01 --update-name ~/SFARI/liftOverPlink/plinkrecodingfile.txt --hwe 0.000001 --geno 0.05 --mind 0.05 --make-bed --out ParisexomemergedQC
./plink --bfile ./PARISexpress_files/Imputed/Parisexpress_merged2 --maf 0.01 --update-name ~/SFARI/liftOverPlink/plinkrecodingfile.txt --hwe 0.000001 --geno 0.05 --mind 0.05 --make-bed --out ParisexpressmergedQC
./plink --bfile ./LEAP_files/Imputed/LEAP_merged2 --maf 0.01 --update-name ~/SFARI/liftOverPlink/plinkrecodingfile.txt --hwe 0.000001 --geno 0.05 --mind 0.05 --make-bed --out LEAPmergedQC
```

## Update the fam files

```{R}
library(data.table)
library(tidyr)

LEAP_fam = fread("LEAPmergedQC.fam")
Express_fam = fread("ParisexpressmergedQC.fam")
Exome_fam = fread("ParisexomemergedQC.fam")

Exome_old = fread("/mnt/b2/home4/arc/vw260/LEAP_PARIS/PARISexome_files/QC2output.fam")
Express_old = fread("/mnt/b2/home4/arc/vw260/LEAP_PARIS/QC2output.fam")
LEAP_old = fread("/mnt/b2/home4/arc/vw260/LEAP_PARIS/LEAP_files/QC2output.fam")

LEAP_fam = LEAP_fam %>% separate(V2, into = c('FID', 'IID'), sep = 6)
Express_fam = Express_fam %>% separate(V2, into = c('FID', 'IID'), sep = 6)
Exome_fam = Exome_fam %>% separate(V2, into = c('FID', 'IID'), sep = 6)

foo <- data.frame(do.call('rbind', strsplit(as.character(LEAP_fam$V2),'_',fixed=TRUE)))
LEAP_fam = cbind(LEAP_fam, foo)

foo <- data.frame(do.call('rbind', strsplit(as.character(Express_fam$V2),'_',fixed=TRUE)))
Express_fam = cbind(Express_fam, foo)

foo <- data.frame(do.call('rbind', strsplit(as.character(Exome_fam$V2),'_',fixed=TRUE)))
Exome_fam = cbind(Exome_fam, foo)

Express_all = merge(Express_fam, Express_old, by.x = "X2", by.y = "V2")
Exome_all = merge(Exome_fam, Exome_old, by.x = "X2", by.y = "V2")
LEAP_all = merge(LEAP_fam, LEAP_old, by.x = "X2", by.y = "V2")

Express_updatenames = Express_all[,c("V1.x", "V2", "X1", "X2")]
Exome_updatenames = Exome_all[,c("V1.x", "V2", "X1", "X2")]
LEAP_updatenames = LEAP_all[,c("V1.x", "V2", "X1", "X2")]

write.table(Express_updatenames, file = "Express_updatenames.txt", row.names = F, col.names = F, quote = F)
write.table(Exome_updatenames, file = "Exome_updatenames.txt", row.names = F, col.names = F, quote = F)
write.table(LEAP_updatenames, file = "LEAP_updatenames.txt", row.names = F, col.names = F, quote = F)

Express_updatesex = Express_all[,c("X1", "X2", "V5.y")]
Exome_updatesex = Exome_all[,c("X1", "X2", "V5.y")]
LEAP_updatesex = LEAP_all[,c("X1", "X2", "V5.y")]

Express_updatepheno = Express_all[,c("X1", "X2", "V6.y")]
Exome_updatepheno = Exome_all[,c("X1", "X2", "V6.y")]
LEAP_updatepheno = LEAP_all[,c("X1", "X2", "V6.y")]

Express_updateparents = Express_all[,c("X1", "X2", "V3.y", "V4.y")]
Exome_updateparents = Exome_all[,c("X1", "X2", "V3.y", "V4.y")]

write.table(Express_updatesex, file = "Express_updatesex.txt", row.names = F, col.names = F, quote = F)
write.table(Exome_updatesex, file = "Exome_updatesex.txt", row.names = F, col.names = F, quote = F)
write.table(LEAP_updatesex, file = "LEAP_updatesex.txt", row.names = F, col.names = F, quote = F)

write.table(Express_updatepheno, file = "Express_updatepheno.txt", row.names = F, col.names = F, quote = F)
write.table(Exome_updatepheno, file = "Exome_updatepheno.txt", row.names = F, col.names = F, quote = F)
write.table(LEAP_updatepheno, file = "LEAP_updatepheno.txt", row.names = F, col.names = F, quote = F)

write.table(Express_updateparents, file = "Express_updateparents.txt", row.names = F, col.names = F, quote = F)
write.table(Exome_updateparents, file = "Exome_updateparents.txt", row.names = F, col.names = F, quote = F)

LEAP_pheno = fread("./LEAP_files/LEAP_Clinical-Info_forVarun_PSC2_short (1)")
LEAP_fam = fread("LEAPmergedQC2.fam")

LEAP_updatefam = LEAP_pheno[,c("Barcodes", "Barcodes", "FID_PSC2", "IID_PSC2")]
LEAP_updateparents = LEAP_pheno[,c("FID_PSC2", "IID_PSC2", "PSC2_father", "PSC2_mother")]

write.table(LEAP_updatefam, file = "LEAPupdatefam.txt", row.names = F, col.names = F, quote = F)
write.table(LEAP_updateparents, file = "LEAPupdateparents.txt", row.names = F, col.names = F, quote = F)
```

Now do what you need to do in Plink

```bash

./plink --bfile ParisexpressmergedQC --update-ids Express_updatenames.txt --make-bed  --out ParisexpressmergedQC2
./plink --bfile ParisexomemergedQC --update-ids Exome_updatenames.txt --make-bed  --out ParisexomemergedQC2
./plink --bfile LEAPmergedQC --update-ids LEAP_updatenames.txt --make-bed  --out LEAPmergedQC2


./plink --bfile LEAPmergedQC2  --update-sex LEAP_updatesex.txt --pheno LEAP_updatepheno.txt --make-bed  --out LEAPmergedQC2
 ./plink --bfile ParisexomemergedQC2 --update-parents Exome_updateparents.txt --update-sex Exome_updatesex.txt --pheno Exome_updatepheno.txt --make-bed  --out ParisexomemergedQC2
 ./plink --bfile ParisexpressmergedQC2 --update-parents Express_updateparents.txt --update-sex Express_updatesex.txt --pheno Express_updatepheno.txt --make-bed  --out ParisexpressmergedQC2
 ./plink --bfile LEAPmergedQC2 --update-ids LEAPupdatefam.txt --make-bed --out LEAPmergedQC3 
 ./plink --bfile LEAPmergedQC3 --update-parents LEAPupdateparents.txt --make-bed --out LEAPmergedQC3
```

## Create files with the pheno
```{R}
setwd("/mnt/b2/home4/arc/vw260/LEAP_PARIS/")
library(data.table)
paris_pheno = fread("PARIS_all_pheno.txt")
Parisexpress_fam = fread("ParisexpressmergedQC2.fam")
Parisexome_fam = fread("ParisexomemergedQC2.fam")
paris_pheno$array = ifelse(paris_pheno$IID %in% Parisexpress_fam$V2, "Express", "NA")
paris_pheno$array = ifelse(paris_pheno$IID %in% Parisexome_fam$V2, "Exome", paris_pheno$array)


LEAP_pheno = fread("./LEAP_files/LEAP_Clinical-Info_forVarun_PSC2_short (1)")
LEAP_fam = fread("LEAPmergedQC3.fam")

LEAP_pheno_2 = subset(LEAP_pheno, RBS_total > 0 & RBS_total < 100)
LEAP_keep_forGRM = LEAP_pheno_2[,c(1:2)]
write.table(LEAP_keep_forGRM, file = "LEAP_keep_forGRM.txt", row.names = F, col.names = F, quote = F

            
paris_express_pheno = subset(paris_pheno, array == "Express")
paris_express_keep = subset(paris_express_pheno, RBSR_NONCOMPLET == "0")
write.table(paris_express_keep[,c("FID", "IID")], file = "PARISexpress_keep_forGRM.txt", row.names = F, col.names = F, quote = F)


paris_exome_pheno = subset(paris_pheno, array == "Exome")
paris_exome_keep = subset(paris_exome_pheno, RBSR_NONCOMPLET == "0")
write.table(paris_exome_keep[,c("FID", "IID")], file = "PARISexome_keep_forGRM.txt", row.names = F, col.names = F, quote = F)

```

Now update in plink, and merge the files
```bash
./plink --bfile LEAPmergedQC3 --keep LEAP_keep_forGRM.txt --make-bed --out LEAPwithpheno
./plink --bfile ParisexpressmergedQC2 --keep PARISexpress_keep_forGRM.txt --make-bed --out Parisexpresswithpheno
./plink --bfile ParisexomemergedQC2 --keep PARISexome_keep_forGRM.txt --make-bed --out Parisexomewithpheno

./plink --bfile LEAPwithpheno --merge-list mergeparisfileswithpheno.txt --make-bed --out parisallwithpheno

./plink --bfile LEAPwithpheno --exclude parisallwithpheno-merge.missnp --make-bed --out LEAPwithpheno
./plink --bfile Parisexpresswithpheno --exclude parisallwithpheno-merge.missnp --make-bed --out Parisexpresswithpheno
./plink --bfile Parisexomewithpheno --exclude parisallwithpheno-merge.missnp --make-bed --out Parisexomewithpheno

./plink --bfile LEAPwithpheno --merge-list mergeparisfileswithpheno.txt --make-bed --out parisallwithpheno

```

## Generate unrelated individuals and conduct PCA
```bash
./gcta64 --bfile ~/LEAP_PARIS/parisallwithpheno --maf 0.01 --make-grm --out PARISLEAPgrm --thread-num 10
./gcta64 --grm PARISLEAPgrm --grm-cutoff 0.05 --make-grm --out PARISLEAPgrmunrelated

./plink --bfile parisallwithpheno --keep ~/GCTA/PARISLEAPgrmunrelated.grm.id --make-bed --out parisallwithphenounrelated

./plink --bfile ~/SFARI/liftOverPlink/hapmap3_hg19_eurfoundersonly --bmerge parisallwithphenounrelated --make-bed --out parisleaphapmap3forpc
./plink --bfile ~/SFARI/liftOverPlink/hapmap3_hg19_eurfoundersonly --exclude parisleaphapmap3forpc-merge.missnp --make-bed --out Hapmap3_formerging
./plink --bfile parisallwithphenounrelated  --exclude parisleaphapmap3forpc-merge.missnp --make-bed --out parisleapfileformerging
./plink --bfile Hapmap3_formerging --bmerge parisleapfileformerging --make-bed --out parisleaphapmap3forpc

./plink --bfile parisleaphapmap3forpc --geno 0.1 --hwe 0.000001 --make-bed --out parisleaphapmap3forpc

./plink --bfile parisleaphapmap3forpc --maf 0.05 --indep-pairwise 100 50 0.2 --out parisleaphapmap3forpcpruned 

./plink --bfile parisleaphapmap3forpc --exclude parisleaphapmap3forpcpruned.prune.out --pca --out parisleap_pcaall
```

Now read the PCs into R and conduct 

```{R}
library(data.table)
library(plyr)
pc = fread("parisleap_pcaall.eigenvec")
selectedRows <- pc[grep("NA", pc$V2), ] #all the hapmapsamples have FID with NA

#Calculate the mean and SD of PC1 and PC2 based on the hapmapsamples
meanV1 = mean(selectedRows$V3)
sdV1 = sd(selectedRows$V3)
meanV2 = mean(selectedRows$V4)
sdV2 = sd(selectedRows$V4)

pc$ZPC1 = (abs(pc$V3 - meanV1))/sdV1
pc$ZPC2 = (abs(pc$V4 - meanV2))/sdV2

selectedRows2 <- pc[!grep("NA", pc$V2), ] #Now restrict it to the SSC samples
PCOK = subset(selectedRows2, ZPC1 < 5 & ZPC2 < 5) #include only samples that are less than 5 SDs away from the mean

PCOK = PCOK[,c(1:2)]

write.table(PCOK, file = "keepfile.txt", row.names = F, col.names = F, quote = F)
```

Next, keep only individuals whose PCs are fine
```bash
./plink --bfile parisallwithphenounrelated --keep keepfile.txt --make-bed --out parisallphenounrelatedpcok 
```

## Merging files and conducting PGS
# Please also see the R script in the attached file.

```{R}
#Polygenic scores in the LEAP and the EU_AIMS cohort
library(data.table)
library(tidyr)
library(dplyr)
library(gtools)

PGS = fread("~/ALSPAC/PRSice2results/PARIS_LEAPSQprsice.all.score", header = T)
LEAP_phenotype = fread("~/LEAP_PARIS/LEAP_files/LEAP_Clinical-Info_forVarun_PSC2_short (1)")
IQ_LEAP = LEAP_phenotype[,c("IID_PSC2", "piq", "fsiq", "viq")]

IQ_LEAP <- IQ_LEAP %>%
  # separate column0 into three columns based on the location of the pipe ,
  tidyr::separate("piq", c('piq1','piq2'), sep ='\\,', remove = FALSE, fill="right") %>%
  # make a newcolumn based on values in column0. Take only letters following the last pipe ,
  dplyr::mutate(newcolumn=sub(pattern=".*\\,([0-9]*$)", replacement="\\1", "piq")) %>%
  # Make the NA's blanks
  dplyr::mutate_all(funs(ifelse(is.na(.),"",.)))      


IQ_LEAP <- IQ_LEAP %>%
  tidyr::separate("viq", c('viq1','viq2'), sep ='\\,', remove = FALSE, fill="right") %>%
  dplyr::mutate(newcolumn=sub(pattern=".*\\,([0-9]*$)", replacement="\\1", "viq")) %>%
  dplyr::mutate_all(funs(ifelse(is.na(.),"",.)))      

IQ_LEAP <- IQ_LEAP %>%
  tidyr::separate("fsiq", c('fsiq1','fsiq2'), sep ='\\,', remove = FALSE, fill="right") %>%
  dplyr::mutate(newcolumn=sub(pattern=".*\\,([0-9]*$)", replacement="\\1", "fsiq")) %>%
  dplyr::mutate_all(funs(ifelse(is.na(.),"",.)))      

IQ_LEAP2 = IQ_LEAP[,c("fsiq1", "piq1", "viq1")]
LEAP_phenotype = cbind(LEAP_phenotype, IQ_LEAP2)

setnames(LEAP_phenotype, 1, "FID")
setnames(LEAP_phenotype, 2, "IID")
setnames(LEAP_phenotype, 3, "Father")
setnames(LEAP_phenotype, 4, "Mother")

LEAP_phenotype2 = LEAP_phenotype[,c("FID", "IID", "Father", "Mother", "sex", "Pheno", "age", "RBS_total",
                                    "fsiq1", "piq1", "viq1", "group")]
LEAP_phenotype2$Category = "LEAP"
LEAP_phenotype2$age = (LEAP_phenotype2$age)/365

PARIS_phenotype = fread("~/LEAP_PARIS/PARIS_all_pheno.txt")
setnames(PARIS_phenotype, "FatherId", "Father")
setnames(PARIS_phenotype, "MotherId", "Mother")
setnames(PARIS_phenotype, "POP_TYP", "group")
setnames(PARIS_phenotype, "ASD", "Pheno")
setnames(PARIS_phenotype, "RBSR_6_TOTAL", "RBS_total")
setnames(PARIS_phenotype, "DIAG_AGE", "age")
setnames(PARIS_phenotype, "Verb_IQ", "viq1")
setnames(PARIS_phenotype, "NVerb_IQ", "piq1")
setnames(PARIS_phenotype, "TOTAL_IQ", "fsiq1")

PARIS_phenotype$Category = "PARIS"

fulldata = smartbind(PARIS_phenotype, LEAP_phenotype2)

fulldata$piq1[fulldata$piq1 == "999"] <- NA
fulldata$fsiq1[fulldata$fsiq1 == "999"] <- NA
fulldata$viq1[fulldata$viq1 == "999"] <- NA
fulldata$RBS_total[fulldata$RBS_total == "999"] <- NA

merged = merge(PGS, fulldata, by = "IID")
pca = fread("parisleap_pcaall.eigenvec")
setnames(pca, "V2", "IID")
merged = merge(merged, pca, by = "IID")
summary(lm(scale(RBS_total) ~ scale(`1.000000`) + as.factor(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) , data = merged))
```
For meta-analysis please see the script uploaded

## Combining CHOP, LEAP, and SSC, creating unrelated individuals and conducting bivariate GCTA-GREML

```{bash}
./plink --bfile CHOPQCunrelated --bmerge ~/LEAP_PARIS/parisallphenounrelatedpcok --make-bed --out CHOPLEAP
./plink --bfile ~/LEAP_PARIS/parisallphenounrelatedpcok --exclude CHOPLEAP-merge.missnp --make-bed --out LEAPformerge
./plink --bfile CHOPQCunrelated --exclude CHOPLEAP-merge.missnp --make-bed --out CHOPformerge
./plink --bfile CHOPformerge --bmerge LEAPformerge --make-bed --out CHOPLEAP

./plink --bfile CHOPLEAP --bmerge ~/SFARI/liftOverPlink/files_imputed/SFARImergedallcasesonly --make-bed --out CHOPLEAPSFARI

./gcta64 --bfile ~/AGRE_data/CHOP/imputed_plinkfile/CHOPLEAPSFARI --make-grm --out allgrm --thread-num 10 --maf 0.01
./gcta64 --grm Allgrm --grm-cutoff 0.05 --make-grm --out Allgrmunrelated --thread-num 20
```
