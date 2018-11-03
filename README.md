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
