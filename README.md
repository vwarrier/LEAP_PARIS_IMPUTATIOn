# LEAP_PARIS_IMPUTATION

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




## PARIS
