#Polygenic scores in the LEAP and the Paris cohort
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
pca = fread("~/LEAP_PARIS/parisleap_pcaall.eigenvec")
setnames(pca, "V2", "IID")
merged = merge(merged, pca, by = "IID")
summary(lm(scale(RBS_total) ~ scale(`1.000000`) + as.factor(sex) + scale(V3) + scale(V4) + scale(V5) + scale(V6) + scale(V7) + scale(V8) + scale(V9) + scale(V10) + scale(V11) + scale(V12) , data = merged))
