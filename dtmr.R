# 0: Env setting
# 1: Load exposure data (eQTL/pQTL)
# 2: MR on T2DM and HbA1c
# 3: MR on PD (discovery: IPDGC excluding UKB)
# 4: MR on PD (replication: FinnGen + UKB)
# 5: Additional MR Analysis
## 5-1: MR with Sex-specific IV using sex-specific pQTL results
## 5-2: MR with DM/HbA1c-modifying IV
## 5-3: MR on REM-sleep Behavior Disorder
# 6: Sensitivity analysis

### 0: Environment setting =====================================================
library(TwoSampleMR); library(MRPRESSO); library(RadialMR); library(dbparser); library(XML); library(ggplot2); library(tidyr); library(data.table); library(dplyr) 
setwd("r2_0.1")


### 1: Load exposure data ======================================================
# 1-1: Load eQTL data from the eQTLGen consortium ==============================
eqtl <- read.csv("../../../eQTLGen/eQTLGen.dpp4i.noPalin.maf.csv", header=TRUE)

# Convert Z-score to beta and SE
eqtl$beta <- eqtl$Zscore / (2 * sqrt(eqtl$MAF * (1 - eqtl$MAF) * (eqtl$NrSamples + eqtl$Zscore^2)))
eqtl$se <- 1 / (2 * sqrt(eqtl$MAF * (1 - eqtl$MAF) * (eqtl$NrSamples + eqtl$Zscore^2)))
eqtl$r2 <- 2*eqtl$beta^2*eqtl$MAF*(1-eqtl$MAF)

# Select and rename columns
eqtl <- eqtl[,c("SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","beta","se","gene_id","gene","MAF","NrSamples","Pvalue","r2")]
names(eqtl) <- c("SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","gene_id","gene","eaf.exposure","samplesize.exposure","pval.exposure","r2.exposure")
eqtl$"exposure" <- "DPP4"
eqtl$"id.exposure" <- "eQTLGen"
dim(eqtl) #177, 15

# Within +/-200 kb of the DPP4 gene (GRCh37: 162848751-162931052, +-200kb: 162648751-163131052)
eqtl <- eqtl[eqtl$pos.exposure >= 162848751-200000 & eqtl$pos.exposure <= 162931052+200000, ]
#eqtl <- eqtl[eqtl$pval.exposure <= 0.00001,]
dim(eqtl) #163, 15
eqtl <- eqtl[eqtl$pval.exposure <= 0.05/nrow(eqtl),]
dim(eqtl) #166, 15

# Clumping
eqtl_clumped <- clump_data(eqtl, clump_r2=0.1); nrow(eqtl_clumped) #8; 10


# 2-2: Load cis-pQTL data from the UK Biobank ==================================
pqtl <- read.csv("../../../pQTL/cispqtl.dpp4.csv",header=TRUE) # +-200kb region selected, palindromic SNPs removed

# Rename columns
pqtl$pval.exposure <- 10^-pqtl$LOG10P
pqtl <- pqtl[,c("SNP","CHROM","POS","ALLELE0","ALLELE1","A1FREQ",
                "N","BETA","SE","pval.exposure")]
names(pqtl) <- c("SNP","chr.exposure","pos.exposure","other_allele.exposure","effect_allele.exposure","eaf.exposure",
                  "samplesize.exposure","beta.exposure","se.exposure","pval.exposure")
pqtl$exposure <- "DPP4"
pqtl$id.exposure <- "cis-pQTL"
pqtl$r2.exposure <- (2*pqtl$eaf.exposure*(1-pqtl$eaf.exposure)*pqtl$beta.exposure^2)/
  ((2*pqtl$eaf.exposure*(1-pqtl$eaf.exposure)*pqtl$beta.exposure^2)+(2*pqtl$eaf.exposure*(1-pqtl$eaf.exposure)*pqtl$samplesize.exposure*pqtl$se.exposure^2))

# Get SNPs with significant p-value
#pqtl <- pqtl[pqtl$pval.exposure<=0.00001,]; nrow(pqtl) #230
#nrow(pqtl[pqtl$pval.exposure<=0.05,]) #485
#nrow(pqtl[pqtl$pval.exposure<=0.01,]) #443
pqtl <- pqtl[pqtl$pval.exposure<=0.05/nrow(pqtl),]; nrow(pqtl) #249

# Clumping
pqtl_clumped <- clump_data(pqtl, clump_r2=0.1); nrow(pqtl_clumped) #25; bon29


### 2. MR analyses on T2DM and HbA1c ===========================================
# 2-0: Make empty lists to save the results =========================================
dm_mr_list <- list()
dm_pleio_list <- list()
dm_hetero_list <- list()
dm_single_list <- list()

hb_mr_list <- list()
hb_pleio_list <- list()
hb_hetero_list <- list()
hb_single_list <- list()

# 2-1: Using eQTL data as exposure =============================================
## T2DM
dm_all2 <- read.table("../../../gwas/DIAMANTE-EUR.sumstat.txt", sep=" ", header=T) 
dm_all2 <- dm_all2[,c("rsID","chromosome.b37.","position.b37.","Fixed.effects_beta","Fixed.effects_SE",
                    "effect_allele","other_allele","effect_allele_frequency","Fixed.effects_p.value")]
names(dm_all2) <- c("SNP","chr.outcome","pos.outcome","beta.outcome","se.outcome",
                   "effect_allele.outcome","other_allele.outcome","eaf.outcome","pval.outcome")
dm_all2$id.outcome <- "T2DM (DIAMANTE)"; dm_all2$outcome <- "T2DM (DIAMANTE)"
dm_all2$samplesize.outcome <- 80154 + 853816

eqtl_dm <- harmonise_data (eqtl_clumped, eqtl_dm, action=3)
eqtl_dm <- eqtl_dm[eqtl_dm$mr_keep==TRUE,]

dm_mr_list[['eQTL']] <- mr(eqtl_dm, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=eqtl_dm)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(eqtl_dm$id.exposure), id.outcome = unique(eqtl_dm$id.outcome), outcome = unique(eqtl_dm$outcome), exposure = unique(eqtl_dm$exposure),
  method = method.presso, nsnp = nrow(eqtl_dm), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
dm_mr_list[['eQTL']] <- rbind(dm_mr_list$eQTL, df.presso)
dm_hetero_list[['eQTL']] <- mr_heterogeneity(eqtl_dm)
dm_pleio_list[['eQTL']] <- mr_pleiotropy_test(eqtl_dm)
dm_single_list[['eQTL']] <- mr_singlesnp(eqtl_dm, all_method = c("mr_ivw"))

dm_hetero_list #0.1712581

## HbA1c
hb_all2 <- read.table("../../gwas/MAGIC1000G_HbA1c_EUR.tsv", sep="\t", header=T) 
hb_all2 <- hb_all2[,c("variant","chromosome","base_pair_location","beta","standard_error",
                      "effect_allele","other_allele","effect_allele_frequency","p_value","sample_size")]
names(hb_all2) <- c("SNP","chr.outcome","pos.outcome","beta.outcome","se.outcome",
                    "effect_allele.outcome","other_allele.outcome","eaf.outcome","pval.outcome","samplesize.outcome")
hb_all2$id.outcome <- "HbA1c (MAGIC)"; hb_all2$outcome <- "HbA1c (MAGIC)"
eqtl_hb <- harmonise_data (eqtl_clumped, hb_all2, action=3)
eqtl_hb <- eqtl_hb[eqtl_hb$mr_keep==TRUE,]

hb_mr_list[['eQTL']] <- mr(eqtl_hb, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=eqtl_hb)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(eqtl_hb$id.exposure), id.outcome = unique(eqtl_hb$id.outcome), outcome = unique(eqtl_hb$outcome), exposure = unique(eqtl_hb$exposure),
  method = method.presso, nsnp = nrow(eqtl_hb), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
hb_mr_list[['eQTL']] <- rbind(hb_mr_list$eQTL, df.presso)
hb_hetero_list[['eQTL']] <- mr_heterogeneity(eqtl_hb)
hb_pleio_list[['eQTL']] <- mr_pleiotropy_test(eqtl_hb)
hb_single_list[['eQTL']] <- mr_singlesnp(eqtl_hb, all_method = c("mr_ivw"))

hb_hetero_list #eQTL: 0.6001620


# 2-2: Using cis-pQTL data as exposure =========================================
## T2DM
pqtl_dm <- harmonise_data (pqtl_clumped, dm_all2, action=2)
pqtl_dm <- pqtl_dm[pqtl_dm$mr_keep==TRUE,]

dm_mr_list[['pQTL']] <- mr(pqtl_dm, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_dm)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_dm$id.exposure), id.outcome = unique(pqtl_dm$id.outcome), outcome = unique(pqtl_dm$outcome), exposure = unique(pqtl_dm$exposure),
  method = method.presso, nsnp = nrow(pqtl_dm), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
dm_mr_list[['pQTL']] <- rbind(dm_mr_list$pQTL, df.presso)
dm_hetero_list[['pQTL']] <- mr_heterogeneity(pqtl_dm)
dm_pleio_list[['pQTL']] <- mr_pleiotropy_test(pqtl_dm)
dm_single_list[['pQTL']] <- mr_singlesnp(pqtl_dm, all_method = c("mr_ivw"))

dm_hetero_list
# pQTL: 0.08135981

# male
dmM <- read.table("../../../gwas/DIAGRAM.Morris2012.SexSpecific.2016JUL05/DIAGRAM.Morris2012.males.txt", sep="\t", header=T)
dmM$beta <- log(dmM$OR.EFFECT_ALLELE.)
dmM$se <- (log(dmM$OR_95U) - log(dmM$OR_95L)) / (2 * 1.96)
dmM <- dmM[,c("ID","CHROMOSOME","POSITION","beta","se",
              "EFFECT_ALLELE","OTHER_ALLELE","P_VALUE")]
names(dmM) <- c("SNP","chr.outcome","pos.outcome","beta.outcome","se.outcome",
                "effect_allele.outcome","other_allele.outcome","pval.outcome")
dmM$samplesize.outcome <- 20219 + 54604
dmM$eaf.outcome <- 0.5
dmM$id.outcome <- "dmM"
dmM$outcome <- "dmM"

pqtl_dmM <- harmonise_data (pqtl_clumped, dmM, action=3)
pqtl_dmM <- pqtl_dmM[pqtl_dmM$mr_keep==TRUE,]

dm_mr_list[['pQTL_dmM']] <- mr(pqtl_dmM, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_dmM)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_dmM$id.exposure), id.outcome = unique(pqtl_dmM$id.outcome), outcome = unique(pqtl_dmM$outcome), exposure = unique(pqtl_dmM$exposure),
  method = method.presso, nsnp = nrow(pqtl_dmM), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
dm_mr_list[['pQTL_dmM']] <- rbind(dm_mr_list$pQTL_dmM, df.presso)
dm_hetero_list[['pQTL_dmM']] <- mr_heterogeneity(pqtl_dmM)
dm_pleio_list[['pQTL_dmM']] <- mr_pleiotropy_test(pqtl_dmM)
dm_single_list[['pQTL_dmM']] <- mr_singlesnp(pqtl_dmM, all_method = c("mr_ivw"))

dm_hetero_list

## HbA1c
pqtl_hb <- harmonise_data (pqtl_clumped, hb_all2, action=2)
pqtl_hb <- pqtl_hb[pqtl_hb$mr_keep==TRUE,]

hb_mr_list[['pQTL']] <- mr(pqtl_hb, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_hb)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_hb$id.exposure), id.outcome = unique(pqtl_hb$id.outcome), outcome = unique(pqtl_hb$outcome), exposure = unique(pqtl_hb$exposure),
  method = method.presso, nsnp = nrow(pqtl_hb), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
hb_mr_list[['pQTL']] <- rbind(hb_mr_list$pQTL, df.presso)
hb_hetero_list[['pQTL']] <- mr_heterogeneity(pqtl_hb)
hb_pleio_list[['pQTL']] <- mr_pleiotropy_test(pqtl_hb)
hb_single_list[['pQTL']] <- mr_singlesnp(pqtl_hb, all_method = c("mr_ivw"))

hb_hetero_list #pQTL: 0.2326634

write.csv(do.call(rbind, dm_mr_list), file="dm_mr.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, dm_pleio_list), file="dm_pleio.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, dm_hetero_list), file="dm_hetero.csv", row.names = FALSE, na = "")

write.csv(do.call(rbind, hb_mr_list), file="hb_mr.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, hb_pleio_list), file="hb_pleio.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, hb_hetero_list), file="hb_hetero.csv", row.names = FALSE, na = "")


### 3: MR on PD ================================================================
# 3-0: Make a list to save the results =========================================
pd_mr_list <- list()
pd_presso_list <- list()
pd_pleio_list <- list()
pd_hetero_list <- list()

pdM_mr_list <- list()
pdM_presso_list <- list()
pdM_pleio_list <- list()
pdM_hetero_list <- list()

pdF_mr_list <- list()
pdF_presso_list <- list()
pdF_pleio_list <- list()
pdF_hetero_list <- list()

# 3-1: eQTL on PD ==============================================================
## All
pd <- read.table("../../../gwas/GCST009324.h.tsv",sep='\t',header=T)
pd <- pd[c("rsid","chromosome","base_pair_location","effect_allele","other_allele",
           "effect_allele_frequency","beta","standard_error","p_value")]
names(pd) <- c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome",
               "eaf.outcome","beta.outcome","se.outcome","pval.outcome")
pd$outcome <- "PD"
pd$id.outcome <- "GCST009324"
pd$samplesize.outcome <- 37688 + 981372

eqtl_in_pd <- eqtl %>%
  filter(SNP %in% pd$SNP)
eqtl_in_pd_clumped <- clump_data(eqtl_in_pd, clump_r2=0.1); nrow(eqtl_in_pd_clumped)
eqtl_pd_dat <- pd %>%
  filter(SNP %in% eqtl_in_pd_clumped$SNP)
nrow(eqtl_pd_dat) #8
eqtl_pd <- harmonise_data (eqtl_in_pd_clumped, eqtl_pd_dat, action = 2)
eqtl_pd <- eqtl_pd[eqtl_pd$mr_keep==TRUE,]

# MR
pd_mr_list[['eQTL']] <- mr(eqtl_pd, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=eqtl_pd)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(eqtl_pd$id.exposure), id.outcome = unique(eqtl_pd$id.outcome), outcome = unique(eqtl_pd$outcome), exposure = unique(eqtl_pd$exposure),
  method = method.presso, nsnp = nrow(eqtl_pd), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pd_mr_list[['eQTL']] <- rbind(pd_mr_list$eQTL, df.presso)
pd_hetero_list[['eQTL']] <- mr_heterogeneity(eqtl_pd)
pd_pleio_list[['eQTL']] <- mr_pleiotropy_test(eqtl_pd)

# check
pd_hetero_list #0.2991976


## Male
pdM <- read.table("../../../gwas/Sumstats_sexSp/MALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz",sep='\t',header=T)
pdM <- pdM[c("ID","Allele1","Allele2","Freq1","Effect","StdErr","P.value")]
names(pdM) <- c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome",
                    "beta.outcome","se.outcome","pval.outcome")
pdM$effect_allele.outcome <- toupper(pdM$effect_allele.outcome)
pdM$other_allele.outcome <- toupper(pdM$other_allele.outcome)
pdM$outcome <- "PD_male"
pdM$id.outcome <- "PD_male"
#pdM$samplesize.outcome <- 13020 + 89660
pdM$samplesize.outcome <- 12054 + 11999
#https://github.com/neurogenetics/Autosomal-sex-differences-PDv2

eqtl_pdM <- harmonise_data (eqtl_clumped, pdM, action=3); nrow(eqtl_pdM) #7

#eqtl_pdM <- pdM[pdM$SNP %in% eqtl_clumped$SNP,]; nrow(eqtl_pdM) #9
#rs12619850 (chr2:162814823) not in pdM, no proxy SNPs were present (r2>0.7)
#eqtl_in_pdM <- eqtl %>%
#  filter(SNP %in% pdM$SNP)
#eqtl_in_pdM_clumped <- clump_data(eqtl_in_pdM, clump_r2=0.1); nrow(eqtl_in_pdM_clumped)
#eqtl_pdM <- pdM %>%
#  filter(SNP %in% eqtl_in_pdM_clumped$SNP)
#eqtl_pdM <- harmonise_data (eqtl_in_pdM_clumped, eqtl_pdM, action=2); nrow(eqtl_pdM)
eqtl_pdM <- eqtl_pdM[eqtl_pdM$mr_keep==TRUE,]


# MR
pdM_mr_list[['eQTL']] <- mr(eqtl_pdM, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=eqtl_pdM)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(eqtl_pdM$id.exposure), id.outcome = unique(eqtl_pdM$id.outcome), outcome = unique(eqtl_pdM$outcome), exposure = unique(eqtl_pdM$exposure),
  method = method.presso, nsnp = nrow(eqtl_pdM), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdM_mr_list[['eQTL']] <- rbind(pdM_mr_list$eQTL, df.presso)
pdM_hetero_list[['eQTL']] <- mr_heterogeneity(eqtl_pdM)
pdM_pleio_list[['eQTL']] <- mr_pleiotropy_test(eqtl_pdM)

# check
pdM_hetero_list # eQTL: 0.4340809

## Female
pdF <- read.table("../../../gwas/Sumstats_sexSp/FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz",sep='\t',header=T)
pdF <- pdF[c("ID","Allele1","Allele2","Freq1","Effect","StdErr","P.value")]
names(pdF) <- c("SNP","effect_allele.outcome","other_allele.outcome","eaf.outcome","beta.outcome","se.outcome","pval.outcome")
pdF$effect_allele.outcome <- toupper(pdF$effect_allele.outcome)
pdF$other_allele.outcome <- toupper(pdF$other_allele.outcome)
pdF$outcome <- "PD_female"
pdF$id.outcome <- "PD_female"
#pdF$samplesize.outcome <- 7947 + 90662
pdF$samplesize.outcome <- 7384 + 12389

eqtl_pdF <- harmonise_data (eqtl_clumped, pdF, action=3); nrow(eqtl_pdF) #7

#eqtl_pdF <- pdF[pdF$SNP %in% eqtl_clumped$SNP,]; nrow(eqtl_pdF) #10
#eqtl_in_pdF <- eqtl %>%
#  filter(SNP %in% pdF$SNP)
#eqtl_in_pdF_clumped <- clump_data(eqtl_in_pdF, clump_r2=0.1); nrow(eqtl_in_pdF_clumped)
#eqtl_pdF <- pdF %>%
#  filter(SNP %in% eqtl_in_pdF_clumped$SNP)
#eqtl_pdF <- harmonise_data (eqtl_in_pdF_clumped, eqtl_pdF, action=2)
eqtl_pdF <- eqtl_pdF[eqtl_pdF$mr_keep==TRUE,]

# MR
pdF_mr_list[['eQTL']] <- mr(eqtl_pdF, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=eqtl_pdF)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(eqtl_pdF$id.exposure), id.outcome = unique(eqtl_pdF$id.outcome), outcome = unique(eqtl_pdF$outcome), exposure = unique(eqtl_pdF$exposure),
  method = method.presso, nsnp = nrow(eqtl_pdF), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdF_mr_list[['eQTL']] <- rbind(pdF_mr_list$eQTL, df.presso)
pdF_hetero_list[['eQTL']] <- mr_heterogeneity(eqtl_pdF)
pdF_pleio_list[['eQTL']] <- mr_pleiotropy_test(eqtl_pdF)

# check
pdF_hetero_list #0.7266726


# 3-2: cis-pQTL on PD ==========================================================
## All
pqtl_in_pd <- pqtl %>%
  filter(SNP %in% pd$SNP)
pqtl_in_pd_clumped <- clump_data(pqtl_in_pd, clump_r2=0.1); nrow(pqtl_in_pd_clumped)
pqtl_pd <- pd %>%
  filter(SNP %in% pqtl_in_pd_clumped$SNP)

pqtl_pd <- harmonise_data (pqtl_in_pd_clumped, pqtl_pd, action = 2); nrow(pqtl_pd) #25
pqtl_pd <- pqtl_pd[pqtl_pd$mr_keep==TRUE,]

# MR
pd_mr_list[['pQTL']] <- mr(pqtl_pd, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_pd)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_pd$id.exposure), id.outcome = unique(pqtl_pd$id.outcome), outcome = unique(pqtl_pd$outcome), exposure = unique(pqtl_pd$exposure),
  method = method.presso, nsnp = nrow(pqtl_pd), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pd_mr_list[['pQTL']] <- rbind(pd_mr_list$pQTL, df.presso)
pd_hetero_list[['pQTL']] <- mr_heterogeneity(pqtl_pd)
pd_pleio_list[['pQTL']] <- mr_pleiotropy_test(pqtl_pd)

# check
pd_hetero_list #0.2604085


## Male
pqtl_in_pdM <- pqtl %>%
  filter(SNP %in% pdM$SNP)
pqtl_in_pdM_clumped <- clump_data(pqtl_in_pdM, clump_r2=0.1); nrow(pqtl_in_pdM_clumped)
pqtl_pdM <- pdM %>%
  filter(SNP %in% pqtl_in_pdM_clumped$SNP)

pqtl_pdM <- harmonise_data (pqtl_in_pdM_clumped, pqtl_pdM, action=2); nrow(pqtl_pdM) #22
pqtl_pdM <- pqtl_pdM[pqtl_pdM$mr_keep==TRUE,]

# MR
pdM_mr_list[['pQTL']] <- mr(pqtl_pdM, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_pdM)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_pdM$id.exposure), id.outcome = unique(pqtl_pdM$id.outcome), outcome = unique(pqtl_pdM$outcome), exposure = unique(pqtl_pdM$exposure),
  method = method.presso, nsnp = nrow(pqtl_pdM), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdM_mr_list[['pQTL']] <- rbind(pdM_mr_list$pQTL, df.presso)
pdM_hetero_list[['pQTL']] <- mr_heterogeneity(pqtl_pdM)
pdM_pleio_list[['pQTL']] <- mr_pleiotropy_test(pqtl_pdM)

#check
pdM_hetero_list # pQTL: 0.5044280


## Female
pqtl_in_pdF <- pqtl %>%
  filter(SNP %in% pdF$SNP)
pqtl_in_pdF_clumped <- clump_data(pqtl_in_pdF, clump_r2=0.1); nrow(pqtl_in_pdF_clumped)
pqtl_pdF <- pdF %>%
  filter(SNP %in% pqtl_in_pdF_clumped$SNP)

pqtl_pdF <- harmonise_data (pqtl_in_pdF_clumped, pqtl_pdF, action=2); nrow(pqtl_pdF) #22
pqtl_pdF <- pqtl_pdF[pqtl_pdF$mr_keep==TRUE,]

# MR
pdF_mr_list[['pQTL']] <- mr(pqtl_pdF, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_pdF)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_pdF$id.exposure), id.outcome = unique(pqtl_pdF$id.outcome), outcome = unique(pqtl_pdF$outcome), exposure = unique(pqtl_pdF$exposure),
  method = method.presso, nsnp = nrow(pqtl_pdF), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdF_mr_list[['pQTL']] <- rbind(pdF_mr_list$pQTL, df.presso)
pdF_hetero_list[['pQTL']] <- mr_heterogeneity(pqtl_pdF)
pdF_pleio_list[['pQTL']] <- mr_pleiotropy_test(pqtl_pdF)

#check
pdF_hetero_list #0.8174754

write.csv(do.call(rbind, pd_mr_list), file="pd_mr.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, pd_pleio_list), file="pd_pleio.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, pd_hetero_list), file="pd_hetero.csv", row.names = FALSE, na = "")

write.csv(do.call(rbind, pdM_mr_list), file="pdM_mr.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, pdM_pleio_list), file="pdM_pleio.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, pdM_hetero_list), file="pdM_hetero.csv", row.names = FALSE, na = "")

write.csv(do.call(rbind, pdF_mr_list), file="pdF_mr.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, pdF_pleio_list), file="pdF_pleio.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, pdF_hetero_list), file="pdF_hetero.csv", row.names = FALSE, na = "")


### 4: MR on FinnGen + UKB (FU) ================================================
# 4-0: Make a list to save the results =========================================
fu_mr_list <- list()
fu_pleio_list <- list()
fu_hetero_list <- list()

# 4-1: eQTL on FinnGen + UKB (FU) ==============================================
fu <- read.table("../../../meta_analysis_ukbb_summary_stats_finngen_R12_G6_PARKINSON_meta_out.dpp4_1Mb.tsv",sep='\t',header=T)
fu_dat <- fu[c("rsid","CHR","POS","ALT","REF","all_meta_N",
               "all_inv_var_meta_beta","all_inv_var_meta_sebeta","all_inv_var_meta_p")]
names(fu_dat) <- c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome","samplesize.outcome",
                   "beta.outcome","se.outcome","pval.outcome")
fu_dat$outcome <- "PD"
fu_dat$id.outcome <- "PD_FinnGen+UKB"
fu_dat$eaf.outcome <- 0.5 # assign arbitraily as palindromic SNPs were already removed

qtl_fu <- fu_dat %>%
  filter(SNP %in% c(eqtl_clumped$SNP,pqtl_clumped$SNP))

eqtl_fu <- harmonise_data (eqtl_clumped, qtl_fu, action=2); nrow(eqtl_fu) #8
eqtl_fu <- eqtl_fu[eqtl_fu$mr_keep==TRUE,]
fu_mr_list[['eQTL']] <- mr(eqtl_fu, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=eqtl_fu)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(eqtl_fu$id.exposure), id.outcome = unique(eqtl_fu$id.outcome), outcome = unique(eqtl_fu$outcome), exposure = unique(eqtl_fu$exposure),
  method = method.presso, nsnp = nrow(eqtl_fu), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
fu_mr_list[['eQTL']] <- rbind(fu_mr_list$eQTL, df.presso)
fu_hetero_list[['eQTL']] <- mr_heterogeneity(eqtl_fu)
fu_pleio_list[['eQTL']] <- mr_pleiotropy_test(eqtl_fu)

fu_hetero_list #0.3108805

# 4-2: cis-pQTL on FinnGen + UKB (FU) ==========================================
pqtl_fu <- harmonise_data (pqtl_clumped, qtl_fu, action=2); nrow(pqtl_fu) #24
pqtl_fu <- pqtl_fu[pqtl_fu$mr_keep==TRUE,]

fu_mr_list[['pQTL']] <- mr(pqtl_fu, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_fu)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_fu$id.exposure), id.outcome = unique(pqtl_fu$id.outcome), outcome = unique(pqtl_fu$outcome), exposure = unique(pqtl_fu$exposure),
  method = method.presso, nsnp = nrow(pqtl_fu), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
fu_mr_list[['pQTL']] <- rbind(fu_mr_list$pQTL, df.presso)
fu_hetero_list[['pQTL']] <- mr_heterogeneity(pqtl_fu)
fu_pleio_list[['pQTL']] <- mr_pleiotropy_test(pqtl_fu)

fu_hetero_list #0.10313389

write.csv(do.call(rbind, fu_mr_list), file="FinnGen+UKB_mr.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, fu_pleio_list), file="FinnGen+UKB_pleio.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, fu_hetero_list), file="FinnGen+UKB_hetero.csv", row.names = FALSE, na = "")


### 5: Additional MR Analysis
## 5-1: MR with Sex-specific IV using sex-specific pQTL results
# Load sex-specific cis-pQTL data from the UK Biobank 
# Male
pqtl_male <- read.table("../../../pqtl_dpp4_male.dpp4.glm.linear", header=F)
colnames(pqtl_male) <- c("chr", "pos", "SNP", "other_allele.exposure", "effect_allele.exposure", "A1.exposure", "TEST", 
                         "samplesize.exposure", "beta.exposure", "se.exposure", "T_STAT", "pval.exposure", "ERRCODE")
pqtl_male <- pqtl_male %>% select(-ERRCODE)
pqtl_male <- pqtl_male %>% filter(!is.na(pval.exposure))
dim(pqtl_male) #112 12

# Remove palindromic SNPs
pqtl_male <- pqtl_male %>%
  filter(!(other_allele.exposure == "A" & effect_allele.exposure == "T") &
           !(other_allele.exposure == "T" & effect_allele.exposure == "A") &
           !(other_allele.exposure == "C" & effect_allele.exposure == "G") &
           !(other_allele.exposure == "G" & effect_allele.exposure == "C"))
dim(pqtl_male) #106 12

# Get AF
pqtl_frq_male <- read.table("../../../ukb_dpp4_frq_male.frq", header=TRUE, stringsAsFactors=FALSE)

# Merge pQTL and AF
pqtl_male <- pqtl_male %>%
  left_join(pqtl_frq_male, by = "SNP")
pqtl_male <- pqtl_male %>%
  mutate(eaf.exposure = ifelse(effect_allele.exposure == A1, MAF, 1 - MAF))
dim(pqtl_male) #106 18

pqtl_male <- pqtl_male[pqtl_male$pval.exposure<=0.05/nrow(pqtl_male),]; nrow(pqtl_male) #33

pqtl_male$exposure <- "DPP4"
pqtl_male$id.exposure <- "cis-pQTL_male"
pqtl_male$r2.exposure <- (2*pqtl_male$eaf.exposure*(1-pqtl_male$eaf.exposure)*pqtl_male$beta.exposure^2)/
  ((2*pqtl_male$eaf.exposure*(1-pqtl_male$eaf.exposure)*pqtl_male$beta.exposure^2)+(2*pqtl_male$eaf.exposure*(1-pqtl_male$eaf.exposure)*pqtl_male$samplesize.exposure*pqtl_male$se.exposure^2))

# Clumping
pqtl_male_in_pdM <- pqtl_male %>%
  filter(SNP %in% pdM$SNP)
pqtl_clumped_male <- clump_data(pqtl_male_in_pdM, clump_r2=0.1)
nrow(pqtl_clumped_male) #13

# MR on PD
pqtl_male_pdM <- harmonise_data(pqtl_clumped_male, pdM, action=3)
pdM_mr_list[['pQTL_male']] <- mr(pqtl_male_pdM, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_male_pdM)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_male_pdM$id.exposure), id.outcome = unique(pqtl_male_pdM$id.outcome), outcome = unique(pqtl_male_pdM$outcome), exposure = unique(pqtl_male_pdM$exposure),
  method = method.presso, nsnp = nrow(pqtl_male_pdM), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdM_mr_list[['pQTL_male']] <- rbind(pdM_mr_list$pQTL_male, df.presso)
pdM_hetero_list[['pQTL_male']] <- mr_heterogeneity(pqtl_male_pdM)


## Female
pqtl_female <- read.table("../../../pqtl_dpp4_female.dpp4.glm.linear", header=F)
colnames(pqtl_female) <- c("chr", "pos", "SNP", "other_allele.exposure", "effect_allele.exposure", "A1.exposure", "TEST", 
                           "samplesize.exposure", "beta.exposure", "se.exposure", "T_STAT", "pval.exposure", "ERRCODE")
pqtl_female <- pqtl_female %>% select(-ERRCODE)
pqtl_female <- pqtl_female %>% filter(!is.na(pval.exposure))
#dim(pqtl_female) #112 12

# Remove palindromic SNPs
pqtl_female <- pqtl_female %>%
  filter(!(other_allele.exposure == "A" & effect_allele.exposure == "T") &
           !(other_allele.exposure == "T" & effect_allele.exposure == "A") &
           !(other_allele.exposure == "C" & effect_allele.exposure == "G") &
           !(other_allele.exposure == "G" & effect_allele.exposure == "C"))
#dim(pqtl_female) #106 12

# Get AF
pqtl_frq_female <- read.table("../../../ukb_dpp4_frq_female.frq", header=TRUE, stringsAsFactors=FALSE)

# Merge pQTL and AF
pqtl_female <- pqtl_female %>%
  left_join(pqtl_frq_female, by = "SNP")
pqtl_female <- pqtl_female %>%
  mutate(eaf.exposure = ifelse(effect_allele.exposure == A1, MAF, 1 - MAF))
dim(pqtl_female) #20 20 

pqtl_female <- pqtl_female[pqtl_female$pval.exposure<=0.05/nrow(pqtl_female),]; nrow(pqtl_female) #28

pqtl_female$exposure <- "DPP4"
pqtl_female$id.exposure <- "cis-pQTL_female"
pqtl_female$r2.exposure <- (2*pqtl_female$eaf.exposure*(1-pqtl_female$eaf.exposure)*pqtl_female$beta.exposure^2)/
  ((2*pqtl_female$eaf.exposure*(1-pqtl_female$eaf.exposure)*pqtl_female$beta.exposure^2)+(2*pqtl_female$eaf.exposure*(1-pqtl_female$eaf.exposure)*pqtl_female$samplesize.exposure*pqtl_female$se.exposure^2))

# Clumping
pqtl_female_in_pdM <- pqtl_female %>%
  filter(SNP %in% pdF$SNP)
pqtl_clumped_female <- clump_data(pqtl_female_in_pdM, clump_r2=0.1)
nrow(pqtl_clumped_female) #10

# MR on PD
pqtl_female_pdF <- harmonise_data(pqtl_clumped_female, pdF, action=2)
pdF_mr_list[['pQTL_female']] <- mr(pqtl_female_pdF, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_female_pdF)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_female_pdF$id.exposure), id.outcome = unique(pqtl_female_pdF$id.outcome), outcome = unique(pqtl_female_pdF$outcome), exposure = unique(pqtl_female_pdF$exposure),
  method = method.presso, nsnp = nrow(pqtl_female_pdF), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdF_mr_list[['pQTL_female']] <- rbind(pdF_mr_list$pQTL_female, df.presso)
pdF_hetero_list[['pQTL_female']] <- mr_heterogeneity(pqtl_female_pdF)

write.csv(do.call(rbind, pdM_mr_list), file="pdM_mr.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, pdM_hetero_list), file="pdM_hetero.csv", row.names = FALSE, na = "")

write.csv(do.call(rbind, pdF_mr_list), file="pdF_mr.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, pdF_hetero_list), file="pdF_hetero.csv", row.names = FALSE, na = "")


## 5-2: MR with DM/HbA1c-modifying IV ==========================================
## DM > PD
dm_all <- read.table("../../../gwas/DIAMANTE-EUR.sumstat.txt", sep=" ", header=T) 
dm_all <- dm_all[,c("rsID","chromosome.b37.","position.b37.","Fixed.effects_beta","Fixed.effects_SE",
                    "effect_allele","other_allele","effect_allele_frequency","Fixed.effects_p.value")]
names(dm_all) <- c("SNP","chr.exposure","pos.exposure","beta.exposure","se.exposure",
                   "effect_allele.exposure","other_allele.exposure","eaf.exposure","pval.exposure")
dm_all$id.exposure <- "T2DM (DIAMANTE)"; dm_all$exposure <- "T2DM (DIAMANTE)"
dm_all$samplesize.exposure <- 80154 + 853816
dm_all_clumped <- clump_data(dm_all, clump_r2 = 0.001, clump_p1=5e-8)
dim(dm_all_clumped) #185, 12
dm_all_pd <- harmonise_data (dm_all_clumped, pd, action=3)
dm_all_pd <- dm_all_pd[dm_all_pd$mr_keep==TRUE,]

# MR
pd_mr_list[['dm']] <- mr(dm_all_pd, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=dm_all_pd)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(dm_all_pd$id.exposure), id.outcome = unique(dm_all_pd$id.outcome), outcome = unique(dm_all_pd$outcome), exposure = unique(dm_all_pd$exposure),
  method = method.presso, nsnp = nrow(dm_all_pd), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pd_mr_list[['dm']] <- rbind(pd_mr_list$dm, df.presso)


# Male
dmM <- read.table("../../../gwas/DIAGRAM.Morris2012.SexSpecific.2016JUL05/DIAGRAM.Morris2012.males.txt", sep="\t", header=T)
dmM$beta <- log(dmM$OR.EFFECT_ALLELE.)
dmM$se <- (log(dmM$OR_95U) - log(dmM$OR_95L)) / (2 * 1.96)
dmM <- dmM[,c("ID","CHROMOSOME","POSITION","beta","se",
              "EFFECT_ALLELE","OTHER_ALLELE","P_VALUE")]
names(dmM) <- c("SNP","chr.exposure","pos.exposure","beta.exposure","se.exposure",
                "effect_allele.exposure","other_allele.exposure","pval.exposure")
dmM$samplesize.exposure <- 20219 + 54604
dmM$eaf.exposure <- 0.5 #arbitrarily assiged as the panlindromic SNPs (AT or CG) were already removed

dmM_clumped <- clump_data(dmM, clump_r2 = 0.001, clump_p1=5e-8)
dim(dmM_clumped) #24, 11
dmM_clumped$exposure <- "Male T2DM"
dmM_clumped$id.exposure <- "Male T2DM"
#dmM_clumped$exposure <- "Male T2DM (E4)"
#dmM_clumped$id.exposure <- "Male T2DM (E4)"
dmM_pdM <- harmonise_data (dmM_clumped, pdM, action=3)
dmM_pdM <- dmM_pdM[dmM_pdM$mr_keep==TRUE,]
pdM_mr_list[['dmM']] <- mr(dmM_pdM, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=dmM_pdM)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(dmM_pdM$id.exposure), id.outcome = unique(dmM_pdM$id.outcome), outcome = unique(dmM_pdM$outcome), exposure = unique(dmM_pdM$exposure),
  method = method.presso, nsnp = nrow(dmM_pdM), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdM_mr_list[['dmM']] <- rbind(pdM_mr_list$dmM, df.presso)


# female
dmF <- read.table("../../../gwas/DIAGRAM.Morris2012.SexSpecific.2016JUL05/DIAGRAM.Morris2012.females.txt", sep="\t", header=T)
dmF$beta <- log(dmF$OR.EFFECT_ALLELE)
dmF$se <- (log(dmF$OR_95U) - log(dmF$OR_95L)) / (2 * 1.96)
dmF <- dmF[,c("ID","CHROMOSOME","POSITION","beta","se",
              "EFFECT_ALLELE","OTHER_ALLELE","P_VALUE")]
names(dmF) <- c("SNP","chr.exposure","pos.exposure","beta.exposure","se.exposure",
                "effect_allele.exposure","other_allele.exposure","pval.exposure")
dmF$samplesize.exposure <- 14621 + 60377
dmF$eaf.exposure <- 0.5 #arbitrarily assiged as the panlindromic SNPs (AT or CG) were already removed

dmF_clumped <- clump_data(dmF, clump_r2 = 0.001, clump_p1=5e-8)
dim(dmF_clumped) #14, 11
dmF_clumped$exposure <- "Female T2DM"
dmF_clumped$id.exposure <- "Female T2DM"
#dmF_clumped$exposure <- "Female T2DM (E4)"
#dmF_clumped$id.exposure <- "Female T2DM (E4)"
dmF_pdF <- harmonise_data (dmF_clumped, pdF, action=3)
dmF_pdF <- dmF_pdF[dmF_pdF$mr_keep==TRUE,]
pdF_mr_list[['dmF']] <- mr(dmF_pdF, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=dmF_pdF)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(dmF_pdF$id.exposure), id.outcome = unique(dmF_pdF$id.outcome), outcome = unique(dmF_pdF$outcome), exposure = unique(dmF_pdF$exposure),
  method = method.presso, nsnp = nrow(dmF_pdF), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdF_mr_list[['dmF']] <- rbind(pdF_mr_list$dmF, df.presso)


### DM risk-modifying
## eQTL
# MR with all SNPs in eQTL data
eqtl_dmSig <- harmonise_data(eqtl_clumped, dm_all2, action=3)
eqtl_dmSig <- eqtl_dmSig[eqtl_dmSig$mr_keep==TRUE,]

# Extract significantly dm-modifying variants
eqtl_dmSig_SNP <- mr_singlesnp(eqtl_dm, all_method = c("mr_ivw")) %>%
  filter(SNP != "All - Inverse variance weighted" & p < 0.05) %>%
  select(SNP)
eqtl_dmSig <- eqtl_clumped %>%
  filter(SNP %in% eqtl_dmSig_SNP$SNP); dim(eqtl_dmSig) #2
eqtl_dmSig_pd <- harmonise_data (eqtl_dmSig, pd, action=2)

# MR
eqtl_dmSig_pd <- eqtl_dmSig_pd[eqtl_dmSig_pd$mr_keep==TRUE,]
eqtl_dmSig_pd$id.exposure <- "eQTLGen_dmSig"
pd_mr_list[['eQTL_dmSig']] <- mr(eqtl_dmSig_pd, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))

# Male
eqtl_dmSig_pdM <- harmonise_data (eqtl_dmSig, pdM, action=2) #from eqtl_dm
eqtl_dmSig_pdM <- eqtl_dmSig_pdM[eqtl_dmSig_pdM$mr_keep==TRUE,]
eqtl_dmSig_pdM$id.exposure <- "eQTLGen_dmSig"
pdM_mr_list[['eQTL_dmSig']] <- mr(eqtl_dmSig_pdM, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
pdM_mr_list[['eQTL_dmSig']] <- mr(eqtl_dmSig_pdM) #from eqtl_dm

# Female
eqtl_dmSig_pdF <- harmonise_data (eqtl_dmSig, pdF, action=2) #from eqtl_dm
eqtl_dmSig_pdF <- eqtl_dmSig_pdF[eqtl_dmSig_pdF$mr_keep==TRUE,]
eqtl_dmSig_pdF$id.exposure <- "eQTLGen_dmSig"
pdF_mr_list[['eQTL_dmSig']] <- mr(eqtl_dmSig_pdF, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))


## pQTL ######
# MR with all SNPs in pQTL data
pqtl_dmSig <- extract_outcome_data(pqtl$SNP, outcomes="ebi-a-GCST90018926", proxies=TRUE)
pqtl_dmSig <- harmonise_data (pqtl, pqtl_dmSig, action=2)
pqtl_dmSig <- pqtl_dmSig[pqtl_dmSig$mr_keep==TRUE,]

# Extract significantly dm-modifying variants
pqtl_dmSig_SNP <- mr_singlesnp(pqtl_dmSig, all_method = c("mr_ivw")) %>%
  filter(SNP != "All - Inverse variance weighted" & p < 0.05) %>%
  select(SNP)
pqtl_dmSig <- pqtl %>%
  filter(SNP %in% pqtl_dmSig_SNP$SNP); dim(pqtl_dmSig) #118, 13
pqtl_dmSig_clumped <- clump_data(pqtl_dmSig, clump_r2=0.1); nrow(pqtl_dmSig_clumped) #10
pqtl_dmSig_pd <- harmonise_data (pqtl_dmSig_clumped, pd, action=2)

# From pqtl_dm
pqtl_dmSig_SNP <- mr_singlesnp(pqtl_dm, all_method = c("mr_ivw")) %>%
  filter(SNP != "All - Inverse variance weighted" & p < 0.05) %>%
  select(SNP)
pqtl_dmSig <- pqtl_clumped %>%
  filter(SNP %in% pqtl_dmSig_SNP$SNP); dim(pqtl_dmSig) #7
pqtl_dmSig_pd <- harmonise_data (pqtl_dmSig, pd, action=2)

# MR
pqtl_dmSig_pd <- pqtl_dmSig_pd[pqtl_dmSig_pd$mr_keep==TRUE,]
pqtl_dmSig_pd$id.exposure <- "pQTLGen_dmSig"
pd_mr_list[['pQTL_dmSig']] <- mr(pqtl_dmSig_pd, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_dmSig_pd)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_dmSig_pd$id.exposure), id.outcome = unique(pqtl_dmSig_pd$id.outcome), outcome = unique(pqtl_dmSig_pd$outcome), exposure = unique(pqtl_dmSig_pd$exposure),
  method = method.presso, nsnp = nrow(pqtl_dmSig_pd), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pd_mr_list[['pQTL_dmSig']] <- rbind(pd_mr_list$pQTL_dmSig, df.presso)

# Male
#pqtl_dmSig_pdM <- harmonise_data (pqtl_dmSig_clumped, pdM, action=2)
pqtl_dmSig_pdM <- harmonise_data (pqtl_dmSig, pdM, action=2) #from pqtl_dm
pqtl_dmSig_pdM <- pqtl_dmSig_pdM[pqtl_dmSig_pdM$mr_keep==TRUE,]
pqtl_dmSig_pdM$id.exposure <- "cis-pQTL_dmSig"
pdM_mr_list[['pQTL_dmSig']] <- mr(pqtl_dmSig_pdM, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_dmSig_pdM)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_dmSig_pdM$id.exposure), id.outcome = unique(pqtl_dmSig_pdM$id.outcome), outcome = unique(pqtl_dmSig_pdM$outcome), exposure = unique(pqtl_dmSig_pdM$exposure),
  method = method.presso, nsnp = nrow(pqtl_dmSig_pdM), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdM_mr_list[['pQTL_dmSig']] <- rbind(pdM_mr_list$pQTL_dmSig, df.presso)

# Female
#pqtl_dmSig_pdF <- harmonise_data (pqtl_dmSig_clumped, pdF, action=2)
pqtl_dmSig_pdF <- harmonise_data (pqtl_dmSig, pdF, action=2) #from pqtl_dm
pqtl_dmSig_pdF <- pqtl_dmSig_pdF[pqtl_dmSig_pdF$mr_keep==TRUE,]
pqtl_dmSig_pdF$id.exposure <- "cis-pQTL_dmSig"
pdF_mr_list[['pQTL_dmSig']] <- mr(pqtl_dmSig_pdF, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_dmSig_pdF)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_dmSig_pdF$id.exposure), id.outcome = unique(pqtl_dmSig_pdF$id.outcome), outcome = unique(pqtl_dmSig_pdF$outcome), exposure = unique(pqtl_dmSig_pdF$exposure),
  method = method.presso, nsnp = nrow(pqtl_dmSig_pdF), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdF_mr_list[['pQTL_dmSig']] <- rbind(pdF_mr_list$pQTL_dmSig, df.presso)


## Insignificant SNPs with T2DM
# eQTL
# Extract Insignificantly dm-modifying variants
eqtl_dmInsig_SNP <- mr_singlesnp(eqtl_dm, all_method = c("mr_ivw")) %>%
  filter(SNP != "All - Inverse variance weighted" & p > 0.05) %>%
  select(SNP)
eqtl_dmInsig <- eqtl_clumped %>%
  filter(SNP %in% eqtl_dmInsig_SNP$SNP); dim(eqtl_dmInsig) #8
eqtl_dmInsig_pd <- harmonise_data (eqtl_dmInsig, pd, action=2)

# MR
eqtl_dmInsig_pd <- eqtl_dmInsig_pd[eqtl_dmInsig_pd$mr_keep==TRUE,]
eqtl_dmInsig_pd$id.exposure <- "eQTLGen_dmInsig"
pd_mr_list[['eQTL_dmInsig']] <- mr(eqtl_dmInsig_pd, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=eqtl_dmInsig_pd)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(eqtl_dmInsig_pd$id.exposure), id.outcome = unique(eqtl_dmInsig_pd$id.outcome), outcome = unique(eqtl_dmInsig_pd$outcome), exposure = unique(eqtl_dmInsig_pd$exposure),
  method = method.presso, nsnp = nrow(eqtl_dmInsig_pd), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pd_mr_list[['eQTL_dmInsig']] <- rbind(pd_mr_list$eQTL_dmInsig, df.presso)

# Male
eqtl_dmInsig_pdM <- harmonise_data (eqtl_dmInsig, pdM, action=2) #from eqtl_dm
eqtl_dmInsig_pdM <- eqtl_dmInsig_pdM[eqtl_dmInsig_pdM$mr_keep==TRUE,]
eqtl_dmInsig_pdM$id.exposure <- "eQTLGen_dmInsig"
pdM_mr_list[['eQTL_dmInsig']] <- mr(eqtl_dmInsig_pdM, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
# presso: not enough IVs
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=eqtl_dmInsig_pdM)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(eqtl_dmInsig_pdM$id.exposure), id.outcome = unique(eqtl_dmInsig_pdM$id.outcome), outcome = unique(eqtl_dmInsig_pdM$outcome), exposure = unique(eqtl_dmInsig_pdM$exposure),
  method = method.presso, nsnp = nrow(eqtl_dmInsig_pdM), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdM_mr_list[['eQTL_dmInsig']] <- rbind(pdM_mr_list$eQTL_dmInsig, df.presso)


# Female
eqtl_dmInsig_pdF <- harmonise_data (eqtl_dmInsig, pdF, action=2) #from eqtl_dm
eqtl_dmInsig_pdF <- eqtl_dmInsig_pdF[eqtl_dmInsig_pdF$mr_keep==TRUE,]
eqtl_dmInsig_pdF$id.exposure <- "eQTLGen_dmInsig"
pdF_mr_list[['eQTL_dmInsig']] <- mr(eqtl_dmInsig_pdF, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=eqtl_dmInsig_pdF)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(eqtl_dmInsig_pdF$id.exposure), id.outcome = unique(eqtl_dmInsig_pdF$id.outcome), outcome = unique(eqtl_dmInsig_pdF$outcome), exposure = unique(eqtl_dmInsig_pdF$exposure),
  method = method.presso, nsnp = nrow(eqtl_dmInsig_pdF), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdF_mr_list[['eQTL_dmInsig']] <- rbind(pdF_mr_list$eQTL_dmInsig, df.presso)


## pQTL ######
# Extract Insignificantly dm-modifying variants
pqtl_dmInsig_SNP <- mr_singlesnp(pqtl_dm, all_method = c("mr_ivw")) %>%
  filter(SNP != "All - Inverse variance weighted" & p > 0.05) %>%
  select(SNP)
pqtl_dmInsig <- pqtl_clumped %>%
  filter(SNP %in% pqtl_dmInsig_SNP$SNP); dim(pqtl_dmInsig) #34
pqtl_dmInsig_pd <- harmonise_data (pqtl_dmInsig, pd, action=2)

# MR
pqtl_dmInsig_pd <- pqtl_dmInsig_pd[pqtl_dmInsig_pd$mr_keep==TRUE,]
pqtl_dmInsig_pd$id.exposure <- "pQTLGen_dmInsig"
pd_mr_list[['pQTL_dmInsig']] <- mr(pqtl_dmInsig_pd, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_dmInsig_pd)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_dmInsig_pd$id.exposure), id.outcome = unique(pqtl_dmInsig_pd$id.outcome), outcome = unique(pqtl_dmInsig_pd$outcome), exposure = unique(pqtl_dmInsig_pd$exposure),
  method = method.presso, nsnp = nrow(pqtl_dmInsig_pd), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pd_mr_list[['pQTL_dmInsig']] <- rbind(pd_mr_list$pQTL_dmInsig, df.presso)

# Male
pqtl_dmInsig_pdM <- harmonise_data (pqtl_dmInsig, pdM, action=2) #from pqtl_dm
pqtl_dmInsig_pdM <- pqtl_dmInsig_pdM[pqtl_dmInsig_pdM$mr_keep==TRUE,]
pqtl_dmInsig_pdM$id.exposure <- "cis-pQTL_dmInsig"
pdM_mr_list[['pQTL_dmInsig']] <- mr(pqtl_dmInsig_pdM, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_dmInsig_pdM)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_dmInsig_pdM$id.exposure), id.outcome = unique(pqtl_dmInsig_pdM$id.outcome), outcome = unique(pqtl_dmInsig_pdM$outcome), exposure = unique(pqtl_dmInsig_pdM$exposure),
  method = method.presso, nsnp = nrow(pqtl_dmInsig_pdM), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdM_mr_list[['pQTL_dmInsig']] <- rbind(pdM_mr_list$pQTL_dmInsig, df.presso)

# Female
pqtl_dmInsig_pdF <- harmonise_data (pqtl_dmInsig, pdF, action=2) #from pqtl_dm
pqtl_dmInsig_pdF <- pqtl_dmInsig_pdF[pqtl_dmInsig_pdF$mr_keep==TRUE,]
pqtl_dmInsig_pdF$id.exposure <- "cis-pQTL_dmInsig"
pdF_mr_list[['pQTL_dmInsig']] <- mr(pqtl_dmInsig_pdF, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_dmInsig_pdF)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_dmInsig_pdF$id.exposure), id.outcome = unique(pqtl_dmInsig_pdF$id.outcome), outcome = unique(pqtl_dmInsig_pdF$outcome), exposure = unique(pqtl_dmInsig_pdF$exposure),
  method = method.presso, nsnp = nrow(pqtl_dmInsig_pdF), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
pdF_mr_list[['pQTL_dmInsig']] <- rbind(pdF_mr_list$pQTL_dmInsig, df.presso)


## 5-3: MR on REM-sleep Behavior Disorder
rbd_mr_list <- list(); rbd_hetero_list <- list()
rbd <- read.table("../../../gwas/GCST90204200.h.tsv",sep='\t',header=T)
rbd <- rbd[c("rsid","chromosome","base_pair_location","effect_allele","other_allele",
             "effect_allele_frequency","beta","standard_error","p_value")]
names(rbd) <- c("SNP","chr","pos","effect_allele.outcome","other_allele.outcome",
                "eaf.outcome","beta.outcome","se.outcome","pval.outcome")
rbd$outcome <- "REM-sleep Behavior Disorder"
rbd$id.outcome <- "GCST90204200"

eqtl_rbd <- harmonise_data (eqtl_clumped, rbd, action=2)
eqtl_rbd <- eqtl_rbd[eqtl_rbd$mr_keep==TRUE,]
rbd_mr_list[['eQTL']] <- mr(eqtl_rbd, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=eqtl_rbd)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(eqtl_rbd$id.exposure), id.outcome = unique(eqtl_rbd$id.outcome), outcome = unique(eqtl_rbd$outcome), exposure = unique(eqtl_rbd$exposure),
  method = method.presso, nsnp = nrow(eqtl_rbd), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
rbd_mr_list[['eQTL']] <- rbind(rbd_mr_list$eQTL, df.presso)
rbd_hetero_list[['eQTL']] <- mr_heterogeneity(eqtl_rbd)


pqtl_rbd <- harmonise_data (pqtl_clumped, rbd, action=3)
pqtl_rbd <- pqtl_rbd[pqtl_rbd$mr_keep==TRUE,]
rbd_mr_list[['pQTL']] <- mr(pqtl_rbd, method_list=c("mr_ivw","mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
presso <- mr_presso(BetaOutcome="beta.outcome",BetaExposure="beta.exposure",SdOutcome="se.outcome",
                    SdExposure="se.exposure",SignifThreshold=0.05,OUTLIERtest=TRUE,DISTORTIONtest=TRUE,data=pqtl_rbd)$`Main MR results`[1, ]
if (!all(is.na(presso[2, c("Causal Estimate", "Sd", "P-value")]))) {
  res.presso <- presso[2, ]
  method.presso <- "MR-PRESSO"
} else {
  res.presso <- presso[1, ]
  method.presso <- "MR-PRESSO"
}
df.presso <- data.frame(
  id.exposure = unique(pqtl_rbd$id.exposure), id.outcome = unique(pqtl_rbd$id.outcome), outcome = unique(pqtl_rbd$outcome), exposure = unique(pqtl_rbd$exposure),
  method = method.presso, nsnp = nrow(pqtl_rbd), b = res.presso$`Causal Estimate`, se = res.presso$Sd, pval = res.presso$`P-value`
)
rbd_mr_list[['pQTL']] <- rbind(rbd_mr_list$pQTL, df.presso)
rbd_hetero_list[['pQTL']] <- mr_heterogeneity(pqtl_rbd)

write.csv(do.call(rbind, rbd_mr_list), file="RBD_mr.csv", row.names = FALSE, na = "")
write.csv(do.call(rbind, rbd_hetero_list), file="RBD_hetero.csv", row.names = FALSE, na = "")



### 6: Sensitivity Analysis and Plot
## All
p1 <- mr_scatter_plot(pd_mr_list$eQTL, eqtl_pd)
ggsave(p1[[1]], file = "scatter_eqtl_pd.png", width=7, height=7)
res_single <- mr_singlesnp(eqtl_pd, all_method = c("mr_ivw"))
p2 <- mr_forest_plot(res_single)
ggsave(p2[[1]], file = "single_eqtl_pd.png", width=7, height=7)
res_loo <- mr_leaveoneout(eqtl_pd)
p3 <- mr_leaveoneout_plot(res_loo)
ggsave(p3[[1]], file = "loo_eqtl_pd.png", width=7, height=7)

p1 <- mr_scatter_plot(pd_mr_list$pQTL, pqtl_pd)
ggsave(p1[[1]], file = "scatter_pqtl_pd.png", width=7, height=7)
res_single <- mr_singlesnp(pqtl_pd, all_method = c("mr_ivw"))
p2 <- mr_forest_plot(res_single)
ggsave(p2[[1]], file = "single_pqtl_pd.png", width=7, height=7)
res_loo <- mr_leaveoneout(pqtl_pd)
p3 <- mr_leaveoneout_plot(res_loo)
ggsave(p3[[1]], file = "loo_pqtl_pd.png", width=7, height=7)


## Male
p1 <- mr_scatter_plot(pdM_mr_list$eQTL, eqtl_pdM)
ggsave(p1[[1]], file = "scatter_eqtl_pdM.png", width=7, height=7)
res_single <- mr_singlesnp(eqtl_pdM, all_method = c("mr_ivw"))
p2 <- mr_forest_plot(res_single)
ggsave(p2[[1]], file = "single_eqtl_pdM.png", width=7, height=7)
res_loo <- mr_leaveoneout(eqtl_pdM)
p3 <- mr_leaveoneout_plot(res_loo)
ggsave(p3[[1]], file = "loo_eqtl_pdM.png", width=7, height=7)

p1 <- mr_scatter_plot(pdM_mr_list$pQTL, pqtl_pdM)
ggsave(p1[[1]], file = "scatter_pqtl_pdM.png", width=7, height=7)
res_single <- mr_singlesnp(pqtl_pdM, all_method = c("mr_ivw"))
p2 <- mr_forest_plot(res_single)
ggsave(p2[[1]], file = "single_pqtl_pdM.png", width=7, height=7)
res_loo <- mr_leaveoneout(pqtl_pdM)
p3 <- mr_leaveoneout_plot(res_loo)
ggsave(p3[[1]], file = "loo_pqtl_pdM.png", width=7, height=7)


## Female
p1 <- mr_scatter_plot(pdF_mr_list$eQTL, eqtl_pdF)
ggsave(p1[[1]], file = "scatter_eqtl_pdF.png", width=7, height=7)
res_single <- mr_singlesnp(eqtl_pdF, all_method = c("mr_ivw"))
p2 <- mr_forest_plot(res_single)
ggsave(p2[[1]], file = "single_eqtl_pdF.png", width=7, height=7)
res_loo <- mr_leaveoneout(eqtl_pdF)
p3 <- mr_leaveoneout_plot(res_loo)
ggsave(p3[[1]], file = "loo_eqtl_pdF.png", width=7, height=7)

p1 <- mr_scatter_plot(pdF_mr_list$pQTL, pqtl_pdF)
ggsave(p1[[1]], file = "scatter_pqtl_pdF.png", width=7, height=7)
res_single <- mr_singlesnp(pqtl_pdF, all_method = c("mr_ivw"))
p2 <- mr_forest_plot(res_single)
ggsave(p2[[1]], file = "single_pqtl_pdF.png", width=7, height=7)
res_loo <- mr_leaveoneout(pqtl_pdF)
p3 <- mr_leaveoneout_plot(res_loo)
ggsave(p3[[1]], file = "loo_pqtl_pdF.png", width=7, height=7)


# Directionality test
steiger_list <- list()
steiger_list[['eQTL_pd']] <- directionality_test(eqtl_pd)
steiger_list[['eQTL_pdM']] <- directionality_test(eqtl_pdM)
steiger_list[['eQTL_pdF']] <- directionality_test(eqtl_pdF)

steiger_list[['pQTL_pd']] <- directionality_test(pqtl_pd)
steiger_list[['pQTL_pdM']] <- directionality_test(pqtl_pdM)
steiger_list[['pQTL_pdF']] <- directionality_test(pqtl_pdF)

write.csv(do.call(rbind, steiger_list), file="steiger.csv", row.names = FALSE, na = "")
