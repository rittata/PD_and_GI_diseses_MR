install.packages("dplyr")
install.packages("tibble")
install.packages("psych")
install.packages("plyr")
install.packages("markdown")
install.packages("googledrive")
install.packages("R.utils")
install.packages("ggplot2")
install.packages("gtable")
install.packages("utils")
install.packages("devtools")
install.packages("remotes")
devtools::install_github("explodecomputer/genetics.binaRies")
devtools::install_github("rondolab/MR-PRESSO")
remotes::install_github("MRCIEU/TwoSampleMR")

require(ggplot2)
library(dplyr)
library(data.table)
library(readr)
library(utils)
library(googledrive)
library(MRPRESSO)
library(TwoSampleMR)
drive_deauth()



### Loading exposure and outcome datasets. See dataset descriptions and download links in README.

# Outcome datasets (PD, B12)
drive_download(file = "https://drive.google.com/file/d/1dlc-FAHxe74nC-AnWWx0bfxjwhBPWhpz/view?usp=drive_link") # Vi.B12.fastGWA
drive_download(file = "https://drive.google.com/file/d/1NamM1FQjkeWuDTH-TrICo8yzVgI_J2rv/view?usp=drive_link") # META_no23_yesUKBB.txt
drive_download(file = "https://drive.google.com/file/d/1jSk0KsP_iWNJ-etq33R-VLjjhw1dkFwe/view?usp=drive_link") # IPDGC_sum_stats_no_UKB.txt
drive_download(file = "https://drive.google.com/file/d/1uSMfwM54N0dtEMCLjWxBc1XfKsR2eBwd/view?usp=drive_link") # IPDGC_AAO_GWAS_sumstats_april_2018_rsid.txt
PD_yesUKBB <- fread("META_no23_yesUKBB.txt")
PD_noUKBB <- fread("IPDGC_sum_stats_no_UKB.txt")
B12 <- fread("Vi.B12.fastGWA")
PD_AAO <- fread("IPDGC_AAO_GWAS_sumstats_april_2018_rsid.txt")
setnames(PD_AAO, c("CHR", "POS", "rsid", "A1", "A2", "Freq1", "FreqSE", "MinFreq", "MaxFreq",
                            "Effect","StdErr","P-value","Direction", "HetISq","HetChiSq","HetDf","HetPVal"))

# Outcome datasets (PD phenotypes and allele reference)
CI_base_not_merged_with_ref <- fread("/kaggle/input/project-summer-19-06/project_summer/base_DEMENTIA.txt/base_DEMENTIA.txt")
CI_surv_not_merged_with_ref <- fread("/kaggle/input/project-summer-19-06/project_summer/surv_DEMENTIA.txt/surv_DEMENTIA.txt")
HY_not_merged_with_ref <- fread("/kaggle/input/project-summer-19-06/project_summer/cont_HY.txt/cont_HY.txt")
MMSE_not_merged_with_ref <- fread("/kaggle/input/project-summer-19-06/project_summer/cont_MMSE.txt/cont_MMSE.txt")
UPDRS3_not_merged_with_ref <- fread("/kaggle/input/project-summer-19-06/project_summer/cont_UPDRS3_scaled.txt/cont_UPDRS3_scaled.txt")
reference <- fread("/kaggle/input/project-summer-19-06/project_summer/reference.txt") # refer
CI <- merge(CI_not_merged_with_ref, reference, by="SNP")
CI_surv <- merge(CI_surv_not_merged_with_ref, reference, by="SNP")
HY <- merge(PD_HY_not_merged_with_ref, reference, by="SNP")
MMSE <- merge(MMSE_not_merged_with_ref, reference, by="SNP")
UPDRS3 <- merge(UPDRS3_not_merged_with_ref, reference, by="SNP")

# Exposure datasets
DivD <- read_tsv("https://cnsgenomics.com/data/wu_et_al_2023_cg/DivD_EUR")
PUD <- read_tsv("https://cnsgenomics.com/data/wu_et_al_2021_nc/1_PUD_summary")
Panc_Alc_Chr <- fread("/kaggle/input/project-summer-19-06/project_summer/Panc_1_bin.tsv")
Panc_nonAlc_Chr <- fread("/kaggle/input/project-summer-19-06/project_summer/Panc_2_bin.tsv")
system("curl -o Panc_4_fss_bin.tsv.gz 'http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90255001-GCST90256000/GCST90255375/GCST90255375_buildGRCh37.tsv.gz'")
Panc_Ac <- fread("Panc_4_fss_bin.tsv.gz")



### Preparing exposure datasets for analysis: filtering by p-value threshold, setting the number of cases and controls
### according to articles and clumping with "relaxed" settings.

# Determining the path to the binary plink file and loading LD reference panels
plink_path <- genetics.binaRies::get_plink_binary()
drive_download(file = "https://drive.google.com/file/d/1r6mKJJicNbzQy_QRIBUror2gUu160YvO/view?usp=sharing", overwrite = TRUE)
dir.create("EUR")
unzip("/kaggle/working/EUR.zip", exdir = "EUR")
eur_path <- '/kaggle/working/EUR'

Panc_Alc_Chr_filtered <- Panc_Alc_Chr %>%
    filter(Panc_Alc_Chr$pval < 5e-8) %>%
    mutate(ncase = 1959,  ncontrol = 6040, samplesize = 7999)
Panc_Alc_Chr_filtered_formated <- format_data(data.frame(Panc_Alc_Chr_filtered), type = "exposure", min_pval = 1e-200)
Panc_Alc_Chr_clumped <- clump_data(Panc_Alc_Chr_filtered_formated, clump_kb = 10000, clump_r2 = 0.2, clump_p1 = 1, clump_p2 = 1, pop = "EUR",
                             bfile = eur_path, plink_bin = plink_path)

Panc_nonAlc_Chr_filtered <- Panc_nonAlc_Chr %>%
    filter(Panc_nonAlc_Chr$pval < 5e-8) %>%
    mutate(ncase = 584,  ncontrol = 6040, samplesize = 6624)
Panc_nonAlc_Chr_filtered_formated <- format_data(data.frame(Panc_nonAlc_Chr_filtered), type = "exposure", min_pval = 1e-200)
Panc_nonAlc_Chr_clumped <- clump_data(Panc_nonAlc_Chr_filtered_formated, clump_kb = 10000, clump_r2 = 0.2, clump_p1 = 1, clump_p2 = 1, pop = "EUR",
                             bfile = eur_path, plink_bin = plink_path)

Panc_Ac_filtered <- Panc_Ac %>%
    filter(Panc_Ac$p_value < 5e-8) %>%
    mutate(ncase = 10630,  ncontrol = 844679, samplesize = 855309)
Panc_Ac_filtered_formated <- format_data(data.frame(Panc_Ac_filtered), type = "exposure", snp_col = "variant_id",
                                  beta_col = "beta", se_col = "standard_error", eaf_col = "effect_allele_frequency",
                                  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                                  pval_col = "p_value", min_pval = 1e-200, chr_col = "chromosome",
                                  pos_col = "base_pair_location")
Panc_Ac_clumped <- clump_data(Panc_Ac_filtered_formated, clump_kb = 10000, clump_r2 = 0.2, clump_p1 = 1, clump_p2 = 1, pop = "EUR",
                                 bfile = eur_path, plink_bin = plink_path)

PUD_filtered <- PUD %>%
    filter(PUD$p < 5e-8) %>%
    mutate(ncase = 16666,  ncontrol = 439661)
PUD_filtered_formated <- format_data(data.frame(PUD_filtered), type = "exposure", snp_col = "SNP", beta_col = "b",
                                     se_col = "se", eaf_col = "freq", effect_allele_col = "A1", other_allele_col = "A2",
                                     pval_col = "p", samplesize_col = "N", min_pval = 1e-200)
PUD_clumped <- clump_data(PUD_filtered_formated, clump_kb = 10000, clump_r2 = 0.2, clump_p1 = 1, clump_p2 = 1, pop = "EUR",
                                bfile = eur_path, plink_bin = plink_path)

DivD_filtered <- DivD %>%
    filter(DivD$p < 5e-8) %>%
    mutate(ncase = 78399,  ncontrol = 645973)
DivD_filtered_formated <- format_data(data.frame(DivD_filtered), type = "exposure", snp_col = "SNP", beta_col = "b",
                                 se_col = "se", eaf_col = "freq", effect_allele_col = "A1",other_allele_col = "A2",
                                 pval_col = "p", samplesize_col = "N", min_pval = 1e-200)
DivD_clumped <- clump_data(DivD_filtered_formated, clump_kb = 10000, clump_r2 = 0.2, clump_p1 = 1, clump_p2 = 1, pop = "EUR",
                           bfile = eur_path, plink_bin = plink_path)



### Preparing datasets for analysis: renaming column names to the default names used in the TwoSampleMR package.

PD_yesUKBB_renamed_for_MR <- PD_yesUKBB %>%
    rename(effect_allele = A1, other_allele = A2, eaf = freq, beta = b, pval = p, ncase = N_cases, ncontrol = N_controls) %>%
    mutate(chr = sub("chr", "", sub(":.*", "", chrpos)), pos = as.numeric(sub(".*:", "", chrpos)))

PD_noUKBB_renamed_for_MR <- PD_noUKBB %>%
    rename(effect_allele = Allele1, other_allele = Allele2, eaf = Freq1, beta = Effect, pval = 'P-value', se = StdErr) %>%
    mutate(chr = sub("chr", "", sub(":.*", "", MarkerName)), pos = as.numeric(sub(".*:", "", MarkerName)))

PD_AAO_renamed_for_MR <- PD_AAO %>%
    rename(chr = CHR, pos = POS, SNP = rsid, effect_allele = A1, other_allele = A2, eaf = Freq1, beta = Effect, pval = 'P-value',
           se = StdErr) %>%
    mutate(ncase = 16502,  ncontrol = 17996, samplesize = 34498)

B12_renamed_for_MR <- B12 %>%
    rename(chr = chromosome, pos = base_pair_location, SNP = variant_id, effect_allele = A1, other_allele = A2, eaf = AF1,
           beta = BETA, pval = p_value, se = SE, samplesize = N)

CI_base_renamed_for_MR <- CI_base %>%
    rename(chrpos = SNP, SNP = RSID, chr = CHR, pos = START, effect_allele = ALT, other_allele = REF, eaf = MAF, beta = BETA,
           pval = P, se = SE, samplesize = N, gene = NearGENE)

CI_surv_renamed_for_MR <- CI_surv %>%
    rename(chrpos = SNP, SNP = RSID, chr = CHR, pos = START, effect_allele = ALT, other_allele = REF, eaf = MAF, beta = BETA,
           pval = P, se = SE, samplesize = N, gene = NearGENE)

HY_renamed_for_MR <- HY %>%
    rename(chrpos = SNP, SNP = RSID, chr = CHR, pos = START, effect_allele = ALT, other_allele = REF, eaf = MAF, beta = BETA,
           pval = P, se = SE, samplesize = N, gene = NearGENE)

MMSE_renamed_for_MR <- MMSE %>%
    rename(chrpos = SNP, SNP = RSID, chr = CHR, pos = START, effect_allele = ALT, other_allele = REF, eaf = MAF, beta = BETA,
           pval = P, se = SE, samplesize = N, gene = NearGENE)

UPDRS3_renamed_for_MR <- UPDRS3 %>%
    rename(chrpos = SNP, SNP = RSID, chr = CHR, pos = START, effect_allele = ALT, other_allele = REF, eaf = MAF, beta = BETA,
           pval = P, se = SE, samplesize = N, gene = NearGENE)



### MR analysis using TwoSampleMR and MR-PRESSO packages.

dir.create("MR_results")
setwd("/kaggle/working/MR_results")
exp <-list(Panc_Alc_Chr_clumped, Panc_nonAlc_Chr_clumped, Panc_Ac_clumped, DivD_clumped, PUD_clumped)
out <- list(PD_yesUKBB_renamed_for_MR, PD_noUKBB_renamed_for_MR, CI_base_renamed_for_MR, CI_surv_renamed_for_MR, HY_renamed_for_MR,
            MMSE_renamed_for_MR, UPDRS3_renamed_for_MR, PD_AAO_renamed_for_MR, B12_renamed_for_MR)
exp_names <- c('Panc_Alc_Chr', 'Panc_nonAlc_Chr', 'Panc_Ac', 'DivD', 'PUD')
out_names <- c('PD_yesUKBB', 'PD_noUKBB', 'CI_base', 'CI_surv', 'HY', 'MMSE', 'UPDRS3', 'PD_AAO', 'B12')

# Setting the disease prevalence to calculate R-squared. See README for links
prevalance_exp <- c(0.0001993, 0.0001772, 0.000348, 0.000012, 0.0014)
prevalance_out <- c(0.02802, 0.02802, 0.12, 0.12, NA, NA, NA, NA, NA)


for(i in 1:length(exp)){
   for(j in 1:length(out)){

      folder <- paste(exp_names[i], " VS ", out_names[j], sep="")
      unlink(paste("/kaggle/working/", folder, sep=""), recursive = TRUE)
      unlink(paste("/kaggle/working/", folder, ".zip", sep=""), recursive = TRUE)
      dir.create(folder)

      exp_data <- exp[[i]]
      out_data <- format_data(data.frame(out[[j]]), type = "outcome", snps = exp_data$SNP)
      harm_data <- harmonise_data(exposure_dat=exp_data, outcome_dat=out_data, action=2)
      harm_data <-subset(harm_data, harm_data$eaf.exposure!="NA")

      # R-squared calculation for exposure and outcome. get_r_from_lor() was used for binary traits,
      # get_r_from_pn() for continuous traits.
      harm_data$prevalence.exposure <- prevalance_exp[i]
      harm_data$r.exposure <- get_r_from_lor(
          lor = harm_data$beta.exposure,
          af = harm_data$eaf.exposure,
          ncase = harm_data$ncase.exposure,
          ncontrol = harm_data$ncontrol.exposure,
          prevalence = harm_data$prevalence.exposure,
          model = "logit",
          correction = FALSE)
      harm_data$units.exposure <-"log odds"
      harm_data$rsq.exposure <- harm_data$r.exposure ** 2
      if (j %in% c(1,2,3,4)){
          harm_data$prevalence.outcome <- prevalance_out[j]
          harm_data$r.outcome <- get_r_from_lor(
              lor = harm_data$beta.outcome,
              af = harm_data$eaf.outcome,
              ncase = harm_data$ncase.outcome,
              ncontrol = harm_data$ncontrol.outcome,
              prevalence = harm_data$prevalence.outcome,
              model = "logit",
              correction = FALSE)
          harm_data$units.outcome <-"log odds"}
      else {harm_data$r.outcome <- get_r_from_pn(harm_data$pval.outcome, harm_data$samplesize.outcome)
           harm_data$units.outcome <-"SD units"}
      harm_data$rsq.outcome <- harm_data$r.outcome ** 2

      # Using the MR-PRESSO method for harmonised datasets with more than 2 SNPs
      if(nrow(harm_data) > length(harm_data$beta.exposure) + 2) {
          presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
                              SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = harm_data,
                              NbDistribution = 1000,  SignifThreshold = 0.05)
          capture.output(print(presso), file = paste(folder, "presso.txt", sep = "/"))}

      # Saving MR results (with calculated F-statistics) in a table
      mr_res <- mr(harm_data)
      mr_res_with_OR <- generate_odds_ratios(mr_res)
      mr_res_with_OR$R2_calc <- mean(harm_data$rsq.exposure)
      n <- mean(harm_data$samplesize.exposure)
      k <- nrow(subset(harm_data, harm_data$ambiguous == FALSE))
      R2 <- mean(harm_data$rsq.exposure)
      mr_res_with_OR$F <- (R2*(n-1-k))/((1-R2)*k)
      write.table(mr_res_with_OR, file=file.path(folder, "mr_res_with_OR.txt"), row.names=FALSE, quote=FALSE, sep = "\t")

      # Create a report with graphs and save sensitivity analysis results in a table for harmonised datasets with more than 2 SNPs
      if(nrow(harm_data) >= 2){
          pleiotropy <- mr_pleiotropy_test(harm_data)
          pleiotropy$method <- 'MR Egger'
          heterogeneity <- mr_heterogeneity(harm_data)
          steiger_direction <- directionality_test(harm_data)
          write.table(pleiotropy, file=file.path(folder, "pleiotropy.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
          write.table(heterogeneity, file=file.path(folder, "heterogeneity.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
          write.table(steiger_direction, file=file.path(folder, "steiger_direction.txt"), row.names=FALSE, quote=FALSE, sep = "\t")
          mr_report(harm_data, study = folder,output_path = folder)
          res_single <- mr_singlesnp(harm_data)
          p5 <- mr_forest_plot(res_single)
          p5[[1]] }

      zip(zipfile = paste(folder, ".zip", sep=""), files = folder)
      gc()
  }
}



### Collecting all tabulated MR results into one table ("MR_results_table.tsv").

all_data <- list()
folder_names <- list.dirs(full.names = FALSE, recursive = FALSE)

for(i in 1:length(folder_names)){
      folder <- folder_names[[i]]
      mr_res_with_OR <- read.table(paste0(folder, "/mr_res_with_OR.txt"), sep="\t", header = TRUE)
      mr_res_with_OR$MR_name <- folder
      if (file.exists(paste0(folder, "/pleiotropy.txt"))) {
          pleiotropy <- read.table(paste0(folder, "/pleiotropy.txt"), header = TRUE, sep = "\t")
          heterogeneity <- read.table(paste0(folder, "/heterogeneity.txt"), header = TRUE, sep = "\t")
          steiger_direction <- read.table(paste0(folder, "/steiger_direction.txt"), header = TRUE, sep = "\t")
          mr_res_with_OR <- mr_res_with_OR %>%
               full_join(pleiotropy, by = c("id.exposure", "id.outcome", "method")) %>%
               full_join(heterogeneity, by = c("id.exposure", "id.outcome", "method")) %>%
               full_join(steiger_direction, by = c("id.exposure", "id.outcome"))
          all_data_1[[i]] <- mr_res_with_OR}
      else {all_data_2[[i]] <- mr_res_with_OR}}

MR_results_table <- do.call(rbind, all_data)
write.table(MR_results_table, "MR_results_table.tsv", sep="\t", row.names = FALSE)

zip(zipfile = "MR_results.zip", '/kaggle/working/MR_results')
