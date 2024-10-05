# check results by
# clu@emma ~/git/3D_cardiac_GWAS/data/48K/gwas_aha_t64/mr_2bp$ cat DBPieu-b-39/*res.txt |grep Inverse

########################################################################################
# This script is to run Mendelian Randomisation 
# 
# Exposure: SBP or DBP
#
# Outcomes: 16 AHA LV traits and 3 global LV traits
# 
# /Volumes/acasis/cardiac/48K/gwas_aha_t64/gwas_raw/40kahat64_sexagebmibsa2bp_*_chrmerged.regenie.gz
# /gwas_raw/AHA_40k_t14_2bp_full_", gtrait, "_chrmerged.regenie.gz"
########################################################################################


# setup
########################################################################################
library(TwoSampleMR)
library(R.utils)
library(ggplot2)
library(MRInstruments)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

exposure_name="DCM"

DCM_outcfilename = "/Volumes/acasis/cardiac/summary_stats/HERMES_HNDC_meta_analysis_2023/FORMAT-METAL_Pheno5_EUR.tsv"
HCM_outcfilename = "/Volumes/acasis/cardiac/summary_stats/HCM_meta_analysis_2023/gwama_sumstat/hcm.gwama.txt.gz"


root_dir = "/Volumes/acasis/cardiac/48K/gwas_aha_t64/"

working_dir = paste0(
    root_dir, "/mr_2bp/Exposure_", exposure_name, '_Outcome_LV-AHA'
)
dir.create(working_dir)
setwd(working_dir)

log_file = "test.log"


# Read keyname from arguments
# args <- commandArgs(trailingOnly = TRUE)
# keyname <- args[1]
########################################################################################


# Get instruments 
# exposure_all <- read_exposure_data(
#           filename = HCM_outcfilename,
#           clump = FALSE,
#           sep = "\t",
#         snp_col = "rsid",
#         beta_col = "beta",
#         se_col = "se",
#         effect_allele_col = "effect_allele",
#         other_allele_col = "noneffect_allele",
#         eaf_col = "eaf",
#         pval_col = "pvalue",
#         samplesize_col = "n_samples"
#       )
exposure_all <- read_exposure_data(
            filename = DCM_outcfilename,
            sep = "\t",
            snp_col = "rsID",
            beta_col = "A1_beta",
            se_col = "se",
            effect_allele_col = "A1",
            other_allele_col = "A2",
            eaf_col = "A1_freq",
            pval_col = "pval",
            samplesize_col = "N_total"
        )

exposure_all = subset(exposure_all, pval.exposure < 5e-8)

exposure_dat = clump_data(
        exposure_all,
        clump_kb = 10000,
        clump_r2 = 0.001,
        clump_p1 = 1,
        clump_p2 = 1,
        pop = "EUR",
        bfile = '/Users/clu/work/genomic_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC',
        plink_bin = '/Users/clu/softwoare/bin/plink'
        )

# select exposure data where pval.exposure < 5e-8
exposure_dat = subset(exposure_dat, pval.exposure < 5e-8)
exposure_dat$exposure = exposure_name

# function to run MR on AHA outcomes
run_MR_AHA_outcome <- function(keyname, i, exposure_dat, log_file){
    filename = paste0(
        root_dir, "/gwas_raw/40kahat64_sexagebmibsa2bp_", keyname, "_AHA_", i, "_chrmerged.regenie.gz"
    )

  # outcome data
  # CHROM	GENPOS	ID	ALLELE0	ALLELE1	A1FREQ	INFO	N	TEST	BETA	SE	CHISQ	LOG10P	EXTRA
  outcome_dat <- read_outcome_data(
      filename = filename,
      snps = exposure_dat$SNP,
      sep = "\t",
      snp_col = "ID",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      eaf_col = "A1FREQ",
      pval_col = "LOG10P",
      log_pval = TRUE,
      samplesize_col = "N",
      chr_col = "CHROM",
      pos_col = "GENPOS"
  )

  outcome_dat$outcome = paste0(keyname, "_AHA_", i)
        
  # Harmonise the exposure and outcome data and output log to file
  dat <- harmonise_data(exposure_dat, outcome_dat)

  # log size of outcome data
  write(paste0(paste0(keyname, "_AHA_", i), ", ", nrow(outcome_dat), ", ", nrow(dat)), log_file, append = TRUE)

  # Perform MR
  res <- mr(dat)

  # Save results
  outnametxt = paste0("40k_", keyname, "_AHA_", i, "_mr_", exposure_name, "_res.txt")
  write.table(res, outnametxt, sep = "\t", row.names = FALSE)
  # Save plot as png, high res
  outnamepng = paste0("40k_", keyname, "_AHA_", i, "_mr_", exposure_name, "_plot.png")
  g = mr_scatter_plot(res, dat)
  ggsave(g[[1]], filename = outnamepng, width = 6, height = 5, units = "in", dpi = 300)

  # done with this one
  print(paste0("Done: 40k_", keyname, "_AHA_", i, "_mr_", exposure_name, "_res.txt"))
}
# run it on all AHA outcomes
# 2bp
keyname = 'WT'
for (i in 1:16){
    tryCatch(
    {
      run_MR_AHA_outcome(keyname, i, exposure_dat, log_file)
    },
    error = function(e) {
        print(e)
    }
    )
}

keyname="Ecc"
for (i in 1:16){
    tryCatch(
    {
        run_MR_AHA_outcome(keyname, i, exposure_dat, log_file)
    },
    error = function(e) {
        print(e)
    }
    )
}

keyname="Err"
for (i in 1:16){
    tryCatch(
    {
        run_MR_AHA_outcome(keyname, i, exposure_dat, log_file)
    },
    error = function(e) {
        print(e)
    }
    )
}

# function to run MR on global outcomes
run_MR_global_outcome <- function(trait, filename, exposure_dat, log_file){

  # outcome data
  # CHROM	GENPOS	ID	ALLELE0	ALLELE1	A1FREQ	INFO	N	TEST	BETA	SE	CHISQ	LOG10P	EXTRA
  outcome_dat <- read_outcome_data(
      filename = filename,
      snps = exposure_dat$SNP,
      sep = "\t",
      snp_col = "ID",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      eaf_col = "A1FREQ",
      pval_col = "LOG10P",
      log_pval = TRUE,
      samplesize_col = "N",
      chr_col = "CHROM",
      pos_col = "GENPOS"
  )

  outcome_dat$outcome = trait
        
  # Harmonise the exposure and outcome data and output log to file
  dat <- harmonise_data(exposure_dat, outcome_dat)

  # log size of outcome data
  write(paste0(trait, ", ", nrow(outcome_dat), ", ", nrow(dat)), log_file, append = TRUE)

  # Perform MR
  res <- mr(dat)

  # Save results
  outnametxt = paste0("40k_", trait, "_mr_", exposure_name, "_res.txt")
  write.table(res, outnametxt, sep = "\t", row.names = FALSE)
  # Save plot as png, high res
  outnamepng = paste0("40k_", trait, "_mr_", exposure_name, "_plot.png")
  g = mr_scatter_plot(res, dat)
  ggsave(g[[1]], filename = outnamepng, width = 6, height = 5, units = "in", dpi = 300)

  # done with this one
  print(paste0("Done: 40k_", trait, "_mr_", exposure_name, "_res.txt"))
}

# run it on all global outcomes
for (gtrait in c("WT_Global", "Ecc_Global", "Err_Global")){
    tryCatch(
    {
        filename = paste0(
            root_dir, "/gwas_raw/AHA_40k_t14_2bp_full_", gtrait, "_chrmerged.regenie.gz"
        )
        run_MR_global_outcome(gtrait, filename, exposure_dat, log_file)
    },
    error = function(e) {
        print(e)
    }
    )
}

