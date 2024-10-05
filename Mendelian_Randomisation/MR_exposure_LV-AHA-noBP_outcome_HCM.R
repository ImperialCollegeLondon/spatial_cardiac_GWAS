# check results by
# clu@emma ~/git/3D_cardiac_GWAS/data/48K/gwas_aha_t64/mr_2bp$ cat DBPieu-b-39/*res.txt |grep Inverse

########################################################################################
# This script is to run Mendelian Randomisation 
# 
# Exposure: 16 AHA LV traits and 3 global LV traits
# 
# Outcomes: HCM
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

root_dir = "/Volumes/acasis/cardiac/48K/gwas_aha_t64/"

aha_gwas_dir = paste0(root_dir, "gwas_raw_nobp/")

working_dir = paste0(
    root_dir, "/mr_nobp/Exposure_LV_Outcome_CM"
)
dir.create(working_dir)
setwd(working_dir)

log_file = "test.log"

# Read keyname from arguments
# args <- commandArgs(trailingOnly = TRUE)
# keyname <- args[1]
########################################################################################


DCM_outcfilename = "/Volumes/acasis/cardiac/summary_stats/HERMES_HNDC_meta_analysis_2023/FORMAT-METAL_Pheno5_EUR.tsv"
HCM_outcfilename = "/Volumes/acasis/cardiac/summary_stats/HCM_meta_analysis_2023/gwama_sumstat/hcm.gwama.txt.gz"

# Get instruments 

run_MR_AHA_exposure <- function(filename, keyname, i, log_file){
  # catch exception
    tryCatch(
    {
      # Get instruments

      # CHROM	GENPOS	ID	ALLELE0	ALLELE1	A1FREQ	INFO	N	TEST	BETA	SE	CHISQ	LOG10P	EXTRA
      exposure_all <- read_exposure_data(
          filename = filename,
          clump = FALSE,
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
      exposure_dat$exposure = paste0(keyname, "_AHA_", i)

      # Get effects of instruments on outcome HCM
      HCM_outcome_dat <- read_outcome_data(
        snps = exposure_dat$SNP,
        filename = HCM_outcfilename,
        sep = "\t",
        snp_col = "rsid",
        beta_col = "beta",
        se_col = "se",
        effect_allele_col = "effect_allele",
        other_allele_col = "noneffect_allele",
        eaf_col = "eaf",
        pval_col = "pvalue",
        samplesize_col = "n_samples"
      )
      HCM_outcome_dat$outcome = "HCM"
      
      # Harmonise the exposure and outcome data
      dat <- harmonise_data(exposure_dat, HCM_outcome_dat)

      # Perform MR
      res <- mr(dat)

      # Save results
      outnametxt = paste0("40k_", keyname, "_AHA_", i, "_mr_HCM_res.txt")
      write.table(res, outnametxt, sep = "\t", row.names = FALSE)
      # Save plot as png, high res
      outnamepng = paste0("40k_", keyname, "_AHA_", i, "_mr_HCM_plot.png")
      g = mr_scatter_plot(res, dat)
      ggsave(g[[1]], filename = outnamepng, width = 6, height = 5, units = "in", dpi = 300)
      print(paste0("Done: 40k_", keyname, "_AHA_", i, "_mr_HCM_res.txt"))

      # Get effects of instruments on outcome DCM
        DCM_outcome_dat <- read_outcome_data(
            snps = exposure_dat$SNP,
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
        DCM_outcome_dat$outcome = "DCM"

        # Harmonise the exposure and outcome data
        dat <- harmonise_data(exposure_dat, DCM_outcome_dat)

        # Perform MR
        res <- mr(dat)

        # Save results
        outnametxt = paste0("40k_", keyname, "_AHA_", i, "_mr_DCM_res.txt")
        write.table(res, outnametxt, sep = "\t", row.names = FALSE)
        # Save plot as png, high res
        outnamepng = paste0("40k_", keyname, "_AHA_", i, "_mr_DCM_plot.png")
        g = mr_scatter_plot(res, dat)
        ggsave(g[[1]], filename = outnamepng, width = 6, height = 5, units = "in", dpi = 300)
        print(paste0("Done: 40k_", keyname, "_AHA_", i, "_mr_DCM_res.txt"))

      return (exposure_dat)
    },
    error = function(e) {
        print(e)
    }
    )
}

keyname = 'WT'
for (i in 1:16){
    filename = paste0(aha_gwas_dir, "/40kahat64_sexagebmibsa_", keyname, "_AHA_", i, "_chrmerged.regenie.gz")
    run_MR_AHA_exposure(filename, keyname, i, log_file)
}

keyname="Ecc"
for (i in 1:16){
    filename = paste0(aha_gwas_dir, "/40kahat64_sexagebmibsa_", keyname, "_AHA_", i, "_chrmerged.regenie.gz")
    run_MR_AHA_exposure(filename, keyname, i, log_file)
}

keyname="Err"
for (i in 1:16){
    filename = paste0(aha_gwas_dir, "/40kahat64_sexagebmibsa_", keyname, "_AHA_", i, "_chrmerged.regenie.gz")
    run_MR_AHA_exposure(filename, keyname, i, log_file)
}

# global outcomes

for (gtrait in c("WT_Global", "Ecc_Global", "Err_Global")){
    tryCatch(
    {
        filename = paste0(
            aha_gwas_dir, "/AHA_40k_t14_nobp_full_", gtrait, "_chrmerged.regenie.gz"
        )
        run_MR_AHA_exposure(filename, gtrait, i, log_file)
    },
    error = function(e) {
        print(e)
    }
    )
}

