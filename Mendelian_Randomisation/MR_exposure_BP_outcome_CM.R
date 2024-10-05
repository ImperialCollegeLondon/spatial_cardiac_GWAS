########################################################################################
# Run this as a sanity check to see if the MR results are consistent with the previous MR results
# SBP to HCM: #SNP: 446, IVW: beta 0.0297, p 1.8e-07; compared to Tardos2023: #SNP: 760, beta 0.027, p 6.9E-09 (SI Table 16)
# SBP to DCM: #SNP: 439, IVW: beta 0.0217, p 4.9E-11; Tardos2023 don't have this data.
# DBP to HCM: #SNP: 447, IVW: beta 0.046, p 3.1E-05; compared to Tardos2023: #SNP: 785, beta 0.041, p 2.1E-06 (SI Table 16) 
# DBP to DCM: #SNP: 440, IVW: beta 0.031, p 1.1E-07; Tardos2023 don't have this data.
########################################################################################


library(TwoSampleMR)
library(R.utils)
library(ggplot2)
library(MRInstruments)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

root_dir = "/Volumes/acasis/cardiac/summary_stats/MR/"

working_dir = paste0(
    root_dir, "/Exposure_BP_Outcome_CM"
)
dir.create(working_dir)
setwd(working_dir)

log_file = "test.log"

# Outcomes
DCM_outcfilename = "/Volumes/acasis/cardiac/summary_stats/HERMES_HNDC_meta_analysis_2023/FORMAT-METAL_Pheno5_EUR.tsv"
HCM_outcfilename = "/Volumes/acasis/cardiac/summary_stats/HCM_meta_analysis_2023/gwama_sumstat/hcm.gwama.txt.gz"


# Get instruments 
# exposure_gwas_key = "ieu-b-38"
# exposure_name="SBP"

exposure_gwas_key = "ieu-b-39"
exposure_name="DBP"

exposure_dat = extract_instruments(
    outcomes = exposure_gwas_key,
    p1 = 5e-8, # significance threshold for instruments, default 5e-8
    p2 = 1, # p-value threshold for clumping, default 5e-8
    r2 = 0.001, # clumping r2 cutoff, default 0.001
    clump = TRUE,
)
# log the size of instruments
write(paste0(exposure_gwas_key, " has number of instruments: ", length(exposure_dat$SNP)), log_file, append = TRUE)

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
outname = paste0("Exposure_", exposure_name, "_Outcome_HCM_res")
write.table(res, paste0(outname, ".txt"), sep = "\t", row.names = FALSE)
# Save plot as png, high res
g = mr_scatter_plot(res, dat)
ggsave(g[[1]], filename = paste0(outname, ".png"), width = 6, height = 5, units = "in", dpi = 300)
print(paste0("Done: ", outname))

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
outname = paste0("Exposure_", exposure_name, "_Outcome_DCM_res")
write.table(res, paste0(outname, ".txt"), sep = "\t", row.names = FALSE)
# Save plot as png, high res
g = mr_scatter_plot(res, dat)
ggsave(g[[1]], filename = paste0(outname, ".png"), width = 6, height = 5, units = "in", dpi = 300)
print(paste0("Done: ", outname))

