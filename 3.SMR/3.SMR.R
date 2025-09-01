
####### SMR analysis ----
system(paste(
  "./smr-1.3.1", 
  "--bfile", "./g1000_eur/g1000_eur", 
  "--gwas-summary", "./vit.txt", 
  "--beqtl-summary", "./eqtl/Whole_Blood.lite", 
  "--out", "./eQTLGen", 
  "--thread-num", "10"
))

system(paste(
  "./smr-1.3.1", 
  "--bfile", "./g1000_eur/g1000_eur", 
  "--gwas-summary", "./vit.txt", 
  "--beqtl-summary", "./gtex/Skin_Not_Sun_Exposed_Suprapubic.lite", 
  "--out", "./Not_Sun_Exposed", 
  "--thread-num", "10"
))

system(paste(
  "./smr-1.3.1", 
  "--bfile", "./g1000_eur/g1000_eur", 
  "--gwas-summary", "./vit.txt", 
  "--beqtl-summary", "./gtex/Skin_Sun_Exposed_Lower_leg.lite.besd", 
  "--out", "./Sun_Exposed", 
  "--thread-num", "10"
))

system(paste(
  "./smr-1.3.1", 
  "--bfile", "./g1000_eur/g1000_eur", 
  "--gwas-summary", "./vit.txt", 
  "--beqtl-summary", "./gtex/Cells_Cultured_fibroblasts.lite", 
  "--out", "./fibroblasts", 
  "--thread-num", "10"
))

system(paste(
  "./smr-1.3.1", 
  "--bfile", "./g1000_eur/g1000_eur", 
  "--gwas-summary", "./vit.txt", 
  "--beqtl-summary", "./mQTL/multiple", 
  "--out", "./mQTL", 
  "--thread-num", "10"
))