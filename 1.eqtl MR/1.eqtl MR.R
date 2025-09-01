### 1. Obtain gene eQTL data ----
# Gene list and related information
{
library(dplyr)
library(TwoSampleMR)
library(readxl)
library(data.table)
library(ieugwasr)
library(parallel)
library(readxl)
library(CMplot)
library(openxlsx)
library(ggsci)
library(ggplot2)
library(ggiraph)
library(ggrepel)
}
# config.R
# ... 所有其他路径定义 ...
# Load drug genelist file
load("./drug_genelist.rda") 
# Load eaf file (allele frequency file for imputation)
load("./eQTLGen_eaf.rda")
# Load eQTL data (from eqgen website, FDR < 0.05)
load("./eQTLGen_eqtl.rda")
# Merge eQTL data and eaf file
eqtl_total = merge(eqtl, eaf, by.x = 'SNP', by.y = 'SNP')

# Align alleles and adjust eaf
eqtl_total <- eqtl_total %>%
  mutate(
    AlleleA_new = ifelse(AssessedAllele == AlleleA & OtherAllele == AlleleB, OtherAllele, AlleleA),
    AlleleB_new = ifelse(AssessedAllele == AlleleA & OtherAllele == AlleleB, AssessedAllele, AlleleB),
    AlleleB_all_new = ifelse(AssessedAllele == AlleleA & OtherAllele == AlleleB, 1 - AlleleB_all, AlleleB_all)
  )

eqtl_total$AlleleA = eqtl_total$AlleleA_new
eqtl_total$AlleleB = eqtl_total$AlleleB_new
eqtl_total$AlleleB_all = eqtl_total$AlleleB_all_new

# Rename gene symbol column
colnames(drug)[3] = 'GeneSymbol'

# Extract SNPs located within ±5000 bp of genes in the list
result <- eqtl_total %>%
  inner_join(drug, by = 'GeneSymbol') %>%
  filter(SNPPos >= start_b37 - 5000 & SNPPos <= end_b37 + 5000)

# Save filtered result
write.csv(result, 'filtered_eqtls.csv', row.names = FALSE)

# Summary statistics
table(result$GeneSymbol)
drug_filter = unique_genes <- unique(result$GeneSymbol)
drug_filter = data.frame(GeneSymbol = drug_filter)

# Calculate beta and standard error (SE)
result$beta <- result$Zscore / sqrt(2 * result$AlleleB_all * (1 - result$AlleleB_all) * (result$NrSamples + result$Zscore^2)) 
result$se   <- 1 / sqrt(2 * result$AlleleB_all * (1 - result$AlleleB_all) * (result$NrSamples + result$Zscore^2))

# Calculate F statistic
result_filter = result[, c(1,2,5,6,9,14,15,23,38,39)]
result_filter$F <- (result_filter$beta / result_filter$se)^2
write.csv(result_filter, 'result_filterF.csv', row.names = FALSE)

# Write SNPs to separate files by gene
save_snps_for_genes <- function(drug_filter, result_filter) {
  # Loop through each row in drug_filter
  for (i in 1:nrow(drug_filter)) {
    # Get current gene symbol
    current_gene <- drug_filter$GeneSymbol[i]
    
    # Filter SNPs belonging to the current gene
    snps_subset <- result_filter %>%
      filter(GeneSymbol == current_gene)
    
    # If non-empty, save as a CSV file
    if (nrow(snps_subset) > 0) {
      # Create output directory
      dir.create(file.path("./drugnew", current_gene), recursive = TRUE, showWarnings = FALSE)
      
      # Build file name
      file_name <- file.path("./drugnew", current_gene, paste0(current_gene, "_snps.csv"))
      
      # Write to CSV
      write.csv(snps_subset, file_name, row.names = FALSE)
    }
  }
}


### 2. Discovery set Mendelian Randomization ----

# Get all subfolder paths in the "drug" directory (each contains SNPs for one gene, >2000 gene folders)
drug_path <- "./drugnew"
folder_paths <- list.dirs(drug_path, full.names = TRUE, recursive = FALSE)

# Load outcome data
outcome <- fread("./vit.csv")

# Set output directory for non-significant results
no_significant_dir <- "./MR_result/discovery/no_significant"
if (!dir.exists(no_significant_dir)) {
  dir.create(no_significant_dir, recursive = TRUE)
}

# Define function to process a single gene folder
process_gene_folder <- function(folder) {
  gene_name <- basename(folder)
  cat("Processing gene:", gene_name, "\n")  # Display gene name in console
  file_path <- file.path(folder, paste0(gene_name, "_snps.csv"))
  
  # Load SNP data for the gene
  tryCatch({
    gene_data <- fread(file_path)
    
    # Format exposure data
    formatted_exposure_data <- format_data(
      gene_data,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      pval_col = "Pvalue",
      effect_allele_col = "AssessedAllele",
      other_allele_col = "OtherAllele"
    )
    
    # Keep only clumped SNPs
    formatted_exposure_data <- clump_data(formatted_exposure_data,
                                          clump_kb = 10000,
                                          clump_r2 = 0.2,
                                          clump_p1 = 0.99)
    
    # Merge exposure and outcome data
    exposure_outcome <- merge(formatted_exposure_data, outcome, by.x = "SNP", by.y = "SNP")
    
    # If no SNPs left after merging
    if (nrow(exposure_outcome) == 0) {
      gene_output_dir <- file.path("./MR_result/discovery/no_snps", gene_name)
      if (!dir.exists(gene_output_dir)) {
        dir.create(gene_output_dir, recursive = TRUE)
      }
      write.csv(data.frame(Message = "No SNPs left after clumping"), 
                file.path(gene_output_dir, paste0(gene_name, "_no_snps.csv")), row.names = FALSE)
      return(NULL)
    }
    
    # Format outcome data
    formatted_outcome_data <- format_data(
      exposure_outcome,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "SE",
      pval_col = "P",
      effect_allele_col = "A1",
      other_allele_col = "A2"
    )
    
    # Harmonise exposure and outcome data
    harmonised_data <- harmonise_data(formatted_exposure_data, formatted_outcome_data)
    
    # Calculate F statistic and filter weak instruments
    harmonised_data$F <- (harmonised_data$beta.exposure / harmonised_data$se.exposure)^2
    harmonised_data <- subset(harmonised_data, F > 10)
    
    # If no valid SNPs remain after filtering
    if (nrow(harmonised_data) == 0) {
      gene_output_dir <- file.path("./MR_result/discovery/no_valid_snps", gene_name)
      if (!dir.exists(gene_output_dir)) {
        dir.create(gene_output_dir, recursive = TRUE)
      }
      write.csv(data.frame(Message = "No valid SNPs after filtering"), 
                file.path(gene_output_dir, paste0(gene_name, "_no_valid_snps.csv")), row.names = FALSE)
      return(NULL)
    }
    
    # Perform MR depending on number of SNPs
    if (nrow(harmonised_data) == 1) {
      mr_results <- mr(harmonised_data)
    } else if (nrow(harmonised_data) > 1) {
      heterogeneity_results <- mr_heterogeneity(harmonised_data)
      het_pval = heterogeneity_results[heterogeneity_results$method == "Inverse variance weighted", "Q_pval"]
      if (het_pval < 0.05) {
        mr_results <- mr(harmonised_data, method_list = c("mr_ivw_mre"))
      } else {
        mr_results <- mr(harmonised_data, method_list = c("mr_ivw_fe"))
      }
    }
    
    # Convert to odds ratios
    OR <- generate_odds_ratios(mr_results)
    significant_results <- subset(OR, OR$pval < 0.05 / 2573)  # Bonferroni correction
    
    
    if (nrow(significant_results) > 0) {
      # Output directory for significant results
      gene_output_dir <- file.path("./MR_result/discovery/vit_eqtl_result", gene_name)
      if (!dir.exists(gene_output_dir)) {
        dir.create(gene_output_dir, recursive = TRUE)
      }
      write.csv(OR, file.path(gene_output_dir, paste0(gene_name, "_MR_results.csv")), row.names = FALSE)
      
      # Only output heterogeneity and pleiotropy tests if significant results exist
      if (nrow(significant_results) > 0) {
        if (exists("heterogeneity_results")) {
          write.csv(heterogeneity_results, file.path(gene_output_dir, paste0(gene_name, "_heterogeneity_results.csv")), row.names = FALSE)
        }
        
        pleio_result <- mr_pleiotropy_test(harmonised_data)
        if (exists("pleio_result")) {
          pleio_filename_suffix <- if (pleio_result$pval < 0.05) {"_pleio_results_fail.csv"} else {"_pleio_results.csv"}
          pleio_filepath <- file.path(gene_output_dir, paste0(gene_name, pleio_filename_suffix))
          write.csv(pleio_result, pleio_filepath, row.names = FALSE)
        }
      }
    } else {
      # Save to non-significant results folder
      write.csv(OR, file.path(no_significant_dir, paste0(gene_name, "_MR_results.csv")), row.names = FALSE)
    }
    
  }, error = function(e) {
    message("Error processing ", gene_name, ": ", e$message)
  })
}

# Parallel processing of each gene folder
mclapply(folder_paths, process_gene_folder, mc.cores = 30)

### 3. Validation set Mendelian Randomization analysis ----

# Get all subfolder paths under the "drug" directory
drug_path <- "./drugnew"
folder_paths <- list.dirs(drug_path, full.names = TRUE, recursive = FALSE)

# Select specific candidate genes
selected_genes <- c("CASP7", "CASP8", "CCDC3", "CD81", "CD226","CNP","CSK","CTSS","EHMT2","ERBB3",
                    "FEN1","FYN","HRAS","MPI","PPP3CA","PTPRC","SLC22A1","SRM","TBXAS1","ULK4")
selected_folder_paths <- folder_paths[basename(folder_paths) %in% selected_genes]


# Load outcome dataset (FinnGen vitiligo GWAS)
outcome <- fread("./finngen_R10_L12_VITILIGO.gz")

# Set directory for non-significant results
no_significant_dir <- "./MR_result/replication/no_significant"
if (!dir.exists(no_significant_dir)) {
  dir.create(no_significant_dir, recursive = TRUE)
}

# Define function to process a single gene folder
process_gene_folder <- function(folder) {
  gene_name <- basename(folder)
  cat("Processing gene:", gene_name, "\n")  # Display gene name in console
  file_path <- file.path(folder, paste0(gene_name, "_snps.csv"))
  
  # Load SNP data
  tryCatch({
    gene_data <- fread(file_path)
    
    # Format exposure data
    formatted_exposure_data <- format_data(
      gene_data,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      pval_col = "Pvalue",
      effect_allele_col = "AssessedAllele",
      other_allele_col = "OtherAllele"
    )
    
    # Perform LD clumping
    formatted_exposure_data <- clump_data(formatted_exposure_data,
                                          clump_kb = 10000,
                                          clump_r2 = 0.2,
                                          clump_p1 = 0.99)

    
    # Merge exposure and outcome data
    exposure_outcome <- merge(formatted_exposure_data, outcome, by.x = "SNP", by.y = "rsids")
    
    # If no SNPs remain
    if (nrow(exposure_outcome) == 0) {
      gene_output_dir <- file.path("./MR_result/replication/no_snps", gene_name)
      if (!dir.exists(gene_output_dir)) {
        dir.create(gene_output_dir, recursive = TRUE)
      }
      write.csv(data.frame(Message = "No SNPs left after clumping"), 
                file.path(gene_output_dir, paste0(gene_name, "_no_snps.csv")), row.names = FALSE)
      return(NULL)
    }
    
    # Format outcome data
    formatted_outcome_data <- format_data(
      exposure_outcome,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "sebeta",
      pval_col = "pval",
      effect_allele_col = "alt",
      other_allele_col = "ref"
    )
    
    # Harmonise exposure and outcome data
    harmonised_data <- harmonise_data(formatted_exposure_data, formatted_outcome_data)
    
    # Calculate F statistic and filter weak instruments
    harmonised_data$F <- (harmonised_data$beta.exposure / harmonised_data$se.exposure)^2
    harmonised_data <- subset(harmonised_data, F > 10)
    
    # If no valid SNPs remain
    if (nrow(harmonised_data) == 0) {
      gene_output_dir <- file.path("./MR_result/replication/no_valid_snps", gene_name)
      if (!dir.exists(gene_output_dir)) {
        dir.create(gene_output_dir, recursive = TRUE)
      }
      write.csv(data.frame(Message = "No valid SNPs after filtering"), 
                file.path(gene_output_dir, paste0(gene_name, "_no_valid_snps.csv")), row.names = FALSE)
      return(NULL)
    }
    
    # Perform MR analysis depending on number of SNPs
    if (nrow(harmonised_data) == 1) {
      mr_results <- mr(harmonised_data)
    } else if (nrow(harmonised_data) > 1) {
      heterogeneity_results <- mr_heterogeneity(harmonised_data)
      het_pval = heterogeneity_results[heterogeneity_results$method == "Inverse variance weighted", "Q_pval"]
      if (het_pval < 0.05) {
        mr_results <- mr(harmonised_data, method_list = c("mr_ivw_mre"))
      } else {
        mr_results <- mr(harmonised_data, method_list = c("mr_ivw_fe"))
      }
    }
    
    # Convert to odds ratios
    OR <- generate_odds_ratios(mr_results)
    significant_results <- subset(OR, OR$pval < 0.05)
    
    if (nrow(significant_results) > 0) {
      # Output directory for significant results
      gene_output_dir <- file.path("./MR_result/replication/vit_eqtl_result", gene_name)
      if (!dir.exists(gene_output_dir)) {
        dir.create(gene_output_dir, recursive = TRUE)
      }
      write.csv(OR, file.path(gene_output_dir, paste0(gene_name, "_MR_results.csv")), row.names = FALSE)
      
      # Only output heterogeneity and pleiotropy tests if significant
      if (exists("heterogeneity_results")) {
        write.csv(heterogeneity_results, file.path(gene_output_dir, paste0(gene_name, "_heterogeneity_results.csv")), row.names = FALSE)
      }
      
      pleio_result <- mr_pleiotropy_test(harmonised_data)
      if (exists("pleio_result")) {
        pleio_filename_suffix <- if (pleio_result$pval < 0.05) {"_pleio_results_fail.csv"} else {"_pleio_results.csv"}
        pleio_filepath <- file.path(gene_output_dir, paste0(gene_name, pleio_filename_suffix))
        write.csv(pleio_result, pleio_filepath, row.names = FALSE)
      }
    } else {
      # Save to non-significant folder
      write.csv(OR, file.path(no_significant_dir, paste0(gene_name, "_MR_results.csv")), row.names = FALSE)
    }
    
  }, error = function(e) {
    message("Error processing ", gene_name, ": ", e$message)
  })
}

# Run MR analysis for selected genes in parallel
mclapply(selected_folder_paths, process_gene_folder, mc.cores = 20)


####4.discovery MR 曼哈顿图
all_data <- data.frame()
file_paths <- list.files(
  path = "./MR_result/discovery/vit_eqtl_result", 
  pattern = "MR_results\\.csv$", 
  full.names = TRUE
)

for (file_path in file_paths) {
  df <- read.csv(file_path)
  full_file_name <- tools::file_path_sans_ext(basename(file_path))
  file_name <- sub("_MR.*", "", full_file_name)
  df$Gene <- file_name
  all_data <- rbind(all_data, df)
}

all_data1 <- data.frame()
file_paths <- list.files(
  path = "./MR_result/discovery/no_significant",
  pattern = "_MR_results\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

# 使用循环依次读取每个文件，并添加处理后的文件名作为新列
for (file_path in file_paths) {
  # 读取文件
  df <- read.csv(file_path)
  
  # 获取文件名（不包含扩展名）
  full_file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # 提取_MR前的字符作为新列
  file_name <- sub("_MR.*", "", full_file_name)
  
  # 添加文件名为新列
  df$Gene <- file_name
  
  # 将当前数据框与所有数据合并
  all_data1 <- rbind(all_data1, df)
}

head(all_data1)
all=rbind(all_data,all_data1)

all <- all %>% select(Gene, everything())

data <- read_excel("./NIHMS80906-supplement-Table_S1.xlsx", sheet = 1)
all_updated <- all %>%
  left_join(data %>% select(hgnc_names, chr_b37), by = c("Gene" = "hgnc_names"))
all_updated <- all_updated %>%
  left_join(data %>% select(hgnc_names, start_b37), by = c("Gene" = "hgnc_names"))

all=all_updated
save(all,file="./Manhattan_plot_data.rda")
cutoff_pvalue <- 0.05/2573

dat_plot <- dat_plot %>%
  mutate(cutoff_color = ifelse(pval < cutoff_pvalue, as.character(chr_b37), NA))


ggplot(dat_plot, aes(x=chr_b37, y=-log10(pval))) +
  geom_jitter(aes(color=cutoff_color), size=3, width=0.3, height=0, show.legend=FALSE) +
  scale_color_manual(values= pal_d3("category20", alpha=0.6)(20)) +
  geom_label_repel(aes(label=labels_choose, color=cutoff_color),
                   size=3, box.padding=2, point.padding=0.5, segment.size=0.5, show.legend=FALSE, max.overlaps=100) +
  geom_hline(yintercept=-log10(cutoff_pvalue), colour="#990000", linetype="dashed") +
  annotate("text", x=22, y=-log10(cutoff_pvalue) + 0.3,
           label=paste0("-log10(pval) = ", round(-log10(cutoff_pvalue), digits=2)), size=3) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 12)) +
  theme_minimal() +
  theme(panel.grid=element_blank(), axis.line=element_line(colour='black', linewidth=0.6),
        plot.title=element_text(size=13, hjust=0.5), plot.subtitle=element_text(size=9, hjust=0.5),
        axis.text.x=element_text(size=13, angle=0, hjust=1),
        axis.text.y=element_text(size=13, angle=0, hjust=1),
        axis.title=element_text(size=13)) +
  labs(y="-log10(Pval)", x="", title="MR results with the discovery database")
dat_plot=dat
# 增加注释数据
dat_plot$OR <- round(dat$or,digits = 2)
dat_plot$LCI <- round(dat$or_lci95,digits = 2)
dat_plot$UCI <- round(dat$or_uci95,digits = 2)
dat_plot$labels <- ifelse((dat_plot$pval < cutoff_pvalue) ,
                          paste0(dat_plot$Gene,
                                 "\n OR: ",dat_plot$OR,
                                 "; 95%CI: [",dat_plot$LCI,"-",dat_plot$UCI,
                                 "]" ) ,"")
dat_plot$name <- paste0(dat_plot$Gene)
dat_plot$labels_choose <- ifelse(dat_plot$pval<cutoff_pvalue,dat_plot$labels,"")

dat_plot$chr_b37 <- factor(dat_plot$chr_b37,levels = c(1:22,"X"))
set.seed(124)
ggplot(dat_plot, aes(x=chr_b37, y=-log10(pval))) +
  geom_jitter(aes(color=chr_b37),size=3,
              width = 0.3,height = 0, show.legend = c(color = FALSE) )+
  scale_shape_manual(values = c(1,16))+
  scale_size_manual(values = c(2,3,4))+
  scale_alpha_continuous(range = c(0,0.8))+
  scale_color_manual(values=c(pal_d3("category20")(20),pal_d3("category20")(3))) +
  
  geom_label_repel( aes(label=labels_choose,colour =chr_b37),
                    size=4,
                    box.padding = 2,
                    point.padding = 0.5,
                    segment.size = 0.5,
                    show.legend = F,
                    max.overlaps = 100)+
  geom_hline(aes(yintercept=-log10(cutoff_pvalue)), colour="#990000", linetype="dashed")+
  annotate("text", x=22-0.5, y=-log10(cutoff_pvalue)+0.3,
           label=paste0("-log10(pval) = ",round(-log10(cutoff_pvalue),digits = 2)),size=4)+
  scale_y_continuous(expand = c(0,0),limits = c(0,12)) +
  
  theme(panel.grid = element_blank(),

        legend.key = element_blank(),
        axis.line = element_line(colour = 'black', linewidth = 0.6),
        panel.background = element_blank(),
        plot.title = element_text(size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5),
        axis.text.x  = element_text(size = 13, color = 'black',angle = 0,hjust = 1),
        axis.text.y  = element_text(size = 13, color = 'black',angle = 0,hjust = 1),
        
        axis.title = element_text(size = 13, color = 'black'))  +
  labs(y="-log10(Pval)" ,x="",title = "MR results with the discovery database")
ggsave("Figure2.discovery_Manhattan_plot.tif",width=16,height=12)




