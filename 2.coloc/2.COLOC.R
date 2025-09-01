### COLOC ----
{
library(data.table)
library(dplyr)
library(gwasglue)
library(coloc)
library(gwasvcf)
library(TwoSampleMR)
library(locuscomparer)
library(ggplot2)
}
gwas2 <-  data.table::fread("./vit.csv",data.table = F,fill = T)
gwas1 <- data.table::fread("./CTSS.vcf.gz.csv",data.table = F,fill = T)
i <- 1
gwas1_c <- gwas1[gwas1$chr == i,]
gwas1_c <- gwas1_c[gwas1_c$pos > 150702672- 300000 &  gwas1_c$pos < 150738433+ 300000,]
gwas2_c <- gwas2[gwas2$SNP %in% gwas1_c$SNP,]
gwas2_c<-gwas2_c[,c("SNP","A1","A2","MAF","P","beta","SE")]
colnames(gwas2_c)=c("SNP","effect allele","other allele","eaf","pval","beta","se")
dat_merge <- merge(gwas1_c,gwas2_c,by=c("SNP"),suffixes = c("_gwas1","_gwas2"))
dat_merge=dat_merge[!duplicated(dat_merge$SNP),]
dat_merge$eaf_gwas1   <- ifelse(dat_merge$eaf_gwas1=='0',0.5,dat_merge$eaf_gwas1)
dat_merge$eaf_gwas2   <- ifelse(dat_merge$eaf_gwas2=='0',0.5,dat_merge$eaf_gwas2)
dat_merge <- dat_merge[rowSums(is.infinite(dat_merge)) == 0, ]
dat_merge<- dat_merge[!is.infinite(dat_merge$beta_gwas1), ]
gwas1_form2 <- list(snp = dat_merge$SNP,
                    position = dat_merge$pos,
                    beta = dat_merge$beta_gwas1,
                    varbeta = dat_merge$se_gwas1 ^ 2,
                    MAF = dat_merge$eaf_gwas1,
                    type = "quant",
                    N = dat_merge$samplesize
                    
)
check_dataset(gwas1_form2)
gwas2_form1<- list(snp = dat_merge$SNP,
                   position =dat_merge$pos,
                   beta = dat_merge$beta_gwas2,
                   varbeta = dat_merge$se_gwas2^2,
                   type = "cc")
check_dataset(gwas2_form1)
my.res1 <- coloc.abf(dataset1 =gwas1_form2 ,dataset2 = gwas2_form1)
subset(my.res1$results,SNP.PP.H4>0.8)
need_result=my.res1$results %>% filter(SNP.PP.H4 > 0.8)
gwas_fn=dat_merge[,c('SNP','pval_gwas2')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas2)
eqtl_fn=dat_merge[,c('SNP','pval_gwas1')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas1)
pdf(file="CTSS.rs2867300.pdf",width=8,height=6,onefile = FALSE)
print(locuscompare(in_fn1=gwas_fn,
                   in_fn2=eqtl_fn,
                   title1='vitiligo GWAS',
                   title2='CTSS eQTL',
                   snp = "rs2867300"))


dev.off()
############
gwas1 <- data.table::fread("./EHMT2.vcf.gz.csv",data.table = F,fill = T)
i <- 6
gwas1_c <- gwas1[gwas1$chr == i,]
gwas1_c <- gwas1_c[gwas1_c$pos > 31847536- 300000 &  gwas1_c$pos < 31865464+ 300000,]
gwas2_c <- gwas2[gwas2$SNP %in% gwas1_c$SNP,]
gwas2_c<-gwas2_c[,c("SNP","A1","A2","MAF","P","beta","SE")]
colnames(gwas2_c)=c("SNP","effect allele","other allele","eaf","pval","beta","se")
dat_merge <- merge(gwas1_c,gwas2_c,by=c("SNP"),suffixes = c("_gwas1","_gwas2"))
dat_merge=dat_merge[!duplicated(dat_merge$SNP),]
dat_merge$eaf_gwas1   <- ifelse(dat_merge$eaf_gwas1=='0',0.5,dat_merge$eaf_gwas1)
dat_merge$eaf_gwas2   <- ifelse(dat_merge$eaf_gwas2=='0',0.5,dat_merge$eaf_gwas2)
dat_merge <- dat_merge[rowSums(is.infinite(dat_merge)) == 0, ]
dat_merge<- dat_merge[!is.infinite(dat_merge$beta_gwas1), ]
gwas1_form2 <- list(snp = dat_merge$SNP,
                    position = dat_merge$pos,
                    beta = dat_merge$beta_gwas1,
                    varbeta = dat_merge$se_gwas1 ^ 2,
                    MAF = dat_merge$eaf_gwas1,
                    type = "quant",
                    N = dat_merge$samplesize
                    
)
check_dataset(gwas1_form2)
gwas2_form1<- list(snp = dat_merge$SNP,
                   position =dat_merge$pos,
                   beta = dat_merge$beta_gwas2,
                   varbeta = dat_merge$se_gwas2^2,
                   type = "cc")
check_dataset(gwas2_form1)
my.res1 <- coloc.abf(dataset1 =gwas1_form2 ,dataset2 = gwas2_form1)
subset(my.res1$results,SNP.PP.H4>0.8)
need_result=my.res1$results %>% filter(SNP.PP.H4 > 0.8)
gwas_fn=dat_merge[,c('SNP','pval_gwas2')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas2)
eqtl_fn=dat_merge[,c('SNP','pval_gwas1')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas1)
pdf(file="EHMT2.rs386480.pdf",width=8,height=6,onefile = FALSE)
print(locuscompare(in_fn1=gwas_fn,
                   in_fn2=eqtl_fn,
                   title1='vitiligo GWAS',
                   title2='EHMT2 eQTL',
                   snp = "rs386480"))


dev.off()
############
gwas1 <- data.table::fread("./FYN.vcf.gz.csv",data.table = F,fill = T)
i <- 6
gwas1_c <- gwas1[gwas1$chr == i,]
gwas1_c <- gwas1_c[gwas1_c$pos > 111981535- 300000 &  gwas1_c$pos < 112194655+ 300000,]
gwas2_c <- gwas2[gwas2$SNP %in% gwas1_c$SNP,]
gwas2_c<-gwas2_c[,c("SNP","A1","A2","MAF","P","beta","SE")]
colnames(gwas2_c)=c("SNP","effect allele","other allele","eaf","pval","beta","se")
dat_merge <- merge(gwas1_c,gwas2_c,by=c("SNP"),suffixes = c("_gwas1","_gwas2"))
dat_merge=dat_merge[!duplicated(dat_merge$SNP),]
dat_merge$eaf_gwas1   <- ifelse(dat_merge$eaf_gwas1=='0',0.5,dat_merge$eaf_gwas1)
dat_merge$eaf_gwas2   <- ifelse(dat_merge$eaf_gwas2=='0',0.5,dat_merge$eaf_gwas2)
dat_merge <- dat_merge[rowSums(is.infinite(dat_merge)) == 0, ]
dat_merge<- dat_merge[!is.infinite(dat_merge$beta_gwas1), ]
gwas1_form2 <- list(snp = dat_merge$SNP,
                    position = dat_merge$pos,
                    beta = dat_merge$beta_gwas1,
                    varbeta = dat_merge$se_gwas1 ^ 2,
                    MAF = dat_merge$eaf_gwas1,
                    type = "quant",
                    N = dat_merge$samplesize
                    
)
check_dataset(gwas1_form2)
gwas2_form1<- list(snp = dat_merge$SNP,
                   position =dat_merge$pos,
                   beta = dat_merge$beta_gwas2,
                   varbeta = dat_merge$se_gwas2^2,
                   type = "cc")
check_dataset(gwas2_form1)
my.res1 <- coloc.abf(dataset1 =gwas1_form2 ,dataset2 = gwas2_form1)
subset(my.res1$results,SNP.PP.H4>0.8)
need_result=my.res1$results %>% filter(SNP.PP.H4 > 0.8)
gwas_fn=dat_merge[,c('SNP','pval_gwas2')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas2)
eqtl_fn=dat_merge[,c('SNP','pval_gwas1')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas1)
pdf(file="FYN.rs1474466.pdf",width=8,height=6,onefile = FALSE)
print(locuscompare(in_fn1=gwas_fn,
                   in_fn2=eqtl_fn,
                   title1='vitiligo GWAS',
                   title2='FYN eQTL',
                   snp = "rs1474466"))


dev.off()
##################
gwas1 <- data.table::fread("./ERBB3.vcf.gz.csv",data.table = F,fill = T)
i <- 12
gwas1_c <- gwas1[gwas1$chr == i,]
gwas1_c <- gwas1_c[gwas1_c$pos > 56473641- 300000 &  gwas1_c$pos < 56497289+ 300000,]
gwas2_c <- gwas2[gwas2$SNP %in% gwas1_c$SNP,]
gwas2_c<-gwas2_c[,c("SNP","A1","A2","MAF","P","beta","SE")]
colnames(gwas2_c)=c("SNP","effect allele","other allele","eaf","pval","beta","se")
dat_merge <- merge(gwas1_c,gwas2_c,by=c("SNP"),suffixes = c("_gwas1","_gwas2"))
dat_merge=dat_merge[!duplicated(dat_merge$SNP),]
dat_merge$eaf_gwas1   <- ifelse(dat_merge$eaf_gwas1=='0',0.5,dat_merge$eaf_gwas1)
dat_merge$eaf_gwas2   <- ifelse(dat_merge$eaf_gwas2=='0',0.5,dat_merge$eaf_gwas2)
dat_merge <- dat_merge[rowSums(is.infinite(dat_merge)) == 0, ]
dat_merge<- dat_merge[!is.infinite(dat_merge$beta_gwas1), ]
gwas1_form2 <- list(snp = dat_merge$SNP,
                    position = dat_merge$pos,
                    beta = dat_merge$beta_gwas1,
                    varbeta = dat_merge$se_gwas1 ^ 2,
                    MAF = dat_merge$eaf_gwas1,
                    type = "quant",
                    N = dat_merge$samplesize
                    # sdY = 1.2
)
check_dataset(gwas1_form2)
gwas2_form1<- list(snp = dat_merge$SNP,
                   position =dat_merge$pos,
                   beta = dat_merge$beta_gwas2,
                   varbeta = dat_merge$se_gwas2^2,
                   type = "cc")
check_dataset(gwas2_form1)
my.res1 <- coloc.abf(dataset1 =gwas1_form2 ,dataset2 = gwas2_form1)
subset(my.res1$results,SNP.PP.H4>0.8)
need_result=my.res1$results %>% filter(SNP.PP.H4 > 0.8)

gwas_fn=dat_merge[,c('SNP','pval_gwas2')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas2)
eqtl_fn=dat_merge[,c('SNP','pval_gwas1')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas1)
pdf(file="ERBB3.rs11171739.pdf",width=8,height=6,onefile = FALSE)
print(locuscompare(in_fn1=gwas_fn,
                   in_fn2=eqtl_fn,
                   title1='vitiligo GWAS',
                   title2='ERBB3 eQTL',
                   snp = "rs11171739"))


dev.off()
##################
gwas1 <- data.table::fread("./PPP3CA.vcf.gz.csv",data.table = F,fill = T)
i <- 4
gwas1_c <- gwas1[gwas1$chr == i,]
gwas1_c <- gwas1_c[gwas1_c$pos > 101944566- 300000 &  gwas1_c$pos < 102269435+ 300000,]
gwas2_c <- gwas2[gwas2$SNP %in% gwas1_c$SNP,]
gwas2_c<-gwas2_c[,c("SNP","A1","A2","MAF","P","beta","SE")]
colnames(gwas2_c)=c("SNP","effect allele","other allele","eaf","pval","beta","se")
dat_merge <- merge(gwas1_c,gwas2_c,by=c("SNP"),suffixes = c("_gwas1","_gwas2"))
dat_merge=dat_merge[!duplicated(dat_merge$SNP),]
dat_merge$eaf_gwas1   <- ifelse(dat_merge$eaf_gwas1=='0',0.5,dat_merge$eaf_gwas1)
dat_merge$eaf_gwas2   <- ifelse(dat_merge$eaf_gwas2=='0',0.5,dat_merge$eaf_gwas2)
dat_merge <- dat_merge[rowSums(is.infinite(dat_merge)) == 0, ]
dat_merge<- dat_merge[!is.infinite(dat_merge$beta_gwas1), ]
gwas1_form2 <- list(snp = dat_merge$SNP,
                    position = dat_merge$pos,
                    beta = dat_merge$beta_gwas1,
                    varbeta = dat_merge$se_gwas1 ^ 2,
                    MAF = dat_merge$eaf_gwas1,
                    type = "quant",
                    N = dat_merge$samplesize
                    
)
check_dataset(gwas1_form2)

gwas2_form1<- list(snp = dat_merge$SNP,
                   position =dat_merge$pos,
                   beta = dat_merge$beta_gwas2,
                   varbeta = dat_merge$se_gwas2^2,
                   type = "cc")
check_dataset(gwas2_form1)
my.res1 <- coloc.abf(dataset1 =gwas1_form2 ,dataset2 = gwas2_form1)

subset(my.res1$results,SNP.PP.H4>0.8)

need_result=my.res1$results %>% filter(SNP.PP.H4 > 0.8)

gwas_fn=dat_merge[,c('SNP','pval_gwas2')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas2)
eqtl_fn=dat_merge[,c('SNP','pval_gwas1')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas1)
pdf(file="PPP3CA.rs2850368.pdf",width=8,height=6,onefile = FALSE)
print(locuscompare(in_fn1=gwas_fn,
                   in_fn2=eqtl_fn,
                   title1='vitiligo GWAS',
                   title2='PPP3CA eQTL',
                   snp = "rs2850368"))


dev.off()
##################
gwas1 <- data.table::fread("./TBXAS1.vcf.gz.csv",data.table = F,fill = T)
i <- 7
gwas1_c <- gwas1[gwas1$chr == i,]
gwas1_c <- gwas1_c[gwas1_c$pos > 139476850- 300000 &  gwas1_c$pos < 139720125+ 300000,]
gwas2_c <- gwas2[gwas2$SNP %in% gwas1_c$SNP,]
gwas2_c<-gwas2_c[,c("SNP","A1","A2","MAF","P","beta","SE")]
colnames(gwas2_c)=c("SNP","effect allele","other allele","eaf","pval","beta","se")
dat_merge <- merge(gwas1_c,gwas2_c,by=c("SNP"),suffixes = c("_gwas1","_gwas2"))
dat_merge=dat_merge[!duplicated(dat_merge$SNP),]
dat_merge$eaf_gwas1   <- ifelse(dat_merge$eaf_gwas1=='0',0.5,dat_merge$eaf_gwas1)
dat_merge$eaf_gwas2   <- ifelse(dat_merge$eaf_gwas2=='0',0.5,dat_merge$eaf_gwas2)
dat_merge <- dat_merge[rowSums(is.infinite(dat_merge)) == 0, ]
dat_merge<- dat_merge[!is.infinite(dat_merge$beta_gwas1), ]
gwas1_form2 <- list(snp = dat_merge$SNP,
                    position = dat_merge$pos,
                    beta = dat_merge$beta_gwas1,
                    varbeta = dat_merge$se_gwas1 ^ 2,
                    MAF = dat_merge$eaf_gwas1,
                    type = "quant",
                    N = dat_merge$samplesize
                    
)
check_dataset(gwas1_form2)
gwas2_form1<- list(snp = dat_merge$SNP,
                   position =dat_merge$pos,
                   beta = dat_merge$beta_gwas2,
                   varbeta = dat_merge$se_gwas2^2,
                   type = "cc")
check_dataset(gwas2_form1)
my.res1 <- coloc.abf(dataset1 =gwas1_form2 ,dataset2 = gwas2_form1)
subset(my.res1$results,SNP.PP.H4>0.8)
need_result=my.res1$results %>% filter(SNP.PP.H4 > 0.8)
gwas_fn=dat_merge[,c('SNP','pval_gwas2')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas2)
eqtl_fn=dat_merge[,c('SNP','pval_gwas1')]%>%
  dplyr::rename(rsid=SNP,pval=pval_gwas1)
pdf(file="TBXAS1.rs2286197.pdf",width=8,height=6,onefile = FALSE)
print(locuscompare(in_fn1=gwas_fn,
                   in_fn2=eqtl_fn,
                   title1='vitiligo GWAS',
                   title2='TBXAS1 eQTL',
                   snp = "rs2286197"))


dev.off()
