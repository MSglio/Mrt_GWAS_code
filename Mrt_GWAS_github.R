###  all the data is from GitHub: https://github.com/MSglio/Mrt_GWAS_code

## required libraries
library(GridLMM) 
library(dplyr)
library(ggplot2)
library(rptR)
library(readxl)
require(data.table)
library(lme4)
library(devtools)
library(ggpubr)
library(readr)
require(rrBLUP)

## Mrt phenotyping data  -----------------------------------------------

Mrt_ALL=read_excel("../Mrt_PAPER_R.xlsx")

Mrt_ALL=Mrt_ALL[,-1]
pheno_keep=Mrt_ALL$isotype
pheno_keep = cbind(Mrt_ALL$isotype, Mrt_ALL$isotype, 1)
fwrite(pheno_keep, file="pheno.keep", col.names=F, sep=" ")

Mrt_ALL=melt(as.data.table(Mrt_ALL), 1:2)

miss=is.na(Mrt_ALL$value)
Mrt_ALL$value=as.numeric(Mrt_ALL$value)

X=is.na(Mrt_ALL$value)&!miss

Mrt_ALL[X,]

Mrt_ALL$obs=tstrsplit(Mrt_ALL$variable, "_")[[1]]
Mrt_ALL$bloc=tstrsplit(Mrt_ALL$variable, "_")[[2]]
Mrt_ALL$rep=tstrsplit(Mrt_ALL$variable, "_")[[3]]
Mrt_ALL$bloceffect=paste(Mrt_ALL$bloc, Mrt_ALL$obs)


#### MODEL - blup model --------------------------------------------------------------

fit_all= glmer(value~(1|bloc)+(1|isotype), Mrt_ALL, family="poisson")

pheno_blups=as.data.frame(ranef(fit_all)$isotype)

## check distrution of random effects
re <- ranef(fit_all)$isotype
qqnorm(re$`(Intercept)`)
qqline(re$`(Intercept)`)

names(pheno_blups)="Mrt_blup"
qplot(pheno_blups$Mrt_blup)

pheno_blups$isotype=row.names(pheno_blups)

#### built genotype matrix with hard filtered SNPs from CENDR vcf processed with PLINK ------------------------------------------------------

### with the restricted set of non Pacific 113 strains
map = fread("../GWAS_CENDR2020_hardfilterisotype_LD099_ALL.map", header = F)
keep = fread("../GWAS_CENDR2020_hardfilterisotype_LD099_ALL.prune.in", header = F)
geno = fread("../GWAS_CENDR2020_hardfilterisotype_LD099_ALL.ped", header = F)

### with the set of 132 strains
#map = fread("../GWAS_CENDR2020_hardfilter_132isotype_LD099_ALL.map", header = F)
#keep = fread("../GWAS_CENDR2020_hardfilter_132isotype_LD099_ALL.prune.in", header = F)
#geno = fread("../GWAS_CENDR2020_hardfilter_132isotype_LD099_ALL.ped", header = F)

head(geno[,1:8])
unique(geno$V1)
samples=geno$V1
x = as.matrix(geno[,-(1:6)])
x = x[,seq(1, ncol(x), 2)]
geno = x[,map$V2 %in% keep$V1]


# GWAS code ---------------------------------------------------------------

SNPs = map[map$V2 %in% keep$V1,c(1,4)]
names(SNPs) = c("chr", "pos")

rownames(geno)=samples
SNPs$chr=as.character(SNPs$chr)
class(SNPs$chr)
SNPs$chr[SNPs$chr==23] <- "X"
SNPs$chr <- factor(SNPs$chr, levels = c( "I", "II", "III", "IV", "V", "X", "MtDNA"))

## 
pheno=pheno_blups

# line up pheno and geno strains
pheno = subset(pheno_blups, isotype %in% samples)
pheno = pheno[order(pheno$isotype),]
geno = geno[samples %in% pheno$isotype,]
geno = geno[order(rownames(geno)),]

# create additive relationship matrix
A = A.mat(geno)

# fit linear mixed model at each snp --------------------------------------

fit1 = GridLMM_GWAS(Mrt_blup~(1|isotype), test_formula=Mrt_blup~1, reduced_formula = Mrt_blup~1, data=pheno, X=geno, X_ID="isotype", relmat = list(isotype=A))

summary(fit1)

SNPs_hf = cbind(SNPs, beta=fit1$results$beta.2, p= -log10(fit1$results$p_value_REML))

# significance threshold - permutations ------------------------------------------------------------

permuts = function(phen, X, nperm=100, np=4 ) {
  library(parallel)
  unlist(mclapply(1:nperm, mc.cores=np, function(i) {
    phen$isotype = factor(phen$isotype, labels=sample(unique(phen$isotype)))
    fit = GridLMM_GWAS(Mrt_blup~(1|isotype), test_formula=Mrt_blup~1, reduced_formula = Mrt_blup~1, data=phen, X=X, X_ID="isotype", relmat = list(isotype=A))
    fit$results$p_value_REML[which.min(fit$results$p_value_REML)]
  }))
}

## adjust the number of nperm to generate the threshold 
#null=permuts(pheno, geno, nperm=1000)
quantile(null, 0.2)
threshold = -log10(quantile(null, 0.2))

## GWAS with hard-filtered snp from isotypes 113 strains, nperm = 1000
#threshold = 5.521903

PAPER_Fig1_E = ggplot(SNPs_hf, aes(pos/1e6, p, factor=chr)) +
  geom_point(size=0.5) +
  facet_grid(.~chr, scale='free') +
  geom_hline(yintercept=threshold, color = "darkred" ) +
  scale_y_continuous(name="-log10(p-value)", limits=c(0, 8)) +
  #ggtitle("Genome Wide Association Study on a subset of 113 wild isolates") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("SNP position along the chromosome (Mb)") +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7)) +
  theme(legend.title = element_text(colour="blue", size=10, face="bold")) +
  theme(axis.text = element_text(size = 8))+
  theme(axis.title = element_text(size = 8))+
  theme(plot.title = element_text(size = 10))
PAPER_Fig1_E


###### FINE MAPPING - GWAS code with ALL the snps region mostly associated on chr III (4 to 7 Mb) ---------------------------------------------

require(data.table)

## with all imputed variants from CENDR 113 strains
map = fread("../impute_isotype_maf005.map", header = F)
geno = fread("../impute_isotype_maf005.ped", header = F)

subset_snps_sf=map[,c(1,4)]
names(subset_snps_sf) = c("chr", "pos")

head(geno[,1:8])
samples=geno$V1
x = as.matrix(geno[,-(1:6)])
x = x[,seq(1, ncol(x), 2)]
geno = x[,map$V2 %in% map$V2]


##### GWAS code
SNPs_sf=subset_snps_sf
dim(subset_snps_sf)
rownames(geno)=samples

## define phenotype
pheno = pheno_blups

# line up phenotype and genotype matrix
pheno = subset(pheno, isotype %in% samples)
pheno = pheno[order(pheno$isotype),]
class(pheno$Mrt_blup)
pheno$Mrt_blup=as.numeric(as.character(pheno$Mrt_blup))

geno = geno[samples %in% pheno$isotype,]
geno = geno[order(rownames(geno)),]


# create additive relationship matrix

A = A.mat(geno)

### model
fit_subset = GridLMM_GWAS(Mrt_blup~(1|isotype), test_formula=Mrt_blup~1, reduced_formula = Mrt_blup~1, data=pheno, X=geno, X_ID="isotype", relmat = list(isotype=A))

SNPs_sf = cbind(SNPs_sf, beta=fit_subset$results$beta.2, p= -log10(fit_subset$results$p_value_REML))

SNPs_sf$chr=as.character(SNPs_sf$chr)
class(SNPs_sf$chr)
SNPs_sf$chr[SNPs_sf$chr==23] <- "X"

SNPs_sf$chr <- factor(SNPs_sf$chr, levels = c( "I", "II", "III", "IV", "V", "X", "MtDNA"))


plot_fit = ggplot(SNPs_sf, aes(pos/1e6, p, factor=chr)) +
  #levels 
  geom_point() +
  facet_grid(.~chr, scale='free') +
  #geom_hline(yintercept=threshold, color = "darkred" ) +
  scale_y_continuous(name="p", limits=c(0, 15))
plot_fit
