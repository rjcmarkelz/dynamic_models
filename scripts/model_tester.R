# Y ~ Mu + Q1 + Q2 + Q3 + Q3:Q2 + error

library(qtl)
library(qtlbim)
#make sure directory is set
setwd("~/git.repos/brassica_meta_analysis/Cleaned_data/")

field_traits <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
                         phefile="all_traits_RQTL.csv", genotypes=c("AA","BB"), na.strings = c("-", "NA"))
head(field_traits)
plot(field_traits)

brassica_traits_qb <- field_traits
brassica_traits_qb <- qb.genoprob(brassica_traits_qb, step=2, stepwidth = "variable")
summary(brassica_traits_qb)

?qb.mcmc
brassica_traits_qb_leaflength <- qb.mcmc(brassica_traits_qb, pheno.col = 2, rancov = 1, seed = 1616, epistasis = TRUE)
#3000 iterations
#3000 iterations
plot(brassica_traits_qb_leaflength)

?qb.hpdone
other_HDI <- qb.hpdone(brassica_traits_qb_leaflength, effects = "estimate")
str(other_HDI)
other_HDI
hist(brassica_traits_qb$pheno$LeafLnUN)
ph_names <- names(brassica_traits_qb$pheno)
ph_names

qbModeltest <- qb.model(brassica_traits_qb, epistatis = TRUE, main.nqtl = 5, intcov = 1,
	                    interval = NULL)
qbModeltest


