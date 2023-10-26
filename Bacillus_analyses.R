setwd("Bacillus/Final_for_paper/Analyses/")

### Loading packages
library(vcfR) # version 1.12.0
library(adegenet) # version 2.1.3
library(MASS) # version 7.3-54
library(phytools) # version 0.7-80
library(ridgeline)


### Custom functions 
'%!in%' <- Negate('%in%')

prop.het <- function(vector) return(sum(substr(vector, 1, 1) != substr(vector, 3, 3), na.rm=T)/sum(is.na(vector)==F)) 
# Returns the proportion of sites that are heterozygous in a genotype file obtained by
# extract.gt(read.vcfr(VCF_FILE), element = "GT")

count.na <- function(x) sum(is.na(x)) # counts the number of NA in a vector
count.not.na <-function(x) sum(is.na(x) == F) # counts the number of non-missing values in a vector

standardize <- function(vector) return(vector / mean(vector, na.rm=T)) # standardizes values of a vector (used here for depth)
# by dividing each value to the average of non-missing values in the vector. This results in a mean of 1.

#######################################
### Reading and formatting the data ###
#######################################

scf_lengths <- read.csv("Brsri_v3_scf_lengths.csv", h=T, sep=";") # This is used for homemade Manhattan plots
scf_lengths$end <- scf_lengths$start + scf_lengths$length # defines the position of each chromosome

sicily.vcf <- read.vcfR("Sicily_filtered_gatk_an23_ql20_snps_dp08_mc3_miss50.recode.vcf") # Reads the main vcf file with all field-caught individuals
sicily.gt <- extract.gt(sicily.vcf, element = "GT", convertNA = T) # extracts the genotype value
sicily.dp <- extract.gt(sicily.vcf, element = "DP", as.numeric = T, convertNA = T) # extracts depth value for each genotype
sicily_std.dp <- apply(sicily.dp, 2, standardize) # "standardizes" depth values for each individual, i.e. the average depth value for each individual becomes 1.
sicily.gl <- vcfR2genlight(sicily.vcf) # converts into genlight format (adegenet)
# excluded 5564 sites because they were triallelic

meta <- read.csv("TableS1.csv", sep=",", h=T)

tree.not.pruned <- read.tree("../mtDNA/All_filtered_with_refs_aligned.phylip_rerooted.newick") # Reads the mtDNA phylogeny 
# retain only tips corresponding to individuals included in sicily.vcf (removing references and individuals with too much missing RAD data)
tree <- drop.tip(phy = tree.not.pruned, tip = tree.not.pruned$tip.label[tree.not.pruned$tip.label %in% meta$ID == F], trim.internal = T)

#######################
### mtDNA phylogeny ###
#######################

## Attributing mtDNA haplotypes to individuals
# Using getDescendants() written by Liam Revell
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr <- vector()
  daughters <- tree$edge[which(tree$edge[,1] == node), 2]
  curr <- c(curr,daughters)
  w <- which(daughters >= length(tree$tip))
  if(length(w) > 0) for(i in 1:length(w))
    curr <- getDescendants(tree, daughters[w[i]], curr)
  return(curr)
}

plot(tree.not.pruned, cex=0.25)
add.scale.bar()
nodelabels(cex=0.4, frame = "none")

## See Figure S3
atticus.mt <- tree.not.pruned$tip.label[getDescendants(tree.not.pruned, node = 1129)]
atticus.mt <- atticus.mt[-which(is.na(atticus.mt))]
ggrandii.mt <- tree.not.pruned$tip.label[getDescendants(tree.not.pruned, node = 1063)]
ggrandii.mt <- ggrandii.mt[-which(is.na(ggrandii.mt))]
gbenazzii.mt <- tree.not.pruned$tip.label[getDescendants(tree.not.pruned, node = 998)]
gbenazzii.mt <- gbenazzii.mt[-which(is.na(gbenazzii.mt))]
gmaretimi.mt <- tree.not.pruned$tip.label[getDescendants(tree.not.pruned, node = 1050)]
gmaretimi.mt <- gmaretimi.mt[-which(is.na(gmaretimi.mt))]
allrossius.mt <- tree.not.pruned$tip.label[getDescendants(tree.not.pruned, node = 576)]
allrossius.mt <- allrossius.mt[-which(is.na(allrossius.mt))]
rrossius.mt <- tree.not.pruned$tip.label[getDescendants(tree.not.pruned, node = 900)]
rrossius.mt <- rrossius.mt[-which(is.na(rrossius.mt))]
rredtenbacheri.mt <- allrossius.mt[which(allrossius.mt %!in% rrossius.mt)]

meta$col.mt <- "grey"
meta$col.mt[which(meta$ID %in% atticus.mt)] <- "#71BF44"
meta$col.mt[which(meta$ID %in% ggrandii.mt)] <- "#BE1E2D"
meta$col.mt[which(meta$ID %in% gbenazzii.mt)] <- "#F7941D"
meta$col.mt[which(meta$ID %in% gmaretimi.mt)] <- "#F9ED32"
meta$col.mt[which(meta$ID %in% rrossius.mt)] <- "#2E368F"
meta$col.mt[which(meta$ID %in% rredtenbacheri.mt)] <- "#0599CE"


#############################
### First MDS on all loci ###
#############################

# Compute Euclidian distances between all individuals
# Warning: triallelic sites have been excluded already
dist_sicily <- dist(sicily.gl)

MDS_sicily <- isoMDS(d = dist_sicily, k = 10, maxit = 100) # Run the MDS. The isoMDS function comes from the package MASS

par(mfrow=c(2,1))
plot(MDS_sicily$points[,2] ~ MDS_sicily$points[,1], las=1, pch=19, col=meta$col.mt, xlab="Dimension 1", ylab="Dimension 2", main="MDS - all markers")
plot(MDS_sicily$points[,3] ~ MDS_sicily$points[,1], las=1, pch=19, col=meta$col.mt, xlab="Dimension 1", ylab="Dimension 3", main="MDS - all markers")

### assigning species identity based on MDS results
meta$sp.nucl <- "unknown"
meta$sp.nucl[which(MDS_sicily$points[,1] < (-200))] <- "rredtenbacheri"
meta$sp.nucl[which(MDS_sicily$points[,2] < (-200))] <- "gbenazzii"
meta$sp.nucl[which(MDS_sicily$points[,2] > 120)] <- "ggrandii"
meta$sp.nucl[which(MDS_sicily$points[,2] < (-100) & MDS_sicily$points[,2] > (-200) & MDS_sicily$points[,1] > 100)] <- "gmaretimi"
meta$sp.nucl[which(MDS_sicily$points[,3] > 250)] <- "atticus"
meta$sp.nucl[which(MDS_sicily$points[,2] > 50 & MDS_sicily$points[,2] < 120 & MDS_sicily$points[,1] > (-100) & MDS_sicily$points[,1] < 100  & MDS_sicily$points[,3]< 50)] <- "SEhybrids"
# SEhybrids contains both B. whitei and hybridogenetic lineage 1. See below for identification.
meta$sp.nucl[which(MDS_sicily$points[,2] > 50 & MDS_sicily$points[,2] < 120 & MDS_sicily$points[,1] > (-100) & MDS_sicily$points[,1] < 100  & MDS_sicily$points[,3]> 50)] <- "lynceorum"
meta$sp.nucl[which(MDS_sicily$points[,2] < (-100) & MDS_sicily$points[,1] < 0)] <- "hybridoNW"

# Assigning colours for plotting
meta$col.nucl <- "grey"
meta$col.nucl[which(meta$sp.nucl == "gbenazzii")] <- "#F7941D"
meta$col.nucl[which(meta$sp.nucl == "gmaretimi")] <- "#F9ED32"
meta$col.nucl[which(meta$sp.nucl == "ggrandii")] <- "#BE1E2D"
meta$col.nucl[which(meta$sp.nucl == "hybridoNW")] <- "#F47F72"
meta$col.nucl[which(meta$sp.nucl == "lynceorum")] <- "#8B5E3C"
meta$col.nucl[which(meta$sp.nucl == "SEhybrids")] <- "#92278F"
meta$col.nucl[which(meta$sp.nucl == "atticus")] <- "#71BF44"
meta$col.nucl[which(meta$sp.nucl == "rredtenbacheri")] <- "#0599CE"
meta$col.nucl[which(meta$sp.nucl == "rrossius")] <- "#2E368F"


#############################
### whitei / hybrido 1 ID ###
#############################

eggs_K2 <- read.table("eggs_cleaned.Q.2")
eggs_meta <- read.table("popmap_eggs.txt")
eggs.vcf <- read.vcfR("cleaned.vcf")
names <- dimnames(eggs.vcf@gt)[[2]][-1]
row.names(eggs_K2) <- eggs_meta$V1
barplot(t(as.matrix(eggs_K2[which(eggs_meta$V2=="-"),])), col=c("#0697CD", "#BC202E"), las=2, border=NA, space=0)


mothers.vcf <- read.vcfR("females_plus_SEhybrids_DP5_miss75_mac20_Q20.recode.vcf")
mothers.gt <- extract.gt(mothers.vcf, "GT")
mothers.gl <- vcfR2genlight(mothers.vcf)

meta_all <- read.csv("../metadata.csv", sep=";", h=T)
metamoms <- data.frame(ID = dimnames(mothers.gt)[[2]])
metamoms$pch.year <- 1
metamoms$pch.year[which(metamoms$ID %in% meta_all$ID[which(meta_all$year == 2018)])] <- 19
metamoms$pch.year[which(metamoms$ID %in% meta_all$ID[which(meta_all$year == 2020)])] <- 15

dist_moms <- dist(mothers.gl)
MDS_moms <- isoMDS(d = dist_moms, k = 10, maxit = 100)

# This is part of Figure 2C
plot(MDS_moms$points[,2] ~ MDS_moms$points[,1], pch=metamoms$pch.year, las=1, xlab="Dimension 1", ylab="Dimension 2", main="")
text(labels = dimnames(mothers.gt)[[2]] , x=MDS_moms$points[,1], y=MDS_moms$points[,2], cex=0.5)

meta$sp.nucl[which(meta$sp.nucl == "SEhybrids" & meta$ID %in% rownames(MDS_moms$points)[which(MDS_moms$points[,2] < 5)])] <- "whitei"
meta$sp.nucl[which(meta$sp.nucl == "SEhybrids" & meta$ID %in% rownames(MDS_moms$points)[which(MDS_moms$points[,2] > 5)])] <- "hybridoSE"

meta$col.nucl[which(meta$sp.nucl == "hybridoSE")] <- "#EC008C"

######################################
### Identifying the sex chromosome ###
######################################

# This is a very unelegant way to do it, but it worked.

# Defining a list of the IDs of females ("F") and males ("M") of each (sub)species of interest.
F_rred <- meta$ID[meta$sp.nucl == "rredtenbacheri" & meta$sex=="F"]
M_rred <- meta$ID[meta$sp.nucl == "rredtenbacheri" & meta$sex=="M"]

F_ggra <- meta$ID[meta$sp.nucl == "ggrandii" & meta$sex=="F"]
M_ggra <- meta$ID[meta$sp.nucl == "ggrandii" & meta$sex=="M"]

F_gben <- meta$ID[which(meta$sp.nucl == "gbenazzii" & meta$sex=="F")]
M_gben <- meta$ID[which(meta$sp.nucl == "gbenazzii" & meta$sex=="M")]

# making a subset vcf for eahc species
rred.vcf <- sicily.vcf[,c(1, which(colnames(sicily.vcf@gt) %in% c(M_rred, F_rred)))]
rred.gt <- extract.gt(rred.vcf, element = "GT")

# extracting depth and removing depth info at ungenotyped sites (i.e. because depth < threshold)
rred.dp.all <- extract.gt(rred.vcf, element = "DP", as.numeric = T)
rred.dp <- rred.dp.all
rred.dp[which(is.na(rred.gt))] <- NA

ggra.vcf <- sicily.vcf[,c(1, which(colnames(sicily.vcf@gt) %in% c(M_ggra, F_ggra)))]
ggra.gt <- extract.gt(ggra.vcf, element = "GT")
ggra.dp.all <- extract.gt(ggra.vcf, element = "DP", as.numeric = T)
ggra.dp <- ggra.dp.all
ggra.dp[which(is.na(ggra.gt))] <- NA

gben.vcf <- sicily.vcf[,c(1, which(colnames(sicily.vcf@gt) %in% c(M_gben, F_gben)))]
gben.gt <- extract.gt(gben.vcf, element = "GT")
gben.dp.all <- extract.gt(gben.vcf, element = "DP", as.numeric = T)
gben.dp <- gben.dp.all
gben.dp[which(is.na(gben.gt))] <- NA


### rossius redtenbacheri

# Making a dataframe with one line per locus to store chromosome, plotting colour and male/female depth ratio information
loci_rred <- data.frame(pos.raw = as.numeric(rred.vcf@fix[,2]), chr = substr(rred.vcf@fix[,1], 13, 18), col2="#00000030")
loci_rred$pos <- loci_rred$pos.raw
loci_rred$pos[which(loci_rred$chr=="1")] <- loci_rred$pos.raw[which(loci_rred$chr=="1")] + scf_lengths$start[which(scf_lengths$scf=="1")]
loci_rred$pos[which(loci_rred$chr=="2")] <- loci_rred$pos.raw[which(loci_rred$chr=="2")] + scf_lengths$start[which(scf_lengths$scf=="2")]
loci_rred$pos[which(loci_rred$chr=="3")] <- loci_rred$pos.raw[which(loci_rred$chr=="3")] + scf_lengths$start[which(scf_lengths$scf=="3")]
loci_rred$pos[which(loci_rred$chr=="4_1")] <- loci_rred$pos.raw[which(loci_rred$chr=="4_1")] + scf_lengths$start[which(scf_lengths$scf=="4_1")]
loci_rred$pos[which(loci_rred$chr=="4_2")] <- loci_rred$pos.raw[which(loci_rred$chr=="4_2")] + scf_lengths$start[which(scf_lengths$scf=="4_2")]
loci_rred$pos[which(loci_rred$chr=="5")] <- loci_rred$pos.raw[which(loci_rred$chr=="5")] + scf_lengths$start[which(scf_lengths$scf=="5")]
loci_rred$pos[which(loci_rred$chr=="6")] <- loci_rred$pos.raw[which(loci_rred$chr=="6")] + scf_lengths$start[which(scf_lengths$scf=="6")]
loci_rred$pos[which(loci_rred$chr=="7")] <- loci_rred$pos.raw[which(loci_rred$chr=="7")] + scf_lengths$start[which(scf_lengths$scf=="7")]
loci_rred$pos[which(loci_rred$chr=="8")] <- loci_rred$pos.raw[which(loci_rred$chr=="8")] + scf_lengths$start[which(scf_lengths$scf=="8")]
loci_rred$pos[which(loci_rred$chr=="9_1")] <- loci_rred$pos.raw[which(loci_rred$chr=="9_1")] + scf_lengths$start[which(scf_lengths$scf=="9_1")]
loci_rred$pos[which(loci_rred$chr=="9_2")] <- loci_rred$pos.raw[which(loci_rred$chr=="9_2")] + scf_lengths$start[which(scf_lengths$scf=="9_2")]
loci_rred$pos[which(loci_rred$chr=="10")] <- loci_rred$pos.raw[which(loci_rred$chr=="10")] + scf_lengths$start[which(scf_lengths$scf=="10")]
loci_rred$pos[which(loci_rred$chr=="11")] <- loci_rred$pos.raw[which(loci_rred$chr=="11")] + scf_lengths$start[which(scf_lengths$scf=="11")]
loci_rred$pos[which(loci_rred$chr=="12")] <- loci_rred$pos.raw[which(loci_rred$chr=="12")] + scf_lengths$start[which(scf_lengths$scf=="12")]
loci_rred$pos[which(loci_rred$chr=="13")] <- loci_rred$pos.raw[which(loci_rred$chr=="13")] + scf_lengths$start[which(scf_lengths$scf=="13")]
loci_rred$pos[which(loci_rred$chr=="14")] <- loci_rred$pos.raw[which(loci_rred$chr=="14")] + scf_lengths$start[which(scf_lengths$scf=="14")]
loci_rred$pos[which(loci_rred$chr=="15")] <- loci_rred$pos.raw[which(loci_rred$chr=="15")] + scf_lengths$start[which(scf_lengths$scf=="15")]
loci_rred$pos[which(loci_rred$chr=="16")] <- loci_rred$pos.raw[which(loci_rred$chr=="16")] + scf_lengths$start[which(scf_lengths$scf=="16")]
loci_rred$pos[which(loci_rred$chr=="17")] <- loci_rred$pos.raw[which(loci_rred$chr=="17")] + scf_lengths$start[which(scf_lengths$scf=="17")]
loci_rred$pos[which(loci_rred$chr=="18")] <- loci_rred$pos.raw[which(loci_rred$chr=="18")] + scf_lengths$start[which(scf_lengths$scf=="18")]
loci_rred$pos[which(as.numeric(loci_rred$chr)>18)] <- NA

loci_rred$col2[which(loci_rred$chr %in% c("2", "4_1", "4_2", "6", "8", "10", "12", "14", "16", "18"))] <- "#0599CE30"

# Standardize depth values to remove individual differences in overall depth
rred.rel_dp <- apply(rred.dp, 2, standardize)

# compute per-locus sex-specific average standardized depth
loci_rred$M.meanDP <- apply(rred.rel_dp[,which(dimnames(rred.rel_dp)[2][[1]] %in% M_rred)], 1, mean, na.rm=T)
loci_rred$F.meanDP <- apply(rred.rel_dp[,which(dimnames(rred.rel_dp)[2][[1]] %in% F_rred)], 1, mean, na.rm=T)

# compute male / female depth ratio
loci_rred$DPratio <- loci_rred$M.meanDP / loci_rred$F.meanDP


par(mfrow=c(1,1))

hist(loci_rred$DPratio, breaks=100, xlim=c(0, 3), border=NA, ylim=c(0,10000), main="", xlab="M/F depth ratio", las=1)

#Along the genome
par(mfrow=c(3,1))
plot(loci_rred$M.meanDP ~ loci_rred$pos, col=loci_rred$col2, las=1, xlab="Position along the genome", ylab="Average standardized depth in males", pch=20)
plot(loci_rred$F.meanDP ~ loci_rred$pos, col=loci_rred$col2, las=1, xlab="Position along the genome", ylab="Average standardized depth in females", pch=20)
plot(loci_rred$DPratio ~ loci_rred$pos, col=loci_rred$col2, las=1, xlab="Position along the genome", ylab="Male/female standardized depth ratio", pch=20)
abline(h=1, lty=2)
abline(h=0.5, lty=3)



### grandii grandii ###

loci_ggra <- data.frame(pos.raw = as.numeric(ggra.vcf@fix[,2]), chr = substr(ggra.vcf@fix[,1], 13, 18), col="#00000030")
loci_ggra$pos <- loci$pos.raw
loci_ggra$pos[which(loci_ggra$chr=="1")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="1")] + scf_lengths$start[which(scf_lengths$scf=="1")]
loci_ggra$pos[which(loci_ggra$chr=="2")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="2")] + scf_lengths$start[which(scf_lengths$scf=="2")]
loci_ggra$pos[which(loci_ggra$chr=="3")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="3")] + scf_lengths$start[which(scf_lengths$scf=="3")]
loci_ggra$pos[which(loci_ggra$chr=="4_1")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="4_1")] + scf_lengths$start[which(scf_lengths$scf=="4_1")]
loci_ggra$pos[which(loci_ggra$chr=="4_2")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="4_2")] + scf_lengths$start[which(scf_lengths$scf=="4_2")]
loci_ggra$pos[which(loci_ggra$chr=="5")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="5")] + scf_lengths$start[which(scf_lengths$scf=="5")]
loci_ggra$pos[which(loci_ggra$chr=="6")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="6")] + scf_lengths$start[which(scf_lengths$scf=="6")]
loci_ggra$pos[which(loci_ggra$chr=="7")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="7")] + scf_lengths$start[which(scf_lengths$scf=="7")]
loci_ggra$pos[which(loci_ggra$chr=="8")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="8")] + scf_lengths$start[which(scf_lengths$scf=="8")]
loci_ggra$pos[which(loci_ggra$chr=="9_1")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="9_1")] + scf_lengths$start[which(scf_lengths$scf=="9_1")]
loci_ggra$pos[which(loci_ggra$chr=="9_2")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="9_2")] + scf_lengths$start[which(scf_lengths$scf=="9_2")]
loci_ggra$pos[which(loci_ggra$chr=="10")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="10")] + scf_lengths$start[which(scf_lengths$scf=="10")]
loci_ggra$pos[which(loci_ggra$chr=="11")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="11")] + scf_lengths$start[which(scf_lengths$scf=="11")]
loci_ggra$pos[which(loci_ggra$chr=="12")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="12")] + scf_lengths$start[which(scf_lengths$scf=="12")]
loci_ggra$pos[which(loci_ggra$chr=="13")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="13")] + scf_lengths$start[which(scf_lengths$scf=="13")]
loci_ggra$pos[which(loci_ggra$chr=="14")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="14")] + scf_lengths$start[which(scf_lengths$scf=="14")]
loci_ggra$pos[which(loci_ggra$chr=="15")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="15")] + scf_lengths$start[which(scf_lengths$scf=="15")]
loci_ggra$pos[which(loci_ggra$chr=="16")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="16")] + scf_lengths$start[which(scf_lengths$scf=="16")]
loci_ggra$pos[which(loci_ggra$chr=="17")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="17")] + scf_lengths$start[which(scf_lengths$scf=="17")]
loci_ggra$pos[which(loci_ggra$chr=="18")] <- loci_ggra$pos.raw[which(loci_ggra$chr=="18")] + scf_lengths$start[which(scf_lengths$scf=="18")]
loci_ggra$pos[which(as.numeric(loci_ggra$chr)>18)] <- NA

loci_ggra$col[which(loci_ggra$chr %in% c("2", "4_1", "4_2", "6", "8", "10", "12", "14", "16", "18"))] <- "#BE1E2D30"


# Standardize depth values to remove individual differences in overall depth
ggra.rel_dp <- apply(ggra.dp, 2, standardize)

# compute per-locus sex-specific average standardized depth
loci_ggra$M.meanDP <- apply(ggra.rel_dp[,which(dimnames(ggra.rel_dp)[2][[1]] %in% M_ggra)], 1, mean, na.rm=T)
loci_ggra$F.meanDP <- apply(ggra.rel_dp[,which(dimnames(ggra.rel_dp)[2][[1]] %in% F_ggra)], 1, mean, na.rm=T)

# compute male / female depth ratio
loci_ggra$DPratio <- loci_ggra$M.meanDP / loci_ggra$F.meanDP

par(mfrow=c(1,1))
hist(loci_ggra$DPratio, breaks=100, xlim=c(0, 3), border=NA, ylim=c(0,20000), main="", xlab="M/F depth ratio", las=1)

#Along the genome
par(mfrow=c(3,1))
plot(loci_ggra$M.meanDP ~ loci_ggra$pos, col=loci_ggra$col, las=1, xlab="Position along the genome", ylab="Average standardized depth in males", pch=20, ylim=c(0,8))
plot(loci_ggra$F.meanDP ~ loci_ggra$pos, col=loci_ggra$col, las=1, xlab="Position along the genome", ylab="Average standardized depth in females", pch=20, ylim=c(0,8))
plot(loci_ggra$DPratio ~ loci_ggra$pos, col=loci_ggra$col, las=1, xlab="Position along the genome", ylab="Male/female standardized depth ratio", pch=20)
abline(h=1, lty=2)
abline(h=0.5, lty=3)


### grandii benazzii ###

loci_gben <- data.frame(pos.raw = as.numeric(gben.vcf@fix[,2]), chr = substr(gben.vcf@fix[,1], 13, 18), col="#00000030")
loci_gben$pos <- loci$pos.raw
loci_gben$pos[which(loci_gben$chr=="1")] <- loci_gben$pos.raw[which(loci_gben$chr=="1")] + scf_lengths$start[which(scf_lengths$scf=="1")]
loci_gben$pos[which(loci_gben$chr=="2")] <- loci_gben$pos.raw[which(loci_gben$chr=="2")] + scf_lengths$start[which(scf_lengths$scf=="2")]
loci_gben$pos[which(loci_gben$chr=="3")] <- loci_gben$pos.raw[which(loci_gben$chr=="3")] + scf_lengths$start[which(scf_lengths$scf=="3")]
loci_gben$pos[which(loci_gben$chr=="4_1")] <- loci_gben$pos.raw[which(loci_gben$chr=="4_1")] + scf_lengths$start[which(scf_lengths$scf=="4_1")]
loci_gben$pos[which(loci_gben$chr=="4_2")] <- loci_gben$pos.raw[which(loci_gben$chr=="4_2")] + scf_lengths$start[which(scf_lengths$scf=="4_2")]
loci_gben$pos[which(loci_gben$chr=="5")] <- loci_gben$pos.raw[which(loci_gben$chr=="5")] + scf_lengths$start[which(scf_lengths$scf=="5")]
loci_gben$pos[which(loci_gben$chr=="6")] <- loci_gben$pos.raw[which(loci_gben$chr=="6")] + scf_lengths$start[which(scf_lengths$scf=="6")]
loci_gben$pos[which(loci_gben$chr=="7")] <- loci_gben$pos.raw[which(loci_gben$chr=="7")] + scf_lengths$start[which(scf_lengths$scf=="7")]
loci_gben$pos[which(loci_gben$chr=="8")] <- loci_gben$pos.raw[which(loci_gben$chr=="8")] + scf_lengths$start[which(scf_lengths$scf=="8")]
loci_gben$pos[which(loci_gben$chr=="9_1")] <- loci_gben$pos.raw[which(loci_gben$chr=="9_1")] + scf_lengths$start[which(scf_lengths$scf=="9_1")]
loci_gben$pos[which(loci_gben$chr=="9_2")] <- loci_gben$pos.raw[which(loci_gben$chr=="9_2")] + scf_lengths$start[which(scf_lengths$scf=="9_2")]
loci_gben$pos[which(loci_gben$chr=="10")] <- loci_gben$pos.raw[which(loci_gben$chr=="10")] + scf_lengths$start[which(scf_lengths$scf=="10")]
loci_gben$pos[which(loci_gben$chr=="11")] <- loci_gben$pos.raw[which(loci_gben$chr=="11")] + scf_lengths$start[which(scf_lengths$scf=="11")]
loci_gben$pos[which(loci_gben$chr=="12")] <- loci_gben$pos.raw[which(loci_gben$chr=="12")] + scf_lengths$start[which(scf_lengths$scf=="12")]
loci_gben$pos[which(loci_gben$chr=="13")] <- loci_gben$pos.raw[which(loci_gben$chr=="13")] + scf_lengths$start[which(scf_lengths$scf=="13")]
loci_gben$pos[which(loci_gben$chr=="14")] <- loci_gben$pos.raw[which(loci_gben$chr=="14")] + scf_lengths$start[which(scf_lengths$scf=="14")]
loci_gben$pos[which(loci_gben$chr=="15")] <- loci_gben$pos.raw[which(loci_gben$chr=="15")] + scf_lengths$start[which(scf_lengths$scf=="15")]
loci_gben$pos[which(loci_gben$chr=="16")] <- loci_gben$pos.raw[which(loci_gben$chr=="16")] + scf_lengths$start[which(scf_lengths$scf=="16")]
loci_gben$pos[which(loci_gben$chr=="17")] <- loci_gben$pos.raw[which(loci_gben$chr=="17")] + scf_lengths$start[which(scf_lengths$scf=="17")]
loci_gben$pos[which(loci_gben$chr=="18")] <- loci_gben$pos.raw[which(loci_gben$chr=="18")] + scf_lengths$start[which(scf_lengths$scf=="18")]
loci_gben$pos[which(as.numeric(loci_gben$chr)>18)] <- NA

loci_gben$col[which(loci_gben$chr %in% c("2", "4_1", "4_2", "6", "8", "10", "12", "14", "16", "18"))] <- "#F7941D30"


# Standardize depth values to remove individual differences in overall depth
gben.rel_dp <- apply(gben.dp, 2, standardize)

# compute per-locus sex-specific average standardized depth
loci_gben$M.meanDP <- apply(gben.rel_dp[,which(dimnames(gben.rel_dp)[2][[1]] %in% M_gben)], 1, mean, na.rm=T)
loci_gben$F.meanDP <- apply(gben.rel_dp[,which(dimnames(gben.rel_dp)[2][[1]] %in% F_gben)], 1, mean, na.rm=T)

# compute male / female depth ratio
loci_gben$DPratio <- loci_gben$M.meanDP / loci_gben$F.meanDP

par(mfrow=c(1,1))
hist(loci_gben$DPratio, breaks=120, xlim=c(0, 3), border=NA, ylim=c(0,6000), main="", xlab="M/F depth ratio", las=1)

#Along the genome
par(mfrow=c(3,1))
plot(loci_gben$M.meanDP ~ loci_gben$pos, col=loci_gben$col, las=1, xlab="Position along the genome", ylab="Average standardized depth in males", pch=20, ylim=c(0,8))
plot(loci_gben$F.meanDP ~ loci_gben$pos, col=loci_gben$col, las=1, xlab="Position along the genome", ylab="Average standardized depth in females", pch=20, ylim=c(0,8))
plot(loci_gben$DPratio ~ loci_gben$pos, col=loci_gben$col, las=1, xlab="Position along the genome", ylab="Male/female standardized depth ratio", pch=20)
abline(h=1, lty=2)
abline(h=0.5, lty=3)




loci_rred$posMb <- loci_rred$pos / 1000000
loci_ggra$posMb <- loci_ggra$pos / 1000000
loci_gben$posMb <- loci_gben$pos / 1000000

loci_rred$col2<- "#00000030"
loci_rred$col2[which(loci_rred$chr %in% c("4_1", "4_2", "6", "8", "10", "12", "14", "16", "18"))] <- "#7f7f7f30"
loci_rred$col2[which(loci_rred$chr == "2")] <- "#0599CE30"

loci_ggra$col2<- "#00000030"
loci_ggra$col2[which(loci_ggra$chr %in% c("4_1", "4_2", "6", "8", "10", "12", "14", "16", "18"))] <- "#7f7f7f30"
loci_ggra$col2[which(loci_ggra$chr == "2")] <- "#BE1E2D30"

loci_gben$col2<- "#00000030"
loci_gben$col2[which(loci_gben$chr %in% c("4_1", "4_2", "6", "8", "10", "12", "14", "16", "18"))] <- "#7f7f7f30"
loci_gben$col2[which(loci_gben$chr == "2")] <- "#F7941D30"

### Figure S1 ###

par(mfrow=c(3,1), mai=c(0.6, 0.75, 0.2, 0.02))
plot(loci_rred$DPratio ~ loci_rred$posMb, col=loci_rred$col2, las=1, xlab="", ylab="", pch=20, main="B. rossius redtenbacheri")
plot(loci_ggra$DPratio ~ loci_ggra$posMb, col=loci_ggra$col2, las=1, xlab="", ylab="Male/female standardized depth ratio", pch=20, main="B. grandii grandii")
plot(loci_gben$DPratio ~ loci_gben$posMb, col=loci_gben$col2, las=1, xlab="Position along the genome [Mbases]", ylab="", pch=20, main="B. grandii benazzii")

##############################
### Sexing all individuals ###
##############################

sicily_auto_std.dp <- sicily_std.dp[which(sicily.vcf@fix[,1] != "Brsri_v3_scf2"),]

sicily_X_std.dp <- sicily_std.dp[which(sicily.vcf@fix[,1] == "Brsri_v3_scf2"),]

meta$avg.auto.dp <- apply(sicily_auto_std.dp, 2, mean, na.rm=T)
meta$avg.X.dp <- apply(sicily_X_std.dp, 2, mean, na.rm=T)
meta$X_auto_ratio <- meta$avg.X.dp / meta$avg.auto.dp

par(mfrow=c(1,1))
hist(meta$X_auto_ratio[meta$sex=="M"], las=1, xlab="Depth ratio (X / autosomes)", main="", xlim=c(0.5, 1.3), ylim=c(0,80), breaks=10, border=F, col="#1C75BC")
par(new=T)
hist(meta$X_auto_ratio[meta$sex=="F"], las=1, xlab="", main="", xlim=c(0.5, 1.3), ylim=c(0,80), ylab="", breaks=20, border=F, col="#F79A70")
par(new=T)
hist(meta$X_auto_ratio[is.na(meta$sex)], las=1, xlab="", main="", xlim=c(0.5, 1.3), ylim=c(0,80), ylab="", breaks=40, border=F, col="#000000")

max(meta$X_auto_ratio[meta$sex=="M"], na.rm=T) # 0.7254365
min(meta$X_auto_ratio[meta$sex=="F"], na.rm=T) # 0.9068253

meta$sex.genetic <- NA
meta$sex.genetic[meta$X_auto_ratio <= max(meta$X_auto_ratio[meta$sex=="M"], na.rm=T)] <- "M"
meta$sex.genetic[meta$X_auto_ratio >= min(meta$X_auto_ratio[meta$sex=="F"], na.rm=T)] <- "F"



#####################
### Figure 2A & B ###
#####################

# We repeat the analysis but exclude the sex chromosome (scaffold 2)
sicily_auto.vcf <- sicily.vcf[which(sicily.vcf@fix[,1] != "Brsri_v3_scf2"),]
sicily_auto.gt <- extract.gt(sicily_auto.vcf, element = "GT", convertNA = T)
sicily_auto.gl <- vcfR2genlight(sicily_auto.vcf) # excludes 5244 sites because they are triallelic

dist_sicily_auto <- dist(sicily_auto.gl)
MDS_sicily_auto <- isoMDS(d = dist_sicily_auto, k = 10, maxit = 100)

#This is Figure 2A and 2B
par(mfrow=c(1,2))
plot(MDS_sicily_auto$points[,2] ~ MDS_sicily_auto$points[,1], las=1, col=meta$col.mt, pch=20, xlab="Dimension 1", ylab="Dimension 2", main="MDS - autosomes")
plot(MDS_sicily_auto$points[,3] ~ MDS_sicily_auto$points[,1], las=1, col=meta$col.mt, pch=20, xlab="Dimension 1", ylab="Dimension 3", main="MDS - autosomes")

#################
### admixture ###
#################

# Selecting loci that are genotyped in all subspecies
# We retain loci if they are genotyped in at least 3 individuals from each subspecies:
# rred, gben, ggra, att. not in gmar because there are too few individuals

N_genotyped_in_rred <- apply(sicily.gt[,which(dimnames(sicily.gt)[[2]] %in% meta$ID[which(meta$sp.nucl1 == "rredtenbacheri")])], 1, count.not.na)
N_genotyped_in_ggra <- apply(sicily.gt[,which(dimnames(sicily.gt)[[2]] %in% meta$ID[which(meta$sp.nucl1 == "ggrandii")])], 1, count.not.na)
N_genotyped_in_gben <- apply(sicily.gt[,which(dimnames(sicily.gt)[[2]] %in% meta$ID[which(meta$sp.nucl1 == "gbenazzii")])], 1, count.not.na)
N_genotyped_in_gmar <- apply(sicily.gt[,which(dimnames(sicily.gt)[[2]] %in% meta$ID[which(meta$sp.nucl1 == "gmaretimi")])], 1, count.not.na)
N_genotyped_in_att <- apply(sicily.gt[,which(dimnames(sicily.gt)[[2]] %in% meta$ID[which(meta$sp.nucl1 == "atticus")])], 1, count.not.na)

sum(N_genotyped_in_rred > 2 & N_genotyped_in_ggra > 2 & N_genotyped_in_gben > 2 & N_genotyped_in_att > 2 )

# Making a list of positions genotyped in at least 3 individuals for each (sub)species for Figure 2E
genotyped_in_all <- sicily.vcf@fix[which(N_genotyped_in_rred > 2 & N_genotyped_in_ggra > 2 & N_genotyped_in_gben > 2 & N_genotyped_in_att > 2),1:2]
write.table(genotyped_in_all, "positions_genotyped_in_all_sspp.txt", quote = F, col.names = F, row.names = F)

### See README.txt for admixture command

# Making a list of positions genotyped in at least 3 individuals for all parental (sub)species of each hybrid lineage for Figure 3
genotyped_in_SE <- sicily.vcf@fix[which(N_genotyped_in_rred > 2 & N_genotyped_in_ggra > 2),1:2]
write.table(genotyped_in_SE, "positions_genotyped_in_SE.txt", quote = F, col.names = F, row.names = F)

genotyped_in_NW <- sicily.vcf@fix[which(N_genotyped_in_rred > 2 & N_genotyped_in_gben > 2),1:2]
write.table(genotyped_in_NW, "positions_genotyped_in_NW.txt", quote = F, col.names = F, row.names = F)

genotyped_in_lynceorums_ancestors <- sicily.vcf@fix[which(N_genotyped_in_rred > 2 & N_genotyped_in_ggra > 2 & N_genotyped_in_att > 2),1:2]
write.table(genotyped_in_lynceorums_ancestors, "positions_genotyped_in_lynceorum_ancestors.txt", quote = F, col.names = F, row.names = F)


# Summarizing Cross-validation error across runs
CV <- read.table("CV_summary_across_rep.txt", h=F)
CV <- CV[order(CV$V1),]
dimnames(CV)[[1]] <- 1:12

par(mfrow=c(1,1))plot(NULL, xlim=range(CV$V1), ylim=range(CV[,2:dim(CV)[2]]), las=1, ylab="Cross-Validation error", xlab = "K")
for (i in 2:dim(CV)[2]) {
  points(CV[,i] ~ CV$V1, type="b", pch=20, col=rgb(0, 0, 0, 0.2))
}

CV.mean <- apply(CV[,-1], 1, mean)
CV.sd <- apply(CV[,-1], 1, sd)
plot(CV.mean, pch=20, cex=0.75, col=NULL, las=1)
for (i in 1:20) segments(x0=i, x1=i, y0=CV.mean[i] - 1.96*CV.sd[i], y1=CV.mean[i] + 1.96*CV.sd[i], lwd=2)
which(CV.mean == min(CV.mean))

# the best replicate was replicate 10
k5_auto <- read.table("admixture_auto_Sicily/rep10/cleaned.5.Q")
dimnames(k5_auto)[[1]] <- meta$ID

### Estimating heterozygosity
meta$het <- apply(sicily.gt, 2, prop.het)


### Figure 2D, 2E and 2F
par(mfrow=c(3,1), mar=c(0, 2, 0.1, 0.1))
barplot(meta$het[match(tree$tip.label, meta$ID)], las=2, border=NA, cex.names = 0)
barplot(t(as.matrix(k5_auto[match(tree$tip.label, meta$ID),])), las=2, border=NA, cex.names = 0, space=0, col=c("#F7941D", "#0599CE", "#F9ED32", "#71BF44", "#BE1E2D"))
plot.phylo(tree, direction="upwards", align.tip.label = T, cex=0.2)
add.scale.bar()
# Individual with identical mtDNA haplotypes were then manually swapped in Illustrator
# to group similar nuclear genotypes together for display purposes.

################################################
### Per-chromosome admixture (Figure 3 & S4) ###
################################################

order_SE <- c(which(read.table("per_chromosome_admixture_SE/chr1/rep1/cleaned.2.Q")$V2 < 0.01), which(read.table("per_chromosome_admixture_SE/chr1/rep1/cleaned.2.Q")$V2 > 0.1 & read.table("per_chromosome_admixture_SE/chr1/rep1/cleaned.2.Q")$V2 < 0.9 & read.table("per_chromosome_admixture_SE/chr2/rep1/cleaned.2.Q")$V2 < 0.1), which(read.table("per_chromosome_admixture_SE/chr1/rep1/cleaned.2.Q")$V2 > 0.1 & read.table("per_chromosome_admixture_SE/chr1/rep1/cleaned.2.Q")$V2 < 0.9 & read.table("per_chromosome_admixture_SE/chr2/rep1/cleaned.2.Q")$V2 > 0.1), which(read.table("per_chromosome_admixture_SE/chr1/rep1/cleaned.2.Q")$V2 > 0.9))
par(mfrow=c(18,1), mar=c(0, 0.1, 0.05, 0.1))
for (i in 1:18){
  x <- read.table(file = paste("per_chromosome_admixture_SE/chr", i, "/rep7/cleaned.2.Q", sep=""), h=F)[,c(2,1)]
  barplot(t(as.matrix(x[order_SE,])), las=2, border=NA, cex.names = 0.1, space=0, col=c("#BE1E2D", "#0599CE"))
}

order_NW <- c(which(read.table("per_chromosome_admixture_NW/chr1/rep1/cleaned.2.Q")$V2 > 0.9), which(read.table("per_chromosome_admixture_NW/chr1/rep1/cleaned.2.Q")$V2 > 0.1 & read.table("per_chromosome_admixture_NW/chr1/rep1/cleaned.2.Q")$V2 < 0.9 & read.table("per_chromosome_admixture_NW/chr2/rep1/cleaned.2.Q")$V1 < 0.1), which(read.table("per_chromosome_admixture_NW/chr1/rep1/cleaned.2.Q")$V2 > 0.1 & read.table("per_chromosome_admixture_NW/chr1/rep1/cleaned.2.Q")$V2 < 0.9 & read.table("per_chromosome_admixture_NW/chr2/rep1/cleaned.2.Q")$V1 > 0.1), which(read.table("per_chromosome_admixture_NW/chr1/rep1/cleaned.2.Q")$V2 < 0.01))
par(mfrow=c(18,1), mar=c(0, 0.1, 0.05, 0.1))
for (i in 1:18){
  x <- read.table(file = paste("per_chromosome_admixture_NW/chr", i, "/rep1/cleaned.2.Q", sep=""), h=F)
  barplot(t(as.matrix(x[order_NW,])), las=2, border=NA, cex.names = 0.1, space=0, col=c("#F7941D", "#0599CE"))
}

#order_lynceorum <- c(which(read.table("per_chromosome_admixture_lynceorum/chr1/rep1/cleaned.3.Q")$V1 > 0.9), which(read.table("per_chromosome_admixture_lynceorum/chr1/rep1/cleaned.3.Q")$V1 > 0.1 & read.table("per_chromosome_admixture_lynceorum/chr1/rep1/cleaned.3.Q")$V2 > 0.1 & read.table("per_chromosome_admixture_lynceorum/chr2/rep1/cleaned.3.Q")$V3 > 0.1), which(read.table("per_chromosome_admixture_lynceorum/chr1/rep1/cleaned.3.Q")$V2 > 0.9), which(read.table("per_chromosome_admixture_lynceorum/chr1/rep1/cleaned.3.Q")$V3 > 0.9))
order_lynceorum <- c(which(read.table("per_chromosome_admixture_lynceorum/chr2/rep1/cleaned.3.Q")$V1 > 0.2 & read.table("per_chromosome_admixture_lynceorum/chr2/rep1/cleaned.3.Q")$V2 > 0.2 & read.table("per_chromosome_admixture_lynceorum/chr2/rep1/cleaned.3.Q")$V3 > 0.2), which(read.table("per_chromosome_admixture_lynceorum/chr2/rep1/cleaned.3.Q")$V2 > 0.4 & read.table("per_chromosome_admixture_lynceorum/chr2/rep1/cleaned.3.Q")$V2 < 0.6))
par(mfrow=c(18,1), mar=c(0, 0.1, 0.05, 0.1))
for (i in 1:18){
  x <- read.table(file = paste("per_chromosome_admixture_lynceorum/chr", i, "/rep1/cleaned.3.Q", sep=""), h=F)
  barplot(t(as.matrix(x[order_lynceorum,])), las=2, border=NA, cex.names = 0.1, space=0, col=c("#0599CE", "#71BF44", "#BE1E2D"))
}


### Figure 3: X chromosome specific
par(mfrow=c(1,1))
admX_NW <- read.table("per_chromosome_admixture_NW/chr2/rep7/cleaned.2.Q", h=F)[,c(2,1)]
dimnames(admX_NW)[[1]] <- NW_ind_list
barplot(t(as.matrix(admX_NW[order_NW,])), las=2, border=NA, cex.names = 0.2, space=0, col=c("#F7941D", "#0599CE"))

admX_SE <- read.table("per_chromosome_admixture_SE/chr2/rep7/cleaned.2.Q", h=F)[,c(2,1)]
dimnames(admX_SE)[[1]] <- SE_ind_list
barplot(t(as.matrix(admX_SE[order_SE,])), las=2, border=NA, cex.names = 0.2, space=0, col=c("#BE1E2D", "#0599CE"))

admX_lynceorum <- read.table("per_chromosome_admixture_lynceorum/chr2/rep1/cleaned.3.Q", h=F)
dimnames(admX_lynceorum)[[1]] <- lynceorum_ind_list
barplot(t(as.matrix(admX_lynceorum[order_lynceorum,])), las=2, border=NA, cex.names = 1, space=0, col=c("#0599CE", "#71BF44", "#BE1E2D"))

##############################################
### Evolution of heterozygosity in hybrids ###
##############################################

# Reading filtered datasets with loci genotyped in both parental species.
SE.vcf <- read.vcfR("SE_shared_gatk_an2_ql20_snps_dp08_mc3_miss50.recode.vcf")
SE.gt <- extract.gt(SE.vcf, element = "GT", convertNA = T)
SE_auto.gt <- extract.gt(SE.vcf[which(SE.vcf@fix[,1] != "Brsri_v3_scf2"),], element = "GT", convertNA = T)
SE_ind_list <- dimnames(SE.gt)[[2]]

NW.vcf <- read.vcfR("NW_shared_gatk_an2_ql20_snps_dp08_mc3_miss50.recode.vcf")
NW.gt <- extract.gt(NW.vcf, element = "GT", convertNA = T)
NW_auto.gt <- extract.gt(NW.vcf[which(NW.vcf@fix[,1] != "Brsri_v3_scf2"),], element = "GT", convertNA = T)
NW_ind_list <- dimnames(NW.gt)[[2]]

# Identifying column indexes of individuals of each lineage
ggra_in_SE <- which(colnames(SE.gt) %in% meta$ID[which(meta$sp.nucl1 == "ggrandii")])
rred_in_SE <- which(colnames(SE.gt) %in% meta$ID[which(meta$sp.nucl1 == "rredtenbacheri")])
SEhybrids_in_SE <- which(colnames(SE.gt) %in% meta$ID[which(meta$sp.nucl1 == "SEhybrids")])

whitei_in_SE <- which(colnames(SE.gt) %in% meta$ID[which(meta$sp.nucl == "whitei")])
hybridoSE_in_SE <- which(colnames(SE.gt) %in% meta$ID[which(meta$sp.nucl == "hybridoSE")])

gben_in_NW <- which(colnames(NW.gt) %in% meta$ID[which(meta$sp.nucl1 == "gbenazzii")])
rred_in_NW <- which(colnames(NW.gt) %in% meta$ID[which(meta$sp.nucl1 == "rredtenbacheri")])
hybridoNW_in_NW <- which(colnames(NW.gt) %in% meta$ID[which(meta$sp.nucl1 == "hybridoNW")])



# This function returns 1 if two groups of individuals (provided) are fixed for opposite alleles;
# it returns 0 otherwise.
# gt.line is one line of a genotype file (obtained with extract.gt(vcf, "GT"))
# one is a vector containing the positions of the columns of a species of interest
# two is a vector containing the positions of the columns of the species for comparison
# fixed.opposite() returns 1 if alternative alleles are fixed in ones and in two, 0 otherwise.

fixed.opposite <- function(gt.line, one, two) {
  fixed.opposite <- 0
  if (sum(!is.na(gt.line[one]))>0 & sum(!is.na(gt.line[two]))>0) {
    ones <- unique(gt.line[one][which(!is.na(gt.line[one]))])
    if (sum(ones %in% c("0/0", "1/1")) == 1) {
      if (sum(!is.na(gt.line[two]))>0) {
        twos <- unique(gt.line[two][which(!is.na(gt.line[two]))])
        if (length(unique(twos))==1 & length(unique(ones))==1 ) {
          if (sum(twos %in% c("0/0", "1/1") == 1) & ones != twos) fixed.opposite <- 1
        }
      }
    }
  }
  return(fixed.opposite)
}

# Checking how many loci are fixed for opposite alleles in species pairs in each dataset
# And as a proportion to have a minimum proportion of possible heterozygosity
sum(apply(SE_auto.gt, 1, fixed.opposite, ggra_in_SE, rred_in_SE))
sum(apply(SE_auto.gt, 1, fixed.opposite, ggra_in_SE, rred_in_SE)) / dim(SE.gt)[[1]]
sum(apply(NW_auto.gt, 1, fixed.opposite, gben_in_NW, rred_in_NW))
sum(apply(NW_auto.gt, 1, fixed.opposite, gben_in_NW, rred_in_NW)) / dim(NW.gt)[[1]]


# This function computes allele frequencies for one locus.
# gt.line is one line of a genotype file (obtained with extract.gt(vcf, "GT"))
# This mimics the structure of the .frq outputs of vcftools (option --freq)
# except "POS" that I don't care about.
# It was more convenient to do it here than import everything from vcftools
get.allele.frq <- function(gt.line) {
  N_CHR <- (2* sum(is.na(gt.line)==F))
  FRQ_REF <- ((2*sum(gt.line == "0/0", na.rm = T) + sum(gt.line == "0/1", na.rm = T)) / N_CHR)
  return(FRQ_REF)
}

get.N.CHR <- function(gt.line) {
  N_CHR <- (2* sum(is.na(gt.line)==F))
  return(N_CHR)
}

# Computing allele frequencies in all lineages for each locus
SE.frq <- data.frame(N_CHR = apply(SE_auto.gt, 1, get.N.CHR), FRQ_REF = apply(SE_auto.gt, 1, get.allele.frq), FRQ_ALT = (1-apply(SE_auto.gt, 1, get.allele.frq)))
ggra_in_SE.frq <- data.frame(N_CHR = apply(SE_auto.gt[,ggra_in_SE], 1, get.N.CHR), FRQ_REF = apply(SE_auto.gt[,ggra_in_SE], 1, get.allele.frq), FRQ_ALT = (1-apply(SE_auto.gt[,ggra_in_SE], 1, get.allele.frq)))
rred_in_SE.frq <- data.frame(N_CHR = apply(SE_auto.gt[,rred_in_SE], 1, get.N.CHR), FRQ_REF = apply(SE_auto.gt[,rred_in_SE], 1, get.allele.frq), FRQ_ALT = (1-apply(SE_auto.gt[,rred_in_SE], 1, get.allele.frq)))

NW.frq <- data.frame(N_CHR = apply(NW_auto.gt, 1, get.N.CHR), FRQ_REF = apply(NW_auto.gt, 1, get.allele.frq), FRQ_ALT = (1-apply(NW_auto.gt, 1, get.allele.frq)))
gben_in_NW.frq <- data.frame(N_CHR = apply(NW_auto.gt[,gben_in_NW], 1, get.N.CHR), FRQ_REF = apply(NW_auto.gt[,gben_in_NW], 1, get.allele.frq), FRQ_ALT = (1-apply(NW_auto.gt[,gben_in_NW], 1, get.allele.frq)))
rred_in_NW.frq <- data.frame(N_CHR = apply(NW_auto.gt[,rred_in_NW], 1, get.N.CHR), FRQ_REF = apply(NW_auto.gt[,rred_in_NW], 1, get.allele.frq), FRQ_ALT = (1-apply(NW_auto.gt[,rred_in_NW], 1, get.allele.frq)))


# This function computes the number of loci that are genotyped 
# in at least two individuals from one groupe and none of another
missing.in.one <- function(gt.line, one, two) {
  missing.here <- 0
  if (sum(!is.na(gt.line[one]))==0) {
    if (sum(!is.na(gt.line[two]))>=2) missing.here <- 1
  }
  return(missing.here)
}

sum(apply(SE_auto.gt, 1, missing.in.one, ggra_in_SE, rred_in_SE))
sum(apply(NW_auto.gt, 1, missing.in.one, gben_in_NW, rred_in_NW))



# Simulated expected heterozygosity in de novo hybrids between the ancestral species
# For this I create genotypes by sampling alleles in the allele frequencies of the parental species
# and I look at the resulting relative heterozygosity
simulate.heterozygosity <- function(frq.A, frq.B, nrep) {
  simulated.het <- numeric(0)
  for(i in 1:nrep) {
    het.per.locus <- logical(0)
    for (l in (1:length(frq.A$N_CHR))[typed.loci <- which(frq.A$N_CHR > 0 & frq.B$N_CHR > 0)]) {
      genoA <- sample(x=c(0,1), size=1, prob=c(frq.A$FRQ_REF[l], frq.A$FRQ_ALT[l]))
      genoB <- sample(x=c(0,1), size=1, prob=c(frq.B$FRQ_REF[l], frq.B$FRQ_ALT[l]))
      het.per.locus <- c(het.per.locus, (genoA + genoB) == 1)
    }
    simulated.het <- c(simulated.het, sum(het.per.locus)/length(typed.loci))
  }
  return(simulated.het)
}

simulated.heterozygosity.rred.ggra <- simulate.heterozygosity(rred_in_SE.frq, ggra_in_SE.frq, 1000)
simulated.heterozygosity.rred.gben <- simulate.heterozygosity(rred_in_NW.frq, gben_in_NW.frq, 1000)

### Observed heterozygosity with the same datasets
### as used for the simulated values
SE.het <- apply(SE_auto.gt, 2, prop.het)
NW.het <- apply(NW_auto.gt, 2, prop.het)


### Figure 4A: Observed vs expected heterozygosity values

par(mfrow=(c(3,1)), mai=c(0.6,0.5,0.05,0.05))
boxplot(SE.het[whitei_in_SE], horizontal = T, ylim=c(0.69, 0.86), border = "#92278f", outline = F, col="#92278f77")
points(jitter(rep(1, length(whitei_in_SE)), factor=4) ~ SE.het[whitei_in_SE], pch=16)

abline(v=min(simulated.heterozygosity.rred.ggra))
abline(v=max(simulated.heterozygosity.rred.ggra))

boxplot(SE.het[hybridoSE_in_SE], horizontal = T, ylim=c(0.69, 0.86), border = "#EC008C", outline = F, col="#EC008C77")
points(jitter(rep(1, length(hybridoSE_in_SE)), factor=4) ~ SE.het[hybridoSE_in_SE], pch=16)

abline(v=min(simulated.heterozygosity.rred.ggra))
abline(v=max(simulated.heterozygosity.rred.ggra))
#hist(simulated.heterozygosity.rred.ggra, xlim=c(0.69, 0.86), las=1, xlab="", ylab="Frequency", border=NA, col="grey")

### Figure 4B: same with NW
boxplot(NW.het[hybridoNW_in_NW], horizontal = T, ylim=c(0.69, 0.86), border = "#f47f72", outline = F, col="#f47f7277")
points(jitter(rep(1, length(hybridoNW_in_NW)), factor=4) ~ NW.het[hybridoNW_in_NW], pch=16)

abline(v=min(simulated.heterozygosity.rred.gben))
abline(v=max(simulated.heterozygosity.rred.gben))


#############################################################
### Where do hybrids loose heterozygosity ? -- per window ###
#############################################################

# The same ugly home-made Manhattan plot method as for the sex chromosome.
chr1 <- data.frame(begin.pos = seq(from=1, to=scf_lengths$end[which(scf_lengths$scf=="1")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="1")], by=500000))
chr1$chr <- "1"
chr2 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="2")], to=scf_lengths$end[which(scf_lengths$scf=="2")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="2")], by=500000))
chr2$chr <- "2"
chr3 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="3")], to=scf_lengths$end[which(scf_lengths$scf=="3")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="3")], by=500000))
chr3$chr <- "3"
chr4_1 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="4_1")], to=scf_lengths$end[which(scf_lengths$scf=="4_1")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="4_1")], by=500000))
chr4_1$chr <- "4_1"
chr4_2 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="4_2")], to=scf_lengths$end[which(scf_lengths$scf=="4_2")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="4_2")], by=500000))
chr4_2$chr <- "4_2"
chr5 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="5")], to=scf_lengths$end[which(scf_lengths$scf=="5")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="5")], by=500000))
chr5$chr <- "5"
chr6 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="6")], to=scf_lengths$end[which(scf_lengths$scf=="6")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="6")], by=500000))
chr6$chr <- "6"
chr7 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="7")], to=scf_lengths$end[which(scf_lengths$scf=="7")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="7")], by=500000))
chr7$chr <- "7"
chr8 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="8")], to=scf_lengths$end[which(scf_lengths$scf=="8")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="8")], by=500000))
chr8$chr <- "8"
chr9_1 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="9_1")], to=scf_lengths$end[which(scf_lengths$scf=="9_1")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="9_1")], by=500000))
chr9_1$chr <- "9_1"
chr9_2 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="9_2")], to=scf_lengths$end[which(scf_lengths$scf=="9_2")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="9_2")], by=500000))
chr9_2$chr <- "9_2"
chr10 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="10")], to=scf_lengths$end[which(scf_lengths$scf=="10")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="10")], by=500000))
chr10$chr <- "10"
chr11 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="11")], to=scf_lengths$end[which(scf_lengths$scf=="11")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="11")], by=500000))
chr11$chr <- "11"
chr12 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="12")], to=scf_lengths$end[which(scf_lengths$scf=="12")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="12")], by=500000))
chr12$chr <- "12"
chr13 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="13")], to=scf_lengths$end[which(scf_lengths$scf=="13")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="13")], by=500000))
chr13$chr <- "13"
chr14 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="14")], to=scf_lengths$end[which(scf_lengths$scf=="14")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="14")], by=500000))
chr14$chr <- "14"
chr15 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="15")], to=scf_lengths$end[which(scf_lengths$scf=="15")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="15")], by=500000))
chr15$chr <- "15"
chr16 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="16")], to=scf_lengths$end[which(scf_lengths$scf=="16")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="16")], by=500000))
chr16$chr <- "16"
chr17 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="17")], to=scf_lengths$end[which(scf_lengths$scf=="17")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="17")], by=500000))
chr17$chr <- "17"
chr18 <- data.frame(begin.pos = seq(from=scf_lengths$start[which(scf_lengths$scf=="18")], to=scf_lengths$end[which(scf_lengths$scf=="18")], by=500000), begin = seq(from=1, to=scf_lengths$length[which(scf_lengths$scf=="18")], by=500000))
chr18$chr <- "18"

windows.1M <- rbind(chr1, chr2, chr3, chr4_1, chr4_2, chr5, chr6, chr7, chr8, chr9_1, chr9_2, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18)



### South-Eastern hybrids ###

het_per_window_SE <- matrix(data = NA, nrow = length(windows.1M$begin), ncol = length(SEhybrids_in_SE))
Nloci_per_window_SE <- matrix(data = NA, nrow = length(windows.1M$begin), ncol = length(SEhybrids_in_SE))

# Iterate over all windows. Compute heterozygosity 
for (i in 1:length(windows.1M$begin.pos)) {
  loci <- which(as.numeric(SE.vcf@fix[,2]) >= windows.1M$begin[i] & as.numeric(SE.vcf@fix[,2]) <= (windows.1M$begin[i] + 999999) & substr(SE.vcf@fix[,1], 13, 15) == windows.1M$chr[i])
  if (length(loci) == 1) {
    het_per_window_SE[i,] <- prop.het(SE.gt[loci,SEhybrids_in_SE])
    Nloci_per_window_SE[i,] <- count.not.na(SE.gt[loci,SEhybrids_in_SE])
  } else if (length(loci) > 1) {
    het_per_window_SE[i,] <- apply(SE.gt[loci,SEhybrids_in_SE], 2, prop.het)
    Nloci_per_window_SE[i,] <- apply(SE.gt[loci,SEhybrids_in_SE], 2, count.not.na)
  }
}


# Simulated expected heterozygosity in de novo hybrids between the ancestral species per locus
# For this I create genotypes by sampling alleles in the allele frequencies of the parental species
# and I look at the resulting relative heterozygosity
simulate.heterozygosity.individual <- function(frq.A, frq.B, nrep) {
  simulated.het <- matrix(data = NA, nrow = length(frq.A$N_CHR), ncol = nrep)
  het.per.locus <- logical(0)
  for (l in (1:length(frq.A$N_CHR))[typed.loci <- which(frq.A$N_CHR > 0 & frq.B$N_CHR > 0)]) {
    genoA <- sample(x=c(0,1), size=nrep, prob=c(frq.A$FRQ_REF[l], frq.A$FRQ_ALT[l]), replace = T)
    genoB <- sample(x=c(0,1), size=nrep, prob=c(frq.B$FRQ_REF[l], frq.B$FRQ_ALT[l]), replace = T)
    simulated.het[l,] <- (genoA + genoB) == 1
  }
  return(simulated.het)
}

# Simulate heterozygosity in hybrids for 1000 individuals (i.e. replicates) at each locus,
# based on allelic frequencies in the parental populations.
simulated.heterozygosity.individual.rred.ggra <- simulate.heterozygosity.individual(rred_in_SE.frq, ggra_in_SE.frq, 1000)
# average across individuals.
average.simulated.heterozygosity.rred.ggra <- apply(simulated.heterozygosity.individual.rred.ggra, 1, sum, na.rm=T)/dim(simulated.heterozygosity.individual.rred.ggra)[2]


# average heterozygosity (already averaged across replicates) per window
sim_het_per_window_SE <- numeric(length = length(average.simulated.heterozygosity.rred.ggra))
for (i in 1:length(windows.1M$begin.pos)) {
  loci <- which(as.numeric(SE.vcf@fix[,2]) >= windows.1M$begin[i] & as.numeric(SE.vcf@fix[,2]) <= (windows.1M$begin[i] + 999999) & substr(SE.vcf@fix[,1], 13, 15) == windows.1M$chr[i])
  sim_het_per_window_SE[i] <- mean(average.simulated.heterozygosity.rred.ggra[loci])
}

### North-Western hybrids ###

het_per_window_NW <- matrix(data = NA, nrow = length(windows.1M$begin), ncol = length(hybridoNW_in_NW))
Nloci_per_window_NW <- matrix(data = NA, nrow = length(windows.1M$begin), ncol = length(hybridoNW_in_NW))

# I'm estimating average heterozygosity across all SNPs in each window for each individual
# I'm also recording how many SNPs were genotyped per window in each individual
for (i in 1:length(windows.1M$begin.pos)) {
  loci <- which(as.numeric(NW.vcf@fix[,2]) >= windows.1M$begin[i] & as.numeric(NW.vcf@fix[,2]) <= (windows.1M$begin[i] + 999999) & substr(NW.vcf@fix[,1], 13, 15) == windows.1M$chr[i])
  if (length(loci) == 1) {
    het_per_window_NW[i,] <- prop.het(NW.gt[loci,hybridoNW_in_NW])
    Nloci_per_window_NW[i,] <- count.not.na(NW.gt[loci,hybridoNW_in_NW])
  } else if (length(loci) > 1) {
    het_per_window_NW[i,] <- apply(NW.gt[loci,hybridoNW_in_NW], 2, prop.het)
    Nloci_per_window_NW[i,] <- apply(NW.gt[loci,hybridoNW_in_NW], 2, count.not.na)
  }
}


simulated.heterozygosity.individual.rred.gben <- simulate.heterozygosity.individual(rred_in_NW.frq, gben_in_NW.frq, 1000)
average.simulated.heterozygosity.rred.gben <- apply(simulated.heterozygosity.individual.rred.gben, 1, sum, na.rm=T)/dim(simulated.heterozygosity.individual.rred.gben)[2]

sim_het_per_window_NW <- numeric(length = length(average.simulated.heterozygosity.rred.gben))

for (i in 1:length(windows.1M$begin.pos)) {
  loci <- which(as.numeric(NW.vcf@fix[,2]) >= windows.1M$begin[i] & as.numeric(NW.vcf@fix[,2]) <= (windows.1M$begin[i] + 999999) & substr(NW.vcf@fix[,1], 13, 15) == windows.1M$chr[i])
  sim_het_per_window_NW[i] <- mean(average.simulated.heterozygosity.rred.gben[loci])
}


par(mfrow=c(5,1), mai=c(0.5, 0.5, 0.25, 0.1))
for (scf in 1:20) {
  plot(NULL, xlim=c(0, scf_lengths$length[scf]), ylim=c(0,1), las=1, xlab="", ylab="", main = scf_lengths$scf[scf])
  for (i in 1:length(hybridoNW_in_NW)) {
    points(het_per_window_NW[which(windows.1M$chr==scf_lengths$scf[scf]),i] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])], col=rgb(0,0,0,0.1), type="l")
  }
  points(sim_het_per_window_NW[which(windows.1M$chr==scf_lengths$scf[scf])] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])], col="#F47F72", type="l", lwd=2)
}



### Figure S6: Expected and observed heterozygosity along all chromosomes

for (scf in 1:20) {
  svg(filename = paste("../Figures/heterozygosity_along_genome_whitei_hybridoSE_NW_chr", scf_lengths$scf[scf], ".svg", sep=""), height = 20, width = scf_lengths$length[scf] / 5000000)
  par(mfrow=c(5,1), mai=c(0.5, 0.8, 0.2, 0.1))
  barplot(apply(Nloci_per_window_SE[which(windows.1M$chr==scf_lengths$scf[scf]),], 1, mean, na.rm=T), border=NA, ylim=c(0,105), las=1)
  plot(NULL, xlim=c(0, scf_lengths$length[scf]), ylim=c(0,1), las=1, xlab="", ylab="", main = "", xaxt="n")
  for (i in which(colnames(SE.gt)[SEhybrids_in_SE] %in% meta$ID[which(meta$sp.nucl == "whitei")])) {
    points(het_per_window_SE[which(windows.1M$chr==scf_lengths$scf[scf]), i] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])], col="#92278F30", type="l")
  }
  points(sim_het_per_window_SE[which(windows.1M$chr==scf_lengths$scf[scf])] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])], col=1, type="l", lwd=2)
  
  plot(NULL, xlim=c(0, scf_lengths$length[scf]), ylim=c(0,1), las=1, xlab="", ylab="", main = "", xaxt="n")
  for (i in which(colnames(SE.gt)[SEhybrids_in_SE] %in% meta$ID[which(meta$sp.nucl == "hybridoSE")])) {
    points(het_per_window_SE[which(windows.1M$chr==scf_lengths$scf[scf]), i] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])], col="#EC008C30", type="l")
  }
  points(sim_het_per_window_SE[which(windows.1M$chr==scf_lengths$scf[scf])] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])], col=1, type="l", lwd=2)
  
  plot(NULL, xlim=c(0, scf_lengths$length[scf]), ylim=c(0,1), las=1, xlab="", ylab="", main = "", xaxt="n")
  for (i in 1:length(hybridoNW_in_NW)) {
    points(het_per_window_NW[which(windows.1M$chr==scf_lengths$scf[scf]),i] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])], col="#F47F7230", type="l")
  }
  points(sim_het_per_window_NW[which(windows.1M$chr==scf_lengths$scf[scf])] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])], col=1, type="l", lwd=2)
  barplot(apply(Nloci_per_window_NW[which(windows.1M$chr==scf_lengths$scf[scf]),], 1, mean, na.rm=T), border=NA, ylim=c(0,105), las=1)
  dev.off()
}


### Figure 3C: first 50 Mb (= 100 windows) of chromosome 1

svg(filename = "../Figures/Figure_3_C_het_all_50Mb_scf1.svg", height = 8, width = 14)

scf=1
par(mfrow=c(5,1), mai=c(0.2, 0.4, 0.2, 0.1))
barplot(apply(Nloci_per_window_SE[which(windows.1M$chr==scf_lengths$scf[scf]),], 1, mean, na.rm=T)[1:100], border=NA, ylim=c(0,105), las=1)
plot(NULL, xlim=c(0, 50000000), ylim=c(0,1), las=1, xlab="", ylab="", main = "", xaxt="n")
for (i in which(colnames(SE.gt)[SEhybrids_in_SE] %in% meta$ID[which(meta$sp.nucl == "whitei")])) {
  points(het_per_window_SE[which(windows.1M$chr==scf_lengths$scf[scf])[1:100], i] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])][1:100], col="#92278F30", type="l")
}
points(sim_het_per_window_SE[which(windows.1M$chr==scf_lengths$scf[scf])][1:100] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])][1:100], col=1, type="l", lwd=2)

plot(NULL, xlim=c(0, 50000000), ylim=c(0,1), las=1, xlab="", ylab="", main = "", xaxt="n")
for (i in which(colnames(SE.gt)[SEhybrids_in_SE] %in% meta$ID[which(meta$sp.nucl == "hybridoSE")])) {
  points(het_per_window_SE[which(windows.1M$chr==scf_lengths$scf[scf])[1:100], i] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])][1:100], col="#EC008C30", type="l")
}
points(sim_het_per_window_SE[which(windows.1M$chr==scf_lengths$scf[scf])][1:100] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])][1:100], col=1, type="l", lwd=2)

plot(NULL, xlim=c(0, 50000000), ylim=c(0,1), las=1, xlab="", ylab="", main = "", xaxt="n")
for (i in 1:length(hybridoNW_in_NW)) {
  points(het_per_window_NW[which(windows.1M$chr==scf_lengths$scf[scf])[1:100],i] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])][1:100], col="#F47F7230", type="l")
}
points(sim_het_per_window_NW[which(windows.1M$chr==scf_lengths$scf[scf])][1:100] ~ windows.1M$begin[which(windows.1M$chr==scf_lengths$scf[scf])][1:100], col=1, type="l", lwd=2)
barplot(apply(Nloci_per_window_NW[which(windows.1M$chr==scf_lengths$scf[scf]),], 1, mean, na.rm=T)[1:100], border=NA, ylim=c(0,105), las=1)

dev.off()

#############################################################
### Heterozygosity distribution in the different lineages ###
#############################################################

meta$het <- apply(sicily.gt, 2, prop.het)
meta$het.auto <- apply(sicily_auto.gt, 2, prop.het)

rred <- meta[which(meta$sp.nucl=="rredtenbacheri"),]
table(rred$pop, rred$sex)

meta$sp_forhet <- NA
meta$sp_forhet[meta$sp.nucl=="rredtenbacheri" & meta$pop %in% c("Patti")] <- "02_rredtenbacheri_sex"
meta$sp_forhet[meta$sp.nucl=="rredtenbacheri" & meta$pop %!in% c("Patti", "Cefalu", "Siracusa")] <- "01_rredtenbacheri_asex"
meta$sp_forhet[meta$sp.nucl=="gmaretimi"] <- "03_gmaretimi"
meta$sp_forhet[meta$sp.nucl=="gbenazzii"] <- "04_gbenazzii"
meta$sp_forhet[meta$sp.nucl=="ggrandii"] <- "05_ggrandii"
meta$sp_forhet[meta$sp.nucl=="atticus"] <- "06_atticus"
meta$sp_forhet[meta$sp.nucl=="whitei"] <- "07_whitei"
meta$sp_forhet[meta$sp.nucl=="lynceorum"] <- "08_lynceorum"

### Figure 5
remotes::install_github("R-CoderDotCom/ridgeline@main")
library(ridgeline)

par(mai=c(0.5, 1.8, 0.1, 0.1))
ridgeline(meta$het.auto, meta$sp_forhet, bw=0.01, palette=c("#0599CE", "#0599CE", "#F9ED32", "#F7941D", "#BE1E2D", "#71BF44", "#92278F", "#8B5E3C"))

anova(aov(meta$het.auto[which(meta$sp_forhet %!in% c("hybridoSE", "hybridoNW", "unknown"))]~meta$sp_forhet[which(meta$sp_forhet %!in% c("hybridoSE", "hybridoNW", "unknown"))]))
TukeyHSD(aov(meta$het.auto[which(meta$sp_forhet %!in% c("hybridoSE", "hybridoNW", "unknown"))]~meta$sp_forhet[which(meta$sp_forhet %!in% c("hybridoSE", "hybridoNW", "unknown"))]))

############################################################
### Allelic depth ratio in B. lynceorum with "missing X" ###
############################################################

# Extract depth values
dp <- extract.gt(sicily.vcf, "DP", as.numeric = T)

# Extract allelic depth values for the reference allele
ad.ref <- extract.gt(sicily.vcf, "AD", as.numeric = T)

# Subset allelic depth values for heterozygous sites only
ad.ref.het <- ad.ref
ad.ref.het[which(sicily.gt != "0/1")] <- NA
ad.ref.het[which(is.na(sicily.gt))] <- NA

# Discard sites with overall depth < 20
ad.ref.het[which(sicily.dp < 20)] <- NA

# Do the same for the alternative allele
ad.alt.het <- dp - ad.ref
ad.alt.het[which(sicily.gt!="0/1")] <- NA
ad.alt.het[which(is.na(sicily.gt))] <- NA
ad.alt.het[which(sicily.dp < 20)] <- NA


# Function to apply a sum like pmin or pmax
# Written by Sven Hohenstein here:
# https://stackoverflow.com/questions/13123638/there-is-pmin-and-pmax-each-taking-na-rm-why-no-psum
psum <- function(...,na.rm=FALSE) {
  dat <- do.call(cbind,list(...))
  res <- rowSums(dat, na.rm=na.rm)
  idx_na <- !rowSums(!is.na(dat))
  res[idx_na] <- NA
  res
}

# compute allelic depth ratio (ADR), defined as
# depth at the allele supported by most reads divided by total depth at the locus.
# takes two matrices with allelic depths (integer) as input
get.ad.ratio <- function(ad.ref, ad.alt) return(pmax(ad.ref, ad.alt) / mapply(psum, ad.ref, ad.alt, na.rm=T))


AD.ratio <- get.ad.ratio(as.matrix(ad.ref.het), as.matrix(ad.alt.het))

### Figure S5: Are the 3 lynceorum with X from 2 species diploid or triploid on the X?
# Autosomes
par(mfrow=c(1,2), mai=c(1, 0.1, 0.1, 0.1))
plot(NULL, xlim=c(0.45, 1), ylim=c(0,6), xlab="Allelic depth ratio (heterozygous alleles only)", ylab="", yaxt="n", las=1)
# take other lynceorum (10 at random, too messy if we take them all)
random_lynceorum <- sample(meta$ID[meta$sp.nucl == "lynceorum" & meta$ID %!in% c("E09", "D32", "D41")], 10)
# take other whitei (10 at random, too messy if we take them all)
random_whitei <- sample(meta$ID[meta$sp.nucl == "whitei"], 10)
for (i in random_lynceorum) {
  points(density(AD.ratio[which(sicily.vcf@fix[,1] != "Brsri_v3_scf2"),which(meta$ID==i)], na.rm=T, bw = 0.02) , type="l", lty=1, col="#000000")
}
for (i in random_whitei) {
  points(density(AD.ratio[which(sicily.vcf@fix[,1] != "Brsri_v3_scf2"),which(meta$ID==i)], na.rm=T, bw = 0.02) , type="l", lty=2, col="#000000")
}
# add focus individuals
points(density(AD.ratio[which(sicily.vcf@fix[,1] != "Brsri_v3_scf2"),which(meta$ID=="E09")], na.rm=T, bw = 0.02) , type="l", lty=1, lwd=3, col="#A3195B")
points(density(AD.ratio[which(sicily.vcf@fix[,1] != "Brsri_v3_scf2"),which(meta$ID=="D32")], na.rm=T, bw = 0.02) , type="l", lty=1, lwd=3, col="#DEDC00")
points(density(AD.ratio[which(sicily.vcf@fix[,1] != "Brsri_v3_scf2"),which(meta$ID=="D41")], na.rm=T, bw = 0.02) , type="l", lty=1, lwd=3, col="#00A19A")


# X chromosome
plot(NULL, xlim=c(0.45, 1), ylim=c(0,6), xlab="Allelic depth ratio (heterozygous alleles only)", ylab="", yaxt="n", las=1)
# take other lynceorum
for (i in random_lynceorum) {
  points(density(AD.ratio[which(sicily.vcf@fix[,1] == "Brsri_v3_scf2"),which(meta$ID==i)], na.rm=T, bw = 0.03) , type="l", lty=1, col="#000000")
}
for (i in random_whitei) {
  points(density(AD.ratio[which(sicily.vcf@fix[,1] == "Brsri_v3_scf2"),which(meta$ID==i)], na.rm=T, bw = 0.03) , type="l", lty=2, col="#000000")
}
# add focus individuals
points(density(AD.ratio[which(sicily.vcf@fix[,1] == "Brsri_v3_scf2"),which(meta$ID=="E09")], na.rm=T, bw = 0.03) , type="l", lty=1, lwd=3, col="#A3195B")
points(density(AD.ratio[which(sicily.vcf@fix[,1] == "Brsri_v3_scf2"),which(meta$ID=="D32")], na.rm=T, bw = 0.03) , type="l", lty=1, lwd=3, col="#DEDC00")
points(density(AD.ratio[which(sicily.vcf@fix[,1] == "Brsri_v3_scf2"),which(meta$ID=="D41")], na.rm=T, bw = 0.03) , type="l", lty=1, lwd=3, col="#00A19A")
