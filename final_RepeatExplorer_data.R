library(vegan)

################
# load data ####
################

setwd("~/git_repos/EuphrasiaGS//")

dat <- read.table("data/COMPARATIVE_ANALYSIS_COUNTS1813_25k.csv", header=T)
head(dat, n=20)
rownames(dat) <- paste0("CL", 1:nrow(dat))


dat[188,] # phage DNA spiked in for sequencing
plCl <- c(70, 57, 39, 37, 41, 51, 60, 79, 92,34, 63, 83) # plastid clusters
contam.counts <- colSums(dat[c(188, plCl),])[3:33]

# remove contaminants
dat <- dat[-c(188, plCl),]
lastCol <- ncol(dat)
# per sample reads in clusters
barplot(colSums(dat)[3:lastCol])

# input read numbers from the comparative counts file:
inReads <- as.numeric(strsplit(readLines("data/COMPARATIVE_ANALYSIS_COUNTS1813_25k.csv", n=1),"\t")[[1]][3:33])
inReads <- inReads - contam.counts # subtract contaminant reads
nams <- names(dat)[3:33]
totalReads <- data.frame(names=nams, inReads=inReads) # contaminations removed
head(totalReads)

# normalising by No of input reads
datNorm <- dat[,3:33]/totalReads$inReads
dim(dat)
dim(datNorm)
head(dat)
head(datNorm)
# Data of individuals without cluster/supercluster numbers, before normalisation
datInds <- dat[,3:lastCol]


# Transate E numbers to something more meaningful. A lookup table:
lookUp <- rbind(c('E033', 'MI0'),
                c('E040', 'RO'),
                c('E032', 'RI'),
                c('E031', 'VI0'),
                c('E022', 'MI2'),
                c('E023', 'MI1'),
                c('E024', 'FO2'),
                c('E025', 'FO1'),
                c('E026', 'AR2'),
                c('E027', 'AR1'),
                c('E028', 'AR0'),
                c('E034', 'X1'),
                c('E035', 'MI3'),
                c('E036', 'AR3'),
                c('E037', 'FO4'),
                c('E038', 'FO3'),
                c('E039', 'X2'),
                c('E030', 'AN0'),
                c('E001', 'AN1'),
                c('E021', 'VI1'),
                c('E020', 'TE'),
                c('E019', 'OS'),
                c('E018', 'FO0'),
                c('E017', 'FH'),
                c('E016', 'NE'),
                c('E002', 'AR4'),
                c('E005', 'AR5'),
                c('E003', 'MI4'),
                c('E006', 'MI5'),
                c('E009', 'CU'))
colnames(datInds)
humanLabels <- sapply(colnames(datInds), function(x) {
  if(x %in% lookUp[,1]){
    lookUp[which(lookUp[,1] == x), 2]
  } else x
})

colnames(datInds) <- humanLabels

totalReads[,1] <- sapply(totalReads[,1], function(x) {
  if(x %in% lookUp[,1]){
    lookUp[which(lookUp[,1] == x), 2]
  } else x
})


totalReads <- data.frame(samples=totalReads[,1], readCount=as.numeric(totalReads[,2]), stringsAsFactors = F)
names(datInds) <- humanLabels
humanLabels %in% c("AN0","AN1", "VI0", "VI1", "CU", "RO", "RI")

# Data all set up.

######################################################################################
# Is the genomic proportion if reads in clusters different between ploidy levels? ####
######################################################################################
barplot(colSums(datInds) / totalReads[,2], las=3, main ="Proportion of reads in RE clusters")


# Does total repeat content differ between diploids and tetraploids?
summary(lm((colSums(datInds) / totalReads[,2])[2:31] ~ (humanLabels %in% c("AN0","AN1", "VI0", "VI1", "CU", "RO", "RI"))[2:31]))
# Yes!
mean((((colSums(datInds) / totalReads[,2]))[(humanLabels %in% c("AN0","AN1", "VI0", "VI1", "CU", "RO", "RI"))]))
mean((((colSums(datInds) / totalReads[,2]))[!(humanLabels %in% c("AN0","AN1", "VI0", "VI1", "CU", "RO", "RI"))])[2:24])

colnames(datNorm) <- humanLabels

# Centre and scale so that variance within each clusert is 1
#datIndsNorm <- t(apply(datIndsNorm, 1, function(x) (x - mean(x)) /sd(x)))
# This leads to the samples  clustering by data set. Skimming vs high-cov (different library types)
# Don't scale. Centring is done automatically by prcomp.

datIndsNorm <- datNorm

# Clustering largely by ploidy
heatmap(as.matrix(datIndsNorm[1:100,]), main = "Read counts per species and RE cluster")
heatmap(t(datIndsNorm[1:100,]), main = "Read counts per species and RE cluster")
heatmap(t(datIndsNorm[1:200,]), main = "Read counts per species and RE cluster")

###########
# PCAs ####
###########

# When Bartsia is included, the Barsia-Euphrasia contrast accounts for much variance.
pc1 <- prcomp(t(datIndsNorm[1:100,]))

str(pc1)
summary(pc1)

plot(pc1$x[,1:2], type="n", xlab="PC1, 50%", ylab="PC2, 18%")
text((pc1$x[,1:2]), labels = rownames(pc1$x),
     col= ifelse(rownames(pc1$x) %in% c("VI1","VI0","RO", "AN0", "AN1", "RI", "CU"), "black", "grey"))

# dropping Bartsia
head(datIndsNorm)
dim(datIndsNorm)
pc11 <- prcomp(t(datIndsNorm[1:100,2:31]))
summary(pc11)
str(pc11)

# figure 3b ####
plot(pc11$x[,1:2], xlab="PC1, 34%", ylab="PC2, 25%")
text((pc11$x[,1:2]), labels = rownames(pc11$x),
     col = ifelse(rownames(pc11$x) %in% c("VI1","VI0","RO", "AN0", "AN1", "RI", "CU"), "black", "grey"))


# figure 3c ####
plot(pc11$rotation[,1], main="Elements of eigenvector 1", xlab="Element position", ylab="Value")
abline(h=c(0.1, -0.1), lty=2)
outliers1 <- which(abs(pc11$rotation[,1]) > 0.1) 
text(cbind(outliers1, pc11$rotation[outliers1,1]), labels = names(outliers1))

which(abs(pc11$rotation[,1]) > 0.1)
# alternative:
biplot(pc11)


# figure 3a ####
plot(hclust(dist(t(datIndsNorm[1:100,1:31]),method = "euc")),
     main="Hierarchical clustering by genomic repeats",
     xlab="")


###################################################
# Element types and difference between samples ####
###################################################
anno <- read.table("data/CLUSTER_TABLE_1813_25k.csv", skip = 6, header=T, stringsAsFactors = F)
head(anno)

plCl <- c(70, 57, 39, 37, 41, 51, 60, 79, 92,34, 63, 83)
# CL 188 is phage DNA spiked in at sequencing
anno <- anno[-c(188, plCl),]

# first 100 clusters only
anno <- anno[1:100,]
head(anno)
table(anno$TAREAN_annotation)
anno$Automatic_annotation[startsWith(anno$TAREAN_annotation, "Putative sat")] <- "All/repeat/satellite"

datAnno <- t(t(dat[1:100,3:33])/totalReads$readCount)
head(datAnno)
head(anno)
# anno and data have the same order and clusters:
all(anno$Cluster[1:100] == dat$cluster[1:100])
colnames(datAnno)
colnames(datAnno) <- humanLabels
# a list of repeat types with relative genome abundance in each sample 
annoSums <- tapply(1:100, anno$Automatic_annotation, function(x){
  
  ifelse(length(x) > 1, return(colSums(datAnno[x,])), return(datAnno[x,]))
  
})
length(annoSums)
#write.table(do.call(rbind, annoSums), "annoSumsOut")
levels(as.factor(anno$Automatic_annotation))
barplot(annoSums$`All/repeat/mobile_element/Class_I/LTR/Ty1_copia/Angela`[2:31],
        las=3)

barplot(annoSums$`All/repeat/rDNA/45S_rDNA`, las=3)
barplot(annoSums$`All/repeat/mobile_element/Class_I/LTR/Ty1_copia/Ale`, las=3)
barplot(annoSums$`All/repeat/rDNA/5S_rDNA`, las=3)
barplot(annoSums$`All/repeat/mobile_element/Class_I/LTR/Ty1_copia/Angela`, las=3)
summary(annoSums$`All/repeat/mobile_element/Class_I/LTR/Ty1_copia/Angela`[2:31])
names(annoSums)
length(names(annoSums))
# 19, exclude the general ones
names(annoSums)[c(3:19)]
# only 17 categories now

barplot(colSums(do.call(rbind, annoSums)), las=3)
# All annotated repeat types
barplot((do.call(rbind, annoSums)), las=3)

# DF of repeat type abundance in each sample
type.df <- as.data.frame(t(do.call(rbind, annoSums)))

# differences between Bartsia alpina and all other samples for each repeat type
t.test.results <- apply(type.df, 2, function(x){
  (t.test(x[2:31]-x[1]))
})
# Corrected for multiple testing, they differ in most repeat types
(sapply(t.test.results, function(x) x$p.value*19<0.05))

# Make a DF with genomic proportions of each repeat type for each sample
type.df <- data.frame(code=humanLabels, species=c("Bartsia alpina", "E. anglica", "E. arctica", "E. micrantha",
                                       "E. arctica", "E. micrantha", "E. cuspidata", "E. nemorosa",
                                       "E. fharaidensis", "E. foulaensis", "E. ostendeldii", "E. tetraquetra",
                                       "E. vigursii", "E. micrantha", "E. micrantha", "E. foulaensis", 
                                       "E. foulaensis", "E. arctica", "E. arctica", "E. arctica", "E. anglica",
                                       "E. vigursii", "E. rivularis", "E. micrantha", "E. tetraploid hybrid",
                                       "E. micrantha", "E. arctica", "E. foulaensis", "E. foulaensis",
                                       "E. tetraploid hybrid", "E. rostkoviana"), type.df)
head(type.df)
table(type.df$species)
#write.table(type.df, "repeat_type_proportions.csv", sep="\t", row.names = F)

# anovas

# logical vector, F=tetraploid, T=diploid
dipTet <- ifelse(humanLabels %in% c("AN001", "CUS009", "VI","AN","RO","RI","VI021"), T, F)[2:31]

# Differences between diploids and tetraploids?
cbind(sapply(3:19, function(x){
        (coef(summary(lm(annoSums[[x]][2:31]~dipTet)))[2,4]*17) < 0.05
      }),
levels(as.factor(anno$Automatic_annotation))[3:19])
which(anno$Automatic_annotation=="All/repeat/mobile_element/Class_I/LTR/Ty1_copia/Ale")

# Yes, Ale (copia) and 45S rDNA!
summary(aov(lm(annoSums$`All/repeat/mobile_element/Class_I/LTR/Ty1_copia/Ale`[2:31]~dipTet)))
summary(aov(lm(annoSums$`All/repeat/rDNA/45S_rDNA`[2:31]~dipTet)))
barplot(annoSums$`All/repeat/rDNA/45S_rDNA`[2:31])

# micrantha compared to other tetraploids
which(!(humanLabels %in% c("AN1", "CU", "VI1","AN0","RO","RI","VI0", "Balp")))
startsWith(humanLabels[!(humanLabels %in% c("AN1", "CU", "VI1","AN0","RO","RI","VI0", "Balp"))], "MI")
data.frame(name=names(annoSums)[3:19],
           differentInMic=sapply(3:19, function(x){
             (coef(summary(lm(annoSums[[x]][which(!(humanLabels %in% c("AN1", "CU", "VI1","AN0","RO","RI","VI0", "Balp")))]~startsWith(humanLabels[!(humanLabels %in% c("AN1", "CU", "VI1","AN0","RO","RI","VI0", "Balp"))], "MI"))))[2,4]*17) < 0.05
           })
)
tetsat <- annoSums$`All/repeat/satellite`[which(!(humanLabels %in% c("AN1", "CU", "VI1","AN0","RO","RI","VI0", "Balp")))]

# logical vector of tetraploids: T=E. micrantha, F=other species
micBoolTetSat <- startsWith(humanLabels[!(humanLabels %in% c("AN1", "CU", "VI1","AN0","RO","RI","VI0", "Balp"))], "MI")

# 0.00436 * 17 repeat types > 0.07, not significant
summary((lm(tetsat~micBoolTetSat)))



###########################################
# overlapping pops GS and RE data sets ####
###########################################

#overlap.gs # from GS script:
#      AR0      AR1      AR2      FO3      FO2      FO1      MI0      MI3      MI1 
# 2.310000 2.508255 2.431333 1.998000 2.477444 2.370333 2.353333 2.365222 2.398167
overlap.gs <- c(2.310000, 2.508255, 2.431333, 1.998000, 2.477444, 2.370333, 2.353333, 2.365222, 2.398167)
names(overlap.gs) <- c("AR0", "AR1", "AR2", "FO3", "FO2", "FO1", "MI0", "MI3", "MI1")
overlap.gs
dat
library(MASS)
head(dat)
dat.overlap <- dat[,c("E028","E027","E026","E038","E024","E025","E033","E035","E023")]
dat.overlap <- as.data.frame(t(dat.overlap))
dim(dat.overlap)
dat.overlap[1:9,1:10]

nclst <- 1000 # 100, 200, 500
tests <- apply(dat.overlap[,1:nclst], 2, function(x) {cor.test(overlap.gs, x)})
str(tests[1])
cor.results <- data.frame(t(sapply(tests, function(x) c(p=x$p.value, r=x$estimate)))[,1:2])
min(cor.results$p, na.rm = T) * nclst

plot(cor.results)  
plot(sort(cor.results$p)*nclst)
plot(sort(cor.results$p)*nclst, ylim=c(0,10))
abline(h=0.05)
# There are no low-enough p values to consider any of the 1000 correlations
#  significant. Neither works with 100, 200, or 500 clusters.


######################
# Repeat overview ####
######################

dim(datIndsNorm)
dIN100 <- datIndsNorm[1:101,2:31]
dIN1001 <- t(apply(dIN100[-64,], 1, function(x) x/max(x)))

image(dIN1001)
library(gplots)
heatmap.2(dIN1001, trace="none", dendrogram="col", Rowv=F)

