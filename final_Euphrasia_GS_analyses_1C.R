library(nlme)
library(ggplot2)
library(vegan)
library(geosphere)
library(emmeans)
library(RColorBrewer)
###################
# Cleaned data ####
###################

#setwd("~/Dropbox/manuscripts/2008_Euphrasia_genome_size/data_and_code/")
#setwd("~/github_repos/xyz/")

dat <- read.table("data/210212_GS_cleaned.csv", stringsAsFactors = F, sep=",", header=T)
head(dat)

table(dat$Sample.identity)
table(dat$ploidy)

# lists of diploid and tetraploid non-hybrid taxa
sp2 <- names(table(dat$Sample.identity[dat$hybrid==F & dat$ploidy=="dip"]))
sp4 <- names(table(dat$Sample.identity[dat$hybrid==F & dat$ploidy=="tet"]))



# plot ordered by ploidy
ord <- order(dat$pg1C)
# figure 2b ####
plot(dat$pg1C[ord], type='n')
grid(nx = NA, ny=NULL)
points(dat$pg1C[ord])
# end of fig 

# Justification for singling out tetraploids with large and small GS:
#  the distribution is discontinuous.
hist(dat$pg1C, breaks=100)
abline(v=c(1.07, 1.26))


# labels for legend
legLabs <- c(levels(as.factor(dat$Sample.identity[ord])), "")
# order of labels
legOrd <- c(1, 2, 15, 20, 21, 23, 24, 3:14,16:19, 22)

# Make data.frame with one row per population: sample no, gs mean, min and max ####
popStats <- tapply(dat$pg1C, dat$specPop, function(x){
  c(mean(x), length(x), min(x), max(x))
})
popStats <- as.data.frame(do.call(rbind, popStats))
head(popStats)
names(popStats) <- c("popMean", "popN", "popMin", "popMax")

popStats$taxon <- as.factor(unlist(lapply(strsplit(rownames(popStats), "[.]"),
                                function(x) paste(x[1:2], collapse = "."))))

popStats$ploidy <- "tet"
popStats$ploidy[popStats$taxon %in% sp2] <- "dip"
popStats$ploidy[1] <- "dip" # one hybrid combination is diploid (anglica x nemorosa)

popStatsOrder <- order(popStats$ploidy, as.character(popStats$taxon), popStats$popMean)
popStats[popStatsOrder,]
popStatsOrderM <- order(popStats$popMean)
unique(popStats$taxon)
head(popStats)

# Make a colour palette. Replace yellow, which is hard to see.
plot(1:8, col = brewer.pal(8, "Set1"), pch=19)
pal8 <- brewer.pal(8, "Set1")
pal8[6] <- "#AAAAAA"
plot(1:8, col = pal8, pch=19)
# figure 1a ####
plot(popStats[popStatsOrder,1],
     pch=as.numeric(popStats$taxon)[popStatsOrder]%%7+1,
     col = pal8[as.numeric(popStats$taxon)[popStatsOrder]%%8+1],
     type='n', xlab="Population, ordered by taxon and genome size", ylab="Genome size (pg/1C)")
grid(nx = NA, ny=NULL)
for(i in 1:nrow(popStats)){
  points(c(i, i),
         c(popStats[popStatsOrder[i], 3], popStats[popStatsOrder[i], 4]),
         type='l')
}
points(popStats[popStatsOrder,1],
     pch=as.numeric(popStats$taxon)[popStatsOrder]%%7+1,
     col = pal8[as.numeric(popStats$taxon)[popStatsOrder]%%8+1])

legend("bottomright",
       col=pal8[as.numeric(unique(popStats$taxon[popStatsOrder]))%%8+1],
       pch=as.numeric(unique(popStats$taxon[popStatsOrder]))%%7+1,
       legend=unique(as.character(popStats$taxon)[popStatsOrder]),
       ncol = 2)
### end figure 1a

# per-taxon stats ####
# 1st subsection of results
table(dat$Sample.identity, dat$hybrid)
colSums(apply(table(dat$Sample.identity, dat$hybrid), c(1, 2), as.logical))
table(dat$ploidy)

lmPloyMeans <- lm(dat$pg1C~ dat$ploidy)
summary(lmPloyMeans)
# Mean GS per ploidy level
emmeans(lmPloyMeans, ~ploidy)
ployMeans <- summary(emmeans(lmPloyMeans, ~ploidy))[,2]
# 11% deficiency in tet GS compared to 2* the dipl size:
((2 * ployMeans[1]) - ployMeans[2]) /ployMeans[2]

# array with taxon stats
meanSdObs_spec_ext <- do.call(rbind,tapply(dat$pg1C, dat$Sample.identity, function(x){
  c(mean(x), min(x), max(x), sd(x), length(x),
    sd(x)/length(x), max(x)/min(x), sd(x)/mean(x))
  }
  )
  )
colnames(meanSdObs_spec_ext) <- c("mean","min","max","sd","obs", "se","ratio",
                                  "cv")
head(meanSdObs_spec_ext)

barplot(meanSdObs_spec_ext[,7], las=3, ylab="Ratio between max and min GS")
abline(h=c(1,1.1,1.2))

barplot(meanSdObs_spec_ext[,8], las=3, ylab="Coefficient of variation of GS")

# geographic distance between min and max individuals within taxa
head(dat)
kmDist <- tapply(1:192, dat$Sample.identity, function(x){
  distm(dat[x,][which(dat[x,]$pg1C == min(dat[x,]$pg1C)),c(6,5)],
  dat[x,][which(dat[x,]$pg1C == max(dat[x,]$pg1C)),c(6,5)])[1]/1000
})
# ms mentions those species which were mentioned befor for their large GS ranges
barplot(kmDist, las=3)



# order by taxon mean GS
ord_spec <- order(meanSdObs_spec_ext[,1])

# Summary stats by taxon and pop
meanSdObs_specPop <- do.call(rbind,tapply(dat$pg1C, dat$specPop,
                    function(x) c(mean(x), sd(x), length(x), sd(x)/length(x))))

colnames(meanSdObs_specPop) <- c("mean","sd","obs", "se")
head(meanSdObs_specPop)
ord_specPop <- order(meanSdObs_specPop[,1])


# pop mean GS by size
plot(meanSdObs_specPop[ord_specPop,1])
# Taxon mean GS by size
plot(meanSdObs_spec_ext[ord_spec,1])

rownames(meanSdObs_specPop)
head(meanSdObs_specPop)
meanSdObs_specPop <- data.frame(meanSdObs_specPop, hyb=0, ploy=4, stringsAsFactors = F) 

meanSdObs_specPop$hyb[grep(" x ", rownames(meanSdObs_specPop))] <- 1
meanSdObs_specPop$hyb[grep("hybrid", rownames(meanSdObs_specPop))] <- 1

meanSdObs_specPop$ploy[sapply(rownames(meanSdObs_specPop),
    function(x) paste(strsplit(x,"[.]")[[1]][1:2],collapse='.')) %in% sp2] <- 2
meanSdObs_specPop$ploy[sapply(rownames(meanSdObs_specPop),
    function(x) paste(strsplit(x,"[.]")[[1]][1:2],collapse='.')) %in% sp4] <- 4
# the only hybrid with ploidy 2:
meanSdObs_specPop[1,]$ploy <- 2

head(meanSdObs_specPop)
meanSdObs_specPop <- data.frame(meanSdObs_specPop, specPop=rownames(meanSdObs_specPop), stringsAsFactors = F)
head(dat)
# do all spec pop combinations have the same lat within groups?
all(tapply(dat$Lat, dat$specPop, function(x) all(x[1]==x)))
all(tapply(dat$Long, dat$specPop, function(x) all(x[1]==x)))
# yes, good

lats <- tapply(dat$Lat, dat$specPop, function(x) x[1])
lons <- tapply(dat$Long, dat$specPop, function(x) x[1])
lats <- data.frame(specPop=rownames(lats), lat=lats, stringsAsFactors = F) 
lons <- data.frame(specPop=rownames(lons), long=lons, stringsAsFactors = F) 


# by lat
specPopLat <- merge(meanSdObs_specPop, lats, by="specPop")
head(specPopLat)
specPopLatLon <- merge(specPopLat, lons, by="specPop")
head(specPopLatLon)
head(dat)
datMerged <- merge(dat, specPopLatLon, by="specPop", all = T)



# max geographic distance between samples of taxon
specPopLatLonTax <- specPopLatLon

specPopLatLonTax$tax <- sapply(specPopLatLon$specPop, function(x){
  a <- strsplit(x, "[.]")[[1]]
  b <- paste0(a[1:(length(a)-1)], collapse = "")
})
tapply(1:nrow(specPopLatLonTax), specPopLatLonTax$tax, function(x) length(x))
maxDists <- tapply(1:nrow(specPopLatLonTax), specPopLatLonTax$tax, function(x){
  nr <- nrow(specPopLatLonTax[x,]) 
  if(nr == 1) {0}
  else{
  dists <- numeric(0)
    for(i in 1:(nr-1)){
    for(j in (i+1):nr){
      dists <- append(dists, distm(specPopLatLonTax[x,9:8][i,], specPopLatLonTax[x,9:8][j,]))
    }
    }

  max(dists)/1000
  }
      

})

barplot(maxDists, las=3)

# INspect on Google and add Country information
specPopLatLon$country <- "Eng"
specPopLatLon[specPopLatLon$long < -6.5,]
specPopLatLon[specPopLatLon$long < -6.5,]$country <- "RIr"
specPopLatLon[specPopLatLon$lat > 55.228,] 
specPopLatLon[specPopLatLon$lat > 55.228,]$country <- "Sco"
specPopLatLon[specPopLatLon$lat < 50,]$country <- "CIs"
specPopLatLon[specPopLatLon$lat > 55.1 & specPopLatLon$long < -3 ,]$country <- "Sco"

specPopLatLon[specPopLatLon$lat > 51.33 &
                specPopLatLon$lat < 53.3 &
                specPopLatLon$long < -3.1 &
                specPopLatLon$long > -5,]$country <- "Wal"
plot(lat~long, data=specPopLatLon, col=factor(country))


# STATS MODELS ####

# anovas ####

head(dat)

plot(dat$pg1C)
partitioningGS1 <- aov(pg1C~ploidy+Sample.identity+Population_number, data=dat)
partitioningGS1nh <- aov(pg1C~ploidy+Sample.identity+Population_number, data=dat[dat$hybrid==F,])
partitioningGS2 <- aov(pg1C~ploidy+Sample.identity, data=dat)
partitioningGS3 <- aov(pg1C~ploidy+Population_number+Sample.identity, data=dat)
partitioningGS4 <- aov(pg1C~ploidy, data=dat)

partitioningLat1 <- aov(Lat~ploidy, data=dat)

s1 <- summary(aov(partitioningGS1))
summary(aov(partitioningGS1)) # sample.id (= taxon) before population (which is nested within sample.id)
summary(aov(partitioningGS1nh))

summary(aov(partitioningGS2))
summary(aov(partitioningGS3)) # population and sample.id the wrong way round
summary(aov(partitioningGS4))
summary(aov(partitioningLat1))


par(mfrow=c(2,2))
plot(partitioningGS1)
plot(partitioningGS1nh)
plot(partitioningGS2)
plot(partitioningGS3)
plot(partitioningGS4)
plot(partitioningLat1)
par(mfrow=c(1,1))


# ratio of Sums of Squares (varience explained by taxon and population)
a <- summary(partitioningGS1)
b <- summary(partitioningGS1nh)
a[[1]]$`Sum Sq`
b[[1]]$`Sum Sq`

# ratio of variance explained by 'population' and 'species'
a[[1]]$`Sum Sq`[3]/a[[1]]$`Sum Sq`[2] # smaller when hybrids are included
# 3.006604
b[[1]]$`Sum Sq`[3]/b[[1]]$`Sum Sq`[2] # larger without hybrids
# 8.001641 ratio is higher, i.e. taxon explains less when hybrids are excluded


# "Taxon" explaines less variance when hybrids are excluded. Do hybrids have
#   less variation (coefficient of variation) in their GS vals?
GScv <- tapply(dat$pg2C, dat$Sample.identity, function(x){sd(x)/mean(x)})
GShyb <- grepl("x",names(GScv))
GSdf <- data.frame(cv=GScv, hyb=GShyb, ploidy=c(2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 2, 2, 4, 2))
summary(lm(GScv~GShyb))
# Yes! But this may be becasue they usually have one population only.

# Differences in cv of GS between ploidy levels? Not significant (no matter
#  whether with or without hybrids):
with(GSdf[GSdf$hyb==F,],{
  lm1 <- lm(cv~ploidy)
  print(summary(lm1))
})
with(GSdf,{
  lm1 <- lm(cv~ploidy)
  print(summary(lm1))
})


# Do arctica and micrantha have different variances?
# compare valriance with an F-test (var.test() in R)

arcMeans <- specPopLatLon$mean[startsWith(specPopLatLon$specPop, "E. arctica.")]
micMeans <- specPopLatLon$mean[startsWith(specPopLatLon$specPop, "E. micrantha.")]
angMeans <- specPopLatLon$mean[startsWith(specPopLatLon$specPop, "E. anglica.")]
mean(arcMeans)
mean(micMeans)
sd(micMeans)
sd(arcMeans)

var.test(arcMeans, micMeans, alternative = "gr")
var.test(angMeans, micMeans, alternative = "gr")
# both samples are not significantly different from normal
shapiro.test(arcMeans)
shapiro.test(micMeans)
# is the variance in GS partitioned by pop in arctica and micrantha?
dim(dat)
head(dat)
summary(aov(pg1C~specPop, data=dat[dat$Sample.identity=="E. micrantha",]))
# not in micrantha
summary(aov(pg1C~specPop, data=dat[dat$Sample.identity=="E. arctica",]))
# yes in arctica
summary(aov(pg1C~specPop, data=dat[dat$Sample.identity=="E. anglica",]))
# yes in anglica
dim(dat[dat$Sample.identity=="E. arctica",])
dim(dat[dat$Sample.identity=="E. micrantha",])
dim(dat[dat$Sample.identity=="E. anglica",])
sort(tapply(dat$Population_number, dat$Sample.identity, function(x) length(unique(x))))
# latitude ####

head(specPopLat)

specPopLat <- data.frame(specPopLat, spec=sapply(specPopLat$specPop,
                function(x) paste(strsplit(x, "[.]")[[1]][1:2], collapse = '.')), stringsAsFactors = F)

# over all ploidy levels
summary(aov(lm(pg1C~Lat,data=dat)))
me1 <- lme(mean~ploy+lat, random=~1|spec,data=specPopLat)
summary(me1)
dim(specPopLat)

# Fig S2 ####
plot(mean~lat, data=specPopLat,
     col = as.numeric(as.factor(specPopLat$spec)),
     pch = as.numeric(as.factor(specPopLat$spec)),
     main="GS by latitude",
     ylab="GS (pg/1C)",
     xlab="Latitude of population")
grid(nx=NA, ny=NULL)
legend("bottomright",
       col = c(1:23,0)[legOrd],
       pch = c(1:23,0)[legOrd],
       legend = legLabs[legOrd],
       ncol = 2)
#### end fig S2
as.numeric(as.factor(specPopLat$spec))
head(specPopLat)
unique(specPopLat$spec)
# tetraploids only
me4s <- lme(mean~lat,random=~1|spec, data=specPopLat[specPopLat$ploy==4,])
summary(me4s)


# diploids only
plot(mean~lat, data=specPopLat[specPopLat$ploy==2,],
     col = as.factor(specPopLat[specPopLat$ploy==2,]$spec), 
     main="Diploid GS by lat coloured by species")
me2s <- lme(mean~lat,random=~1|spec, data=specPopLat[specPopLat$ploy==2,])
summary(me2s)
# no difference with or without taxon included

# E. arctica only
head(specPopLat)

lmArc <- lm(mean~ lat, data=specPopLat[specPopLat$spec=="E. arctica",])
summary(aov(lmArc))
summary((lmArc))

lmMic <- lm(mean~ lat, data=specPopLat[specPopLat$spec=="E. micrantha",])
summary(aov(lmMic))
summary((lmMic))

lmAng <- lm(mean~ lat, data=specPopLat[specPopLat$spec=="E. anglica",])
summary(aov(lmAng))
summary((lmAng))


# fig 2c ####
with(specPopLat[specPopLat$spec %in% c("E. anglica", "E. arctica", "E. micrantha"), ],{
  plot(mean~lat, col = as.factor(spec),
       pch = as.numeric(as.factor(spec)), xlab="Latitude in degrees", ylab="1C in pg",
       type='n')#, ylim=c(1.2,2.7))
  grid()
  points(mean~lat, col = pal8[c(2, 4, 6)][as.factor(spec)],
         pch = c(2, 4, 7)[as.numeric(as.factor(spec))], xlab="Latitude", ylab="2C in pg")
  abline(lmAng, lty=2, col=pal8[2])
  abline(lmArc, col=pal8[4])
  legend("right",pch=c(2, 4, 7), col=pal8[c(2, 4, 6)],legend=c("E. anglica","E. arctica","E. micrantha"))
})
# fig end

par(mfrow=c(1,1))

# country ####
par(mfrow=c(2,2))
head(specPopLatLon)
lmEngScot <- lm(mean~country, data=specPopLatLon[specPopLatLon$country %in% c("Eng","Sco") &
                                                   specPopLatLon$ploy==4,])
summary(lmEngScot)
plot(lmEngScot)
summary(aov(lmEngScot))

lmCountry <- lm(mean~country, data=specPopLatLon[specPopLatLon$ploy==4,])

plot(lmCountry)
summary(aov(lmCountry))
par(mfrow=c(1,1))

# mantel tests ####

# mantel tests per individual
dim(dat)
head(dat)
table(dat$ploidy)
table(dat$Sample.identity)


dips.int.mant <- dat[dat$ploidy=="dip",]
tets.int.mant <- dat[dat$ploidy=="tet",]
arcs.int.mant <- dat[dat$Sample.identity %in% c("E. arctica"),]
mics.int.mant <- dat[dat$Sample.identity=="E. micrantha",]

head(dips.int.mant)
dips.int.mant[,6:5]
dip.diffs.ind <- matrix(0, nrow(dips.int.mant), nrow(dips.int.mant))
dip.dists.ind <- matrix(0, nrow(dips.int.mant), nrow(dips.int.mant))
tet.diffs.ind <- matrix(0, nrow(tets.int.mant), nrow(tets.int.mant))
tet.dists.ind <- matrix(0, nrow(tets.int.mant), nrow(tets.int.mant))
arc.diffs.ind <- matrix(0, nrow(arcs.int.mant), nrow(arcs.int.mant))
arc.dists.ind <- matrix(0, nrow(arcs.int.mant), nrow(arcs.int.mant))
mic.diffs.ind <- matrix(0, nrow(mics.int.mant), nrow(mics.int.mant))
mic.dists.ind <- matrix(0, nrow(mics.int.mant), nrow(mics.int.mant))
for(i in 1:nrow(dips.int.mant)){
  for(j in 1:nrow(dips.int.mant)){
    if(i!=j){
      dip.diffs.ind[i,j] <- abs(dips.int.mant$pg2C[i]-dips.int.mant$pg2C[j])
      dip.dists.ind[i,j] <- distm(dips.int.mant[i,6:5],dips.int.mant[j,6:5],fun=distHaversine) # haversine assumes spherical earth
    }
  }
}
for(i in 1:nrow(tets.int.mant)){
  for(j in 1:nrow(tets.int.mant)){
    if(i!=j){
      tet.diffs.ind[i,j] <- abs(tets.int.mant$pg2C[i]-tets.int.mant$pg2C[j])
      tet.dists.ind[i,j] <- distm(tets.int.mant[i,6:5],tets.int.mant[j,6:5],fun=distHaversine) # haversine assumes spherical earth
    }
  }
}
for(i in 1:nrow(arcs.int.mant)){
  for(j in 1:nrow(arcs.int.mant)){
    if(i!=j){
      arc.diffs.ind[i,j] <- abs(arcs.int.mant$pg2C[i]-arcs.int.mant$pg2C[j])
      arc.dists.ind[i,j] <- distm(arcs.int.mant[i,6:5],arcs.int.mant[j,6:5],fun=distHaversine) # haversine assumes spherical earth
    }
  }
}
for(i in 1:nrow(mics.int.mant)){
  for(j in 1:nrow(mics.int.mant)){
    if(i!=j){
      mic.diffs.ind[i,j] <- abs(mics.int.mant$pg2C[i]-mics.int.mant$pg2C[j])
      mic.dists.ind[i,j] <- distm(mics.int.mant[i,6:5],mics.int.mant[j,6:5],fun=distHaversine) # haversine assumes spherical earth
    }
  }
}


plot(dip.diffs.ind, dip.dists.ind)
plot(tet.diffs.ind, tet.dists.ind)
plot(arc.diffs.ind, arc.dists.ind)


# per population
# distances
dips <- specPopLatLon[specPopLatLon$ploy==2,]
tets <- specPopLatLon[specPopLatLon$ploy==4,]
arcs <- specPopLatLon[sapply(specPopLatLon$specPop, function(x) {
  paste(strsplit(x, "[.]")[[1]][1:2], collapse='.')
}) == "E. arctica",]
mics <- specPopLatLon[sapply(specPopLatLon$specPop, function(x) {
  paste(strsplit(x, "[.]")[[1]][1:2], collapse='.')
})=="E. micrantha",]
mics <- specPopLatLon[sapply(specPopLatLon$specPop, function(x) {
  paste(strsplit(x, "[.]")[[1]][1:2], collapse='.')
}) == "E. micrantha",]




dip.diffs <- matrix(0, nrow(dips), nrow(dips))
dip.dists <- matrix(0, nrow(dips), nrow(dips))
tet.diffs <- matrix(0, nrow(tets), nrow(tets))
tet.dists <- matrix(0, nrow(tets), nrow(tets))
arc.diffs <- matrix(0, nrow(arcs), nrow(arcs))
arc.dists <- matrix(0, nrow(arcs), nrow(arcs))
mic.diffs <- matrix(0, nrow(mics), nrow(mics))
mic.dists <- matrix(0, nrow(mics), nrow(mics))
for(i in 1:nrow(dips)){
  for(j in 1:nrow(dips)){
    if(i!=j){
      dip.diffs[i,j] <- abs(dips$mean[i]-dips$mean[j])
      dip.dists[i,j] <- distm(dips[i,c("long","lat")],dips[j,c("long","lat")],fun=distHaversine) # haversine assumes spherical earth
    }
  }
}
for(i in 1:nrow(tets)){
  for(j in 1:nrow(tets)){
    if(i!=j){
      tet.diffs[i,j] <- abs(tets$mean[i]-tets$mean[j])
      tet.dists[i,j] <- distm(tets[i,c("long","lat")],tets[j,c("long","lat")],fun=distHaversine)
    }
  }
}
for(i in 1:nrow(arcs)){
  for(j in 1:nrow(arcs)){
    if(i!=j){
      arc.diffs[i,j] <- abs(arcs$mean[i]-arcs$mean[j])
      arc.dists[i,j] <- distm(arcs[i,c("long","lat")],arcs[j,c("long","lat")],fun=distHaversine)
    }
  }
}
for(i in 1:nrow(mics)){
  for(j in 1:nrow(mics)){
    if(i!=j){
      mic.diffs[i,j] <- abs(mics$mean[i]-mics$mean[j])
      mic.dists[i,j] <- distm(mics[i,c("long","lat")],mics[j,c("long","lat")],fun=distHaversine) # haversine assumes spherical earth
    }
  }
}
# Mantel results ####
set.seed(123345)
mantel(dip.diffs.ind, dip.dists.ind)
mantel(tet.diffs.ind, tet.dists.ind)
mantel(arc.diffs.ind, arc.dists.ind)
mantel(mic.diffs.ind, mic.dists.ind)

mantel(dip.diffs, dip.dists)
mantel(tet.diffs, tet.dists)
mantel(arc.diffs, arc.dists)
mantel(mic.diffs, mic.dists)


plot(dip.diffs, dip.dists)
plot(tet.diffs, tet.dists)
plot(arc.diffs, arc.dists)
plot(mic.diffs, mic.dists)

# link pops, GS, and RE ####

# overlap.df <- dat[dat$Population_number %in% c('FIH072','FIH092','FIMXB024',
# 'FIH098','FIA105','FIA111','FIA107','BERW1','AT1401','AT1405'),] # including one anglica
overlap.df <- dat[dat$Population_number %in% c('FIH072','FIH092','FIMXB024',
        'FIH098','FIA105','FIA111','FIA107','BERW1','AT1401'),] # tetraploids only


barplot(tapply(overlap.df$pg1C, overlap.df$specPop, mean), las=3)
plot(do.call(rbind, tapply(overlap.df$pg1C, overlap.df$specPop,
                           function(x) c(mean(x), sd(x)/mean(x)))),
     xlab="Mean GS", ylab="CV over populations")
overlap.gs <- tapply(overlap.df$pg2C, overlap.df$specPop, mean)
names(overlap.gs) <- c("AR0", "AR1","AR2","FO3","FO2","FO1","MI0","MI3","MI1")
barplot(overlap.gs)


