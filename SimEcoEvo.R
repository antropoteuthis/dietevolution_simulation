#SIMULATING CAUSAL PATHWAYS FOR DIET DETERMINATION

rm(list=ls())
par(ask=F)
#Load libraries
library(reshape2)
library(dplyr)
library(tidyverse)
library(stringr)
library(vegan)
library(ape)
library(phangorn)
library(phytools)
library(picante)
library(geiger)
library(phylobase)
library(fields)
library(phylosignal)
library(OUwie)
library(geomorph)

#Directory setup
setwd("~/Dropbox/Simulation_EcoEvo")
getwd()

#TREE BUILD

rRNA18S = read.nexus('~/Dropbox/bipartitions_root.boot10')
transcriptome = read.tree('~/Dropbox/Microscopy_Data/ultimatephylo2017.tre')
#plot.phylo(consensus)
#is.rooted(consensus)
#nodelabels()
rRNA18S = reroot(rRNA18S, 76, 0.001)
rRNA18S$tip.label = str_replace_all(rRNA18S$tip.label,'_',' ')
rRNA18S$tip.label[rRNA18S$tip.label=="Cordagalma cordiforme"] = "Cordagalma ordinatum"
rRNA18S$tip.label[rRNA18S$tip.label=="Cordagalma sp "] = "Cordagalma ordinatum"
rRNA18S$tip.label[rRNA18S$tip.label=="Clausophyes ovata"] = "Kephyes ovata"
rRNA18S$tip.label[rRNA18S$tip.label=="Erenna sp"] = "Erenna sirena"
rRNA18S$tip.label[rRNA18S$tip.label=="Apolemia sp"] = "Apolemia lanosa"
rRNA18S$tip.label[rRNA18S$tip.label=="Apolemia sp "] = "Apolemia rubriversa"
rRNA18S$tip.label[rRNA18S$tip.label=="Sphaeronectes gracilis"] = "Sphaeronectes koellikeri"
rRNA18S$tip.label[rRNA18S$tip.label=="Prayidae D27SS7"] = "Desmophyes haematogaster"
rRNA18S$tip.label[rRNA18S$tip.label=="Prayidae D27D2"] = "Craseoa lathetica"
rRNA18S = drop.tip(rRNA18S, 41:42)

transcriptome$tip.label = str_replace_all(transcriptome$tip.label,'_',' ')
transcriptome$tip.label[transcriptome$tip.label=="Cordagalma cordiforme"] = "Cordagalma ordinatum"
transcriptome$tip.label[transcriptome$tip.label=="Cordagalma sp "] = "Cordagalma ordinatum"
transcriptome$tip.label[transcriptome$tip.label=="Clausophyes ovata"] = "Kephyes ovata"
transcriptome$tip.label[transcriptome$tip.label=="Erenna sp"] = "Erenna sirena"
transcriptome$tip.label[transcriptome$tip.label=="Apolemia sp"] = "Apolemia lanosa"
transcriptome$tip.label[transcriptome$tip.label=="Apolemia sp "] = "Apolemia rubriversa"
transcriptome$tip.label[transcriptome$tip.label=="Sphaeronectes gracilis"] = "Sphaeronectes koellikeri"
transcriptome$tip.label[transcriptome$tip.label=="Prayidae D27SS7"] = "Desmophyes haematogaster"
transcriptome$tip.label[transcriptome$tip.label=="Prayidae D27D2"] = "Craseoa lathetica"
transcriptome = drop.tip(transcriptome, 29:37)

multitree = list(rRNA18S,transcriptome)
class(multitree)="multiPhylo"
STree = superTree(multitree)
#plot(STree)
consensus = STree[[2]]
consensus = drop.tip(consensus, "Erenna sirena")
consensus = bind.tip(consensus,tip.label = "Erenna sirena", where=which(consensus$tip.label=="Erenna richardi"), edge.length = 1)
consensus$edge.length[consensus$edge.length==0] = c(0.1,0.1)
consensus$tip.label = str_replace_all(consensus$tip.label,'_',' ')
consensus = reroot(consensus,59,mean(consensus$edge.length))
#plot(consensus)
consensus = drop.tip(consensus, which(consensus$tip.label == "Nectadamas diomedeae" | consensus$tip.label == "Muggiaea atlantica"))
consensus$edge.length[consensus$edge.length<0] <- 0.1
consensus = drop.tip(consensus, 50)
ultratree = chronos(consensus)
plot(ultratree)

##### SIMULATION #####
nmeasured = 10
kind_M = "BM"
nunmeasured = 8
kind_U = "BM"

#Measured Traits
measured = as.data.frame(ultratree$tip.label)

if(kind_M == "BM"){N_BM_m = nmeasured
} else {N_BM_m = round(nmeasured/2)}
if(kind_M == "BM" | kind_M == "mix"){
for(i in 1:N_BM_m){
  VAL = runif(1, min=10, max=100)
  measured = cbind(measured, rTraitCont(ultratree, model="BM", root.value = VAL, sigma = VAL-runif(1, min=1, max=VAL)))
}
}

if(kind_M == "OU"){N_OU_m = nmeasured
} else {N_OU_m = round(nmeasured/2)}
if(kind_M == "OU" | kind_M == "mix"){
for(i in 1:N_OU_m){
  VAL = runif(1, min=10, max=100)
  measured = cbind(measured, rTraitCont(ultratree, model="OU", alpha=runif(1,min=0,max=1), sigma=VAL-runif(1, min=1, max=VAL), root.value = VAL, theta=VAL+runif(1, min=-VAL/2, max=VAL)))
}}
measured = measured[,-1]
names(measured) = paste(rep("MT",ncol(measured)), 1:ncol(measured), sep= "")

#Unmeasured Traits
unmeasured = as.data.frame(ultratree$tip.label)

if(kind_U == "BM"){N_BM_u = nunmeasured
} else {N_BM_u = round(nunmeasured/2)}
if(kind_U == "BM" | kind_U == "mix"){
for(i in 1:N_BM_u){
  VAL = runif(1, min=10, max=100)
  unmeasured = cbind(unmeasured, rTraitCont(ultratree, model="BM", root.value = VAL, sigma=VAL-runif(1, min=1, max=VAL)))
}}

if(kind_U == "OU"){N_OU_u = nunmeasured
} else {N_OU_u = round(nunmeasured/2)}
if(kind_U == "OU" | kind_U == "mix"){
for(i in 1:N_OU_u){
  VAL = runif(1, min=10, max=100)
  unmeasured = cbind(unmeasured, rTraitCont(ultratree, model="OU", alpha=runif(1,min=0,max=1), sigma = VAL-runif(1, min=1, max=VAL), root.value = VAL, theta=VAL+runif(1, min=-VAL/2, max=VAL)))
}}
unmeasured = unmeasured[,-1]
names(unmeasured) = paste(rep("UT",ncol(unmeasured)), 1:ncol(unmeasured), sep= "")

#Habitat
habitat = abs(2^abs(rTraitCont(ultratree, model="BM", root.value = 6, sigma = 3)) + runif(length(ultratree$tip.label),min=0, max=200))
phylosig(ultratree, habitat)

#Prey availability
preyfield = as.data.frame(ultratree$tip.label)
foodtypes = c("Krill", "Copepods", "Ostracods", "Decapods", "Gelata", "Fish", "Chaetognaths")
for(i in 1:length(foodtypes)){
  #preyfield = cbind(preyfield, as.integer(habitat*runif(nrow(preyfield), 0.1, 10))) #uncomment to make it habitat dependet and thus insert residual phylogenetic signal
  preyfield = cbind(preyfield, as.integer(runif(nrow(preyfield), 0.1, 1000)))
}
preyfield = preyfield[,-1]
preyfield = (preyfield/rowSums(preyfield))*100
names(preyfield) = foodtypes
rownames(preyfield) = ultratree$tip.label

#Diet
if(nunmeasured>0){traits = cbind(measured, unmeasured)
} else if(nmeasured<0){traits=unmeasured
} else 
    traits=measured
traits= abs(traits)

selectivity = preyfield
 for(i in 1:length(foodtypes)){
   #traits_i = traits[,sample(1:ncol(traits),2)]
   traits_i = traits[,sample(1:ncol(traits),floor(ncol(traits)/4))]
   if(ncol(traits)<4){traits_i = traits[,sample(1:ncol(traits),1)]}
   for(j in 1:length(ultratree$tip.label))
     #selectivity[j,i] <- (mean(traits_i[,1])/traits_i[j,1])*(mean(traits_i[,2])/traits_i[j,2])
     selectivity[j,i] <- prod(sapply(traits_i, mean)/traits_i[j,])
 }
selectivity = (selectivity/colMeans(selectivity))/10
selectivity = selectivity/rowSums(selectivity)

diet = preyfield*selectivity
diet = 100*(diet/rowSums(diet))
#heatmap(as.matrix(diet))

###### Quality Assesment ######
dietdist <- vegdist(diet, "bray") %>% as.matrix()
selecdist <- vegdist(selectivity, "euc") %>% as.matrix()
traitdist <- vegdist(traits, "euc") %>% as.matrix()
preyfielddist <- vegdist(preyfield,"bray") %>% as.matrix()
mantel(preyfielddist, traitdist)$signif > 0.05
mantel(dietdist, selecdist)$signif < 0.05
mantel(traitdist, selecdist)$signif < 0.05
mantel(dietdist, traitdist)$signif < 0.05
cor.table(cbind(traits,diet))$P[1:ncol(traits),((ncol(traits)+1):(ncol(traits)+length(foodtypes)))] %>% .[.<0.05 & .>0] %>% length()

###### Analysis ######

#Evolution of measured traits
#Test BM, WhiteNoise, OU, EB, Trend, delta, kappa, drift
# for(c in 1:ncol(measured)){
#   C = measured[,c]
#   names(C) = rownames(measured)
#   Ctree = ultratree
#   BMmodel <- fitContinuous(Ctree, C, model="BM")
#   WNmodel <- fitContinuous(Ctree, C, model="white")
#   DRIFTmodel <- fitContinuous(Ctree, C, model="drift")
#   EBmodel <- fitContinuous(Ctree, C, model="EB")
#   OUmodel <- fitContinuous(Ctree, C, model="OU")
#   TRmodel <- fitContinuous(Ctree, C, model="trend")
#   Dmodel <- fitContinuous(Ctree, C, model="delta")
#   Lmodel <- fitContinuous(Ctree, C, model="lambda")
#   Kmodel <- fitContinuous(Ctree, C, model="kappa")
#   print(names(measured)[c])
#   print("WN AICc")
#   print(WNmodel$opt$aicc)
#   print("BM AICc")
#   print(BMmodel$opt$aicc)
#   print("drift AICc")
#   print(DRIFTmodel$opt$aicc)
#   print("EB AICc")
#   print(EBmodel$opt$aicc)
#   print("OU AICc")
#   print(OUmodel$opt$aicc)
#   print("Trend AICc")
#   print(TRmodel$opt$aicc)
#   print("Delta AICc")
#   print(Dmodel$opt$aicc)
#   print("Lambda AICc")
#   print(Lmodel$opt$aicc)
#   print("Kappa AICc")
#   print(Kmodel$opt$aicc)
# }
# }
#THIS STUFF WORKS, ALREADY TESTED. ONLY ODD THING IS THAT OU TRAITS COME OUT AS Delta or Kappa

#Mantel testing
phylodist <- cophenetic(ultratree) %>% as.matrix()
morphdist <- vegdist(measured,"euc") %>% as.matrix()
unmeasureddist <- vegdist(unmeasured,"euc") %>% as.matrix()

mantel(dietdist, morphdist)$signif < 0.05
mantel(dietdist, preyfielddist)$signif < 0.05
mantel(dietdist, phylodist)$signif < 0.05
mantel(morphdist, phylodist)$signif < 0.05

MM_dm = multi.mantel(dietdist, morphdist, nperm=99)
dm_residuals = as.matrix(MM_dm$residuals)
MM_dprey = multi.mantel(dietdist, preyfielddist, nperm=99)
dprey_residuals = as.matrix(MM_dprey$residuals)
MM_dp = multi.mantel(dietdist, phylodist, nperm=99)
dp_residuals = as.matrix(MM_dp$residuals)

       #Unmeasured residuals
UM_Residuals = as.matrix(multi.mantel(dm_residuals, preyfielddist, nperm=99)$residuals)
mantel(UM_Residuals, unmeasureddist)$signif < 0.05
mantel(UM_Residuals, phylodist)$signif < 0.05

#Figuring out number of traits involved
cbind(sqrt(as.numeric(phylodist)), as.numeric(morphdist)) %>% as.data.frame() -> PDEDdists
names(PDEDdists) = c("SqrtPhylogeneticDistance","FunctionalDistance")
PDEDdists %>% lm(SqrtPhylogeneticDistance~FunctionalDistance, data=.) -> LMPDED
ggplot(PDEDdists,aes(x=SqrtPhylogeneticDistance,y=FunctionalDistance))+geom_density2d()+geom_smooth()+geom_jitter(width=0.05,height = 0.05)

cbind(sqrt(as.numeric(phylodist)), as.numeric(dietdist)) %>% as.data.frame() -> PDEDdists
names(PDEDdists) = c("SqrtPhylogeneticDistance","EcologicalDistance")
PDEDdists %>% lm(SqrtPhylogeneticDistance~EcologicalDistance, data=.) -> LMPDED
ggplot(PDEDdists,aes(x=SqrtPhylogeneticDistance,y=EcologicalDistance))+geom_density2d()+geom_smooth()+geom_jitter(width=0.05,height = 0.05)

       #Principal components
UM_PC1 = prcomp(t(unmeasured))$rotation[,1]
UMR_PC1 = prcomp(t(UM_Residuals))$rotation[,1]
phylosig(ultratree, UMR_PC1)
startree = rescale(ultratree, model = "lambda", 0)
StarModel <- fitContinuous(startree, UMR_PC1, model="BM")
BMmodel <- fitContinuous(ultratree, UMR_PC1, model="BM")
WNmodel <- fitContinuous(ultratree, UMR_PC1, model="white")
DRIFTmodel <- fitContinuous(ultratree, UMR_PC1, model="drift")
EBmodel <- fitContinuous(ultratree, UMR_PC1, model="EB")
OUmodel <- fitContinuous(ultratree, UMR_PC1, model="OU")
TRmodel <- fitContinuous(ultratree, UMR_PC1, model="trend")
Dmodel <- fitContinuous(ultratree, UMR_PC1, model="delta")
Lmodel <- fitContinuous(ultratree, UMR_PC1, model="lambda")
Kmodel <- fitContinuous(ultratree, UMR_PC1, model="kappa")
print("WN AICc")
WNmodel$opt$aicc
print("Star AICc")
StarModel$opt$aicc
print("BM AICc")
print(BMmodel$opt$aicc)
print("drift AICc")
print(DRIFTmodel$opt$aicc)
print("EB AICc")
print(EBmodel$opt$aicc)
print("OU AICc")
print(OUmodel$opt$aicc)
print("Trend AICc")
print(TRmodel$opt$aicc)
print("Delta AICc")
print(Dmodel$opt$aicc)
print("Lambda AICc")
print(Lmodel$opt$aicc)
print("Kappa AICc")
print(Kmodel$opt$aicc)

diet_PC1 = prcomp(t(diet))$rotation[,1]
AllMeasured_LM = lm(diet_PC1~measured[,1] + measured[,2] + measured[,3] + PF_PC1)

AMLM_res = AllMeasured_LM$residuals
physignal(AMLM_res, ultratree)


PF_PC1 = prcomp(t(preyfield))$rotation[,1]

phy_PC1 = prcomp(t(phylodist))$rotation[,1]
morph_PC1 = prcomp(t(measured))$rotation[,1]

plot(phy_PC1, morph_PC1)
text(phy_PC1, morph_PC1, names(phy_PC1), cex=0.4, col="blue")
plot(as.numeric(sqrt(phylodist[upper.tri(phylodist)])), as.numeric(morphdist[upper.tri(morphdist)]))


FPDists = as.matrix(vegdist(cbind(phy_PC1, morph_PC1), "euc"))
heatmap(FPDists, symm=T)

