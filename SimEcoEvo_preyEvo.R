#SIMULATING CAUSAL PATHWAYS FOR DIET DETERMINATION

rm(list=ls())
par(ask=F)
#Load libraries
library(tidyverse)
library(vegan)
library(ape)
library(phytools)
library(picante)
library(geiger)
library(phylobase)
library(geomorph)

##### SIMULATION #####
ultratree <- rtree(50) %>% chronos()
ultratree$tip.label <- paste("SP", 1:length(ultratree$tip.label), sep="")

nmeasured = 5
kind_M = "BM"
nunmeasured = 5
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

if(nunmeasured>0){traits = cbind(measured, unmeasured)
} else if(nmeasured<0){traits=unmeasured
} else 
  traits=measured
traits= abs(traits)

#Reference traits
refTs <- cbind(rTraitCont(ultratree, model="OU"), rTraitCont(ultratree, model="BM")) %>% as.data.frame()
names(refTs) = c("OU", "BM")
BMrefDist = vegdist(refTs$BM,"euc") %>% as.matrix()
OUrefDist = vegdist(refTs$OU,"euc") %>% as.matrix()
colnames(OUrefDist) = ultratree$tip.label
colnames(BMrefDist) = ultratree$tip.label
rownames(OUrefDist) = ultratree$tip.label
rownames(BMrefDist) = ultratree$tip.label

#Habitat
habitat = abs(abs(rTraitCont(ultratree, model="BM", root.value = 200, sigma = 200)) + runif(length(ultratree$tip.label),min=200, max=200))
phylosig(ultratree, habitat)
#barplot(habitat[order(habitat, decreasing = T)])
habitat_zones = habitat
habitat_zones[habitat < quantile(habitat)[5] | habitat == quantile(habitat)[5]] <- 4
habitat_zones[habitat < quantile(habitat)[4]] <- 3
habitat_zones[habitat < quantile(habitat)[3]] <- 2
habitat_zones[habitat < quantile(habitat)[2]] <- 1
phylosignal(habitat_zones,ultratree)

#Prey evolution
preytypes = paste(rep("P",10),1:10,sep="")
preytree = chronos(rtree(length(preytypes)))
preytree$tip.label <- preytypes
preytraits = as.data.frame(preytree$tip.label)
for(i in 1:ncol(traits)){
    VAL = mean(traits[,i])
    preytraits = cbind(preytraits, rTraitCont(preytree, model="BM", root.value = VAL, sigma=VAL/4))
  }
preytraits = preytraits[,-1]
names(preytraits) = paste(rep("PT",ncol(preytraits)), 1:ncol(preytraits), sep= "")

refPT = rTraitCont(preytree, model="BM")

#Prey availability
preyfields = list()
for(i in 1:max(habitat_zones)){
  preyfields[[i]] <- as.integer(runif(length(preytypes), 0.1, 1000))
  names(preyfields[[i]]) = preytypes
}
preyfield = preytypes
for(i in 1:length(habitat_zones)){
  preyfield = rbind(preyfield, preyfields[[habitat_zones[i]]]) 
}
preyfield = as.data.frame(apply(preyfield[-1,], 2, as.numeric))
rownames(preyfield) = ultratree$tip.label
preyfield = (preyfield/rowSums(preyfield))*100
physignal(preyfield, ultratree)

#Diet
selectivity = preyfield
 for(i in 1:length(preytypes)){
   for(j in 1:length(ultratree$tip.label))
     SEL_VAL <- abs(traits[j,]-preytraits[i,])
     SEL_VAL = SEL_VAL[1,]/rowMeans(preytraits[i,])
     selectivity[j,i] <- 1/sum(t(SEL_VAL)[,1])[1]
 }
selectivity = (selectivity/colMeans(selectivity))/10
selectivity = selectivity/rowSums(selectivity)

diet = preyfield*selectivity
diet = 100*(diet/rowSums(diet))
#heatmap(as.matrix(diet))

#DISTANCE MATRICES
phylodist <- cophenetic(ultratree) %>% as.matrix()
morphdist <- vegdist(measured,"euc") %>% as.matrix()
unmeasureddist <- vegdist(unmeasured,"euc") %>% as.matrix()
traitdist <- vegdist(traits, "euc") %>% as.matrix()
dietdist <- vegdist(diet, "bray") %>% as.matrix()
selecdist <- vegdist(selectivity, "euc") %>% as.matrix()
preyfielddist <- vegdist(preyfield,"bray") %>% as.matrix()
REV_dietdist <- vegdist(t(diet),"bray") %>% as.matrix()
REV_preyfielddist <- vegdist(t(preyfield),"bray") %>% as.matrix()
REV_selectivitydist <- vegdist(t(selectivity),"euc") %>% as.matrix()
REV_preytraitdist <- vegdist(preytraits,"euc") %>% as.matrix()
preyphylodist <- cophenetic(preytree) %>% as.matrix()

###### Quality Assesment ######
  #These should all be TRUE
mantel.partial(preyfielddist, traitdist, phylodist)$signif > 0.05
mantel(dietdist, selecdist)$signif < 0.05
mantel(traitdist, selecdist)$signif < 0.05
mantel(dietdist, traitdist)$signif < 0.05
mantel(preyphylodist, REV_preytraitdist)$signif < 0.05
mantel(preyphylodist, REV_selectivitydist)$signif < 0.05
mantel(REV_preytraitdist, REV_selectivitydist)$signif < 0.05
mantel(REV_preytraitdist, REV_dietdist)$signif < 0.05

varpart(t(selectivity), preytraits, refPT)$part$frac$Adj.R.squared[2]*100
varpart(selectivity, traits, refTs$BM) $part$frac$Adj.R.squared[1]*100

  #There should be more than 10 correlations between traits and prey types ingested
cor.table(cbind(traits,diet))$P[1:ncol(traits),((ncol(traits)+1):(ncol(traits)+length(preytypes)))] %>% .[.<0.05 & .>0] %>% length()

  #EXPECTED VARIANCE DECOMPOSITION:
paste("Phylogeny in Diet:   ", paste(floor((varpart(diet, refTs$BM, refTs$OU)$part$fract$Adj.R.square[1])*100)), "%", sep="")
paste("Habitat in Preyfield:   ", paste(floor((varpart(preyfield, habitat, refTs$OU)$part$fract$Adj.R.square[1])*100)), "%", sep="")
paste("Phylogeny in measured:   ", paste(floor((varpart(measured, refTs$BM, refTs$OU)$part$fract$Adj.R.square[1])*100)), "%", sep="")
paste("Phylogeny in unmeasured:   ", paste(floor((varpart(unmeasured, refTs$BM, refTs$OU)$part$fract$Adj.R.square[1])*100)), "%", sep="")
paste("Phylogeny in habitat:   ", paste(floor((varpart(habitat, refTs$BM, refTs$OU)$part$fract$Adj.R.square[1])*100)), "%", sep="")
VP_Sel = varpart(selectivity, measured, unmeasured)$part$indfract$Adj.R.square
paste("Measured Selectivity:   ", paste(floor((VP_Sel[1])*100)), "%", sep="")
paste("Uneasured Selectivity:   ", paste(floor((VP_Sel[3])*100)), "%", sep="")
paste("Stochastic Selectivity:   ", paste(floor((VP_Sel[4])*100)), "%", sep="")
VP_ALL = varpart(diet, X = preyfield, habitat, selectivity)$part
paste("   Preyfield:", paste(abs(floor((VP_ALL$indfrac$Adj.R.square)*100)[1]), "%", sep=""), "   Habitat:", paste(abs(floor((VP_ALL$fract$Adj.R.square)*100)[2]), "%", sep=""), "   Selectivity:", paste(abs(floor(VP_ALL$contr1$Adj.R.square[5]*100)), "%", sep=""), "   Non-Evolutionary Residuals:", paste(floor(VP_ALL$indfrac$Adj.R.square[8]*100), "%", sep=""), sep=" ")

###### Analysis ######

#Evolution of measured traits
#Test BM, WhiteNoise, OU, EB, Trend, delta, kappa, drift
# for(c in 1:ncol(measured)){
#   C = measured[,c]
#   names(C) = rownames(measured)
#   ultratree = ultratree
#   BMmodel <- fitContinuous(ultratree, C, model="BM")
#   WNmodel <- fitContinuous(ultratree, C, model="white")
#   DRIFTmodel <- fitContinuous(ultratree, C, model="drift")
#   EBmodel <- fitContinuous(ultratree, C, model="EB")
#   OUmodel <- fitContinuous(ultratree, C, model="OU")
#   TRmodel <- fitContinuous(ultratree, C, model="trend")
#   Dmodel <- fitContinuous(ultratree, C, model="delta")
#   Lmodel <- fitContinuous(ultratree, C, model="lambda")
#   Kmodel <- fitContinuous(ultratree, C, model="kappa")
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
mantel(dietdist, morphdist)$signif < 0.05
mantel(dietdist, preyfielddist)$signif < 0.05
mantel(dietdist, phylodist)$signif < 0.05
mantel(morphdist, phylodist)$signif < 0.05
mantel.partial(dietdist, phylodist, morphdist)$signif 
mantel.partial(dietdist, morphdist, phylodist)$signif 
mantel.partial(dietdist, preyfielddist, phylodist)$signif < 0.05
mantel.partial(dietdist, phylodist, preyfielddist)$signif < 0.05
mantel.partial(dietdist, preyfielddist, morphdist)$signif < 0.05
mantel.partial(dietdist, morphdist, preyfielddist)$signif < 0.05
DM_P = rda(diet, measured, preyfield)
DM_P$CA
DP_M = rda(diet, preyfield, measured)

VP_DPM = varpart(Y = as.dist(dietdist), X =  preyfield, measured)
#Variance explained by:
#Preyfield   |   Measured Morphology
paste(as.integer(VP_DPM$part$fract$Adj.R.square[1:2]*100), c("%","%"), sep="")
#Residuals left:
paste(as.integer(VP_DPM$part$indfract$Adj.R.square[4]*100), "%", sep="")

VP_DHM = varpart(Y = as.dist(dietdist), X =  habitat_zones, measured)
#Variance explained by:
#Habitat   |   Measured Morphology
paste(as.integer(VP_DHM$part$fract$Adj.R.square[1:2]*100), c("%","%"), sep="")
#Residuals left:
paste(as.integer(VP_DHM$part$indfract$Adj.R.square[4]*100), "%", sep="")

VP_DHdM = varpart(Y = as.dist(dietdist), X =  habitat, measured)
#Variance explained by:
#Habitat Depth   |   Measured Morphology
paste(as.integer(VP_DHdM$part$fract$Adj.R.square[1:2]*100), c("%","%"), sep="")
#Residuals left:
paste(as.integer(VP_DHdM$part$indfract$Adj.R.square[4]*100), "%", sep="")

MM_dm = multi.mantel(dietdist, morphdist, nperm=99)
dm_residuals = as.matrix(MM_dm$residuals)
MM_dprey = multi.mantel(dietdist, preyfielddist, nperm=99)
dprey_residuals = as.matrix(MM_dprey$residuals)
MM_dp = multi.mantel(dietdist, phylodist, nperm=99)
dp_residuals = as.matrix(MM_dp$residuals)
MM_mp = multi.mantel(morphdist, phylodist, nperm=99)
mp_residuals = as.matrix(MM_mp$residuals)

       #Unmeasured residuals
UM_Residuals = as.matrix(multi.mantel(dm_residuals, preyfielddist, nperm=99)$residuals)
mantel(UM_Residuals, unmeasureddist)$signif < 0.05
mantel(UM_Residuals, phylodist)$signif < 0.05
UM_ran_Residuals = as.matrix(multi.mantel(UM_Residuals, phylodist, nperm=99)$residuals)
UM_evo_Residuals = as.matrix(multi.mantel(UM_Residuals, phylodist, nperm=99)$fitted.values)

Phy_PCOA <- pcoa(phylodist)
PCOAPCAPhy <- prcomp(t(Phy_PCOA$vectors))$rotation[,1:8]

#NON-Distance based approach
VP_DPM_NonD = varpart(diet, X = preyfield, measured)
VP_DPM_NonD$part %>% str()

#Figuring out number of traits involved
cbind(sqrt(as.numeric(phylodist)), as.numeric(morphdist)) %>% as.data.frame() -> PDEDdists
names(PDEDdists) = c("SqrtPhylogeneticDistance","FunctionalDistance")
PDEDdists %>% lm(SqrtPhylogeneticDistance~FunctionalDistance, data=.) -> LMPDED
ggplot(PDEDdists,aes(x=SqrtPhylogeneticDistance,y=FunctionalDistance))+geom_density2d()+geom_smooth()+geom_jitter(width=0.05,height = 0.05)

cbind(sqrt(as.numeric(phylodist)), as.numeric(dietdist)) %>% as.data.frame() -> PDEDdists
names(PDEDdists) = c("SqrtPhylogeneticDistance","EcologicalDistance")
PDEDdists %>% lm(SqrtPhylogeneticDistance~EcologicalDistance, data=.) -> LMPDED
ggplot(PDEDdists,aes(x=SqrtPhylogeneticDistance,y=EcologicalDistance))+geom_density2d()+geom_smooth()+geom_jitter(width=0.05,height = 0.05)

    #Principal Coordinate Analysis
diet_PCOA <- pcoa(dietdist)
diet_PCOA$values$Broken_stick %>% barplot(.[order(.)])
PF_PCOA <- pcoa(preyfielddist)
PF_PCOA$values$Broken_stick %>% barplot(.[order(.)])
Phy_PCOA$values$Broken_stick %>% barplot(.[order(.)])
UMR_PCOA <- pcoa(UM_Residuals)
UMR_PCOA$values$Broken_stick %>% barplot(.[order(.)])
plot(cmdscale(UM_Residuals), cmdscale(phylodist))
Morph_PCOA <- pcoa(morphdist)
Morph_PCOA$values$Broken_stick %>% barplot(.[order(.)])
nmeasured

  #Identify the evolution behind morphological traits from the PCos of morphological dissimilarity square matrix WORKS
# kind_M
# startree = rescale(ultratree, model = "lambda", 0)
# for(c in 1:ncol(Morph_PCOA$vectors)){
#   C = Morph_PCOA$vectors[,c]
#   names(C) = ultratree$tip.label
#   StarModel <- fitContinuous(ultratree, C, model="BM")
#   BMmodel <- fitContinuous(ultratree, C, model="BM")
#   WNmodel <- fitContinuous(ultratree, C, model="white")
#   DRIFTmodel <- fitContinuous(ultratree, C, model="drift")
#   EBmodel <- fitContinuous(ultratree, C, model="EB")
#   OUmodel <- fitContinuous(ultratree, C, model="OU")
#   TRmodel <- fitContinuous(ultratree, C, model="trend")
#   Dmodel <- fitContinuous(ultratree, C, model="delta")
#   Lmodel <- fitContinuous(ultratree, C, model="lambda")
#   Kmodel <- fitContinuous(ultratree, C, model="kappa")
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

       #Principal components
kind_U
UMR_PC1 = prcomp(t(UM_Residuals))$rotation[,1]
phylosig(ultratree, UMR_PC1)
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
PF_PC1 = prcomp(t(preyfield))$rotation[,1]
phy_PC1 = prcomp(t(phylodist))$rotation[,1]
morph_PC1 = prcomp(t(measured))$rotation[,1]

phylosig(ultratree, lm(lm(diet_PC1 ~ PF_PC1)$residuals ~ morph_PC1)$residuals)

plot(phy_PC1, morph_PC1)
text(phy_PC1, morph_PC1, names(phy_PC1), cex=0.4, col="blue")
plot(as.numeric(sqrt(phylodist[upper.tri(phylodist)])), as.numeric(morphdist[upper.tri(morphdist)]))


FPDists = as.matrix(vegdist(cbind(phy_PC1, morph_PC1), "euc"))
heatmap(FPDists, symm=T)

