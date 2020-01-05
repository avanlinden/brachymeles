###### Brachymeles PGLS body elongation analysis 
###### Abby Vander Linden 
###### 8/7/18, edited 1/11/19, edited again 11/6/19, 
###### cleaned up and took the swears out before sharing on 1/5/20

### Load libraries =======================

library(phytools)
library(ape)
library(geiger)
library(caper)
library(tidyverse)
library(phangorn)


### Load tree(s) ===========================

#load nexus file with ML and Bayesian trees
silerPhy <- read.nexus(file = "Siler_Brachymeles_Evolution_2011.nex")

#inspect trees

plot(ladderize(silerPhy$ML), cex = 0.6)
plot(silerPhy$BI, cex = 0.6)

# We'll work with the maximum likelihood tree for our analysis


### Load body shape data =============================

skinkData <- read.csv(file = "skink_measurements.csv")

#format as tibble
skinkDat <- as_tibble(skinkData)

#check
skinkDat

#reformat species names to match tree tip labels with underscore

skinkDat <- skinkDat %>% 
  mutate(specNum = as.character(specNum)) %>% 
  mutate(species = as.character(species)) %>% 
  mutate(species = str_replace(species, " ", "_"))

# write cleaned data to file if you want
write.csv(skinkDat, file = "cleaned_Brachymeles_data.csv")

#create a tibble with the species names as they appear in the measurement data (we will check this against tree names later)
datSpec <- skinkDat %>% 
  dplyr::select(species)


### Adjusting tree tips ==================

#make a vector of the tree tips
evoTips <- silerPhy$ML$tip.label

# make that vector into a tibble for easy string manipulation, I'm not here to do stuff with strings in base R 
# (it will get mad that you made a vector into a tibble but it will still do it; if it ever stops working, replace "as_tibble" with "enframe")
treeSpec <- evoTips %>% 
  as_tibble() %>% 
  dplyr::rename(species = value) 

#remove the spec numbers from the tip labels so it's just species_names
treeSpecEdit <- treeSpec %>% 
  separate(species, c("genus", "spec", "subsp", "specNum"), "_") %>% #warning messages here are ok -- some tips have a subspecies and some don't -- treeSpecEdit tibble should be created anyway 
  unite(gen_spec, genus, spec, sep = "_") %>% 
  mutate(specNum_suffix = ifelse(is.na(specNum), subsp, specNum)) %>% 
  dplyr::select(species = gen_spec, specNum = specNum_suffix)
 
#this should now have the genus_species name in the same format as the data, plus the species number as a separate column 
treeSpecEdit

#append these nicely formatted species names to the tip labels of the ML tree
skinkTree <- silerPhy$ML
skinkTree$tip.label <- treeSpecEdit$species

#now the tip labels are in our format
plot(skinkTree, cex = 0.5)

#CHECK -- are all of our sample skinks in the tree?
treeSpecEdit %>% 
  filter(species %in% c(datSpec$species)) %>% 
  dplyr::select(species) %>% 
  unique() #ok so for whatever reason only 8 are in the tree (spelling, synonym issues)

#there is no B. vermis because it's not in the tree -- no tissue samples exist and restricted access to island; we have to drop B. vermis from our sample data because we have no phylo info for it

## start editing some tree tip names to correspond to our specimen names (resolve synonyms, misspellings, etc)

# vectors of species with names we need to change using specimen numbers
samad <- c("KU310825", "KU311216", "KU311220")
cobos <- c("KU324020", "KU324022")
ilocandia <- c("KU308004", "KU307967")
boulengeri_boulengeri <- c("KU323409", "KU307752", "KU307753")

# create a tibble with a new column that we will edit to have the correct species names
treeSpecRev <- treeSpecEdit %>% 
  mutate(specRevised = species)

# create a new tibble with the "old" species names for samad, cobos, etc replaced with our versions
treeSpecOk <- treeSpecRev %>%
  filter(specNum %in% samad | specNum %in% cobos | specNum %in% ilocandia | specNum %in% boulengeri_boulengeri) %>% 
  mutate(specRevised = str_replace(specRevised, "bonitae", "ilocandia")) %>% 
  mutate(specRevised = str_replace(specRevised, "samarensis", "bicolandia")) %>% 
  mutate(specRevised = str_replace(specRevised, "gracilis", "samad")) %>% 
  mutate(specRevised = str_replace(specRevised, "boulengeri", "boulengerib")) %>% 
  right_join(treeSpecEdit)

#check all of our edited tree tips -- if they didn't need changed, specRevised should have NA, if not they should have the new genus_species name
print(treeSpecOk, n = Inf)

# grab a tibble of the correct species names for all tips
treeTipsRev <- treeSpecOk %>% 
  mutate(specRevised = ifelse(is.na(specRevised), species, specRevised)) %>%
  dplyr::select(specRevised)

#yes creating all these new data objects is clunky, I know, it's old code, I'm sorry, it still works

#check concordance between tree tip labels and species labels in our dataset
treeTipsRev %>% 
  filter(specRevised %in% c(datSpec$species)) %>% 
  unique() ## aight, we got everyone except vermis, we good

#rename the tree tip labels to our revised versions
skinkTree$tip.label <- treeTipsRev$specRevised

#prune the tree to just the skinks in our dataset
dropTips <- treeTipsRev %>% 
  filter(!specRevised %in% c(datSpec$species))

skinkTreePruned <- drop.tip(skinkTree, dropTips$specRevised)

tree <- skinkTreePruned

#check our gorgeous pruned tree
plot(skinkTreePruned, cex = 0.6)
plot(tree)

### Combining individual tips =====================

# for PGLS, we need to combine tips of the same species into one branch

# I used Liam Revell's method from his phytools blog

# define some variables: 

taxa <- as.vector(unique(skinkTreePruned$tip.label))
nodes<-1:tree$Nnode+Ntip(tree)
subtrees<-list()

for(i in 1:tree$Nnode) subtrees[[i]]<-extract.clade(tree,nodes[i])

names(subtrees)<-nodes ## all subtrees

all.foos <- c()

all.foos<-nodes[sapply(subtrees,function(x) all(x$tip.label %in% taxa))] #this is actully just all the nodes

foo.mrcas<-monoNodes[sapply(monoNodes,function(x,tree,y)
  !Ancestors(tree,x,"parent")%in%y,
  tree=tree,y=monoNodes)]

w.foos<-which(tree$tip.label %in% taxa)

#plot tree with node labels so I can note down the nodes of each species (nope, this is not fancy)
plot(tree, cex = 0.6)
nodelabels(cex = 0.5)

#nodes of species monophyly (by visual inspection aka "looking" with my "eyes"):

monoNodes <- c(42, 43, 46, 47, 51, 54, 57, 61, 66, 68, 71, 74)

## rename tips uniquely
tree$tip.label[w.foos] <- paste(tree$tip.label[w.foos],
                                replicate(length(w.foos), paste(sample(letters, 6),
                                                                collapse =
                                                                  "")), sep = "_")
collapsed <- tree

## iterate over all MRCAs of foo clades
for(i in 1:length(foo.mrcas)){
  M<-matchNodes(tree,collapsed)
  nn<-M[which(M[,1]==foo.mrcas[i]),2]
  dd<-Descendants(collapsed,nn)[[1]]
  h<-sapply(dd,nodeheight,tree=collapsed)
  collapsed$tip.label[dd[1]]<-paste(foo.mrcas[i])
  ind<-which(collapsed$edge[,2]==dd[1])
  collapsed$edge.length[ind]<-collapsed$edge.length[ind]+
    mean(h)-h[1]
  if(length(dd)>1)
    collapsed<-drop.tip(collapsed,
                        collapsed$tip.label[dd[2:length(dd)]])
}

plot(collapsed, cex = 0.6) #it worked!

#put the species names back as the tip labels

newTree <- collapsed
newTree$tip.label <- c("Brachymeles_tridactylus", "Brachymeles_ilocandia", "Brachymeles_bicolandia", "Brachymeles_boulengerib", "Brachymeles_talinis", "Brachymeles_kadwa", "Brachymeles_muntingkamay", "Brachymeles_bicolor", "Brachymeles_samad", "Brachymeles_pathfinderi", "Brachymeles_makusog", "Brachymeles_minimus")

plot(newTree, cex = 0.6)

### This tree is THE TREE that is ready to use! ================

theTree <- newTree
plot(theTree)

# write nexus file of this tree and save for later
write.tree(theTree, file = "pruned_renamed_Brachymeles.nex")

### PGLS ======================================

### To re-perform the PGLS analyses, run the following code.
### Make sure you have loaded the libraries at the top of this file, especially tidyverse and caper. 

#If needed, load the CORRECT tree file:
theTree <- read.tree(file = "pruned_renamed_Brachymeles.nex")

# Load the CLEANED and CORRECT data if needed:
skinkDat <- as_tibble(read.csv("cleaned_Brachymeles_data.csv"))

#tidy data and reformat species names if you reloaded
skinkData <- skinkDat %>% 
  dplyr::select(-X) %>% 
  mutate(species = str_replace(species, " ", "_")) %>% 
  as.data.frame()

#create the comparative data object
#there will be a message about data being dropped -- this is because B. vermis is not in our tree
skinkComp <- comparative.data(phy = theTree, data = skinkData, names.col = species, vcv = TRUE, warn.dropped = TRUE)

### PGLS models for area vs body size ==============

head <- pgls(log10(headArea) ~ log10(SVLDig), data = skinkComp, lambda = "ML")
summary(head)
plot(pgls.profile(head, which = "lambda")) #this is a plot of the likelihood surface of the estimates for lambda, it's not important for the outcome of the model but it can be useful to check if you want to see how the model is behaving

mid <- pgls(log10(midDorsalArea) ~ log10(SVLDig), data = skinkComp, lambda = "ML")
summary(mid)

tail <- pgls(log10(midTailArea) ~ log10(SVLDig), data = skinkComp, lambda = "ML")
summary(tail)

# store all pgls area results in a table

row_vars = c("Int_est","Int_p", "slope_est", "slope_p", "lambda_est", "lambda_ci_low", "lambda_ci_high", "adj_rsq")

headAparam <- c(summary(head)$coeff[1,1], 
                summary(head)$coeff[1,4],
                summary(head)$coeff[2,1],
                summary(head)$coeff[2,4],
                head$param.CI$lambda[1]$opt,
                head$param.CI$lambda$ci.val[1],
                head$param.CI$lambda$ci.val[2],
                summary(head)$adj.r.squared
)

midAparam <- c(summary(mid)$coeff[1,1], 
               summary(mid)$coeff[1,4],
               summary(mid)$coeff[2,1],
               summary(mid)$coeff[2,4],
               mid$param.CI$lambda[1]$opt,
               mid$param.CI$lambda$ci.val[1],
               mid$param.CI$lambda$ci.val[2],
               summary(mid)$adj.r.squared
)

tailAparam <- c(summary(tail)$coeff[1,1], 
                summary(tail)$coeff[1,4],
                summary(tail)$coeff[2,1],
                summary(tail)$coeff[2,4],
                tail$param.CI$lambda[1]$opt,
                tail$param.CI$lambda$ci.val[1],
                tail$param.CI$lambda$ci.val[2],
                summary(tail)$adj.r.squared
)

pglsAreaResults <- tibble(row_vars, headAparam, midAparam, tailAparam)

write.csv(pglsAreaResults, file = "skink_area_pgls_results.csv")


### PGLS models for CIRCUMFERENCE vs SVL ==================

headC <- pgls(log10(headCirc) ~ log10(SVLDig), data = skinkComp, lambda = "ML")
summary(headC)
plot(pgls.profile(headC, which = "lambda")) #this is a plot of the likelihood surface of the estimates for lambda, it's not important for the outcome of the model

midC <- pgls(log10(midDorsalCirc) ~ log10(SVLDig), data = skinkComp, lambda = "ML")
summary(midC)

tailC <- pgls(log10(midTailCirc) ~ log10(SVLDig), data = skinkComp, lambda = "ML")
summary(tailC) 

#gather and store pgls results: intercept estimate and p-value, slope (SVL) estimate and p-value, lambda estimate, lambda 95% CI intervals, adjusted R-squared)

row_vars = c("Int_est","Int_p", "slope_est", "slope_p", "lambda_est", "lambda_ci_low", "lambda_ci_high", "adj_rsq")


headCparam <- c(summary(headC)$coeff[1,1], 
                summary(headC)$coeff[1,4],
                summary(headC)$coeff[2,1],
                summary(headC)$coeff[2,4],
                headC$param.CI$lambda[1]$opt,
                headC$param.CI$lambda$ci.val[1],
                headC$param.CI$lambda$ci.val[2],
                summary(headC)$adj.r.squared
                )

midCparam <- c(summary(midC)$coeff[1,1], 
                summary(midC)$coeff[1,4],
                summary(midC)$coeff[2,1],
                summary(midC)$coeff[2,4],
               midC$param.CI$lambda[1]$opt,
               midC$param.CI$lambda$ci.val[1],
                midC$param.CI$lambda$ci.val[2],
                summary(midC)$adj.r.squared
)

tailCparam <- c(summary(tailC)$coeff[1,1], 
               summary(tailC)$coeff[1,4],
               summary(tailC)$coeff[2,1],
               summary(tailC)$coeff[2,4],
               tailC$param.CI$lambda[1]$opt,
               tailC$param.CI$lambda$ci.val[1],
               tailC$param.CI$lambda$ci.val[2],
               summary(tailC)$adj.r.squared
)

pglsCircResults <- tibble(row_vars, headCparam, midCparam, tailCparam)

write.csv(pglsCircResults, file = "skink_circumference_pgls_results.csv")


### Phylomorphospaces for Area ======================

#use this code to make phylomorphospace plots of variable relationships

skinkDatPrune <- skinkData[-4,] #drop vermis

#head area and svl
headMat <- cbind(log10(skinkDatPrune$SVLDig), log10(skinkDatPrune$headArea))
rownames(headMat) <- skinkDatPrune$species
colnames(headMat) <- c("SVL", "headArea")

phylomorphospace(theTree, headMat, label = "horizontal", fsize = 0.7)

#middorsal area and svl
midMat <- cbind(log10(skinkDatPrune$SVLDig), log10(skinkDatPrune$midDorsalArea))
rownames(midMat) <- skinkDatPrune$species
colnames(midMat) <- c("SVL", "midDorsalArea")

phylomorphospace(theTree, midMat, label = "horizontal", fsize = 0.7)

#tail area and svl
tailMat <- cbind(log10(skinkDatPrune$SVLDig), log10(skinkDatPrune$midTailArea))
rownames(tailMat) <- skinkDatPrune$species
colnames(tailMat) <- c("SVL", "headArea")

phylomorphospace(theTree, tailMat, label = "horizontal", fsize = 0.7)

### Phylomorphospaces for Circumference ======================

#head circ and svl
headCMat <- cbind(log10(skinkDatPrune$SVLDig), log10(skinkDatPrune$headCirc))
rownames(headCMat) <- skinkDatPrune$species
colnames(headCMat) <- c("SVL", "headCirc")

phylomorphospace(theTree, headCMat, label = "horizontal", fsize = 0.7)

#middorsal circ and svl
midCMat <- cbind(log10(skinkDatPrune$SVLDig), log10(skinkDatPrune$midDorsalCirc))
rownames(midCMat) <- skinkDatPrune$species
colnames(midCMat) <- c("SVL", "midDorsalCirc")

phylomorphospace(theTree, midCMat, label = "horizontal", fsize = 0.7)

#tail circ and svl
tailCMat <- cbind(log10(skinkDatPrune$SVLDig), log10(skinkDatPrune$midtailirc))
rownames(tailCMat) <- skinkDatPrune$species
colnames(tailCMat) <- c("SVL", "tailcirc")

phylomorphospace(theTree, tailCMat, label = "horizontal", fsize = 0.7)


