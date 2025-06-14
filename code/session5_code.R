
## Session 5 -  Assigning Individuals to Populations

# follow the links (ctrl click) to find the accompanying session material
# https://green-striped-gecko.github.io/iccb2025/session06.html


## ------------------------------------------------------------------------------------------------------------------------
library(dartRverse)



## ----------------------------------------------------------------------------------------------------------------------------------------------------
# We will first set the verbosity globally to level 3
gl.set.verbosity(3)



## ----------------------------------------------------------------------------------------------------------------------------------------------------
# Read in the data set for the worked example
gl <- readRDS("./data/assignment_example1.Rdata")
# Familiarize yourself with its contents
gl  
nLoc(gl)
nInd(gl)
nPop(gl)
# Display a list of populations and sample sizes
table(pop(gl))


## Analysis 1: Assignment by genotype likelihood ---------------------------------------------------------------------------------------------------------------------
gen.result<-gl.assign.on.genotype(gl, unknown="AA011731", nmin=10)



## Analysis 2: Assignment by Private Alleles -------------------------------------------------------------------------------------------------------------------
pa.result <- gl.assign.pa(gl, unknown="AA011731", nmin=10, alpha=0.05)



## Analysis 3: Assignment by PCA --------------------------------------------------------------------------------------------------------------------
pca_pa_result <-gl.assign.pca(pa.result, unknown="AA011731")


## Analysis 4: Assignment by Mahalanobis Distances --------------------------------------------------------------------------------------------------------------------------------
mahal_result <- gl.assign.mahal(pa.result,unknown="AA011731")



# Exercise ----------------------------------------------------------------------------------------------------------------------------------
# The data
gl
# The unknown
Unknown = "AA046092"
# Preliminaries
popNames(gl)

gl2 <- gl.keep.pop(gl, pop.list=c("EmsubBamuAli", "EmsubFlyGuka", "EmsubFlyJikw",
                                  "EmsubJardine", "EmsubKerema", "EmsubMorehead"))

# Knock yourself out


