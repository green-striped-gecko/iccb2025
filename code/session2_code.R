## Session 2 - Estimating Effective Population Size and Key Stats

# follow the links (ctrl click) to find the accompanying session material
# https://green-striped-gecko.github.io/iccb2025/session02.html

## ----------------------------------------------------------------------------------
library(dartRverse)


## ---------------------------------------------------------------------------------------------------------------
gls <- possums.gl[c(1:5,31:35),1:7]  #small data set 



## ---------------------------------------------------------------------------------------------------------------
#by hand number of alleles
as.matrix(gls[1:5,])
#count the number of alleles per locus in each popu

 

nas1 <- gl.report.allelerich(gls)
nas1$`Allelic Richness per population`

nas2 <- gl.report.diversity(gls, table = "D")




## ---------------------------------------------------------------------------------------------------------------
gg <- gl.report.nall(gls, simlevels = 1:10, reps = 20, ncores = 1)  #change the number of cores if you have more available


## ---------------------------------------------------------------------------------------------------------------
gl.report.heterozygosity(gls, method = "ind")
gl.report.heterozygosity(gls, method = "pop")




## ---------------------------------------------------------------------------------------------------------------
#create a sample data sets (two populations, once with 10 individuals and one with 5 individuals)

gls1 <- possums.gl[c(1:10, 31:35),]

#cr
subfun <- function(x) { 
  xx <- gl.subsample.ind(x, n = 5, replace = TRUE, verbose = 0)
  out <- gl.report.heterozygosity(xx, method = "pop", verbose = 0)
  return(out$Ho)
  
  }

subfun(gls1)
subfun(gls1)
#ignore the warning
res <- sapply(1:50, function(x) subfun(gls1) ) 
rownames(res) <- popNames(gls1)
boxplot(t(res))
summary(t(res))




## ---------------------------------------------------------------------------------------------------------------
gg <- gl.report.heterozygosity(gls)

gg





## ---------------------------------------------------------------------------------------------------------------
#does not make much sense (sample size too low)
gl.report.hwe(gls,subset = "each",min_sample_size = 1 )

# only 10 loci...
gl.report.hwe(possums.gl[1:90,],subset = "each")


## ---------------------------------------------------------------------------------------------------------------
gl.report.heterozygosity(gls)




## ---------------------------------------------------------------------------------------------------------------
gl.map.interactive(possums.gl[1:120,])
gl.report.pa(possums.gl[1:120,] )

gl.fixed.diff(possums.gl[1:120,])  


## ---------------------------------------------------------------------------------------------------------------
#simulate some relatedness data

glsim <- gl.sim.Neconst(ninds = 50, nlocs = 1000)
Amat <- gl.grm(glsim)
### F is diag(A)-1
# centered around 1
hist(diag(Amat)-1)

#proportion of shared alleles  IBD between individuals
relA <- Amat
diag(relA)<-NA 
hist(relA)






## ---------------------------------------------------------------------------------------------------------------
#create an 10 offsprings from indivivdual 1 and 2
offsprings <- gl.sim.offspring(glsim[1,], glsim[2,], noffpermother = 10, popname = "off12")
glsimoff <- rbind(glsim, offsprings)

Amat <- gl.grm(glsimoff)



## ---------------------------------------------------------------------------------------------------------------
# matrix between offsprings and parents

Apo <- Amat[c(1,2,51:60),c(1,2,51:60) ]
round(Apo,2)

## only related to parents and each other
relpo <- Apo
diag(relpo)<-NA
hist(relpo, main = "Relatedness between offsprings and parents", xlab = "Relatedness", ylab = "Frequency")




## ---------------------------------------------------------------------------------------------------------------

binary_dir <- './binaries/'
#install emibd9 from Wang 2022 in the temporary directory
dir <- gl.download.binary("emibd9", out.dir = binary_dir)



## ---------------------------------------------------------------------------------------------------------------

#ignore the warnings...
kinmat <- gl.run.EMIBD9(glsimoff, Inbreed = 1, emibd9.path = binary_dir)



## ---------------------------------------------------------------------------------------------------------------

### F is diag(kinmat$rel)-1
hist((2*diag(kinmat$rel) - 1), main = "Kinship", xlab = "Kinship", ylab = "Frequency" )


#kinship values between parents and offsprings and between offsprings =0.25
kinpo <- kinmat$rel[c(1,2,51:60),c(1,2,51:60) ]
diag(kinpo)<- NA
hist(kinpo, main = "Kinship between offsprings and parents", xlab = "Kinship", ylab = "Frequency" )

#variance is different. 
#Parent off spring variance is zero.





## ---------------------------------------------------------------------------------------------------------------
#compare the two matrices
Arel <- Amat[lower.tri(Amat, diag = FALSE)]
kinrel <- kinmat$rel[lower.tri(kinmat$rel, diag = FALSE)]
plot(Arel, kinrel, 
     xlab = "Relatedness (gl.grm)", 
     ylab = "Kinship (EMIBD9)", 
     main = "Comparison of F values from gl.grm and EMIBD9",
     pch = 19, col = "blue")


## ---------------------------------------------------------------------------------------------------------------
#install the Neestimator package

dir <- dartRverse::gl.download.binary("neestimator", out.dir = binary_dir)


dir_Ne <- './binaries/Neestimator/'

## ---------------------------------------------------------------------------------------------------------------
#simulate a population of 50 individuals with 1000 loci

sim50 <- gl.sim.Neconst(ninds = 50, nlocs = 3000)
#need to reproduce a bit (to get the population to be in Hardy-Weinberg equilibrium)

sim50 <- gl.sim.offspring(sim50, sim50, noffpermother = 1) #ideal population
sim50 <- gl.sim.offspring(sim50, sim50, noffpermother = 1) #ideal population
sim50 <- gl.sim.offspring(sim50, sim50, noffpermother = 1) #ideal population
sim50 <- gl.sim.offspring(sim50, sim50, noffpermother = 1) #ideal population
sim50 <- gl.sim.offspring(sim50, sim50, noffpermother = 1) #ideal population


gg <-gl.LDNe(sim50, neest.path = dir_Ne, mating = "random", critical = c(0.1,0.05))
         




## ---------------------------------------------------------------------------------------------------------------
#install the Neestimator package

dir <- dartRverse::gl.download.binary("epos", out.dir = binary_dir)

## ---------------------------------------------------------------------------------------------------------------
gl.sfs(sim50)


## ---------------------------------------------------------------------------------------------------------------
#parameters for epos
# mutation rate
u =1e-8  #from simulation

#To get L we need estimate the Length of the genome in base pairs.
# watterson estimate for a sample of 50 is 4.499
L <- sum(gl.sfs(sim50)[-1])/( 4 * 50 *1e-8 * 4.499)
  

gepos <- gl.run.epos(sim50, epos.path = dir, L=L, u=u, boot=10, minbinsize = 2)


## ---------------------------------------------------------------------------------------------------------------
library(dartRverse)
north <- readRDS("./data/TympoNorth.rds")


#filters 
north2 <- gl.filter.callrate(north, method = "loc", threshold = 0.95)
north3 <- gl.filter.rdepth(north2, lower = 10, upper=40)
north4 <- gl.filter.monomorphs(north3)
north5 <- gl.filter.callrate(north4, threshold = 0.9, method="ind")

#checks
north5
gl.sfs(north5)
nInd(north5)
nLoc(north5)


### takes 10 minutes or so....
#curNe <- gl.LDNe(north5, neest.path = dir, critical = 0.05)  #-> ~40


#for dart data (times 75), 69+
L <- nLoc(north5)*75*200
mu = 16.17e-9


# takes a minute or so
gg <- gl.run.epos(north5, epos.path = dir,L = L, u = mu, method = "greedy", depth = 2, boot=5)

gg$plot

