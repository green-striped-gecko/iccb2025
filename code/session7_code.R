## Session 3 - Identifying Population Structure

# follow the links (ctrl click) to find the accompanying session material
# https://green-striped-gecko.github.io/iccb2025/session03.html


## --------------------------------------------
library(dartRverse)
library(ggplot2)
library(tidyr)


## --------------------------------------------------------------------------
rfbe <- readRDS("./data/rfbe.rds")
gl.report.basics(rfbe)


## --------------------------------------------------------------------------
source <- gl.keep.pop(rfbe, pop.list="PJTub1.2.3")
nInd(source)
source



## --------------------------------------------------------------------------
#remove all missing data
source2 <- gl.filter.callrate(source, method="loc", threshold=1)

#impute
#source2 <- gl.impute(source, method="random")




## --------------------------------------------------------------------------
#create a new genlight object based on allele frequencies from source2
transfer <- gl.subsample.ind(source2, n=10, replace = FALSE, verbose = 0)
res<- mean(gl.He(transfer))   #mean heterozygosity 
res




## --------------------------------------------------------------------------
for (gen in 2:21){
  #cloning snails
  transfer <- gl.sim.offspring(transfer, transfer, noffpermother = 1)
  res[gen] <- mean(gl.He(transfer))
  
}
  
plot(res, type="b", xlab="Generation", ylab="Expected heterozygosity", 
     main="Simulation of expected heterozygosity over time")




## --------------------------------------------------------------------------

simHe <- function(x, nInd=10, ngens=20){
  #remove all missing data

  
  #create a new genlight object based on allele frequencies from source2
  transfer <- gl.subsample.ind(x, n=nInd, replace = FALSE, verbose = 0)
  
  res <- mean(gl.He(transfer))   #mean heterozygosity 
  
  
  for (gen in 2:ngens)
  {
    #cloning snails
    transfer <- gl.sim.offspring(transfer, transfer, noffpermother = 1)
    res[gen] <- mean(gl.He(transfer))
    
  }
  
  return(res)
}
#test the function

out <- simHe(source2, nInd = 30, ngens = 20)
plot(out, type="b", xlab="Generation", ylab="Expected heterozygosity",
     main="Simulation of expected heterozygosity over time")



## --------------------------------------------------------------------------

nreps <- 10
out <- lapply(1:nreps, function(x) simHe(source2, nInd = 10, ngens = 20))



## --------------------------------------------------------------------------
outmat <- data.frame(do.call(cbind, out))

matplot(outmat, type="b", xlab="Generation", ylab="Expected heterozygosity",
        main="Simulation of expected heterozygosity over time", col=1:nreps, pch=1:nreps)





## --------------------------------------------------------------------------
library(ggplot2)
outmat2 <- pivot_longer(outmat, cols=everything(), names_to="replicate", values_to="He")
outmat2$generation <- rep(1:20, each=nreps)
ggplot(outmat2, aes(x=generation, y=He, group=generation)) +
  geom_boxplot() +
  labs(x="Generation", y="Expected heterozygosity",
       title="Simulation of expected heterozygosity over time") +
  theme_minimal() +
  scale_x_continuous(breaks=1:20)



## --------------------------------------------------------------------------

#use only 1000 loci (due ot speed)
glsim1 <- gl.sim.ind(source2[,1:1000], n = 100)

gl.smearplot(glsim1)




## --------------------------------------------------------------------------

glsim1 <-  gl.sim.ind(source2, n = 10, popname = "pop1")
glsim2 <-  gl.sim.ind(source2, n = 10, popname = "pop2")
fst <- NA
for (i in 1:20) {
  
 glsim1 <- gl.sim.ind(glsim1, n=10, popname = "pop1")
 glsim2 <- gl.sim.ind(glsim2, n=10, popname = "pop2")
 gg <- rbind(glsim1, glsim2)
fst[i] <- gl.fst.pop(gg, verbose = 0)[2,1]
}


plot(fst, type="b", xlab="Generation", ylab="Fst", 
     main="Simulation of Fst over time")



## --------------------------------------------------------------------------
simalllos <- function(x, ngens=30, nind=100, nloc=1000, 
                      surv=0.8, repro=3, K=100){
  #create a population of individuals
  dummy <- gl.sim.ind(x, n = nind, popname = "pop1")
  dummy <- dummy[,1:nloc]
    #allocate sex
    
    pop(dummy) <- sample(c("M","F"), nInd(dummy),replace = TRUE)
  
  
  for (gen in 1:ngens){

  
    #offsprings
    
    offsprings <- gl.sim.offspring(fathers = dummy[pop="M"], mothers = dummy[pop="F"],
                                   noffpermother = repro)
    pop(offsprings) <- sample(c("M","F"), nInd(offsprings),replace = TRUE)
    
    #death of adults
    survival <- rbinom(nInd(dummy), size = 1, prob = surv) #survival
    dummy <- dummy[survival == 1, ] #keep only survivors
    
    #combine population
    dummy <- rbind(dummy, offsprings)
    
    #check population size > K
    if (nInd(dummy) > K) {
      #randomly remove individuals to keep population size at K
      remove <- sample(1:nInd(dummy), nInd(dummy) - K, replace = FALSE)
      dummy <- dummy[-remove, ]
    }
    #check alleles lost

  }
  alleleslost <- sum(colMeans(as.matrix(dummy)) %% 2 ==0)
  
  return(alleleslost)
}
#test 
start <- gl.filter.monomorphs(source2)

out <- simalllos(start, ngens=30, nind=20, nloc=1000, surv=0.8, repro=3, K=100)
out



## --------------------------------------------------------------------------
paras <- expand.grid(ngens=1:10, ninds=20, nlocs=1000, surv=c(0.5, 0.8, 0.9), repro=3, K=50)

res <- sapply(1:nrow(paras), function(i) {
  simalllos(start, ngens=paras$ngens[i], nind=paras$ninds[i], nloc=paras$nlocs[i],
            surv=paras$surv[i], repro=paras$repro[i], K=paras$K[i])
})



paras$alleleslost <- res

ggplot(paras, aes(x=ngens, y=alleleslost, color=factor(ninds), group=factor(ninds))) +
  geom_line() +
  labs(x="Generations", y="Alleles lost", color="Population size") +
  scale_color_brewer(palette="Set1") +
  facet_wrap(~surv)


## --------------------------------------------------------------------------
table(pop(rfbe))

#lets keep the first 9 populations
pops <- popNames(rfbe)[1:9]
rfbe2 <- gl.keep.pop(rfbe, pop.list=pops, verbose = 0)

#and impute missing data

rfbe3 <-  gl.impute(rfbe2,method="neighbour" )



## --------------------------------------------------------------------------
panel <- gl.select.panel(rfbe3, method="pic", nl=50, verbose = 0)
res <- gl.check.panel(panel, rfbe3, parameter="Fst")
out <- summary(lm(res[,1] ~ res[,2]))$r.squared
out


## --------------------------------------------------------------------------

xx <- rfbe3
xx <- xx[order(pop(xx)),]


for ( i in 2:10) {
  pops <- seppop(xx)
  #panmictiv pops, but completely seperated
  pops <- lapply(pops, function(yy)   {
    dummy <- gl.sim.offspring(yy,yy, noffpermother = 3)
    dummy <- gl.sample(dummy, nInd(yy), replace = FALSE, verbose = 0)
  })
  pops <- do.call(rbind, pops)  
  pop(pops)<- pop(xx)
  
  res <- gl.check.panel(x =panel, pops, parameter="Fst")
out[i] <- summary(lm(res[,1] ~ res[,2]))$r.squared
xx<- pops
  
cat(paste("Generation", i, "R-squared:", out[i], "\n"))
flush.console()  
  
}


plot(out, type="b", xlab="Generation", ylab="R-squared", 
     main="Performance of SNP panel over time")



