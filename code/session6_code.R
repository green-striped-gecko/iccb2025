## ----warning=FALSE, message=FALSE-------------------------------------------------------------------------------
library(dartRverse)
library(webexercises)
library(knitr)


## ---------------------------------------------------------------------------------------------------------------

rfbe <- readRDS("./data/rfbe.rds")
gl.report.basics(rfbe)



## ---------------------------------------------------------------------------------------------------------------
tt <- table(pop(rfbe))
pop20 <- names(tt)[tt>10]

rfbe20 <- gl.keep.pop(rfbe, pop.list=pop20)
kable(table(pop(rfbe20)))



## ---------------------------------------------------------------------------------------------------------------
rfbe20_1 <- gl.filter.allna(rfbe20, by.pop = T)
rfbe20_2 <- gl.filter.callrate(rfbe20_1, threshold=0.99)
rfbe20_3 <- gl.filter.maf(rfbe20_2, threshold = 5, by.pop = FALSE)
nLoc(rfbe20_3)


## ---------------------------------------------------------------------------------------------------------------

index <- nchar(as.character(rfbe20_3@other$loc.metrics$TrimmedSequence))>29
rfbe20_4 <- rfbe20_3[, index]



## ---------------------------------------------------------------------------------------------------------------
#takes a while to run, hence load the pre run data set
#rfbe20_5 <- gl.blast(rfbe20_4,ref_genome = "d:/bernd/r/Elise_pansnp/final.genome.scf.fasta", task = "blastn", number_of_threads = 10 )
#saveRDS(rfbe20_5, "d:/bernd/r/Elise_pansnp/rfbe20_5_blast.rds")

rfbe20_5 <- readRDS("./data/rfbe20_5_blast.rds")



## ---------------------------------------------------------------------------------------------------------------
#bitscore >=100
index <- rfbe20_5@other$loc.metrics$bitscore >= 100
index <- ifelse(is.na(index), FALSE, index)
rfbe20_6 <- rfbe20_5[,index]



## ---------------------------------------------------------------------------------------------------------------
#takes a while to run
panel <- gl.select.panel(rfbe20_6, method="dapc", nl = 50)
nLoc(panel)


## ---------------------------------------------------------------------------------------------------------------
outdapc <- gl.check.panel(panel, rfbe20_6, parameter = "Fst")


## ---------------------------------------------------------------------------------------------------------------
panel_random <- gl.select.panel(rfbe20_6, method="random", nl = 50)
outrandom <- gl.check.panel(panel_random, rfbe20_6, parameter = "Fst")

