# Session 1 Importing and Handling Genomic Data

# follow the links (ctrl click) to find the accompanying session material
# https://green-striped-gecko.github.io/iccb2025/session01.html

# library -----------------------------------------------------------------------------------------------
library(dartRverse)


# dartR Fundamentals -----------------------
# https://green-striped-gecko.github.io/iccb2025/session01.html#dartr-fundamentals

## data ---------------------------------------------------------------------------------------------------------------
# dartR comes with a built in test dataset called testset.gl. 
# We first examine the contents of testset.gl
testset.gl


## interrogate ---------------------------------------------------------------------------------------------------------------
# Next let us copy the contents of testset.gl to another genlight object called gl.
gl <- testset.gl
# Use adegenet accessors to interrogate the genlight object further
nInd(gl)
nLoc(gl)
nPop(gl)
popNames(gl)
indNames(gl)
locNames(gl)


## individuals per pop --------------------------------------------------------------------------------------------------------------------
# So if you want to tablulate the number of individuals in each population, use
table(pop(gl))


## matrix of snp data---------------------------------------------------------------------------------------------------------------
as.matrix(gl)[1:7,1:5]


## Again1 ---------------------------------------------------------------------------------------------------------------
# Copy the original testset.gl
gl <- testset.gl
# Set the global verbosity to 3 – this will result in some detailed comments as you run each # script
gl.set.verbosity(3)


## reports -----------------------------------------------------------------
# https://green-striped-gecko.github.io/iccb2025/session01.html#reporting
# Now lets use a report function to report the call rate for genlight object gl
gl.report.callrate(gl)

gl.report.callrate(gl,method="ind")

gl.report.reproducibility(gl)



## filtering ---------------------------------------------------------------
# http://localhost:6901/session01.html#filtering


# First look at a report to decide a threshold
gl.report.callrate(gl,method="ind")

# Then filter using that threshold
gl <- gl.filter.callrate(gl, method="ind", threshold=0.80)

# Use a smear plot for a visual assessment of the effectiveness of filtering
gl <- testset.gl
gl.smearplot(gl)

# Filter on callrate by locus, then on individual
gl <- gl.filter.callrate(gl,verbose=0)
gl <- gl.filter.callrate(gl, method= "ind", threshold=0.80, verbose=0)

# Examine the smearplot again
gl.smearplot(gl)


# sex linked markers ------------------------------------------------------
## load data--------------------------------------------------------------
LBP                   # Explore the dataset
LBP@n.loc             # Number of SNPs
length(LBP@ind.names) # Number of individuals


head(LBP@other$ind.metrics)

## function help page ----------------------------------------------
help(gl.report.sexlinked)

## run report ----------------------------------------------------------------------------------------------------------------------------------
out <- gl.report.sexlinked(LBP, system = "xy")

out

## filte sex-linked ---------------------------------------------------------------------------------------------------------
new.gl <- gl.drop.sexlinked(LBP, system = "xy")
new.gl
new.gl@n.loc



# load data ---------------------------------------------------------------

# how to load a genlight object
gl <- gl.read.dart(filename="./data/sample_data_2row.csv", 
                   ind.metafile="./data/sample_metadata.csv")



