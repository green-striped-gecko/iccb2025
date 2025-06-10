

## library -----------------------------------------------------------------------------------------------
library(dartRverse)


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


# reports -----------------------------------------------------------------

# Now lets use a report function to report the call rate for genlight object gl
gl.report.callrate(gl)

gl.report.callrate(gl,method="ind")

gl.report.reproducibility(gl)

# First look at a report to decide a threshold
gl.report.callrate(gl,method="ind")

# Then filter using that threshold
gl <- gl.filter.callrate(gl, method="ind", threshold=0.80)

# Use a smear plot for a visual assessment of the effectiveness of filtering
gl <- testset.gl
gl.smearplot(gl)



## more filtering -----------------------------------------------------------------------------------------------------------
# Filter on callrate by locus, then on individual
gl <- gl.filter.callrate(gl,verbose=0)
gl <- gl.filter.callrate(gl, method= "ind", threshold=0.80, verbose=0)

# Examine the smearplot again
gl.smearplot(gl)


## ----output=FALSE---------------------------------------------------------------------------------------------------------------
# how to load a genlight object
gl <- gl.read.dart(filename="./data/sample_data_2Row.csv", 
                   ind.metafile="./data/sample_metadata.csv")

