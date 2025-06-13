## Session 8 - Sex Linked Markers

# follow the links (ctrl click) to find the accompanying session material
# https://green-striped-gecko.github.io/iccb2025/session08.html


## -----------------------------------------------------------------------------------------------------------------------
library(dartRverse)
library(dartR.sexlinked)


## ----------------------------------------------------------------------------------------------------------------------------------------------------
data("EYR")
EYR                   # Explore the dataset
EYR@n.loc             # Number of SNPs
length(EYR@ind.names) # Number of individuals


## ----------------------------------------------------------------------------------------------------------------------------------------------------
EYR@other$ind.metrics


## ----------------------------------------------------------------------------------------------------------------------------------------------------
out <- gl.report.sexlinked(EYR, system = "zw")


## ----------------------------------------------------------------------------------------------------------------------------------------------------
EYR@other$ind.metrics[!(EYR@other$ind.metrics$sex %in% c("M", "F")), ]


## ----------------------------------------------------------------------------------------------------------------------------------------------------
EYR_sexLinked <- gl.keep.sexlinked(EYR, system = "zw") # save sex-linked loci
inferred.sexes <- gl.infer.sex(gl_sexlinked = EYR_sexLinked, 
                               system = "zw", seed = 124) # use sex-linked loci


## ----------------------------------------------------------------------------------------------------------------------------------------------------
head(inferred.sexes, 10) 


## ----------------------------------------------------------------------------------------------------------------------
# Filter for read depth
gl.report.rdepth(EYR)
EYR.sloppy <- gl.filter.rdepth(EYR, lower = 3, upper = 11, verbose = 0)

# Filter for loci call rate
gl.report.callrate(EYR.sloppy, method = "loc")
EYR.sloppy <- gl.filter.callrate(EYR.sloppy, method = "loc",  threshold = 0.75, verbose = 0, recalc = TRUE)

# Filter for individual call rate
gl.report.callrate(EYR.sloppy, method = "ind")
EYR.sloppy <- gl.filter.callrate(EYR.sloppy, method = "ind", threshold = 0.65, verbose = 0, recalc = TRUE)

# Filter for MAC (= 3)
gl.report.maf(EYR.sloppy)
EYR.sloppy <- gl.filter.maf(EYR.sloppy, threshold = 3, verbose = 0, recalc = TRUE)


## ------------------------------------------------------------------------------------------------------------------------
# Filter for sex-linked loci
EYR.correct <- gl.drop.sexlinked(EYR, system = "zw")  

# Filter for read depth
EYR.correct <- gl.filter.rdepth(EYR.correct, lower = 3, upper = 11, verbose = 0)

# Filter for loci call rate
gl.report.callrate(EYR.correct, method = "loc")
EYR.correct <- dartR.base::gl.filter.callrate(EYR.correct, method = "loc",  threshold = 0.75, verbose = 0, recalc = TRUE)

# Filter for individual call rate
gl.report.callrate(EYR.correct, method = "ind")
EYR.correct <- gl.filter.callrate(EYR.correct, method = "ind", threshold = 0.65, verbose = 0, recalc = TRUE)

# Filter for MAC (= 3)
gl.report.maf(EYR.correct)
EYR.correct <- dartR.base::gl.filter.maf(EYR.correct, threshold = 3, verbose = 0, recalc = TRUE)


## --------------------------------------------------------------------------------------------------------
# Sloppy
PCA.sloppy <- gl.pcoa(EYR.sloppy, verbose = 0)
pcplot_sloppy <- gl.pcoa.plot(PCA.sloppy, EYR.sloppy, xaxis = 1, yaxis = 2)

# Correct
PCA.correct <- gl.pcoa(EYR.correct, verbose = 0)
pcplot_correct <- gl.pcoa.plot(PCA.correct, EYR.correct, xaxis = 1, yaxis = 2)


## --------------------------------------------------------------------------------------------------------------------------------------
basic.sloppy  <- utils.basic.stats(EYR.sloppy)
basic.correct <- utils.basic.stats(EYR.correct)
basic.sloppy$overall
basic.correct$overall


## ---------------------------------------------------------------------------------------------------------------------------------------
gl.fst.pop(EYR.sloppy, verbose = 0)
gl.fst.pop(EYR.correct, verbose = 0)

