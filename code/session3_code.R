## Session 3 - Identifying Population Structure

# follow the links (ctrl click) to find the accompanying session material
# https://green-striped-gecko.github.io/iccb2025/session03.html


## ----------------------------------------------------------------------------------
library(dartRverse)


## --------------------------------------------------------------------------------------------------------------

table(pop(possums.gl)) #check the individuals and the populations

gl.map.interactive(possums.gl)



## --------------------------------------------------------------------------------------------------------------
# Undertake a PCA on the raw data
pc <- gl.pcoa(possums.gl, verbose = 3)



## -----------------------------------------------------------------------------
# Plot the first two dimensions of the PCA
pc_a1a2 <- gl.pcoa.plot(glPca = pc, x = possums.gl)

# Plot the first and third dimensions of the PCA
pc_a1a3 <- gl.pcoa.plot(glPca = pc, x = possums.gl, xaxis = 1, yaxis = 3)

# Plot the first three dimensions of the PCA
pc_a1a3 <- gl.pcoa.plot(glPca = pc, x = possums.gl, xaxis = 1, yaxis = 2, zaxis = 3)



## --------------------------------------------------------------------------------------------------------------
# Select only the data from one cluster in the primary PCA
temp <- gl.drop.pop(x = possums.gl, pop.list = c('D', 'A', 'E', 'F', 'H', 'G'))
# Plot the first two dimensions of the secondary PCA
pc <- gl.pcoa(temp, verbose = 3)
pc_plot <- gl.pcoa.plot(glPca =  pc, x = temp)



## --------------------------------------------------------------------------------------------------------------

structure_file <- ifelse('structure.exe' %in% list.files('./binaries/'), 
                         './binaries/structure.exe', './binaries/console/structure')
srnoad <- gl.run.structure(possums.gl, k.range = 2:7, num.k.rep = 2, 
                           exec = structure_file,plot.out = FALSE,
                           burnin=500, numreps=1000, 
                           noadmix=TRUE)




## --------------------------------------------------------------------------------------------------------------
ev <- gl.evanno(srnoad)


## --------------------------------------------------------------------------------------------------------------
qmatnoad <- gl.plot.structure(srnoad, K=3:5)
head(qmatnoad[[1]])


## ----structure map---------------------------------------------------------------------------------------------
gm <- gl.map.structure(qmat = qmatnoad, x = possums.gl,K=5, scalex=1, scaley=0.5 )



## --------------------------------------------------------------------------------------------------------------
platypus.gl

gl.map.interactive(platypus.gl)


## ------------------------------------------------------------------------------------------
my_fast <- gl.run.faststructure(platypus.gl,
                                k.range = 2:4,
                                num.k.rep = 1,
                                exec = "./binaries/fastStructure",
                                exec.plink = "./binaries/", output = tempdir())

gl.plot.faststructure(sr = my_fast,k.range = 3, border_ind = 0)



## --------------------------------------------------------------------------------
my_snmf <- gl.run.snmf(possums.gl, minK = 2, maxK = 7, rep = 2, regularization = 10)

gl.plot.snmf(snmf_result = my_snmf,plot.K = 3, border_ind = 0)



## --------------------------------------------------------------------------------------------------------------
my_popcluster <- gl.run.popcluster(x = possums.gl, minK = 2, maxK = 7, rep = 2,
                                   popcluster.path = './binaries/')

my_plot_popcluster <- gl.plot.popcluster(my_popcluster,plot.K = 3)

gl.map.popcluster(x = possums.gl, qmat =  my_plot_popcluster)


