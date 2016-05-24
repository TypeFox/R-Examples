## ------------------------------------------------------------------------
#Load package for calculation of multivariate richness & sensitivity analysis tools
library(multirich) 

#Set up labels
sp.lbl = sprintf("sp%s",seq(1,15,1))
com.lbl = c("pool","com1","com2","com3")
tr.lbl = c("tr1","tr2")

#Set up traits and species x trait matrix
tr1 = c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5)
tr2 = c(1,2,3,4,5,1,2,3,4,1,2,3,1,2,1)
in.mat = matrix(c(tr1,tr2),ncol = 2, dimnames = list(sp.lbl,tr.lbl))

#Set up communities
pool = rep(1,15)
com1 = c(1,0,0,0,1,0,0,0,0,0,0,1,0,0,1)
com2 = c(1,1,0,0,0,1,1,0,0,0,0,0,0,0,0)
com3 = c(1,0,0,0,0,0,1,0,0,0,0,0,0,0,0)
in.com = matrix(c(pool,com1,com2,com3),nrow = 4,byrow = T,dimnames = list(com.lbl,sp.lbl))

tr1.breaks = tr2.breaks = get.breaks(1,5)
breaks = list(tr1.breaks, tr2.breaks)
out.pdf = "none" #Specifying a file here will save the result to file.

# Traitspace here is less than whole - because of a "known" trade-off between the two traits
in.traitspaces = c(3,6,10,15) 

# The results object can be used to do further manipulations of the output data
results = sensitivity.analysis(in.mat, in.com, breaks, out.pdf, in.traitspaces = in.traitspaces)


