### R code from vignette source 'glpk-gmpl-intro.Rnw'

###################################################
### code chunk number 1: glpk-gmpl-intro.Rnw:23-24
###################################################
library(glpkAPI)


###################################################
### code chunk number 2: glpk-gmpl-intro.Rnw:49-56
###################################################
mip <- initProbGLPK()
setProbNameGLPK(mip, "transport")
trans <- mplAllocWkspGLPK()
result <- mplReadModelGLPK(trans,
	system.file("extdata", "transport.mod", package = "glpkAPI"), skip=0)
result <- mplGenerateGLPK(trans)
result <- mplBuildProbGLPK(trans, mip)


###################################################
### code chunk number 3: glpk-gmpl-intro.Rnw:69-76
###################################################
numrows <- getNumRowsGLPK(mip)
numrows


for (i in 1:numrows){
	print(getRowNameGLPK(mip, i))
}


###################################################
### code chunk number 4: glpk-gmpl-intro.Rnw:81-88
###################################################
numcols <- getNumColsGLPK(mip)
numcols

for (j in 1:numcols){
	print(getColNameGLPK(mip, j))
}
print(getNumNnzGLPK(mip))


###################################################
### code chunk number 5: glpk-gmpl-intro.Rnw:93-95
###################################################
return <- solveSimplexGLPK(mip)
return <- mplPostsolveGLPK(trans, mip, GLP_MIP);


###################################################
### code chunk number 6: glpk-gmpl-intro.Rnw:100-104
###################################################
for (i in 1:numrows){
print(getRowNameGLPK(mip, i))
print(getRowPrimGLPK(mip, i))
}


###################################################
### code chunk number 7: glpk-gmpl-intro.Rnw:109-113
###################################################
for (j in 1:numcols){
print(getColNameGLPK(mip, j))
print(getColPrimGLPK(mip, j))
}


###################################################
### code chunk number 8: glpk-gmpl-intro.Rnw:118-120
###################################################
mplFreeWkspGLPK(trans)
delProbGLPK(mip)


###################################################
### code chunk number 9: glpk-gmpl-intro.Rnw:128-136
###################################################
print ("USING API")
canneries <- c("Seattle", "San-Diego")
capacity <- c(350, 600)
markets <- c("New-York", "Chicago", "Topeka")
demand <- c(325, 300, 275)
distance <- c(2.5, 2.5, 1.7, 1.8, 1.8, 1.4)
dim(distance) <- c(2, 3)
freight <- 90


###################################################
### code chunk number 10: glpk-gmpl-intro.Rnw:141-145
###################################################
lpi <- initProbGLPK()
setProbNameGLPK(lpi, "cannery API")
setObjNameGLPK(lpi, "Total Cost")
setObjDirGLPK(lpi, GLP_MIN)


###################################################
### code chunk number 11: glpk-gmpl-intro.Rnw:150-164
###################################################
numlinks <- length(distance)
nummarkets <- length(markets)
numcanneries <- length(canneries)
addColsGLPK(lpi, numlinks)
for (i in 1:numcanneries){
	cannerystartrow <- (i-1) * nummarkets
	for (j in 1:nummarkets){
		colname <-toString(c(canneries[i], markets[j]))
		transcost <- distance[i, j]*freight/1000
		setColNameGLPK(lpi, cannerystartrow+j, colname)
		setColBndGLPK(lpi, cannerystartrow+j, GLP_LO, 0.0, 0.0)
		setObjCoefsGLPK(lpi, cannerystartrow+j, transcost)
	}
}


###################################################
### code chunk number 12: glpk-gmpl-intro.Rnw:169-181
###################################################
numcanneries <- length(canneries)
nummarkets <- length(markets)
addRowsGLPK(lpi, numcanneries+nummarkets+1)
setRowsNamesGLPK(lpi, 1, getObjNameGLPK(lpi))
for (i in 1:numcanneries){
	setRowsNamesGLPK(lpi, i+1, toString(c("Supply", canneries[i])))
	setRowBndGLPK(lpi, i+1, GLP_UP, 0, capacity[i])
}
for (j in 1:nummarkets){
	setRowsNamesGLPK(lpi, numcanneries+j+1, toString(c("Demand", markets[j])))
	setRowBndGLPK(lpi, numcanneries+j+1, GLP_LO, demand[j], 0)
}


###################################################
### code chunk number 13: glpk-gmpl-intro.Rnw:185-216
###################################################
# create variables to hold the constraint information
ia <- numeric()
ja <- numeric()
ar <- numeric()

# add in objective coefficients

for (i in 1:numcols){
	ia[i] <- 1
	ja[i] <- i
	ar[i] <- getObjCoefGLPK(lpi, i)
}

for (i in 1:numcanneries){
	#supply constraints
	cannerysupplyrow = numcols + (i-1)*nummarkets
	for (j in 1:nummarkets){
		ia[cannerysupplyrow+j] <- (i+1)
		ja[cannerysupplyrow+j] <- (i-1)+numcanneries *(j-1)+1
		ar[cannerysupplyrow+j] <- 1
	}
	#demand constraints
	marketdemandrow = numcols+numcanneries * nummarkets
	for (j in 1:nummarkets){
		colnum <- (i-1)*nummarkets+j
		ia[marketdemandrow + colnum] <- numcanneries+j+1
		ja[marketdemandrow + colnum] <- colnum
		ar[marketdemandrow + colnum] <- 1
	}
}
loadMatrixGLPK(lpi, length(ia), ia, ja, ar)


###################################################
### code chunk number 14: glpk-gmpl-intro.Rnw:221-232
###################################################
numrows <- getNumRowsGLPK(lpi)
numrows
numcols <- getNumColsGLPK(lpi)
numcols
for (i in 1:numrows){
	print(getRowNameGLPK(lpi, i))
}
for (j in 1:numcols){
	print(getColNameGLPK(lpi, j))
}
print(getNumNnzGLPK(lpi))


###################################################
### code chunk number 15: glpk-gmpl-intro.Rnw:237-247
###################################################
solveSimplexGLPK(lpi)

for (i in 1:numrows){
print(getRowNameGLPK(lpi, i))
print(getRowPrimGLPK(lpi, i))
}
for (j in 1:numcols){
print(getColNameGLPK(lpi, j))
print(getColPrimGLPK(lpi, j))  
}


###################################################
### code chunk number 16: glpk-gmpl-intro.Rnw:252-253
###################################################
printSolGLPK(lpi, "transout.api")


###################################################
### code chunk number 17: glpk-gmpl-intro.Rnw:262-269
###################################################
cindex <- createIndexGLPK(lpi)
new_york_row = findRowGLPK(lpi, "Demand, New-York") 
topeka_row = findRowGLPK(lpi, "Demand, Topeka")
new_york_row
topeka_row
setRowBndGLPK(lpi, new_york_row, GLP_LO, 300, 0)
setRowBndGLPK(lpi, topeka_row, GLP_LO, 300, 0)


###################################################
### code chunk number 18: glpk-gmpl-intro.Rnw:274-287
###################################################
solveSimplexGLPK(lpi)

for (i in 1:numrows){
print(getRowNameGLPK(lpi, i))
print(getRowPrimGLPK(lpi, i))
print(getRowDualGLPK(lpi, i))
}
for (j in 1:numcols){
print(getColNameGLPK(lpi, j))
print(getColPrimGLPK(lpi, j))
print(getColDualGLPK(lpi,j))
print(getObjCoefGLPK(lpi, j))
}


###################################################
### code chunk number 19: glpk-gmpl-intro.Rnw:292-293
###################################################
delProbGLPK(lpi)


