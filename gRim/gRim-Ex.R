pkgname <- "gRim"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('gRim')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CGstats")
### * CGstats

flush(stderr()); flush(stdout())

### Name: CGstats
### Title: Mean, covariance and counts for grouped data
### Aliases: CGstats CGstats.data.frame CGstats_internal print.CGstats
### Keywords: utilities

### ** Examples


data(milkcomp)

CGstats(milkcomp)
CGstats(milkcomp, c(1,2))
CGstats(milkcomp, c("lactime","treat"))
CGstats(milkcomp, c(3,4))
CGstats(milkcomp, c("fat","protein"))

CGstats(milkcomp, c(2,3,4), simplify=FALSE)
CGstats(milkcomp, c(2,3,4), homogeneous=FALSE)
CGstats(milkcomp, c(2,3,4), simplify=FALSE, homogeneous=FALSE)





cleanEx()
nameEx("ciTest_df")
### * ciTest_df

flush(stderr()); flush(stdout())

### Name: ciTest_df
### Title: Test for conditional independence in a dataframe
### Aliases: ciTest_df
### Keywords: htest

### ** Examples

data(milkcomp1)
ciTest(milkcomp1, set=~tre+fat+pro)
ciTest_df(milkcomp1, set=~tre+fat+pro)



cleanEx()
nameEx("ciTest_generic2")
### * ciTest_generic2

flush(stderr()); flush(stdout())

### Name: ciTest
### Title: Generic function for conditional independence test
### Aliases: ciTest ciTest.data.frame ciTest.table ciTest.list print.citest
###   summary.citest
### Keywords: htest

### ** Examples


## contingency table:
data(reinis)
## dataframe with only numeric variables:
data(carcass)
## dataframe with numeric variables and factors:
data(milkcomp1)

ciTest(cov.wt(carcass, method='ML'), set=~Fat11+Meat11+Fat12)
ciTest(reinis, set=~smo+phy+sys)
ciTest(milkcomp1, set=~tre+fat+pro)




cleanEx()
nameEx("ciTest_mvn")
### * ciTest_mvn

flush(stderr()); flush(stdout())

### Name: ciTest_mvn
### Title: Test for conditional independence in the multivariate normal
###   distribution
### Aliases: ciTest_mvn
### Keywords: htest

### ** Examples

data(carcass)
ciTest(cov.wt(carcass, method='ML'), set=~Fat11+Meat11+Fat12)
ciTest_mvn(cov.wt(carcass, method='ML'), set=~Fat11+Meat11+Fat12)



cleanEx()
nameEx("ciTest_ordinal")
### * ciTest_ordinal

flush(stderr()); flush(stdout())

### Name: ciTest_ordinal
### Title: A function to compute Monte Carlo and asymptotic tests of
###   conditional independence for ordinal and/or nominal variables.
### Aliases: ciTest_ordinal
### Keywords: htest

### ** Examples

library(gRim)
data(dumping, package="gRbase")

ciTest_ordinal(dumping, c(2,1,3), stat="jt", N=1000)
ciTest_ordinal(dumping, c("Operation","Symptom","Centre"), stat="jt", N=1000)
ciTest_ordinal(dumping, ~ Operation + Symptom + Centre, stat="jt", N=1000)

data(reinis)
ciTest_ordinal(reinis, c(1,3,4:6),N=1000)

# If data is a dataframe
dd     <- as.data.frame(dumping)
ncells <- prod(dim(dumping))
ff     <- dd$Freq
idx    <- unlist(mapply(function(i,n) rep(i,n),1:ncells,ff))
dumpDF <- dd[idx, 1:3]
rownames(dumpDF) <- 1:NROW(dumpDF)

ciTest_ordinal(dumpDF, c(2,1,3), stat="jt", N=1000)
ciTest_ordinal(dumpDF, c("Operation","Symptom","Centre"), stat="jt", N=1000)
ciTest_ordinal(dumpDF, ~ Operation + Symptom + Centre, stat="jt", N=1000)




cleanEx()
nameEx("ciTest_table")
### * ciTest_table

flush(stderr()); flush(stdout())

### Name: ciTest_table
### Title: Test for conditional independence in a contingency table
### Aliases: ciTest_table
### Keywords: htest

### ** Examples

data(reinis)
ciTest(reinis, set=~smo+phy+sys)
ciTest_table(reinis, set=~smo+phy+sys)



cleanEx()
nameEx("cmod")
### * cmod

flush(stderr()); flush(stdout())

### Name: cmod
### Title: Graphical Gaussian model
### Aliases: cmod
### Keywords: models

### ** Examples

## Graphical Gaussian model
data(carcass)
cm1<-cmod(~.^., carcass)

## Stepwise selection based on BIC
cm2<-backward(cm1,k=log(nrow(carcass)))

## Stepwise selection with fixed edges
cm3<-backward(cm1,k=log(nrow(carcass)),
 fixinMAT=matrix(c("LeanMeat","Meat11","Meat12","Meat13","LeanMeat","Fat11","Fat12","Fat13"),
 ncol=2))



cleanEx()
nameEx("dModel-class")
### * dModel-class

flush(stderr()); flush(stdout())

### Name: dModel-class
### Title: Class '"dModel"'
### Aliases: dModel-class cModel-class mModel-class
### Keywords: classes

### ** Examples

showClass("dModel")



cleanEx()
nameEx("dmod")
### * dmod

flush(stderr()); flush(stdout())

### Name: dmod
### Title: Log-linear model
### Aliases: dmod print.dModel
### Keywords: models

### ** Examples


## Graphical log-linear model
data(reinis)
dm1<-dmod(~.^., reinis)
dm2<-backward(dm1, k=2)
dm3<-backward(dm1, k=2, fixin=list(c("family","phys","systol")))
## At most 3-factor interactions
dm1<-dmod(~.^., data=reinis,interactions=3)




cleanEx()
nameEx("effloglin")
### * effloglin

flush(stderr()); flush(stdout())

### Name: effloglin
### Title: Fitting Log-Linear Models by Message Passing
### Aliases: effloglin
### Keywords: models

### ** Examples

data(reinis)
glist <-list(c("smoke", "mental"), c("mental", "phys"), c("phys", "systol"
), c("systol", "smoke"))

stab <- lapply(glist, function(gg) tableMargin(reinis, gg))
fv3 <- effloglin(stab, glist, print=FALSE)





cleanEx()
nameEx("getEdges")
### * getEdges

flush(stderr()); flush(stdout())

### Name: getEdges
### Title: Find edges in a graph and edges not in an undirected graph.
### Aliases: getEdges getEdges.list getEdges.graphNEL getEdges.matrix
###   getInEdges getOutEdges getInEdgesMAT getOutEdgesMAT
### Keywords: utilities

### ** Examples

gg     <- ug(~a:b:d+a:c:d+c:e)
glist  <- getCliques(gg)
adjmat <- as.adjMAT(gg)

#### On a glist
getEdges(glist)
getEdges(glist,type="decomposable")
# Deleting (a,d) would create a 4-cycle

getEdges(glist, ingraph=FALSE)
getEdges(glist,type="decomposable", ingraph=FALSE)
# Adding (e,b) would create a 4-cycle

#### On a graphNEL
getEdges(gg)
getEdges(gg,type="decomposable")
# Deleting (a,d) would create a 4-cycle

getEdges(gg, ingraph=FALSE)
getEdges(gg,type="decomposable", ingraph=FALSE)
# Adding (e,b) would create a 4-cycle

#### On an adjacency matrix
getEdges(adjmat)
getEdges(adjmat,type="decomposable")
# Deleting (a,d) would create a 4-cycle

getEdges(adjmat, ingraph=FALSE)
getEdges(adjmat,type="decomposable", ingraph=FALSE)
# Adding (e,b) would create a 4-cycle


## Marked graphs; vertices a,b are discrete; c,d are continuous
UG <- ug(~a:b:c+b:c:d)
disc <- c("a","b")
getEdges(UG)
getEdges(UG, discrete=disc)
## Above: same results; there are 5 edges in the graph

getEdges(UG, type="decomposable")
## Above: 4 edges can be removed and will give a decomposable graph
##(only removing the edge (b,c) would give a non-decomposable model)

getEdges(UG, type="decomposable", discrete=c("a","b"))
## Above: 3 edges can be removed and will give a strongly decomposable
## graph. Removing (b,c) would create a 4--cycle and removing (a,b)
## would create a forbidden path; a path with only continuous vertices
## between two discrete vertices.




cleanEx()
nameEx("ggmfit")
### * ggmfit

flush(stderr()); flush(stdout())

### Name: ggmfit
### Title: Iterative proportional fitting of graphical Gaussian model
### Aliases: ggmfit ggmfitr
### Keywords: multivariate models

### ** Examples

## Fitting "butterfly model" to mathmark data
## Notice that the output from the two fitting functions is not
## entirely identical.
data(math)
ddd <- cov.wt(math, method="ML")
glist <- list(c("al","st","an"), c("me","ve","al"))
ggmfit (ddd$cov, ddd$n.obs, glist)
ggmfitr(ddd$cov, ddd$n.obs, glist)





cleanEx()
nameEx("iModel-stepwise")
### * iModel-stepwise

flush(stderr()); flush(stdout())

### Name: stepwise.iModel; backward; forward
### Title: Stepwise model selection in (graphical) interaction models
### Aliases: stepwise.iModel backward forward
### Keywords: models

### ** Examples

data(reinis)
## The saturated model
m1 <- dmod(~.^., data=reinis)
m2 <- stepwise(m1)
m2




cleanEx()
nameEx("loglinDim")
### * loglinDim

flush(stderr()); flush(stdout())

### Name: loglinDim
### Title: Return the dimension of a log-linear model
### Aliases: loglinGenDim loglinDecDim
### Keywords: models

### ** Examples


## glist contains variable names and tableinfo is a named vector:
loglinGenDim(list(c("a","b"),c("b","c")), c(a=4,b=7,c=6))

## glist contains variable names and tableinfo is not named:
loglinGenDim(list(c(1,2),c(2,3)), c(4,7,6))

## For decomposable models:
loglinDecDim(list(c("a","b"),c("b","c")), c(a=4,b=7,c=6),adjust=FALSE)



cleanEx()
nameEx("mmod")
### * mmod

flush(stderr()); flush(stdout())

### Name: mmod
### Title: Mixed interaction model.
### Aliases: mmod coef.mModel coefficients.mModel print.mModel
###   summary.mModel mmod_dimension
### Keywords: models

### ** Examples

### FIXME: To be written



cleanEx()
nameEx("modify_glist")
### * modify_glist

flush(stderr()); flush(stdout())

### Name: modify_glist
### Title: Modify generating class for a graphical/hierarchical model
### Aliases: modify_glist
### Keywords: utilities

### ** Examples

glist <- list(c(1,2,3),c(2,3,4))

## Add edges
modify_glist(glist, items=list(add.edge=c(1,4)))
modify_glist(glist, items=list(add.edge=~1:4))

## Add terms
modify_glist(glist, items=list(add.term=c(1,4)))
modify_glist(glist, items=list(add.term=~1:4))

## Notice: Only the first term is added as the second is already 
## in the model.
modify_glist(glist, items=list(add.term=list(c(1,4),c(1,3))))
modify_glist(glist, items=list(add.term=~1:4+1:3))

## Notice: Operations are carried out in the order given in the
## items list and hence we get different results: 
modify_glist(glist, items=list(drop.edge=c(1,4), add.edge=c(1,4)))
modify_glist(glist, items=list(add.edge=c(1,4), drop.edge=c(1,4)))



cleanEx()
nameEx("testEdges")
### * testEdges

flush(stderr()); flush(stdout())

### Name: testInEdges; testOutEdges
### Title: Test edges in graphical models with p-value/AIC value
### Aliases: testInEdges testOutEdges testEdges testEdges.iModel
### Keywords: models htest

### ** Examples

data(math)
cm1 <- cmod(~me:ve+ve:al+al:an, data=math)
testInEdges(cm1, getEdges(cm1$glist))
testOutEdges(cm1, getEdges(cm1$glist, ingraph=FALSE))



cleanEx()
nameEx("testadd")
### * testadd

flush(stderr()); flush(stdout())

### Name: testadd
### Title: Test addition of edge to graphical model
### Aliases: testadd testadd.iModel print.testadd testadd.mModel
### Keywords: models htest

### ** Examples

## ## ## testadd
## ## ## 

## ## Discrete model
## ## 
data(reinis)
## A decomposable model
##
mf <- ~smoke:phys:mental+smoke:systol:mental
object <- dmod(mf, data=reinis)
testadd(object,c("systol","phys"))


## A non-decomposable model
##
mf <- ~smoke:phys+phys:mental+smoke:systol+systol:mental
object <- dmod(mf, data=reinis)
testadd(object,c("phys","systol"))


## ## Continuous model
## ## 
data(math)
## A decomposable model
##
mf <- ~me:ve:al+al:an
object <- cmod(mf, data=math)
testadd(object,c("me","an"))

## A non-decomposable model
##
mf <- ~me:ve+ve:al+al:an+an:me
object <- cmod(mf, data=math)
testadd(object,c("me","al"))



cleanEx()
nameEx("testdelete")
### * testdelete

flush(stderr()); flush(stdout())

### Name: testdelete
### Title: Test deletion of edge from an interaction model
### Aliases: testdelete testdelete.iModel print.testdelete
###   testdelete.mModel
### Keywords: models htest

### ** Examples

## ## ## testdelete
## ## ## 

## ## Discrete model
## ## 
data(reinis)
## A decomposable model
##
mf <- ~smoke:phys:mental+smoke:systol:mental
object <- dmod(mf, data=reinis)

testdelete(object,c("phys","mental"))
testdelete(object,c("smoke","mental"))
#testdelete(object,c("systol","phys"))


## A non-decomposable model
##
mf <- ~smoke:phys+phys:mental+smoke:systol+systol:mental
object <- dmod(mf, data=reinis)

testdelete(object,c("phys","mental"))
#testdelete(object,c("systol","phys"))
#testdelete(object,c("smoke","mental"))


## ## Continuous model
## ## 
data(math)
## A decomposable model
##
mf <- ~me:ve:al+me:al:an
object <- cmod(mf, data=math)

testdelete(object,c("ve","al"))
testdelete(object,c("me","al"))




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
