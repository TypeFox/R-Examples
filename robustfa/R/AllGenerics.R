## show, print, summary, predict, plot were defined in packages.
## However, they are not S4 generic functions. After setGeneric, they are S4 generic functions.
## setGeneric("show") # standardGeneric for "show" defined from package "methods"
setGeneric("print") # standardGeneric for "print" defined from package "base"
setGeneric("summary") # standardGeneric for "summary" defined from package "base"
setGeneric("predict") # standardGeneric for "predict" defined from package "stats"
setGeneric("plot") # standardGeneric for "plot" defined from package "graphics"    

if(!isGeneric("getQuan"))
	setGeneric("getQuan", function(obj) standardGeneric("getQuan"))
## returns the number of observations used in the computation of the FA (n for classic)

if(!isGeneric("getCenter"))
	setGeneric("getCenter", function(obj) standardGeneric("getCenter")) 

if(!isGeneric("getLoadings"))
	setGeneric("getLoadings", function(obj) standardGeneric("getLoadings")) 

if(!isGeneric("getEigenvalues"))
	setGeneric("getEigenvalues", function(obj) standardGeneric("getEigenvalues")) 

if(!isGeneric("getSdev"))
	setGeneric("getSdev", function(obj) standardGeneric("getSdev")) 

if(!isGeneric("getScores"))
	setGeneric("getScores", function(obj) standardGeneric("getScores")) 

setGeneric("getFa", function(obj) standardGeneric("getFa")) 
## return a factanal() compatible object to use the available standard plots 
## (i.e. a list with loadings, uniquenesses, covariance, correlation, usedMatrix, scores, scoresMethod, sdev)