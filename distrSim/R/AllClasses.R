################################
##
## Class: Dataclass
##
################################
setClass("SeqDataFrames", representation(data = "list"),
          prototype = list(data.frame(1)),
          validity = function(object){
               len <- length(object@data)
               if (len > 1)
                  { if (!all(unlist(lapply(object@data, is.data.frame))))
                        stop("all elements must be data frames")
                    f <- function(y) {list(ncol(y), names(y))}
                    g <- function(y) identical(f(y), f(object@data[[1]]))
                    if (!all(unlist(lapply(object@data, g))))
                        stop("all elements must have the same column structure")     
                  }
              return(TRUE) }
    )


################################
##
## Some Class Unions
##
################################


setClassUnion("ArrayorNULLorVectororDataframeorSeqDataFrames",c("array", "NULL",
               "vector", "data.frame", "SeqDataFrames"))
setClassUnion("DataframeorSeqDataFrames",c("data.frame", "SeqDataFrames"))
setClassUnion("ArrayorNULLorVector",c("array", "NULL","vector"))
setClassUnion("MatrixorNULLorVector",c("matrix", "NULL","vector"))

################################
##
## Class: Dataclass
##
################################


.pkgv <- as.character(
         if(all(is.na(packageDescription("distrSim"))))
             packageVersion("distr") else packageVersion("distrSim")
         )

setClass("Dataclass",
         representation(filename = "vectororNULL",
#old:                        Data = "vectororNULL",
                        Data = "ArrayorNULLorVectororDataframeorSeqDataFrames",
                        obsDim ="numeric",   ### new v.1.8
                        runs = "numeric",
                        samplesize = "numeric",
                        name = "character", ### new v.1.8
                        version = "character" ### new v.1.8
                        ),
         prototype=list(filename = "Data-set", Data = NULL, 
                        runs = 1, samplesize = 1, 
                        obsDim = 1, 
                        version = .pkgv,
                        name = "Data-Set"))

################################
##
## Class: Simulation
##
################################

### changed from version 1.8 on:
## ith observation in ith line of datamatrix/array
## jth item/dimension of each observation in jth column of datamatrix/array
## kth run/time of each observation in kth slide of datamatrix/array

## ++old
## +ith run in ith line of datamatrix
## +jth samples of each run in jth column of datamatrix


setClass("Simulation",
         representation("Dataclass",
                        seed = "list",
 ##new 03-10-06:
                        distribution = "Distribution"
 ###old:        distribution = "UnivariateDistribution"
         ),
         contains="Dataclass")


################################
##
## Class: Contsimulation
##
################################

setClass("Contsimulation",
         representation("Dataclass",
                        ind = "MatrixorNULLorVector",
                        rate = "numeric",
 ##new 03-10-06:
                        Data.id = "ArrayorNULLorVector",
                        Data.c =  "ArrayorNULLorVector",
                        distribution.c = "Distribution",
                        distribution.id = "Distribution",
 ###old:        Data.id = "vectororNULL",
 ###old:        Data.c = "vectororNULL",
 ###old:        distribution.c = "UnivariateDistribution"
 ###old:        distribution.id = "UnivariateDistribution"
                        seed = "list"),
         contains = "Dataclass")            
