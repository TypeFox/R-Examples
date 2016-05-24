 # # # # # # # # # # # # # # # # # #
# # #     Class  definition     # # #
 # # #       ParKmlShape       # # #
  # # # # # # # # # # # # # # # # #

.ParKmlShape.validity <- function(object){
    return(TRUE)
}


setClass(
    Class="ParKmlShape",
    representation=representation(
        aggregationMethod="character",
        shuffle="logical",
        sampleSize="integer",
        methodHclust="character",
        maxIter="integer"
   ),
    prototype=prototype(
        aggregationMethod=character(),
        shuffle=logical(),
        sampleSize=integer(),
        methodHclust=character(),
        maxIter=integer()
   ),
    validity=.ParKmlShape.validity
)


cat("### Constructor ###\n")
parKmlShape <- function(aggregationMethod="all",shuffle=TRUE,sampleSize=128,methodHclust="average",maxIter=100){
    return(new("ParKmlShape",aggregationMethod=aggregationMethod,shuffle=shuffle,sampleSize=as.integer(sampleSize),methodHclust=methodHclust,maxIter=as.integer(maxIter)))
}


cat("### Show ###\n")
.ParKmlShape.show <- function(object){
     cat("   ~~~ Class:",class(object),"~~~ ")
    cat("\n ~ aggregationMethod : ",object@aggregationMethod)
    cat("\n ~ shuffle : ",object@shuffle)
    cat("\n ~ sampleSize : ",object@sampleSize)
    cat("\n ~ methodHclust : ",object@methodHclust)
    cat("\n ~ maxIter : ",object@maxIter)
    cat("\n")
    return(invisible())
}
setMethod(f="show",signature="ParKmlShape",definition=.ParKmlShape.show)


cat("### Getteur ###\n")
.ParKmlShape.get <- function(x,i,j,drop){
    switch(EXPR=i,
       "aggregationMethod"={return(x@aggregationMethod)},
       "shuffle"={return(x@shuffle)},
       "sampleSize"={return(x@sampleSize)},
       "methodHclust"={return(x@methodHclust)},
       "maxIter"={return(x@maxIter)},
       stop("[ParKmlShape:get] ",i," is not a 'ParKmlShape' slot")
    )
    return(invisible())
}
setMethod(f="[",signature="ParKmlShape",definition=.ParKmlShape.get)


cat("### Setteur ###\n")
.ParKmlShape.set <- function(x,i,j,value){
    switch(EXPR=i,
       "aggregationMethod"={x@aggregationMethod<-value},
       "shuffle"={x@shuffle<-value},
       "sampleSize"={x@sampleSize<-value},
       "methodHclust"={x@methodHclust<-value},
       "maxIter"={x@maxIter<-value},
       stop("[ParKmlShape:set] ",i," is not a 'ParKmlShape' slot")
    )
    validObject(x)
    return(x)
}
setMethod(f="[<-",signature="ParKmlShape",definition=.ParKmlShape.set)
