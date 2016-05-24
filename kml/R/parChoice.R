cat("\n####################################################################
########################## Class parChoice #########################
############################# Creation #############################
####################################################################\n")

.ParChoice.validity <- function(object){
#    cat("**** validity ParChoice <empty> ****")
    return(TRUE)
}


setClass(
    Class="ParChoice",
    representation=representation(
         toDo="character",
         xy="numeric",
         nbTime="numeric",
         critMatrix="matrix",
         selectedPart="list",
         toPlot="character",
         styleTrajRank="numeric",
         styleMeanRank="numeric",
         critRank="numeric",
         critSorted="logical",
         cex="numeric",
         pchPeriod="numeric",
         yLegend="numeric"
    ),
    prototype=prototype(
         toDo=character(),
         xy=numeric(),
         nbTime=numeric(),
         critMatrix=matrix(),
         selectedPart=list(),
         toPlot=character(),
         styleTrajRank=numeric(),
         styleMeanRank=numeric(),
         critRank=numeric(),
         critSorted=logical(),
         cex=numeric(),
         pchPeriod=numeric(),
         yLegend=numeric()
    ),
    validity=.ParChoice.validity
)


cat("### Constructor ###\n")
parChoice <- function(xy,nbTime,critMatrix,selectedPart=list(),toPlot="both",styleTrajRank=1,styleMeanRank=1,
                      critRank=1,critSorted=TRUE,cex=1.2,pchPeriod=1,yLegend=-0.12){
    return(new("ParChoice",toDo="",xy=xy,nbTime=nbTime,critMatrix=critMatrix,selectedPart=selectedPart,toPlot=toPlot,styleTrajRank=styleTrajRank,styleMeanRank=styleMeanRank,
        critRank=critRank,critSorted=critSorted,cex=cex,pchPeriod=pchPeriod,yLegend=yLegend))
}



cat("### Show ###\n")
.ParChoice.show <- function(object){
     cat("   ~~~ Class: ParChoice ~~~ ")
     cat("\n ~ toDo          : ",object@toDo)
     cat("\n ~ xy            : ",object@xy)
     cat("\n ~ nbTime        : ",object@nbTime)
     cat("\n ~ critMatrix\n");print(object@critMatrix)
     cat("\n ~ selectedPart  : ");print(object@selectedPart)
     cat(" ~ toPlot        : ",object@toPlot)
     cat("\n ~ styleTrajRank : ",object@styleTrajRank)
     cat("\n ~ styleMeanRank : ",object@styleMeanRank)
     cat("\n ~ critRank      : ",object@critRank)
     cat("\n ~ critSorted    : ",object@critSorted)
     cat("\n ~ cex           : ",object@cex)
     cat("\n ~ pchPeriod     : ",object@pchPeriod)
     cat("\n ~ yLegend       : ",object@yLegend)
     cat("\n")
    return(invisible())
}
setMethod(f="show",signature="ParChoice",definition=.ParChoice.show)


cat("### Getteur ###\n")
.ParChoice.get <- function(x,i,j,drop){
    switch(EXPR=i,
        "toDo"={return(x@toDo)},
        "xy"={return(x@xy)},
        "nbTime"={return(x@nbTime)},
        "critMatrix"={return(x@critMatrix)},
        "selectedPart"={return(x@selectedPart)},
        "toPlot"={return(x@toPlot)},
        "styleTrajRank"={return(x@styleTrajRank)},
        "styleMeanRank"={return(x@styleMeanRank)},
#        "critPossible"={return(x@critPossible)},
        "critRank"={return(x@critRank)},
        "critSorted"={return(x@critSorted)},
        "cex"={return(x@cex)},
        "pchPeriod"={return(x@pchPeriod)},
        "yLegend"={return(x@yLegend)},
        stop("[ParChoice:get] ",i," is not a 'ParChoice' slot")
    )
    return(invisible())
}
setMethod(f="[",signature="ParChoice",definition=.ParChoice.get)

cat("### Setteur ###\n")
 .ParChoice.set <- function(x,i,j,value){
    switch(EXPR=i,
        "toDo"={x@toDo<-value},
        "xy"={x@xy<-value},
        "nbTime"={x@xy<-value},
        "critMatrix"={x@critMatrix<-value},
        "selectedPart"={x@selectedPart<-value},
        "toPlot"={x@toPlot<-value},
        "styleTrajRank"={x@styleTrajRank<-value},
        "styleMeanRank"={x@styleMeanRank<-value},
 #       "critPossible"={x@critPossible<-value},
        "critRank"={x@critRank<-value},
        "critSorted"={x@critSorted<-value},
        "cex"={x@cex<-value},
        "pchPeriod"={x@pchPeriod<-value},
        "yLegend"={x@yLegend<-value},
        stop("[ParChoice:set] ",i," is not a 'ParChoice' slot")
     )
    validObject(x)
    return(x)
}
setMethod(f="[<-",signature="ParChoice",definition=.ParChoice.set)

cat("\n--------------------------------------------------------------------
------------------------ Fin Test ParChoice ------------------------
--------------------------------------------------------------------\n")
