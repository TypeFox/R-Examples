.ParKml.validity <- function(object){
#    cat("**** validity ParKml <empty> ****\n")
    return(TRUE)
}

setClass(
    Class="ParKml",
    representation=representation(
        saveFreq="numeric",
        maxIt="numeric",
        imputationMethod="character",
        distanceName="character",
        power="numeric",
        distance="function",
        centerMethod="function",
        startingCond="character",
 #       distanceStartingCond="function",
        nbCriterion="numeric",
        scale="logical"
    ),
    prototype=prototype(
        saveFreq=numeric(),
        maxIt=numeric(),
        imputationMethod=character(),
        distanceName=character(),
        power=numeric(),
        distance=function(){},
        centerMethod=function(){},
        startingCond=character(),
#        distanceStartingCond=function(){},
        nbCriterion=numeric(),
        scale=logical()
    ),
    validity=.ParKml.validity
)

parKml <- function(saveFreq,maxIt,imputationMethod,distanceName,power,distance,
                   centerMethod,startingCond,nbCriterion,scale){
    if(distanceName %in%DISTANCE_METHODS){
        eval(parse(text=paste("distance <- function(x,y){dist(rbind(x,y),method='",distanceName,"',p=",power,")}",sep="")))
    }else{}
    new("ParKml",saveFreq=saveFreq,maxIt=maxIt,imputationMethod=imputationMethod,
        distanceName=distanceName,power=power,distance=distance,
        centerMethod=centerMethod,startingCond=startingCond,#distanceStartingCond=distanceStartingCond,
        nbCriterion=nbCriterion,scale=scale)
}

parALGO <- function(saveFreq=100,maxIt=200,imputationMethod="copyMean",
                   distanceName="euclidean",power=2,distance=function(){},
                   centerMethod=meanNA,startingCond="nearlyAll",#distanceStartingCond=function(x,y)dist(rbind(x,y)),
                   nbCriterion=1000,scale=TRUE){
    if(distanceName %in%DISTANCE_METHODS){
        eval(parse(text=paste("distance <- function(x,y){dist(rbind(x,y),method='",distanceName,"',p=",power,")}",sep="")))
    }else{}
    new("ParKml",saveFreq=saveFreq,maxIt=maxIt,imputationMethod=imputationMethod,
        distanceName=distanceName,power=power,distance=distance,
        centerMethod=centerMethod,startingCond=startingCond,#distanceStartingCond=distanceStartingCond,
        nbCriterion=nbCriterion,scale=scale)
}


setMethod("[","ParKml",
    function(x,i,j,drop){
        switch(EXPR=i,
            "saveFreq"={return(x@saveFreq)},
            "maxIt"={return(x@maxIt)},
            "imputationMethod"={return(x@imputationMethod)},
            "distanceName"={return(x@distanceName)},
            "power"={return(x@power)},
            "distance"={return(x@distance)},
            "centerMethod"={return(x@centerMethod)},
            "startingCond"={return(x@startingCond)},
#            "distanceStartingCond"={return(x@distanceStartingCond)},
            "nbCriterion"={return(x@nbCriterion)},
            "scale"={return(x@scale)},
            stop("[ParKml:get]: there is not such a slot in ParWindows")
        )
    }
)

setMethod(f="[<-",signature="ParKml",
    definition=function(x,i,j,value){
        switch(EXPR=i,
            "saveFreq"={x@saveFreq <- value},
            "maxIt"={x@maxIt <- value},
            "imputationMethod"={x@imputationMethod <- value},
            "distanceName"={
                x@distanceName <- value
                if(x@distanceName %in%DISTANCE_METHODS){
                    eval(parse(text=paste("x@distance <- function(x,y){dist(rbind(x,y),method='",x@distanceName,"',p=",x@power,")}",sep="")))
                }else{}
            },
            "power"={x@power <- value},
            "distance"={
                x@distance <- value
                x@distanceName <- "User defined"
            },
            "centerMethod"={x@centerMethod <- value},
            "startingCond"={x@startingCond <- value},
#            "distanceStartingCond"={x@distanceStartingCond <- value},
            "nbCriterion"={x@nbCriterion <- value},
            "scale"={x@scale <- value},
            stop("[ParKml:set]: there is not such a slot in ParWindows")
        )
        validObject(x)
        return(x)
    }
)


cat("### Method : 'show' for ParKml ###\n")
.ParKml.show <- function(object){
    cat("   ~~~ Class: ParKml ~~~ ")
    cat("\n ~ saveFreq             :",object@saveFreq)
    cat("\n ~ maxIt                :",object@maxIt)
    cat("\n ~ imputationMethod     :",object@imputationMethod)
    cat("\n ~ distanceName         :",object@distanceName)
    cat("\n ~ power                :",object@power)
    cat("\n ~ distance             : ");print(object@distance)
    cat(" ~ centerMethod         : ");print(object@centerMethod)
    cat(" ~ startingCond         :",object@startingCond)
#    cat("\n ~ distanceStartingCond : ");print(object@distanceStartingCond)
    cat("\n ~ nbCriterion          :",object@nbCriterion,"\n")
    cat("\n ~ scale                :",object@scale,"\n")
    return(invisible(object))
}
setMethod(f="show",signature="ParKml",definition=.ParKml.show)
