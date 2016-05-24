# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Routines for storing values in 'spfrontier' call's specific environment
# Stack of 'spfrontier' calls is allowed
#


#
# Package-specific environment
#
.SpfrontierEnv <- new.env(hash = TRUE)

#
# Returns estimator-specific environment
#
getEstimatorEnvir <- function(){
    res <- NULL
    if (exists("estimatorEnvir", envir = .SpfrontierEnv)){
        res <- get("estimatorEnvir", envir = .SpfrontierEnv)
    }
    return(res)
}

#
# Creates a clone a given environment
#
cloneEnvir <- function(envir){
    return(as.environment(as.list(envir, all.names=TRUE)))
}

#
# Initializes a estimator-specific environment, putting all arguments into it
#
initEnvir <- function(...){
    envir <- getEstimatorEnvir()
    newenv <- new.env()
    
    # if an environment is already set (this is a nested call of estimator), 
    # save it for later restoration
    if (!is.null(envir)){
        en <- cloneEnvir(envir)
        assign("parentEnvir", en, envir = newenv)
    }
    argg <- list(...)
    for (a in names(argg)){
        if (!is.null(argg[[a]])) {
            assign(a, argg[[a]], envir=newenv)
        }
    }
    assign("estimatorEnvir", newenv, envir = .SpfrontierEnv)
}

#
# Finalizes a estimator-specific environment
#
finalizeEnvir <- function(){
    parentEnv <- envirGet("parentEnvir")
    
    # if a parent environment is already set (this is a nested call of estimator), 
    # restore it
    if(!is.null(parent.env)){
        assign("estimatorEnvir", parentEnv, envir = .SpfrontierEnv)
    }else{
        rm("estimatorEnvir",envir = .SpfrontierEnv)
    }
}

#
# Returns an object from the current estimator-specific environment
#
envirGet <- function(name,...){
    res <- NULL
    if (exists(name, envir = getEstimatorEnvir(),...)){
        res <- get(name, envir = getEstimatorEnvir(),...)
    }
    return(res)
}

#
# Assigns an object to a name in the current estimator-specific environment
#
envirAssign <- function(name,obj,...){
    assign(name, obj, envir = getEstimatorEnvir(),...)
}

#
# Counts a number of calls to a specified 'name' within the current estimator
#
envirCounter <- function(name,...){
    suffix <- ".count"
    cName <- paste(name,suffix,sep="")
    counter <- envirGet(cName,...)
    if (!is.null(counter)){
        counter <- counter + 1
    }else{
        counter <- 1
    }
    assign(cName, counter, envir = getEstimatorEnvir(),...)
    return(counter)
}