infiniteSearchStep <- function(oldResult=NULL, ..., start.designs) {
    if (is.null(oldResult)) {
        result <- searchCrossOverDesign(start.designs=start.designs, ...)
        return(result)
    }
    result <- searchCrossOverDesign(start.designs=list(getDesign(oldResult)), ...)
    return(combineResults(oldResult, result))
}

infiniteSearch <- function(...) {
    oldResult <- NULL
    assign(".search", TRUE, envir=Crossover.env)
    while(get(".search", envir=Crossover.env)) {
        result <- tryCatch({        
            infiniteSearchStep(oldResult=oldResult, ...)
        }, interrupt = function(interrupt) {            
            assign(".search", FALSE, envir=Crossover.env)    
        }, finally = {
            print("Another loop.")
        })
        if (!is.logical(result)) {
            if (is.null(oldResult)) {
                oldResult <- result
            } else {
                oldResult <- combineResults(oldResult, result)
            }
            print(result)
            gc()
        }
    }
    return(result)
}

combineResults <- function(x, y, save.history=FALSE) {    
    return(new("CrossoverSearchResult", 
               design=y@design, 
               startDesigns=ifelse(save.history, c(x@startDesigns, y@startDesigns), list()), 
               eff=ifelse(save.history, c(x@eff, y@eff), list()),
               search=x@search, 
               model=x@model, 
               time=c(x@time, y@time), 
               misc=ifelse(save.history, c(x@misc, y@misc), list())))
    
}

getMaxEffPerRun <- function(x) {
    return(unlist(lapply(x@eff, function(x) {max(x[!is.na(x)])})))
}

getMinEffPerRun <- function(x) {
    return(unlist(lapply(x@eff, function(x) {min(x[!is.na(x)])})))
}