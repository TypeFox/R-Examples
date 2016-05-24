levenshtein.neighbors <-
function(xsource, targets){
    results<-list()
    distances<-levenshtein.distance(xsource, targets)
    for (distance in min(distances):max(distances)){
        results[distance]=list(names(which(distances==distance)))
    }
    return(results)
}

