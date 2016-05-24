levenshtein.damerau.neighbors <-
function(xsource, targets){
    results<-list()
    distances<-levenshtein.damerau.distance(xsource, targets)
    for (distance in sort(unique(distances))){
        results[distance]=list(names(which(distances==distance)))
    }
    return(results)
}

