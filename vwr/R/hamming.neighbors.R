hamming.neighbors <-
function(source, targets){
    results<-list()
    distances<-hamming.distance(source, targets)
    for (distance in sort(unique(distances))){
        results[distance]=list(names(which(distances==distance)))
    }
    return(results)
}

