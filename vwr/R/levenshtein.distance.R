levenshtein.distance <-
function(xsource, targets){
    distances<-stringdist(xsource, targets, method='lv')
    names(distances)<-targets
    return(distances)
}

