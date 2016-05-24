levenshtein.damerau.distance <-
function(xsource, targets){
    distances<-stringdist(xsource, targets, method='osa')
    names(distances)<-targets
    return(distances)
}