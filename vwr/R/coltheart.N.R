coltheart.N <-
function(sources, targets, distance=1, method='hamming', parallel=FALSE){
    if(typeof(sources)!='character'){
        sources<-as.character(sources)
    }
    distance.function<-which.distance.function(method)
    do.one<-function(xsource, targets) return(sum(distance.function(xsource, targets)==distance, na.rm=TRUE))
    
    if (parallel==TRUE){
        ncores=detectCores(logical = FALSE)
        results<-unlist(mclapply(sources, do.one, targets, mc.cores=ncores))
        names(results)<-sources
    }
    else{
        results<-unlist(lapply(sources,do.one,targets))
        names(results)<-sources
    }
    return(results)
}

