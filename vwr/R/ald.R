ald <-
function(sources, targets, n, method='levenshtein', parallel=FALSE){
  
  if (method=="levenshtein"){
    distance.function=levenshtein.distance
  }
  if (method=="levenshtein.damerau"){
    distance.function=levenshtein.damerau.distance
  }
  
  if (typeof(sources)!='character'){
        sources<-as.character(sources)
    }

    if (typeof(targets)!='character'){
        targets<-as.character(targets)
    }
    
    target.lengths<-nchar(targets)
        
    do.one<-function(source, targets, target.lengths, n){
        distances<-rep(NA,length(targets))    
        min.distances<-abs(nchar(source)-target.lengths)
        unique.distances<-sort(unique(min.distances))
        for(current.distance in unique.distances){
            indexes<-which(min.distances==current.distance)
            distances[indexes]<-distance.function(source,targets[indexes])
            dfreq<-tabulate(distances)
            cutoff.index<-max(which(cumsum(dfreq)<n)+1,1)
            if (cutoff.index<current.distance+1){
                dfreq<-dfreq[1:cutoff.index]
                break
            }
        }
        return(mean(rep(seq_along(dfreq),dfreq)[1:n]))
    }
    if (parallel==TRUE){
        ncores = detectCores(logical = FALSE)
        results<-unlist(mclapply(sources, do.one, targets, target.lengths, n, mc.cores=ncores))
        names(results)<-sources
        }
    else{
        results<-sapply(sources, do.one, targets, target.lengths, n)
    }
    return(results)
}

