ldknn <-
function(stimuli, types, reference, k=1, method='levenshtein', parallel=FALSE){
    distance.function<-which.distance.function(method)
    n<-length(stimuli)
    do.one<-function(index, stimuli, k){
        distances<-distance.function(stimuli[index],stimuli[1:index-1])
        unique.distances<-sort(unique(distances))
        minimum.distance<-unique.distances[k]
        indexes<-which(distances<=minimum.distance)
        distribution<-types[1:index-1][indexes]
        probability<-sum(distribution==reference)/length(distribution)
        return(probability)
    }
    if(parallel==TRUE){
        ncores=detectCores(logical = FALSE)
        probabilities<-unlist(mclapply(2:n, do.one, stimuli, k, mc.cores=ncores))
    }
    else{
        probabilities<-unlist(lapply(2:n, do.one, stimuli, k))
    }
    probabilities<-c(0.5,probabilities)
    # print(probabilities)
    d<-data.frame(stimulus=stimuli,type=types,p=probabilities)
    odds<-ldknn.odds(d$type,d$p,reference)
    results<-list(data=d, reference.level=reference, k=k, odds=odds)
    class(results)<-'ldknn.run'
    return(results)
}

