tiger.cluster <- function(answer, som.init="sample", som.topol="hexa", maxc = 15, som.dim = c(20,20)){

    answer$som.dim <- som.dim
    answer$som.init <- som.init
    answer$som.topol <- som.topol
    answer$maxc <- maxc

    if(answer$use.som){
        if(length(answer$modelled) - answer$window.size < som.dim[1]*som.dim[2]){
            stop("SOM dimensions too large for length of time series")
        }

#        require(som)
        cat("Calculating SOM takes a while\n")
        mysom.rer <- som(as.matrix(answer$measures.uniform[!answer$na.rows,]),
            xdim=answer$som.dim[1], ydim=answer$som.dim[2],
            init=answer$som.init, topol=answer$som.topol)
        answer$som <- mysom.rer
        values <- mysom.rer$code
    } else {
        values <- answer$measures.uniform[!answer$na.rows,]
    }


    #Bestimmen der Anzahl Cluster
#    require(e1071)
    validity <- rep(NA, maxc)
    all.cluster <- list()
    for(centers in 2:maxc){
        cat("Result for", centers, "clusters\n")
        cluster<-cmeans(x=values, centers=centers, method="cmeans", m=2)
        #print(cluster)
        validity[centers] <- validityIndex(cluster , values)
        all.cluster[[centers]] <- cluster
    }
    
    answer$validityMeasure <- validity
    answer$cluster.assignment <- all.cluster

    if(answer$has.peaks){
       answer <- tiger.som.peaks(answer)
    }

    return(answer)

}
