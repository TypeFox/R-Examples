validityIndex<-function(cclust, values, verbose=FALSE){
         #the goal function of the classes (each member as close to
         #the center as possible
         mysum <- 0
         # The minimum distance between centers
         minDenom <-suppressWarnings(min())
         for(cluster in 1:length(cclust$size)){
              #Zähler aus Choi 2007
              mysum <- mysum + sum(cclust$membership[,cluster]*
                                     apply(values, 1, FUN=eD, y=cclust$centers[cluster,])
                                  )
              #Nenner aus Choi 2007 (alle mit allen vergleichen
              if(cluster <length(cclust$size)){
              for (cluster2 in (cluster+1):length(cclust$size)){
                   current <- eD(cclust$centers[cluster,] ,cclust$centers[cluster2,])
                   if(current < minDenom) minDenom <- current
              }
              }

              
         }
         if(verbose){
            cat(paste("Within cluster:", signif(mysum/NROW(values)), "between clusters:", signif(minDenom), sep="\t"), "\n")
         }
         return(mysum/NROW(values)/minDenom)
}

