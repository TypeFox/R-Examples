index.G2 <- function(d,cl){
  cn <- max(cl)
  n <- length(cl)
  dmat <- as.matrix(d)
  diameter <- average.distance <- median.distance <- separation <-
    average.toother <- 
    cluster.size <- within.dist <- between.dist <- numeric(0)
  separation.matrix <- matrix(0,ncol=cn,nrow=cn)
  di <- list()
  for (i in 1:cn){
    cluster.size[i] <- sum(cl==i)
    #print(i)
    #print(cl==i)
    #print(dmat[cl==i,cl==i])
    di <- as.dist(dmat[cl==i,cl==i])
    within.dist <- c(within.dist,di)
    #diameter[i] <- max(di)
    average.distance[i] <- mean(di)
    median.distance[i] <- median(di)
    bv <- numeric(0)
    for (j in 1:cn){
      if (j!=i){
        sij <- dmat[cl==i,cl==j]
        bv <- c(bv,sij)
        if (i<j){
          separation.matrix[i,j] <- separation.matrix[j,i] <- min(sij)
          between.dist <- c(between.dist,sij)
        }
      }
    }   
   }
   nwithin<-length(within.dist)
   nbetween<-length(between.dist)
   .C("fng2",as.double(within.dist),as.integer(nwithin),as.double(between.dist),as.integer(nbetween),wynik=double(1),PACKAGE="clusterSim")$wynik[1]
}
