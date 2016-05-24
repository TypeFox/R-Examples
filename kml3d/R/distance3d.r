#####################################################
################### Distance Traj ###################
#####################################################
dist3d <- function(x,y,method="euclidian",power=2){
    return(dist(rbind(c(x),c(y)),method=method,p=power))
#    return(sqrt(sum((c(x)-c(y))^2)))
}

#dist3dCoef <- function(x,y,coef){
#    distance <- 0
#    for(i in 1:ncol(x)){
#        distance <- distance+coef[i]*dist(rbind(x[,i],y[,i]))
#    }
#    return(distance/sum(coef))
#    return(sqrt(sum((c(x)-c(y))^2)))
#}


matDist3d <- function(data,distance=dist3d){
   nbRow <- dim(data)[1]
   matDistance <- matrix(0,nbRow,nbRow)
   duree <- system.time(matDistance[1,2] <- matDistance[1,2] <- distance(data[1,,],data[2,,]))[1]*(nbRow^2)/2
   if(duree>10){cat("Approximative time of distance matrix computation: ",duree,"\n")}else{}
   for(i in 1:(nbRow-1)){
       for(j in (i+1):nbRow){
           matDistance[i,j] <- matDistance[j,i] <- distance(data[i,,],data[j,,])
       }
   }
   return(matDistance)
}
