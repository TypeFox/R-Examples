## the canonical gamma MME -- but the cMLE is lower variance and just as fast!
gammaMME <- function(x) { 
 return(c(shape=(mean(as.matrix(x),na.rm=T)/apply(as.matrix(x),2,sd,na.rm=T))^2,
          scale=apply(as.matrix(x),2,var,na.rm=T)/mean(as.matrix(x),na.rm=T)))
} 
