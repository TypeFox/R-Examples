preprocessing <- function(data,stand,percAccomm,mahal=TRUE){
  
  if(stand == TRUE){
    data1 <- scale(data,center=sapply(data,mean),scale=sapply(data,sd))   
  }else{
    data1 <- data  
  }
  
  if(percAccomm != 1){
   if(mahal == TRUE){   
    Sx <- cov(data1)
    D2 <- mahalanobis(data1,colMeans(data1), Sx)
    indivYes <- which(D2 <= qchisq(percAccomm, df=dim(data1)[2]))
    indivNo <- which(D2 > qchisq(percAccomm, df=dim(data1)[2]))
    perc <- (length(indivYes) / dim(data1)[1]) * 100
    data1 <- data1[indivYes,]
   }else{
     if(ncol(data) <= 3){
      appr <- FALSE
     }else{
       appr <- TRUE
     }   
      dt = c()
      for(i in 1 : nrow (data1)){
       dt[i] <- depth(data1[i,], data1, approx=appr) 
      }
     num <- sum(dt == min(dt))
     indivYes <- which(dt != min(dt)) 
     indivNo <- which(dt == min(dt)) 
     perc <- (length(indivYes) / dim(data1)[1]) * 100
     data1 <- data1[indivYes,]
   }
  }else{
    data1 <- data1
  }

 if(percAccomm != 1){  
  print(paste("The percentage of accommodation is exactly ", round(perc,2), "%",sep="")) 
 }

 if(percAccomm != 1){  
  return(list(data=data1,indivYes=indivYes,indivNo=indivNo))
 }
 else{
  return(list(data=data1)) 
 }  
}