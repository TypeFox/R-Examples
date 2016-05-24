library(msr)

x =  cbind(rep(0:10, 11), rep(0:10, each=11))/10-0.55
y =  exp( - rowSums(x^2))

#build Morse-Smale complex
ms <- msc.nn(y=y, x=x, pLevel=0.1, knn = 15)


#check partition assignments
for(pId in 1:length(ms$level[[1]]$partitionSize)){
   tmp <- x[ms$partition == pId, ]
   for(d in 1:2){
     t1 <- sum(tmp[,d] < 0 )== 0
     t2 <- sum(tmp[,d] < 0) == ms$level[[1]]$partitionSize[pId]  
     if( (t1 || t2) == F){
       stop("Error in Morse-Smale computation")
     }
   }
}


