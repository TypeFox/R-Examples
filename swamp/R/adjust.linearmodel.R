"adjust.linearmodel"<-function (g, o.batches, robust.lm=F, small.memory=F){              
       
       if(class(g)!="matrix"){stop("g is not a matrix")}
       if(nrow(data.frame(o.batches))!=ncol(g)){stop("o.batches has not the same number of samples as ncol(g)")}
       if(!is.data.frame(o.batches)){
       o.batches<-data.frame(o.batches,row.names=colnames(g))}
         naing<-sum(apply(g,2,function(x) sum(is.na(x))))
         gc()
         if (naing>0){stop("No NAs in g are allowed for adjust.linearmodel")}
    nas<-apply(o.batches,1,function(x) sum(is.na(x)))
    index<-which(nas==0)
    gafter<-matrix(nrow=nrow(g),ncol=ncol(g),dimnames=dimnames(g))
    if (length(which(nas>0))>0){gafter[,which(nas>0)]<-g[,which(nas>0)]
                                warning(paste("Samples",toString(which(nas>0)),"are not adjusted because of NAs in o.batches"))}    
    if(robust.lm==F){
      if(small.memory==F){
    fit1 <- lm(t(g) ~ ., o.batches) #full input, residuals come without NA samples
    gafter[,index] = t(fit1$residuals)} 
       else{
          for(i in 1:nrow(g)){
          gafter[i,index] <- lm(g[i,] ~ ., o.batches)$residuals} 
          }
       }           
          else {
          require(MASS)
          for(i in 1:nrow(g)){
          gafter[i,index] <- rlm(g[i,] ~ ., o.batches)$residuals}
          }
    rowm<-rowMeans(g[,index])      
    gc()
    gafter[,index]<-gafter[,index]+rowm
    gc()
    return(gafter)
} 
