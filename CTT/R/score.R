`score` <-
function(items,key,output.scored=FALSE,ID=NA,rel=FALSE){ 
  t<- as.vector(ID)                                          
  t<- table(ID)  
  if(any(t>1)){ for(i in 1:length(ID)){
                   for(j in 1:nrow(subset(t,t>1))){
                   if(ID[i]==(rownames(subset(t,t>1)))[j])
                   {ID[i]<- paste(ID[i],"/",i)}}}
        warning("Duplicate ID exists; the duplicate ID has been renamed and retained in the calculation")
               }
  
   if(!missing(ID)){
     if(length(ID)==nrow(items)) rownames(items) <- ID
        else warning("The length of ID vector does not match the sample size.")}
               
   if(missing(key)){
    warning("No key provided, assuming pre-scored data.")
	
    scored <- apply(items,2, function(XXX){
                                	if(! is.numeric(XXX)) XXX <- as.numeric(XXX)
									XXX
									})
   } 
  else {
    if(length(key)==ncol(items)) scored <- t(apply(items,1,function(X){ifelse(X==(key),1,0)}))
    else stop("Number of items is not equal to the length of key.")
  }
 scores <- rowSums(scored)
 names(scores)<-paste("P",c(seq(1:nrow(items))),sep="")
 if(!rel==FALSE)reli<-reliability(scored)
 if(output.scored==FALSE & rel==FALSE) out<-list(score=scores)
 if(output.scored==FALSE & rel==TRUE)out<-list(score=scores,reliability=reli) 
 if(output.scored==TRUE & rel==FALSE)out<-list(score=scores,scored=scored)
 if(output.scored==TRUE & rel==TRUE) out<- list(score=scores,reliability=reli, scored=scored)
 out
}

