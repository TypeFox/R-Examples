#'@export

prune.DStree <- function(tree,data,gamma=2,which,...){
  
  fitr <- tree
  class(fitr) <-"rpart"
  
  if (missing(which)){
    which <- "BRIER"
  }
  if(missing(data)){
  time <- as.numeric(fitr$y[,1])
  states <- fitr$y[,2]
  }else{
  time <- as.numeric(data[,tree$names[1]])
  states <- data[,tree$names[2]]+1
  }
  
  npred <- length(time)
  z <- sort(unique(time))
  
  
  cptable<-fitr$cptable[-1,]
  cp<-cptable[,1]
  n.cp <- length(cp)
  n.lev <- length(fitr$parms[[1]]) 
  
  states2 <- abs(states-3)
  KM<-survfit(Surv(as.numeric(time),states2)~1,data.frame(time,states2))$surv
  DEV<-CRIT<-BRIER<-rep(NA,n.cp)
  fitr$frame$yval <- as.integer(rownames(fitr$frame))
  z <- min(time):round(median(time))
  if(min(z)==0){
    z<-z+1
    time<-time+1
  }
  
  
  for(i in 1:n.cp){
    prunedfitr <- prune.rpart(fitr,cp=cp[i])
    Pred <- predict(prunedfitr,data,type="matrix")
    nodes <- predict(prunedfitr,data,type="vector")
    N <- by(cbind(time,states),nodes,FUN=computeN,lev=fitr$parms[[1]])
    unique.nodes <- attributes(N)$dimnames$INDICES
    S <- subset(prunedfitr$frame,var=="<leaf>")
    S.ord <- S[match(unique.nodes,rownames(S)),]$yval2
    k <-length(unique.nodes)
    lik <- rep(0,k)
    
    for (j in 1:k){
      lik[j] <- lik(ncens=unlist(N[[j]][1]),
                                nuncens = unlist(N[[j]][2]),pi=S.ord[j,1:(n.lev)],S=S.ord[j,(n.lev+1):(2*n.lev)])    
    }
    
    
    CRIT[i] <- -2*sum(lik)+gamma*k*n.lev
    DEV[i] <- -2*sum(lik)
   
    BRIER[i]<- sum(sapply(z,brier,S=Pred[,(n.lev+1):(2*n.lev)],
                time=time,states=states,KM=KM,npred=npred))/length(z)
    
    M <- cbind(CRIT,DEV,BRIER)
  }
  ind <- which.min(M[,which])
  fitr$frame$yval <- tree$frame$yval
  prunedfitr <- prune.rpart(fitr,cp=cp[ind])
  prunedfit <- prunedfitr
  class(prunedfit)<-"DStree"
  
  return(list(nsplit=as.numeric(cptable[,2]),CRIT=CRIT,DEV=DEV,BRIER=BRIER,prunedfit=prunedfit))
  
}


