cal.Dvine <- function(obj,val) {
  if(!class(obj$Dvine)=="penDvine") stop("obj$Dvine is not from class penDvine")
  Dvine.save <- obj$Dvine
  doMC <- TRUE
  registerDoParallel()
  K <- obj$K
  S <- obj$S
  N <-  obj$N
  base <- obj$base
                      
  Dvine <- list()
  SSi <-  list()
  for ( i in 1:length(S))
    {
      vine.knot <-  list(j1=S[i],j2=NULL,D=NULL, v=NULL,U=val[,S[i]] )
      SSi <- append(SSi, list(vine.knot))
    }
  Dvine <- append(Dvine, list(SSi))

  level <-  1
  SS <-  Dvine[[level]]
  nSS <- length(SS)
  while (nSS>1)
    {
      SSi <-  list()
      for (i in 2:nSS)
        {
          S1 <- c(SS[[i-1]]$j1, SS[[i-1]]$j2,SS[[i-1]]$D)
          S2 <- c(SS[[i]]$j1, SS[[i]]$j2,SS[[i]]$D)
          index1 <- rep(TRUE,length(S1))
          index2 <- rep(TRUE,length(S2))
          S3 <- c()
          for (j in 1:length(S1))
            {
              indexi <- S1[j]==S2
              if (sum(indexi)>0)
                {
                  S3 <- c(S3, S1[j])  
                  index1[j] <- FALSE
                  index2[indexi] <- FALSE                
                 }
            }
          if(!is.null(S3)) S3 <- sort(S3)
          S1 <- S1[index1]
          S2 <- S2[index2]
          vine.knot <-  list(j1=S1,j2=S2,D=S3,U=NULL)
          SSi <- append(SSi, list(vine.knot))
        }
      Dvine <- append(Dvine, list(SSi))
      level <- level + 1
    SS <-  Dvine[[level]]
      nSS <- length(SS)
    }
  
  level <-  2
  Tree.l <- Dvine[[level]]
  
  Index.basis.D <- matrix(NA,(K+1)^2,2)
  Index.basis.D[,1] <- rep(seq(1,(K+1)),(K+1))
  Index.basis.D[,2] <- sort(Index.basis.D[,1])
 
log.like <- AIC <- 0
  
if(doMC){
  Tree.l.temp <- foreach(i=1:length(Tree.l),.combine=list, .multicombine=TRUE) %dopar%   {
    index.j1 <- Tree.l[[i]]$j1
    index.j2 <- Tree.l[[i]]$j2
    U.hat <- val[,c(index.j1,index.j2)]
    #U <- matrix(plot.paircopula(Dvine.save[[level]][[i]]$cop,val=U.hat,int=FALSE,marg=FALSE,contour=FALSE,plot=FALSE)[,3],ncol=1)
    U <- eval.paircopula(val[,c(index.j1,index.j2)],K=K,int=FALSE,Index.basis.D,ck.val=Dvine.save[[level]][[i]]$v,base)
    ll <- sum(sapply(U,log))
    list(j1=index.j1,j2=index.j2,D=NULL,v=Dvine.save[[level]][[i]]$v,U.hat=U.hat,U=U,cop=Dvine.save[[level]][[i]]$cop,log.like=ll) 
  }
}
else {
  Tree.l.temp <- list()
  for(i in 1:length(Tree.l)) { 
    index.j1 <- Tree.l[[i]]$j1
    index.j2 <- Tree.l[[i]]$j2
    U.hat <- val[,c(index.j1,index.j2)]
    U <- eval.paircopula(val[,c(index.j1,index.j2)],K=K,int=FALSE,Index.basis.D,ck.val=Dvine.save[[level]][[i]]$v,base)
    ll <- sum(sapply(U,log))
    Tree.l.temp[[i]] <- list(j1=index.j1,j2=index.j2,D=NULL,v=Dvine.save[[level]][[i]]$v,U.hat=U.hat,U=U,cop=Dvine.save[[level]][[i]]$cop,log.like=ll) 
  }
}

  for(j in 1:length(Tree.l.temp)) {
    log.like <- log.like+Tree.l.temp[[j]]$log.like
  }

  Dvine[[level]] <- Tree.l.temp

  while( length(Dvine[[level]])>2) # loop over the levels of the Dvine. We stop if the tree in the Dvine has just one knot
  {
    level <-  level+1
    Tree.l <- Dvine[[level]] # current tree in vine
    Tree.l1 <- Dvine[[level-1]] # previous tree in vine, one level up
    if(doMC){
      Tree.l.temp <- foreach(i=1:length(Tree.l),.combine=list, .multicombine=TRUE) %dopar% {
        U.hat <- c()
        index <- list(c(Tree.l[[i]]$j1,Tree.l[[i]]$D),c(Tree.l[[i]]$j2,Tree.l[[i]]$D))
        for(j in 1:2)
          {
            indexi <- index[[j]]
            index.ancestor <- c()
            for (m in 1:length(Tree.l1))
              {
                index.ancestor <- c(index.ancestor, all( sort(indexi)==sort(c(Tree.l1[[m]]$j1,Tree.l1[[m]]$j2,Tree.l1[[m]]$D))))
              }
            ancestor.knot <- Tree.l1[index.ancestor][[1]] # Ancestor knot of j-th Element in Index set, i.e. knot with indices {j1,D}  # for j=1 and {j2,D} for j=2, respectively.
            if(j==1) U.hat <-cbind(U.hat,cond.cop(data=ancestor.knot$U.hat,coef=ancestor.knot$v,K=K,diff="u2",Index.basis.D,base))
            if(j==2) U.hat <-cbind(U.hat,cond.cop(data=ancestor.knot$U.hat,coef=ancestor.knot$v,K=K,diff="u1",Index.basis.D,base))          
          }
        #U  <- matrix(plot.paircopula(Dvine.save[[level]][[i]]$cop,val=U.hat,int=FALSE,marg=FALSE,contour=FALSE,plot=FALSE)[,3],ncol=1)
        U <- eval.paircopula(val=U.hat,K=K,int=FALSE,Index.basis.D,ck.val=Dvine.save[[level]][[i]]$v,base)
        ll <- sum(sapply(U,log))
        list(j1=Tree.l[[i]]$j1,j2=Tree.l[[i]]$j2,D=Tree.l[[i]]$D,U.hat=U.hat,U=U,v=Dvine.save[[level]][[i]]$v,cop=Dvine.save[[level]][[i]]$cop,log.like=ll)
      }
    }
  else {
    level <-  level+1
    Tree.l <- Dvine[[level]] # current tree in vine
    Tree.l1 <- Dvine[[level-1]] # previous tree in vine, one level up
    Tree.l.temp <- list()
    for (i in 1:length(Tree.l)) {
      U.hat <- c()
      index <- list(c(Tree.l[[i]]$j1,Tree.l[[i]]$D),c(Tree.l[[i]]$j2,Tree.l[[i]]$D))
      for(j in 1:2)
        {
          indexi <- index[[j]]
          index.ancestor <- c()
          for (m in 1:length(Tree.l1))
            {
              index.ancestor <- c(index.ancestor, all( sort(indexi)==sort(c(Tree.l1[[m]]$j1,Tree.l1[[m]]$j2,Tree.l1[[m]]$D))))
            }
          ancestor.knot <- Tree.l1[index.ancestor][[1]] # Ancestor knot of j-th Element in Index set, i.e. knot with indices {j1,D}
                                        # for j=1 and {j2,D} for j=2, respectively.
          if(j==1) U.hat <-cbind(U.hat,cond.cop(data=ancestor.knot$U.hat,coef=ancestor.knot$v,K=K,diff="u2",Index.basis.D,base))
          if(j==2) U.hat <-cbind(U.hat,cond.cop(data=ancestor.knot$U.hat,coef=ancestor.knot$v,K=K,diff="u1",Index.basis.D,base))          
          }
      U <- eval.paircopula(val=U.hat,K=K,int=FALSE,Index.basis.D,ck.val=Dvine.save[[level]][[i]]$v,base)
      ll <- sum(sapply(U,log))
      Tree.l.temp[[i]] <- list(j1=Tree.l[[i]]$j1,j2=Tree.l[[i]]$j2,D=Tree.l[[i]]$D,U.hat=U.hat,U=U,v=Dvine.save[[level]][[i]]$v,cop=Dvine.save[[level]][[i]]$cop,log.like=ll)
      }
    }
    
    if(length(Dvine[[level]])==1)  {
      Tree.l.temp <- list(Tree.l.temp)
    }
    Dvine[[level]] <- Tree.l.temp
    
    for(i in 1:length(Tree.l.temp)) {
      log.like <- log.like+Tree.l.temp[[i]]$log.like
    }
  }
  level <-  level+1
  Tree.l <- Dvine[[level]] # current tree in vine
  Tree.l1 <- Dvine[[level-1]] # previous tree in vine, one level up
  i <- 1
  U.hat <- c()
  index <- list(c(Tree.l[[i]]$j1,Tree.l[[i]]$D),c(Tree.l[[i]]$j2,Tree.l[[i]]$D))
  for(j in 1:2)
    {
      indexi <- index[[j]]
      index.ancestor <- c()
      for (m in 1:length(Tree.l1))
        {
          index.ancestor <- c(index.ancestor, all( sort(indexi)==sort(c(Tree.l1[[m]]$j1,Tree.l1[[m]]$j2,Tree.l1[[m]]$D))))
        }
      ancestor.knot <- Tree.l1[index.ancestor][[1]] # Ancestor knot of j-th Element in Index set, i.e. knot with indices {j1,D}
      if(j==1) U.hat <-cbind(U.hat,cond.cop(data=ancestor.knot$U.hat,coef=ancestor.knot$v,K=K,diff="u2",Index.basis.D,base))
      if(j==2) U.hat <-cbind(U.hat,cond.cop(data=ancestor.knot$U.hat,coef=ancestor.knot$v,K=K,diff="u1",Index.basis.D,base))
    }
  U <- eval.paircopula(val=U.hat,K=K,int=FALSE,Index.basis.D,ck.val=Dvine.save[[level]][[i]]$v,base)
  ll <- sum(sapply(U,log))
  Tree.l.temp <- list(j1=Tree.l[[i]]$j1,j2=Tree.l[[i]]$j2,D=Tree.l[[i]]$D,U=U,v=Dvine.save[[level]][[i]]$v,cop=Dvine.save[[level]][[i]]$cop,log.like=ll)
  Dvine[[level]][[1]] <- Tree.l.temp
  log.like <- log.like+ll
  U.result <- rep(1,dim(val)[1])
  for(i in 2:length(S))  {
    for(j in 1:length(Dvine[[i]]))  {
      U.result <- U.result*Dvine[[i]][[j]]$U
    }
  }
  return(list(cal=U.result,log.like=log.like))
}
