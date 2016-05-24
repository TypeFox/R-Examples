A.dep.tests <- function(Xmat,choice=1,d=0,m=d,freqname="",type="text") {

  if (!(choice %in% 1:2)) stop("choice should be 1 or 2")
  
  Xd <- function(x,d) {
    n <- length(x)
    V1 <- NULL
    for (i in 1:d) V1 <- cbind(V1,x[i:(n-d+i)])
    V1
  }
  
  
  if (choice==1) {
    if (is.table(Xmat)) {
      d <- length(dim(Xmat))
      tmat <- sapply(dimnames(Xmat),length)
    } else {
      if (freqname=="") {
        d <- ncol(Xmat)
        tmat <- sapply(lapply(Xmat,MARGIN=2,FUN=unique),length)
        V1 <- Xmat
      } else {
        Xmat <- eval(parse(text=paste("xtabs(formula=",freqname,"~.,data=Xmat)")))
        d <- length(dim(Xmat))
        tmat <- sapply(dimnames(Xmat),length)
      }
    }
  }
  
  if (choice==2) {
    if (d==0) stop("You should enter a value for d")
    if (m>d) stop("You should enter a value for m which is <=d")
    V1 <- Xd(Xmat,d) # la fonction Xd est ci-dessus
    #Xmat <- table(as.data.frame(V1))
    tmat <- sapply(lapply(data.frame(V1),MARGIN=2,FUN=unique),length)
  }
  
                                        
  RES <- as.list(1:(m-1))
  TheTAs <- RES

  if (choice==1) {
  
    for (cardA in 2:m) {
      RES[[cardA-1]] <- as.matrix(combn(d,cardA))
      TheTAs[[cardA-1]] <- as.data.frame(t(as.matrix(RES[[cardA-1]][1,])))
      dimnames(RES[[cardA-1]]) <- list(NULL,apply(RES[[cardA-1]],2,f<-function(...){paste(...,collapse=",")}))
      names(TheTAs[[cardA-1]]) <- dimnames(RES[[cardA-1]])[[2]]
    }
  
    ThefAs <- TheTAs
    ThepvalAs <- TheTAs
  
    for (cardA in 2:m) {
      bs <- choose(d,cardA)
      for (j in 1:bs) {
        
        A <- RES[[cardA-1]][,j]
        if (is.table(Xmat)) X2A <- summary(margin.table(Xmat, A))$statistic else X2A <- summary(table(as.data.frame(V1[,A])))$statistic
        
        somme <- 0
        if (cardA>2) {
          for (cardB in 2:(cardA-1)) {
            somme <- somme + sum(TheTAs[[cardB-1]][list(NULL,apply(as.matrix(combn(A,cardB)),2,f<-function(...){paste(...,collapse=",")}))[[2]]])
          }
        }
        
        TA <- X2A-somme
        TheTAs[[cardA-1]][paste(as.character(A),collapse=",")] <- TA
        
        dA <- length(A)
        I <- rep(0,dA)
        
        for (k in 1:dA) I[k] <- tmat[A[k]]
        fA <- prod(as.vector(I-1))
        ThefAs[[cardA-1]][paste(as.character(A),collapse=",")] <- fA
        ThepvalAs[[cardA-1]][paste(as.character(A),collapse=",")] <- pchisq(TA,df=fA,lower.tail=FALSE)
        
      }
      
    }

  } else { # choice==2

    for (cardA in 2:m) {
      RES[[cardA-1]] <- rbind(rep(1,choose(d-1,cardA-1)),as.matrix(combn(2:d,cardA-1))) # tous les ensembles A avec 1
      TheTAs[[cardA-1]] <- as.data.frame(t(as.matrix(RES[[cardA-1]][1,])))
      dimnames(RES[[cardA-1]]) <- list(NULL,apply(RES[[cardA-1]],2,f<-function(...){paste(...,collapse=",")}))
      names(TheTAs[[cardA-1]]) <- dimnames(RES[[cardA-1]])[[2]]
    }
    
    ThefAs <- TheTAs
    ThepvalAs <- TheTAs
  
    for (cardA in 2:m) {
      bs <- choose(d-1,cardA-1)
      for (j in 1:bs) {
        
        A <- RES[[cardA-1]][,j]
        if (is.table(Xmat)) X2A <- summary(margin.table(Xmat, A))$statistic else X2A <- summary(table(as.data.frame(V1[,A])))$statistic
     
        somme <- 0
        X2Awo1 <- 0
        if (cardA>2) {
        if (is.table(Xmat)) X2Awo1 <- summary(margin.table(Xmat, A[-1]))$statistic else X2Awo1 <- summary(table(as.data.frame(V1[,A[-1]])))$statistic
          for (cardB in 2:(cardA-1)) {
            somme <- somme + sum(TheTAs[[cardB-1]][list(NULL,apply(rbind(rep(1,choose(cardA-1,cardB-1)),as.matrix(combn(A[-1],cardB-1))),2,f<-function(...){paste(...,collapse=",")}))[[2]]])
          }
        }

        TA <- X2A-X2Awo1-somme
        TheTAs[[cardA-1]][paste(as.character(A),collapse=",")] <- TA
        
        dA <- length(A)
        I <- rep(0,dA)
        
        for (k in 1:dA) I[k] <- tmat[A[k]]
        fA <- prod(as.vector(I-1))
        ThefAs[[cardA-1]][paste(as.character(A),collapse=",")] <- fA
        ThepvalAs[[cardA-1]][paste(as.character(A),collapse=",")] <- pchisq(TA,df=fA,lower.tail=FALSE)
        
      }
      
    }



  }
    

  names(TheTAs) <- paste("|A|=",2:m,sep="")
  names(ThefAs) <- names(TheTAs)
  names(ThepvalAs) <- names(TheTAs)
  
  res <- list(TA=TheTAs,fA=ThefAs,pvalA=ThepvalAs)
  
  
  X <- c()
  for (cardA in 2:m) {
    X <- rbind(X,eval(parse(text=paste("t(rbind(res$TA$`|A|=",cardA,"`,res$fA$`|A|=",cardA,"`,res$pvalA$`|A|=",cardA,"`))",sep=""))))
  }
  
  colnames(X) <- c("TA","fA","pvalA")
  
# Test of Mutual independence
  M <- sum(X[,1])
  fM <- sum(X[,2])
  pvalM <- pchisq(M,df=fM,lower.tail=FALSE)
  
  X2tmp <- rbind(c("","",""),round(X,2),c("","",""),round(c(M,fM,pvalM),2))
#  Anames <- sapply(strsplit(rownames(X),"",fixed=TRUE),FUN=function(x) paste(x,collapse=","))
  Anames <- rownames(X)
  if (choice==1) rownames(X2tmp) <- c("A",Anames,"","X2") else rownames(X2tmp) <- c("A",Anames,"","Y2")


  if (type=="html") {
    require(xtable)
    sortie <- cbind(subset=rownames(X),round(X,2))
    sortie <- rbind(sortie,c("",sortie[nrow(sortie),2:4]))
    if (choice==1) {sortie[nrow(sortie)-1,] <- c("<B>all</B>","<B>X<sup>2</sup></B>","<B>total</B>","<B>p.value</B>")}
    if (choice==2) {sortie[nrow(sortie)-1,] <- c("<B>all</B>","<B>Y<sup>2</sup></B>","<B>total</B>","<B>p.value</B>")}
    print(xtable(data.frame(sortie)),type="html",include.rownames=FALSE)
  }
  else {
    print(data.frame(X2tmp))
  }
  cat("\n")
  
  if (choice == 1) {
    return(invisible(c(res,list(X=X,X2=M,f=fM,pval=pvalM))))
  } else {
    return(invisible(c(res,list(X=X,Y2=M,f=fM,pval=pvalM))))
  }
  
}



