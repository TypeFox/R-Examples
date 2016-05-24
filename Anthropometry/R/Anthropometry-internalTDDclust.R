assign("last.warning", NULL, envir = baseenv())

is.empty <- function(x, mode=NULL){
  if (is.null(mode)) mode <- class(x)
  identical(vector(mode,1),c(x,vector(class(x),1)))
}

rowweight <- function(X,a){ #Product function. 
  xm <- X * a 
  return(xm) 
} 

Distvec <- function(Xmatr,y){ #L1 distance from y to all data points. 
  Dv <- apply(Xmatr,1,l1norm,b=y) 
  return(Dv) 
} 

l1norm <- function(a,b){ #L1 distance from a to b. 
  l1 <- sqrt(sum((a-b)^2)) 
  return(l1) 
} 

unitvec <- function(a,b){ #Unit vector from a to b. 
  
  l1 <- sqrt(sum((a-b)^2)) 
  
  if (l1 > 10^(-10)){ 
    uv <- (b-a)/l1 
  }
  
  if (l1 < 10^(-10)){ 
    uv <- matrix(0,1,length(a)) 
  } 
  
  return(uv) 
} 


initW2 <- function(Xmatr,Nvec){ #Initial y for Weiszfeld algorithm. 
  #Xmatr = cell data matrix, Nmat = multiplicities. 
  dimX <- dim(Xmatr) 
  mm <- apply(Xmatr,2,mean) 
  ystart <- mm 
  return(ystart) 
} 


Weisziteradj <- function(Xmatr,Nvec,ystart,B){ #Iterates the Weiszfeld step B times # starting value ystart. 
  
  y <- ystart 
  
  for(k in (1:B)){ 
    W <- Weiszfeldadj(Xmatr,Nvec,y) 
    y <- W$y 
  } 
  
  DDy <- W$DDy 
  
  return(list(y=y,DDy=DDy))
}


Weiszfeldadj <- function(Xmatr,Nvec,ystart){ 
  #Xmatr distinct data matrix, Nvec vector of multiplicities, ystart is the starting point. 
  #Also gives data depth of new y. 
  l <- rbind(Xmatr,ystart)
  y <- ystart 
  
  lcorr <- cor(t(l)) 
  diag(lcorr) <- 0 
  uc <- upper.tri(lcorr) 
  lcorr[uc==0] <- 0 
  lcorr[abs(lcorr-1) < 10^(-10)] <- 1 
  lcorr[abs(lcorr-1) > 10^(-10)] <- 0 
  ccl <- apply(lcorr,1,sum) #Default for y weight. 
  
  ny <- 0 
  
  if (max(ccl)==1){ #y is in the set of observations. 
    cpp <-ccl[1:(length(ccl)-1)] 
    ny <- Nvec[cpp==1] 
  } 
  
  dvec <- Distvec(Xmatr,y) #Distance from Xmatr elements to y. 
  
  dind <- seq(1,length(dvec)) 
  
  dinds <- dind[abs(dvec) > 10^(-10)]  
  
  wvec <- matrix(0,1,dim(Xmatr)[1]) 
  
  wvec[dinds] <- (Nvec[dinds] / dvec[dinds]) 
  
  wvec <- wvec / sum(wvec) 
  
  Ttilde <- apply(apply(Xmatr,2,rowweight,a=wvec),2,sum) 
  
  unitvecs <- t(apply(Xmatr,1,unitvec,b=y)) 
  
  wunitvecs <- apply(unitvecs,2,rowweight,a=Nvec) 
  
  ry <- sqrt(sum(apply(wunitvecs,2,sum)^2)) 
  
  ynew <- Ttilde * max(0,(1-ny/ry)) + y * min(1,ny/ry) 
  
  unitvecs <- t(apply(Xmatr,1,unitvec,b=ynew)) 
  
  ebar <- sqrt(sum(apply(apply(unitvecs,2,rowweight,a=Nvec/sum(Nvec)),2,sum)^2)) 
  
  l <- rbind(ynew,Xmatr) 
  lcorr <- cor(t(l)) 
  diag(lcorr) <- 0 
  uc<-upper.tri(lcorr) 
  lcorr[uc==0] <- 0 
  lcorr[abs(lcorr-1) < 10^(-10)] <- 1 
  lcorr[abs(lcorr-1) > 10^(-10)] <- 0 
  ccl <- apply(lcorr,1,sum) 
  
  ny<-0 
  
  if(max(ccl)==1){ 
    ny<-Nvec[ccl==1] 
  } 
  if(ny==0){ 
    DDy <- 1 - ebar
  } 
  if(ny>0) { 
    DDy <- 1 - max(0,(ebar-ny/sum(Nvec))) 
    #DDy <- 1 - ebar
  } 
  y<-ynew 
  
  return(list(y=y,DDy=DDy))
} 


##TO RUN THE PROGRAM## 
#yr <- initW2(Xmatr,Nvecr) 
#W <- Weisziteradj(Xmatr,Nvecr,yr,B) #Usually enough to run B=5, can check by DDy to see that it is equal to 1. 
#Y <- W$Y #This is the multivariate median. 
#Computing depth of an observations with respect to a cluster. 

DDapply <- function(Xv,Xmat) { #Computes the depths of a vector of observations Xv with respect to
#a cloud of observations Xmat. 
  Nvec <- rep(1,dim(Xmat)[1]) 
  dd <- DDfcnadj(Xmat,Nvec,Xv) 
  return(dd) 
} 

DDfcnadj <- function(Xmatr,Nvec,y){ #Computes the data depth of point y in cloud Xmatr. 
 
  aux <- apply(Xmatr,1,unitvec,b=y)
  if(class(aux) == "list"){
    unitvecs <- t(matrix(unlist(aux),dim(Xmatr)[1],dim(Xmatr)[2])) 
  }else{
    unitvecs <- t(aux) 
  }
  
  #unitvec<-function(a,b){ #Unit vector from a to b. 
  #l1<-sqrt(sum((a-b)^2)) 
  #if (l1 > 10^(-10)){ #10^(-10) 
  # uv <- (b-a)/l1 
  #}
  # if (l1 < 10^(-10)){ 
  #  uv <- matrix(0,1,length(a)) 
  # } 
  # return(uv) 
  #} 
  
  gvv = matrix(as.numeric(unitvecs), dim(unitvecs)[1], dim(unitvecs)[2])
 
  a <- Nvec/sum(Nvec)  
  ebar <- sqrt( sum ( apply( rowweight(gvv,a) ,2,sum )^2))   
  
  #rowweight<-function(X,a){ #Product function 
  # xm<-X*a 
  # return(xm) 
  #}
  
  l <- rbind(y,Xmatr) 
  
  lcorr <- cor(t(l)) 
  
  diag(lcorr) <- 0 
  
  uc <- upper.tri(lcorr) 
  
  lcorr[uc==0] <- 0 
  lcorr[abs(lcorr-1) < 10^(-10)] <- 1 
  lcorr[abs(lcorr-1) > 10^(-10)] <- 0 
  ccl <- apply(lcorr,1,sum) 
  
  ny <- 0 
  
  if(max(ccl)==1){ 
    ny <- sum(Nvec[ccl==1])
  } 
  
  if(ny == 0){ 
    DDy <- 1-ebar
  }
  
  if(ny > 0){ 
    DDy <- 1-max(0,(ebar-ny/sum(Nvec))) 
    #DDy <- 1-ebar 
  } 
  return(DDy) 
} 

DDapplyloo <- function(Xv,Xmat){ 
  
  dl <- apply(Xmat,1,l1norm,b=Xv) 
  
  Xmatu <- Xmat[dl>0,] 
  
  Nvec <- rep(1,dim(Xmatu)[1]) 
  
  dd <- DDfcnadj(Xmatu,Nvec,Xv) 
  
  return(dd) 
} 

##TO RUN## 
#x is a vector, X is a group of data,N is the vector of multiplicities.
#Dx <- apply(x,1,DDapply,Xmat=X,Nvec=N) #Depth of a single observation x. 
#Dx <- DDfcnadj(Xmat,Nvec,x) 


DDcalc2 <- function(DDi,NN,K,Km,norm=1){ 
  
  n <- NN[1,] 
  nvec <- matrix(0,K,1) 
  ovec <- matrix(0,length(n),1) 
  svec <- matrix(0,length(n),1) 
  indvec <- seq(1,length(n)) 
  normvec <- matrix(0,K,1) 
  
  for(kz in (1:K)){ 
    nvec[kz] <- length(n[n==kz]) 
    Vin <- DDi[n==kz,kz] 
    normvec[kz] <- mean(Vin) 
  } 
  if(norm==0){ 
    normvec <- normvec/normvec
  } 
  DDb <- matrix(0,length(n),1) 
  DDw <- DDb 
  DD <- DDb 
  Vt <- 0 
  VZD <- DDi 
  
  for(ml in (1:length(n))){ 
    if(ml==1){ 
      Vt <- DDi[ml,c(NN[2:Km,ml])] / normvec[c(NN[2:Km,ml])]
      
      VZD[ml,] <- DDi[ml,] / normvec[c(NN[,ml])] 
    }
    if(ml>1){ 
      Vt <- c(Vt,DDi[ml,c(NN[2:Km,ml])] / normvec[c(NN[2:Km,ml])]) 
      
      VZD[ml,] <- DDi[ml,] / normvec[c(NN[,ml])] 
    }
  } 
  TT <- Vt[rev(sort.list(Vt))[length(NN[1,])+1]] #VZD2<-VZD #VZD[VZD<-0 
  
  for(kt in (1:length(n))){ 
    Nv <- NN[2:Km,kt]
    
    Vv <- sum(VZD[kt,Nv]) 
    DDw[kt] <- DDi[kt,NN[1,kt]] / normvec[NN[1,kt]] 
    DDb[kt] <- VZD[kt,NN[2,kt]] / normvec[c(NN[2,kt])] 
    ti <- DDi[kt,NN[1,kt]] / normvec[NN[1,kt]]
    DD[kt] <- (ti-Vv) 
  } 
  
  return(list(DDw=DDw,DD=DD))
} 


pamsil <- function(X,pc,K) { #Silwidth calculation. 
  
  sil <- matrix(0,length(pc),1) 
  
  Kvec <- seq(min(pc):max(pc)) 
  
  Nuvec <- pc 
  
  ccl <- rep(0,length(pc)) 
  
  clvec <- matrix(0,length(pc),K) 
  
  clvec[,1] <- pc 
  
  clvec2 <- clvec 
  
  indvec <- seq(1:length(Nuvec)) 
  
  for(kt in (1:length(pc))){ 
    kk <- Nuvec[kt] 
    
    kc <- setdiff(indvec,kt) 
    
    Nt <- Nuvec[kc] 
    
    ai <- mean(Distvec(X[kc[Nt==kk],],X[kt,])) 
    
    Ks <- setdiff(Kvec,kk) 
    
    for(kz in (1:length(Ks))){ 
      if(kz==1){ 
        bi <- mean(Distvec(X[Nuvec==Ks[kz],],X[kt,])) 
        
        ccl[kt] <- Ks[kz] 
        
        clvec2[kt,2] <- bi 
        
        clvec[kt,2] <- Ks[kz] 
      } 
      if(kz>1){ 
        nnd <- mean(Distvec(X[Nuvec==Ks[kz],],X[kt,])) 
        if(bi>nnd){ 
          ccl[kt] <- Ks[kz] 
        } 
        bi <- min(bi,nnd) 
        clvec2[kt,kz+1] <- nnd 
        clvec[kt,kz+1] <- Ks[kz] 
      } 
    } 
    
    sili <- (bi-ai) / max(bi,ai) 
    
    sil[kt] <- sili 
  } 
  clvec <- clvec[,2:K] 
  clvec2 <- clvec2[,2:K] 
  
  if(K>2){ 
    Sm <- apply(clvec2,1,sort.list) 
    clvec3 <- clvec 
    
    for(kk in (1:length(Nuvec))){ 
      clvec3[kk,] <- clvec[kk,Sm[,kk]] 
    } 
    clvec <- clvec3 
  } 

  return(list(sil=sil,ccl=ccl,clvec=clvec)) 
} 


NNDDVQA1 <- function(X,pcc,lambda,norm){

  #X is the data matrix. 
  #pcc particion actual.
  #MM is the matrix of multivariate medians. 
  #T is current temperature. 
  
  K <- dim(pcc$med)[1] 
  Km <- 2 
  Kmat <- matrix(0,dim(X)[1],K) 
  
  for(kk in (1:K)){ 
    Kmat[,kk]<-apply(X,1,l1norm,b=pcc$med[kk,]) 
  } 
  
  Smat <- apply(Kmat,1,sort.list) 
  
  Kmat <- apply(Kmat,1,sort) 
  NNuse <- Smat           
  
  Nvec <- rep(1,dim(X)[1]) 
  DDi <- matrix(0,dim(X)[1],K) 
  
  for(kz in (1:K)){ 
    for(ky in (1:dim(X)[1])){ 
      Xmatr <- X[NNuse[1,] == kz,]
      
      Nvecr <- Nvec[NNuse[1,] == kz] 
      
      DDi[ky,kz] <- DDfcnadj(Xmatr,Nvecr,X[ky,]) 
    } 
  } 
  
  DD <- DDcalc2(DDi,NNuse,K,Km,norm) 
    
  pcc <- pamsil(X,NNuse[1,],K)   
  
  Kmata <- pcc$sil * (1 - lambda) + lambda * DD$DD 
  
  NN <- NNuse 
  NN[2:K,] <- t(pcc$clvec) 
  
  Nuvec <- NNuse[1,] 
  
  Kmatb <- 0 
  prs <- 0 
  stopnow <- 0 
  
  return(list(NN=NN,DDi=DDi,DD=DD,pcc=pcc,Kmata=Kmata,Nuvec=Nuvec,Kmatb=Kmatb,prs=prs,stopnow=stopnow)) 
} 


NNDDVQE <- function(X,MM,prT,lambda,NNold,Km,Th,Trimm,kl,NNold_old,verbose){ 
  
  #X is the data matrix. 
  #MM is the matrix of multivariate medians. 
  #Th is current temperature. 
  
  change <- c()
  iter_while <- 0
  improv <- c()
  norm <- 0 
  K <- dim(MM)[1] 
  
  if(K <= 3){ 
    Km <- 2 
  } 
    
  aux <- sort(NNold$Kmata)
  
  tri <- Trimm*(dim(X)[1])
  if(tri < 1){
   tri1 <- ceiling(tri)  
  }else{
    tri1 <- round(tri)
  }
  
  phi <- aux[abs(tri1)]
  trimm <- which(NNold$Kmata <= phi)
  if(verbose){
   cat("Trimmed women: ")
   print(trimm)
   cat("\n")
  } 
  X <- X[-trimm,]
  
  sil <- as.matrix(NNold$pcc$sil[-trimm,])
  ccl <- NNold$pcc$ccl[-trimm]
  clvec <- as.matrix(NNold$pcc$clvec[-trimm,])
  pcc_aux <- list(sil=sil,ccl=ccl,clvec=clvec)  
  
  pcc <- pcc_aux
  NNuse <- NNold$NN[,-trimm]
  NNuse[2:K,] <- t(pcc$clvec) 
  
  Nvec <- rep(1,dim(X)[1])
  DDi <- NNold$DDi[-trimm,] 
  DDw <- as.matrix(NNold$DD$DDw[-trimm,])
  DD <- as.matrix(NNold$DD$DD[-trimm,])
  DD <- list(DDw=DDw,DD=DD)
  
  Kmata <- pcc$sil * (1 - lambda) + lambda * DD$DD 
  
  if((kl-1) == 0){
    if(verbose){
     cat("Initial value of the partition:")
     print(mean(Kmata))
     cat("\n")
    } 
  }else{
    if(verbose){
     cat("Initial value of the partition:")
     print(NNold$Cost)
     cat("\n")
    } 
  } 
  
  Th <- as.numeric(quantile(Kmata,0.25))
  
  Kmat0 <- Kmata 
  Kmatb <- Kmata
  Kmatb[Kmata > Th] <- 1 
  Kmatb[Kmata <= Th] <- 0 
  
  NNnew <- NNuse 
  NNusea <- NNuse 
  indvec <- seq(1,dim(X)[1]) 
  
  Iuse <- indvec[Kmatb==0] #Elements that can be moved. 
  Iuse <- setdiff(Iuse,trimm)
  if(verbose){
   cat("Set S of observations that can be relocated:")
   print(Iuse)
   cat("\n")
  }
  ct <- 0 
  nostop <- 0 
  nb <- min(50,length(Iuse)) 
  
  if(verbose){
   cat("While loop starts: \n")
  } 
  while(length(Iuse)>0 & nostop<=nb){ 
    iter_while <- iter_while + 1
    if(length(Iuse)==1){ 
      iuse <- Iuse
      if(verbose){
       cat("Random subset E of S:")
       print(iuse)
      } 
    } 
    if(length(Iuse)>1){ 
      E <- round(runif(1,min=1,max=min(10,length(Iuse)))) 
      iuse <- sample(Iuse,E) 
      if(verbose){
       cat("Random subset E of S:")
       print(iuse)
      } 
    } 
    for(zz in (1:length(iuse))){ 
      NNusea[1,iuse[zz]] <- NNuse[2,iuse[zz]] 
      NNusea[2,iuse[zz]] <- NNuse[1,iuse[zz]] 
    } 
    DDia <- DDi-DDi 
    
    for(kz in (1:K)){ 
      Xmatr <- X[NNusea[1,]==kz,] 
      Nvecr <- Nvec[NNusea[1,]==kz] 
      DDia[,kz] <- apply(X,1,DDapply,Xmat=Xmatr) 
    } 
    
    DDa <- DDcalc2(DDia,NNusea,K,Km,norm) 
    pcca <- pamsil(X,NNusea[1,],K) 
    if(verbose){
     cat("Candidate partition:")
     print(table(NNusea[1,]))
     cat("\n")
    } 
    Kmataa <- pcca$sil * (1 - lambda) + lambda * DDa$DD #Value of new partition. 
    
    if((kl-1) == 0){
      Del <- mean(Kmataa) - mean(Kmata) 
    }else{
      Del <- mean(Kmataa) - NNold$Cost
    }
    
    if(Del == "NaN"){
     stop("This partition cannot be evaluated \n")  
     break
    }
    
    if((kl-1) == 0){
      if(verbose){
       cat("Value of the previous partition C(I_1^{K}):")
       print(mean(Kmata))
      } 
    }else{
      if(verbose){
       cat("Value of the previous partition C(I_1^{K}):")
       print(NNold$Cost)
      } 
    } 
    
    if(verbose){
     cat("Value of the new partition C(tilde{I}_1^{K}):")
     print(mean(Kmataa))
     cat("Current Delta value:")
     print(Del)
     cat("\n")
    }
    #Stopping criterion:
    if(Del > 0){
      if((Del / mean(Kmataa)) < 0.01){
        if(verbose){
         cat("STOPPING CRITERION FULFILLED")
         cat("\n")
        }
         change[iter_while] <- 0
        break         
      }
    } 
    
    pru <- runif(1)
    
    if(prT>0){
      PRT <- exp(Del/prT) 
      if(verbose){
       cat("Probability of mis-allocation:")
       print(PRT)
       cat("\n")
      }  
    } 
    if(prT==0){ 
      PRT <- 0 
      if(verbose){
       cat("Probability of mis-allocation:")
       print(PRT)
       cat("\n")
      } 
    } 

    if(Del>0 | (Del<=0 & pru==PRT)){ 
      NNuse <- NNusea
      Kmata <- Kmataa
      if(verbose){
       cat("The partition is accepted.\n")
       cat("New value of the partition:")
       print(mean(Kmata))
       cat("Partition:\n")
       print(table(NNuse[1,]))
      }
       NNold$Cost <- mean(Kmata)
      if(verbose){
       cat("\n")
      }
       improv[kl] <- kl - 1 
       change[iter_while] <- 1
    }else{
      NNold$Cost <- NNold$Cost  
      if(verbose){
       cat("The partition is rejected.\n")
       cat("Current value of the partition:")
       print(NNold$Cost)
       cat("\n")
      } 
      change[iter_while] <- 0
    } 
    nostop <- nostop+1 
  }
  if(verbose){
   cat("While loop ends: \n")
  }
  if(max(change) == 1){
    if(ct==0){ 
      stopnow <- 1 
    } 
    if(ct>0){ 
      stopnow <- 0 
    } 
    prs <- ct/length(NNuse[1,]) 
    NNnew <- NNuse 
    NN <- NNnew 
    Nuvec <- NNnew[1,] 
    
    DDia <- matrix(0,nrow=dim(X)[1],ncol=K)
    for(kz in (1:K)){ 
      Xmatr <- X[NNuse[1,]==kz,] 
      Nvecr <- Nvec[NNuse[1,]==kz] 
      DDia[,kz]<-apply(X,1,DDapply,Xmat=Xmatr) 
    } 
    
    DDa <- DDcalc2(DDia,NNuse,K,Km,norm) 
    pcc <- pamsil(X,NNuse[1,],K) 
    
    return(list(Iuse=Iuse,ct=ct,NN=NN,DDi=DDia,DD=DDa,pcc=pcc,Kmata=Kmata,Kmat0=Kmat0,Nuvec=Nuvec,Kmatb=Kmatb,
                prs=prs,stopnow=stopnow,Cost=NNold$Cost,trimmed=trimm,improv=improv))
  }else if(max(change) == 0 | is.null(change)){
    NNold_old <- NNold_old
    if((kl-1) == 0){
      NNold_old$Cost <- mean(Kmata)
    }else{
      NNold_old$Cost <- NNold$Cost
     } 
    NNold_old$ct <- ct
    NNold_old$trimmed <- trimm
    NNold_old$improv <- improv
    return(NNold_old)
   } 
} 


NNDDVQEstart <- function(X,MM,prT,lambda,NNold,Km,Th,Trimm,kl,NNold_old,verbose){ 
  #X is the data matrix. 
  #MM is the matrix of multivariate medians. 
  #Th is current temperature. 
   
  change <- c()
  iter_while <- 0
  improv <- c()
  norm <- 0 
  K <- dim(MM)[1] 
  
  if(K < 3){ 
    Km <- 2 
  } 
  
  
  aux <- sort(NNold$Kmata)
  
  tri <- Trimm*(dim(X)[1])
  if(tri < 1){
    tri1 <- ceiling(tri)  
  }else{
    tri1 <- round(tri)
  }
  
  phi <- aux[abs(tri1)]
  trimm <- which(NNold$Kmata <= phi)
  if(verbose){
   cat("Trimmed women: ")
   print(trimm)
   cat("\n")
  }
   X <- X[-trimm,]
  
  sil <- as.matrix(NNold$pcc$sil[-trimm,])
  ccl <- NNold$pcc$ccl[-trimm]
  clvec <- as.matrix(NNold$pcc$clvec[-trimm,])
  pcc_aux <- list(sil=sil,ccl=ccl,clvec=clvec)  
  
  pcc <- pcc_aux
  NNuse <- NNold$NN[,-trimm]
  NNuse[2:K,] <- t(pcc$clvec) 

  Nvec <- rep(1,dim(X)[1])
  DDi <- NNold$DDi[-trimm,] 
  DDw <- as.matrix(NNold$DD$DDw[-trimm,])
  DD <- as.matrix(NNold$DD$DD[-trimm,])
  DD <- list(DDw=DDw,DD=DD)
  
  Kmata <- pcc$sil * (1 - lambda) + lambda * DD$DD 
  
  if((kl-1) == 0){
    if(verbose){
     cat("Initial value of the partition:")
     print(mean(Kmata))
     cat("\n")
    } 
  }else{
    if(verbose){
     cat("Initial value of the partition:")
     print(NNold$Cost)
     cat("\n")
    }
  } 
  
  Th <- as.numeric(quantile(Kmata,0.25))
  
  Kmat0 <- Kmata 
  Kmatb <- Kmata
  Kmatb[Kmata > Th] <- 1 
  Kmatb[Kmata <= Th] <- 0 
  
  NNnew <- NNuse 
  NNusea <- NNuse 
  indvec <- seq(1,dim(X)[1]) 
  
  Iuse <- indvec[Kmatb==0] #Elements that can be moved. 
  Iuse <- setdiff(Iuse,trimm)
  if(verbose){
   cat("Set S of observations that can be relocated:")
   print(Iuse)
   cat("\n")
  }
  ct <- 0 
  fir <- 0 
  nostop <- 0 
  nb <- min(50,length(Iuse)) 
  
  if(verbose){
   cat("While loop starts:\n")
  }
  while(length(Iuse)>0 & nostop<=nb){
    iter_while <- iter_while + 1
    if(length(Iuse)==1){ 
      iuse <- Iuse 
      if(verbose){
       cat("Random subset E of S:")
       print(iuse)
      } 
    } 
    if(length(Iuse)>1){ 
      E <- round(runif(1,min=1,max=min(5,length(Iuse)))) 
      iuse <- sample(Iuse,E) 
      if(verbose){
       cat("Random subset E of S:")
       print(iuse)
      }  
     } 
    for(zz in (1:length(iuse))){ 
      NNusea[1,iuse[zz]] <- NNuse[2,iuse[zz]]
      NNusea[2,iuse[zz]] <- NNuse[1,iuse[zz]] 
    } 
    DDia <- DDi-DDi 
    
    for(kz in (1:K)){ 
      Xmatr <- X[NNusea[1,]==kz,] 
      Nvecr <- Nvec[NNusea[1,]==kz] 
      DDia[,kz]<-apply(X,1,DDapply,Xmat=Xmatr) 
    } 
    
    DDa <- DDcalc2(DDia,NNusea,K,Km,norm) 
    pcca <- pamsil(X,NNusea[1,],K) 
    if(verbose){
     cat("Candidate partition:")
     print(table(NNusea[1,]))
     cat("\n")
    }
     Kmataa <- pcca$sil * (1 - lambda) + lambda * DDa$DD 
    
    if((kl-1) == 0){
      Del <- mean(Kmataa) - mean(Kmata)
    }else{
      Del <- mean(Kmataa) - NNold$Cost
    }
    
    if(Del == "NaN"){
      stop("This partition cannot be evaluated \n")  
      break
    }
    
    if((kl-1) == 0){
      if(verbose){
       cat("Value of the previous partition C(I_1^{K}):")
       print(mean(Kmata))
      }  
    }else{
      if(verbose){
       cat("Value of the previous partition C(I_1^{K}):")
       print(NNold$Cost)
      } 
    } 
    
    if(verbose){
     cat("Value of the new partition C(tilde{I}_1^{K}):")
     print(mean(Kmataa))
     cat("Current Delta value:")
     print(Del)
     cat("\n")
    }
    #Stopping criterion:
    if(Del > 0){
      if((Del / mean(Kmataa)) < 0.01){
        if(verbose){
         cat("STOPPING CRITERION FULFILLED")
         cat("\n")
        }
        change[iter_while] <- 0
        break         
      }
    }
    
    pru <- runif(1) 
    if(fir==0 & prT>0){ 
      prT <- -abs(Th)/log(prT) 
      T0 <- prT 
      PRT <- exp(Del/prT) 
      fir <- 1 
      if(verbose){
       cat("Probability of mis-allocation:")
       print(PRT)
       cat("\n")
      }    
    } 
    if(fir!=0 & prT>0){ 
      PRT <- exp(Del/prT) 
      if(verbose){
       cat("Probability of mis-allocation:")
       print(PRT)
       cat("\n")
      } 
    } 
    if(prT==0){ 
      PRT <- 0 
      if(verbose){
       cat("Probability of mis-allocation:")
       print(PRT)
       cat("\n")
      } 
    } 
    if(fir==0){ 
      PRT <- -abs(Del)/log(prT) 
      if(verbose){
       cat("Probability of mis-allocation:")
       print(PRT)
       cat("\n")
      }
      T0 <- prT 
      fir <- 1 
    } 
    
    if(verbose){
     cat("Probability of mis-allocation:")
     print(PRT)
     cat("\n")
    }
    if(Del>0 | (Del<=0 & pru==PRT)){ 
      NNuse <- NNusea
      Kmata <- Kmataa
      if(verbose){
       cat("The partition is accepted.\n")
       cat("New value of the partition:")
       print(mean(Kmata))
       cat("Partition:\n")
       print(table(NNuse[1,]))
      }
       NNold$Cost <- mean(Kmata)
      if(verbose){
       cat("\n")
      }
       improv[kl] <- kl - 1 
      change[iter_while] <- 1
    }else{
      NNold$Cost <- NNold$Cost  
      if(verbose){
       cat("The partition is rejected.\n")
       cat("Current value of the partition:")
       print(NNold$Cost)
       cat("\n")
      }
       change[iter_while] <- 0
    } 
    nostop <- nostop+1 
  }
  if(verbose){
   cat("While loop ends: \n")  
  }
  if(max(change) == 1){
    if(ct==0){ 
      stopnow <- 1 
    } 
    if(ct>0){ 
      stopnow <- 0 
    } 
    prs <- ct/length(NNuse[1,]) 
    NNnew <- NNuse 
    NN <- NNnew 
    Nuvec <- NNnew[1,] 
    
    DDia <- matrix(0,nrow=dim(X)[1],ncol=K)
    for(kz in (1:K)){ 
      Xmatr <- X[NNuse[1,]==kz,] 
      Nvecr <- Nvec[NNuse[1,]==kz] 
      DDia[,kz]<-apply(X,1,DDapply,Xmat=Xmatr) 
    } 
    
    DDa <- DDcalc2(DDia,NNuse,K,Km,norm) 
    pcc <- pamsil(X,NNuse[1,],K) 
    
    return(list(Iuse=Iuse,ct=ct,NN=NN,DDi=DDia,DD=DDa,pcc=pcc,Kmata=Kmata,Kmat0=Kmat0,Nuvec=Nuvec,Kmatb=Kmatb,
                prs=prs,stopnow=stopnow,Cost=NNold$Cost,trimmed=trimm,improv=improv,T0=T0))
  }else if(max(change) == 0 | is.null(change)){
    NNold_old <- NNold_old
    if((kl-1) == 0){
      NNold_old$Cost <- mean(Kmata)
    }else{
      NNold_old$Cost <- NNold$Cost
    } 
    NNold_old$ct <- ct
    NNold_old$trimmed <- trimm
    NNold_old$improv <- improv
    NNold_old$T0 <- T0
    return(NNold_old)
  }
}
