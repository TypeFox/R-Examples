#########################################################################################################
##########
##########              FONCTIONS C++ CALL IN R
##########
#########################################################################################################


ClassiSeg <- function(geno, grille, Kmax){
	nRow <- length(geno)
	nGrille <- length(grille)
      	### Appel C
	A <- .C("ClassiSeg",as.double(geno), as.integer(nRow), as.integer(Kmax),  res1=double(Kmax*nRow), res2=integer(Kmax*nRow), as.integer(nGrille), moyennes=as.double(grille), PACKAGE="cghseg")
	A$res1 <- matrix(A$res1, nrow=Kmax, byrow=TRUE)
	A$res2 <- matrix(A$res2, nrow=Kmax, byrow=TRUE)
	n <- ncol(A$res1)
	res3 <- matrix(NA, nrow=nrow(A$res2), ncol=nrow(A$res2))
	res3[1, 1] <- 0
	for(i in 2: nrow(A$res2)){
		res3[i, i-1] <- A$res2[i, n]
		for(k in 1:(i-1)){
		res3[i, i-1-k] <- A$res2[i-k, res3[i, i-k]]
		}		
	}
	diag(res3) <- ncol(A$res1)
	return(list(J.est=A$res1[,Kmax],t.est=res3))
	} 
### logP = donner -logProba du melange
### variance = variance du modèle
ClassiSeg_2 <- function(geno, grille, Kmax, logP, variance){
	nRow <- length(geno)
	nGrille <- length(grille)
### Appel C
	A <- .C("ClassiSegProba",as.double(geno), as.integer(nRow), as.integer(Kmax),  res1=double(Kmax*nRow), 
			res2=integer(Kmax*nRow), as.integer(nGrille), moyennes=as.double(grille), 
			logP=as.double(logP), variance=as.double(variance), PACKAGE="cghseg")
	A$res1 <- matrix(A$res1, nrow=Kmax, byrow=TRUE)
	A$res2 <- matrix(A$res2, nrow=Kmax, byrow=TRUE)
	n <- ncol(A$res1)
	res3 <- matrix(NA, nrow=nrow(A$res2), ncol=nrow(A$res2))
	res3[1, 1] <- 0
	for(i in 2: nrow(A$res2)){
		res3[i, i-1] <- A$res2[i, n]
		for(k in 1:(i-1)){
			res3[i, i-1-k] <- A$res2[i-k, res3[i, i-k]]
		}		
	}
	diag(res3) <- ncol(A$res1)
	return(list(J.est=A$res1[,Kmax],t.est=res3))
} 

ClassiSeg_ <- function(geno, grille, Kmax){
	nRow <- length(geno)
	nGrille <- length(grille)
      	### Appel C
	A <- .C("ClassiSeg",as.double(geno), as.integer(nRow), as.integer(Kmax),  res1=double(Kmax*nRow), res2=integer(Kmax*nRow), as.integer(nGrille), moyennes=as.double(grille), PACKAGE="cghseg")
	A$res1 <- matrix(A$res1, nrow=Kmax, byrow=TRUE)
	A$res2 <- matrix(A$res2, nrow=Kmax, byrow=TRUE)
	return(A)
	} 

ClassiSeg_2_ <- function(geno, grille, Kmax, logP, variance){
	nRow <- length(geno)
	nGrille <- length(grille)
### Appel C
	A <- .C("ClassiSegProba",as.double(geno), as.integer(nRow), as.integer(Kmax),  res1=double(Kmax*nRow), 
			res2=integer(Kmax*nRow), as.integer(nGrille), moyennes=as.double(grille), 
			logP=as.double(logP), variance=as.double(variance), PACKAGE="cghseg")
	A$res1 <- matrix(A$res1, nrow=Kmax, byrow=TRUE)
	A$res2 <- matrix(A$res2, nrow=Kmax, byrow=TRUE)
	return(A)
} 
colibriR_c <- function(signalBruite, Kmax, mini=min(signalBruite), maxi=max(signalBruite)){
	n <- length(signalBruite)
    A <- .C("colibriR_c", signal=as.double(signalBruite), n=as.integer(n), Kmax=as.integer(Kmax),   min=as.double(mini), max=as.double(maxi), path=integer(Kmax*n), cost=double(Kmax)
	, PACKAGE="cghseg")
    A$path <- matrix(A$path, nrow=Kmax, byrow=T)
    A$cost <- A$cost + sum(signalBruite^2)
    return(A);	
} 

meanRuptR_c <- function(Ym, rupt, k){
	A <- .C("meanRuptR_c", data=as.double(Ym), position=as.integer(rupt), k=as.integer(k), res=double(k), PACKAGE="cghseg")	
	return(A$res)
}

meansqRuptR_c <- function(Ym, rupt, k){
	A <- .C("meansqRuptR_c", data=as.double(Ym), position=as.integer(rupt), k=as.integer(k), res=double(k), PACKAGE="cghseg")	
	return(A$res)
}


retour <- function(path, i){
   chaine <- integer(i)
   chaine[i] <- ncol(path)
   for(j in i:2)  chaine[j-1] <- path[j, chaine[j]]
   return(chaine)
}

segmeanCO <- function(geno, Kmax){
  out        = colibriR_c(geno, Kmax)
  t.est      = t(sapply(1:Kmax,FUN=function(k){c(retour(out$path, k),rep(0,Kmax-k))}))
  t.est[1,1] = length(geno)
  res        = list(J.est  = out$cost,t.est  = t.est)
  return(res)
}


Estep <- function(logdensity,phi) {
	storage.mode(logdensity)<-"double";
	storage.mode(phi)<-"double";
	.Call("sc_estep",logdensity,as.integer(dim(logdensity)[1]),as.integer(dim(logdensity)[2]),phi);
}

hybrid <- function(x,P,Kmax,lmin=1,lmax=length(x),vh=TRUE,fast) {
  
  checkoptions = TRUE
  if ((vh==FALSE) & (lmin<=1)){
    checkoptions = FALSE
    cat("Error in hybrid : lmin must be greater than 2 when vh=FALSE","\n")     
  }
  
  if (Kmax>floor(length(x)/lmin)){
    checkoptions = FALSE
    cat("Error in hybrid : Kmax must be lower than [length(x)/lmin]","\n")
  }
  
  if (Kmax<floor(length(x)/lmax)){
    checkoptions = FALSE
    cat("Error in hybrid : Kmax must be greater than [length(x)/lmax]","\n")
  }
  
  if (P>Kmax){
    checkoptions = FALSE
    cat("Error in hybrid : the number of groups must be lower than the number of segments")
  }
   if (sum(is.na(x))>0){
     checkoptions = FALSE
     cat("Error in hybrid : the data must not contain missing value","\n")
   }
  
  if (checkoptions==TRUE){
    storage.mode(x)<-"double"
    res         = .Call("sc_hybrid",x,as.integer(P),as.integer(Kmax),as.integer(lmin),as.integer(lmax),as.logical(vh),as.logical(fast))   
    tmp         = res$Linc
    if (sum(tmp==0)>0){
      cat("The algorithm has faced convergence problems for P =", P,"\n")
      cat("and for configurations with", which(tmp==0)," segments","\n")
      tmp[tmp==0] = -Inf
      res$Linc    = tmp
    }
  }
  invisible(res)	
}

logdens <- function(x,P,phi) {
  storage.mode(x)<-"double";
  storage.mode(phi)<-"double";
  .Call("sc_logdens",x,as.integer(P),phi);
}

segmeanGR <- function(geno, Kmax){
  nRow       = length(geno)
  A          = .C("LinProgDyn",as.double(geno), as.integer(nRow), as.integer(Kmax),res1=double(Kmax*nRow), res2=integer(Kmax*nRow))
  A$res1     = matrix(A$res1, nrow=Kmax, byrow=TRUE)
  A$res2     = matrix(A$res2, nrow=Kmax, byrow=TRUE)
  n          = ncol(A$res1)
  res3       = matrix(NA, nrow=nrow(A$res2), ncol=nrow(A$res2))
  res3[1, 1] = 0  
  for(i in c(2: nrow(A$res2))){
    res3[i,(i-1)] = A$res2[i, n]
    for(k in c(1:(i-1))){
      res3[i, (i-1-k)] =  A$res2[i-k, res3[i, (i-k)]]
    }    
  }
  diag(res3) = ncol(A$res1)  
  res        = list()
  res$J.est  = A$res1[, ncol(A$res1)]
  res$t.est  = res3
  return(res)
}


segmixtGR <- function(geno,P,Kmax,phi){
  moyennes   = phi[1:P]
  variances  = phi[(P+1):(2*P)]^2
  prop       = phi[(2*P+1):(3*P)]  
  nRow       = length(geno)
  A          = .C("LinProgDynMelange",as.double(geno), as.integer(nRow), as.integer(Kmax),
    res1=double(Kmax*nRow), res2=integer(Kmax*nRow), as.integer(length(prop)), as.double(moyennes), 
    as.double(variances), as.double(prop))
  A$res1     = matrix(A$res1, nrow=Kmax, byrow=TRUE)
  A$res2     = matrix(A$res2, nrow=Kmax, byrow=TRUE)
  n          = ncol(A$res1)
  res3       = matrix(NA, nrow=nrow(A$res2), ncol=nrow(A$res2))
  res3[1, 1] = 0  
  for(i in c(2: nrow(A$res2))){
    res3[i,(i-1)] = A$res2[i, n]
    for(k in c(1:(i-1))){
      res3[i, (i-1-k)] =  A$res2[i-k, res3[i, (i-k)]]
    }    
  }
  diag(res3) = ncol(A$res1)  
  res        = list()
  res$J.est  = A$res1[, ncol(A$res1)]
  res$t.est  = res3
  return(res)
}



segmean <- function(x,Kmax,lmin=1,lmax=length(x),vh=TRUE) {
  
  checkoptions = TRUE
  if ((vh==FALSE) & (lmin<=1)){
    checkoptions = FALSE
    cat("Error in segmean : lmin must be greater than 2 when vh=FALSE","\n")     
  }
  if (Kmax> floor(length(x)/lmin)){
    checkoptions = FALSE
    cat("Error in segmean : Kmax must be lower than [length(x)/lmin]","\n")
  }
  if (Kmax<floor(length(x)/lmax)){
    checkoptions = FALSE
    cat("Error in segmean : Kmax must be greater than [length(x)/lmax]","\n")
  }
  
  if (sum(is.na(x))>0){
    checkoptions = FALSE
    cat("Error in segmean : the data must not contain missing values","\n")
  }
  
  if (checkoptions == TRUE){
    storage.mode(x)<-"double"
    .Call("sc_segmean",x,as.integer(lmin),as.integer(lmax),as.integer(Kmax),as.logical(vh))
  }
}

segmixt <- function(x,P,Kmax,phi,lmin=1,lmax=length(x)) {
  
  checkoptions = TRUE
  if (Kmax> floor(length(x)/lmin)){
    cat("Error in segmixt : Kmax must be lower than [length(x)/lmin]","\n")
    checkoptions = FALSE
  }
  
  if (Kmax<floor(length(x)/lmax)){
    checkoptions = FALSE
    cat("Error in segmixt : Kmax must be greater than [length(x)/lmax]","\n")
  }
  
  if (P>Kmax){
    checkoptions = FALSE
    cat("Error in segmixt : the number of groups must be lower than the number of segments","\n")
  }
  if (sum(is.na(x))>0){
    checkoptions = FALSE
    cat("Error in segmixt : the data must not contain missing values","\n")
  }
  if (checkoptions == TRUE){
    storage.mode(x)<-"double"
    storage.mode(phi)<-"double"
    .Call("sc_segmixt",x,as.integer(lmin),as.integer(lmax),as.integer(Kmax),phi, as.integer(P))
  }
}


segibp <- function(Contrast,Kseq,multiKmax){
  Contrast  = as.vector(Contrast$J.est)
  Kseq      = as.integer(Kseq)
  multiKmax = as.integer(multiKmax)
  checkoptions = TRUE
  if (checkoptions == TRUE){
    storage.mode(Contrast)  <-"double";
    storage.mode(Kseq)      <-"integer";
    storage.mode(multiKmax) <-"integer";
    .Call("sc_segibp",Contrast,Kseq,multiKmax)
  }
}


EMalgo <- function(x,phi,rupt,P,vh=TRUE){
  checkoptions = TRUE
  K = dim(rupt)[1]
  
  if (P>K){
    checkoptions = FALSE
    cat("Error in EMalgo : the number of groups must be lower than the number of segments","\n")
  }
  if (checkoptions == TRUE){
    storage.mode(x)<-"double"
    storage.mode(phi) <-"double"
    storage.mode(rupt) <- "double"
    .Call("sc_EMalgo",x,phi,rupt,as.integer(K),as.integer(P),as.logical(vh))
  }
  
}


compactEMalgo <- function(xk,x2k,phi,nk,P,vh=TRUE){
  checkoptions = TRUE
  K = length(xk)
  
  if (P>K){
    checkoptions = FALSE
    cat("Error in EMalgo : the number of groups must be lower than the number of segments","\n")
  }
  if (checkoptions == TRUE){
    storage.mode(xk)<-"double"
    storage.mode(x2k)<-"double"
    storage.mode(nk)<-"double"
    storage.mode(phi) <-"double"
     .Call("sc_compactEMalgo",xk,x2k,phi,nk,as.integer(K),as.integer(P),as.logical(vh))
  }
  
}

quicklvinc <- function(xk,x2k,phi,nk,P,vh=TRUE){
  checkoptions = TRUE
  K = length(xk)
  
  if (P>K){
    checkoptions = FALSE
    cat("Error in EMalgo : the number of groups must be lower than the number of segments","\n")
  }
  if (checkoptions == TRUE){
    storage.mode(xk)<-"double"
    storage.mode(x2k)<-"double"
    storage.mode(nk)<-"double"
    storage.mode(phi) <-"double"
     .Call("sc_Lvinc",xk,x2k,phi,nk,as.integer(K),as.integer(P),as.logical(vh))
  }
  
}


compactEMinit <- function(xk,x2k,nk,P,R_OMP_NUM_THREADS, vh=TRUE){
  checkoptions = TRUE
  K = length(xk)
  if (P>K){
    checkoptions = FALSE
    cat("Error in EMinit : the number of groups must be lower than the number of segments","\n")
  }
  if (checkoptions == TRUE){
    storage.mode(xk)<-"double"
    storage.mode(x2k)<-"double"
    storage.mode(nk)<-"double"    
    .Call("sc_compactEMinit",xk,x2k,nk,as.integer(K),as.integer(P),as.integer(R_OMP_NUM_THREADS),as.logical(vh))
  }
  
}

##compactEMinitReg <- function(xk,x2k,yk,y2k,xyk,nk,P,R_OMP_NUM_THREADS, vh=TRUE){
##  checkoptions = TRUE
##  K = length(xk)
##  if (P>K){
##    checkoptions = FALSE
##    cat("Error in EMinit : the number of groups must be lower than the number of segments","\n")
##  }
##  if (checkoptions == TRUE){
##    storage.mode(xk)<-"double"
##    storage.mode(x2k)<-"double"
##    storage.mode(yk)<-"double"
##    storage.mode(y2k)<-"double"
##    storage.mode(xyk)<-"double"
##    storage.mode(nk)<-"double"    
##    .Call("sc_compactEMinit",xk,x2k,yk,y2k,xyk,nk,as.integer(K),as.integer(P),as.integer(R_OMP_NUM_THREADS),as.logical(vh))
##  }  
##}




EMinit <- function(x,rupt,P,vh=TRUE){
  checkoptions = TRUE
  K = dim(rupt)[1]
  if (P>K){
    checkoptions = FALSE
    cat("Error in EMinit : the number of groups must be lower than the number of segments","\n")
  }
  if (checkoptions == TRUE){
    storage.mode(x)<-"double"
    storage.mode(rupt) <- "double"
    .Call("sc_EMinit",x,rupt,as.integer(K),as.integer(P),as.logical(vh))
  }
  
}


