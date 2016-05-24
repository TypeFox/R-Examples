Random.Start <- function(k){
  qr.Q(qr(matrix(rnorm(k*k),k)))
  }

NormalizingWeight <- function(A, normalize=FALSE){
 if ("function" == mode(normalize)) normalize <- normalize(A)
 if (is.logical(normalize)){
    if (normalize) normalize <- sqrt(rowSums(A^2))
    else return(array(1, dim(A)))
    }
 if (is.vector(normalize)) 
    {if(nrow(A) != length(normalize))
        stop("normalize length wrong in NormalizingWeight")
     return(array(normalize, dim(A)))
    }
 stop("normalize argument not recognized in NormalizingWeight")
}

print.GPArotation <- function (x, digits=3, Table=FALSE, ...){
   cat(if(x$orthogonal)"Orthogonal" else "Oblique")
   cat(" rotation method", x$method)
   cat((if(!x$convergence)" NOT" ), "converged.\n")

   cat("Loadings:\n")
   print(x$loadings, digits=digits)

   cat("\nRotating matrix:\n")
   print(t(solve(x$Th)), digits=digits)

   if(!x$orthogonal){
     cat("\nPhi:\n")
     print(x$Phi, digits=digits)
     }

   if(Table){
     cat("\nIteration table:\n")
     print(x$Table, digits=digits)
     }
   invisible(x)
   }

summary.GPArotation <- function (object, ...){ 
   r <- list(loadings=object$loadings, method=object$method, orthogonal=object$orthogonal, 
          convergence=object$convergence, iters=nrow(object$Table))
   class(r) <- "summary.GPArotation"
   r
   }

print.summary.GPArotation <- function (x, digits=3, ...){
   cat(if(x$orthogonal)"Orthogonal" else "Oblique")
   cat(" rotation method", x$method)
   if(!x$convergence) cat(" NOT" )
   cat(" converged in ", x$iter, " iterations.\n")

   cat("Loadings:\n")
   print(x$loadings, digits=digits)
  }

GPForth <- function(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, 
		    method="varimax", methodArgs=NULL){
 if((!is.logical(normalize)) || normalize) {
     W <- NormalizingWeight(A, normalize=normalize)
     normalize <- TRUE
     A <- A/W
     }
 if(1 >= ncol(A)) stop("rotation does not make sense for single factor models.")
 al <- 1
 L <- A %*% Tmat
 #Method <- get(paste("vgQ",method,sep="."))
 #VgQ <- Method(L, ...)
 Method <- paste("vgQ",method,sep=".")
 VgQ <- do.call(Method, append(list(L), methodArgs))
 G <- crossprod(A,VgQ$Gq)
 f <- VgQ$f
 Table <- NULL
 #set initial value for the unusual case of an exact initial solution 
 VgQt <- do.call(Method, append(list(L), methodArgs))   
 for (iter in 0:maxit){
   M <- crossprod(Tmat,G)
   S <- (M + t(M))/2
   Gp <- G - Tmat %*% S
   s <- sqrt(sum(diag(crossprod(Gp))))
   Table <- rbind(Table, c(iter, f, log10(s), al))
   if (s < eps)  break
   al <- 2*al
   for (i in 0:10){
     X <- Tmat - al * Gp
     UDV <- svd(X)
     Tmatt <- UDV$u %*% t(UDV$v)
     L <- A %*% Tmatt
     #VgQt <- Method(L, ...)
     VgQt <- do.call(Method, append(list(L), methodArgs))
     if (VgQt$f < (f - 0.5*s^2*al)) break
     al <- al/2
     }
   Tmat <- Tmatt
   f <- VgQt$f
   G <- crossprod(A,VgQt$Gq)
   }
 convergence <- (s < eps)
 if ((iter == maxit) & !convergence)
     warning("convergence not obtained in GPForth. ", maxit, " iterations used.")
 if(normalize) L <- L * W
 dimnames(L) <- dimnames(A)
 r <- list(loadings=L, Th=Tmat, Table=Table, 
        method=VgQ$Method, orthogonal=TRUE, convergence=convergence, Gq=VgQt$Gq)
 class(r) <- "GPArotation"
 r
}

GPFoblq <- function(A, Tmat=diag(ncol(A)), normalize=FALSE, eps=1e-5, maxit=1000, 
		    method="quartimin", methodArgs=NULL){
 if(1 >= ncol(A)) stop("rotation does not make sense for single factor models.")
 if((!is.logical(normalize)) || normalize) {
     W <- NormalizingWeight(A, normalize=normalize)
     normalize <- TRUE
     A <- A/W
     }
 al <- 1
 L <- A %*% t(solve(Tmat))
 #Method <- get(paste("vgQ",method,sep="."))
 #VgQ <- Method(L, ...)
 Method <- paste("vgQ",method,sep=".")
 VgQ <- do.call(Method, append(list(L), methodArgs))
 G <- -t(t(L) %*% VgQ$Gq %*% solve(Tmat))
 f <- VgQ$f
 Table <- NULL
 #Table <- c(-1,f,log10(sqrt(sum(diag(crossprod(G))))),al)
 #set initial value for the unusual case of an exact initial solution 
 VgQt <- do.call(Method, append(list(L), methodArgs))   
 for (iter in 0:maxit){
   Gp <- G - Tmat %*% diag(c(rep(1,nrow(G)) %*% (Tmat*G)))
   s <- sqrt(sum(diag(crossprod(Gp))))
   Table <- rbind(Table,c(iter,f,log10(s),al))
   if (s < eps) break
   al <- 2*al
   for (i in 0:10){
     X <- Tmat - al*Gp
     v <- 1/sqrt(c(rep(1,nrow(X)) %*% X^2))
     Tmatt <- X %*% diag(v)
     L <- A %*% t(solve(Tmatt))
     #VgQt <- Method(L, ...)
     VgQt <- do.call(Method, append(list(L), methodArgs))
     improvement <- f - VgQt$f 
     if (improvement >  0.5*s^2*al) break
     al <- al/2
     }
   Tmat <- Tmatt
   f <- VgQt$f
   G <- -t(t(L) %*% VgQt$Gq %*% solve(Tmatt))
   }
 convergence <- (s < eps)
 if ((iter == maxit) & !convergence)
     warning("convergence not obtained in GPFoblq. ", maxit, " iterations used.")
 if(normalize) L <- L * W
 dimnames(L) <- dimnames(A)
 # N.B. renaming Lh to loadings in specificific rotations 
 #   uses fact that  Lh is first.
 r <- list(loadings=L, Phi=t(Tmat) %*% Tmat, Th=Tmat, Table=Table,
      method=VgQ$Method, orthogonal=FALSE, convergence=convergence, Gq=VgQt$Gq)
 class(r) <- "GPArotation"
 r
}

#######################


oblimin <- function(L, Tmat=diag(ncol(L)), gam=0, normalize=FALSE, eps=1e-5, maxit=1000){
   GPFoblq(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit, 
                method="oblimin",  methodArgs=list(gam=gam))
   }

vgQ.oblimin <- function(L, gam=0){
  X <- L^2 %*% (!diag(TRUE,ncol(L))) 
  if (0 != gam) {
     p <- nrow(L)
     X <- (diag(1,p) - matrix(gam/p,p,p)) %*% X
     }
  list(Gq=L*X,
       f=sum(L^2 * X)/4,
       Method=  if (gam == 0)  "Oblimin Quartimin" else
		if (gam == .5) "Oblimin Biquartimin" else
		if (gam == 1)  "Oblimin Covarimin"   else
		         paste("Oblimin g=", gam,sep="")  )
}

# original
# vgQ.oblimin <- function(L, gam=0){
#   Method <- paste("Oblimin g=",gam,sep="")
#   if (gam == 0) Method <- "Oblimin Quartimin"
#   if (gam == .5) Method <- "Oblimin Biquartimin"
#   if (gam == 1) Method <- "Oblimin Covarimin"
#   k <- ncol(L)
#   p <- nrow(L)
#   N <- matrix(1,k,k)-diag(k)
#   f <- sum(L^2 * (diag(p)-gam*matrix(1/p,p,p)) %*% L^2 %*% N)/4
#   Gq <- L * ((diag(p)-gam*matrix(1/p,p,p)) %*% L^2 %*% N)
#   return(list(Gq=Gq,f=f,Method=Method))
# }
##vgQ.oblimin(FA2)$f - vgQ.origoblimin(FA2)$f
#vgQ.oblimin(FA2)$Gq - vgQ.origoblimin(FA2)$Gq


quartimin <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000){
   GPFoblq(L, Tmat=Tmat, method="quartimin", normalize=normalize, eps=eps, maxit=maxit)
   }

vgQ.quartimin <- function(L){
  X <-  L^2 %*% (!diag(TRUE,ncol(L))) 
  list(Gq= L*X,
       f= sum(L^2 * X)/4,  
       Method=  "Quartimin" )
  }

#original
#vgQ.quartimin <- function(L){
#  Method="Quartimin"
#  L2 <- L^2
#  k <- ncol(L)
#  M <- matrix(1,k,k)-diag(k)
#  f <- sum(L2 * (L2 %*% M))/4
#  Gq <- L * (L2 %*% M)
#  return(list(Gq=Gq,f=f,Method=Method))
#} 


targetT <- function(L, Tmat=diag(ncol(L)), Target=NULL, normalize=FALSE, eps=1e-5, maxit=1000) {
   if(is.null(Target)) stop("argument Target must be specified.")
   GPForth(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
           method="target", methodArgs=list(Target=Target))
   }

targetQ <- function(L, Tmat=diag(ncol(L)), Target=NULL, normalize=FALSE, eps=1e-5, maxit=1000) {
   if(is.null(Target)) stop("argument Target must be specified.")
   GPFoblq(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
           method="target", methodArgs=list(Target=Target))
   }

vgQ.target <- function(L, Target=NULL){
   if(is.null(Target)) stop("argument Target must be specified.")
   #   e.g.  Target <- matrix(c(rep(NA,4),rep(0,8),rep(NA,4)),8) 
   #  approximates  Michael Brown approach
   Gq <-  2 * (L - Target)
   Gq[is.na(Gq)] <- 0 #missing elements in target do not affect the first derivative 
   list(Gq=Gq,
        f=sum((L-Target)^2, na.rm=TRUE),  
	Method="Target rotation")  #The target rotation ? option in Michael Browne's algorithm should be NA
   }

pstT <- function(L, Tmat=diag(ncol(L)), W=NULL, Target=NULL, normalize=FALSE, eps=1e-5, maxit=1000) {
   if(is.null(W))      stop("argument W must be specified.")
   if(is.null(Target)) stop("argument Target must be specified.")
   GPForth(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit, 
	    method="pst", methodArgs=list(W=W, Target=Target))
   }

pstQ <- function(L, Tmat=diag(ncol(L)), W=NULL, Target=NULL, normalize=FALSE, eps=1e-5, maxit=1000) {
   if(is.null(W))      stop("argument W must be specified.")
   if(is.null(Target)) stop("argument Target must be specified.")
   GPFoblq(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit, 
	  method="pst", methodArgs=list(W=W, Target=Target))
   }

vgQ.pst <- function(L, W=NULL, Target=NULL){
   if(is.null(W))      stop("argument W must be specified.")
   if(is.null(Target)) stop("argument Target must be specified.")
   # Needs weight matrix W with 1's at specified values, 0 otherwise
   # e.g. W = matrix(c(rep(1,4),rep(0,8),rep(1,4)),8). 
   # When W has only 1's this is procrustes rotation
   # Needs a Target matrix Target with hypothesized factor loadings.
   # e.g. Target = matrix(0,8,2)
   Btilde <- W * Target
   list(Gq= 2*(W*L-Btilde), 
        f = sum((W*L-Btilde)^2),
        Method="Partially specified target")
}

oblimax <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000){
   GPFoblq(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
            method="oblimax")
   }


vgQ.oblimax <- function(L){
  list(Gq= -(4*L^3/(sum(L^4))-4*L/(sum(L^2))),
       f= -(log(sum(L^4))-2*log(sum(L^2))),
       Method="oblimax")
}

entropy <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000) {
   GPForth(L, Tmat=Tmat, method="entropy", normalize=normalize, eps=eps, maxit=maxit)
   }

vgQ.entropy <- function(L){
  list(Gq= -(L*log(L^2 + (L^2==0)) + L),
       f= -sum(L^2*log(L^2 + (L^2==0)))/2, 
       Method="Minimum entropy")
}

quartimax <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000) {
   GPForth(L, Tmat=Tmat, method="quartimax", normalize=normalize, eps=eps, maxit=maxit)
   }

vgQ.quartimax <- function(L){
  list(Gq= -L^3,
       f= -sum(diag(crossprod(L^2)))/4, 
       Method="Quartimax")
}

Varimax <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000) {
   GPForth(L, Tmat=Tmat, method="varimax", normalize=normalize, eps=eps, maxit=maxit)
   }

vgQ.varimax <- function(L){
  QL <- sweep(L^2,2,colMeans(L^2),"-")
  list(Gq= -L * QL,
       f= -sqrt(sum(diag(crossprod(QL))))^2/4, 
       Method="varimax")
}

simplimax <- function(L, Tmat=diag(ncol(L)), k=nrow(L), normalize=FALSE, eps=1e-5, maxit=1000) {
   GPFoblq(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
         method="simplimax", methodArgs=list(k=k))
   }

vgQ.simplimax <- function(L, k=nrow(L)){
  # k: Number of close to zero loadings
  Imat <- sign(L^2 <= sort(L^2)[k])
  list(Gq= 2*Imat*L,
       f= sum(Imat*L^2), 
       Method="Simplimax")
}

bentlerT <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000) {
   GPForth(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit, 
	     method="bentler")
   }

bentlerQ <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000) {
   GPFoblq(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
           method="bentler")
   }

vgQ.bentler <- function(L){
  L2 <- L^2
  M <- crossprod(L2)
  D <- diag(diag(M))
  list(Gq= -L * (L2 %*% (solve(M)-solve(D))),
       f= -(log(det(M))-log(det(D)))/4,
       Method="Bentler's criterion")
}

tandemI <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000) {
   GPForth(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
            method="tandemI")
   }

#vgQ.tandemI <- function(L){  # Tandem Criterion, Comrey, 1967.
#  Method <- "Tandem I"
#  LL <- (L %*% t(L))
#  LL2 <- LL^2
#  f <- -sum(diag(crossprod(L^2, LL2 %*% L^2)))
#  Gq1 <- 4 * L *(LL2 %*% L^2)
#  Gq2 <- 4 * (LL * (L^2 %*% t(L^2))) %*% L
#  Gq <- -Gq1 - Gq2 
#  return(list(Gq=Gq,f=f,Method=Method))
#}

vgQ.tandemI <- function(L){  # Tandem Criterion, Comrey, 1967.
  LL <- (L %*% t(L))
  LL2 <- LL^2
  Gq1 <- 4 * L *(LL2 %*% L^2)
  Gq2 <- 4 * (LL * (L^2 %*% t(L^2))) %*% L
  Gq <- -Gq1 - Gq2 
  list(Gq=Gq,
       f= -sum(diag(crossprod(L^2, LL2 %*% L^2))), 
       Method="Tandem I")
  }

tandemII <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000) {
   GPForth(L, Tmat=Tmat, method="tandemII", normalize=normalize, eps=eps, maxit=maxit)
   }


vgQ.tandemII <- function(L){  # Tandem Criterion, Comrey, 1967.
  LL <- (L %*% t(L))
  LL2 <- LL^2
  f <- sum(diag(crossprod(L^2, (1-LL2) %*% L^2)))
  Gq1 <- 4 * L *((1-LL2) %*% L^2)
  Gq2 <- 4 * (LL * (L^2 %*% t(L^2))) %*% L
  Gq <- Gq1 - Gq2 
  list(Gq=Gq,
       f=f, 
       Method="Tandem II")
  }

geominT <- function(L, Tmat=diag(ncol(L)), delta=.01, normalize=FALSE, eps=1e-5, maxit=1000){
   GPForth(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
           method="geomin", methodArgs=list(delta=delta))
   }

geominQ <- function(L, Tmat=diag(ncol(L)), delta=.01, normalize=FALSE, eps=1e-5, maxit=1000){
   GPFoblq(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
            method="geomin", methodArgs=list(delta=delta))
   }

vgQ.geomin <- function(L, delta=.01){
  k <- ncol(L)
  p <- nrow(L)
  L2 <- L^2 + delta
  pro <- exp(rowSums(log(L2))/k) 
  list(Gq=(2/k)*(L/L2)*matrix(rep(pro,k),p),
       f= sum(pro), 
       Method="Geomin")
  }

cfT <- function(L, Tmat=diag(ncol(L)), kappa=0, normalize=FALSE, eps=1e-5, maxit=1000) {
   GPForth(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
             method="cf", methodArgs=list(kappa=kappa))
   }

cfQ <- function(L, Tmat=diag(ncol(L)), kappa=0, normalize=FALSE, eps=1e-5, maxit=1000) {
   GPFoblq(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
             method="cf", methodArgs=list(kappa=kappa))
   }

vgQ.cf <- function(L, kappa=0){
  k <- ncol(L)
  p <- nrow(L)
  # kappa <- 0 # quartimax 
  # kappa <- 1/p # varimax
  # kappa <- m/(2*p) # equamax
  # kappa <- (m-1)/(p+m-2) # parsimax
  # kappa <- 1 # factor parsimony
  N <- matrix(1,k,k)-diag(k)
  M <- matrix(1,p,p)-diag(p)
  L2 <- L^2
  f1 <- (1-kappa)*sum(diag(crossprod(L2,L2 %*% N)))/4
  f2 <- kappa*sum(diag(crossprod(L2,M %*% L2)))/4
  list(Gq= (1-kappa) * L * (L2 %*% N) + kappa * L * (M %*% L2),
       f= f1 + f2,
       Method=paste("Crawford-Ferguson:k=",kappa,sep=""))
}

infomaxT <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000) {
   GPForth(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
            method="infomax")
   }

infomaxQ <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000) {
   GPFoblq(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
            method="infomax")
   }

vgQ.infomax <- function(L){
  Method <- "Infomax"
  k <- ncol(L)
  p <- nrow(L)
  S <- L^2
  s <- sum(S)
  s1 <- rowSums(S)
  s2 <- colSums(S)
  E <- S/s
  e1 <- s1/s
  e2 <- s2/s
  Q0 <- sum(-E * log(E))
  Q1 <- sum(-e1 * log(e1))
  Q2 <- sum(-e2 * log(e2))
  f <- log(k) + Q0 - Q1 - Q2
  H <- -(log(E) + 1)
  alpha <- sum(S * H)/s^2
  G0 <- H/s - alpha * matrix(1, p, k)
  h1 <- -(log(e1) + 1)
  alpha1 <- s1 %*% h1/s^2
  G1 <- matrix(rep(h1,k), p)/s - as.vector(alpha1) * matrix(1, p, k)
  h2 <- -(log(e2) + 1)
  alpha2 <- h2 %*% s2/s^2
  G2 <- matrix(rep(h2,p), ncol=k, byrow=T)/s - as.vector(alpha2) * matrix(1, p, k)
  Gq <- 2 * L * (G0 - G1 - G2)
  list(Gq=Gq,f=f,Method=Method)
}

mccammon <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, eps=1e-5, maxit=1000) {
   GPForth(L, Tmat=Tmat, method="mccammon", normalize=normalize, eps=eps, maxit=maxit)
   }

vgQ.mccammon <- function(L){
  Method <- "McCammon entropy"
  k <- ncol(L)
  p <- nrow(L)
  S <- L^2
  M <- matrix(1,p,p)
  s2 <- colSums(S)
  P <- S / matrix(rep(s2,p),ncol=k,byrow=T)
  Q1 <- -sum(P * log(P))
  H <- -(log(P) + 1)
  R <- M %*% S
  G1 <- H/R - M %*% (S*H/R^2)
  s <- sum(S)
  p2 <- s2/s
  Q2 <- -sum(p2 * log(p2))
  h <- -(log(p2) + 1)
  alpha <- h %*% p2
  G2 <- rep(1,p) %*% t(h)/s - as.vector(alpha)*matrix(1,p,k)
  Gq <- 2*L*(G1/Q1 - G2/Q2)
  Q <- log(Q1) - log(Q2)
  list(Gq=Gq, f=Q, Method=Method)
}


bifactorT <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, 
     eps=1e-5, maxit=1000){
  #adapted from Jennrich and Bentler 2011. code provided by William Revelle	
    GPForth(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
            method="bifactor")
    }

bifactorQ <- function(L, Tmat=diag(ncol(L)), normalize=FALSE, 
     eps=1e-5, maxit=1000){
  #the oblique case
  #adapted from Jennrich and Bentler 2011. code provided by William Revelle	
    GPFoblq(L, Tmat=Tmat, normalize=normalize, eps=eps, maxit=maxit,
            method="bifactor")
    }


vgQ.bifactor <- function(L) {
  # code provided by William Revelle	
  D <- function(L) {
    L2 <- L^2
    L2N <- L2 %*% ! diag(NCOL(L))
    list(f=sum(L2 * L2N),
         Gq=4 * L * L2N)
  }
  lvg <- D(L[,-1, drop=FALSE])
  G <- lvg$Gq
  G <-cbind(G[,1],G)
  G[,1] <- 0
  list(f=lvg$f,
       Gq=G)
}


# promax is already defined in the stats (previously mva) package
# 
#GPromax <- function(A,pow=3){
# method <- "Promax"
# # Initial rotation: Standardized Varimax
# require(statsa)
# xx <- promax(A,pow)
# Lh <- xx$loadings
# Th <- xx$rotmat
# orthogonal <- F
# Table <- NULL
#return(list(loadings=Lh,Th=Th,Table=NULL,method,orthogonal=orthogonal))
#}


