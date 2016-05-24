## reconstruct.R (2014-10-24)

##   Ancestral Character Estimation

## Copyright 2014 Manuela Royer-Carenzi, Didier Gilles

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

#renvoie la racine d'arbre
racine <- function(arbre) {
    Ntip(arbre) + 1
}

# renvoie une liste dont la premiere composante est l'arbre renumerote
# de telle sorte que l'index d'un enfant est superieur a celui de son pere,
# la seconde compopsante est la fonction de l'index initial vers le second,
# et la troisieme son inverse
# (attention probleme pour l'image de 0 mise a l'image du max)
#
renumeroteArbre <- function(arbre) {
  m <- Ntip(arbre) + Nnode(arbre)
  v<-numeric(m)
  t<-numeric(m)
  stack<-numeric(m)
  istack<-1
  stack[istack]<-racine(arbre)
  codeI<-1
  codeL<-Nnode(arbre)+1
  while(istack>0){
    cour<-stack[istack]
    istack<-istack-1
    l <- which(arbre$edge[, 1] == cour)
    if(length(l)>0){
      v[cour] <- codeI
      t[codeI] <- cour
      codeI <- codeI+1
      for(i in 1:length(l)) {
	istack<-istack+1
	stack[istack] <- arbre$edge[l[i], 2]
      }
    } else {
      v[cour] <- codeL
      t[codeL] <- cour
      codeL <- codeL+1
    }
  }
  arbrebis<-arbre
#renumeroter les noms
  for(i in 1:Nedge(arbre)) {
    arbrebis$edge[i,1] <- v[arbre$edge[i,1]]
    arbrebis$edge[i,2] <- v[arbre$edge[i,2]]
  }
  l <- list(arbre = arbrebis, cod = v, dec = t)
  l
}

#calcule la matrice C selon le modele BM
#
calculeC <- function(arbre) {
  m <- Ntip(arbre) + Nnode(arbre)
  C <- matrix(0,nrow=m,ncol=m)
  for(i in 1:(m)) {
    l <- which(arbre$edge[, 2] == i)
    if(length(l)>0){
      for(j in 1:(m)) {
	C[j,i] <- C[j, arbre$edge[l[1], 1]]
      }
    }
    C[i,i]<-1;
  }
  t(C)
}





#########################
#calcul Variance
#########################

getSumSquare <- function(value, arbre) {
 sum <- 0.
 for(eu in 1:Nedge(arbre))
   sum <- sum + (value[arbre$edge[eu,2]]-value[arbre$edge[eu,1]])^2/arbre$edge.length[eu]
 sum
}


getMLHessian <- function(value, arbre) {
   sumSqu <- getSumSquare(value, arbre)
   nI <- Nnode(arbre)
   nT <- length(arbre$tip.label)
   nE <- nI+nT-1
   sizeH<-nI+1
   hessian <- matrix(0., nrow=sizeH, ncol=sizeH)
   var <- sumSqu/nE
   sd <- sqrt(var)
   hessian[1,1] <- -nE/(2*var^2)+sumSqu/var^3
   for(i in 1:nI) {
     	child <- which(arbre$edge[, 1] == nT+i)
  	if(length(child)>0) {
      		for(j in 1:length(child)) {
     		hessian[1,i+1] <- hessian[1,i+1]-(value[arbre$edge[child[j],2]]-value[nT+i])/arbre$edge.length[child[j]]
      		hessian[i+1,i+1] <- hessian[i+1,i+1]+1./arbre$edge.length[child[j]]
       		if(arbre$edge[child[j],2]>nT)
		  hessian[i+1,arbre$edge[child[j],2]-nT+1] <- -1./(var*arbre$edge.length[child[j]])
      		}
     	}
      	anc <- which(arbre$edge[, 2] == nT+i)
     	if(length(anc)>0) {
     	 	for(j in 1:length(anc)) {
       		hessian[1,i+1] <- hessian[1,i+1]+(value[nT+i]-value[arbre$edge[anc[j],1]])/arbre$edge.length[anc[j]]
     		hessian[i+1,i+1] <- hessian[i+1,i+1]+1./arbre$edge.length[anc[j]]
       		hessian[i+1,arbre$edge[anc[j],1]-nT+1] <- -1./(var*arbre$edge.length[anc[j]])
      		}
     	}
   hessian[1,i+1] <- -hessian[1,i+1]/sd^2
   hessian[i+1,1] <- hessian[1,i+1]
   hessian[i+1,i+1] <- hessian[i+1,i+1]/var
   }
   hessian
}

getREMLHessian <- function(value, arbre, sigma2) {
   nI <- Nnode(arbre)
   nT <- length(arbre$tip.label)
   sizeH<-nI
   hessian <- matrix(0., nrow=sizeH, ncol=sizeH)
   for(i in 1:nI) {
     	child <- which(arbre$edge[, 1] == nT+i)
  	if(length(child)>0) {
      		for(j in 1:length(child)) {
      		hessian[i,i] <- hessian[i,i]+1./arbre$edge.length[child[j]]
       		if(arbre$edge[child[j],2]>nT)
		  hessian[i,arbre$edge[child[j],2]-nT] <- -1./(sigma2*arbre$edge.length[child[j]])
      		}
     	}
      	anc <- which(arbre$edge[, 2] == nT+i)
     	if(length(anc)>0) {
     	 	for(j in 1:length(anc)) {
      		hessian[i,i] <- hessian[i,i]+1./arbre$edge.length[anc[j]]
       		hessian[i,arbre$edge[anc[j],1]-nT] <- -1./(sigma2*arbre$edge.length[anc[j]])
      		}
     	}
      hessian[i,i] <- hessian[i,i]/sigma2
   }
   hessian
}



reconstruct <- function (x, phyInit, method = "ML", CI = TRUE) {
 if (!inherits(phyInit, "phylo"))
  stop("object \"phy\" is not of class \"phylo\"")
 if (is.null(phyInit$edge.length))
  stop("tree has no branch lengths")
 nN <- phyInit$Nnode
 nT <- length(x)
#renumber tree
 transf <- renumeroteArbre(phyInit)
 phy <- transf$arbre
 vY <- x
 for (iF in 1:nT) {
  indT <- transf$cod[iF] - nN
  vY[indT] <- x[iF]
  names(vY)[indT]  <- names(x)[iF]
 }
 LongBranche <- phy$edge.length
 Fils <- phy$edge[,2]
 Tau <- phy$edge.length[order(phy$edge[,2])]
 vJ <- rep(1, length=nT)
 sigmaAcc <- Tau
 V2Acc <- diag(sigmaAcc)
 CABM <- calculeC(phy)
 C <- CABM[2:(nN+nT), 2:(nN+nT)]
 sigmaNoeuds <- C %*% V2Acc %*% t(C)
 sigma11 <- sigmaNoeuds[(1:(nN-1)), (1:(nN-1))]
 sigma22 <- sigmaNoeuds[(nN:(nN+nT-1)) , (nN:(nN+nT-1)) ]
 sigma12 <- sigmaNoeuds[(1:(nN-1)), (nN:(nN+nT-1)) ]
 sigma21 <- sigmaNoeuds[(nN:(nN+nT-1)) , (1:(nN-1))]
 GM <- solve( t(vJ) %*% solve(sigma22) %*% vJ ) %*% t(vJ) %*% solve(sigma22) %*% vY
 ZA <- rep(GM, length=nN-1)+ sigma12 %*% solve(sigma22) %*% (vY-rep(GM, length=nT))
 Intern <- c(GM, ZA)
 ValueTmp <- c(Intern, vY)
 Value <- ValueTmp
 for(iF in 1:(nN+nT)) {
  Value[transf$dec[iF]] <- ValueTmp[iF]
 }
 switch(method,
  ML = {
   Hessian <- getMLHessian(Value, phyInit)
   se <- sqrt(diag(solve(Hessian)))
   se <- se[-1]
   tmp <- se*qt(0.025, nN)
  },
  REML={
   minusLogLik <- function(sig2) {
    if (sig2 < 0) return(1e+100)
    V <- sig2 * vcv(phyInit)
    distval <- stats::mahalanobis(x, center = mu, cov = V)
    logdet <- sum(log(eigen(V, symmetric = TRUE, only.values = TRUE)$values))
    (nT * log(2 * pi) + logdet + distval)/2
   }
   mu <- rep(GM, nT)
   out <- nlm(minusLogLik, 1, hessian = FALSE)
   sigma2 <- out$estimate
   Hessian <- getREMLHessian(Value, phyInit, sigma2)
   se <- sqrt(diag(solve(Hessian)))
   tmp <- se*qt(0.025, nN)
  },
  GLS = {
   Hessian <- (sigma11 - sigma12 %*% solve(sigma22) %*% t(sigma12))
   seTmp <- c(0, sqrt(diag(Hessian)))
   se <- seTmp
   for(iF in 1:nN) {
    se[iF] <- seTmp[transf$cod[iF+nT]]
   }
   tmp <- se*qnorm(0.025)
  }
 )
 InternOP <- Intern
 for (iF in 1:nN) {
  InternOP[iF] <- Value[iF+nT]
 }
 CI95 <- cbind(lower=InternOP+tmp, upper=InternOP-tmp)
if (CI==TRUE)
 list(ace=InternOP, CI95=CI95)
else
  list(ace=InternOP)
}








