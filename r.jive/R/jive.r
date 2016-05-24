jive <- function (data, rankJ=1, rankA=rep(1,length(data)), method="perm", dnames=names(data), conv="default", maxiter=1000, scale=TRUE, orthIndiv=TRUE, est=TRUE, showProgress=TRUE) {
   # Get the number of data sets
  
   l <- length(data)

   # Calculate the dimensions of the data (for BIC)
   # n is the total number of values in each reduced individual matrix (=nd[k])
   # d is the number of rows of each dataset
   n <- c()
   d <- c()
   for (i in 1:length(data)) {
      n[i] <- nrow(data[[i]]) * ncol(data[[i]])
      d[i] <- nrow(data[[i]])
   }

   # Impute missing values using SVDmiss
   for(i in 1:l) {
      temp <- SVDmiss(data[[i]], ncomp=ncol(data[[i]]))[[1]]
      data[[i]] <- temp$u %*% diag(x=temp$d) %*% t(temp$v)
   }

   # Center and scale individual data sets by frobenius norm
   for (i in 1:l) {
      data[[i]] <- data[[i]] - matrix(rep(apply(data[[i]],1,mean,na.rm=T),ncol(data[[i]])),nrow=nrow(data[[i]]))
      if (scale) { data[[i]] <- data[[i]] / norm(data[[i]], type="f")*sqrt(sum(n))}
   }
   
   if(conv=="default"){conv = 10^(-6)*norm(do.call(rbind, data),type="f")}

   # Centered and scaled (but not reduced) data
   orig <- data

   # Call JIVE algorithm for different rank methods
   if (method=="given") {
      # Reduce matrix to speed computations
      # Should make reduction a separate function
      if (est) {
         u <- list()
         for(i in 1:l) {
            if(nrow(data[[i]]) > ncol(data[[i]])) {
               temp <- svdwrapper(data[[i]], nu=ncol(data[[i]]), nv=ncol(data[[i]]))
               data[[i]] <- diag(x=temp$d[1:ncol(data[[1]])], nrow=ncol(data[[1]])) %*% t(temp$v[,1:ncol(data[[1]])])
               u[[i]] <- temp$u
            } else {
               # Set u[[i]] to identity matrix
               u[[i]] <- diag(1, nrow(data[[i]]))
            }
         }
      }
      if (showProgress) { cat("Running JIVE algorithm for ranks:\njoint rank:", rankJ,", individual ranks:", rankA, "\n") }
      temp <- jive.iter(data, rankJ, rankA, conv=conv, maxiter=maxiter, orthIndiv=orthIndiv, showProgress=showProgress)
      joint <- temp$joint
      individual <- temp$individual
      if (est) {
         for (i in 1:l) {
            joint[[i]] <- u[[i]] %*% joint[[i]]
            individual[[i]] <- u[[i]] %*% individual[[i]]
         }
      }
   } else if (method=="perm") {
      temp <- jive.perm(data, est=est, conv=conv, maxiter=maxiter, orthIndiv=orthIndiv, showProgress=showProgress)
      joint <- temp$joint
      individual <- temp$individual
      rankJ <- temp$rankJ
      rankA <- temp$rankA
      converged <- temp$converged
   } else if (method=="bic") {
      # Reduce matrix to speed computations
      if (est) {
         u <- list()
         for(i in 1:l) {
            if(nrow(data[[i]]) > ncol(data[[i]])) {
               temp <- svdwrapper(data[[i]], nu=ncol(data[[i]]), nv=ncol(data[[i]]))
               data[[i]] <- diag(x=temp$d[1:ncol(data[[1]])], nrow=ncol(data[[1]])) %*% t(temp$v[,1:ncol(data[[1]])])
               u[[i]] <- temp$u
            } else {
               # Set u[[i]] to identity matrix
               u[[i]] <- diag(1, nrow(data[[i]]))
            }
         }
      }
      temp <- bic.jive(data, n, d, conv=conv, maxiter=maxiter, orthIndiv=orthIndiv, showProgress=showProgress)
      joint <- temp$joint
      individual <- temp$individual
      rankJ <- temp$rankJ
      rankA <- temp$rankA
      bic.table <- temp$bic.table
      if (est) {
         for (i in 1:l) {
            joint[[i]] <- u[[i]] %*% joint[[i]]
            individual[[i]] <- u[[i]] %*% individual[[i]]
         }
      }
   }

   # Add names to data
   if (is.null(dnames)) {
      for (i in 1:l) {
         names(orig)[[i]] <- paste("Source",i,sep="_")
      }
   } else {
      names(orig) <- dnames
   }

   result <- list(data=orig, joint=joint, individual=individual, rankJ=rankJ, rankA=rankA, method=method)
   if (method=="bic") { 
      result$bic.table <- bic.table
   }
   if (method=="perm") { 
      result$converged <- converged
   }
   class(result) <- "jive"
   return(result)
}


# Run one iteration of the JIVE algorithm
jive.iter <- function (data, rankJ=1, rankA=rep(1,length(data)), conv=0.000001, maxiter=1000, orthIndiv=TRUE, showProgress=TRUE) {
   # Get the number of data sets
   l <- length(data)

   # Initialize A matrix
   A <- list()
   for (i in 1:l) {
      A[[i]] <- matrix(0, nrow(data[[i]]), ncol(data[[i]]))
   }

   Xtot <- do.call(rbind, data)
   Jtot <- matrix(-1, nrow(Xtot), ncol(Xtot))
   Atot <- matrix(-1, nrow(Xtot), ncol(Xtot))

   # Keep track of number of iterations and convergence criteria, respectively
   nrun = 0
   converged = F

   if (orthIndiv) { Vind <- list() }

   # Iterate to find A and J
   while (nrun < maxiter & !converged) {
      # Store previous iteration
      Jlast <- Jtot
      Alast <- Atot

     timeJ <- system.time({
      # For fixed A, calculate J
      if (rankJ > 0) { # Only change J if rankJ is non-zero
         temp <- Xtot - Atot
         s <- svdwrapper(temp, nu=rankJ, nv=rankJ)
         Jtot <- s$u[,1:rankJ] %*% diag(x=s$d[1:rankJ], nrow=rankJ) %*% t(s$v[,1:rankJ])
         V <- s$v[,1:rankJ]
      } else {
         Jtot <- matrix(0, nrow(Xtot), ncol(Xtot))
         V <- matrix(0,ncol(Xtot),rankJ)
      }
      temp <- Jtot
      J <- list()
      for (i in 1:l) {
         J[[i]] <- temp[1:nrow(data[[i]]),]
         temp <- temp[-(1:nrow(data[[i]])),]
      }
     })

     timeA <- system.time({
      # For fixed J, calculate A
      A <- list()
      for (i in 1:l) {
         # If rank is 0, set to 0 matrix.  Otherwise, estimate.
         if (rankA[i] > 0) {
            temp <- (data[[i]] - J[[i]]) %*% (diag(ncol(Xtot)) - V %*% t(V))
            if (orthIndiv & nrun > 0) {
               for (j in (1:l)[-i]) {
                  temp <- temp %*% (diag(ncol(Xtot)) - Vind[[j]] %*% t(Vind[[j]]))
               }
            }
            s <- svdwrapper(temp, nu=rankA[i], nv=rankA[i])
            if (orthIndiv) { Vind[[i]] <- s$v[,1:rankA[i]] }
            A[[i]] <- s$u[,1:rankA[i]] %*% diag(x=s$d[1:rankA[i]], nrow=rankA[i]) %*% t(s$v[,1:rankA[i]])
         } else {
            A[[i]] <- matrix(0, nrow(data[[i]]), ncol(data[[i]]))
            if (orthIndiv) { Vind[[i]] <- matrix(0,ncol(Xtot),rankA[[i]]) }
         }
      }
     })

      if (orthIndiv & nrun == 0) {
         for (i in 1:l) {
            for (j in (1:l)[-i]) {
               A[[i]] <- A[[i]] %*% (diag(ncol(Xtot)) - Vind[[j]] %*% t(Vind[[j]]))
            }}
         for(i in 1:l){     ########re-estimate the Vind's
            if (rankA[i] > 0) {
             s <- svdwrapper(A[[i]], nu=rankA[i], nv=rankA[i])
             Vind[[i]] <- s$v[,1:rankA[i]]
            }
         }
      }

      Atot <- do.call(rbind, A)

      if (norm(Jtot - Jlast, type="f") <= conv & norm(Atot - Alast, type="f") <= conv) {
         converged <- T
      }
      nrun = nrun + 1
   }

   # Export data as list (split for each data source)
   # Find breaks between data sources
   #d_red <- c()
   #for (i in 1:length(data)) { d_red[i] <- nrow(data[[i]]) }
   #dlow <- c(1, cumsum(d_red)[-length(data)]+1)
   #dhi <- cumsum(d_red)

   # Split results into lists
   #joint <- list()
   #individual <- list()
   #for (i in 1:l) {
   #   joint[[i]] <- Jtot[dlow[i]:dhi[i],]
   #   individual[[i]] <- Atot[dlow[i]:dhi[i],]
   #}
   
   if (showProgress) {
    if (converged) { 
      cat(paste("JIVE algorithm converged after ", nrun, " iterations.\n")) 
    } else {
      cat(paste("JIVE algorithm did not converge after ", nrun, " iterations.\n")) 
    }
   }

   return(list(data=data, joint=J, individual=A, rankJ, rankA, method="given"))
}

# Function to simplify p.jive calculation
pjsum <- function (dim,rank) {
   s<-0
   for(i in 1:length(dim)) {
      s=s+ifelse(rank[i]==0,0,sum(dim[i]:(dim[i]-rank[i]+1)))
   }
   return(s)
}

# Determine ranks with bic
bic.jive <- function (data, n=unlist(lapply(data,ncol))*unlist(lapply(data,nrow)), d=unlist(lapply(data,nrow)), conv=.000001, maxiter=1000, orthIndiv=TRUE, showProgress=TRUE) {
   # Get the number of data sets
   l <- length(data)

   nc <- ncol(data[[1]])

   # Set coefficient of p.jive
   lambda <- log(sum(n))

   # sse will store sum of squared error for each individual structure matrix
   sse <- c()
   # Find breaks between data sources
   #d_red <- c()
   #for (i in 1:length(data)) { d_red[i] <- nrow(data[[i]]) }
   #dlow <- c(1, cumsum(d_red)[-length(data)]+1)
   #dhi <- cumsum(d_red)
   Xtot <- do.call(rbind, data)

   bic.improve <- T
   rankJ <- 0
   rankA <- rep(0,l)
   if (showProgress) { cat("Running JIVE algorithm for ranks:\njoint rank:", rankJ,", individual ranks:", rankA, "\n") }
   current <- jive.iter(data, rankJ, rankA, conv, maxiter, orthIndiv, showProgress=showProgress)
    # Calculate BIC
    for (k in 1:length(data)) { sse[k] <- norm(data[[k]] - current$joint[[k]] - current$individual[[k]], type="f")^2 }
    p.jive <- 0
    current.bic <-  sum(n*log(sse/n)) + p.jive * lambda
   bic.table <- c(rankJ, rankA, current.bic)

   while (bic.improve) {
      bic.improve <- F
      temp <- list()
      if (showProgress) { cat("Running JIVE algorithm for ranks:\njoint rank:", rankJ+1,", individual ranks:", rankA, "\n") }
      temp[[1]] <- jive.iter(data, rankJ+1, rankA, conv, maxiter, orthIndiv, showProgress=showProgress)
       # Calculate BIC
       for (k in 1:length(data)) { sse[k] <- norm(data[[k]] - temp[[1]]$joint[[k]] - temp[[1]]$individual[[k]], type="f")^2 }
       p.jive <- sum(sum(d):(sum(d)-(rankJ+1)+1)) + sum(nc:(nc-(rankJ+1)+1)) + pjsum(d, rankA) + pjsum(rep(nc,length(data))-(rankJ+1), rankA)
       bic <- sum(n*log(sse/n)) + p.jive * lambda
      bicvec <- bic
      bic.table <- rbind(bic.table, c(rankJ+1, rankA, bic))
      for (i in 1:l) {
        tempR <- rankA
        tempR[i] <- tempR[i] + 1
        if (tempR[i] < min(n,nrow(data[[i]]))) {
         if (showProgress) { cat("Running JIVE algorithm for ranks:\njoint rank:", rankJ,", individual ranks:", tempR, "\n") }
         temp[[i+1]] <- jive.iter(data, rankJ, tempR, conv, maxiter, orthIndiv, showProgress=showProgress)
          # Calculate BIC
          for (k in 1:length(data)) { sse[k] <- norm(data[[k]] - temp[[i+1]]$joint[[k]] - temp[[i+1]]$individual[[k]], type="f")^2 }
          p.jive <- ifelse(rankJ==0,0,sum(sum(d):(sum(d)-rankJ+1)) + sum(nc:(nc-rankJ+1))) + pjsum(d, tempR) + pjsum(rep(nc,length(data))-rankJ, tempR)
          bic <- sum(n*log(sse/n)) + p.jive * lambda
        } else {
         bic <- NA
        }
         bicvec <- c(bicvec, bic)
         bic.table <- rbind(bic.table, c(rankJ, tempR, bic))
      }
      lowest.bic <- temp[[which.min(bicvec)]]
      if (min(bicvec,na.rm=T) < current.bic) {
         bic.improve <- T
         current <- lowest.bic
         current.bic <- min(bicvec,na.rm=T)
         if (which.min(bicvec)==1) {
            rankJ <- rankJ + 1 
         } else {
            rankA[which.min(bicvec)-1] <- rankA[which.min(bicvec)-1] + 1 
         }
      }
   }
   return(list(data=data, joint=current$joint, individual=current$individual, rankJ=rankJ, rankA=rankA, bic.table))
}


# Determine ranks by a permutation test
jive.perm <- function (data, nperms=100, alpha=0.05, est=TRUE, conv=0.000001, maxiter=1000, orthIndiv=TRUE, showProgress=TRUE) {
   nrun <- 0

   # Find breaks between data scources
   #d_red <- c()
   #for (i in 1:length(data)) { d_red[i] <- nrow(data[[i]]) }
   #dlow <- c(1, cumsum(d_red)[-length(data)]+1)
   #dhi <- cumsum(d_red)

  # Initialize J and A matrices
  Jperp <- list()
  Aperp <- list()
  for (i in 1:length(data)) {
     Jperp[[i]] <- matrix(0, nrow=nrow(data[[i]]), ncol=ncol(data[[i]]))
     Aperp[[i]] <- matrix(0, nrow=nrow(data[[i]]), ncol=ncol(data[[i]]))
  }
  last <- rep(-2, length(data)+1)
  current <- rep(-1, length(data)+1)
  while (!isTRUE(all.equal(last, current)) & nrun < 10) {
   last <- current
   if (showProgress) {
      if (nrun == 0) {
         cat("Estimating  joint and individual ranks via permutation...\n")
      } else {
         cat("Re-estimating  joint and individual ranks via permutation...\n")
      }
   }

   # Permute columns of joint structure
   full <- do.call(rbind, data) - do.call(rbind, Aperp)
   n <- ncol(full)
   actual <- svdwrapper(full,nu=0,nv=0)$d
   # Each row of perms will be the singular values of a single permutation
   # The ith column of perms will be the ith singular value from all permutations
   perms <- matrix(NA, nperms, min(n,sum(unlist(lapply(data,nrow)))))
   for (i in 1:nperms) {
      temp <- list()
      for (j in 1:length(data)) {
         temp[[j]] <- data[[j]][, sample(1:n, n, replace=F)]
      }
      full <- do.call(rbind, temp)
      perms[i,] <- svdwrapper(full,nu=0,nv=0)$d
   }
   rankJ <- 0
   for (i in 1:n) {
      if (actual[i] > quantile(perms[,i], 1-alpha)) { 
         rankJ <- rankJ +  1 
      } else {
         # Don't want to count any more rows after one is not greater than the threshold
         break
      }
   }
   # Force joint rank to be non-decreasing
   rankJ <- max(rankJ,last[1])

   # Permute columns of each row of individual structures
   rankA <- c()
   for (i in 1:length(data)) {
      ind <- data[[i]] - Jperp[[i]]
      actual <- svdwrapper(ind,nu=0,nv=0)$d
      perms <- matrix(NA, nperms, min(n,nrow(data[[i]])))
      for (k in 1:nperms) {
         perm <- t(data[[i]])
         pind <- order(c(col(perm)), runif(length(perm))) 
         perm <- matrix(perm[pind], nrow=nrow(ind), ncol=n, byrow=TRUE)

         perms[k,] <- svdwrapper(perm,nu=0,nv=0)$d
      }
      rankA[i] <- 0
      for (j in 1:n) {
         if (actual[j] > quantile(perms[,j], 1-alpha)) { 
            rankA[i] <- rankA[i] +  1 
         } else {
            # Don't want to count any more rows after one is not greater than the threshold
            break
         }
      }
   }

   dataR <- list()
   if (est) {
      u <- list()
      for(i in 1:length(data)) {
         if(nrow(data[[i]]) > ncol(data[[i]])) {
            temp <- svdwrapper(data[[i]], nu=ncol(data[[i]]), nv=ncol(data[[i]]))
            dataR[[i]] <- diag(x=temp$d[1:ncol(data[[1]])], nrow=ncol(data[[1]])) %*% t(temp$v[,1:ncol(data[[1]])])
            u[[i]] <- temp$u
         } else {
            # Set u[[i]] to identity matrix
            u[[i]] <- diag(1, nrow(data[[i]]))
            dataR[[i]] <- data[[i]]
         }
      }
   } else {
      dataR <- data
   }
   if (showProgress) { cat("Running JIVE algorithm for ranks:\njoint rank:", rankJ,", individual ranks:", rankA, "\n") }
   tempjive <- jive.iter(dataR, rankJ, rankA, conv=conv, maxiter=maxiter, orthIndiv=orthIndiv, showProgress=showProgress)
   Jperp=tempjive$joint
   Aperp=tempjive$individual
   if (est) {
      for (i in 1:length(data)) {
         Jperp[[i]] <- u[[i]] %*% Jperp[[i]]
         Aperp[[i]] <- u[[i]] %*% Aperp[[i]]
      }
   }
   current <- c(rankJ, rankA)
   nrun <- nrun + 1
  }
  converged <- ifelse(nrun == 10,F,T)
  return(list(data=data, joint=Jperp, individual=Aperp, rankJ=rankJ, rankA=rankA, converged=converged))
}


summary.jive <- function (object, ...) {
   method <- object$method
   # Display chosen ranks
   Rank <- c(object$rankJ, object$rankA)
   Source <- c("Joint", names(object$data))
   # Show variance explained by each component
   var.table <- c()
   for (i in 1:length(object$data)) {
      ssj <- norm(object$joint[[i]], type="f")^2
      ssi <- norm(object$individual[[i]], type="f")^2
      sse <- norm(object$data[[i]] - object$joint[[i]] - object$individual[[i]], type="f")^2
      sst <- norm(object$data[[i]], type="f")^2
      var.table <- cbind(var.table, round(c(ssj/sst, ssi/sst, sse/sst),3))
   }
   colnames(var.table) <- c(names(object$data))
   rownames(var.table) <- c("Joint", "Individual", "Residual")
   result <- list(method, cbind(Source, Rank), var.table)
   names(result) <- c("Method", "Ranks", "Variance")
   result
}


print.jive <- function (x, ...) {
   # Show variance explained by each component
   var.table <- c()
   for (i in 1:length(x$data)) {
      ssj <- norm(x$joint[[i]], type="f")^2
      ssi <- norm(x$individual[[i]], type="f")^2
      sse <- norm(x$data[[i]] - x$joint[[i]] - x$individual[[i]], type="f")^2
      sst <- norm(x$data[[i]], type="f")^2
      var.table <- cbind(var.table, round(c(ssj/sst, ssi/sst, sse/sst),3))
   }
   colnames(var.table) <- c(names(x$data))
   rownames(var.table) <- c("Joint", "Individual", "Residual")
   var.table
}


plot.jive <- function (x, type="var", ...) {
   if ("var" %in% type) {
      showVarExplained(x, ...)
   }
   if ("heat" %in% type) {
      showHeatmaps(x, ...)
   }
   if ("pca" %in% type) {
      showPCA(x,...)
   }
}

svdwrapper = function( x, nu, nv, verbose=F ){
  # workaround by Art Owen to avoid LAPACK errors
  # See https://stat.ethz.ch/pipermail/r-help/2007-October/143508.html  
  gotit = F
  try( {svdx = svd(x,nu,nv); gotit=T}, silent = !verbose )
  if( gotit )return(svdx)
  try( {svdtx = svd(t(x),nv,nu); gotit=T}, silent = !verbose )
  if( !gotit )stop("Error: svd(x) and svd(t(x)) both failed to converge.")
#  if( verbose )print("svd(x) failed but svd(t(x)) worked.")
  temp    = svdtx$u
  svdtx$u = svdtx$v
  svdtx$v = temp
  svdtx
}




