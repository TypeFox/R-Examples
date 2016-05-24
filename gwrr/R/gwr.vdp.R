gwr.vdp <-
function(form, locs, data, phi, kernel="exp", sel.ci=30, sel.vdp=0.5){
   # Parse variables in formula to pass to function
   rhs <- as.character(form)[3]
   rhs.v <- strsplit(rhs, " + ", fixed=TRUE)   # Returns a list with 1 first element, unknown 2 elements
   n.l <- length(rhs.v[[1]])   # get number of x variables

   # Create design matrix
   db <- data
   N <- dim(db)[1]
   X <- rep(1,N)   # Assume intercept for now
   for(i in 1:n.l) X <- cbind(X, db[,rhs.v[[1]][i]])
   
   # Calculate pairwise distances
   library(fields)   
   S <- rdist(locs)   # Assume Euclidean distance is appropriate for now

   if (kernel == "exp") W <- w.exp(phi, S)
   if (kernel == "gauss") W <- w.gauss(phi,S)
   
   N <- dim(X)[1]
   p <- dim(X)[2]
   vdp.k <- array(0, dim=c(N,p))      # condition index
   vdp.pi <- array(0, dim=c(N,p,p))   # VDPs

   for (i in 1:N){
      W.i <- diag(W[i,])
      W.i.h <- W.i^.5
      W.X <- W.i.h %*% X
      vdp <- local.vdp(W.X, N)
      vdp.k[i,] <- vdp$k.X
      vdp.pi[i,,] <- vdp$pi.ij
   }

   # Return largest column of condition index and associated VDPs for now
   vdp.k[,p]
   vdp.pi[,p,]
   #params <- list(vdp.k[,p], vdp.pi[,p,])
   #names(params) <- c("condition", "vdp")
   #params
   
   # Flag shared components of large size or large condition index.
   TF <- vdp.pi >= sel.vdp
   TF.sum <- apply(TF,c(1,2),sum)   # row sums are now in a row vector
   TF.sum.TF <- TF.sum > 1
   flag.k <- vdp.k[,p] >= sel.ci
   flag.pi <- TF.sum.TF[,p]
   flag.k.pi <- flag.k == TRUE & flag.pi == TRUE
   params <- list(vdp.k[,p], vdp.pi[,p,], flag.k, flag.pi, flag.k.pi)
   names(params) <- c("condition", "vdp", "flag.cond", "flag.vdp", "flag.cond.vdp")
   params
}

