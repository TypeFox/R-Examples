
# examples of callback functions

# cgraph: a callback function to check whether a graph is connected, and
# optionally, to find minimum distances; assumes m is 1s and 0s, showing
# adjacency; see http://heather.cs.ucdavis.edu/parprocbook

# arguments:
#    ev:  the ev in matpow()
#    cbinit:  when TRUE, initialization will be done
#    mindist:  if TRUE, the matrix of min distances will be computed
# return values (connected, dists) placed in ev

# currently set up for ordinary R matrices or gputools

# it is recommended that matpow() be called with squaring = TRUE

cgraph <- function(ev,cbinit=FALSE,mindist=FALSE) {
   if (cbinit) {
      diag(ev$m[,]) <- 1
      diag(ev$prod1[,]) <- 1
      ev$k <- nrow(ev$m) - 1
      ev$dists <- ev$m
      ev$connected <- FALSE
      if (mindist && ev$squaring) 
         stop("squaring cannot be used with the mindist option")
      return()
   }
   # requires the matrix class to have an all() method
   if (ev$i %% 2 == 1) prd <- ev$prod2 else prd <- ev$prod1
   if (all(prd > 0)) {
      ev$connected <- TRUE
      ev$stop <- TRUE
   }
   if (mindist) {
      tmp <- prd[,] > 0
      ev$dists[tmp & ev$dists == 0] <- ev$i + 1
   }
}

# eig:  callback to find the principal eigenvalue of a matrix, using the
# power method

# arguments:
#    ev:  the ev in matpow()
#    cbinit:  when TRUE, initialization will be done
#    x:  optional initial guess for the principal eigenvector
#    eps:  epsilon, to check convergence
# return values placed in copy of ev; look at ev$x for the eigenvector,
# ev$i for the number of iterations, and ev$stop to see if
# convergence was reached

# currently set up for ordinary R matrices or gputools

eig <- function(ev,cbinit=FALSE,x=NULL,eps=1e-08) {
   if (cbinit) {
      if (is.null(x)) 
         x <- rep(1,nrow(ev$m))
      ev$x <- x
      ev$oldx <- x
      ev$eps <- eps
      return()
   }
   ev$x <- ev$x / normvec(ev$x)
   cmd <- ev$genmulcmd("ev$m","ev$x","ev$x")
   eval(parse(text=cmd))
   diff <- normvec(ev$x - ev$oldx)
   if (diff / sum(abs(ev$x)) < ev$eps) ev$stop <- TRUE
   ev$oldx <- ev$x
}

# mc:  callback to find the long-run distribution of an aperiodic,
# discrete-time Markov chain

mc <- function(ev,cbinit=FALSE,eps=1e-08) {
   if (cbinit) {
      return()
   }
   diff <- norm(ev$prod1 - ev$prod2)
   if (ev$i %% 2 == 1) prd <- ev$prod2 else prd <- ev$prod1
   if (diff / norm(prd) < eps) {
      ev$stop <- TRUE
      # long-run distribution vector
      ev$pivec <- colMeans(prd)
   }
}

# mexp:  compute the exponential of a matrix, used in solutions to
# differential equations, control theory and continuous-time Markov
# chains

mexp <- function(ev,cbinit=FALSE,eps=1e-08) {
   if (cbinit) {
      if (ev$squaring) stop("squaring not allowed with mexp")
      # ev$esum will be the current sum of the series
      ev$esum <- diag(nrow(ev$m)) + ev$m
      ev$esumold <- ev$esum
      return()
   }
   if (ev$i %% 2 == 1) prd <- ev$prod2 else prd <- ev$prod1
   ev$esum <- ev$esum + (1/factorial(ev$i+1)) * prd
   diff <- norm(ev$esum - ev$esumold)
   if (diff / norm(ev$esum) < eps) {
      ev$stop <- TRUE
   }
   ev$esumold <- ev$esum
}
