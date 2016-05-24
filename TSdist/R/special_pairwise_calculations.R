# # # ADDITIONAL FUNCTIONS FOR TRAIN/TEST DISTANCE MATRIX CALCULATION# # # 

# These functions are small modifications of some general functions of the
# package TSclust. The modifications are done to enable pairwise
# dissimilarity calculations between series of two different databases.

PairwiseDistances1 <- function(X, distfun, ...) {
  if (is.matrix(X)) {
  n <- dim(X)[1]
  }
  if (is.list(X)) {
  n <- length(X)
  }
  distances <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
        d <-   distfun(X, Y=NULL, i, j, ...)
        distances[i,j] <- d
        distances[j,i] <- d
    }
}
return(distances)
}

PairwiseDistances2 <- function(X, Y, distfun, ...) {
  if (is.matrix(X)) {
  n <- dim(X)[1]
  m <- dim(Y)[1]
  }
  if (is.list(X)) {
  n <- length(X)
  m <- length(Y)
  }
  distances <- matrix(0, n, m)
  for (i in 1:n) {
    for (j in 1:m) {
        d <-   distfun(X, Y, i, j, ...)
        distances[i,j] <- d
    }
  }
  return(distances)
}

# Function that converts matrix to list
MatrixToList <- function(X) {
  aux <- X
  X <- list()
  for (i in 1:nrow(aux)) {
    X[[i]] <- aux[i,]
  }
  names(X) <- rownames(aux)
  return(X)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Here we define the pairwise comparison of some special distances for 
# one or two databases.

# # AR.PIC distance between i-th and j-th series of X (or X and Y).
PairwiseARPicDistance <- function(X, Y=NULL, i, j, order.x=NULL, order.y=NULL, 
                                    permissive=TRUE) {

  
  # If data is given as a matrix convert into list
  if (! is.list(X)) {X <- MatrixToList(X)}
  
  # Show initial errors, which will not allow distance matrix calculation.
  options(show.error.messages = TRUE)
  if (! is.null(order.x) && dim(order.x)[1] != length(X)) {
  stop("The length of order.x must be equal to the number of 
       series in the database.")}
  
  # Case with only one database
  if (is.null(Y)) {

  # Do not show errors arising in pairwise calculations. Simply return NA.
    options(show.error.messages = FALSE)
    d <- diss.AR.PIC(X[[i]], X[[j]], order.x[i, ], order.x[j, ], permissive)
  
  # Case of two databases (for TRAIN/TEST)
  } else {
    
    # If data is given as a matrix convert into list
    if (! is.list(Y)) {Y <- MatrixToList(Y)}
  
    if (! is.null(order.y) && dim(order.y)[1] != length(Y)) {
      stop("The length of order.y must be equal to the number of 
       series in the database Y.")}
    

      # Do not show errors arising in pairwise calculations. Simply return NA.
      options(show.error.messages = FALSE)
      d <- diss.AR.PIC(X[[i]], Y[[j]], order.x[i, ], order.y[j, ], permissive)
      # Calculate pairwise distances between series in database X and 
      # series in database Y.
    }
  options(show.error.messages = TRUE)
  return(d) 
}

# AR.LPC.CEPS dissimilarity measure between i-th and j-th series of 
# X (or X and Y).

PairwiseARLPCCepsDistance <- function(X, Y=NULL, i, j, k=50, order.x=NULL, 
                                         order.y=NULL, seasonal.x=NULL, 
                                         seasonal.y=NULL, permissive=TRUE) {
  
  # If data is given as a matrix convert into list
  if (! is.list(X)) {X <- MatrixToList(X)}
  
  # Show initial errors, which will not allow distance matrix calculation.
  options(show.error.messages = TRUE)
  
  if (! is.null(order.x) && dim(order.x)[1] != length(X)) 
    stop("The number of rows of order.x must be equal to the number 
         of series in the database X.")
  
  # Case with only one database
  if (is.null(Y)) {
    

    
    # Calculate pairwise distances for database X.
        if (is.null(seasonal.x)) {
          seasonal.x[[i]] <-  list(order=c(0, 0, 0), period=NA)
          seasonal.x[[j]] <-  list(order=c(0, 0, 0), period=NA)
        }
        # Do not show errors arising in pairwise calculations. Simply return NA.
        options(show.error.messages = FALSE)
        d <-   diss.AR.LPC.CEPS(X[[i]], X[[j]], k, order.x[i, ], order.x[j, ], 
                                seasonal.x[[i]], seasonal.x[[j]], permissive)      
    # For TRAIN/TEST case
    } else {
      
      # If data is given as a matrix convert into list
      if (! is.list(Y)) {Y <- MatrixToList(Y)}
      
      if (! is.null(order.y) && dim(order.y)[1] != length(Y)) {
        stop("The length of order.y must be equal to the number of 
       series in the database Y.")
      }
      if (is.null(seasonal.x)) {
      seasonal.x[[i]] <-  list(order=c(0, 0, 0), period=NA)
      }
      if (is.null(seasonal.y)) {
      seasonal.y[[j]] <-  list(order=c(0, 0, 0), period=NA)
      }
     
      
      # Do not show errors arising in pairwise calculations. Simply return NA.
      options(show.error.messages = FALSE)
      d <-   diss.AR.LPC.CEPS(X[[i]], Y[[i]], k, order.x[i, ], order.y[j, ], seasonal.x[[i]], seasonal.y[[j]],permissive) 
      }
  options(show.error.messages = TRUE)
  return(d) 
}

# PRED dissimilarity measure for datasets.

PairwisePredDistance <- function(X, Y=NULL, h, B=500, logarithms.x=NULL, logarithms.y=NULL, differences.x=NULL, differences.y=NULL, plot=FALSE) {
 
  # If data is given as a matrix convert into list
  if (! is.list(X)) {X <- MatrixToList(X)}
  n1 <- length(X)
  
  if (! is.null(logarithms.x) && length(logarithms.x) != n1) {
    stop("The length of logarithms.x must be equal to the number of series in X.")}
  
  if (! is.null(differences.x) && length(differences.x) != n1) {
    stop("The length of differences.x must be equal to the number of series in X.")}
  
  if (is.null(logarithms.x)) {
    logarithms.x <- rep(FALSE, n1)
  }
  if (is.null(differences.x)) {
    differences.x <- rep(0, n1)
  }
  
  # Case of only one database
  if (is.null(Y)) {
    
  # Calculate all the individual densities by applying diss.PRED to database X.
    individual.dens1 <- list()
    ii <- 1
    while (ii < n1) {
      dP <- diss.PRED(X[[ii]], X[[ii + 1]], h , B, logarithms.x[ii], logarithms.x[ii + 1], differences.x[ii], differences.x[ii+1], FALSE )
      individual.dens1[[ii]] <- list(dens=dP$dens.x, bw=dP$bw.x)
      individual.dens1[[ii + 1]] <- list(dens=dP$dens.y, bw=dP$bw.y)
      ii = ii + 2
      if (ii == n1) {ii <- ii -1}
    }
    
  densities <- list()
  distances <- matrix(0, n1, n1)
  rownames(distances) <- names(X)
  
  for (i in 1:(n1 - 1) ) {
    for (j in (i + 1):n1 ) {
      # Calculate integrated L1 distance between the forecasts.
      distance <- L1dist(individual.dens1[[i]]$dens, individual.dens1[[j]]$dens, individual.dens1[[i]]$bw, individual.dens1[[j]]$bw )
      distances[i, j] <- distance
      distances[j, i] <- distance
    }
  } 
  
  # Case of two databases (X and Y)
  } else {
    
    # If data is given as a matrix convert into list
    if (! is.list(Y)) {Y <- MatrixToList(Y)}
    n2 <- length(Y)
    
    if (! is.null(logarithms.y) && length(logarithms.y) != n2) {
      stop("The length of logarithms.y must be equal to the number of series in Y.")}
    
    if (! is.null(differences.y) && length(differences.y) != n2) {
      stop("The length of differences.y must be equal to the number of series in Y.")}
    
    if (is.null(logarithms.y)) {
      logarithms.y <- rep(FALSE, n2)
    }
    if (is.null(differences.y)) {
      differences.y <- rep(0, n2)
    }
    densities <- list()
    
    # Calculate all the individual densities by applying diss.PRED to database X and Y, together.
    individual.dens1 <- list()
    individual.dens2 <- list()
    ii <- 1
    XY <- append(X, Y)
    logarithms.xy <- append(logarithms.x, logarithms.y)
    differences.xy <- append(differences.x, differences.y)
    
    while (ii < n1+n2) {
      dP <- diss.PRED(XY[[ii]], XY[[ii + 1]], h , B, logarithms.xy[ii], logarithms.xy[ii+1] , differences.xy[ii], differences.xy[ii + 1], FALSE )
      individual.dens1[[ii]] <- list(dens=dP$dens.x, bw=dP$bw.x)
      individual.dens1[[ii + 1]] <- list(dens=dP$dens.y, bw=dP$bw.y)
      ii <- ii + 2
      if (ii == n1+n2) {ii <- ii - 1}
    }
    
    densities <- list()
    distances <- matrix(0, n1, n2)
    rownames(distances) <- names(X)
    colnames(distances) <- names(Y)
    
    for (i in 1:n1 ) {
      for (j in 1:n2 ) {
        # Calculate integrated L1 distance between the forecasts.
        distance <- L1dist(individual.dens1[[i]]$dens, individual.dens1[[n1+j]]$dens, individual.dens1[[i]]$bw, individual.dens1[[n1+j]]$bw )
        distances[i, j] <- distance
      }
    } 
    
  }
  
  return(distances)

}


# # SPEC.LLR dissimilarity measure for dataset X (or X and Y).

PairwiseSpecLLRDistance <- function(X, Y=NULL, method="DLS", 
                              alpha=0.5, plot=FALSE, n=NULL) {
  
  # If data is given as matrix, convert to list
  if (! is.list(X)) {X <- MatrixToList(X)}
  if (! is.null(Y) && ! is.list(Y)) {Y <- MatrixToList(Y)}
  
  if (is.null(n)) {
    n <- length(X[[1]])
  }
  
  
  if (n > 0) {
    interpfun <- NULL
    type <-  (pmatch(method, c("DLS", "LK" )))
    if (is.na(type)) {
      stop(paste("Unknown method", method))
    } else if (type == 1) {
      interpfun <- interp.SPEC.LS
    }
    else if (type == 2) {
      interpfun <- interp.W.LK
    }
    d <- PairwiseInterpSpecDistance(X, Y, n, interpfun, integrate.divergenceW, alpha)
    
    } else {
    if (is.null(Y)) {
    d <- dist(X, method="TSDistances", distance="spec.llr", alpha=alpha,  plot=FALSE, n=n)
    } else {
    d <- dist(X, Y, method="TSDistances", distance="spec.llr", alpha=alpha,  plot=FALSE, n=n)
    }
  }
  options(show.error.messages = TRUE)
  return(d)
}


# # SPEC.GLK dissimilarity measure for dataset X (or datasets X and Y).

PairwiseSpecGLKDistance <- function(X, Y=NULL, plot=FALSE) {
  
  # If data is given as matrix, convert to list
  if (! is.list(X)) {X <- MatrixToList(X)}
  if (! is.null(Y) && ! is.list(Y)) {Y <- MatrixToList(Y)}
  
  PairwiseInterpSpecDistance(X, Y, floor(length(X[[1]])/2), interp.SPEC.GLK, integrate.GLK)
}

# # SPEC.ISD dissimilarity measure for datasets.

PairwiseSpecISDDistance<- function(X, Y=NULL, plot=FALSE, n=NULL) {

  # If data is given as matrix, convert to list
  if (! is.list(X)) {X <- MatrixToList(X)}
  if (! is.null(Y) && ! is.list(Y)) {Y <- MatrixToList(Y)}
  
  
  if (is.null(n)) {
    n <- length(X[[1]])
  }
  
  
  if (n > 0) {
  d <- PairwiseInterpSpecDistance(X, Y, n, interp.SPEC.LOGLIKELIHOOD , integrate.ISD)
  }else {
    if (is.null(Y)) {
    d <- dist(X, method="TSDistances", distance= "spec.isd", plot=FALSE, n=n)
    } else {
    d <- dist(X, Y, method="TSDistances", distance= "spec.isd", plot=FALSE, n=n) 
    }
  }
  return(d)
}



# Adaptation of function multidiss.interp.SPEC of the TSclust package
# to possibilite the calculation of distances between two different 
# databases.
PairwiseInterpSpecDistance <- function(X, Y=NULL, n, interpfun, 
                                         integrationfun, ...) {
  l1 <- length(X)
  # get the interpolated values
  interps1 <- lapply(X, interpfun, n)
  base <- interps1[[1]]$x
  # For only one database
  if (is.null(Y)) {
  dists <- matrix(0, l1, l1) 
  for (i in 1:(l1-1)) {
    for (j in (i+1):l1) {
      d <- integrationfun(base, interps1[[i]]$y , interps1[[j]]$y, ...)
      dists[i,j] <- d
      dists[j,i] <- d
    }
  }
  # For two databases, TRAIN/TEST framework
  } else {
  l2 <- length(Y)
  dists <- matrix(0, l1, l2) 
  interps2 <- lapply(Y, interpfun, n)
  # # calc the function with the interpolated values
  for (i in 1:l1) {
    for (j in 1:l2) {
      d <- integrationfun(base, interps1[[i]]$y , interps2[[j]]$y, ...)
      dists[i,j] <- d
    }
  }
  }
  dists
}


# # # CDM distance (does not give the same result if we use dist.)

PairwiseCDMDistance <- function(X, Y=NULL, i, j, ...) {
  
  # If data is given as matrix, convert to list
  if (! is.list(X)) {X <- MatrixToList(X)}
  if (! is.null(Y) && ! is.list(Y)) {Y <- MatrixToList(Y)}
  
  if (is.null(Y)) {
    as.numeric(diss.CDM(X[[i]], X[[j]], ...))
  } else {
    as.numeric(diss.CDM(X[[i]], Y[[j]], ...))
  }
}

# # # NCD distance (does not give the same result if we use dist.)

PairwiseNCDDistance <- function(X, Y=NULL, i, j, ...) {
  
  # If data is given as matrix, convert to list
  if (! is.list(X)) {X <- MatrixToList(X)}
  if (! is.null(Y) && ! is.list(Y)) {Y <- MatrixToList(Y)}
  
  if (is.null(Y)) {
    as.numeric(diss.NCD(X[[i]], X[[j]], ...))
  } else {
    as.numeric(diss.NCD(X[[i]], Y[[j]], ...))
  }
}

# Special function to calculate Frechet distance for databases. 
# This function does not work with dist.

PairwiseFrechetDistance <- function(X, Y=NULL, i, j, ...) {
  
  # If data is given as matrix, convert to list
  if (! is.list(X)) {
    X <- MatrixToList(X)
  }
  if (! is.null(Y) && ! is.list(Y)) {
    Y <- MatrixToList(Y)
  }
  if (is.null(Y)) {
    as.numeric(FrechetDistance(X[[i]], X[[j]], ...))
  } else {
    as.numeric(FrechetDistance(X[[i]], Y[[j]], ...))
  }
}


###ADDITIONAL FUNCTIONS FOR TRAIN/TEST DISTANCE MATRIX CALCULATION for
## PDC distance

divergenceMatrix2 <- function(codebooks1, codebooks2, divergence)
{
  l1 <- dim(codebooks1)[1]
  l2 <- dim(codebooks2)[1]
  mt <- matrix(rep(0, l1 * l2), l1, )
  for (i in 1:l1)
  {
    for (j in 1:l2)
    {
      mt[i, j] <- divergence(codebooks1[i, ], codebooks2[j, ])
    }
  }
  return(mt);
}


PDCDist2 <- function(X, Y, m=NULL, t=NULL, divergence=symmetricAlphaDivergence)
{
  if (is.null(t) | is.null(m)) {
    ent <- entropyHeuristic(cbind(X,Y))
    
    if (is.null(m)) {
      m <- ent$m;
    }
    
    if (is.null(t)) {
      t <- ent$t;
    }
  }
  
  codebooks1 <- ConvertMatrix(X,m,t);
  codebooks2 <- ConvertMatrix(Y,m,t);
  D <- divergenceMatrix2( codebooks1,  codebooks2, divergence );
  
  
  
  return(D);
}

#Internal function used by pdc to calculate the permutation distributions from 
#the series
ConvertMatrix <-function(X, m, td){
  return( t(apply(X, 2, codebook, m=m, t=td)) )
}