#' @title k prototypes clustering
#' @description Computes k prototypes clustering for mixed type data.
#' 
#' @details The algorithm like k means iteratively recomputes cluster prototypes and reassigns clusters.
#' Clusters are assigned using \eqn{d(x,y) =  d_{euclid}(x,y) + \lambda d_{simple\,matching}(x,y)}.
#' Cluster prototypes are computed as cluster means for numeric variables and modes for factors 
#' (cf. Huang, 1998).
#' 
#' 
#' @rdname kproto
#' @export 
kproto <- function (x, ...) 
  UseMethod("kproto")

#' @aliases kproto kproto.default 
#' 
#' @param x Data frame with both mumerics and factors.
#' @param k Either the number of clusters or a vector specifying indices of initial prototypes.
#' @param lambda Parameter > 0 to trade off between Euclidean distance of numeric variables 
#' and simple matching coefficient between categorical variables.
#' @param iter.max Maximum number of iterations if no convergence before.
#' @param nstart If > 1 repetetive computations with random initializations are computed and the result with minimum tot.dist is returned.
#' @param \dots Currently not used.
#
#' @return \code{\link{kmeans}} like object of class \code{kproto}:
#' @return \item{cluster}{Vector of cluster memberships.}
#' @return \item{centers}{Data frame of cluster prototypes.}
#' @return \item{lambda}{Distance parameter lambda.}
#' @return \item{size}{Vector of cluster sizes.}
#' @return \item{withinss}{Vector of summed distances to the cluster prototype per cluster.}
#' @return \item{tot.withinss}{Target function: sum of all distances to clsuter prototype.}
#' @return \item{dists}{Matrix with distances of observations to all cluster prototypes.}
#' @return \item{iter}{Prespecified maximum number of iterations.}
#' @return \item{trace}{List with two elements (vectors) tracing the iteration process: 
#' \code{tot.dists} and \code{moved} number of observations over all iterations.}
#'   
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 100
#' prb <- 0.9
#' muk <- 1.5 
#' clusid <- rep(1:4, each = n)
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' # apply k prototyps
#' kpres <- kproto(x, 4)
#' clprofiles(kpres, x)
#' 
#' # in real world  clusters are often not as clear cut
#' # by variation of lambda the emphasize is shifted towards factor / numeric variables    
#' kpres <- kproto(x, 2)
#' clprofiles(kpres, x)
#' 
#' kpres <- kproto(x, 2, lambda = 0.1)
#' clprofiles(kpres, x)
#' 
#' kpres <- kproto(x, 2, lambda = 25)
#' clprofiles(kpres, x)
#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' @references Z.Huang (1998): Extensions to the k-Means Algorithm for Clustering Large Data Sets with Categorical Variables, 
#' Data Mining and Knowledge Discovery 2, 283-304.
#' 
#' @rdname kproto
#' 
#' @method kproto default
#' @export 
kproto.default <- function(x, k, lambda = NULL,  iter.max = 100, nstart=1, ...){
  
  # initial error checks
  if(!is.data.frame(x)) stop("x should be a data frame!")
  if(ncol(x) < 2) stop("For clustering x should contain at least two variables!")
  if(iter.max < 1 | nstart < 1) stop("iter.max and nstart must not be specified < 1!")
  if(!is.null(lambda)) if(lambda < 0) stop("lambda must be specified >= 0 !")
  
  # check for numeric and factor variables
  numvars <- sapply(x, is.numeric)
  anynum <- any(numvars)
  catvars <- sapply(x, is.factor)
  anyfact <- any(catvars)
  if(!anynum) cat("\n No numeric variables in x! \n\n")
  if(!anyfact) cat("\n No factor variables in x! \n\n")
  
  # automatic calculation of lambda
  if(is.null(lambda)){
    if(anynum & anyfact){
      vnum <- mean(sapply(x[,numvars], var))
      vcat <- mean(sapply(x[,catvars], function(z) return(1-sum((table(z)/length(z))^2))))
      if (vnum == 0){
        warning("All numerical variables have zero variance.")
        anynum <- FALSE
      } 
      if (vcat == 0){
        warning("All categorical variables have zero variance.")
        anycat <- FALSE
      } 
      if(anynum & anyfact) {lambda <- vnum/vcat; cat("Estimated lambda:", lambda, "\n\n")}  
      else lambda <- 1   
    }
  } 
  
  # initialize prototypes
  if (length(k) == 1){
    ids <- sample(nrow(x), k)
    protos <- x[ids,]
  }
  if (length(k) > 1){
    ids <- k
    k <- length(ids)
    if(length(unique(ids)) != length(ids)) stop("If k is specified as a vector it should contain different indices!")
    if(any(ids<1)|any(ids>nrow(x))) stop("If k is specified as a vector all elements must be valid indices of x!")
    #check for integer
    protos <- x[ids,]
  } 
  
  # check for any equal prototypes and reduce cluster number in case of occurence
  keep.protos <- rep(TRUE,k)
  for(l in 1:(k-1)){
    for(m in (l+1):k){
      d1 <- sum((protos[l,numvars]-protos[m,numvars])^2) # euclidean for numerics
      d2 <- sum(protos[l,catvars] != protos[m,catvars]) # wtd simple matching for categorics 
      if((d1+d2) == 0) keep.protos[m] <- FALSE 
    }
  }
  if(!all(keep.protos)){
    protos <- protos[keep.protos,]
    k <- sum(keep.protos)
    cat("Equal prototyps merged. Cluster number reduced to:", k, "\n\n")      
  }
  
  rm(ids)
  
  # initialize clusters
  clusters  <- numeric(nrow(x)) 
  tot.dists <- NULL
  moved   <- NULL
  
  iter <- 1
  while(iter < iter.max){
    
    # compute distances 
    nrows <- nrow(x)
    dists <- matrix(NA, nrow=nrows, ncol = k)
    for(i in 1:k){
      #a0 <- proc.time()[3]      
      #d1 <- apply(x[,numvars],1, function(z) sum((z-protos[i,numvars])^2)) # euclidean for numerics
      d1 <- (x[,numvars] - matrix(rep(as.numeric(protos[i,numvars]), nrows), nrow=nrows, byrow=T))^2
      d1 <- rowSums(d1)
      #a1 <- proc.time()[3]      
      #d2 <- lambda * apply(x[,catvars],1, function(z) sum((z != protos[i,catvars]))) # wtd simple matching for categorics 
      d2 <- sapply(which(catvars), function(j) return(x[,j] != rep(protos[i,j], nrows)) )
      d2 <- lambda * rowSums(d2)
      #a2 <- proc.time()[3]      
      dists[,i] <- d1 + d2
      #cat(a1-a0, a2-a1, "\n")
    }
    
    # assign clusters 
    old.clusters  <- clusters
    # clusters      <- apply(dists, 1, function(z) which.min(z))
    clusters      <- apply(dists, 1, function(z) {a <- which.min(z); if (length(a)>1) a <- sample(a,1); return(a)}) # sample in case of multiple minima
    size          <- table(clusters)  
    min.dists     <- apply(cbind(clusters, dists), 1, function(z) z[z[1]+1])
    within        <- as.numeric(by(min.dists, clusters, sum))
    tot.within    <- sum(within)
    # prevent from empty classes
    #tot.within    <- numeric(k)
    #totw.list     <- by(min.dists, clusters, sum) 
    #tot.within[names(totw.list)] <- as.numeric(totw.list)
    
    # ...check for empty clusters and eventually reduce number of prototypes    
    if (length(size) < k){
      k <- length(size)
      protos <- protos[1:length(size),]  
      cat("Empty clusters occur. Cluster number reduced to:", k, "\n\n")
      
    }
    
    # trace
    tot.dists <- c(tot.dists, sum(tot.within))      
    moved <- c(moved, sum(clusters != old.clusters))
    
    # compute new prototypes
    for(i in 1:k){
      if(size[i] > 0){
        protos[i, numvars] <- sapply(x[clusters==i, numvars], mean)
        protos[i, catvars] <- sapply(x[clusters==i, catvars], function(z) levels(z)[which.max(table(z))])        
      }  
    }    
    
    # check for any equal prototypes and reduce cluster number in case of occurence
    keep.protos <- rep(TRUE,k)
    for(l in 1:(k-1)){
      for(m in (l+1):k){
        d1 <- sum((protos[l,numvars]-protos[m,numvars])^2) # euclidean for numerics
        d2 <- sum(protos[l,catvars] != protos[m,catvars]) # wtd simple matching for categorics 
        if((d1+d2) == 0) keep.protos[m] <- FALSE 
      }
    }
    if(!all(keep.protos)){
      protos <- protos[keep.protos,]
      k <- sum(keep.protos)
      cat("Equal prototyps merged. Cluster number reduced to:", k, "\n\n")      
    }
    
    # add stopping rules
    if(moved[length(moved)] ==  0) break
    if(k == 1) break
    
    #cat("iter", iter, "moved", moved[length(moved)], "tot.dists",tot.dists[length(tot.dists)],"\n" )      
    iter <- iter+1  
  }
  
  # create result: 
  res <- list(cluster = clusters,  
              centers = protos, 
              lambda = lambda, 
              size = size,
              withinss = within,
              tot.withinss = tot.within,   
              dists = dists, 
              iter = iter, 
              trace = list(tot.dists = tot.dists, moved = moved))
  
  # loop: if nstart > 1:
  if(nstart > 1)
    for(j in 2:nstart){
      res.new <- kproto(x=x, k=k, lambda = lambda,  iter.max = iter.max, nstart=1)
      if(res.new$tot.withinss < res$tot.withinss) res <- res.new
    }  
  
  class(res) <- "kproto"
  return(res)
}



#' @title k prototypes clustering
#'
#' @description Predicts k prototypes cluster memberships and distances for new data.
#' 
#' @details The algorithm like k means iteratively recomputes cluster prototypes and reassigns clusters.
#' Clusters are assigned using \eqn{d(x,y) =  d_{euclid}(x,y) + \lambda d_{simple\,matching}(x,y)}.
#' Cluster prototypes are computed as cluster means for numeric variables and modes for factors 
#' (cf. Huang, 1998).
#' 
#' @param object Object resulting from a call of resulting \code{kproto}.
#' @param newdata New data frame (of same structure) where cluster memberships are to be predicted.
#' @param \dots Currently not used.
#
#' @return \code{\link{kmeans}} like object of class \code{kproto}:
#' @return \item{cluster}{Vector of cluster memberships.}
#' @return \item{dists}{Matrix with distances of observations to all cluster prototypes.}
# 
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 100
#' prb <- 0.9
#' muk <- 1.5 
#' clusid <- rep(1:4, each = n)
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' # apply k prototyps
#' kpres <- kproto(x, 4)
#' predicted.clusters <- predict(kpres, x) 
#' 
#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' @rdname predict.kproto
#' @export
predict.kproto <-function(object, newdata, ...){
  lambda <- object$lambda
  k <- length(object$size)
  protos <- object$centers
  x <- newdata
  
  # check for numeric and factor variables  
  numvars <- sapply(x, is.numeric)
  catvars <- sapply(x, is.factor)
  
  if(!all(names(object$centers) == names(x))) warning("Variable names of newdata do not match!")
  if(ncol(object$centers) != ncol(x)) stop("Coloumn number of newdata does not match object!")
  if(!all(sapply(object$centers[,numvars], is.numeric))) stop("Numeric variables of object and newdata do not match!")
  if(!all(sapply(object$centers[,catvars], is.factor))) stop("Factor variables of object and newdata do not match!")
  
  # compute distances 
  dists <- matrix(NA, nrow=nrow(x), ncol = k)
  nrows <- nrow(x)
  for(i in 1:k){
    #a0 <- proc.time()[3]      
    #d1 <- apply(x[,numvars],1, function(z) sum((z-protos[i,numvars])^2)) # euclidean for numerics
    d1 <- (x[,numvars] - matrix(rep(as.numeric(protos[i,numvars]), nrows), nrow=nrows, byrow=T))^2
    d1 <- rowSums(d1)
    #a1 <- proc.time()[3]      
    #d2 <- lambda * apply(x[,catvars],1, function(z) sum((z != protos[i,catvars]))) # wtd simple matching for categorics 
    d2 <- sapply(which(catvars), function(j) return(x[,j] != rep(protos[i,j], nrows)) )
    d2 <- lambda * rowSums(d2)
    #a2 <- proc.time()[3]      
    dists[,i] <- d1 + d2
    #cat(a1-a0, a2-a1, "\n")
  }
  
  # assign clusters
  clusters      <- apply(dists, 1, function(z) {a <- which.min(z); if (length(a)>1) a <- sample(a,1); return(a)}) # sample in case of multiple minima
  
  res <- list(cluster = clusters, dists = dists)
  return(res)
}


#' @title profile k prototypes clustering
#'
#' @description Visualization of k prototypes clustering result for cluster interpretation.
#' 
#' @details For numerical variables boxplots and for factor variables barplots of each cluster are generated.
#' 
#' @param object Object resulting from a call of resulting \code{kproto}. Also other \code{kmeans} like objects with \code{object$cluster} and \code{object$size} are possible. 
#' @param x Original data.
#' @param vars Vector of either coloumn indices or variable names.
#
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 100
#' prb <- 0.9
#' muk <- 1.5 
#' clusid <- rep(1:4, each = n)
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' # apply k prototyps
#' kpres <- kproto(x, 4)
#' clprofiles(kpres, x)
#' 
#' # in real world  clusters are often not as clear cut
#' # by variation of lambda the emphasize is shifted towards factor / numeric variables    
#' kpres <- kproto(x, 2)
#' clprofiles(kpres, x)

#' 
#' kpres <- kproto(x, 2, lambda = 0.1)
#' clprofiles(kpres, x)
#' 
#' kpres <- kproto(x, 2, lambda = 25)
#' clprofiles(kpres, x)

#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' @rdname clprofiles
#' @importFrom graphics barplot
#' @importFrom graphics boxplot
#' @importFrom graphics legend
#' @importFrom graphics par
#' 
#' @export
clprofiles <- function(object, x, vars = NULL){
  if(length(object$cluster) != nrow(x)) stop("Size of x does not match cluster result!") 
  if(is.null(vars)) vars <- 1:ncol(x)
  if(!is.numeric(vars)) vars <- sapply(vars, function(z) return(which(colnames(x)==z)))
  if(length(vars) < 1) stop("Specified variable names do not match x!")
  
  #clusids <- as.numeric(names(object$size)) # object size not named for kmeans
  clusids <- sort(unique(object$cluster)) 
  
  par(ask=TRUE)
  for(i in vars){
    if(is.numeric(x[,i])){
      boxplot(x[,i]~object$cluster, col = clusids, main = colnames(x)[i])
      legend("topright", legend=clusids, fill = clusids)
    } 
    if(is.factor(x[,i])){
      tab <- table(x[,i], object$cluster)
      for(j in 1:length(object$size)) tab[,j] <- tab[,j]/object$size[j]
      barplot(t(tab), beside = TRUE, main = colnames(x)[i])
    } 
  } 
  invisible()
}


#' @title compares variance of all variables
#'
#' @description Investigation of variances to specify lambda for k prototypes clustering .
#' 
#' @details Variance of numeric variables and \eqn{1-\sum_i p_i^2} of categorical variance is computed.
#' 
#' @param x Original data.
#' 
#' @return \item{lambda}{Ratio of averages over all numeric/factor variables is returned.}
#' 
#' @examples
#' # generate toy data with factors and numerics
#' 
#' n   <- 100
#' prb <- 0.9
#' muk <- 1.5 
#' clusid <- rep(1:4, each = n)
#' 
#' x1 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x1 <- c(x1, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x1 <- as.factor(x1)
#' 
#' x2 <- sample(c("A","B"), 2*n, replace = TRUE, prob = c(prb, 1-prb))
#' x2 <- c(x2, sample(c("A","B"), 2*n, replace = TRUE, prob = c(1-prb, prb)))
#' x2 <- as.factor(x2)
#' 
#' x3 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' x4 <- c(rnorm(n, mean = -muk), rnorm(n, mean = muk), rnorm(n, mean = -muk), rnorm(n, mean = muk))
#' 
#' x <- data.frame(x1,x2,x3,x4)
#' 
#' a <- lambdaest(x)
#' res <- kproto(x, 2, lambda = a)
#' 
#' @author \email{gero.szepannek@@web.de}
#' 
#' @rdname lambdaest
#' 
#' @importFrom stats var
#' 
#' @export
lambdaest <- function(x){
  # initial error checks
  if(!is.data.frame(x)) stop("x should be a data frame!")
  
  # check for numeric and factor variables
  numvars <- sapply(x, is.numeric)
  anynum <- any(numvars)
  catvars <- sapply(x, is.factor)
  anyfact <- any(catvars)
  if(!anynum) cat("\n No numeric variables in x! \n\n")
  if(!anyfact) cat("\n No factor variables in x! \n\n")
  
  if(anynum) vnum <- sapply(x[,numvars], var)
  if(anyfact) vcat <- sapply(x[,catvars], function(z) return(1-sum((table(z)/length(z))^2)))
  if (mean(vnum) == 0){
    warning("All numerical variables have zero variance.\n
            No meaninful estimation for lambda.\n
            Rather use kmodes{klaR} instead of kprotos().")
    anynum <- FALSE
  } 
  if (mean(vcat) == 0){
    warning("All categorical variables have zero variance.\n
            No meaninful estimation for lambda!\n
            Rather use kmeans() instead of kprotos().")
    anyfact <- FALSE
  } 
  cat("Numeric variances:\n")
  print(vnum)
  cat("Average numeric variance:", mean(vnum), "\n\n")
  
  cat("Categorical variances:\n")
  print(vcat)
  cat("Average categorical variance:", mean(vcat), "\n\n")
  
  if(anynum & anyfact) {
    lambda <- mean(vnum)/mean(vcat); cat("Estimated lambda:", lambda, "\n\n")  
    return(lambda)
  }
  if(!(anynum & anyfact)) invisible()
}

#' @export
print.kproto <- function(x, ...){  
  cat("Numeric predictors:", sum(sapply(x$centers, is.numeric)), "\n")
  cat("Categorical predictors:", sum(sapply(x$centers, is.factor)), "\n")
  cat("Lambda:", x$lambda, "\n\n")
  
  cat("Number of Clusters:", length(x$size), "\n")
  cat("Cluster sizes:", x$size, "\n")
  cat("Within cluster error:", x$withinss, "\n\n")
  
  cat("Cluster prototypes:\n")
  print(x$centers)
}
