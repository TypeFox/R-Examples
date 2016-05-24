#' Perform imputation of a data frame using k-NN.
#' 
#' Perform imputation of missing data in a data frame using the k-Nearest Neighbour algorithm.
#' For discrete variables we use the mode, for continuous variables the median value is instead taken.
#' 
#' @name knn.impute
#' @rdname knn.impute
#' 
#' @param data a data frame
#' @param k number of neighbours to be used; for categorical variables 
#'        the mode of the neighbours is used, for continuous variables 
#'        the median value is used instead. Default: 10.
#' @param cat.var vector containing the indices of the variables to be 
#'        considered as categorical. Default: all variables.
#' @param to.impute vector indicating which rows of the dataset are to be imputed. 
#'        Default: impute all rows.
#' @param using vector indicating which rows of the dataset are to be used 
#'        to search for neighbours. Default: use all rows.
#'      
#' @return imputed data frame.
#' 
#' @export knn.impute
knn.impute <- function( data, k = 10, cat.var = 1:ncol(data), 
	to.impute = 1:nrow(data), using = 1:nrow(data) )
{
  n.var <- dim(data)[2]
  num.var <- setdiff(1:n.var,cat.var)
  
  use.data <- data[using,,drop=FALSE] # retain dimensions even for one row
  imp.data <- data[to.impute,,drop=FALSE] 
  
  storage.mode(use.data) <- "double"
  storage.mode(imp.data) <- "double"
  
  use.num.var.max <- apply(use.data[,num.var],2,max,na.rm=TRUE)
  use.num.var.min <- apply(use.data[,num.var],2,min,na.rm=TRUE)
  
  na.cases <- (rowSums(is.na(imp.data)) > 0)
   
  neigh <- rep(0,k)
  storage.mode(num.var) <- "integer"
  
  for( i in 1:nrow(imp.data) )
  {
    if( na.cases[i] )
    {
      num.var.max <- pmax(use.num.var.max, imp.data[i,num.var], na.rm=TRUE)
      num.var.min <- pmin(use.num.var.min, imp.data[i,num.var], na.rm=TRUE)
      num.var.range <- num.var.max - num.var.min
      
      d <- .Call( "heom_dist", imp.data[i,], use.data, 
                  num.var, num.var.range, PACKAGE = "bnstruct" )
      s <- sort(d, index.return=TRUE)$ix
      
      for( j in which(is.na(imp.data[i,])) )
      {
        # find the k closest neighbours with nonmissing value for j
        ind.neigh <- 1
        # if to.impute[i] is in _using_, it is the first of the list
        if (to.impute[i] %in% using)
          ind.s <- 2
        else
          ind.s <- 1
        while( ind.neigh < k + 1)
        {
          if( !is.na(use.data[s[ind.s],j]) )
          {
            neigh[ind.neigh] <- use.data[s[ind.s],j]
            ind.neigh <- ind.neigh + 1
          }
          ind.s <- ind.s + 1
        }
        
        # impute from the neighbours
        if( j %in% cat.var )
          imp.data[i,j] <- stat.mode.ord(neigh)
        else
          imp.data[i,j] <- median(neigh)      
      }
    }
  }
  
  return(imp.data)
}

# Heterogeneous Euclidean Overlap Metric between 
# vector x and the rows of matrix Xr
# heom.dist <- function( x, Xr, num.var, num.var.range ) 
# {
#   dXr <- dim(Xr)
#   rx <- matrix(rep(x,dXr[1]),dXr[1],dXr[2],byrow=TRUE)
#   cat.var <- setdiff(1:dXr[2],num.var)
#   
#   # compute distance matrix (numeric variables)
#   num.dist <- (rx[,num.var] - Xr[,num.var]) / num.var.range
#   num.dist[is.na(num.dist)] <- 1
#   
#   # compute distance matrix (categorical variables)
#   cat.dist <- rx[,cat.var] != Xr[,cat.var]
#   cat.dist[is.na(cat.dist)] <- 1
#   
#   dist = 0
#   if ( length(cat.var) < 2 )
#     dist <- dist + sum(cat.dist^2)
#   else
#     dist <- dist + rowSums(cat.dist^2)
#   if( length(num.var) < 2 )
#     dist <- dist + sum(num.dist^2)
#   else
#     dist <- dist + rowSums(num.dist^2)
#   dist <- sqrt(dist)
#   
#   return(dist)
# }

# statistical mode, for ordered categorical vectors. 
# returns the first mode in case of multiple modes
stat.mode.ord <- function(x) 
{
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# statistical mode, for categorical vectors. 
# returns a vector in case of multiple modes
stat.mode <- function(x) 
{
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[which(tab == max(tab))]
}


#' tune the parameter k of the knn algorithm used in imputation.
#' 
#' @param data a data frame
#' @param cat.var vector containing the categorical variables
#' @param k.min minimum value for k
#' @param k.max maximum value for k
#' @param frac.miss fraction of missing values to add
#' @param n.iter number of iterations for each k
#' @param seed random seed
#' 
#' @return matrix of error distributions
#' 
#' @export tune.knn.impute
tune.knn.impute <- function( data, cat.var = 1:ncol(data), k.min = 1, k.max = 20, frac.miss = 0.1, n.iter = 20, seed = 0 )
{
  n.var <- dim(data)[2]
  num.var <- setdiff(1:n.var,cat.var)
  num.var.range <- diff(apply(data[,num.var],2,range,na.rm = T,finite=T))
  
  storage.mode(data) <- "double"
  storage.mode(num.var) <- "integer"
  
  # results: error distribution for each k
  err.dist <- matrix(0,n.iter,k.max-k.min+1)
  err.dist.ref <- rep(0,n.iter)
  
  set.seed(seed)
  
  # probability of a new missing value for each variable
  miss.p.var <- frac.miss * colSums(is.na(data)) / (nrow(data) - colSums(is.na(data)))  
  for ( it in 1:n.iter )
  {
    cat("iter:",it,"\n")
    err.all <- matrix(0,0,k.max-k.min+1)
    err.ref = c()
    miss.data <- data
    # add missingness
    for ( var in 1:ncol(data) )
    {
      nomiss <- which(!is.na(data[,var]))
      miss.data[ nomiss[ runif(length(nomiss)) < miss.p.var[var] ], var] <- NA
    }
    
    new.na.cases <- (rowSums(is.na(miss.data) &! is.na(data)) > 0)
    neigh <- rep(0,k.max)
    
    # compute neighbours and impute
    for( i in 1:nrow(miss.data) )
    {
      if( new.na.cases[i] )
      {
        d <- .Call( "heom_dist", miss.data[i,], miss.data, 
                    num.var, num.var.range, PACKAGE = "bnstruct" )
        s <- sort(d, index.return=TRUE)$ix
        
        for( j in which(is.na(miss.data[i,]) &!is.na(data[i,])) )
        {
          # find the k.max closest neighbours with nonmissing value for j
          ind.neigh <- 1
          ind.s <- 2
          while( ind.neigh < k.max + 1)
          {
            if( !is.na(miss.data[s[ind.s],j]) )
            {
              neigh[ind.neigh] <- miss.data[s[ind.s],j]
              ind.neigh <- ind.neigh + 1
            }
            ind.s <- ind.s + 1
          }
          
          # compute imputation error with k.min to k.max neighbours,
          # with HEOM distance
          err.k <- rep(0,k.max-k.min+1)
          for( k in k.min:k.max )
          {
            if( j %in% cat.var )
              err.k[k-k.min+1] <- (stat.mode.ord(neigh[1:k]) != data[i,j])
            else
              err.k[k-k.min+1] <- abs(median(neigh[1:k]) - data[i,j]) /
              num.var.range[num.var==j]
          }
          if( j %in% cat.var )
            err.ref <- c(err.ref,stat.mode.ord(miss.data[
              !is.na(miss.data[,j]),j]) != data[i,j])
          else err.ref <- c(err.ref,
                            abs(median(miss.data[,j],na.rm=T) - data[i,j]) / 
                              num.var.range[num.var==j])
          # add to results
          err.all <- rbind(err.all, err.k)
        }
      }
    }
    err.dist[it,] <- colMeans(err.all)
  }
  return(err.dist)
}
