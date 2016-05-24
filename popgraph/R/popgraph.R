#' Default constructor for a \code{PopulationGraph} object
#' 
#' This function is the default constructor for \code{PopulationGraph} objects.  Mechanistically,
#'  a \code{PopulationGraph} is a class that has some meta data and an \code{igraph} object 
#'  within it.
#' @param x An object of type \code{matrix} that holds the data to be analyzed.
#' @param groups A factor indicating population membership of each row in \code{x}
#' @param alpha The significance level to test edge retention (default = 0.05).
#' @param tol A measure of tolerance for the retention of multivariate data columns (default sdev=1e-4)
#' @return An object of type \code{popgraph}.
#' @export
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
popgraph <- function( x, groups, alpha=0.05, tol=1.0e-4 ) {
  
  if( missing(x) )
    stop("You must use a matrix to pass data to this function, object to create a 'PopulationGraph'" )
  if( is(x,"data.frame") ) 
    stop("The data passed to popgraph() needs to be a numeric matrix. If you are using gstudio, convert your data first using to_mv().")
  if( missing( groups) )
    stop("You need to specify which 'groups' the nodes will represent.")

  # make sure they are only factors with samples
  groups <- factor(as.character(groups))
    
  
  # sort genos in order of Population
  #x <- x[ with(x,order(x[[groups]])), ]
  
  # throw warning of some groups are small
  t <- table(groups)
  if( any( t < 4 ) ) {
    popnames <- paste( names(which(t<4)), collapse=", ")
    warning( paste( "You have strata (",popnames,") that have fewer than 4 individuals. This anlaysis needs to have a good estimate of within stratum variance."))
  }

  # make the graph
  N <- length(groups)
  K <- length(levels(groups))
  
  critVal <- qchisq( 1.0-alpha,1 )
  EdgeStr <- matrix(0,K,K)
  
  # rotate all data and keep stuff that is not invariant.
  pcfit <- prcomp(x,retx=T)
  mv <-pcfit$x[,pcfit$sdev > tol] 
  
  P <- ncol(mv)
  Pop.priors <- as.numeric(table(groups)/N)
  K <- length( Pop.priors)
  pop.means <- tapply( mv, list( rep(groups,P),col(mv)), mean)
  sigma.w <- sqrt( diag( var( mv - pop.means[groups,])))
  
  if( any( sigma.w < tol ) )
    warning( paste( "Dropped rotated genetic variable '",
                    paste((as.character(1:P)[sigma.w<tol]),collapse=","),
                    "' from the analysis due to constancy within groups",
                    sep="") )  
  
  scaling <- diag(1/sigma.w,,P)
  fac <- 1/(N-K)
  X <- sqrt(fac) * (mv - pop.means[groups,]) %*% scaling
  X.s <- svd(X, nu = 0)
  rank <- sum(X.s$d > tol)
  if(rank < P) 
    warning( paste( (P-rank),
                    " variables are collinear and being dropped from the discriminant rotation.",
                    sep=""))
  scaling <- scaling %*% X.s$v[, 1:rank] %*% diag(1/X.s$d[1:rank],,rank)
  mu <- colSums(Pop.priors %*% pop.means)
  X <- sqrt((N * Pop.priors)/(K-1)) * scale(pop.means, center=mu, scale=FALSE) %*% scaling
  X.s <- svd(X, nu = 0)
  rank <- sum(X.s$d > tol * X.s$d[1L])
  scaling <- scaling %*% X.s$v[, 1L:rank]
  
  means <- colMeans(pop.means)
  LDValues <- scale( mv, center=means, scale=FALSE ) %*% scaling
  
  allLD <- centroid_distance( LDValues, groups )
  allSD <- centroid_variance( LDValues, groups )
  
#  Partitions <- Partition(data.frame(groups,LDValues),stratum="groups") 
#  allLD <- lapply( Partitions, function(x) { y <- x[2:(dim(x)[2])]; return( apply(y,2,sum) ) } ) 
#  allSD <- lapply( Partitions, function(x) { y <- x[2:(dim(x)[2])]; return( apply(y,2,var) ) } )
  D <- matrix(0.0,nrow=K,ncol=K)
  
  for(i in seq(1,K)) 
    for(j in seq(i,K)){
      if( i != j ) {
        p1 <- unlist( allLD[i,])
        p2 <- unlist( allLD[j,])
        D[i,j] <- D[j,i] <- sqrt( sum( (p1-p2)^2 ))        
      }
    }
  rownames(D) <- colnames(D) <- rownames(allLD)
  totMean <- mean(D)
  colMean <- colMeans(D)
  colMeanMatrix <- matrix(colMean,K,K,byrow=T)
  rowMeanMatrix <- matrix(colMean,K,K,byrow=F)
  C <- -0.5 * (D - colMeanMatrix - rowMeanMatrix + totMean)
  R <- matrix(1,K,K)
  
  for( i in 1:K ) for( j in 1:K ) if(i!=j)
    R[i,j] = C[i,j]/sqrt( C[i,i] * C[j,j] )
  
  # one MASS import
  SRI <- matrix(1,K,K)
  RI <- MASS::ginv(R)
  EED <- matrix(0,K,K)
  
  for(i in seq(1,K)) for(j in seq(1,K)) if(i!=j) 
    SRI[i,j] <- -1*RI[i,j]/sqrt( RI[i,i]*RI[j,j])
  
  SRI <- 1-SRI^2
  SRI[ SRI < 0 ] <- 0
  EED <- -N *log( SRI ) 
  
# EdgeStr <- -0.5 * log(SRI)
#  for(i in seq(1,K)) for(j in seq(1,K)) if(i!=j) {
#    if( SRI[i,j]^2 > 1 ) {
#      EED[i,j] = -N * log(1-SRI[i,j]^2)
#      EdgeStr[i,j] = -0.5 * log(1-SRI[i,j]^2)
#    }   
#  }  
  
  D[ EED<=critVal ] <- 0
  
  graph <- graph.adjacency(D,mode="undirected",weighted=TRUE,diag=FALSE)
  V(graph)$name <- row.names(D)  
  
#  popSD <- unlist( lapply( allSD, sum ) )
#  popSD <- scale( popSD, center=min(popSD), scale=TRUE) * 5 + 5    
  popSD <- scale( allSD, center=min(allSD), scale=TRUE) * 5 + 5
  V(graph)$size <- popSD
  


  class( graph ) <- c("igraph", "popgraph")
  
  return( graph )
}








