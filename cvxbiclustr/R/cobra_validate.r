#' Get Subgroup Means
#'
#' \code{get_subgroup_means_full} computes the subgroup means.
#'
#' @param X Data matrix
#' @param clusters_row Row cluster object
#' @param clusters_col Column cluster object
#' @export
#' @seealso \code{\link[cvxbiclustr]{get_subgroup_means}}
get_subgroup_means_full <- function(X,clusters_row,clusters_col) {
  Y <- X
  num_clusters_row <- length(clusters_row$size)
  num_clusters_col <- length(clusters_col$size)
  M <- matrix(NA,num_clusters_row,num_clusters_col)
  for (i in 1:num_clusters_row) {
    ix_row <- which(clusters_row$cluster==i)
    for (j in 1:num_clusters_col) {
      ix_col <- which(clusters_col$cluster==j)
      M[i,j] <- mean(Y[ix_row,ix_col],na.rm=TRUE)
    }
  }
  return(M)
}

#' Get Subgroup Means
#'
#' \code{get_subgroup_means} computes the subgroup means on the training set for validation.
#'
#' @param X Data matrix
#' @param Theta validation index set - vector
#' @param clusters_row Row cluster object
#' @param clusters_col Column cluster object
#' @export
#' @seealso \code{\link[cvxbiclustr]{get_subgroup_means_full}}
get_subgroup_means <- function(X,Theta,clusters_row,clusters_col) {
  Y <- X
  Y[Theta] <- NA
  num_clusters_row <- length(clusters_row$size)
  num_clusters_col <- length(clusters_col$size)
  M <- matrix(NA,num_clusters_row,num_clusters_col)
  for (i in 1:num_clusters_row) {
    ix_row <- which(clusters_row$cluster==i)
    for (j in 1:num_clusters_col) {
      ix_col <- which(clusters_col$cluster==j)
      M[i,j] <- mean(Y[ix_row,ix_col],na.rm=TRUE)
    }
  }
  return(M)
}

#' Select Validation Set for a Matrix
#'
#' \code{get_validation} selects a random sample of matrix indices for regularization parameter selection.
#'
#' @param p Number of rows
#' @param n Number of columns
#' @param fraction Fraction of entries for hold out
#' @param seed Seed for random number generator
#' @export
#' @examples
#' n <- 5
#' p <- 4
#' fraction <- 0.1
#' Theta <- get_validation(p,n,fraction)
#'
#' M <- matrix(rnorm(n*p),p,n)
#' YT <- t(M)
#' YT[Theta$ThetaV] <- NA
#' Y <- t(YT)
#' Theta
get_validation <- function(p,n,fraction=0.1,seed=12345) {
  set.seed(seed)
  ix <- sample(x=1:(n*p),size=round(fraction*n*p),replace=FALSE)
  #  return(Theta=ix)
  rows <- ceiling(ix/n)
  cols <- ix - n*(rows-1)
  return(list(ThetaM=cbind(rows,cols),ThetaV=ix))
}

#' Perform validation to select regularization parameter
#'
#' \code{cobra_validate} performs an MM algorithm wrapper to do parameter selection.
#'
#' @param X Data matrix
#' @param E_row Edge-incidence matrix for row graph
#' @param E_col Edge-incidence matrix for column graph
#' @param w_row Vector of weights for row graph
#' @param w_col Vector of weights for column graph
#' @param gamma Regularization parameter for shrinkage
#' @param Lambda_row Initial guess of row Langrage multipliers
#' @param Lambda_col Initial guess of column Langrage multipliers
#' @param fraction Fraction of entries for hold out
#' @param max_iter Maximum number of iterations
#' @param tol Stopping criterion
#' @param max_iter_inner Maximum number of inner cobra iterations
#' @param tol_inner Stopping criterion for inner cobra loop
#' @useDynLib cvxbiclustr
#' @export
#' @examples
#' ## Create bicluster path
#' ## Example: Lung
#' X <- lung
#' X <- X - mean(X)
#' X <- X/norm(X,'f')
#'
#' ## Create annotation for heatmap
#' types <- colnames(lung)
#' ty <- as.numeric(factor(types))
#' cols <- rainbow(4)
#' YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
#' hmcols <- colorRampPalette(YlGnBu5)(256)
#'
#' ## Construct weights and edge-incidence matrices
#' phi <- 0.5; k <- 5
#' wts <- gkn_weights(X,phi=phi,k_row=k,k_col=k)
#' w_row <- wts$w_row
#' w_col <- wts$w_col
#' E_row <- wts$E_row
#' E_col <- wts$E_col
#'
#' ## Connected Components of Row and Column Graphs
#' wts$nRowComp
#' wts$nColComp
#'
#' #### Initialize path parameters and structures
#' nGamma <- 7
#' gammaSeq <- 10**seq(0,1,length.out=nGamma)
#'
#' ## Generate solution path
#' sol <- cobra_validate(X,E_row,E_col,w_row,w_col,gammaSeq,fraction=0.01)
#'
#' ## Plot validation error
#' verr <- sol$validation_error
#' plot(verr)
#'
#' ## Heatmap of data smoothed at the model selected to minimize validation error
#' ix <- which.min(verr)
#' groups_row <- sol$groups_row[[ix]]
#' groups_col <- sol$groups_col[[ix]]
#' M <- biclust_smooth(X,groups_row,groups_col)
#' heatmap(M,col=hmcols,labRow=NA,labCol=NA,ColSideCol=cols[ty])
cobra_validate <- function(X,E_row,E_col,w_row,w_col,gamma,
                           Lambda_row=matrix(0,n,nrow(E_row)),
                           Lambda_col=matrix(0,p,nrow(E_col)),
                           fraction=0.1,
                           max_iter=1e2,tol=1e-3,max_iter_inner=1e3,tol_inner=1e-4) {
  n <- ncol(X); p <- nrow(X)
  ThetaOut <- get_validation(p,n,fraction)
  ThetaM <- ThetaOut$ThetaM
  ThetaV <- ThetaOut$ThetaV
  U_last <- U <- matrix(0,p,n)
  nGamma <- length(gamma)
  UHx <- vector(mode='list',length=nGamma)
  VrHx <- vector(mode='list',length=nGamma)
  VcHx <- vector(mode='list',length=nGamma)

  ## Initialize variables for computing validation error
  validation_error <- double(nGamma)
  acc_count <- integer(nGamma)
  groups_row <- groups_col <- vector(mode='list',length=nGamma)

  XT <- t(X)
  Lambda_row <- matrix(0,n,nrow(E_row))
  Lambda_col <- matrix(0,p,nrow(E_col))
  for (ig in 1:nGamma) {
    sol <- cobra_pod(X, Lambda_row, Lambda_col, E_row, E_col, gamma[ig]*w_row, gamma[ig]*w_col, ThetaV,
                    max_iter=max_iter, tol=tol, max_iter_inner=max_iter_inner, tol_inner=tol_inner)
    UHx[[ig]] <- sol$U
    VcHx[[ig]] <- sol$V_col
    VrHx[[ig]] <- sol$V_row

    ## Compute Error on Validation Set
    clusters_row <- find_clusters(create_adjacency(sol$V_row,E_row))
    clusters_col <- find_clusters(create_adjacency(sol$V_col,E_col))
    groups_row[[ig]] <- clusters_row
    groups_col[[ig]] <- clusters_col
    MM <- get_subgroup_means(X,ThetaM[,1]+p*(ThetaM[,2]-1),clusters_row,clusters_col)
    MM[which(is.nan(MM))] <- 0
    ixi <- clusters_row$cluster[ThetaM[,1]]
    ixj <- clusters_col$cluster[ThetaM[,2]]
    errors <- double(nrow(ThetaM))
    for (i in 1:nrow(ThetaM)) {
      errors[i] <- MM[ixi[i],ixj[i]] - XT[ThetaV[i]]
    }
    validation_error[ig] <- sqrt(sum(errors**2))

    print(paste('*****','Completed gamma =',ig,'*****'))
  }
  return(list(U=UHx,V_row=VrHx,V_col=VcHx,ThetaM=ThetaM,ThetaV=ThetaV,
              groups_row=groups_row,groups_col=groups_col,
              validation_error=validation_error))
}

#' Bicluster Smooth
#'
#' \code{biclust_smooth} computes the bicluster estimates given row and column partitions
#'
#' @param X Data matrix
#' @param clusters_row Row cluster object
#' @param clusters_col Column cluster object
#' @export
#' @seealso \code{\link[cvxbiclustr]{get_subgroup_means_full}}
biclust_smooth <- function(X,clusters_row,clusters_col) {
  p <- nrow(X); n <- ncol(X)
  Y <- matrix(NA,p,n)
  M <- get_subgroup_means_full(X,clusters_row,clusters_col)
  num_clusters_row <- length(clusters_row$size)
  num_clusters_col <- length(clusters_col$size)
  for (i in 1:num_clusters_row) {
    ixi <- which(clusters_row$cluster == i)
    for (j in 1:num_clusters_col) {
      ixj <- which(clusters_col$cluster == j)
      Y[ixi,ixj] <- M[i,j]
    }
  }
  return(Y)
}
