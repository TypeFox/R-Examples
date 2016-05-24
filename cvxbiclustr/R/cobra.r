#' Convex Biclustering via Dykstra-like Proximal Algorithm
#'
#' \code{cobra_internal} is an R wrapper around C code.
#'
#' @param X Data matrix
#' @param Lambda_row Initial guess of row Langrage multipliers
#' @param Lambda_col Initial guess of column Langrage multipliers
#' @param E_row Edge-incidence matrix for row graph
#' @param E_col Edge-incidence matrix for column graph
#' @param w_row Vector of weights for row graph
#' @param w_col Vector of weights for column graph
#' @param gamma Regularization parameter for shrinkage
#' @param max_iter Maximum number of AMA iterations
#' @param tol stopping tolerance
#' @useDynLib cvxbiclustr
cobra_internal <- function(X, Lambda_row,Lambda_col,E_row,E_col,w_row,w_col,gamma,
                           max_iter=1e2, tol=1e-3) {
  ## Get matrix dimensions
  m_row <- as.integer(nrow(E_row))
  m_col <- as.integer(nrow(E_col))
  n <- as.integer(ncol(X)); p <- as.integer(nrow(X))

  ## Unpack sparse matrix information for E_row
  column_ptr_row <- as.integer(E_row@p)
  values_row <- as.double(E_row@x)
  row_indices_row <- as.integer(E_row@i)

  ## Unpack sparse matrix information for E_col
  column_ptr_col <- as.integer(E_col@p)
  values_col <- as.double(E_col@x)
  row_indices_col <- as.integer(E_col@i)

  ## Cast types
  XT <- t(X)
  storage.mode(XT) <- "double"

  LambdaT_row <- t(Lambda_row)
  storage.mode(LambdaT_row) <- "double"
  LambdaT_temp_row <- matrix(0,m_row,n)
  storage.mode(LambdaT_temp_row) <- "double"
  LambdaT_old_row <- matrix(0,m_row,n)
  storage.mode(LambdaT_old_row) <- "double"
  dLambdaT_row <- matrix(0,m_row,n)
  storage.mode(dLambdaT_row) <- "double"
  gLambdaT_row <- matrix(0,m_row,n)
  storage.mode(gLambdaT_row) <- "double"
  gLambdaT_old_row <- matrix(0,m_row,n)
  storage.mode(gLambdaT_old_row) <- "double"

  LambdaT_col <- t(Lambda_col)
  storage.mode(LambdaT_col) <- "double"
  LambdaT_temp_col <- matrix(0,m_col,p)
  storage.mode(LambdaT_temp_col) <- "double"
  LambdaT_old_col <- matrix(0,m_col,p)
  storage.mode(LambdaT_old_col) <- "double"
  dLambdaT_col <- matrix(0,m_col,p)
  storage.mode(dLambdaT_col) <- "double"
  gLambdaT_col <- matrix(0,m_col,p)
  storage.mode(gLambdaT_col) <- "double"
  gLambdaT_old_col <- matrix(0,m_col,p)
  storage.mode(gLambdaT_old_col) <- "double"

  VT_row <- matrix(0,m_row,n)
  storage.mode(VT_row) <- "double"
  VT_col <- matrix(0,m_col,p)
  storage.mode(VT_col) <- "double"

  UT <- matrix(0,n,p)
  storage.mode(UT) <- "double"

  YT <- matrix(0,p,n)
  storage.mode(YT) <- "double"

  PT <- matrix(0,n,p)
  storage.mode(PT) <- "double"

  QT <- matrix(0,p,n)
  storage.mode(QT) <- "double"

  w_row <- as.double(gamma*w_row)
  w_col <- as.double(gamma*w_col)
  nu_row <- as.double(1/nrow(X))
  nu_col <- as.double(1/ncol(X))

  primal_row <- double(max_iter)
  dual_row <- double(max_iter)
  primal_col <- double(max_iter)
  dual_col <- double(max_iter)
  max_iter <- as.integer(max_iter)
  tol <- as.double(tol)

  #  void test_convex_bicluster_dlpa(double *xt,
  #                                  double *lambdat_row,
  #                                  double *lambdat_temp_row, double *lambdat_old_row, double *dlambdat_row,
  #                                  double *glambdat_row, double *glambdat_old_row,
  #                                  double *lambdat_col,
  #                                  double *lambdat_temp_col, double *lambdat_old_col, double *dlambdat_col,
  #                                  double *glambdat_col, double *glambdat_old_col,
  #                                  double *vt_row, double *vt_col,
  #                                  double *ut, double *yt, double *pt, double *qt,
  #                                  int *column_ptr_row, int *row_indices_row, double *values_row,
  #                                  int *column_ptr_col, int *row_indices_col, double *values_col,
  #                                  int *m_row, int *m_col, int *n, int *p,
  #                                  double *w_row, double *w_col,
  #                                  double *nu_row, double *nu_col,
  #                                  double *primal_row, double *dual_row,
  #                                  double *primal_col, double *dual_col,
  #                                  int *max_iter, int *iter, double *tol) {
  sol <- .C('test_convex_bicluster_dlpa',XT=XT,
            LambdaT_row=LambdaT_row,
            LambdaT_temp_row=LambdaT_temp_row,LambdaT_old_row=LambdaT_old_row,dLambdaT_row=dLambdaT_row,
            gLambdaT_row=gLambdaT_row,gLambdaT_old_row=gLambdaT_old_row,
            LambdaT_col=LambdaT_col,
            LambdaT_temp_col=LambdaT_temp_col,LambdaT_old_col=LambdaT_old_col,dLambdaT_col=dLambdaT_col,
            gLambdaT_col=gLambdaT_col,gLambdaT_old_col=gLambdaT_old_col,
            VT_row=VT_row,VT_col=VT_col,
            UT=UT,YT=YT,PT=PT,QT=QT,
            column_ptr_row=column_ptr_row,row_indices_row=row_indices_row,values_row=values_row,
            column_ptr_col=column_ptr_col,row_indices_col=row_indices_col,values_col=values_col,
            m_row=m_row,m_col=m_col,n=n,p=p,
            w_row=w_row,w_col=w_col,
            nu_row=nu_row,nu_col=nu_col,
            primal_row=primal_row,dual_row=dual_row,
            primal_col=primal_col,dual_col=dual_col,
            max_iter=max_iter,iter=integer(1),tol=tol)
  return(list(UT=sol$UT,YT=sol$YT,LambdaT_row=sol$LambdaT_row,VT_row=sol$VT_row,
              LambdaT_col=sol$LambdaT_col,VT_col=sol$VT_col,
              nu_row=sol$nu_row,nu_col=sol$nu_col,
              primal_row=sol$primal_row[1:sol$iter],dual_row=sol$dual_row[1:sol$iter],
              primal_col=sol$primal_col[1:sol$iter],dual_col=sol$dual_col[1:sol$iter],
              iter=sol$iter))
}

#' Convex Biclustering Algorithm
#'
#' \code{cobra} computes a convex biclustering path via Dykstra-like Proximal Algorithm
#'
#' @param X The data matrix to be clustered. The rows are the features, and the columns are the samples.
#' @param E_row Edge-incidence matrix for row graph
#' @param E_col Edge-incidence matrix for column graph
#' @param w_row Vector of weights for row graph
#' @param w_col Vector of weights for column graph
#' @param gamma A sequence of regularization parameters for row and column shrinkage
#' @param max_iter Maximum number of iterations
#' @param tol Stopping criterion
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
#' nGamma <- 5
#' gammaSeq <- 10**seq(0,3,length.out=nGamma)
#'
#' ## Generate solution path
#' sol <- cobra_validate(X,E_row,E_col,w_row,w_col,gammaSeq)
#'
#' ix <- 4
#' heatmap(sol$U[[ix]],col=hmcols,labRow=NA,labCol=NA,ColSideCol=cols[ty])
cobra <- function(X,E_row,E_col,w_row,w_col,gamma,max_iter=1e2,tol=1e-3) {
  call <- match.call()
  nGamma <- length(gamma)

  ## Get matrix dimensions
  m_row <- as.integer(nrow(E_row))
  m_col <- as.integer(nrow(E_col))
  n <- as.integer(ncol(X)); p <- as.integer(nrow(X))

  ## Check matrix dimensions
  check_E_row <- ncol(E_row) == p
  check_E_col <- ncol(E_col) == n
  check <- check_E_row && check_E_col

  ## Check weights vectors have correct length and are positive
  check_w_row <- (length(w_row) == m_row) && all(w_row > 0)
  check_w_col <- (length(w_col) == m_col) && all(w_col > 0)
  check <- check && check_w_row && check_w_col

  ## Check that m_row <= p*(p-1)/2 and m_col <= n*(n-1)/2
  check_m <- (m_row <= p*(p-1)/2) && (m_col <= n*(n-1)/2)
  check <- check && check_m

  if ( !check )
    stop("Dimensions of E_row, E_col, w_row, w_col, and X are not consistent.")

  ## Initialize dual variables
  Lambda_row <- matrix(0,n,m_row)
  Lambda_col <- matrix(0,p,m_col)
  list_U <- vector(mode="list",length=nGamma)
  list_V_col <- vector(mode="list",length=nGamma)
  list_V_row <- vector(mode="list",length=nGamma)
  list_Lambda_col <- vector(mode="list",length=nGamma)
  list_Lambda_row <- vector(mode="list",length=nGamma)
  list_cluster_col <- vector(mode="list",length=nGamma)
  list_size_col <- vector(mode="list",length=nGamma)
  list_cluster_row <- vector(mode="list",length=nGamma)
  list_size_row <- vector(mode="list",length=nGamma)

  iter_vec <- integer(nGamma)
  print("gamma  | row gap         col gap       max gap           ")
  print("---------------------------------------------------------")
  for (ig in 1:nGamma) {
    gam <- as.double(gamma[ig])
    sol <- cobra_internal(X,Lambda_row,Lambda_col,E_row,E_col,w_row,w_col,gam)

    ## Homotopy
    list_Lambda_col[[ig]] <- Lambda_col <- t(sol$LambdaT_col)
    list_Lambda_row[[ig]] <- Lambda_row <- t(sol$LambdaT_row)

    list_U[[ig]] <- t(sol$U)
    list_V_row[[ig]] <- t(sol$VT_row)
    list_V_col[[ig]] <- t(sol$VT_col)

    row_gap <- sol$primal_row[sol$iter]-sol$dual_row[sol$iter]
    col_gap <- sol$primal_col[sol$iter]-sol$dual_col[sol$iter]
    print(sprintf("%5d  | %5f        %5f      %7f", ig, signif(row_gap,4),
                  signif(col_gap,4), signif(max(row_gap,col_gap),4)))
  }
  return(list(U=list_U,V_row=list_V_row,V_col=list_V_col,
              Lambda_col=list_Lambda_col,Lambda_row=list_Lambda_row))
}
