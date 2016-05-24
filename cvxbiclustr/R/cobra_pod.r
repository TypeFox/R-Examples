#' MM algorithm for Convex Biclustering with Missing Data
#'
#' \code{cobra_pod} performs convex biclustering on incomplete data matrices using an MM algorithm.
#'
#' @param X The data matrix to be clustered. The rows are the features, and the columns are the samples.
#' @param Lambda_row Initial guess of row Langrage multipliers
#' @param Lambda_col Initial guess of column Langrage multipliers
#' @param E_row Edge-incidence matrix for row graph
#' @param E_col Edge-incidence matrix for column graph
#' @param w_row Vector of weights for row graph
#' @param w_col Vector of weights for column graph
#' @param Theta A vector of missing indices - row major order
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
#' ## Generate random initial dual variables
#' set.seed(12345)
#' n <- ncol(X); p <- nrow(X)
#' m_row <- nrow(E_row); m_col <- nrow(E_col)
#' Lambda_row <- matrix(rnorm(n*m_row),n,m_row)
#' Lambda_col <- matrix(rnorm(p*m_col),p,m_col)
#'
#' #### Initialize path parameters and structures
#' gam <- 200
#'
#' ## Create random mask
#' nMissing <- floor(0.1*n*p)
#' Theta <- sample(1:(n*p), nMissing, replace=FALSE)
#'
#' sol <- cobra_pod(X,Lambda_row,Lambda_col,E_row,E_col,gam*w_row,gam*w_col,Theta)
#'
#' heatmap(sol$U,col=hmcols,labRow=NA,labCol=NA,ColSideCol=cols[ty])
cobra_pod <- function(X, Lambda_row, Lambda_col, E_row, E_col, w_row, w_col, Theta,
                        max_iter=1e2, tol=1e-3, max_iter_inner=1e3, tol_inner=1e-4) {

  ## Get matrix dimensions
  m_row <- as.integer(nrow(E_row))
  m_col <- as.integer(nrow(E_col))
  n <- as.integer(ncol(X)); p <- as.integer(nrow(X))

  ## Check matrix dimensions
  check_Lambda_row <- (nrow(Lambda_row) == n) && (ncol(Lambda_row) == m_row)
  check_Lambda_col <- (nrow(Lambda_col) == p) && (ncol(Lambda_col) == m_col)
  check_E_row <- ncol(E_row) == p
  check_E_col <- ncol(E_col) == n
  check <- check_Lambda_row && check_Lambda_col && check_E_row && check_E_col

  ## Check weights vectors have correct length and are positive
  check_w_row <- (length(w_row) == m_row) && all(w_row > 0)
  check_w_col <- (length(w_col) == m_col) && all(w_col > 0)
  check <- check && check_w_row && check_w_col

  ## Check that m_row <= p*(p-1)/2 and m_col <= n*(n-1)/2
  check_m <- (m_row <= p*(p-1)/2) && (m_col <= n*(n-1)/2)
  check <- check && check_m

  if ( !check )
    stop("Lambda_row, Lambda_col, E_row, E_col, and X do not have consistent dimensions.")
  call <- match.call()

  ## Cast types
  XT <- t(X)
  storage.mode(XT) <- "double"

  UT <- t(X)
  storage.mode(UT) <- "double"

  LambdaT_row <- t(Lambda_row)
  storage.mode(LambdaT_row) <- "double"

  LambdaT_col <- t(Lambda_col)
  storage.mode(LambdaT_col) <- "double"

  VT_row <- matrix(0,m_row,n)
  storage.mode(VT_row) <- "double"
  VT_col <- matrix(0,m_col,p)
  storage.mode(VT_col) <- "double"

  Theta <- as.integer(Theta-1)
  nMissing <- as.integer(length(Theta))

  ## Unpack sparse matrix information for E_row
  column_ptr_row <- as.integer(E_row@p)
  values_row <- as.double(E_row@x)
  row_indices_row <- as.integer(E_row@i)

  ## Unpack sparse matrix information for E_col
  column_ptr_col <- as.integer(E_col@p)
  values_col <- as.double(E_col@x)
  row_indices_col <- as.integer(E_col@i)

  w_row <- as.double(w_row)
  w_col <- as.double(w_col)
  nu_row <- as.double(1/nrow(X))
  nu_col <- as.double(1/ncol(X))

  max_iter <- as.integer(max_iter)
  tol <- as.double(tol)
  max_iter_inner <- as.integer(max_iter_inner)
  tol_inner <- as.double(tol_inner)

  #  void test_convex_bicluster_impute(double *mt, double *ut,
  #                                    double *lambdat_row, double *lambdat_col,
  #                                    double *vt_row, double *vt_col,
  #                                    int *column_ptr_row, int *row_indices_row, double *values_row,
  #                                    int *column_ptr_col, int *row_indices_col, double *values_col,
  #                                    int *m_row, int *m_col, int *n, int *p,
  #                                    int *Theta, int *nMissing,
  #                                    double *w_row, double *w_col,
  #                                    double *nu_row, double *nu_col,
  #                                    double *mm_loss,
  #                                    int *max_iter, int *iter, double *tol,
  #                                    int *max_iter_inner, double *tol_inner) {
  sol <- .C('test_convex_bicluster_impute',MT=XT,UT=UT,
            LambdaT_row=LambdaT_row,
            LambdaT_col=LambdaT_col,
            VT_row=VT_row,VT_col=VT_col,
            column_ptr_row=column_ptr_row,row_indices_row=row_indices_row,values_row=values_row,
            column_ptr_col=column_ptr_col,row_indices_col=row_indices_col,values_col=values_col,
            m_row=m_row,m_col=m_col,n=n,p=p,
            Theta=Theta,nMissing=nMissing,
            w_row=w_row,w_col=w_col,
            nu_row=nu_row,nu_col=nu_col,
            mm_loss=double(max_iter),
            max_iter=max_iter,iter=integer(1),tol=tol,
            max_iter_inner=max_iter_inner,tol_inner=tol_inner)
  return(list(U=t(sol$UT),Lambda_row=t(sol$LambdaT_row),V_row=t(sol$VT_row),
              Lambda_col=t(sol$LambdaT_col),V_col=t(sol$VT_col),
              nu_row=sol$nu_row,nu_col=sol$nu_col,
              mm_loss=sol$mm_loss[1:sol$iter],
              call=call))
}
