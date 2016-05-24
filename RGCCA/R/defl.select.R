#' The function defl.select() computes residual matrices \eqn{\mathbf{X}_{1,h+1}, \ldots, \mathbf{X}_{J,h+1}}. These 
#' residual matrices are determined according to the following formula: \eqn{\mathbf{X}_{j,h+1} = \mathbf{X}_{jh} - \mathbf{y}_{jh} \mathbf{p}_{jh}^t}.
#' @param yy  A matrix that contains the SGCCA block components of each block: \eqn{\mathbf{y}_{1h}, \ldots, \mathbf{y}_{Jh}}
#' @param rr  A list that contains the residual matrices \eqn{\mathbf{X}_{1h}, \ldots, \mathbf{X}_{Jh}} 
#' @param nncomp A \eqn{1 \times J} vector that contains the number of components to compute for each block.
#' @param nn  A \eqn{1 \times J} vector that contains the numbers of already computed components for each block
#' @param nbloc Number of blocks.
#' @return \item{resdefl}{A list of \eqn{J} elements that contains \eqn{\mathbf{X}_{1,h+1}, \ldots, \mathbf{X}_{J,h+1}}.}
#' @return \item{pdefl}{A list of \eqn{J} elements that contains \eqn{\mathbf{p}_{1h}, \ldots, \mathbf{p}_{Jh}}.}
#' @title deflation function
#' @export defl.select

defl.select <- function(yy,rr,nncomp,nn,nbloc) {
  resdefl <- NULL
  pdefl <- NULL
  for (q in 1:nbloc) {
    if ( nn <= nncomp[q] ) {
       defltmp <- deflation(rr[[q]],yy[ , q])
       resdefl[[q]] <- defltmp$R
       pdefl[[q]]   <- defltmp$p
    } else {
        resdefl[[q]] <- rr[[q]]
        pdefl[[q]]   <-	rep(0,NCOL(rr[[q]]))
    }
  }
  return(list(resdefl=resdefl,pdefl=pdefl))
}
