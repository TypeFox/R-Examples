#' Alternating Least Squares with Groups for Constrained Dual Scaling
#' 
#' Alternating least-sqaures for estimating row and column scores in constrained
#' dual scaling, where different groups are allowed for.
#' 
#' @param a A \code{2n}-vector of row scores.
#' @param m Integer; the number of items.
#' @param q Integer; the rating scale from \code{1:q}.
#' @param G An indicator matrix of size \code{n} by \code{K}. 
#' @param Fr.cent The centred F_r matrix.
#' @param eps The numerical tolerance level for the loss.
#' @param maxit Integer; the maximum number of iterations allowed.
#' @param Mmat Matrix of spline basis functions.
#' @param info.level Integer controlling the amount of information printed.
#' @param const The constant part of the loss function.
#' @param K The number of latent classes.
#' @param n The number of samples.
#' @param tol tolerance \code{tol} passed to \code{\link{lsei}} of the
#' \pkg{limSolve} package
#' @keywords multivariate
#' @export group.ALS
group.ALS <-
function (a, m, q, G, Fr.cent, eps = 1e-1, maxit = 50, Mmat, 
          info.level = 2, const, K, n, tol) {
time1 <- proc.time()[3]
Vmat <- matrix(NA, nrow = m + q - 1, ncol = K)
alphamat <- matrix(NA, ncol = K, nrow = 4)
crit <- crit.scaled <- rep(NA, maxit)
iter <- 1
a.upd <- a

while(iter <= maxit)
	{
	awts <-  sqrt(colSums(matrix(a.upd*a.upd,ncol=K,nrow=2*n)*rbind(G,G)))
	Vmat <-  apply(rbind(G,G), 2, function(x, a, Fmat) 2*crossprod(Fmat,matrix(a*x,ncol=1))/
                      (sqrt(as.numeric(crossprod(a,a*x)))*(m + q - 2)), 
                    a = a.upd, Fmat = Fr.cent)
    ### Check dimensions of G and a !!!!!
  
#     apply(G, 2, function(x,a,Fmat) 2*t(Fmat)%*%matrix(a*x,ncol=1)/(sqrt(sum(a*a*x))*(m + q - 2)), 
#                   a = a.upd, Fmat = Fr.cent)
  # dim Vmat: (m + q - 1) x K
  
	### Update bstar
	b1.hat <- rowSums(matrix(apply(rbind(awts, Vmat[1:m,]), 2, function(x) x[1]*x[-1]), ncol = K, nrow = m))/crossprod(a.upd)
#     sum(a.upd*a.upd)
  # rowSums matrix: m x K
  
	### Update alphamat
  # might add try()? PCS 25/10/2013
	alphamat <- apply(rbind(awts, matrix(Vmat[-(1:m), ], ncol = K)), 2, function(x, Mmat, tol) 
		limSolve::lsei(A = x[1] * Mmat / max(awts), B = x[-1] / max(awts), G = cbind(0, diag(3)), H = rep(0, 3), 
            type = 2, tol = tol)$X, Mmat = Mmat, tol = tol)
  # dim alphamat: 4 x K
	
	### update a
	bkmat <- apply(alphamat, 2, function(x, b, Mmat) c(b, Mmat %*% matrix(x, ncol = 1)),
	               b = b1.hat, Mmat = Mmat)
	# dim bkmat: m + q - 1 x K
  bprod <- colSums(bkmat * bkmat)
	denom <- rowSums(matrix(1/bprod, byrow = TRUE, nrow = n, ncol = K)*G)
	denom <- c(denom, denom)
  Fr.bk <- Fr.cent %*% bkmat
	a.upd <- 2 * rowSums(rbind(G, G) * Fr.bk) * denom/(m + q - 2)
  
  # rescale Dgk*a to unit length, and rescale bk similarly
#   if(rescale.grp){
#     wts <- apply(rbind(G,G), 2, function(x, a) sd((a*x)[x != 0]), a = a.upd)
#     tmp <- matrix(wts, byrow = T, nrow = 2*n, ncol = K)
#     a.upd <- rowSums(rbind(G,G)*a.upd/tmp)
#     alphamat <- alphamat*matrix(wts, nrow = 4, ncol = K)
#     bkmat <- bkmat*matrix(wts,nrow = m + q - 1, ncol = K)
#     bprod <- colSums(bkmat*bkmat)
#     rm(tmp, wts)
#   }

	### Calculate loss
	awts2 <- colSums(matrix(a.upd * a.upd, ncol = K, nrow = 2 * n) * rbind(G, G))
	last <- sum((matrix(a.upd, ncol = K, nrow = 2 * n) * rbind(G, G)) * (Fr.bk))
	crit[iter] <- const + 0.25 * (m + q - 2) * (m + q - 2) * crossprod(bprod, awts2) - (m + q - 2) * last
	crit.scaled[iter] <- 2 * crit[iter] / (n * (m + q - 1) * (m + q - 2)^2)
	
	if(iter > 1 && abs(crit[iter - 1] - crit[iter]) < eps) {
		if(crit[iter - 1] - crit[iter] < 0 & abs( 1 - crit[iter - 1]/crit[iter]) > sqrt(.Machine$double.eps)) {
      warning("ALS: Not monotone decreasing")
      print(crit[!is.na(crit)], digits = 15)
		}
		if(info.level==2 || info.level == 4) {
      time2 <- proc.time()[3]
      cat("\t\t\t\t ALS loss = ", crit[iter], "\t | \t", round(time2 - time1, 2), "sec's \t | \t", 
          iter, "iter's | \n")
		  }
    break
		}		
	if(iter == maxit) {
		# warning("ALS: Maximum number of iterations reached")
    if(info.level == 2 || info.level == 4)
      {
		  time2 <- proc.time()[3]
		  cat("\t\t\t\t ALS loss = ", crit[iter], "\t | \t", round(time2 - time1,2), "sec's \t | \t", 
            iter, "iter's | *** \n")
      }
		break
		}
	iter <- iter + 1
	}

list(b1 = b1.hat, alpha = alphamat, a = a.upd, bk = bkmat, crit = crit[1:iter],
	crit.scaled = crit.scaled[1:iter], iter = iter, bwts2 = bprod, 
      awts2 = awts2, last = last, Fr.bk = Fr.bk, crit.opt = crit[iter])
}
