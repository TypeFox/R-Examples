## excursions.int.R
##
##   Copyright (C) 2012, 2013, 2014 David Bolin, Finn Lindgren
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.

gaussint <- function(mu,
                     Q.chol,
                     Q,
                     a,
                     b,
                     lim = 0,
                     n.iter = 10000,
                     ind,
                     use.reordering = c("natural","sparsity","limits"),
                     max.size,
                     max.threads=0,
                     seed,
                     LDL=FALSE)
{

  if( missing(Q) && missing(Q.chol))
	  stop('Must specify a precision matrix or its Cholesky factor')

  if(missing(a))
    stop('Must specify lower integration limit')

  if(missing(b))
    stop('Must specify upper integration limit')

  if(!missing(mu))
    mu <- private.as.vector(mu)

  if(!missing(ind))
    ind <- private.as.vector(ind)

  a <- private.as.vector(a)
  b <- private.as.vector(b)

  if(!missing(Q))
    Q <- private.as.Matrix(Q)

  if(!missing(Q.chol))
    Q.chol <- private.as.Matrix(Q.chol)


  n = length(a)
  if(length(b) != n)
    stop('Vectors with integration limits are of different length.')

  use.reordering <- match.arg(use.reordering)

  if(!missing(ind) && !is.null(ind)){
    a[!ind] = -Inf
    b[!ind] = Inf
  }

  if(missing(max.size))
    max.size = n

  reordered = FALSE
  if(!missing(Q.chol) && !is.null(Q.chol)){
    ## Cholesky factor is provided, use that and do not reorder
      L = Q.chol
      if(dim(L)[1] != dim(L)[2]){
        stop("Q.chol is not symmetric")
      } else if(dim(L)[1] != n) {
        stop("Dimensions of Q.chol is different from the length of the integration limits.")
      }
  } else if(!missing(Q) && !is.null(Q)){
    if(dim(Q)[1] != dim(Q)[2]){
      stop("Q is not symmetric")
      } else if(dim(Q)[1] != n) {
        stop("Dimensions of Q is different from the length of the integration limits.")
      }
    ## Cholesky factor is not provided and we are told to reorder
    if(use.reordering == "limits")
    {
      #Reorder by moving nodes with limits -inf ... inf first
      inf.ind = (a==-Inf) & (b==Inf)
      if(sum(inf.ind)>0){
        max.size = min(max.size,n-sum(inf.ind))
        n = length(a)
        cind = rep(1,n)
        cind[inf.ind] = 0
        reo = rep(0,n)
			  out <- .C("reordering",nin = as.integer(n), Mp = as.integer(Q@p),
		                        Mi = as.integer(Q@i), reo = as.integer(reo),
		                        cind = as.integer(cind))
		    reo = out$reo+1
		    Q = Q[reo,reo]
		    reordered = TRUE
		  }

      if(LDL) {
  		  L = suppressWarnings(t(as(Cholesky(Q,perm=FALSE),"Matrix")))
		  } else {
		     L = chol(private.as.spam(Q),pivot=FALSE)
		  }
		} else if(use.reordering == "sparsity"){
		  #Reorder for sparsity, let spam do it...
		  if(LDL) {
  		  L = suppressWarnings(t(as(Cholesky(Q,perm=FALSE),"Matrix")))
		    reo = L@perm
		    ireo[reo] = 1:length(reo)
		  } else {
		     L = chol(private.as.spam(Q))
  		  reo = L@pivot
	  	  ireo = L@invpivot
		  }
		  reordered = TRUE
		} else {
		  #Do not reorder
      if(LDL) {
        L = suppressWarnings(t(as(Cholesky(Q,perm=FALSE),"Matrix")))
		  } else {
		    L = chol(private.as.spam(Q),pivot=FALSE)
		  }
		}
	}

  # Note: If lim > 0 and reorder == TRUE, we should calculate marginal
  # probabilities, see if bound is above lim, and then reorder

  if(!missing(mu) && !is.null(mu)){
    if(length(mu) != n){
      stop("The length of mu is different from the length of the integration limits.")
    }
    a = a - mu
    b = b - mu
  }

  a[a==Inf]  = .Machine$double.xmax
  b[b==Inf]  = .Machine$double.xmax
  a[a==-Inf] = -.Machine$double.xmax
  b[b==-Inf] = -.Machine$double.xmax

	if(reordered == TRUE){
		a = a[reo]
		b = b[reo]
  }

  if(is(L,'spam.chol.NgPeyton')){
     L = as(as(spam::as.dgRMatrix.spam(spam::as.spam(L)), "TsparseMatrix"),"dtCMatrix")
  } else if (is(L, "Matrix")) {
     L <- as(as(L, "CsparseMatrix"), "dtCMatrix")
  } else {
    stop("Unsuported matrix type.")
  }

  if(!missing(seed) && !is.null(seed)){
    seed_provided = 1
    seed.in = seed
  } else {
    seed_provided = 0
    seed.in = as.integer(rep(0,6))
  }

  Pv = Ev = rep(0,dim(L)[1])

  opts = c(n,n.iter,max.size,max.threads,seed_provided)

  out <- .C("shapeInt", Mp = as.integer(L@p), Mi = as.integer(L@i),
              Mv = as.double(L@x), a = as.double(a), b = as.double(b),
              opts = as.integer(opts), lim = as.double(lim),
              Pv = as.double(Pv), Ev = as.double(Ev),seed_in=seed.in)

  if(use.reordering == "limits") {
    out$Pv[1:(dim(L)[1]-max.size)] = out$Pv[dim(L)[1]-max.size+1]
    out$Ev[1:(dim(L)[1]-max.size)] = out$Ev[dim(L)[1]-max.size+1]
  } else if(use.reordering == "sparsity") {
    out$Pv = out$Pv[ireo]
    out$Ev = out$Ev[ireo]
  }
  return(list(Pv = out$Pv, Ev = out$Ev, P = out$Pv[1], E = out$Ev[1]))
}
