## utils.R
##
##   Copyright (C) 2013 David Bolin, Finn Lindgren
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




## calculates variances given either a cholesky factor L in Matrix or spam f
## format, or given a precision matrix Q. If Q is provided, the cholesky factor
## is calculated and the variances are then returned in the same ordering as Q
## If L is provided, the variances are returned in the same ordering as L, even
## if L@invpivot exists.
excursions.variances<-function(L,Q, max.threads=0, use.spam=FALSE)
{

  if(!missing(L) && !is.null(L)){
    ireo = FALSE
    if(is(L,'spam.chol.NgPeyton')){
     L = as(as(spam::as.dgRMatrix.spam(spam::as.spam(L)), "TsparseMatrix"),"dtCMatrix")
   } else {
      if (!is(L, "dtCMatrix"))
       stop("L needs to be in ccs format for now.")
    }
  } else {
    if (requireNamespace("spam", quietly=TRUE)){
      L = chol(private.as.spam(Q))
      ireo = TRUE
      reo = L@invpivot
      L = as(as(spam::as.dgRMatrix.spam(spam::as.spam(L)),
                "TsparseMatrix"),"dtCMatrix")
    } else {
      L = Cholesky(Q,LDL=FALSE)
      ireo = TRUE
      reo <- 1:(L@Dim[1])
      reo[L@perm+1] <- 1:(L@Dim[1])
    }
  }

  out<- .C("Qinv",Rir = as.integer(L@i),
                  Rjc = as.integer(L@p),
                  Rpr = as.double(L@x),
                  variances=double(dim(L)[1]),
                  n = as.integer(dim(L)[1]),
                  n_threads = as.integer(max.threads))

  if(ireo){
    return(out$variances[reo])
  } else {
    return(out$variances)
  }
}


excursions.marginals <- function(type, rho,vars, mu, u, QC = FALSE)
{
  rl = list()
  if(type == "=" || type == "!="){
    if(QC){
      rl$rho_ngu = rho
      rl$rho_l = pnorm(mu-u,sd=sqrt(vars), lower.tail=FALSE)
			rl$rho_u = 1-rl$rho_l
			rl$rho = pmax(rl$rho_u,rl$rho_l)
			rl$rho_ng = pmax(rl$rho_ngu,1-rl$rho_ngu)
    } else {
      if(!missing(rho)){
        rl$rho_u = rho
      } else {
        rl$rho_u = 1 - pnorm(mu-u,sd=sqrt(vars), lower.tail=FALSE)
		  }
		  rl$rho_l = 1-rl$rho_u
		  rl$rho = pmax(rl$rho_u,rl$rho_l)
    }
  } else {
    if(QC) {
      rl$rho_ng = rho
      if(type == ">"){
			  rl$rho = pnorm(mu-u,sd=sqrt(vars))
			} else {
				rl$rho = pnorm(mu-u,sd=sqrt(vars), lower.tail=FALSE)
			}
    } else {
      if(missing(rho)){
       if(type == ">"){
  				rl$rho = pnorm(mu-u,sd=sqrt(vars))
       } else {
  				rl$rho = pnorm(mu-u,sd=sqrt(vars), lower.tail=FALSE)
        }
      } else {
        rl$rho = rho
      }
    }
  }
  return(rl)
}


excursions.permutation <- function(rho, ind, use.camd = TRUE,alpha,Q)
{
  if(!missing(ind) && !is.null(ind))
    rho[!ind] = -1
  n = length(rho)
  v.s = sort(rho,index.return=TRUE)
  reo = v.s$ix
  rho_sort = v.s$x
  ireo = integer(n)
  ireo[reo] = 1:n
  if(use.camd){
    k=0
    i = n-1
    #add nodes to lower bound
    cindr = cind = rep(0,n)
    while(rho_sort[i]> 1 - alpha && i>0){
      cindr[i] = k
			i = i-1
			k=k+1
    }
    if(i>0){
      #reorder nodes below the lower bound for sparsity
      while (i>0) {
				cindr[i] = k;
				i = i-1
			}
      #change back to original ordering
      for (i in 1:n) {
				cind[i] = k - cindr[ireo[i]]
			}
			#call CAMD
			out <- .C("reordering",nin = as.integer(n), Mp = as.integer(Q@p),
		                        Mi = as.integer(Q@i), reo = as.integer(reo),
		                        cind = as.integer(cind))
		  reo = out$reo+1
    }
  }
  return(reo)
}


excursions.setlimits <- function(marg, vars,type,QC,u,mu)
{
	if(QC){
		if (type=="<") {
		  uv = sqrt(vars)*qnorm(pmin(pmax(marg$rho_ng,0),1))
		} else if (type==">"){
			uv = sqrt(vars)*qnorm(pmin(pmax(marg$rho_ng,0),1),lower.tail=FALSE)
		} else if (type=="=" || type == "!="){
			uv = sqrt(vars)*qnorm(pmin(pmax(marg$rho_ngu,0),1),lower.tail=FALSE)
		}
	} else {
		uv = u-mu
	}
	if (type == "=" || type == "!=") {
	  if(QC){
	    a = b = uv
	    a[marg$rho_ngu<=0.5] = -Inf
	    b[marg$rho_ngu>0.5] = Inf
		} else {
		  a = b = uv
		  a[marg$rho_u<=0.5] = -Inf
      b[marg$rho_u>0.5] = Inf
		}
	} else if (type == ">") {
	  a = uv
	  b = rep(Inf,length(mu))
	} else if (type == "<"){
    a = rep(-Inf,length(mu))
    b = uv
  }
  return(list(a=a,b=b))
}



excursions.call <- function(a,b,reo,Q, is.chol = FALSE, lim, K, max.size,n.threads, seed,LDL=FALSE)
{
  if(is.chol == FALSE){
    a.sort = a[reo]
    b.sort = b[reo]
    Q = Q[reo,reo]

    if(LDL){
  		L = suppressWarnings(t(as(Cholesky(Q,perm=FALSE),"Matrix")))
    } else {
      L = suppressWarnings(chol(private.as.spam(Q),pivot = FALSE))
    }

    res = gaussint(Q.chol = L, a = a.sort, b = b.sort, lim = lim,
                                 n.iter = K, max.size = max.size,
                                 max.threads = n.threads, seed = seed)
  } else {
    #assume that everything already is ordered
    res = gaussint(Q = Q, a= a, b = b, lim = lim, n.iter = K,
                    max.size = max.size,
                    max.threads = n.threads,seed = seed,LDL=LDL)
  }
  return(res)
}

private.as.spam <- function(A)
{
  if(is(A,"spam")){
    return(A)
  }
  else if(is(A,"dsyMatrix")){
    return(spam::as.spam.dgCMatrix(as(as.matrix(A),"dgCMatrix")))
  } else {
    return(spam::as.spam.dgCMatrix(as(A,"dgCMatrix")))
  }
}


private.check.integer <- function(v)
{
  if(is.null(v)) {
    stop("Anticipated scalar value, got NULL")
  } else if(!is.null(dim(v))){
    stop("Anticipated scalar value, got matrix")
  } else if(length(v)>1) {
    stop("Anticipated scalar value, got vector")
  }
}

private.as.vector <- function(v)
{
  if(is.null(v) || is.vector(v)) {
    return(v)
  } else {
    if(min(dim(v)>1)) {
      stop("vector has wrong dimensions")
    }
    return(as.vector(v))
  }
  return(c(v))
}

private.as.Matrix <- function(M)
{
  if(is.null(M) || is(M,"Matrix")){
    return(M)
  } else {
    return(as(M,"Matrix"))
  }
}

##
# Distribution function of Gaussian mixture \Sum_k w[k]*N(mu[k],sigma[k]^2)
##
Fmix <- function(x,mu,sd,w) sum(w*pnorm(x,mean=mu,sd=sd))

##
# Quantile function of Gaussian mixture
##
Fmix_inv <- function(p,mu,sd,w,br=c(-1000,1000))
{
   G = function(x) Fmix(x,mu,sd,w) - p
   return(uniroot(G,br)$root)
}

##
# Function for optimization of interval for mixtures
##
fmix.opt <- function(x,
                     alpha,
                     sd,
                     Q.chol,
                     w,
                     mu,
                     limits,
                     verbose,
                     max.threads,
                     ind)
{
  K = dim(mu)[1]
  n = dim(mu)[2]
  q.a = sapply(seq_len(n),function(i) Fmix_inv(x/2,
                                               mu = mu[,i],
                                               sd = sd[,i],
                                               w=w,
                                               br=limits))
  q.b = sapply(seq_len(n),function(i) Fmix_inv(1-x/2,
                                               mu = mu[,i],
                                               sd = sd[,i],
                                               w=w,
                                               br=limits))

  prob = 0
  stopped = 0

  k.seq = sort(w,decreasing=TRUE,index.return=TRUE)$ix

  ki = 1
  for(k in k.seq){
    ws = 0
    if(ki<K){
      ws = sum(w[k.seq[(ki+1):K]])
    }
    ki = ki+1
    lim = (1-alpha - prob - ws)/w[k]
    p = gaussint(mu = mu[k,],
                 Q.chol=Q.chol[[k]],
                 a=q.a,
                 b=q.b,
                 ind = ind,
                 lim=max(0,lim),
                 max.threads=max.threads)
    if(p$P==0){
      stopped = 1
      break
    } else {
      prob = prob + w[k]*p$P
    }
  }


  if(stopped==1){ #too large alpha
    if(prob == 0){
      val = 10*(1+x)
    } else {
      val = 1+(prob-(1-alpha))^2
    }
  } else { #too small x
    val = (prob-(1-alpha))^2
  }

  if(verbose){
    cat("in optimization: ",x," ", prob, " ", val, "\n")
  }

  return(val)

}


fmix.samp.opt <- function(x, alpha,mu, sd, w, limits, samples)
{
  n = dim(mu)[2]
  q.a = sapply(seq_len(n),function(i) Fmix_inv(x/2,
                                               mu = mu[,i],
                                               sd = sd[,i],
                                               w=w,
                                               br=limits))
  q.b = sapply(seq_len(n),function(i) Fmix_inv(1-x/2,
                                               mu = mu[,i],
                                               sd = sd[,i],
                                               w=w,
                                               br=limits))

  cover = sapply(seq_len(dim(samples)[1]), function(i) (sum(samples[i,]>q.b) + sum(samples[i,]<q.a))==0)

  prob = mean(cover)
  val = (prob-(1-alpha))^2
  cat("in optimization: ",x," ", prob, " ", val, "\n")
  return(val)

}

fsamp.opt <- function(x, samples,verbose=FALSE)
{
  q.a = apply(samples,1,quantile,1,probs=c(x/2))
  q.b = apply(samples,1,quantile,1,probs=c(1-x/2))
  prob = mean(apply((samples<q.b)*(samples>q.a),2,prod))
  if(verbose)
    cat("in optimization: ",x," ", prob, "\n")
  return(prob)

}


mix.sample <- function(n.samp = 1, mu,Q.chol,w)
{
  K = length(mu)
  n = length(mu[[1]])
  idx = sample(seq_len(K), n.samp, prob = w, replace = TRUE)
  idx = sort(idx)
  n.idx = numeric(K)
  n.idx[] = 0
  for (i in 1:K) {
    n.idx[i] = sum(idx == i)
  }
  samples = c()
  for (k in 1:K){
    if (n.idx[k] > 0) {
      xx = mu[[k]] + solve(Q.chol[[k]],matrix(rnorm(n.idx[k]*n),n,n.idx[k]))
      samples = rBind(samples,t(as.matrix(xx)))
    }
  }
  return(samples)
}

excursions.rand <- function(n,seed,n.threads=1)
{
  if(!missing(seed) && !is.null(seed)){
    seed_provided = 1
    seed.in = seed
  } else {
    seed_provided = 0
    seed.in = as.integer(rep(0,6))
  }

  x = rep(0,n)
  opt = c(n,n.threads,seed_provided)
  out<- .C("testRand",opt  = as.integer(opt),
                      x = as.double(x),
                      seed_in  = as.integer(seed.in))

  return(out$x)
}



## Turn off all warnings for require(), to allow clean completion of
## examples that require unavailable Suggested packages.
require.nowarnings <- function(package, lib.loc = NULL, character.only = FALSE)
{
  if (!character.only)
    package <- as.character(substitute(package))
  op <- options("warn")
  on.exit(options(op))
  options(warn = -1)
  require(package, lib.loc = lib.loc, quietly = TRUE, character.only = TRUE)
}



excursions.marginals.mc <- function(X,type, rho, mu, u)
{
  rl = list()
  if(type == "=" || type == "!="){
      if(!missing(rho)){
        rl$rho_u = rho
      } else {
        rl$rho_u = 1 - rowMeans(X<u)
		  }
		  rl$rho_l = 1-rl$rho_u
		  rl$rho = pmax(rl$rho_u,rl$rho_l)
  } else {
      if(missing(rho)){
       if(type == ">"){
  				rl$rho = rowMeans(X>u)
       } else {
  				rl$rho = rowMeans(X<u)
        }
      }
    }
  return(rl)
}

mcint <- function(X,
                  a,
                  b,
                  ind)
{

  if(missing(a))
    stop('Must specify lower integration limit')

  if(missing(b))
    stop('Must specify upper integration limit')

  n = length(a)
  if(length(b) != n)
    stop('Vectors with integration limits are of different length.')

  if(!missing(ind) && !is.null(ind)){
    a[!ind] = -Inf
    b[!ind] = Inf
  }

  Pv = rowMeans(apply(apply(apply(a < X & X < b,2,rev),2,cumprod),2,rev))

  #Estimate of MC error, not implemented yet
  Ev <-  rep(0,n)

  return(list(Pv = Pv, Ev = Ev, P = Pv[1], E = Ev[1]))
}


summary.excurobj <- function(object,...)
{
  out <- list()
  class(out) <- "summary.excurobj"
  out$calculation = object$meta$calculation
  out$call = object$meta$call
    if(object$meta$calculation == "excursions"){
      if(object$meta$type == ">"){
        out$computation = "Positive excursion set, E_{u,alpha}^+"
      } else if(object$meta$type == "<"){
         out$computation = "Negative excursion set, E_{u,alpha}^-"
      } else if(object$meta$type == "="){
        out$computation = "Contour credible region, E_{u,alpha}^c"
      } else {
        out$computation = "Contour avoiding set, E_{u,\alpha}"
      }
      out$u = object$meta$level
      out$alpha = object$meta$alpha
      out$F.limit = object$meta$F.limit
      out$method = object$meta$method
    } else if(object$meta$calculation == "simconf"){
      out$computation = "Simultaneous confidence band"
      out$alpha = object$alpha
    } else if(object$meta$calculation == "contourmap"){
      out$computation = "Contour map"
      out$u = object$u
      out$type = object$meta$contourmap.type
      out$F.computed = object$meta$F.computed
      if(out$F.computed)
        out$F.limit = object$meta$F.limit

      if(is.null(object$P0) && is.null(object$P1) && is.null(object$P2) &&
         is.null(object$P0.bound) && is.null(object$P1.bound) &&
         is.null(object$P2.bound)){
      } else {
        out$measures = list()
        if(!is.null(object$P0))
          out$measures$P0 = object$P0

        if(!is.null(object$P1))
          out$measures$P1 = object$P1

       if(!is.null(object$P2))
          out$measures$P2 = object$P2

        if(!is.null(object$P0.bound))
          out$measures$P0.bound = object$P0.bound

        if(!is.null(object$P1.bound))
          out$measures$P1.bound = object$P1.bound

       if(!is.null(object$P2.bound))
          out$measures$P2.bound = object$P2.bound
    }
  }
  return(out)
}


print.summary.excurobj <- function(x,...)
{

  cat("Call: \n")
  print(x$call)
  cat("\nComputation:\n")
  cat(x$computation,"\n\n")

  if(x$calculation == "excursions"){
    cat("Level: u = ")
    cat(x$u,"\n")
    cat("Error probability: alpha = ")
    cat(x$alpha,"\n")
    cat("Limit for excursion function computation: F.limit = ")
    cat(x$F.limit,"\n")
    cat("Method used : ")
    cat(x$method,"\n")
  } else if(x$calculation == "simconf"){
    cat("Error probability: alpha = ")
    cat(x$alpha,"\n")
  } else if(x$calculation == "contourmap"){
    cat("Level: u = ")
    cat(x$u,"\n")
    cat("Type of contour map: ")
    cat(x$type,"\n")
    if(x$F.computed) {
      cat("Contour map function computed\n")
      cat("Limit for excursion function computation: F.limit = ")
      cat(x$F.limit,"\n")
    } else {
      cat("Contour map function not computed\n")
    }
    cat("Quality measures computed : ")
    if(is.null(x$measures)){
      cat("none\n")
    } else {
      for(i in 1:length(x$measures)){
        cat(names(x$measures)[i], " = ", x$measures[[i]],"\n")
      }
    }
  }
}


print.excurobj <- function(x,...) {
  print.summary.excurobj(summary(x))
}
