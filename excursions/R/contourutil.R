## Calculate the contour map function
contourfunction.mc <- function(lp,mu,X,ind, alpha, verbose=FALSE)
{

  if(missing(mu))
    stop('Must supply mu')

  if(missing(lp))
    stop('Must supply level plot')

	F.limit = 1

  if(verbose) cat("set limits\n")

	lim <- excursions.limits(lp=lp,mu=mu,measure=0)

	m.size = length(mu)
  indices = NULL

  if (!missing(ind)) {
    if(is.logical(ind)){
      indices = ind
      m.size = sum(ind)
    } else {
      indices = rep(FALSE,length(mu))
      indices[ind] = TRUE
      m.size = length(ind)
    }
  }
  if(verbose) cat("calculate marginals\n")

  rho <- contourmap.marginals.mc(X=X,lim=lim,ind=indices)


  reo <- excursions.permutation(rho, indices, use.camd = FALSE)

  if(verbose) cat("integrate\n")

  res = mcint(X=X[reo,],a=lim$a[reo],b=lim$b[reo])

  n = length(mu)
  ii = which(res$Pv[1:n] > 0)
  if (length(ii) == 0) i=n+1 else i=min(ii)

  F = Fe  = E = rep(0,n)
  F[reo] = res$Pv
  Fe[reo] = res$Ev

  ireo = NULL
  ireo[reo] = 1:n

  ind.lowF = F < 1-F.limit
  E[F>1-alpha] = 1
  F[ind.lowF] = Fe[ind.lowF] = NA

  M = rep(-1,n)
  for(i in 1:(lp$n.levels+1)){
    M[(lp$G == (i-1)) & (E == 1)] = i-1
  }

  return(list(F=F, Fe=Fe, E=E, M=M, rho=rho))
}


## Calculate the contour map function
contourfunction <- function(lp,mu,Q,vars,ind, alpha, n.iter=10000,
                            F.limit, Q.chol,max.threads=0,
                            seed=seed,verbose=FALSE)
{
	if (!missing(Q.chol) && !is.null(Q.chol)) {
      Q = Q.chol
      is.chol = TRUE
  } else if (!missing(Q) && !is.null(Q)) {
      is.chol = FALSE
  } else {
    stop('Must supply Q or Q.chol')
  }

  if(missing(mu))
    stop('Must supply mu')

  if(missing(lp))
    stop('Must supply level plot')

  if(missing(F.limit) || is.null(F.limit))
    F.limit = 1

  if(!missing(alpha) && !is.null(alpha))
	  F.limit = max(alpha,F.limit)

  if(verbose) cat("set limits\n")
	lim <- excursions.limits(lp=lp,mu=mu,measure=0)


  if (missing(vars)) {
    if(verbose) cat("calculate variances\n")
    if(is.chol) {
      vars <- excursions.variances(L=Q)
    } else {
      vars <- excursions.variances(Q=Q)
    }
  }

	m.size = length(mu)
  indices = NULL

  if (!missing(ind)) {
    if(is.logical(ind)){
      indices = ind
      m.size = sum(ind)
    } else {
      indices = rep(FALSE,length(mu))
      indices[ind] = TRUE
      m.size = length(ind)
    }
  }
  if(verbose) cat("calculate marginals\n")
  rho <- contourmap.marginals(mu=mu,vars=vars,lim=lim,ind=indices)

  lim$a <- lim$a - mu
	lim$b <- lim$b - mu

    if(verbose) cat("calculate permutation\n")
  use.camd = !missing(ind) || alpha < 1
  reo <- excursions.permutation(rho = rho, ind = indices,
                                use.camd = use.camd,alpha = F.limit,Q = Q)

    if(verbose) cat("integrate\n")
  res <- excursions.call(lim$a,lim$b,reo,Q, is.chol = is.chol,
                         1-F.limit, K = n.iter, max.size = m.size,
                         n.threads = max.threads,seed=seed)

  n = length(mu)
  ii = which(res$Pv[1:n] > 0)
  if (length(ii) == 0) i=n+1 else i=min(ii)

  F = Fe  = E = rep(0,n)
  F[reo] = res$Pv
  Fe[reo] = res$Ev

  ireo = NULL
  ireo[reo] = 1:n

  ind.lowF = F < 1-F.limit
  E[F>1-alpha] = 1
  F[ind.lowF] = Fe[ind.lowF] = NA

  M = rep(-1,n)
  for(i in 1:(lp$n.levels+1)){
    M[(lp$G == (i-1)) & (E == 1)] = i-1
  }

  return(list(F=F, Fe=Fe, E=E, M=M, rho=rho))
}

## Calculate marginal probabilities P(lim$a < X < lim$b) for
## when X is N(mu,vars)-distributed
contourmap.marginals <- function(mu,vars,lim,ind)
{
  if(!missing(ind) && !is.null(ind)){
	  marg = rep(0,length(mu))
	  marg[ind] = pnorm(lim$b[ind], mu[ind], sqrt(vars[ind])) -
	              pnorm(lim$a[ind], mu[ind], sqrt(vars[ind]))
	} else {
	  marg = pnorm(lim$b, mu, sqrt(vars)) - pnorm(lim$a, mu, sqrt(vars))
	}
	return(marg)
}

## Calculate marginal probabilities P(lim$a < X < lim$b) for
## when X
contourmap.marginals.mc <- function(X,lim,ind)
{
  if(!missing(ind) && !is.null(ind)){
	  marg = rep(0,length(lim$a))
	  marg[ind] = rowMeans(lim$a[ind] < X[ind,] & X[ind,] < lim$b[ind])
	} else {
	  marg = rowMeans(lim$a < X & X < lim$b)
	}
	return(marg)
}

## Create a levelplot with given levels/number of levels
excursions.levelplot <- function(mu,n.levels,ind,levels,
                                 equal.area=FALSE,
                                 pretty.cm=FALSE)
{
	n = length(mu)
	if(missing(ind)) ind = 1:n
	x.mean = rep(-Inf,n)
	x.mean[ind] = mu[ind]
	r = range(mu[ind])
	if(missing(levels)){
		if(missing(n.levels)){
			stop('Must specify levels or number of levels')
		} else {
		  private.check.integer(n.levels)
		}
		if(equal.area){
			levels <- as.vector(quantile(mu[ind],
								(1:n.levels)/(n.levels+1),type=4))
		} else if(pretty.cm){
      levels <- pretty(mu[ind], n = n.levels)
      n.levels = length(levels)
		} else {
			levels <- seq(from=r[1],to=r[2],
						  length.out = (n.levels+2))[2:(n.levels+1)]
		}
	} else {
		if(!missing(n.levels)){
			if(n.levels != length(levels)){
			  warning('n.levels is not equal to the length of levels')
				n.levels = length(levels)
			}
		} else {
			n.levels = length(levels)
		}
	}

	u.e = NULL
	E = vector("list",n.levels+1)

  if(n.levels == 1){
    l1 = c(levels-diff(r),levels,levels+diff(r))
  } else {
  	l1 = c(2*levels[1]- levels[2],
	         levels,
	         2*levels[n.levels]-levels[n.levels-1])
  }
	for(i in 1:(n.levels+1)) u.e[i] = (l1[i]+l1[i+1])/2

	for(i in 1:(n.levels)){
		E[[i]] = which((l1[i] <= x.mean) & (x.mean < l1[i+1]))
	}
	E[[n.levels+1]] = which((l1[n.levels+1] <= x.mean))
	map = rep(0,n)
	G = rep(-1,n)
	for(i in 1:(n.levels+1)){
	  map[E[[i]]] = u.e[i]
	  G[E[[i]]] = i-1
	}
	return(list(u = levels,
	            n.levels = n.levels,
	            u.e = u.e,
	            E=E,
	            map=map,
	            G=G))
}

## Create a P-optimal levelplot.
## The function will take A LOT of time to run if use.marginals=FALSE.
excursions.opt.levelplot <- function(mu,vars,Q,n.levels, measure=2, use.marginals=TRUE, ind)
{
	if( (measure != 1) && (measure != 2) && (measure != 0))
		stop('only measure 0, 1, or 2 allowed')

	if(missing(ind)) ind=1:length(mu)

	r = range(mu[ind])

	u <- seq(from=r[1],to=r[2],length.out = (n.levels+2))[2:(n.levels+1)]
	if(n.levels==1){
  	l = diff(r)
	} else {
	  l <- diff(u[1:2])
	}

  Q.chol = chol(Q)
  u.add = seq(from=-l, to = l,length.out = 19)
  P.add = NULL
  for(jj in 1:length(u.add)){
    P.add[jj] = -restricted.lim.func(u.add = u.add[jj],
                                     u0 = u,
                                     mu = mu,
                                     vars = vars,
                                     Q.chol = Q.chol,
                                     Q=Q,
                                     measure = measure,
                                     use.marginals = use.marginals,
                                     ind=ind)
  }
  plot(u.add,P.add)
  k = which.max(P.add)
  u.add = u.add[k]
  P.k = P.add[k]

  u <- u.add + u
	lp = excursions.levelplot(mu,levels = u,ind=ind)
	if(measure==2){
		if(use.marginals){
			lp$P2.bound = P.k
		} else {
			lp$P2 = P.k
		}
	} else if(measure==1){
		if(use.marginals){
			lp$P1.bound = P.k
		} else {
			lp$P1 = P.k
		}
	} else if(measure==0){
		if(use.marginals){
			lp$P0.bound = P.k
		} else {
			lp$P0 = P.k
		}
	}
	return(lp)
}
## Internal function for optimization of restricted P-optimal contour map
restricted.lim.func <- function(u.add, u0, mu, vars, Q.chol, Q, measure,
                                use.marginals, ind=ind)
{
	lp = excursions.levelplot(mu=mu,levels = (u.add + u0),ind=ind)
  if(use.marginals){
	  val = -Pmeasure.bound(lp=lp,mu=mu,vars=vars,type=measure,ind=ind)
		cat(u0, ': ', -val, '\n')
	} else {
		val = -Pmeasure(lp,mu=mu,Q=Q,Q.chol=Q.chol,type=measure,
		                ind=ind,vars=vars)
	  cat(u.add, ': ', -val, '\n')
	}
	return(val)
}


## Internal function for optimization of P-optimal contour map
excursions.lim.func <- function(u, mu, vars, Q.chol, Q, measure,
                                use.marginals, ind=ind)
{
	lp = excursions.levelplot(mu,levels = u,ind=ind)
	if( min(u)<=min(mu[ind]) | max(u) >= max(mu[ind])){
	  # levels should be in (min(mu),max(mu))
	  val = 0
	} else if (max(sort(u,index.return=TRUE)$ix - seq_len(length(u)))>0) {
	  # levels should be sorted
	  val = 0
	} else {
	  v = TRUE;
	  if(length(u)>1){
  	  for(i in 1:(length(u)-1)){
	      v = v & (sum(mu[ind]<u[i+1] & mu[ind] > u[i])>0)
	    }
	  }
	  if(!v) {
	    # all sets E should be non-empty
	    val = 0
	  } else {
  		if(use.marginals){
		  	val = -Pmeasure.bound(lp=lp,mu=mu,vars=vars,type=measure,ind=ind)
		    cat(u, ': ', -val, '\n')
		  } else {
			  val = -Pmeasure(lp,mu=mu,Q=Q,Q.chol=Q.chol,type=measure,
			                ind=ind,vars=vars)
	      cat(u, ': ', -val, '\n')
		  }
		}
	}
	return(val)
}

Pmeasure.bound <- function(lp, mu, vars, type, ind=NULL)
{
  limits = excursions.limits(lp,mu,measure=type)
  if(type==0){
    return(mean(contourmap.marginals(mu,vars,limits,ind)[ind]))
  } else {
	  return(min(contourmap.marginals(mu,vars,limits,ind)[ind]))
  }
}

## Function that calculates the P measure for a given contour map.
Pmeasure <- function(lp,mu,Q,Q.chol, ind=NULL,type,vars=vars)
{
  if(type==0){
    res <- contourfunction(lp=lp,mu=mu,Q=Q,vars=vars,ind=ind)
    p = mean(res$F[ind])
  } else {
    if(type == 1 && length(lp$u) == 1){
      return(1)
    }
    limits = excursions.limits(lp=lp,mu=mu,measure=type)
	  res = gaussint(mu = mu, Q=Q, Q.chol = Q.chol, a=limits$a,
	                b=limits$b,ind=ind,use.reordering="limits")
	  p = res$P[1]
  }
	return(p)
}

## Function that calculates the P measure for a given contour map.
Pmeasure.mc <- function(lp,mu,X, ind=NULL,type)
{
  if(type==0){
    res <- contourfunction.mc(lp=lp,X=X,ind=ind)
    p = mean(res$F[ind])
  } else {
    if(type == 1 && length(lp$u) == 1){
      return(1)
    }
    limits = excursions.limits(lp=lp,mu=mu,measure=type)
    res = mcint(X=X, a=limits$a, b=limits$b,ind=ind)
	  p = res$P[1]
  }
	return(p)
}

## Set integration limits for a given measure.
excursions.limits <- function(lp,mu,measure)
{
	n = length(mu)
	n.l = length(lp$u)

	a = rep(-Inf,n)
	b = rep(Inf,n)

	if(measure==2){
		b[lp$E[[1]]] = lp$u.e[2]
		a[lp$E[[n.l+1]]] = lp$u.e[n.l]

		if(n.l>1){
			for(i in 2:n.l){
				a[lp$E[[i]]] = lp$u.e[i-1]
				b[lp$E[[i]]] = lp$u.e[i+1]
			}
		}
	} else if(measure ==1){
		if(n.l>1){
			b[lp$E[[1]]] = lp$u[2]
			a[lp$E[[n.l+1]]] = lp$u[n.l-1]
		}
		if(n.l>2){
			b[lp$E[[2]]] = lp$u[3]
			a[lp$E[[n.l]]] = lp$u[n.l-2]
		}
		if(n.l>3){
			for(i in 3:(n.l-1)){
				a[lp$E[[i]]] = lp$u[i-2]
				b[lp$E[[i]]] = lp$u[i+1]
			}
		}
	} else if(measure==0){

		b[lp$E[[1]]] = lp$u[1]
		a[lp$E[[n.l+1]]] = lp$u[n.l]

		if(n.l>1){
			for(i in 2:n.l){
				a[lp$E[[i]]] = lp$u[i-1]
				b[lp$E[[i]]] = lp$u[i]
			}
		}
	} else {
		stop('Measure must be 0, 1, or 2')
	}
	return(list(a=a,b=b))
}

contourmap.colors <- function(lp,zlim,col,credible.col)
{

  if(missing(zlim) || is.null(zlim)){
    zlim = lp$meta$mu.range
  }

  u.e <- sort(unique(lp$map)) #extracts only used breakpoints
  u.e[1] = max(u.e[1],zlim[1])
  u.e[length(u.e)] = min(u.e[length(u.e)],zlim[2])

  breaks <- seq(zlim[1]-1e-16, zlim[2]+1e-16,
                length.out = length(col)+1)
  cmap = col[cut(c(u.e), breaks)@.Data]
  if(!missing(credible.col) && !is.null(credible.col))
    cmap = c(credible.col,cmap)

  return(cmap)
}



