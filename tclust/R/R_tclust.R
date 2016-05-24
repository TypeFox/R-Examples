

##	get a matrix object out of the sigma - tensor
.ssclmat <- function (x, k) as.matrix (x[,,k])

##	calculates the initial cluster sizes
.getini <- function (K, no.trim)
{
	if (K == 1)
		return (no.trim)

	pi.ini  <- runif(K)
	ni.ini <- sample(x = K, size = no.trim, replace = T, prob = pi.ini / sum (pi.ini))
	return (tabulate(ni.ini, nbins = K))
}

##	calculates the initial cluster assignment and parameter
.InitClusters <- function (X, iter, pa)
{
	for (k in 1:pa$K)
	{
		idx <- sample (1:pa$n, pa$p+1)
		X.ini = X [drop = F,idx,]#sample (1:pa$n, pa$p+1),]	##	selecting observations randomly for the current init - cluster
		iter$center [k,] <- colMeans (X.ini)				#n	(moved)	calculating the center

#n		iter$sigma[,,k] <- cov_fact (X.ini, iter$center [k,], pa$p/(pa$p+1))	#n	using the "new" covariance function, which is numerically more similiar to the C++ - version
		iter$sigma[,,k] <- (pa$p/(pa$p+1))*cov (X.ini)		##	calculating sigma (cc = current cov)
#o		iter $center [k,] <- colMeans (X.ini)				##	--> moved up
	}

	if (pa$equal.weights)									##	if we're considering equal weights, cw is set here **ONLY ONCE AND NEVER CHANGED**
		iter$csize <- rep (pa$no.trim / pa$K, pa$K)
	else
		iter$csize = .getini (pa$K, pa$no.trim)
	iter$cw <- iter$csize / pa$no.trim
									##		if we're considering different weights, calculate them, and they're gonna be recalculated every time in .estimClustPar
	return (iter)
}

.estimClustPar <- function (X, iter, pa)
{				

		for (k in 1:pa$K)
		{
			if (iter$csize[k] > pa$zero.tol)				##	this cluster's size is > 0
			{
				iter$center[k,] = (t(iter$z_ij[,k]) %*% X) / iter$csize[k]
				X.c <- (X - matrix (iter$center[k,], ncol = pa$p, nrow = pa$n, byrow = TRUE))
				iter$sigma[,,k] <- (t(X.c * iter$z_ij[,k]) %*% X.c) / iter$csize[k]
			}
			else											##	this cluster's size has decreased to 0
				iter$sigma[,,k] <- 0						##	do the same as in the non - fuzzy approach
		}

	return (iter)
}

.div_by0_set0 <- function (x, y)								##	a function which divides two elements, setting x/0 := 0
{
	ret <- x / y
	ret [y == 0] <- 0
	ret
}

.calcFuzzy_row <- function (x, pa)							#n	new function calculating a row of iter$Z, corresponding to a single observation
{															##	2DO: instead of <= 0 we might use pa$zero.tol ?!
    idx.max = which.max (x)

	if (x[idx.max] >= 1)
		return (as.numeric (1:pa$K == idx.max))						##	this gives an 0, 0, 1, 0, 0 vector for k = 5, idx.max = 3

	x.log <- -log (x)
	x.log [x <= 0] <- 0												##	avoid -Inf. Note, that x < 1. thus log (x) can never be 0!

	la2 <- matrix (x.log, ncol = pa$K, nrow = pa$K)
	la3 <- t (la2)
	z2_3d <- .div_by0_set0 (la3, la2)								##	this is 2d now actually...

	z2_ij <- colSums (z2_3d^(1/(pa$m-1)))
	z2_ij <- .div_by0_set0 (1, z2_ij)

#o	if (sum (z2_ij) <= 0)
#o		return (as.numeric (1:pa$K == idx.max))						##	this was a wrong solution to an issue discovered in the original version
	if (sum (z2_ij) <= pa$zero.tol)									#n	use zero.tol instead!
		return (rep (1/pa$K, pa$K))									#n	changed on 20110216

	return (z2_ij^pa$m)
}

.findClustAssig <- function (X, iter, pa)
{															##	finds the current cluster assignment based on the given center and sigma
	ll = matrix (NA, pa$n, pa$K)

	for (k in 1:pa$K)										##	calculating the probability (or density ?!) for each observation in each cluster.
		ll[,k] <- iter$cw[k] * .dmnorm(X, iter$center[k,], .ssclmat (iter$sigma,k)) ## dmvnorm could be used here...

	if (pa$fuzzy > 0)
	{
        old.z_ij <- iter$z_ij
		iter$assig <- apply(ll,1,which.max)

		if (pa$fuzzy >= 2)											#n	another possibility of calculating z_ij, which is closer related to the C++ implementation.
			iter$z_ij <- t (apply (ll, 1, .calcFuzzy_row, pa = pa))	#n apply the calculation for an observation (".calcFuzzy_row") for each row of the ll - vector.
		else
		{
	#n        iter$assig <- apply(ll, 1, select.clust, pa = pa)	#n	##	searching the cluster which fits best for each observation

			## to obtain FUZZY values of z_ij
			iter$z_ij <- matrix (0, ncol = pa$K, nrow = pa$n)	##	calculating "absolute z_ij" for !pa$fuzzy, which only contains 0 and 1

#o			z1_ij <-	matrix (0, ncol = pa$K, nrow = pa$n)	##  INITIALIZATION
#o			z2_ij <-	matrix (0, ncol = pa$K, nrow = pa$n)	##	INITIALIZATION
#o			pre.z_ij <-	matrix (0, ncol = pa$K, nrow = pa$n)	##	INITIALIZATION

			log.ll <- pre.z_ij <- z2_ij <- z1_ij <- iter$z_ij	#n  INITIALIZATION (these matrices are all the same at the beginning...)

			assig_ <- cbind (1:pa$n, iter$assig)
			z1_ij[assig_] <- 1
			max.ll <- apply (ll, 1, max)
			max.ll.greater1 <- max.ll >= 1
	#o		max.ll.greater1_ <-matrix (max.ll.greater1, nrow = pa$n, ncol = pa$K, byrow = FALSE)	#	not needed anymore, we use "max.ll.greater1" directly as column indices (e.g. z1_ij[max.ll.greater1, ])
	#o		log.ll <- ll * 0		##	initialized at the beginngin (-> "INITIALIZATION")

	#o		log.ll[ll > 0] <- -log(ll[ll > 0]) * (ll[ll > 0] < 1)
			ll.bzo <- ll > 0 & ll < 1							#n	ll between zero and one
			log.ll[ll.bzo] <- -log(ll[ll.bzo])					#n

	#o		ll2 <- log.ll %*% matrix (rep (diag (pa$K), pa$K), nrow=pa$K)	##	repeat matrix log.ll pa$K-times (columnwise)
	#i		ll2 <- matrix (log.ll, nrow = nrow (ll), ncol = ncol (ll) * pa$K)
			la2 <- array (log.ll, dim = c(pa$n, pa$K, pa$K))				#n	a (n x K x K) array. this is easier to handly

	#o		yy <- matrix (0, nrow = pa$K, ncol = pa$K * pa$K)
	#o		for (ii in 1:pa$K) 
	#o			yy[ii,((ii-1)*pa$K+1):(ii*pa$K)]=rep(1,pa$K)
	#o		ll3 <- log.ll %*% yy
			la3 <- aperm (la2, c (1, 3, 2))									##	an array of same dimension as la2, but the 2nd and 3rd dimensions are swapped (this corresponds to yy) (this is done with function "aperm")

	#o		z2_ij__ <- (((ll3 * (ll2 > 0) / (ll2 + (ll2 == 0)))^(1/(pa$m-1))) %*% t(yy))	##	this has been split into the 3 following "#i" lines

	#i		z2_ij__ <- (ll3 / ll2)^(1/(pa$m-1))
	#i		z2_ij__ [ll2 == 0] <- 0
			z2_3d <- .div_by0_set0 (la3, la2)

	#i		z2_ij__ <- z2_ij__  %*% t(yy)
			z2_ij <- apply (z2_3d^(1/(pa$m-1)), c (1, 3), sum)

	#o		z2_ij <- 1*(z2_ij__  > 0) / (z2_ij__ + (z2_ij__ == 0))
	#i		z2_ij <- .div_by0_set0 (1, z2_ij__)
			z2_ij <- .div_by0_set0 (1, z2_ij)

	#o		sum_z2_ij = cbind (apply (z2_ij, 1, 'sum'))				##	we use rowSums instead
			sum_z2_ij <- rowSums (z2_ij)

	#o?		z2_ij[sum_z2_ij == 0,] <- assig_[sum_z2_ij == 0,]		##	is this really supposed to be assig_? or rather z1_ij ? note that assig_ contains values up to pa$n
	#o		z2_ij[sum_z2_ij == 0,] <- z1_ij[sum_z2_ij == 0,]		#i	i think it's supposed to be like this.. WRONG!!
			z2_ij[sum_z2_ij <= pa$zero.tol,] <- 1/pa$K				#n	this is the proper solution to this issue. changed on 20110216

	#o		pre.z_ij <- z1_ij * max.ll.greater1_ + z2_ij * (max.ll.greater1_ == 0)
			pre.z_ij[max.ll.greater1, ] <- z1_ij[max.ll.greater1, ]			#n
			pre.z_ij[!max.ll.greater1, ] <- z2_ij[!max.ll.greater1, ]		#n

			iter$z_ij <- pre.z_ij^pa$m
		}

        ##PUT TRIMMING in iter$assig and in iter$z_ij
        ll[ll < 0] <- 0
        log.ll_ <- log(ll)

		log.ll_[iter$z_ij == 0] <- 0							#n BUG !? this is a new line.. in order to avoid 0*-Inf division in iter$z_ij * log.ll_.
#o        disc <- apply (iter$z_ij * log.ll_, 1, sum)			#o use rowSums instead

        disc <- rowSums (iter$z_ij * log.ll_)

#o      idx.out <- rank (disc, ties.method = c ("random"))		#o	this is not an idx.out, but a rank - variable. let's use a "real" idx.out - variable
		idx.out <- which (rank (disc, ties.method = c ("random"))<=pa$trim)	#i	now this is a "real" idx.out - vector
#n		idx.out <- select.outliers (pa, disc)					#n	selecting the trimmed observations with special treatment of ties (2come)

#o		iter$assig [idx.out<=pa$trim] <- 0						##	setting outliers to class 0
#o		iter$z_ij[idx.out<=pa$trim,] <- 0

        iter$assig [idx.out] <- 0								#n	setting outliers to class 0
        iter$z_ij[idx.out,] <- 0								#n	using the "real" idx.out - vector

        ##CHECKING EQUALITY IN Z_IJ
        iter$code <- all (old.z_ij == iter$z_ij)				##	setting the code - parameter, signaling whether the assig - array is still the same

		iter$csize <- colSums (iter$z_ij)						##	calculating the cluster size
	}
	else
	{
		old.assig <- iter$assig
		iter$assig <- apply(ll,1,which.max)						#o	searching the cluster which fits best for each observation
#n		iter$assig <- apply(ll, 1, select.clust, pa = pa)		##	searching the cluster which fits best for each observation

		disc <-ll [cbind (1:pa$n, iter$assig)]

#o		idx.out=rank(disc,ties.method = c( "random"))			#o	this is not an idx.out, but a rank - variable. let's use a "real" idx.out - variable
		idx.out <- which (rank (disc, ties.method = c ("random"))<=pa$trim)	#i	now this is a "real" idx.out - vector
#n		idx.out <- select.outliers (pa, disc)					#n	selecting the trimmed observations with special treatment of ties (2come)

#o		iter$assig [idx.out<=pa$trim]=0							#o	setting outliers to class 0
		iter$assig [idx.out] <- 0								#n	using the "real" idx.out - vector	#	setting outliers to class 0

		iter$code <- all (old.assig == iter$assig)				##	setting the code - parameter, signaling whether the assig - array is still the same
		iter$csize <- tabulate (iter$assig, pa$K)
		iter$z_ij<- matrix (0, ncol = pa$K, nrow = pa$n)		##	calculating "absolute z_ij" for !pa$fuzzy, which only contains 0 and 1

#o		iter$z_ij[cbind ((1:pa$n)[idx.out>pa$trim], iter$assig[idx.out>pa$trim])] <- 1		
		iter$z_ij[cbind (1:pa$n, iter$assig)[-idx.out, ]] <- 1	#n	using the "real" idx.out - vector	

#n		idx.max <- cbind (1:pa$n, iter$assig)					##
#n		iter$z_ij[idx.max[,idx.out>pa$trim]] <- 1
	}

	if (!pa$equal.weights)										##	calculate cluster weights (if necessary)
		iter$cw <- iter$csize / sum(iter$csize)

	return (iter)
}

.calcobj <- function (X, iter, pa)
{			
		##	calculates the obj. functions value
#	iter$obj <- ifelse (pa$equal.weights, 0, sum ( (iter$csize *log(iter$cw)  )[iter$csize != 0] ) )
	iter$obj <- ifelse (pa$equal.weights, 0, sum ( (iter$csize *log(iter$cw)  )[iter$csize > pa$zero.tol] ) )

	for (k in 1:pa$K)
	{
		w <- .dmnorm(X,iter$center[k,],.ssclmat (iter$sigma,k))

		if (sum(iter$z_ij[  ,k]) >  pa$zero.tol)
#o			if (sum(iter$z_ij[!w,k]) <= pa$zero.tol)
			if (sum(iter$z_ij[w <= 0, k]) <= pa$zero.tol)
				iter$obj <- iter$obj + sum(iter$z_ij[w > 0,k] * log(w[w > 0]))  
			else 
			{
				iter$obj <- iter$obj -Inf
				return (iter)					#n	we can already return -Inf, no need to do any further calculation
			}
	}
	return (iter)
}

.TreatSingularity <- function (iter, pa) 
{	
	warning ("After trimming, all points in the data set are concentrated in k subspaces.") ##  a single point is a subspace too.

	iter$code <- 2	# indicating the data's concentration in either k subspaces or points.
	return (iter)
}

#################						##	a development-function. This it takes a function argument (f.restr) for specifying the restriction.
##	.tclust.R  ##
#################

.tclust.R <- function (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, f.restr = .restr.diffax, equal.weights = FALSE, fuzzy = FALSE, m = 2, zero.tol = 1e-16,  trace = 0, store.x = TRUE, f.hook.iter, f.hook.model, ...)
{
	if (is.data.frame(x))
		x <- data.matrix(x)
	else if (!is.matrix (x))
		x <- matrix(x, nrow = length(x), ncol = 1, dimnames = list (names(x), deparse(substitute(x))))
	if (!is.numeric (x))
		stop ("parameter x: numeric matrix/vector expected")

	parlist <- list (k = k, alpha = alpha, nstart = nstart, iter.max = iter.max, f.restr = f.restr, equal.weights = equal.weights, fuzzy = fuzzy, m = m, zero.tol = zero.tol, trace = trace, store.x = store.x, ...)

	if (store.x)
		parlist$x <- x

	n <- nrow (x)
	p <- ncol (x)
	no.trim <- floor(n * (1 - alpha))

	# preparing lot's of lists	
	pa <- list (									##	these are variables which all the iterations have in common, and will never change (pa for "params")
		n = n,													##	number of observations
		p = p,													##	number of dimensions
		no.trim = no.trim,										##	number of observations which are considered as to be not outlying
		trimm = n-no.trim,										##	number of observations which are considered as to be outlying
		K = k,													##	number of clusters to be searched for
		equal.weights = equal.weights,							##	wether equal weights shall be assumed for all clusters
		zero.tol = zero.tol,									##	zero tolerance				##	XXXS yet to be implemented
		trace = trace,											##	trace - level giving more information on the iterations
		fuzzy = fuzzy,											##	whether fuzzy clustering shall be performed
		m = m													##  ?? value specifying what ??
	)

	iter <- list (									##	these variables change each iteration - this object is passed to all the functions, modified and returned by them.
		obj = -Inf,												##	current objective value
		assig = array (0, n),									##	cluster assignment
		csize = array (NA, k),									##	cluster sizes
		cw = rep (NA, k),										##	cluster weights 
		sigma = array (NA, c (p, p, k)),						##	cluster's sigmas
		center = array (NA, c(k, p)),							##	cluster's centers
		code = NA,												##	this is a return code supplied by functions like .findClustAssig
		z_ij = matrix (0, nrow = n, ncol = k )					##	z_ij ## -> what was it's old name - "ind" - right?
	)

#o	best.iter <- list (obj = -Inf)								##	empty iter - element, only containing the obj - value
	best.iter <- iter											#n	initializing the best.iter - structure


	for (j in 1:nstart)
	{
		iter <- .InitClusters (x, iter, pa)

		if (!missing (f.hook.iter))
			f.hook.iter (x, .findClustAssig (x, iter, pa), pa, parlist, j, 0)

		lastobj <- -Inf
		for (i in 0:iter.max)
		{
			iter <- f.restr (iter = iter, pa = pa, ...) ##	restricting the clusters' scatter structure

			if (!iter$code)
			{								##	all eigenvalues are zero / a singularity has been detected..
				if (i)
					return (.Parsetclust.Res (x, .TreatSingularity (.calcobj (x, iter, pa), pa), parlist))
				else
					iter$sigma[,,] = diag (pa$p)
            }

			iter <- .findClustAssig (x, iter, pa)				##		finding the cluster assignment on behalf of the current sigma & center information

			if (i && !missing (f.hook.iter))
				f.hook.iter (x, iter, pa, parlist, j, i)

			if (iter$code ||									##		if .findClustAssig returned 1, meaning that the cluster Assignment has not changed
				i == iter.max)									##		or we're in the last concentration step:
				break											##		break the for - loop - we finished this iteration! dont re-estimate cluster parameters this time

			iter <- .estimClustPar (x, iter, pa)					##		estimates the cluster's parameters (cov, center) based on the current iter$assig 

			if (trace >= 2)
			{
				curobj <- .calcobj (x, iter, pa)$obj
				if (curobj < lastobj)
					cat ("obj. function dropped by", lastobj - curobj, "from", lastobj, "to", curobj, "\n")

				lastobj <- curobj
			}
		}

		iter <- .calcobj (x, iter, pa)							##	calculates the obj - value of struct iter
		iter$code = as.numeric (i == iter.max)					##	return code 1 : did not converge
																##				0 : converged

		if (!missing (f.hook.model))
			f.hook.model (x, iter, pa, parlist, j)

		if (j == 1 || iter$obj > best.iter$obj)
			best.iter = iter
	}

	# no more re-calculation of the cluster - assignment necessary, because the algorithm stopped at the right position (after calculating the cluster assignment) 

	return (.Parsetclust.Res (x, best.iter, parlist))
}

####################						##	this function is supposed to act exactly as the final C implementation of tclust
##	.nr_tclust.R  ##
####################
	
.nr_tclust.R <- function (..., restr = c ("eigen", "deter", "prop", "dir.eigen", "dir.deter", "sigma"))
{
	restr <- match.arg (restr, c ("eigen", "deter", "prop", "dir.eigen", "dir.deter", "sigma"))

	if (restr == "eigen")
		return (.tclust.R (f.restr = .restr.diffax, restr.deter = FALSE, ...))
	else if (restr == "deter")
		return (.tclust.R (f.restr = .restr.diffax, restr.deter = TRUE, ...))
	else if (restr == "prop")
		return (.tclust.R (f.restr = .restr.prop, ...))
	else if (restr == "dir.eigen")
		return (.tclust.R (f.restr = .restr.dir, restr.deter = FALSE, ...))
	else if (restr == "dir.deter")	
		return (.tclust.R (f.restr = .restr.dir, restr.deter = TRUE, ...))
	else if (restr == "sigma")
		return (.tclust.R (f.restr = .restr.avgcov, ...))

## execution never comes here. (match.arg checks for unknown/invalid valuesof argument "restr")
}

########################
##	.Parsetclust.Res  ##
########################

.Parsetclust.Res <- function (x, iter, parlist)
{	##	converts the output of function .tclust.R to a "tclust" - object

	idx.clust <- order (iter$csize, decreasing = TRUE)
	idx.nz <- iter$csize[idx.clust] != 0
	idx.clust <- idx.clust [idx.nz]

	id.clust <- 0
	id.clust [1 + idx.clust] <- 1:length (idx.clust)

	int <- list (
		iter.successful = 0,
		iter.converged = 0,
		dim = dim (x)
		)

	ret <- list (
		centers = t (iter$center[idx.clust, , drop = FALSE]), 
		cov = iter$sigma [,, idx.clust, drop = FALSE], 
		cluster = id.clust [iter$assig + 1], 
		par = parlist,
		k = length (idx.clust),
		obj = iter$obj,
		size = iter$csize [idx.clust],
		weights = iter$cw [idx.clust], 
		ret.orig = iter,
		int = int
	)
	class (ret) <- "tclust"
	ret
}
