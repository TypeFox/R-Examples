
.restr.none <- function (iter, pa)
{
	iter$code = TRUE
	return (iter)
}

##	restricts the clusters' covariance matrices by averaging them
.restr.avgcov <- function (iter, pa)
{
	s.all <- matrix (0, pa$p, pa$p)
	for (k in 1:pa$K)
		s.all <- s.all + iter$sigma[,,k] * iter$csize[k]/sum(iter$csize)

	iter$sigma[,,] <- s.all			##

	iter$code = sum(diag (s.all)) > pa$zero.tol
	return (iter)
}

##	restricts the clusters' covariance matrices by restricting it's eigenvalues or determinants without restrics directions 
.restr.diffax <- function (iter, pa, f.restr.eigen = .restr2_eigenv, restr.deter, restr.fact = 12)
{
	if (pa$p == 1)								#	one - dimensional
		if (restr.fact == 1)					#	all variances should be equal -> use the simpler .restr.avgcov function instead
			return (.restr.avgcov (iter, pa))	#
		else									#
			restr.deter <- FALSE				#	else - if p == 1 always use the simpler eigen - restriction 

	if (!missing (restr.deter))					#	evaluating the appropriate restriction function
		f.restr.eigen <- if (restr.deter) .restr2_deter_ else .restr2_eigenv								#n

	u <- array (NA, c(pa$p, pa$p, pa$K))
	d <- array (NA, c(pa$p, pa$K))
	for (k in 1:pa$K)
	{
		ev <- eigen (iter$sigma[,,k])
		u [,,k] <- ev$vectors
		d [,k] <- ev$values
	}

	d [d < 0] <- 0		##	all eigenvalue < 0 are restricted to 0 BECAUSE: min (d) < 0 -> max (d) / min (d) < 0 which would always fit our eigenvalue - criterion!
						##	there is *NO* reason *AT ALL* for touching *ANY* eigenvalue > 0!! (e.g. 0 <= d < 1e-16!)

	d <- f.restr.eigen (d, iter$csize, restr.fact, pa$zero.tol)

						##	checking for singularity in all clusters.
	iter$code = max(d) > pa$zero.tol

	if (!iter$code)
		return (iter)

	for (k in 1:pa$K)	##	re-composing the sigmas
		iter$sigma[,,k] <- u[,,k] %*% diag (d[,k], nrow = pa$p) %*% t(u[,,k])
	return (iter)
}

.diag_matmult_y_x_yt <- function (x, y) diag (y %*% x %*% t(y))			#	a simple helper function for use with "apply"
.sum_diag_matmult_y_x <- function (x, y) sum (diag (y %*% x))			#	a simple helper function for use with "apply"

# Function for restricted principal components which allows restrictions type 0 and 1
.restr.dir <- function (iter, pa, restr.deter = FALSE, restr.fact = 12, n_inic = 10, n_iter1 = 10, n_iter2 = 10, ...)
{
	if (pa$p == 1 ||								# only one dim -> use the simpler .restr.diffax function instead, as directions make no sense here.
		(restr.fact == 1 && !restr.deter))			# all eigenvalues should be equal -> use the simpler .restr.diffax function instead!
		return (.restr.diffax (iter, pa, restr.fact = restr.fact))

													#	evaluating the appropriate restriction function
	f.restr.eigen <- if (restr.deter) .restr2_deter_ else .restr2_eigenv

	for (inic in 1:n_inic)
	{
		#initial solution for u
#o		ee <- matrix(runif(pa$p*pa$p, min=0, max=1),pa$p,pa$p)
		ee	<- matrix (rnorm (pa$p * (pa$p+1)), ncol=pa$p)					#n	as suggested by agustin 20110217
		w <- eigen(t(ee)%*%ee)
		u <- w$vectors

		#iterations
		for (iter1 in 1: n_iter1)
		{
			dd <- apply (iter$sigma, 3, .diag_matmult_y_x_yt , y = t(u))

			dd <- f.restr.eigen (dd, iter$csize, restr.fact, pa$zero.tol)

			u <- .optvectors (iter, pa, u, dd, n.iter = n_iter2, ...)
		}

		res <- 0

		for (i in 1:pa$K)
			res <- res + iter$csize[i] * sum (diag (u %*% diag (1 / dd[,i]) %*% t (u) %*% iter$sigma[,,i]))

										##	no need to store all solutions; only saving the best solution

		if (inic == 1 || res < best$res)	#n	store the first/best solution
			best <- list (u = u, dd = dd, res = res)
	}

	## reconstructing sgima
										
	iter$code <- max(best$dd) > pa$zero.tol

	if (!iter$code)						#	if function failed
		return (iter)
										#	calculate the new sigmas
	for (i in 1:pa$K)					#	using the best found solution
		iter$sigma[,, i] <- best$u %*% diag (best$dd[, i]) %*% t (best$u)

	return(iter)
}

# Function for vector optimization
.optvectors <- function (iter, pa, u, d, n.iter, ovv = 0)
{
	#iterations in each pair of directions

	for (iter2 in 1:n.iter)
	{
#t		for (j in 1:pa$p)
		for (j in 1:(pa$p - 1))
		{
#t			for (m in 1: pa$p)
			for (m in (j+1): pa$p)
			{
				a_mj <- matrix(0, pa$p, pa$p)
				for (i in 1:pa$K)
					a_mj <- a_mj + (d[j,i]   - d[m,i]  ) / (d[j,i]   * d[m,i]  ) * iter$csize[i] * iter$sigma[,,i]		## NOTE - this is actually equal for each iteration of iter2. Thus these values could be calculated in a 

				a2_mj <- t (u[,c(j,m)]) %*% a_mj %*% u[,c(j,m)]
				v <- eigen(a2_mj)
				uu <- u[,c(j,m)] %*% v$vectors

#cat ("values:\n")
#print (v$values)
#cat ("vecs:\n")
#print (v$vectors)
##browser ()

				if (ovv == 0)
				{
##	### NEW:	new part without checking for labels and signs.
					u[,j] <- uu[,1]
					u[,m] <- uu[,2]
				}
				else if (ovv == 1)
				{
##	### 2DELETE:	old part which checked the labels and signs
					#checking labels and signs in order to allow identifiability

					up1 <- up2 <- u
					up1[,j] <- up2[,m] <- uu[,1]
					up1[,m] <- up2[,j] <- uu[,2]


					resp1 <- resp2 <- 0
					for (i in 1:pa$K)
					{
						resp1 <- resp1 + iter$csize[i] * sum (diag (up1 %*% diag (1 / d[,i]) %*% t (up1) %*% iter$sigma[,,i]))	#n	extracted iter$csize[i] out of sum (diag (...))
						resp2 <- resp2 + iter$csize[i] * sum (diag (up2 %*% diag (1 / d[,i]) %*% t (up2) %*% iter$sigma[,,i]))	#n
					}

					less <- (resp1 - resp2) <= sqrt (pa$zero.tol)

					neg1 <- t (uu[, 1 + !less]) %*% u[, j] < 0		#	transforming the logical variable ("less") into the appropriate index
					neg2 <- t (uu[, 1 +  less]) %*% u[, m] < 0		#	

					u <- if (less) up1 else up2						#	storing the proper up1 or up2 matrix in u

					if (neg1) u[, j] <- -u[, j]						#	changing the sign of the proper columns (if necessary)
					if (neg2) u[, m] <- -u[, m]
				}
#	print (u)
#	cat ("\n")
			}
		}
	}
	return(u)
}

# Function for proporcionality case (only restrictions type 1 are allowed )
.restr.prop <- function (iter, pa, restr.fact = 12, n_inic = 10, n_iter1 = 10, n_iter2 = 10, ...)
{
	if (pa$p == 1)							## only one dim -> use the simpler .restr.diffax function instead!
		return (.restr.diffax (iter, pa, restr.fact = restr.fact))
	if (restr.fact == 1)					## all covmats should be exactly equal -> use the simpler .restr.avgcov function instead!
		return (.restr.avgcov (iter, pa))

	for (inic in 1:n_inic)
	{
		#initial solution
#o		ee	<- matrix (runif (pa$p * pa$p, min = 0, max = 1), pa$p, pa$p)
		ee	<- matrix (rnorm (pa$p * (pa$p+1)), ncol=pa$p)					#n	as suggested by agustin 20110217
		w	<- eigen (t (ee) %*% ee)
		u	<- w$vectors
		sf	<- w$values / (prod (w$values)^(1/pa$p))
				##	remark: if w$values[pa$p] is very small we might run in numerical problems here!

		#iterations
		for (iter1 in 1: n_iter1)
		{
			dd <- apply (iter$sigma, 3, .diag_matmult_y_x_yt , y = t(u))
			st <- colSums (dd / sf) / pa$p

			st <- .restr2_eigenv (t (st), iter$csize, restr.fact^(1 / pa$p), pa$zero.tol)

			sf <- colSums	(t (dd) * iter$csize / as.numeric (st))

			sf <- sf / prod(sf)^(1 / pa$p)

			#u estimation

			dd <- sf %*% st							#	st is a 1 x K matrix (from .restr2_eigenv)
			u <- .optvectors(iter, pa, u, dd, n.iter = n_iter2, ...)
		}

		tempSum <- apply (iter$sigma, 3, .sum_diag_matmult_y_x, y = u %*% diag (1/sf) %*% t (u))
		res <- sum (iter$csize * (log (st^pa$p) + tempSum / st))

		if (inic == 1 || res < best$res)
			best <- list (u = u, sf = sf, st = st, res = res)
	}

	iter$code <- max(best$sf) > pa$zero.tol

	if (!iter$code)								#	if function failed
		return (iter)
												# reconstruct the sigma - matrices
	for (i in 1:pa$K)
		iter$sigma[,,i] <- best$st[i] * best$u %*% diag (best$sf) %*% t (best$u)

#a	iter$sigma <- sapply (best$st, get ("*"), sigmaBase, simplify = FALSE) - iter$sigma	#a alternative: think about future representation of sigma as a list!

	return(iter)
}
