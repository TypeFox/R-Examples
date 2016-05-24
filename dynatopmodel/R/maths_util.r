#-------------------------------------------------------------------------------
#  maths utilities including
#  * numerical solution methods
#  * matrix operations
#-------------------------------------------------------------------------------
# given initial conditions x0 at t=t0, obtain approximate solution using Euler's
# method at t=t1 to ODE dx/dt=func(x,t)
solve.euler <- function(x0, t1, func, t0=0, nstep=10, return.all=F, implicit=F, ...)
{
	if(implicit)
	{
		return(# backwards implicit euler scheme using fixed point iteration
			solve.euler.implicit(x0, t1, func, t0, nstep, return.all=return.all, ...)
		)
	}
	nx <- length(x0)
	x <- matrix(NA,nrow=nstep, ncol=nx, byrow=T)
	x[1,]<- x0
	dt <- (t1-t0)/nstep
	for(i in 2:nstep)
	{
		t <- t0 + dt*(i-1)
		x[i,] <- x[i-1,] + dt*func(x[i-1,], t=t, ...)
	}
	if(return.all)
	{
		return(cbind(seq(t0, t1, length.out=nstep), x))
	}
	return(x[nstep,])
}

# backwards implicit euler scheme using fixed point iteration
solve.euler.implicit <- function(q0, t1, dqdt, t0=0, nstep=10, niter=30, eps=1e-6, return.all=F, ...)
{
	nq <- length(q0)
	# qk+1=qk + dt * q'(tk+1, qk+1)
	# first approximation
	q <- q0
	dt <- (t1-t0)/nstep
	t <- t0
	for(istep in 0:nstep)
	{
		unsolved <- 1:nq
		i <- 1
		while(i <= niter & length(unsolved)>0)
		{
			prev <- q[unsolved]
			# qk+1 = qk + h *f(tk+1, qk+1)
			#
			dq.iter <- dt*dqdt(t=t+dt*(istep+1), q=q, ...)
			q[unsolved] <- q[unsolved] + dq.iter[unsolved]
			dq <- ifelse(prev==q[unsolved], 0, abs(prev-q[unsolved])/prev)
			unsolved <- unsolved[-which(dq < eps)]
			i <- i+1
		}
		if(length(unsolved)>0)
		{
			warning(paste0("Non-convergence in implicit solution, elements ", unsolved, " of ", nq, ". Try increased number of iterations (currently ", niter, ")."))
		}
	}

	return(q)

}

#------------------------------------------------------------------------------#
identity.matrix <- function(n)
{
  mat <- matrix(ncol=n,nrow=n)
  mat[] <- 0
  mat[col(mat)==row(mat)]<-1
  return(mat)
}

# Ron E. VanNimwegen wrote:
#   >
#   > I was looking for a similar operator, because R uses scalar products
# > when raising a matrix to a power with "^".  There might be something
# > more elegant, but this little loop function will do what you need for a
# > matrix "mat" raised to a power "pow":
#   >
#   > mp <- function(mat,pow){
#     > ans <- mat
#     > for ( i in 1:(pow-1)){
#       > ans <- mat%*%ans
#       > }
#     > return(ans)
#     > }
# >
#   This function is extremely inefficient for high powers of pow
# [an unhappy choice of variables, because pow is the power
#  function in many languages]
#
# A better method would keep track of the powers of two,
# and would optimize according to it.
#
# For example, in order to compute mat^42, we just have to
# compute mat^2, mat^4, mat^8, mat^16, mat^32, and
# then mat^40 and mat^42, with "only" 7 multiplications, instead
# of 41.
#
# See the technique in...
# http://en.wikipedia.org/wiki/Exponentiation_by_squaring
matrix.power <- function(mat, n)
{
  # test if mat is a square matrix
  # treat n < 0 and n = 0 -- this is left as an exercise
  # trap non-integer n and return an error
  if (n == 1) return(mat)
  result <- diag(1, ncol(mat))
  while (n > 0) {
    if (n %% 2 != 0) {
      result <- result %*% mat
      n <- n - 1
    }
    mat <- mat %*% mat
    n <- n / 2
  }
  return(result)
}

# ensure sum of all rows in matrix is unity
normalise.rows <- function(w)
{
  if(!is.matrix(w)){warning("matrix input expected to normalise.rows"); return(w)}
  if(any(rowSums(w)!=1))
  {
    cat("Normalising rows...\n")

    #     w <-apply(w, MARGIN=1,
    #               FUN=function(x)
    #               {
    #                 sum <- sum(x, na.rm=T)
    #                 res<- ifelse(sum==0, rep(0, length(x)), x/sum)
    #                 return(res)
    #               }
    #     )


    w <- apply(w, MARGIN=1,
               FUN=function(x)
               {
                 sum <- sum(x, na.rm=T)
                 res<- x/sum
                 return(res)
               })
    w[!is.finite(w)]  <- 0
    #  w[1:nrow(w),]<-w[1:nrow(w),]/rowSums(w)
    # transpose to get in rightg row. col order
    w <- t(w)
  }
  # need to transpose to get in right row-col order
  return(w)
}

# graph utilities
# produice an adjacency matrix for the graph (directed=F) or digraph (T)
# i, j element is 1 if vertex i connects to vertex j
# see Foulds (1992) Chapter 6 p.76-78
AdjacencyMatrix <- function(conns,
                            n=max(conns), # max vertexes considered,defaults to max encountered in connection matrix
                            directed=T)
{
  res <- matrix(ncol=n, nrow=n)
  res[]<-0

  # conns a list of from -> to vertex numbers
  ir <-1
  while(ir < nrow(conns))
  {
    from <- conns[ir,1]
    to <- conns[ir,2]
    res[from, to]<-1
    ir <- ir + 1
  }
  # if undirected make symmetric...
  return(res)
}
