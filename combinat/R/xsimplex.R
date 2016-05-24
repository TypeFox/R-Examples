"xsimplex"<-
function(p, n, fun = NULL, simplify = TRUE, ...)
{
#       DATE WRITTEN: 11 February 1992          LAST REVISED:  10 July 1995
#       AUTHOR:  Scott Chasalow
#
#       DESCRIPTION:
#             Generates all points on a {p,n} simplex lattice (i.e. a p-part 
#             composition of n).  Each point is represented as x, a 
#             p-dimensional vector of nonnegative integers that sum to n.
#             If argument "fun" is not null,  applies a function given
#             by the argument to each point.  If simplify is FALSE,  returns 
#             a list; else returns a vector or an array.  "..." are passed 
#             unchanged to function given by argument fun,  if any.
#       EXAMPLE:
#             Compute Multinomial(n = 4, pi = rep(1/3, 3)) p.f.:
#             xsimplex(3, 4, dmnom, prob=1/3) 
#
	if(p < 1 || n < 0) return(if(simplify) numeric(0) else list())
	p1 <- p - 1
	x <- numeric(p)
	x[1] <- n
	nofun <- is.null(fun)
	out <- if(nofun) x else fun(x, ...)
	if(p == 1 || n == 0) {
		return(if(simplify) out else list(out))
	}
	count <- nCm(n + p - 1, n)
	if(simplify) {
		dim.use <- NULL
		if(nofun) {
			if(count > 1)
				dim.use <- c(p, count)
		}
		else {
			d <- dim(out)
			if(count > 1) {
				if(length(d) > 1)
				  dim.use <- c(d, count)
				else if(length(out) > 1)
				  dim.use <- c(length(out), count)
			}
			else if(length(d) > 1)
				dim.use <- d
		}
	}
	out <- vector("list", count)
	target <- 1
	i <- 0
	while(1) {
		i <- i + 1
		out[[i]] <- if(nofun) x else fun(x, ...)
		x[target] <- x[target] - 1
		if(target < p1) {
			target <- target + 1
			x[target] <- 1 + x[p]
			x[p] <- 0
		}
		else {
			x[p] <- x[p] + 1
			while(x[target] == 0) {
				target <- target - 1
				if(target == 0) {
				  i <- i + 1
				  out[[i]] <- if(nofun) x else fun(x, ...)
				  if(simplify) {
				    if(is.null(dim.use))
				      out <- unlist(out)
				    else out <- array(unlist(out), dim.use)
				  }
				  return(out)
				}
			}
		}
	}
}

