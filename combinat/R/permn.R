"permn"<-
function(x, fun = NULL, ...)
{
# DATE WRITTEN: 23 Dec 1997          LAST REVISED:  23 Dec 1997
# AUTHOR:  Scott D. Chasalow (Scott.Chasalow@users.pv.wau.nl)
#
# DESCRIPTION:
#             Generates all permutations of the elements of x, in a minimal-
#	change order. If x is a	positive integer,  returns all permutations
#	of the elements of seq(x). If argument "fun" is not null,  applies
#	a function given by the argument to each point. "..." are passed
#	unchanged to the function given by argument fun, if any.
#
#	Returns a list; each component is either a permutation, or the
#	results of applying fun to a permutation.
#
# REFERENCE:
#	Reingold, E.M., Nievergelt, J., Deo, N. (1977) Combinatorial
#	Algorithms: Theory and Practice. NJ: Prentice-Hall. pg. 170.
#
# SEE ALSO:
#	sample, fact, combn, hcube, xsimplex
#
# EXAMPLE:
#	# Convert output to a matrix of dim c(6, 720)
#	t(array(unlist(permn(6)), dim = c(6, gamma(7))))
#
#	# A check that every element occurs the same number of times in each
#	# position
#	apply(t(array(unlist(permn(6)), dim = c(6, gamma(7)))), 2, tabulate, 
#		nbins = 6)
#
#	# Apply, on the fly, the diff function to every permutation
#	t(array(unlist(permn(6, diff)), dim = c(5, gamma(7))))
#
	if(is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x) x <- seq(
			x)
	n <- length(x)
	nofun <- is.null(fun)
	out <- vector("list", gamma(n + 1))
	p <- ip <- seqn <- 1:n
	d <- rep(-1, n)
	d[1] <- 0
	m <- n + 1
	p <- c(m, p, m)
	i <- 1
	use <-  - c(1, n + 2)
	while(m != 1) {
		out[[i]] <- if(nofun) x[p[use]] else fun(x[p[use]], ...)
		i <- i + 1
		m <- n
		chk <- (p[ip + d + 1] > seqn)
		m <- max(seqn[!chk])
		if(m < n)
			d[(m + 1):n] <-  - d[(m + 1):n]
		index1 <- ip[m] + 1
		index2 <- p[index1] <- p[index1 + d[m]]
		p[index1 + d[m]] <- m
		tmp <- ip[index2]
		ip[index2] <- ip[m]
		ip[m] <- tmp
	}
	out
}

