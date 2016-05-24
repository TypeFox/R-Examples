#' Generate all possible permutations of k elements out of n 
#'
#' This function generates all permutations of n elements taken k at a time.
#' The most efficient way to generate n! permutations of vec 
#' in a Lexicographic order without recursive algorithm in R.
#'
#' @param n number of the whole elements
#' @param k number of elements to permute (default the same as n)
#' @param vec the source vector of length n (default as c(1:n))
#' @return a matrix or a vector, each row is a permutation
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>
#' @examples
#' generate.perms(10, 6)
#' @seealso generate.combn
generate.perms <- function (n, k = n, vec = 1:n)
{
	# Generate all possilbe permutations of k elements out of n 
	# An efficient way to generate n! permutations of vec without recursive 
	# algorithm in a Lexicographic order
	# Author 	: Li Qinglong
	# Email 	: liqinglong0830@163.com
	# Created   : Oct 29, 2013
	if (k > n)
	{
		stop("k should be smaller than or equal to n!")
	}
	if (length(vec) != n)
	{
		stop("The vector should contain exactly n elements!")
	}
	if (n > 1)
	{
		permns <-function(n, vec = 1:n)
		{
			# Generate all possible permutaions of vector contains n elements
			ans = matrix(1, nrow = factorial(n), ncol = n)
			faci = 1
			if (n > 1)
			{
				for (i in 2:n)
				{   
					faci_1 = faci   # i.e. factorail(i-1)
					faci = faci_1 * i # i.e. factorial(i)
					series = 1:i
					head = faci - faci_1 + 1
					for (p in i:1)
					{
						rvec = series[-p]
						ans[head:(head + faci_1 - 1), 2:i] =
							rvec[ans[1:faci_1, 1:(i-1)]]
						head = head - faci_1
					}  
					ans[1:faci, 1] = rep(1:i, each = faci_1)
				}
			}
			ans = matrix(vec[ans], ncol = n)
			return(ans)
		}
		mat = t(as.matrix(combn(vec, k)))
	    ans = apply(mat, 1, permns, n = k)
		ans = matrix(ans, ncol = k)
		return(ans)
	} else return(as.matrix(vec, nrow = 1))
}
