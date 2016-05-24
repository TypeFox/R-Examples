#' Generate all incomplete rankings of k elements out of n complete elements
#'
#' This function generates all the possible incomplete rankings of 
#' k elements out of n complete elements. E.g. c(1, 2, NA) is an 
#' incomplete ranking out with 2 non-missing elements of 3 elements
#'
#' @param n number of complete elements
#' @param k number of non-missing elements in the incomplete ranking
#' @return a matrix each of whose column is a possible incomplete ranking
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>
#' @examples
#' incomplete.rankings(5, 3) 
incomplete.rankings <- function(n, k)
{
	if (n < k)
	{
		stop("n should be larger than k!")
	}

	output_mat = matrix(NA, ncol = factorial(n) / factorial(n - k), nrow = n)
	u = t(generate.perms(k))
	ind = generate.combs(n, k)
	len = dim(u)[2]
	head = 1
	for (i in seq(dim(ind)[2]))
	{
		tail = head + len - 1
		output_mat[ind[, i], head:tail] = u
		head = tail + 1
	}
	return(output_mat)
}