#' Set of compatible rankings
#'
#' This function generates the set of complete permutations compatible
#' with the input_ranking.
#'
#' @param input_ranking a vector, using NAs to stand for the missing ranks
#' @return a matrix or a vector, each column is a complete ranking.
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>
#' @examples
#' u_star = c(2, NA, 3, 4, 1)
#' C_set = compatible.rankings(u_star)
#' is.compatible(C_set, u_star)
compatible.rankings <- function(input_ranking)
{
	# Generate the set of complete permutations compatible
	# with the incomplete_ranking.
	# Using NA to stand for the missing ranks.
	# Input:
	# input_ranking: vector, using NA to stand for the missing ranks.
	# Output:
	# A matrix contains all complete permutations compatible
	# with the incomplete_ranking. Each column is a ranking.
	input_ranking = rank(as.vector(input_ranking), 
		ties.method = "average", na.last = "keep")
	n = length(input_ranking)
	k = sum(!is.na(input_ranking))
	if (k < 2)
	{
		stop("k should be as least 2(k >= 2)")
	}
	na_ind = is.na(input_ranking)

	# Begin to deal with ties
	relative_ranking = input_ranking[!na_ind]
	inv_ind = rank(relative_ranking, ties.method = "random")
	unique_ranks = sort(unique(relative_ranking))
	nGroups = length(unique_ranks)
	sorted_ranking = sort(relative_ranking)
	# Tie pattern
	g = rowSums(matrix(rep(sorted_ranking, nGroups), nrow = nGroups, byrow = T) == unique_ranks)
	# Successive ranks who share the same rank in each group	
	v_list = mapply(seq, to = cumsum(g), length.out = g)
	# All possible permutations in each group
	group_perm_list = mapply(generate.perms, n = g, k = g, vec = v_list)
	# I don't how to explain the next few lines in the comments...hope you can understand
	Mleft = c(1, cumprod(factorial(g[-nGroups])))
	Mright = rev(c(1, cumprod(factorial(g[nGroups:2]))))
	# Using kronecker products to replicate matrix
	construt_mat <- function(A, M1, M2)  return(matrix(1, M1, 1) %x% A %x% matrix(1, M2, 1))
	group_perm_list = mapply(construt_mat, A = group_perm_list, M1 = Mright, M2 = Mleft)
	# All possible relative rankings  
	relative_mat = matrix(unlist(group_perm_list), ncol = k)
	relative_mat = matrix(relative_mat[, inv_ind], ncol = k)
	# End of dealing with ties

	# Begin to deal with NAs
	if (k < n) # incomplete ranking case
	{
		# Permutate the missing ranks
		na_perm = generate.perms(n, n - k)
		# Remaining ranks of 1:n when n-k missing positions are ranked
		remain = t(apply(na_perm, 1, setdiff, x = 1:n))
		nrold = dim(relative_mat)[1]
		nrnew = dim(na_perm)[1]
		na_perm = repmat(na_perm, nrold, 1)
		remain = repmat(remain, nrold, 1)
		# Using kronecker product to do matrix replication
		# Relative ranks of remain, ie the subscript index of remain, that's why I
		# plus k * seq(0, nrnew * nrold - 1) as remain is a (nrnew*nrold) * n matrix
		relative_mat = relative_mat %x% matrix(1, nrnew, 1) +
			k * seq(0, nrnew * nrold - 1)

	    output_mat = matrix(0, ncol = n, nrow = nrold * nrnew)
	    output_mat[, na_ind] = na_perm
	    output_mat[, !na_ind] = matrix(t(remain)[c(relative_mat)], ncol = k)
	} else # complete ranking case
	{
		output_mat = relative_mat
	}
	# End of dealing with NAs
	return(t(output_mat)) 
}