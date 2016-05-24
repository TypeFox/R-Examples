#' Is compatible with the candidate_ranking
#'
#' This function judges whether the complete rangkings are compatible with 
#' an incomplete ranking.
#'
#' @param complete_ranking a matrix or a vector, each column is a complete ranking.
#' @param candidate_ranking a vector, using NA to stand for the missing ranks.
#' @return a vector of TRUEs and FALSEs 
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>
#' @examples
#' u_star = c(2, NA, 3, 4, 1)
#' C_set = compatible.rankings(u_star)
#' is.compatible(C_set, u_star)
is.compatible <- function(complete_ranking, candidate_ranking)
{
	# Whether the complete rangkings are compatible with an incomplete
	# ranking of k of these subjects
	# Input:
	# complete_ranking: a matrix or a vector, each column is a complete
	# ranking.
	# incomplete_ranking: a vector, using NA to stand for the missing ranks.
	na_ind = is.na(candidate_ranking)

	y = rank(candidate_ranking[!na_ind], ties.method = "average")

	x = apply(matrix(complete_ranking[!na_ind], nrow = sum(!na_ind)), 
			2, rank)
	# Inverse index, i.e. for any X, sort(X)[inv_ind] = X
	inv_ind = apply(x, 2, rank, ties.method = "random")
	sorted_ranking = sort(y)
	x = matrix(sorted_ranking[inv_ind], nrow = length(y))
	# return(apply(x, 2, identical, y))
	# more efficient way
	return(colSums(abs(x - y)) == 0)
}
