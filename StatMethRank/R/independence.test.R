#' Nonparametric Rank Tests for Independence
#'
#' This function performs a nonparametric test of ranking data based 
#' on the correlation. This function can be applied to the ranking data
#' with missing ranks and tie ranks.
#'
#' @param X1 a vector, using NA to stand for the missing ranks
#' @param X2 the same as X1 
#' @param method whether the test is based on Spearman correlation or Kendall
#' correlation
#' @return a list of the test statistics
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>
#' @examples
#' Arith = c(14, 18, 23, 26, 27, 30, 40, NA, NA)
#' Lang = c(28, 14, 46, NA, 53, NA, 54, 50, NA)
#' independence.test(Arith, Lang, method = "spearman")
#' independence.test(Arith, Lang, method = "kendall")
#' @references Rank Correlation Methods for Missing Data, Mayer Alvo and Paul 
#' Cablio \cr
#' Nonparametric Rank Tests for Independence in Opinion Surveys, Philip L.H. Yu, 
#' K.F. Lam, and Mayer Alvo
independence.test <- function(X1, X2, method = c("spearman", "kendall"))
{
	# Test for independence extended to incomplete and tied rankings
	X1 = as.vector(X1)
	X2 = as.vector(X2)
	t = length(X1)
	if (t != length(X2))
	{
		stop("The lengths of X1 and X2 are not equal!")
	}
	k1 = sum(!is.na(X1)) 
	k2 = sum(!is.na(X2))
	
	if (k1 == 0 || k2 == 0) 
	{
		stop("No observations in X1 or X2!")
	}
	# Make sure k1 <= k2
	if (k1 > k2)
	# Swap X1 and X2
	{
		temp = X1
		X1 = X2
		X2 = temp
		temp = k1
		k1 = k2
		k2 = temp
	}

	rank1 = rank(X1, ties.method = "average", na.last = "keep")	
	rank2 = rank(X2, ties.method = "average", na.last = "keep")

	# Using right invariant property
	o = which(!is.na(rank1)) # Label of obejects in ranking 1
	o_star = rank2[o] 
	o_star[is.na(o_star)] = (k2 + 1) / 2
	o_mean = mean(o_star)

	ind = !(is.na(X1 + X2))
	k_star = sum(ind)
	
	method = match.arg(method)
    if (method == "spearman")
    {
		As_star = (t + 1)^2 / ((k1 + 1) * (k2 + 1)) * 
			sum((rank1[ind] - (k1 + 1) / 2) * (rank2[ind] - (k2 + 1) / 2))

		# Under H1
		varAs_star = ((t + 1) ^ 2 / ((k1 + 1) * (k2 + 1)))^2 * (1 / (k1 - 1)) * 
		    sum((o_star - o_mean)^2) * sum((rank1[o] - (k1 + 1) / 2)^2)

		# Under H2
		K1 = k1 * (k1 - 1) / (k1 + 1)
		K2 = k2 * (k2 - 1) / (k2 + 1)
		# varAs_star = (t + 1)^4 / (144 * (t - 1)) * K1 * K2

		Cs = t * (t^2 - 1) / 12
		# Genaralized Spearman distance
		ds_star = Cs - As_star
		if (k1 %% 2) # if k1 is odd
		{
			rs = (k1^2 - 1) * (3 * k2 - k1)		
		}
		else # if k1 is even
		{
			rs = k1 * (k1 * (3 * k2 - k1) - 2)
		}
		rs = rs * (t + 1) ^ 2 / (24 * (k1 + 1) * (k2 + 1))
		# ms = Cs - rs
		# Ms = Cs + rs

		# Correlation computation is only right when there is no tied rankings
		# Type a correlation
		Alpha_a = As_star / rs
		# Type b correlation
		Alpha_b = 12 * As_star / ((t + 1)^2 * sqrt(K1 * K2))
		# Standarlized test statistic
		z_stat = As_star / sqrt(varAs_star)
		# Double tail test
		p_value = 2 * (1 - pnorm(abs(z_stat)))
		
		res = list(	Similarity_Measure = As_star, 
					# Spearman_Correlation = Alpha_a, 
					Test_Statistic = z_stat, 
					P_value = p_value,
					Distance = ds_star
					)	
    }
    else if (method == "kendall")
    {
    	a1 = matrix(rep(0, t^2), ncol = t)
		a2 = matrix(rep(0, t^2), ncol = t)
		for (i in 1:(t-1))
		{
			for (j in (i+1):t)
			{
				a1[i, j] =  sum(c(sign(rank1[i] - rank1[j]) * !is.na(rank1[i]) * !is.na(rank1[j]),
					is.na(rank1[j]) * !is.na(rank1[i]) * (1 - 2 * rank1[i] / (k1 + 1)),
					is.na(rank1[i]) * !is.na(rank1[j]) * (2 * rank1[j] / (k1 + 1) - 1)), 
					na.rm = TRUE)

				a2[i, j] =  sum(c(sign(rank2[i] - rank2[j]) * !is.na(rank2[i]) * !is.na(rank2[j]),
					is.na(rank2[j]) * !is.na(rank2[i]) * (1 - 2 * rank2[i] / (k2 + 1)),
					is.na(rank2[i]) * !is.na(rank2[j]) * (2 * rank2[j] / (k2 + 1) - 1)), 
					na.rm = TRUE)
			}
		}
		Ak_star = sum(a1 * a2)

		# Under H1 (16 / t^2) * varAs_star
		varAk_star = (16 / t^2) * ((t + 1) ^ 2 / ((k1 + 1) * (k2 + 1)))^2 * 
			(1 / (k1 - 1)) * sum((o_star - o_mean)^2) * sum((rank1[o] - (k1 + 1) / 2)^2)

		# Under H2
		K1 = k1 * (k1 - 1) / (k1 + 1)
		K2 = k2 * (k2 - 1) / (k2 + 1)
		# K1 = k1 * (k1 - 1) / (k1 + 1) * (1 - sum(g1^3 - g1) / (k1^3 - k1))
		# K2 = k2 * (k2 - 1) / (k2 + 1) * (1 - sum(g2^3 - g2) / (k2^3 - k2))
		# varAk_star = K1 * K2 / (9 * t * (t - 1)) * 
		#	((2 * t + k1 + 3) * (2 * t + k2 + 3) / 2 + (t^2 - k1 - 2) * (t^2 - k2 - 2) / (t - 2))

		Ck = t * (t - 1) / 2
		# Genaralized Spearman distance
		dk_star = Ck - Ak_star
		if (k1 %% 2) # if k1 is odd
		{
			rk = (k1^2 - 1) * (t * (3 * k2 - k1) + k2 * (k1 + 3))		
		}
		else # if k1 is even
		{
			rk = k1 * (3 * k1 * k2 * (t + 1) - (k1^2 + 2) * (t - k2) - 3 * (k2 + 1))
		}
		rk = rk / (6 * (k1 + 1) * (k2 + 1))
		# ms = Cs - rs
		# Ms = Cs + rs

		# Correlation computation is only right when there is no tied rankings
		# Type a correlation
		Alpha_a = Ak_star / rk
		# Type b correlation
		Alpha_b = 6 * Ak_star / sqrt((2 * t + k1 + 3) * (2 * t + k2 + 3) * K1 * K2)
		# Standarlized test statistic
		z_stat = Ak_star / sqrt(varAk_star)
		# Double tail test
		p_value = 2 * (1 - pnorm(abs(z_stat)))
		
		res = list(	Similarity_Measure = Ak_star, 
					# Kendall_Correlation = Alpha_a, 
					Test_Statistic = z_stat, 
					P_value = p_value,
					Distance = dk_star
					)	
    }
	return(res)
}