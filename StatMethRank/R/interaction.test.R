#' Test for Interaction of Ranking Data
#'
#' This function performs a test of interaction of ranking data.
#'
#' @param X a I * J * N array, two-factor design dataset. 
#' X[i, j, n] denotes the response of the n-th replicate in the (i, j) cell .
#' @return a list of test statistic and the p value
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>
#' @examples
#' # Box-cox data revisited
#' # See Box, G. and Cox, D. (1964). An analysis of transformations. 
#' boxcoxdat = array(0, c(3, 4, 4))
#' boxcoxdat[1, , ] = 
#'     matrix(c(
#'         0.31, 0.82, 0.43, 0.45,
#'         0.45, 1.10, 0.45, 0.71,
#'         0.46, 0.88, 0.63, 0.66,
#'         0.43, 0.72, 0.76, 0.62), nrow = 4)
#' boxcoxdat[2, , ] = 
#'     matrix(c(
#'         0.36, 0.92, 0.44, 0.56,
#'         0.29, 0.61, 0.35, 1.02,
#'         0.40, 0.49, 0.31, 0.71,
#'         0.23, 1.24, 0.40, 0.38), nrow = 4)
#' boxcoxdat[3, , ] = 
#'     matrix(c(
#'         0.22, 0.30, 0.23, 0.30,
#'         0.21, 0.37, 0.25, 0.36,
#'         0.18, 0.38, 0.24, 0.31,
#'         0.23, 0.29, 0.22, 0.33), nrow = 4)
#' interaction.test(boxcoxdat)
#' @references A Nonparametric test for interaction in two-way layouts, Xin Gao and Mayer Alvo
interaction.test <- function(X)
{
	# Code by Li Qinglong
	# Last editted: Oct 18, 2013
	# Input:
	# X 	: a I * J * N array, two-factor design dataset. 
	#		Let X[i, j, n] denote the response of the n-th replicate in the (i, j) cell 
    # require(MASS)
    size = dim(X)
    II = size[1] # number of rows
    JJ = size[2] # number of columns
    N = size[3]  # number of levels(sample size)
        
    A = -matrix(1, ncol = II, nrow = II) %x% diag(JJ) / II + 
        diag(II) %x% diag(JJ)
    B = -diag(II) %x% matrix(1, ncol = JJ, nrow = JJ) / JJ + 
        diag(II) %x% diag(JJ)

    Row_rank = array(0, dim = dim(X))
    Col_rank = array(0, dim = dim(X))
    SN = array(0, dim = c(II, JJ)) # sum of the row ranks
    TN = array(0, dim = c(II, JJ)) # sum of the column ranks
    tempmat1 = apply(X, 1, rank)
    tempmat2 = apply(X, 2, rank)
    for (i in 1:N)
    {
        head = (i - 1) * JJ + 1
        tail = head + JJ - 1
        Row_rank[, , i] = t(tempmat1[head:tail, ])
        SN = SN + Row_rank[, , i]

        head = (i - 1) * II + 1
        tail = head + II - 1
        Col_rank[, , i] = tempmat2[head:tail, ]
        TN = TN + Col_rank[, , i]
    }
    SN = as.vector(t(SN)) / (N * JJ + 1) 
    TN = as.vector(t(TN)) / (N * II + 1)
	# Emprical Hajec projection of SN(i, j) onto X(a, b, n)
    Cmat = array(0, dim = c(II, JJ, JJ, N))
    for (i in 1:II)
    {
    	for (j in 1:JJ)
    	{
    		# when b != j
    		for (b in (1:JJ)[-j])
    		{
    			Cmat[i, j, b, ] = -rowSums(matrix(X[i, b, ] >= rep(X[i, j, ], each = N), 
                                                nrow = N))
    		}
    		# when b = j
    		Cmat[i, j, j, ] = rowSums(matrix(X[i, j, ] >= rep(X[i, -j, ], each = N), 
                                             nrow = N))
    	}
    }
    Cmat = Cmat / (N * JJ)
	# Emprical Hajec projection of TN(i, j) onto X(a, b, n)
    Gmat = array(0, dim = c(II, JJ, II, N))
    for (i in 1:II)
    {
    	for (j in 1:JJ)
    	{
    		# when a != i
    		for (a in (1:II)[-i])
    		{
    			Gmat[i, j, a, ] = -rowSums(matrix(X[a, j, ] >= rep(X[i, j, ], each = N), 
                                                nrow = N))
    		}
    		# when a = i
    		Gmat[i, j, i, ] = rowSums(matrix(X[i, j, ] >= rep(X[-i, j, ], each = N), 
                                            nrow = N))
    	}
    }
    Gmat = Gmat / (N * II)
	# Limiting covariance matrices of SN and TN
    sigma1 = array(0, dim = c(II * JJ, II * JJ))
    sigma2 = array(0, dim = c(II * JJ, II * JJ))
    sigma12 = array(0, dim = c(II * JJ, II * JJ))
    for (iRow in 1:(II * JJ))
    {
        i1 = ceiling(iRow / JJ)
        j1 = iRow - (i1 - 1) * JJ
        for (iCol in iRow:(II * JJ))
        {
            i2 = ceiling(iCol / JJ)
            j2 = iCol - (i2 - 1) * JJ
            if (i1 == i2) # a = i1 = i2, otherwise 0
            {
                vec1 = Cmat[i1, j1, , ]
                vec2 = Cmat[i2, j2, , ]                                        
                sigma1[iRow, iCol] = 
                    sum((vec1 - rowMeans(vec1)) * (vec2 - rowMeans(vec2))) / N
                sigma1[iCol, iRow] = sigma1[iRow, iCol]                
            }    
            if (j1 == j2) # b = j1 = j2, otherwise 0
            {
                vec1 = Gmat[i1, j1, , ]
                vec2 = Gmat[i2, j2, , ]
                sigma2[iRow, iCol] = 
                    sum((vec1 - rowMeans(vec1)) * (vec2 - rowMeans(vec2))) / N
                sigma2[iCol, iRow] = sigma2[iRow, iCol]
            }
			# a = i1, b = j2, otherwise 0
            vec1 = Cmat[i1, j1, j2, ]
            vec2 = Gmat[i2, j2, i1, ]
            sigma12[iRow, iCol] = sum((vec1 - mean(vec1)) * (vec2 - mean(vec2))) / N
			sigma12[iCol, iRow] = sigma12[iRow, iCol]
        }        
    }
   
    # Estimated covariance matrix
    Sigma = A %*% sigma1 %*% t(A) +
            B %*% sigma2 %*% t(B) +
            2 * A %*% sigma12 %*% t(B)

    tempvec = A %*% SN + B %*% TN
    # Test statistic, degree of freedom = (I-1)(J-1)
    W = t(tempvec) %*% ginv(Sigma) %*% tempvec / N
	# Under null hypothesis
    p_value = 1 - pchisq(W, df = (II - 1) * (JJ - 1))    
    output = list(Statistic = W, p_value = p_value)
	return(output)
}