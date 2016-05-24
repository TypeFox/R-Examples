#' @name CovHT
#' @title Covariance estimator between two Horvitz - Thompson estimators
#' 
#' @description Computes the covariance estimator between two Horvitz - Thompson estimators of population total from survey data obtained from a single stage sampling design
#' 
#' @usage CovHT(y, x, pikl)
#' @param y A numeric vector of size n containing information about first variable of interest in the sample
#' @param x A numeric vector of size n containing information about second variable of interest in the sample
#' @param pikl A square numeric matrix of dimension n containing first and second order inclusion probabilities for units included in the sample
#' @details Covariance estimator between two Horvitz - Thompson estimators of population total is given by
#'  \deqn{\hat{Cov}(\hat{Y}_{HT}, \hat{X}_{HT}) = \sum_{k \in s}\sum_{l \in s} \frac{\pi_{kl} - \pi_k \pi_l}{\pi_{kl}}\frac{y_k}{\pi_k}\frac{x_l}{\pi_l}}
#' @return A numeric value representing covariance estimator between two Horvitz - Thompson estimators for population total for considered values
#' @references Horvitz, D. G. and Thompson, D. J. (1952)
#'  \emph{A generalization of sampling without replacement from a finite universe.}
#'  Journal of the American Statistical Association, 47, 663 - 685
#'  @references Sarndal, C. E., Swensson, B. and Wretman, J. (1992)
#'  \emph{Model Assisted Survey Sampling}. Springer-Verlag. New York.
#' @seealso \code{\link{HT}} \code{\link{VarHT}}
#' @examples
#' ##########   Example 1   ##########
#' Indicators <- c(1, 2, 3, 4, 5)
#' X <- c(13, 18, 20, 14, 9)
#' Y <- c(2, 0.5, 1.2, 3.3, 2)
#' #Let draw two simple random samples without replacement of size 2
#' s <- sample(Indicators, 2)
#' sX <- X[s]
#' sY <- Y[s]
#' #Now, let calculate the associated probability matrix with first and
#' #second order inclusion probabilities
#' Ps <- matrix(c(0.4,0.2, 0.2,0.4), 2, 2)
#' CovHT(sX, sY, Ps)
#' 
#' ##########   Example 2   ##########
#' data(DatA)
#' attach(DatA)
#' data(PiklA)
#' #Let calculate Horvitz - Thompson estimator for total of variable Clothing in Frame A.
#' HT(Clo, ProbA)
#' #Let calculate Horvitz - Thompson estimator for total of variable Feeding in Frame A.
#' HT(Feed, ProbA)
#' #And now, let compute the covariance between the previous estimators
#' CovHT(Clo, Feed, PiklA)
#' @export
CovHT = function (y, x, pikl)
{
	if (any(is.na(y)))
		stop("There are missing values in y.")
	if (any(is.na(x)))
		stop("There are missing values in x.")
	if (any(is.na(pikl)))
		stop("There are missing values in pikl.")
	if (length(y) != length (x))
		stop("x and y have different sizes.")

  pik<-as.matrix(diag(pikl))
  y <- as.matrix(y)
	x <- as.matrix(x)
  delta <- pikl - (pik %*% t(pik))
  ykxlexp <- (y/pik) %*% t(x/pik)
  A <- (delta) * (ykxlexp)/pikl
  return(sum(A))
}