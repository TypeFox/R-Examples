#' @name VarHT
#' @title Variance estimator of Horvitz - Thompson estimator
#' 
#' @description Computes the variance estimator of Horvitz - Thompson estimator of population total
#' 
#' @usage VarHT(y, pikl)
#' @param y A numeric vector of size n containing information about variable of interest
#' @param pikl A square numeric matrix of dimension n containing first and second order inclusion probabilities for units included in \code{y}
#' @details Variance estimator of Horvitz - Thompson estimator of population total is given by
#'  \deqn{\hat{Var}(\hat{Y}_{HT}) = \sum_{k \in s}\frac{y_k^2}{\pi_k^2}(1 - \pi_k) + \sum_{k \in s}\sum_{l \in s, l \neq k} \frac{y_k y_l}{\pi_k \pi_l} \frac{\pi_{kl} - \pi_k \pi_l}{\pi_{kl}}}
#' @return A numeric value representing variance estimator of Horvitz - Thompson estimator for population total for considered values
#' @references Horvitz, D. G. and Thompson, D. J. (1952)
#'  \emph{A generalization of sampling without replacement from a finite universe.}
#'  Journal of the American Statistical Association, 47, 663 - 685
#' @references Sarndal, C. E., Swensson, B. and Wretman, J. (1992)
#'  \emph{Model Assisted Survey Sampling}. Springer-Verlag. New York.
#' @seealso \code{\link{HT}} \code{\link{CovHT}}
#' @examples
#' ##########   Example 1   ##########
#' U <- c(13, 18, 20, 14, 9)
#' #A simple random sample of size 2 without replacement is drawn from population
#' s <- sample(U, 2)
#' #Horvitz - Thompson estimator of population total is calculated.
#' ps <- c(0.4, 0.4)
#' HT(s, ps)
#' #Now, we calculate variance estimator of the Horvitz - Thompson estimator.
#' Ps <- matrix(c(0.4,0.1, 0.1,0.4), 2 ,2)
#' VarHT(s, Ps)
#' 
#' ##########   Example 2   ##########
#' data(DatA)
#' attach(DatA)
#' data(PiklA)
#'
#' #Let calculate Horvitz - Thompson estimator for total of variable Clothing in Frame A.
#' HT(Clo, ProbA)
#' #And now, let compute the variance of the previous estimator
#' VarHT(Clo, PiklA)
#' 
#' @export
VarHT = function (y, pikl)
{
    if (any(is.na(pikl))) 
        stop("There are missing values in pikl.")
    if (any(is.na(y))) 
        stop("There are missing values in y.")
    if (nrow(pikl) != ncol(pikl))
	stop("pikl is not a square matrix.")
    if (length(y) != nrow(pikl)) 
        stop("y and pikl have different sizes.")

    pik<-as.matrix(diag(pikl))
    y <- as.matrix(y)
    delta <- pikl - (pik %*% t(pik))
    ykylexp <- (y/pik) %*% t(y/pik)
    A <- (delta) * (ykylexp)/pikl
    return(sum(A))
}