#' @name HT
#' @aliases HT
#' @title Horvitz - Thompson estimator
#' 
#' @description Computes the Horvitz - Thompson estimator 
#' 
#' @usage HT(y, pik)
#' @param y A numeric vector of size n containing information about variable of interest
#' @param pik A numeric vector of size n containing first order inclusion probabilities for units included in \code{y}
#' @details Horvitz - Thompson estimator of population total is given by
#'  \deqn{\hat{Y}_{HT} = \sum_{k \in s} \frac{y_k}{\pi_k}}
#' @return A numeric value representing Horvitz - Thompson estimator for population total for considered values
#' @references Horvitz, D. G. and Thompson, D. J. (1952)
#'  \emph{A generalization of sampling without replacement from a finite universe.}
#'  Journal of the American Statistical Association, 47, 663 - 685
#' @examples
#' ##########   Example 1   ##########
#' U <- c(13, 18, 20, 14, 9)
#' #A simple random sample of size 2 without replacement is drawn from population
#' s <- sample(U, 2)
#' ps <- c(0.4, 0.4)
#' HT(s, ps)
#' 
#' ##########   Example 2   ##########
#' data(DatA)
#' attach(DatA)
#' #Let estimate population total for variable Feeding in frame A
#' HT(Feed, ProbA)
#' @seealso \code{\link{VarHT}}
#' @export
HT = function (y, pik)
{	
    if (any(is.na(pik))) 
        stop("There are missing values in pik.")
    if (any(is.na(y))) 
        stop("There are missing values in y.")
    if (length(y) != length(pik)) 
        stop("y and pik have different sizes")

    return (as.numeric(crossprod(y, 1/pik)))
}