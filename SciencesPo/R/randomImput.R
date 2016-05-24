#' @encoding UTF-8
#' @title Simple Random Imputation
#'
#' @description Performs random imputation in a vector containing missing values.
#'
#' @param x a vector whose missing values (\code{NA}) is to be replaced.
#'
#' @details Indeed a very simple but somewhat limited approach is to impute missing values from observed ones chosen at random with replacement (MCAR), assuming that \deqn{p(R|Z_{obs}, Z_{mis}) = p(R|\phi)}. Sampling with replacement is important since it continues to favor values with higher incidence (preserving the MCAR empirical distribution). It  may also be combined with apply for matrix imputation drills, but keep in mind that it is experimental (actually, I wrote this for teaching purposes).
#'
#' @examples
#' x <- c(1,2,NA,4,5,NA)
#' randomImput(x)
#'
#' if (interactive()) {
#' n = 100
#' mat <- matrix(ncol=3, nrow=n)
#' for(i in 1:n){
#' mu = mean(randomImput(x))
#' med = median(randomImput(x))
#' mod = Mode(randomImput(x))
#' mat[i,] <- c(mu, med, mod[1])
#' }
#' print(mat)
#' }
#' @export
`randomImput` <- function(x){
  y <- is.na(x)
  xx <- x[!y]
  x[y] <- sample(x=xx,size=sum(y),replace=TRUE)
  return(x)
}
NULL
