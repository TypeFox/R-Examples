#' Parametric per trait Repeatabilities
#'
#' Estimates the variance in the sample not due to measurment error
#'
#' @note Requires at least two observations per individual
#' @param ID indentity of individuals
#' @param ind.data individual measurments
#' @return vector of repeatabilities
#' @export
#' @references Lessels, C. M., and Boag, P. T. (1987).
#'    Unrepeatable repeatabilities: a common mistake.
#'    The Auk, 2(January), 116-121.
#' @author Guilherme Garcia
#' @examples
#' num.ind = length(iris[,1])
#' ID = rep(1:num.ind, 2)
#' ind.data = rbind(iris[,1:4], iris[,1:4]+array(rnorm(num.ind*4, 0, 0.1), dim(iris[,1:4])))
#' CalcRepeatability(ID, ind.data)
CalcRepeatability <- function (ID, ind.data)
{
  models.list <- apply (ind.data, 2, function (vec){return (lm (vec ~ as.factor(ID)))})
  models.list <- lapply (models.list, anova)
  rep.itself <- function (lm.model){
    msq <- lm.model$'Mean Sq' ## 1 entre, 2 dentro
    s2a <- (msq[1] - msq[2])/2
    output <- s2a / (s2a + msq[2])
    return (output)
  }
  output <- sapply (models.list, rep.itself)
  names (output) <- colnames (ind.data)
  return (output)
}
