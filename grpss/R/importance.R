#' Arrange and visualize the importance of groups
#' @description Arranges and visualizes the importance of groups for the results of
#' \code{\link{grp.criValues}}.
#' @param grp.values The result from \code{\link{grp.criValues}}.
#' @param n Number of important groups to display. The default is the top 10 groups.
#' @param plot A logical value indicating whether to make a barplot of the importance of
#' groups.
#' @details This function arranges the values of screening criterion from the most important
#' to the least important, and then make a barplot to visualize the importance.
#' @return A matrix containing the first \code{n} important group indices and values of
#' screening criterion.
#' @author Debin Qiu, Jeongyoun Ahn
#' @seealso \code{\link{grp.criValues}}
#' @examples
#' library(MASS)
#' n <- 30 # sample size
#' p <- 3  # number of predictors in each group
#' J <- 50 # number of groups
#' group <- rep(1:J,each = 3)  # group indices
#' Sigma <- diag(p*J)  # covariance matrix
#' X <- mvrnorm(n,seq(0,5,length.out = p*J),Sigma)
#' beta <- runif(12,-2,5)  # coefficients
#' y <- X%*%matrix(c(beta,rep(0,p*J-12)),ncol = 1) + rnorm(n)
#'
#' crivalues <- grp.criValues(X,y,group)  # gSIS
#' importance(crivalues, n = 20)
#' @export
#'
#' @importFrom graphics barplot
#' @importFrom grDevices topo.colors


importance <- function(grp.values, n = 10, plot = TRUE) {
  grp.values <- grp.values[order(grp.values[,2], decreasing = TRUE),]
  if (plot)
    barplot(grp.values[1:n,2], col = topo.colors(n),main = "Importance plot",
            las = 2, names.arg = grp.values[1:n,1],xlab = "group indices",
            ylab = "importance")
  return(grp.values[1:n,])
}
