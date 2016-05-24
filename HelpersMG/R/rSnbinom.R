#' rSnbinom returns random numbers for the sum of random variable with negative binomial distributions
#' @title Random generation for the sum of random variable with negative binomial distributions. 
#' @author Marc Girondot
#' @return rSnbinom returns random number
#' @param n number of observations.
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success in each trial. 0 < prob <= 1.
#' @param mu alternative parametrization via mean.
#' @description Random numbers for the sum of random variable with negative binomial distributions.
#' @family Distribution of sum of random variable with negative binomial distributions
#' @examples
#' alpha <- c(2.1, 2.05, 2)
#' mu <- c(10, 30, 20)
#' rep <- 100000
#' distEmpirique <- rSnbinom(n=rep, size=alpha, mu=mu)
#' tabledistEmpirique <- rep(0, 301)
#' names(tabledistEmpirique) <- as.character(0:300)
#' tabledistEmpirique[names(table(distEmpirique))] <- table(distEmpirique)/rep
#' 
#' plot(0:300, dSnbinom(0:300, size=alpha, mu=mu, infinite=1000), type="h", bty="n", 
#'    xlab="x", ylab="Density", ylim=c(0,0.02))
#' plot_add(0:300, tabledistEmpirique, type="l", col="red")
#' legend(x=200, y=0.02, legend=c("Empirical", "Theoretical"), 
#'    text.col=c("red", "black"), bty="n")
#' @export

rSnbinom <- function(n=1, 
                     size=NULL, 
                     prob=NULL, mu=NULL) {
  
  # prob=NULL; mu=NULL; log = FALSE; infinite=10
  
  if (is.null(mu) + is.null(size) + is.null(prob) != 1) stop("Two values among mu, size and prob must be provided")
  
  m <- max(c(length(size), length(prob), length(mu)))
  if (!is.null(mu)) mu <- rep(mu, m)[1:m]
  if (!is.null(size)) size <- rep(size, m)[1:m]
  if (!is.null(prob)) prob <- rep(prob, m)[1:m]
  
  if (is.null(prob)) prob <- size/(size+mu)
  if (is.null(mu)) mu <- size/prob - size
  if (is.null(size))  size  <- (prob * mu) / (1 - prob)  
#  if (length(prob)<length(size)) prob <- rep(prob, length(size))[1:length(size)]
#  if (length(size)<length(prob)) size <- rep(size, length(prob))[1:length(prob)]
  
  m <- matrix(1:length(size), nrow=1)
  
  rl <- apply(m, MARGIN=2, function(i) rnbinom(n, size=size[i], prob=prob[i]))
  
  if (class(rl)=="integer") rl <- as.data.frame(matrix(rl, nrow=1))
  
  return(apply(rl, MARGIN=1, sum))
}