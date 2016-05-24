#'Sampler for the concentration parameter of a Dirichlet process
#'
#'Sampler udating the concentration parameter of a Dirichlet process given
#'the number of observations and a Gamma(\code{a}, \code{b}) prior, following the augmentation
#'strategy of West, and of Escobar and West.
#'
#'@param alpha_old the current value of alpha
#'
#'@param n the number of data points
#'
#'@param K current number of cluster
#'
#'@param a shape hyperparameter of the Gamma prior
#'on the concentration parameter of the Dirichlet Process.
#'Default is \code{0.0001}.
#'
#'@param b scale hyperparameter of the Gamma prior
#'on the concentration parameter of the Dirichlet Process.
#'Default is \code{0.0001}. If \code{0} then the concentration is fixed and this function
#'returns \code{a}.
#'
#'@details A Gamma prior is used.
#'
#'@references M West, Hyperparameter estimation in Dirichlet process mixture models,
#'Technical Report, Duke University, 1992.
#'
#'@references MD Escobar, M West, Bayesian Density Estimation and Inference Using Mixtures
#'\emph{Journal of the American Statistical Association}, 90(430):577-588, 1995.
#'
#'@importFrom stats rbeta rgamma runif
#'
#'@export
#'
#'@examples
#' #Test with a fixed K
#' ####################
#'
#' alpha_init <- 1000
#' N <- 10000
#' #n=500
#' n=10000
#' K <- 80
#' a <- 0.0001
#' b <- a
#' alphas <- numeric(N)
#' alphas[1] <- alpha_init
#' for (i in 2:N){
#'  alphas[i] <- sample_alpha(alpha_old = alphas[i-1], n=n, K=K, a=a, b=b)
#' }
#'
#' postalphas <- alphas[floor(N/2):N]
#' alphaMMSE <- mean(postalphas)
#' alphaMAP <- density(postalphas)$x[which.max(density(postalphas)$y)]
#'
#' expK <- sum(alphaMMSE/(alphaMMSE+0:(n-1)))
#' round(expK)
#'
#'
#'  prioralpha <- data.frame("alpha"=rgamma(n=5000, a,1/b),
#'                          "distribution" =factor(rep("prior",5000),
#'                          levels=c("prior", "posterior")))
#'
#'  library(ggplot2)
#'  p <- (ggplot(prioralpha, aes(x=alpha))
#'        + geom_histogram(aes(y=..density..),
#'                         colour="black", fill="white")
#'        + geom_density(alpha=.2, fill="red")
#'        + ggtitle(paste("Prior distribution on alpha: Gamma(", a,
#'                  ",", b, ")\n", sep=""))
#'       )
#'  p
#'
#' postalpha.df <- data.frame("alpha"=postalphas,
#'                          "distribution" = factor(rep("posterior",length(postalphas)),
#'                          levels=c("prior", "posterior")))
#'  p <- (ggplot(postalpha.df, aes(x=alpha))
#'        + geom_histogram(aes(y=..density..), binwidth=.1,
#'                         colour="black", fill="white")
#'        + geom_density(alpha=.2, fill="blue")
#'        + ggtitle("Posterior distribution of alpha\n")
#'        # Ignore NA values for mean
#'        # Overlay with transparent density plot
#'        + geom_vline(aes(xintercept=mean(alpha, na.rm=TRUE)),
#'                     color="red", linetype="dashed", size=1)
#'      )
#'  p
#'
#'
#'
#'
#'
#'

sample_alpha <- function(alpha_old, n, K, a=0.0001, b=0.0001){

  if(b > 0){
    # Sample scale factor in Dirichlet Process
    x <- stats::rbeta(n=1, shape1=alpha_old + 1, shape2=n)
    temp <- (a+K-1) / (n*(b-log(x)))
    pi <- temp/(1+temp)
    u <- stats::runif(1)
    if (u<pi){
        alpha_new  <-  stats::rgamma(n=1, shape=a + K, scale=1/(b - log(x)))
    } else{
      alpha_new <-  stats::rgamma(n=1, shape=a + K - 1, scale=1/(b - log(x)))
    }
    return(alpha_new)
  }else if(b==0){
    alpha_new=a
  }else{
    print("b cannot be negative")
  }
}

