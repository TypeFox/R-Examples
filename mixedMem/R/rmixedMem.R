#' Simulate Mixed Membership Data
#' 
#' Simulate data from a mixed membership model.
#' 
#' \code{rmixedMem} simulates data from a mixed membership model given the specified parameters and dimensions. The function
#' returns a random sample of observations \code{obs}, context indicators \code{Z}, and group membership scores \code{lambda}.
#'    
#' @param Total the number of individuals in the sample.
#' @param J the number of variables observed on each individual.
#' @param Rj a vector of positive integers of length \code{J} specifying the number of repeated measurements
#'  for each variable.
#' @param Nijr an array of dimension (\code{Total}, \code{J}, \code{max(Rj)}) indicating the number
#'  of ranking levels for each replication. For multinomial and Bernoulli
#'  variables, \code{Nijr}[i,j,r] = 1. For rank variables, \code{Nijr}[i,j,r] indicates the
#'  number of items ranked for each individual.
#' @param K the number of latent sub-populations.
#' @param Vj a vector of length \code{J} specifying the number of possible candidates
#'  for each variable. For a Bernoulli variable \code{Vj}[j] = 1. For a multinomial
#'   or rank variable, \code{Vj}[j] is the number of possible categories/items.
#' @param dist a vector of strings of length \code{J} specifying variable types. Options are
#'  "bernoulli", "multinomial" or "rank" corresponing to the distributions
#'   of the observed variables.
#' @param alpha a positive \code{K}-length vector  which is the parameter for the Dirichlet
#'  distribution of membership scores.
#' @param theta an 3 way array of dimensions (\code{J},\code{K},\code{max(Vj)}) which governs the variable
#'  distributions. Parameter \code{theta}[j,k,] governs the distirbution of responses on variable j for an inidvidually completely in sub-population k.
#'  If the number of items/categories differs across variables, any
#'  unusued portions of \code{theta} should be set to 0.
#'  @param lambda an optional matrix of dimensions (\code{Total}, \code{K}) containing the membership scores for each individual. If the \code{lambda}
#'  argument is not specified, the group membership scores will be automatically sampled from a Dirichlet(\code{alpha})
#' @return \code{rmixedMem} returns a list containing a three items: A matrix of group memberships scores \code{lambda}, 
#' an array of context indicators \code{Z} and an array of observations \code{obs}.
rmixedMem <- function(Total, J, Rj, Nijr, K, Vj, dist, theta, alpha, lambda = NULL)
{
  if(is.null(lambda)){
    lambda <- gtools::rdirichlet(Total, alpha)
  }
  Z <-array(-1, dim = c(Total, J, max(Rj), max(Nijr)))
  obs <- array(-1, dim = c(Total, J, max(Rj), max(Nijr)))
  for(i in 1:Total)
  {
    for(j in 1:J)
    {
      for(r in 1:Rj[j])
      {
        if(dist[j] =="multinomial") {
          sub.pop <- sample.int(K, size = 1, prob = lambda[i,]) # sub-population which governs the response
          obs[i,j,r,1] <- sample.int(Vj[j], size = 1, prob = theta[j,sub.pop,c(1:Vj[j])])-1 # must be in 0:(Vj[j]-1)
          Z[i,j,r,1] <- sub.pop-1
        }
        if(dist[j] == "rank") {
          select = c()
          for(n in 1:Nijr[i,j,r])
          {
            sub.pop <- sample.int(K, size = 1, prob = lambda[i,]) # sub-population which governs the response
            prob = theta[j,sub.pop,c(1:Vj[j])]
            prob[c(select)] <- 0
            select <- c(select, sample.int(Vj[j], size = 1, replace = F, prob = prob))
            obs[i,j,r,n] <- select[n]-1  # must be in 0:(Vj[j]-1)
            Z[i,j,r,n] <- sub.pop-1
          }
        }
        if(dist[j] == "bernoulli"){
          sub.pop <- sample.int(K, size = 1, prob = lambda[i,]) # sub-population which governs the response
          obs[i,j,r,1] <- rbinom(n = 1, size = 1, prob = theta[j,sub.pop,1])
          Z[i,j,r,1] <- sub.pop-1
        }
      }
    }
  }
  
  out <- list(lambda = lambda, Z = Z, obs = obs)
  return(out)
}