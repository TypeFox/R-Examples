emg.nllik <- function(x, mu, sigma, lambda)
{
  #print(paste('nllemg',length(x), mu, sigma, lambda, '=>', demg(x, mu, sigma, lambda)))
  -sum(demg(x, mu, sigma, lambda, log=TRUE))
}
