#<<BEGIN>>
dbern <- function(x,prob=.5,log=FALSE)
#TITLE The Bernoulli Distribution
#NAME bernoulli
#KEYWORDS distribution
#DESCRIPTION
#Density, distribution function, quantile function and random generation
#for the Bernoulli distribution with probability equals to \samp{prob}.
#INPUTS
#{x,q}<<vector of quantiles.>>
#{p}<<vector of probabilities.>>
#{n}<<number of observations. If \samp{length(n) > 1}, the length is taken to be the number required.>>
#[INPUTS]
#{prob}<<vector of probabilities of success of each trial.>>
#{log, log.p}<<logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.>>
#{lower.tail}<<logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.>>
#DETAILS
#These fonctions use the corresponding functions from the \code{\link{binomial}} distribution with argument \samp{size = 1}.
#Thus, 1 is for success, 0 is for failure.
#VALUE
#\samp{dbern} gives the density, \samp{pbern} gives the distribution function,
#\samp{qbern} gives the quantile function, and \samp{rbern} generates random deviates.
#EXAMPLE
#rbern(n=10,prob=.5)
#rbern(n=3,prob=c(0,.5,1))
#SEE ALSO
#\code{\link{Binomial}}

#CREATED 08-02-20
#--------------------------------------------
{
  return(dbinom(x, 1, prob=prob, log = log))}

#<<BEGIN>>
pbern <- function(q,prob=.5,lower.tail = TRUE, log.p = FALSE)
#ISALIAS dbern
#--------------------------------------------
{
  return(pbinom(q, 1, prob=prob, lower.tail = lower.tail, log.p = log.p))
  }

#<<BEGIN>>
qbern <- function(p,prob=.5,lower.tail=TRUE, log.p = FALSE)
#ISALIAS dbern
#--------------------------------------------
{
  return(qbinom(p, 1, prob=prob, lower.tail = lower.tail, log.p = log.p))
  }


#<<BEGIN>>
rbern <- function(n,prob=.5)
#ISALIAS dbern
#--------------------------------------------
{
  return(rbinom(n, 1, prob=prob))
  }
