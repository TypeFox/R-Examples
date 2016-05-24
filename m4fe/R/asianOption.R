"permn" <-function(x, fun = NULL, ...) {
  if(is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == x) x <- seq(x)
  n <- length(x)
  nofun <- is.null(fun)
  out <- vector("list", gamma(n + 1))
  p <- ip <- seqn <- 1:n
  d <- rep(-1, n)
  d[1] <- 0
  m <- n + 1
  p <- c(m, p, m)
  i <- 1
  use <-  - c(1, n + 2)
  while(m != 1) {
    out[[i]] <- if(nofun) x[p[use]] else fun(x[p[use]], ...)
    i <- i + 1
    m <- n
    chk <- (p[ip + d + 1] > seqn)
    m <- max(seqn[!chk])
    if(m < n)
      d[(m + 1):n] <-  - d[(m + 1):n]
    index1 <- ip[m] + 1
    index2 <- p[index1] <- p[index1 + d[m]]
    p[index1 + d[m]] <- m
    tmp <- ip[index2]
    ip[index2] <- ip[m]
    ip[m] <- tmp
  }
  out
}

#' @title Asian Option Price
#' @description Returns the price of an asian option using a binomial tree approach
#' @param S the initial stock price
#' @param K the strike price
#' @param r the risk free (continuously compounded interest rate)
#' @param delta the annual dividend rate
#' @param sigma the volatility
#' @param t the expiration time (default one year)
#' @param call TRUE if option is a call, FALSE is option is a put
#' @param arithmetic TRUE if arithmetic average is used, FALSE if geometric average is used
#' @param price TRUE if average price is used, FALSE if average strike is used
#' @param h the number of subdivisions between 0 and t (default 10)
#' @details Uses a forward tree to compute u and d. p is the risk-neutral probability
#' @examples asianOption(40, 39, 0.05, 0, 0.3, 3/12, call=FALSE, arithmetic=TRUE, price=TRUE, h=3)
#' @export
asianOption <- function(S, K, r, delta, sigma, t = 1, call=TRUE, arithmetic=TRUE, price=TRUE, h=10) {
  u = exp( (r-delta)*h + sigma*sqrt(h))
  d = exp( (r-delta)*h - sigma*sqrt(h))
  p = (exp((r-delta)*(t/h)) - d)/(u-d) 
  
  paths <-numeric(0) # compute the possible (ordered) paths
  for(i in 0:h) { 
    path = c(rep(0,i),rep(1,(h-i)))
    paths = c(paths, unique(permn(path)))
  }
  
  payoffs <- numeric(0) # determine the payoff for each path
  probabilities <- numeric(0) # determine the probability for each path
  averages <- numeric(0)
  
  for(path in paths) {
    previousPrice = S
    probability = 1
    avg <- numeric(0)
    payoff <- numeric(0)
    
    prices <- numeric(0)
    for(num in path) {
      if(num==1) { # u 
        prices = c(prices, previousPrice*u)
        previousPrice = previousPrice * u
        probability = probability * p
      } else { # d
        prices = c(prices, previousPrice*d)
        previousPrice = previousPrice * d
        probability = probability * (1-p)
      }
    }
    
    if(arithmetic==TRUE) {
      avg = mean(prices)
    } else { # geometric average
      avg = exp(mean(log(prices)))
    }
    
    if(price==TRUE) { 
      if(call==TRUE) {
        payoff=max(0, avg-K) 
      } else { # call == FALSE
        payoff = max(0, K-avg)
      }
    } else { # strike option (K=avg, S = prices[h])
      if(call==TRUE) {
        payoff=max(0, prices[h]-avg)
      } else {
        payoff=max(0, avg-prices[h])
      }
    }
    payoffs = c(payoffs, payoff)
    probabilities = c(probabilities, probability)
    averages = c(averages, avg)
  }
  payoffs <<- payoffs
  averages <<- averages
  probabilities <<- probabilities
  
  exp(-r*t)*sum(probabilities*payoffs)
}