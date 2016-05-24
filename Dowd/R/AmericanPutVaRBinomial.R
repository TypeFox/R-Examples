#' Estimates VaR of American vanilla put using binomial tree.
#' 
#' Estimates VaR of American Put Option using binomial tree to price the option
#' and historical method to compute the VaR.
#' 
#' @param amountInvested Total amount paid for the Put Option.
#' @param stockPrice Stock price of underlying stock.
#' @param strike Strike price of the option.
#' @param r Risk-free rate.
#' @param volatility Volatility of the underlying stock.
#' @param maturity Time to maturity of the option in days.
#' @param numberSteps The number of time-steps considered for 
#' the binomial model.
#' @param cl Confidence level for which VaR is computed.
#' @param hp Holding period of the option in days.
#' @return VaR of the American Put Option
#' @references Dowd, Kevin. Measuring Market Risk, Wiley, 2007.
#' 
#' Lyuu, Yuh-Dauh. Financial Engineering & Computation: Principles, 
#' Mathematics, Algorithms, Cambridge University Press, 2002.
#' 
#' @author Dinesh Acharya
#' @examples
#' 
#'    # Market Risk of American Put with given parameters.
#'    AmericanPutVaRBinomial(0.20, 27.2, 25, .16, .05, 60, 20, .95, 30)
#' 
#' @export
AmericanPutVaRBinomial <- function(amountInvested, stockPrice, strike, r,
                                   volatility, maturity, numberSteps, cl, hp){
  s <- stockPrice
  x <- strike
  time <- maturity
  sigma <- volatility
  # convert maturity and holding period unit from days to years
  time <- time/360
  hp <- hp/360
  # calculate the length of each interval, dt, and number of steps to end of
  # holding period 
  n <- numberSteps
  dt <- time/n
  m <- round(hp/dt) # number of steps to end of holding period
  # 
  # calculate movements and associated probabilities
  #
  u <- exp(sigma*sqrt(dt))
  d <- 1/u
  a <- exp(r*dt)
  p <- (a-d)/(u-d)
  jspan <- n*.5
  jspan <- -((jspan>=0)*floor(jspan)+(jspan<0)*ceiling(jspan)) # j-th node 
  # offset number
  ispan <- round(time/dt)%%2 # i-th node offset number
  i <- ispan:(n+ispan) # i-th node numbers
  j <- jspan:(n+jspan) # j-th node numbers
  
  # expand i and j to eliminate for loop
  jex <- matrix(rep(j, each=length(i)), ncol=length(i), byrow=TRUE)
  iex <- matrix(rep(i, each=length(j)), nrow=length(j))
  #
  # asset price at nodes, matrix is flipped so tree appears correct visually
  pr <- t(apply(apply((s*(u^jex)*(d^(iex-jex))), 2, rev), 1, rev))
  # get upper triangle of pr
  lower.tri(pr)
  pr[lower.tri(pr)] <- 0
  #
  # option valuation along tree
  opt <- matrix(0,nrow(pr),ncol(pr))
  opt[,n+1] <- pmax(x-pr[,n+1],0) # determine final option value  from 
  # underlying price
  for(l in seq(n,1,by=-1)){
    k=1:l
    # probable option values discounted back one time step 
    discopt=(p*opt[k,l+1]+(1-p)*opt[k+1,l+1])*exp(-r*dt)
    # option value is max of X - current price or discopt
    opt[,l]=c(pmax(x-pr[1:l,l],discopt),rep(0,n+1-l))
  }
  # initial option price and number of options in portfolio
  initialOptionPrice <- opt[1,1]
  numberOptions <- amountInvested/initialOptionPrice
  # option tree values at end of holding period
  endHpOptionPrice <- opt[,m+1]
  endHpOptionPrice <- endHpOptionPrice[1:(m+1)]
  # returns<-(endHpOptionPrice-initialOptionPrice)/initialOptionPrice
  # option position Profit and Loss
  profitOrLoss <- (endHpOptionPrice-initialOptionPrice)*numberOptions
  #compute VaR
  VaR <- HSVaR(profitOrLoss, cl)
  return(VaR)
}