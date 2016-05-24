################################################################################
#' Function for simulating a time series 
#'
#' Function for simulating a time series and creating a sts-object 
#' As the counts are generated using a negative binomial distribution
#' one also gets the (1-alpha) quantile for each timepoint (can be interpreted
#' as an in-control upperbound for in-control values).
#' The baseline and outbreaks are created as in Noufaily 2012.
################################################################################
# Parameters for the negbin
#' @references An improved algorithm for outbreak detection in multiple surveillance systems, Noufaily, A., Enki, D.G., Farrington, C.P., Garthwaite, P., Andrews, N.J., Charlett, A. (2012), Statistics in Medicine, published online.

###
#' @param theta baseline frequency of reports
#' @param beta time trend
#' @param gamma1 seasonality
#' @param gamma2 seasonality
#' @param m seasonality
#' @param overdispersion overdispersion (size in rnbinom for the parameterization with mean and size)
#' @param delayMax maximal delay in time units
###
# Parameters for the time series
###
#' @param dates dates of the time series
#' @param densityDelay density distribution for the delay

###
# Parameters for outbreaks
###
#' @param sizesOutbreak sizes of all the outbreaks (vector) 
#' @param datesOutbreak dates of all the outbreaks (vector)
#' # alpha
#' @param alpha alpha for getting the (1-alpha) quantile of the negative binomial distribution at each timepoint
#' @examples
#' set.seed(12345)
#' # Time series parameters
#' scenario4 <- c(1.6,0,0.4,0.5,2)
#' theta <- 1.6
#' beta <- 0
#' gamma1 <-0.4 
#' gamma2 <- 0.5
#' overdispersion <- 1
#' m <- 1
#' # Dates 
#' firstDate <- "2006-01-01"
#' lengthT=350
#' dates <- as.Date(firstDate,origin='1970-01-01') + 7 * 0:(lengthT - 1)
#' # Maximal delay in weeks
#' D=10
#' # Dates and sizes of the outbreaks
#' datesOutbreak <- c(as.Date("2008-03-30"),as.Date("2011-09-25",origin="1970-01-01"))
#' sizesOutbreak <- c(2,5)
#' # Delay distribution
#' data("salmAllOnset")
#' in2011 <- which(formatDate(epoch(salmAllOnset), "%G") == 2011)
#' rT2011 <- salmAllOnset@@control$reportingTriangle$n[in2011,]
#' densityDelay <- apply(rT2011,2,sum, na.rm=TRUE)/sum(rT2011, na.rm=TRUE)
#' # alpha for the upperbound
#' alpha <- 0.05
#' # Create the sts with the full time series
#' stsSim <- sts_creation(theta=theta,beta=beta,gamma1=gamma1,gamma2=gamma2,m=m,overdispersion=overdispersion,dates=dates,sizesOutbreak=sizesOutbreak,datesOutbreak=datesOutbreak,delayMax=D,
#'                          densityDelay=densityDelay,alpha=alpha)
#' plot(stsSim)                         
#' @export
sts_creation <- function(theta,beta,gamma1,gamma2,m,overdispersion,dates,
                         sizesOutbreak,datesOutbreak,delayMax,alpha,
                         densityDelay){
  lengthT <- length(dates)
  firstDate=dates[1]
  # Baseline
  observed <- rep(NA,lengthT)
  upperbound <- rep(NA,lengthT)
  state <- logical(length=lengthT)

  for (t in 1:lengthT) {
    if (m==0){season=0}
    if (m==1){season=gamma1*cos(2*pi*t/52)+ gamma2*sin(2*pi*t/52)}
    if (m==2){season=gamma1*cos(2*pi*t/52)+ gamma2*sin(2*pi*t/52)+gamma1*cos(4*pi*t/52)+ gamma2*sin(4*pi*t/52)}
    mu <- exp(theta + beta*t + season)
    observed[t] <- rnbinom(mu=mu,size=overdispersion,n=1)
    upperbound[t] <- qnbinom(mu=mu,size=overdispersion,p=(1-alpha))
    
  }
  
  # Outbreaks
  nOutbreaks <- length(sizesOutbreak)
  if (nOutbreaks>1){
    dens <- lognormDiscrete(Dmax=20,logmu=0,sigma=0.5)
    for (i in 1:nOutbreaks){
      tOutbreak <- which(dates==datesOutbreak[i])
      numberOfCases <- rpois(n=1,lambda=sizesOutbreak[i]*(mu*(1+mu/overdispersion)))
      cases <- rep(0,length(dens))
      if (numberOfCases!=0){
        for (case in 1:numberOfCases){
          t <- sample(x=1:length(dens),size=1,prob=dens)
          cases[t] <- cases[t] + 1
        }
      }
      cases <- cases[cases>0]
      if(sum(cases)>0){
      observed[tOutbreak:(tOutbreak+length(cases)-1)] <- observed[tOutbreak:(tOutbreak+length(cases)-1)] + cases
      state[tOutbreak:(tOutbreak+length(cases)-1)] <- TRUE
      }
    }
    
  }
  observed <- observed[1:lengthT]
  
  # Reporting triangle
  if (!is.null(densityDelay)){
    # use density delay  
    n <- matrix(0, lengthT, delayMax + 1,dimnames=list(as.character(dates),NULL))
    for (t in 1:lengthT){
      if(observed[t]!=0){
        for (case in 1:observed[t]){
          delay <- sample(x=0:delayMax,size=1,prob=densityDelay)
          if (delay > delayMax) {delay <- delayMax}
          n[t, delay + 1] <- n[t, delay + 1] + 1
        }
      }
    }    
  }
  else{
    # Using a poisson as for the outbreaks because it looks good
    n <- matrix(0, lengthT, D + 1,dimnames=list(as.character(dates),NULL))
    for (t in 1:lengthT){
      if(observed[t]!=0){
        for (case in 1:observed[t]){
          delay <- rpois(n=1, lambda=1.5) 
          if (delay > D) {delay <- D}
          n[t, delay + 1] <- n[t, delay + 1] + 1
        }
      }
    }
  }
  # Create the sts
  firstYear <- isoWeekYear(as.Date(firstDate,origin="1970-01-01"))$ISOYear
  firstWeek <- isoWeekYear(as.Date(firstDate,origin="1970-01-01"))$ISOWeek
  newSts <- new("sts", epoch = as.numeric(dates), start = c(2006, 1), upperbound = as.matrix(upperbound),
                freq = 52, observed = observed, state = as.matrix(state), epochAsDate = TRUE)
  newSts@control$reportingTriangle$n <- n
  return(newSts)
}
################################################################################
# FUNCTION FOR DISCRETIZING THE LOG NORM DISTRIBUTION
################################################################################
lognormDiscrete <- function(Dmax=20,logmu=0,sigma=0.5){
  Fd <- plnorm(0:Dmax, meanlog = logmu, sdlog = sigma)
  FdDmax <- plnorm(Dmax, meanlog = logmu, sdlog = sigma)
  
  #Normalize 
  prob <- diff(Fd)/FdDmax
  return(prob)
}

