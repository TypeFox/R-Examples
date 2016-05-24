require(rootSolve)
#
mrchk <- function(XYdt,P1 = 10,P2 = 1,P3 = 0,...) {
  # define header
  names(XYdt) <- c("XC","YC")
  # define nonlinear system and solve it
  FFn <- nls(
    YC ~ P1 * exp(P2 * (XC ^ (0.5)) - P3 * (XC ^ 3)),
    start = list(P1 = P1,P2 = P2,P3 = P3),
    data = XYdt,na.exclude
  )
  # return output from anaysis
  FFn
}

mrgsn <- function(XYdt,...) {
  # define header
  names(XYdt) <- c("XC","YC")
  ##  controll<-nls.control(maxiter=50,
  ##  tol=1e-10, minFactor = 1/1024,
  ##  printEval = FALSE, warnOnly = FALSE)
  # define nonlinear system and solve it
  FFn <- nls(
    YC ~ P1 + P2 * (XC) ^ 0.5 + P3 * XC,
    start = list(P1 = 50,P2 = 1,P3 = 0),
    data = XYdt,na.exclude
  )
  # return output from anaysis
  FFn
}

tello <- function(XYdt,...) {
  # tello's method is highly dependent of guess values to obtain its parameters
  # so the method for calculating a more approximated guess value is described in
  # his article and implemented below
  #
  # calculate derivative of data
  df.sys <- diff(XYdt[,1]) / diff(XYdt[,2])
  nrow.sys <- nrow(XYdt)
  S <- smooth.spline(df.sys)$fit$coef[1:nrow.sys]
  XC <- XYdt[,1]
  # merge data in a dataframe
  coef.est <- LLSRxy(XC,S)
  # define header
  names(coef.est) <- c("XC","S")
  # define nonlinear system to estimate parameters and solve it
  FFnEst <- nls(S ~ (P2 / P1) + XC / P1,
                start = list(P1 = -.1,P2 = .001),
                data = coef.est,na.exclude)
  # obtain 1st and 2nd coefficients
  coefEst <- summary(FFnEst)
  coef.1 <- coefEst$coefficients[1]
  coef.2 <- coefEst$coefficients[2]
  # calculate third coefficient
  x <- mean(XYdt[,1])
  y <- mean(XYdt[,2])
  coef.3 <- log(exp(y) * ((x + coef.2) ^ (-coef.1)))
  # any negative difference between coefficients and mass fraction will result in error,
  # so the guess values are set to be bigger than the minimum value in the dataset
  if (coef.2 < -min(XYdt[,1]))
    coef.2 <- -(min(XYdt[,1]) - 0.001)
  #
  # define header
  names(XYdt) <- c("XC","YC")
  # define nonlinear system and solve it using estimated parameters as guess
  FFn <- nls(
    XC ~ exp((YC - P3) / P1) - P2,
    #YC ~ P1*log(XC + P2) + P3,
    start = list(P1 = coef.1, P2 = coef.2, P3 = coef.3),
    #algorithm = "port",
    #lower = c(P1 = -Inf, P2 = -(min(XYdt[,1])-0.001), P3 = -Inf),
    data = XYdt, na.action = na.exclude
  )
  # return output from anaysis
  FFn
}
