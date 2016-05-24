ss.aipe.reg.coef<-function (Rho2.Y_X = NULL, Rho2.j_X.without.j = NULL, p = NULL, 
                            b.j = NULL, width, which.width = "Full", sigma.Y = 1, sigma.X = 1, 
                            RHO.XX = NULL, Rho.YX = NULL, which.predictor = NULL, Noncentral = FALSE, 
                            alpha.lower = NULL, alpha.upper = NULL, conf.level = 0.95, 
                            degree.of.certainty = NULL, assurance = NULL, certainty = NULL, 
                            Suppress.Statement = FALSE)
{
  if(!requireNamespace("gsl", quietly = TRUE)) stop("The package 'gsl' is needed; please install the package and try again.")
  
  RHO.XX <- as.matrix(RHO.XX)  
  Rho.YX <- as.matrix(Rho.YX)  
  if (!is.null(p)) {			
    if (p == 1) {
      if (!is.null(Rho2.j_X.without.j)) {
        if (Rho2.j_X.without.j != 0) 
          stop("If p=1, how could 'Rho2.j_X.without.j' be nonzero?")
      }
      if (!is.null(RHO.XX)) 
      {                if (dim(RHO.XX) != c(1, 1)) stop("If p=1, 'RHO.XX' should be a 1 x 1 matrix?")
        if (dim(Rho.YX) != c(1, 1)) stop("If p=1, how can 'RHO.XX' not be a 1 x 1 matrix?")
        Rho.YX
      }
    }
  }
  # Thanks to Jan Herman for modificatoin of the function just above here to allow p=1 to work properly. 
  
  
  
  if (!is.null(certainty) & is.null(degree.of.certainty) & 
      is.null(assurance)) 
    degree.of.certainty <- certainty
  if (is.null(assurance) && !is.null(degree.of.certainty) & 
      is.null(certainty)) 
    assurance <- degree.of.certainty
  if (!is.null(assurance) && is.null(degree.of.certainty) & 
      is.null(certainty)) 
    degree.of.certainty <- assurance
  if (!is.null(assurance) && !is.null(degree.of.certainty) && 
      assurance != degree.of.certainty) 
    stop("The arguments 'assurance' and 'degree.of.certainty' must have the same value.")
  if (!is.null(assurance) && !is.null(certainty) && assurance != 
      certainty) 
    stop("The arguments 'assurance' and 'certainty' must have the same value.")
  if (!is.null(degree.of.certainty) && !is.null(certainty) && 
      degree.of.certainty != certainty) 
    stop("The arguments 'degree.of.certainty' and 'certainty' must have the same value.")
  Expected.R2 <- function(Population.R2, N, p) {
    Value <- 1 - ((N - p - 1)/(N - 1)) * (1 - Population.R2) * 
      gsl::hyperg_2F1(1, 1, 0.5 * (N + 1), Population.R2)
    Value <- max(0, Value)
    return(Value)
  }
  if (Noncentral == TRUE & is.null(sigma.X)) 
    sigma.X <- 1
  if (Noncentral == TRUE & is.null(sigma.Y)) 
    sigma.Y <- 1
  if (Noncentral == TRUE & sigma.Y != 1) 
    stop("Since you've specified 'Noncentral=TRUE', all variances should be one (your 'sigma.Y' (i.e., the Y standard deviation) is not one), for a standardized solution.")
  if (Noncentral == TRUE & sigma.X != 1) 
    stop("Since you've specified 'Noncentral=TRUE', all variances should be one (your 'sigma.X' (i.e., the X standard deviation) is not one), for a standardized solution.")
  if (is.null(p) & is.null(RHO.XX)) 
    stop("Since RHO.XX is not specified, you must specify 'p'.")
  if (!is.null(RHO.XX)) {
    if (!(sum(round(RHO.XX, 5) == round(t(RHO.XX), 5)) == 
          dim(as.matrix(RHO.XX))[1] * dim(as.matrix(RHO.XX))[2])) 
      stop("The correlation matrix, 'RHO.XX' should be symmetric.")
  }
  if (is.null(p) & !is.null(RHO.XX)) 
    p <- dim(RHO.XX)[1]
  char.expand(which.width, c("Full", "Lower", "Upper"), nomatch = stop("Problems with 'which.width' specification. You must choose either 'Full', 'Lower', or 'Upper'.", 
                                                                       call. = FALSE))
  if (which.width == "Lower" | which.width == "Upper") 
    stop("At the present time, only the 'which.width' of 'Full' is implemented.", 
         call. = FALSE)
  if (is.null(conf.level)) {
    if (!is.numeric(alpha.lower) | !is.numeric(alpha.upper)) 
      stop("Since 'conf.level' is not specified, you need to correctly specify 'alpha.lower' and 'alpha.upper'.")
    if (alpha.lower < 0 | alpha.lower >= 1 | alpha.upper < 
        0 | alpha.upper >= 1) 
      stop("You have not correctly specified 'alpha.lower' and/or 'alpha.upper'.")
  }
  if (!is.null(conf.level)) {
    if (!is.null(alpha.lower) | !is.null(alpha.upper)) 
      stop("Since 'conf.level' is not specified, you need to correctly specify 'alpha.lower' and 'alpha.upper'.")
    if (!is.null(alpha.lower) | !is.null(alpha.upper)) 
      stop("Since 'conf.level' is specified, do not specifiy 'alpha.lower' and 'alpha.upper'.")
    alpha.lower <- alpha.upper <- (1 - conf.level)/2
  }
  if (!is.null(degree.of.certainty)) {
    if ((degree.of.certainty <= 0) | (degree.of.certainty >= 
                                      1)) 
      stop("The 'degree.of.certainty' must either be NULL or some value greater than .50 and less than 1.", 
           call. = FALSE)
    if (degree.of.certainty <= 0.5) 
      stop("The 'degree.of.certainty' should be > .5 (but less than 1).", 
           call. = FALSE)
  }
  if ((!is.null(Rho2.j_X.without.j) & !is.null(Rho2.Y_X)) & 
      (!is.null(RHO.XX) & !is.null(Rho.YX))) 
    stop("Since 'Rho2.j_X.without.j' and 'Rho2.Y_X' are specified, do not specify 'RHO.XX' or 'Rho.YX' (or vice versa).", 
         call. = FALSE)
  if (!is.null(RHO.XX) & !is.null(Rho.YX)) {
    if (!is.null(Rho2.j_X.without.j) | !is.null(Rho2.Y_X)) 
      stop("Since 'RHO.XX' and 'Rho.YX' are specified, do not specify 'Rho2.Y_X' or 'Rho2.j_X.without.j'.", 
           call. = FALSE)
    Rho2.Y_X <- (Rho.YX %*% solve(RHO.XX) %*% Rho.YX)
    Rho2.j_X.without.j <- 1 - 1/solve(RHO.XX)[which.predictor, 
                                              which.predictor]
    # if-statement added by Jan Herman to cope with the p=1 case
    if (p>1) {
      Rho2.Y_X.without.j <- (Rho.YX[-which.predictor] %*% solve(RHO.XX[-which.predictor, 
                                                                       -which.predictor]) %*% Rho.YX[-which.predictor])
      b.j.tmp <- (solve(RHO.XX) %*% Rho.YX)[which.predictor]
    } else {
      Rho2.Y_X.without.j<-as.matrix(0)   #If p=1 then Rho2.Y_X.without.j will be zero	
    }
    if (!is.null(b.j)) {
      if (round(b.j, 5) != b.j.tmp) 
        stop("The covariance structure implied regression coefficient and 'b.j' are not equal; this is a problem.")
    }
  }
  if (!is.null(Rho2.j_X.without.j) & !is.null(Rho2.Y_X) & (is.null(RHO.XX) | 
                                                           is.null(Rho.YX))) {
    if (is.null(b.j)) 
      stop("Since 'RHO.XX' and 'Rho.YX' are not specified, implying 'Rho2.j_X.without.j' and 'Rho2.Y_X' are specified, 'b.j' must also be specified.", 
           call. = FALSE)
  }
  n0 <- (qnorm(1 - (alpha.lower + alpha.upper)/2)/(width * 
                                                     0.5))^2 * (1 - Rho2.Y_X)/(1 - Rho2.j_X.without.j) * (sigma.Y^2/sigma.X^2) + 
    p + 1
  n1 <- max(ceiling(n0) - 10, 2 * p)
  Diff <- 1
  while (Diff > 0) {
    n1 <- n1 + 1
    E.Rho2.Y_X <- Expected.R2(Population.R2 = Rho2.Y_X, N = n1, 
                              p = p)
    E.Rho2.j_X.without.j <- Expected.R2(Population.R2 = Rho2.j_X.without.j, 
                                        N = n1, p = p)
    CV.i <- (qt(1 - (alpha.lower + alpha.upper)/2, n1 - p - 
                  1))
    SD.i <- sqrt(((1 - E.Rho2.Y_X)/((1 - E.Rho2.j_X.without.j) * 
                                      (n1 - p - 1)))) * (sigma.Y/sigma.X)
    Current.Width <- CV.i * SD.i * 2
    Diff <- Current.Width - width
  }
  N <- n1
  if (!is.null(degree.of.certainty)) {
    if ((degree.of.certainty <= 0) | (degree.of.certainty >= 
                                      1)) 
      stop("The 'degree.of.certainty' must either be NULL or some value greater than zero and less than unity.", 
           call. = FALSE)
    if (degree.of.certainty <= 0.5) 
      stop("The 'degree.of.certainty' should be > .5 (but less than 1).", 
           call. = FALSE)
    E.Rho2.Y_X <- Expected.R2(Population.R2 = Rho2.Y_X, N = N, 
                              p = p)
    E.Rho2.j_X.without.j <- Expected.R2(Population.R2 = Rho2.j_X.without.j, 
                                        N = N, p = p)
    N_M <- (qt(1 - (alpha.lower + alpha.upper)/2, N - p - 
                 1)/(width * 0.5))^2 * ((1 - E.Rho2.Y_X)/(1 - E.Rho2.j_X.without.j)) * 
      (sigma.Y^2/sigma.X^2) * (qchisq(degree.of.certainty, 
                                      N - 1)/(N - p - 1)) + p + 1
    N_M <- ceiling(N_M)
  }
  if (Noncentral == FALSE & is.null(degree.of.certainty)) {
    if (Suppress.Statement == FALSE) 
      print(paste("Necessary sample size such that the expected", 
                  (1 - alpha.lower - alpha.upper) * 100, "confidence interval width is", 
                  width, "is", N))
    return(as.numeric(N))
  }
  if (Noncentral == FALSE & !is.null(degree.of.certainty)) {
    if (Suppress.Statement == FALSE) 
      print(paste("Necessary sample size such that the expected", 
                  (1 - alpha.lower - alpha.upper) * 100, "confidence interval will be no wider than", 
                  width, "with", degree.of.certainty * 100, "certainty is", 
                  N_M))
    return(as.numeric(N_M))
  }
  if (Noncentral == TRUE) {
    if (is.null(b.j)) {
      b.j <- (solve(RHO.XX) %*% Rho.YX)[which.predictor]
      if (is.null(b.j)) 
        stop("b.j must be specified directly or obtained from other combinations of parameters.")
    }
    n2 <- max(N - 6, 2 * p + 1)
    Diff <- 1
    while (Diff > 0) {
      n2 <- n2 + 1
      E.Rho2.Y_X <- Expected.R2(Population.R2 = Rho2.Y_X, 
                                N = n2, p = p)
      E.Rho2.j_X.without.j <- Expected.R2(Population.R2 = Rho2.j_X.without.j, 
                                          N = n2, p = p)
      CI.Result.NC <- ci.reg.coef(b.j = b.j, SE.b.j = NULL, 
                                  s.Y = sigma.Y, s.X = sigma.X, N = n2, p = p, 
                                  R2.Y_X = E.Rho2.Y_X, R2.j_X.without.j = E.Rho2.j_X.without.j, 
                                  conf.level = NULL, R2.Y_X.without.j = NULL, t.value = NULL, 
                                  alpha.lower = alpha.lower, alpha.upper = alpha.upper, 
                                  Noncentral = TRUE, Suppress.Statement = TRUE)
      current.width <- CI.Result.NC$Upper.Limit - CI.Result.NC$Lower.Limit
      Diff <- current.width - width
    }
    N_NC <- n2
    if (Noncentral == TRUE & is.null(degree.of.certainty)) {
      if (Suppress.Statement == FALSE) 
        print(paste("Necessary sample size such that the expected", 
                    (1 - alpha.lower - alpha.upper) * 100, "confidence interval using noncentral methods is", 
                    width, "is", N_NC))
      return(as.numeric(N_NC))
    }
    if (!is.null(degree.of.certainty)) {
      E.Rho2.Y_X <- Expected.R2(Population.R2 = Rho2.Y_X, 
                                N = N_NC, p = p)
      E.Rho2.j_X.without.j <- Expected.R2(Population.R2 = Rho2.j_X.without.j, 
                                          N = N_NC, p = p)
      SE.b.M <- sqrt((((1 - E.Rho2.Y_X) * sigma.Y^2)/(((1 - 
                                                          E.Rho2.j_X.without.j) * (N_NC - p - 1)) * sigma.X^2)))
      CI.R2 <- ci.R2(R2 = E.Rho2.Y_X, df.1 = p, df.2 = N_NC - 
                       p - 1, conf.level = degree.of.certainty, F.value = NULL, 
                     N = NULL, p = NULL, alpha.lower = NULL, alpha.upper = NULL)$Lower
      n3 <- N_NC
      Diff <- 1
      while (Diff > 0) {
        n3 <- n3 + 1
        E.Rho2.Y_X <- Expected.R2(Population.R2 = CI.R2, 
                                  N = n3, p = p)
        E.Rho2.j_X.without.j <- Expected.R2(Population.R2 = Rho2.j_X.without.j, 
                                            N = n3, p = p)
        CI.Result.NC <- ci.reg.coef(b.j = b.j, SE.b.j = NULL, 
                                    s.Y = sigma.Y, s.X = sigma.X, N = n3, p = p, 
                                    R2.Y_X = E.Rho2.Y_X, R2.j_X.without.j = E.Rho2.j_X.without.j, 
                                    conf.level = NULL, R2.Y_X.without.j = NULL, 
                                    t.value = NULL, alpha.lower = alpha.lower, 
                                    alpha.upper = alpha.upper, Noncentral = TRUE, 
                                    Suppress.Statement = TRUE)
        current.width <- CI.Result.NC$Upper.Limit - CI.Result.NC$Lower.Limit
        Diff <- current.width - width
      }
      N.NC_M <- n3
      if (Noncentral == TRUE & !is.null(degree.of.certainty)) {
        if (Suppress.Statement == FALSE) 
          print(paste("Necessary sample size such that the expected", 
                      (1 - alpha.lower - alpha.upper) * 100, "confidence interval using noncentral methods will be no wider than", 
                      width, "with", degree.of.certainty * 100, 
                      "certainty is", N.NC_M, ". Caution, although the method used here seems to work well, it has not been definitively shown to be the optimal method. It will perhaps be best to evaluate the estimated sample size given here with the 'ss.aipe.reg.coef.sensitivity' function."))
        return(as.numeric(N.NC_M))
      }
    }
  }
}
