submean <- function(y, PSU, N, M, Nl, m.weight, n.weight, method = "simple",
                    level = 0.95)
{ 
### input treatment

  method <- match.arg(method, c("simple", "ratio"))

  # population size
  if(missing(N)){
    stop("Invalid input: population size ", sQuote("N"), " must be given.")
  }
  # data and primary sampling units (PSU)
  if (!(is.numeric(y) && is.vector(y))){
    stop("Invalid input: ", sQuote("y"), " has to be a numeric vector.")
  }
  if( !is.factor(PSU) ){
    stop("Invalid input: ", sQuote("PSU"), " has to be a factor.")
  }else{
    m <- nlevels(PSU) # number of sampled primary units
    nl <- table(PSU) # number of samples in each PSU
  }
  if( length(y) != length(PSU) ){
    stop("Invalid input: ", sQuote("y"), " and ", sQuote("PSU")," must have the same length.")
  }
  if( N < length(y)){
    stop("Invalid input: population size ", sQuote("N"), " must be larger than sample size.")
  }

  # define data set
  data <- data.frame(target=as.numeric(y),s1=as.factor(PSU))
  # omit NA values in data
  if(any(is.na(data))){
    data <- na.omit(data)
    warning("Missing values will be ignored in calculations.")
  }
  # each strata need more than one observation
  if(any(table(data$s1) < 2)){
    stop("Invalid input: PSU with only one observation was entered.")
  }

  ## stage group weights
  if(missing(M) == missing(m.weight) || missing(Nl) == missing(n.weight)){  
    stop("Invalid input: wether ", sQuote("M"), " and ", sQuote("Nl"), " or ", sQuote("m.weight"), " and ", sQuote("n.weight"), " must be given.")
  }
  if(!missing(M) && !missing(Nl)){
    if(M < m){
      stop("Invalid input: the number of primary units in the population", sQuote("M")," has to be greater than the number of primary units in the sample ", m , ".")
    }
    if(is.vector(Nl) && m != length(Nl) ){
      stop("Invalid input: the number of PSU given in ", sQuote("Nl")," is not correct according the sample data.")
    }
    # with finite poulation correction
    fpc <- TRUE
    fpc1 <- (M - m) / M
    fpc2 <- (Nl - nl) / Nl
  }else{ 
    M = m / m.weight
    if(is.vector(n.weight) && m != length(n.weight) ){
      stop("Invalid input: the length of ", sQuote("n.weight")," does not equal the number of given PSU.")
    }else{
      Nl = nl / n.weight
    }
    # without finite population correction
    fpc <- FALSE
    fpc1 <- 1
    fpc2 <- 1
  }

  # level
  if( level < 0 || level > 1 ){
    stop("Invalid input: ", sQuote("level")," has to be probability between 0 and 1.")
  }else{
    q <- qnorm((1 + level) / 2) 
  }

### compute mean in each group
  ylbar <- ylsum <- varl <- NA
  # divide data according to primary sample units s1
  splitted <- split(data, data$s1)
  for (i in 1:m){
    ylbar[i] <- mean(splitted[[i]]$target)
    varl[i] <- var(splitted[[i]]$target)
  }

### compute overall mean
  if( method == "simple"){
    Y <- M / m / N * sum(Nl * ylbar)
    VB <- sum( ( Nl * ylbar - 1 / m * sum(Nl * ylbar) ) ^ 2 )
  }
  if( method == "ratio"){
    Y <- sum(Nl * ylbar) / sum(Nl)
    VB <- sum( ( Nl*ylbar - Nl*Y ) ^ 2 )
  }
### compute overall variance and confidence interval
  Var <- (M/N)^2 * (fpc1 / (m*(m-1)) * VB + (1/(m*M)) * sum(Nl^2 * fpc2 * varl / nl))
  CIo <- Y + q * sqrt(Var)
  CIu <- Y - q * sqrt(Var)

### return argument
  ret <- list()
  ret$call <- list(y = y, PSU = PSU, N = N, M = M, fpc = fpc,  method = method, level = level)
  ret$mean <- Y
  ret$se <- sqrt(Var)
  ret$ci <- c(CIu, CIo)
  structure(ret, class = "submean")
}
