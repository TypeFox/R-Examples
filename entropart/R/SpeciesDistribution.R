as.SpeciesDistribution <-
function (x) 
{
  UseMethod("as.SpeciesDistribution")
}


as.SpeciesDistribution.data.frame <-
function (x) 
{
  if (any(x < 0)) stop("Species distribution abundances or probabilities must be positive.")

  # Try to save the names before applying as.numeric
  spNames <- names(x)
  spD <- as.numeric(as.matrix(x))
  if (length(spNames) == length(spD))
    names(spD) <- spNames
  
  class(spD) <- c("SpeciesDistribution", class(spD))
  return(spD)
}


as.SpeciesDistribution.numeric <-
function (x) 
{
  if (any(x < 0)) stop("Species distribution abundances or probabilities must be positive.")
  spD <- x
  class(spD) <- c("SpeciesDistribution", class(spD))
  return(spD)
}


as.SpeciesDistribution.integer <-
function (x) 
{
  return(as.SpeciesDistribution.numeric(x))
}


is.SpeciesDistribution <-
function (x) 
{
  inherits(x, "SpeciesDistribution")
}


as.ProbaVector <-
function (x, Correction = "None", Unveiling = "None", RCorrection = "Chao1", JackOver = FALSE, CEstimator = "ZhangHuang", CheckArguments = TRUE) 
{
  UseMethod("as.ProbaVector")
}


as.ProbaVector.data.frame  <-
function (x, Correction = "None", Unveiling = "None", RCorrection = "Chao1", JackOver = FALSE, CEstimator = "ZhangHuang", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()

  spD <- as.SpeciesDistribution(x)
  
  return(as.ProbaVector.numeric(spD, CheckArguments=FALSE))
}


as.ProbaVector.numeric <-
function (x, Correction = "None", Unveiling = "None", RCorrection = "Chao1", JackOver = FALSE, CEstimator = "ZhangHuang", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Try to save the names before applying as.vector
  spNames <- names(x)
  spD <- as.SpeciesDistribution(as.vector(x))
  if (length(spNames) == length(spD))
    names(spD) <- spNames

  if (Correction == "None") {
    spD <- spD/sum(spD)
  } else {
    # Integer abundances are required
    NsInt <- round(spD)
    if (any(abs(NsInt-spD) > sum(spD)*.Machine$double.eps)) warning("Integer abundance values are required to estimate community probabilities. Abundances have been rounded.")
    
    # Eliminate 0 and calculate elementary statistics
    Ns <- NsInt[NsInt > 0]
    S <- length(Ns)
    N <- sum(Ns)
    Ps <- Ns/N
    C <- Coverage(Ns, Estimator=CEstimator)
    if (Correction == "Chao2015" | Unveiling == "Chao2015") {
      Singletons <- sum(Ns==1)
      Doubletons <- sum(Ns==2)
      if (Doubletons==0) {
        Singletons <- max(Singletons - 1, 0)
        Doubletons <- 1
      }
      Tripletons <- max(sum(Ns==3), 1)
      # 1 minus sample coverage (i.e. Coverage Deficit) of order 2
      CD2 <- Doubletons / choose(N, 2) * ((N-2)*Doubletons / ((N-2)*Doubletons + 3*Tripletons))^2
    }
    
    # Estimate the number of unobserved species
    Sestimate <- ceiling(bcRichness(Ns, Correction=RCorrection, JackOver=JackOver))
    S0 <- Sestimate - S
    
    # Tune the probabilities of observed species
    if (C == 1) {
      # Sample coverage equal to 1: do not tune
      PsTuned <- Ps
    } else {
      PsTuned <- NA
      if (Correction == "ChaoShen") {
        PsTuned <- C*Ps
      }
      if (Correction == "Chao2013") {
        # Single parameter estimation, Chao et al. (2013)
        denominator <- sum(Ps*(1-Ps)^N)
        if (denominator == 0) {
          # N too big, denominator equals 0. Just multiply by C.
          PsTuned <- Ps*C
        } else {
          # General case
          lambda <- (1 - C)/denominator
          PsTuned <- Ps*(1 - lambda*(1-Ps)^N)
        }      
      } 
      if (Correction == "Chao2015")  {
        # Two parameters, Chao et al. (2015). 
        # Code inspired from JADE function DetAbu(), http://esapubs.org/archive/ecol/E096/107/JADE.R
        theta.solve <- function(theta){
          lambda <- (1-C) / sum(Ps * exp(-theta*Ns))
          out <- sum((Ps * (1 - lambda * exp(-theta*Ns)))^2) - sum(choose(Ns,2)/choose(N,2)) + CD2
          abs(out)
        }
        # Estimate theta. Set it to 1 if impossible
        theta <- tryCatch(stats::optimize(theta.solve, c(0,1))$min, error = function(e) {1})
        lambda <- (1-C) / sum(Ps * exp(-theta*Ns))
        PsTuned <- Ps * (1 - lambda * exp(-theta*Ns))
      }
      if (any(is.na(PsTuned))) {
        warning("Correction was not recognized")
        return (NA)
      }
    }
    names(PsTuned) <- names(spD[spD > 0])

    # Unobserved species
    if (S0) {
      if (Unveiling == "None") {
        spD <- PsTuned
      } else {
        Ps0 <- NA
        if (Unveiling == "geom") {
          # Code inspired from JADE function UndAbu(), http://esapubs.org/archive/ecol/E096/107/JADE.R
          if (S0 == 1) {
            # A single unobserved species
            Ps0 <- 1-C
          } else {
            r <- (1-C)^2/CD2
            i <- 1:S0
            beta.solve <- function(beta){
              out <- sum(beta^i)^2 / sum((beta^i)^2) - r
              abs(out)
            }
            beta <-  tryCatch(optimize(beta.solve, lower=(r-1)/(r+1), upper=1, tol=.Machine$double.eps)$min, error = function(e) {(r-1)/(r+1)})
            alpha <- (1-C) / sum(beta^i)
            Ps0 <- alpha * beta^i
            # Sometimes fails when the distribution is very uneven (sometimes r < 1) 
            # Then, go back to the uniform distribution
            if (any(Ps0 <= 0)) Unveiling <- "unif"
          }
        }      
        if (Unveiling == "unif") {
          # Add S0 unobserved species with equal probabilities
          Ps0 <- rep((1-sum(PsTuned))/S0, S0)
        }
        if (any(is.na(Ps0))) {
          warning("Unveiling method was not recognized")
          return(NA)
        } else {
          names(Ps0) <- paste("UnobsSp", 1:(length(Ps0)), sep="")
          spD <- c(PsTuned, Ps0)
        }         
      }
    } else {
      spD <- PsTuned
    }
    spD <- as.SpeciesDistribution(spD)
  }
  class(spD) <- c("ProbaVector", class(spD))
  return(spD)
}


as.ProbaVector.integer <-
function (x, Correction = "None", Unveiling = "None", RCorrection = "Chao1", JackOver = FALSE, CEstimator = "ZhangHuang", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return(as.ProbaVector.numeric(x, Correction=Correction, Unveiling=Unveiling, RCorrection=RCorrection, JackOver=JackOver, CEstimator=CEstimator, CheckArguments=FALSE))
}


is.ProbaVector <-
function (x) 
{
  inherits(x, "ProbaVector")
}


as.AbdVector <-
function (x, Round = TRUE)
{
  UseMethod("as.AbdVector")
}


as.AbdVector.data.frame <-
function (x, Round = TRUE) 
{
  # Try to save the names before applying as.vector
  spNames <- names(x)
  
  if (Round) {
    intx <- as.integer(as.matrix(round(x)))
    spD <- as.SpeciesDistribution(as.vector(intx))
  } else {      
    spD <- as.SpeciesDistribution(as.vector(x))
  }

  # Restore the names
  if (length(spNames) == length(spD))
    names(spD) <- spNames
  
  class(spD) <- c("AbdVector", class(spD))
  return(spD)
}


as.AbdVector.numeric <-
function (x, Round = TRUE) 
{
  # Try to save the names before applyinf as.vector
  spNames <- names(x)

  if (Round) {
    intx <- as.integer(round(x))
    spD <- as.SpeciesDistribution(as.vector(intx))
  } else {      
    spD <- as.SpeciesDistribution(as.vector(x))
  }
  
  # Restore the names
  if (length(spNames) == length(spD))
    names(spD) <- spNames
  
  class(spD) <- c("AbdVector", class(spD))
  return(spD)
}


as.AbdVector.integer <-
function (x, Round = TRUE) 
{
  return(as.AbdVector.numeric(x))
}


is.AbdVector <-
function (x) 
{
  inherits(x, "AbdVector")
}


plot.SpeciesDistribution <-
function(x, ..., Distribution = NULL, 
         type = "b", log = "y", main = NULL, xlab = "Rank", ylab = NULL) 
{
  # Eliminate zeros and sort
  Ns <- sort(x[x > 0], decreasing = TRUE)
  N <- sum(Ns)
  S <- length(Ns)
  
  # Prepare ylab
  if (is.null(ylab)) {
    if (is.ProbaVector(x)) {
      ylab <- "Probability"
    } else {
      ylab <- "Abundance" 
    }
  }
  
  graphics::plot(Ns, type=type, log=log, main=main, xlab=xlab, ylab=ylab, axes=FALSE, ...)
  # x axis ticks must start from 1
  graphics::axis(1, graphics::axTicks(1)+1)
  graphics::axis(2)
  graphics::box()
  
  if (!is.null(Distribution)) {
    if (Distribution == "lnorm") {
      # Fit a lognormal distribution
      mu <- mean(log(Ns))
      sigma <- stats::sd(log(Ns))
      ranks <- S*(1-stats::pnorm(log(Ns), mu, sigma))
      graphics::lines(ranks, Ns, col = "red")
      return(list(mu = mu, sigma = sigma))
    }
    if (Distribution == "geom") {
      lNs <- log(Ns)
      Rank <- 1:S
      reg <- stats::lm(lNs~Rank)
      graphics::lines(Rank, exp(reg$coefficients[1]+reg$coefficients[2]*Rank), col = "red")
      return(list(prob = as.numeric(-reg$coefficients[2])))
    }
    if (Distribution == "lseries") {
      # evaluate alpha
      alpha <- vegan::fisher.alpha(Ns)
      # May (1975) Ecology and evolution of communities, Harvard University Press
      sei <- function(t) exp(-t)/t
      ranks <- vapply(Ns, function(x) {
        n <- x * log(1 + alpha/N)
        f <- stats::integrate(sei, n, Inf)
        fv <- f[["value"]]
        return(alpha * fv)}
      , 0)
      graphics::lines(ranks, Ns, col = "red")
      return(list(alpha = alpha))
    }
    if (Distribution == "bstick") {
      f1 <- sort(cumsum(1/(S:1)), decreasing = TRUE)
      Rank <- 1:S
      Abundance <- N*f1/sum(f1)
      graphics::lines(Rank, Abundance, col = "red")
      return(list(max = max(Abundance)))
    }
    warning("The distribution to fit has not been recognized")
    return(NA)
  }
}
