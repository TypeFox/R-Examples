# this is to get around issues with ggplot2 code not passing package checks
if(getRversion() >= "2.15.1")  utils::globalVariables(c("alternative"))


################################################################################
# E-Step and M-Step
################################################################################

densities.to.utilities <- function(rum.densities, Data, race = FALSE, utilities.per.agent = 20) {
  
  GIBBSCONVERGENCE <- 20 # number of iterations Gibbs sampler must go through before we declare convergence
  
  RUFs <- function(i, x) {
    n <- length(rum.densities[[i]])
    rum.densities[[i]][max(min(round(x * (n - 1)) + 1, n), 1)]
  }
  
  n <- nrow(Data)
  
  get.utilities <- function(row) {      
    missing <- sum(as.numeric(Data[row, ]) == 0)
    if(missing == 0 || race == TRUE) {
      race.true.or.no.missing(row)
    } else {
      partial.rank(row)
    }
    
  }
  
  race.true.or.no.missing <- function(row){
    ranks <- as.numeric(Data[row, ])
    
    # this allows for "race" cases
    ranks <- ranks[ranks > 0]
    m <- length(ranks)
    
    # initializing utilities
    utilities <- c(1, sort(runif(m), decreasing = TRUE), 0)
    
    output <- data.frame()
    
    # run Gibbs sampler to sample utilities
    for(j in 1:(GIBBSCONVERGENCE + 2 * utilities.per.agent)) { # number of iterations of Gibbs sampler
      for(i in scramble(2:m + 1)) # order in which we want to sample the alternatives
        utilities[i] <- uni.slice(utilities[i], function(x) log(RUFs(ranks[i - 1], x)), upper = utilities[i - 1], lower = utilities[i + 1])
      if(j > GIBBSCONVERGENCE & j %% 2 == 0)
        output <- rbind(output, data.frame(alternative = ranks, utilities = utilities[c(-1, -(m+2))], agent = paste0(row, "-", j - GIBBSCONVERGENCE)))
    }
    
    output
  }
  
  partial.rank <- function(row){
    ranks <- as.numeric(Data[row, ])
    
    m <- length(ranks)
    ranks <- ranks[ranks > 0]
    present <- length(ranks)
    which.missing <- (1:m)[!(1:m %in% ranks)]
    
    # initializing utilities
    utilities <- c(1, sort(runif(m), decreasing = TRUE), 0)
    
    output <- data.frame()
    # run Gibbs sampler to sample utilities
    for(j in 1:(GIBBSCONVERGENCE + 2 * utilities.per.agent)) { # number of iterations of Gibbs sampler
      for(i in scramble(1:present + 1))
        utilities[i] <- uni.slice(utilities[i], function(x) log(RUFs(ranks[i - 1], x)), upper = utilities[i - 1], lower = max(utilities[i + 1], utilities[(present+1):m + 1]))
      for(i in (present+1):m + 1)
        utilities[i] <- uni.slice(utilities[i], function(x) log(RUFs(which.missing[i - present - 1], x)), upper = utilities[present + 1], lower = 0)
      if(j > GIBBSCONVERGENCE & j %% 2 == 0)
        output <- rbind(output, data.frame(alternative = c(ranks, which.missing), utilities = utilities[c(-1, -(m+2))], agent = paste0(row, "-", (j - GIBBSCONVERGENCE) / 2 )))
    }
    
    output
  }
  
  ldply(1:n, get.utilities)
}

utilities.to.densities <- function(u, m, bw = 0.025) {
  gety <- function(alternative) {
  	if(sum(u$alternative == alternative) < 2) rep(1, 512)
  	else {density(u$utilities[u$alternative == alternative], from = 0, to = 1, bw = bw)$y}
  }
  llply(1:m, gety)
}
#
################################################################################
# Running Algorithm
################################################################################

#' Nonparametric RUM Estimator
#' 
#' Given rank data (full, top partial, or sub partial), this function returns an
#' inference object that fits nonparametric latent utilties on the rank data.
#' 
#' @param Data full, top partial, or sub partial rank data
#' @param m number of alternatives
#' @param iter number of EM iterations to run
#' @param bw bandwidth, or smoothing parameter for KDE
#' @param utilities.per.agent Number of utility vector samples that we get per
#' agent. More generally gives a more accurate estimate
#' @param race TRUE if data is sub partial, FALSE (default) if not
#' @export
#' @examples
#' data(Data.Test)
#' Estimation.RUM.Nonparametric(Data.Test, m = 5, iter = 3)
Estimation.RUM.Nonparametric <- function(Data, m, iter = 10, bw = 0.025, utilities.per.agent = 20, race = FALSE) {
  t0 <- proc.time()
  x.star <- seq(0, 1, len = 512)
  
  # initializer
  rum.densities <- replicate(m, rep(1, length(x.star)), simplify = FALSE)
  
  for(i in 1:iter) {
    print(paste0("Iteration ", i, "/", iter))
    utilities     <- densities.to.utilities(rum.densities, Data, race = race, utilities.per.agent) # MC E-Step
    rum.densities <- utilities.to.densities(utilities, m, bw = bw) # Variational M-Step
  }
  ordering <- order(-ddply(utilities, .(alternative), summarize, mean = mean(utilities))$mean)
  
  Time <- proc.time() - t0
  
  # reshape utilities to usable format
  utilities.by.agent <- reshape(utilities, v.names = "utilities", timevar = "alternative", idvar = "agent", direction = "wide", sep = "")
  utilities.by.agent$agent <- NULL
  
  # make the columns if they don't exist
  unames <- paste("utilities", 1:m, sep="")
  utilities.by.agent[, unames[!(unames %in% names(utilities.by.agent))]] <- rep(NA, nrow(utilities.by.agent))
  
  # order by alternative
  utilities.by.agent <- as.matrix(utilities.by.agent[, unames])
  
  list(m = m, order = ordering, iter = iter, bw = bw, Time = Time, utilities = utilities, utilities.by.agent = utilities.by.agent, rum.densities = rum.densities)
}

get.nomodel.KL <- function(Data.train, Data.test) {
  Data.train.pairs   <- Breaking(Data.train, method = "full")
  Data.test.pairs    <- Breaking(Data.test, method = "full")
  Data.train.preferences  <- generateC(Data.train.pairs, m)
  Data.test.preferences   <- generateC(Data.test.pairs, m)
  KL(Data.test.preferences, Data.train.preferences)
}

################################################################################
# Generate pairwise matrix
################################################################################

#' Generate pairwise matrix for an NPRUM model
#' 
#' Generates a matrix where entry i, j is the estimated probabiltiy that
#' alternative i beats alternative j
#' 
#' @param Estimate fitted NPRUM object
#' @param bw bandwidth used for generating the pairwise probabilites
#' @export
#' @examples
#' data(Data.Test)
#' Estimate <- Estimation.RUM.Nonparametric(Data.Test, m = 5, iter = 3)
#' generateC.model.Nonparametric(Estimate)
generateC.model.Nonparametric <- function(Estimate, bw = 0.1) {
  sample.ranks <- Generate.NPRUM.Data(Estimate, nrow(Estimate$utilities.by.agent), bw = bw)
  sample.ranks.pairs <- Breaking(sample.ranks, method="full")
  generateC(sample.ranks.pairs, Estimate$m)
}

################################################################################
# Likelihood Function
################################################################################
#' Calculate Likelihood for the nonparametric model
#' 
#' Computes likelihood in the case that we assume no correlation structure
#' 
#' @param Data full, top partial, or subpartial data
#' @param Estimate fitted NPRUM object
#' @param race indicator that the data is from subpartial data
#' @export
#' @examples
#' data(Data.Test)
#' Estimate <- Estimation.RUM.Nonparametric(Data.Test, m = 5, iter = 3)
#' Likelihood.Nonparametric(Data.Test, Estimate)
Likelihood.Nonparametric <- function(Data, Estimate, race = FALSE) {
  n <- nrow(Data)
  m <- ncol(Data)
  rum.densities <- replicate(m, list(NA))
  
  normalize.density <- function(density.vector) {
    n <- length(density.vector)
    density.vector[1] <- density.vector[1] * .5
    density.vector[n] <- density.vector[n] * .5
    multiplier <- (n - 1) / sum(density.vector)
    multiplier * density.vector
  }
  
  for(i in 1:m) rum.densities[[i]] <- normalize.density(Estimate$rum.densities[[i]])
  
  buckets <- length(rum.densities[[1]])
  np.pdf <- function(i, x) rum.densities[[i]][max(min(round(x * (buckets - 1)) + 1, buckets), 1)]
  x.star <- 0:(buckets - 1)/(buckets - 1)
  
  row_num_to_ll <- function(i) {
    #print(i)
    CDF = matrix(1,1,buckets)
    mj <- sum(Data[i, ] > 0)

    if(mj < m & !race) {
      for(jt in setdiff(1:m, Data[i, 1:mj]))
        CDF <- sapply(x.star, function(a) np.pdf(jt, a)) * CDF
    }
    
    for(j in mj:1) {
      PDF <- sapply(x.star, function(a) np.pdf(Data[i, j], a)) * CDF
      CDF <- cumsum(PDF) / buckets
    }
    log(CDF[buckets])

  }
  
  sum(sapply(1:n, row_num_to_ll))
}

################################################################################
# Data generation
################################################################################

#' Generate data from an NPRUM model
#' 
#' This is useful for performing inference tasks for NPRUM
#' 
#' @param Estimate fitted NPRUM object
#' @param n number of agents that we want in our sample
#' @param bw smoothing parameter to use when sampling data
#' @export
#' @examples
#' Data.Tiny <- matrix(c(1, 2, 3, 3, 2, 1, 1, 2, 3), ncol = 3, byrow = TRUE)
#' Estimate <- Estimation.RUM.Nonparametric(Data.Tiny, m = 3, iter = 3)
#' Generate.NPRUM.Data(Estimate, 3, bw = 0.1)
Generate.NPRUM.Data <- function(Estimate, n, bw = 0.1) {
  u <- -Estimate$utilities.by.agent + rnorm(n = length(Estimate$utilities.by.agent), sd = bw)
  if(n < nrow(u)) {
    u <- u[sample(1:nrow(u), n, replace = FALSE), ]
  } else if(n > nrow(u)) {
    u <- u[sample(1:nrow(u), n, replace = TRUE), ]
  }
  t(apply(u, 1, function(x) order(x, na.last = NA)))
}