#' Create a transition probability matrix for the states in the Hidden Markov Model
#' 
#' @details
#' Each row in the transition probability matrix sums to 1. For N regimes
#' we create N(N-1) paramters. 
#' 
#' @param tau parameter that is optimized 
#' @param numRegimes number of regimes to fit to the data
#' 
#' @export
#'  
.TransitionProb <- function(tau, numRegimes) {
  M <- matrix(tau, nrow = numRegimes, byrow = TRUE)
  eM <- exp(M); reM <- rowSums(eM)
  P <- cbind(eM/(1+reM), 1/reM)
  P
}

#' Initialize control paramters
#' 
#' @details
#' This is an internal function for initializing differnt control paramters for
#' the optimizer and setting up bounds for the paramters that need to be 
#' optimized
#' 
#' @param numRegimes number of regimes to fit to the data
#' @param tau lower and upper bounds for the tau paramter required to compute the
#'        transition probability matrix
#' @param theta lower and upper bounds for the theta paramter required to compute the
#'        DCC model for each regime#' 
#' @param ... additional control paramters required by DEoptim. Several control
#'              paramters are defualted if not passed
#'  
.optim.params <- function(numRegimes, tau = c(2,10), theta = c(0,1), ...) {
  
  lower.tau <- rep(tau[1], numRegimes*(numRegimes-1))
  lower.thetas <- rep(theta[1], 2*numRegimes); 
  lower <- c(lower.tau, lower.thetas)
  
  upper.tau <- rep(tau[2], numRegimes*(numRegimes-1))
  upper.thetas <- rep(theta[2], 2*numRegimes); 
  upper <- c(upper.tau, upper.thetas)
  
  add.args <- list(...)
  
  if(!"NP" %in% names(add.args)) 
    add.args[["NP"]] <- 10*(numRegimes*(numRegimes + 1))
  
  NP <- add.args[["NP"]]
  
  if(!"initialpop" %in% names(add.args)) {
    pop1 <- maximinLHS(NP,numRegimes*(numRegimes-1))
    pop1 <- pop1*(upper.tau - lower.tau) + lower.tau
    
    pop2 <- maximinLHS(NP, 2*numRegimes)
    pop2 <- pop2*(upper.thetas - lower.thetas) + lower.thetas
    add.args[["initialpop"]] <- cbind(pop1, pop2)
  }
  
  if(!"parallelType" %in% names(add.args)) 
    add.args[["parallelType"]] <- 1
  
  if(!"itermax" %in% names(add.args)) 
    add.args[["itermax"]] <- 100
  
  if(!"trace" %in% names(add.args)) 
    add.args[["trace"]] <- 20
  
  list(lower = lower, upper = upper, add.args = add.args )
}

#' Calculate the initial value of the correlation
#' 
#' @details
#' This is an internal function for calculating the intial estimate of 
#' correlation. The entire data set is partioned into parts based on the 
#' weight vector and correlation is computed for each part.
#' 
#' @param R xts object of returns data
#' @param numRegimes number of regimes to fit to the data
#' @param w weight vector to select entries from returns
#' 
#' @return Q0 is a list of correlation estimates for each regime
#' 
.initQ <- function(R, numRegimes, w) {

  w <- w/sum(w)
  size <- ceiling(nrow(R)*w) 
  size[numRegimes] <- size[numRegimes] - (sum(size) - nrow(R))
  starts <- 1 + c(0, head(size, -1))
  
  parts <- lapply(1:length(size), function(i) 
    seq(from = starts[i], length.out = size[i]))

  Q0 <- lapply(parts, function(indices) cor(R[indices,]))
  names(Q0) <- 1:numRegimes
  Q0
}


#' Calculates the unconditional state probabilities
#' 
#' @details
#' This is an internal function for calculating the intial steady state
#' probabilities of being in a state.
#' 
#' @param P esitmate of the transition probability matrix
#' 
#' @export
#' 
.initFilterProb <- function(P) {
  e <- eigen(t(P))
  e$vectors[,1]/sum(e$vectors[,1])
}

#' Calculate the log-liklihood using Expectation Maximization
#' 
#' @details
#' This is an internal function. Given a set of paramters it returns the 
#' log-liklihood using Expectation Maximization procedure described in (Lee, 2010).
#' It also returns the final covariance matrix and final probability of being
#' in each state. 
#' 
#' @param taus parameters required to fit a transition probability matrix
#' @param theta1 parameters required to fit a DCC model
#' @param theta2 parameters required to fit a DCC model
#' @param returns time series of returns data
#' @param numRegimes number of regimes to fit to the data
#' @param Q initial value of the correlation
#' @param Q.bar constant unconditional correlation matrix
#' 
#' @export
#'  
.helper.loglik <- function(taus, theta1, theta2, returns, Sigma, 
                           numRegimes, Q, Q.bar) {
  
  rstd <- returns/Sigma
  
  P <- .TransitionProb(taus, numRegimes); T <- nrow(returns)
  filtprob <- .initFilterProb(P); loglik <- 0
  likelihood <- rep(NA, numRegimes) 

  Cov <- replicate(T,list())
  
  prob_t <- matrix(NA, nrow = T, ncol = numRegimes )
  prob_t[1, ] <- filtprob
  
  for(t in 2:T) {
    
    filtprob <- P %*% filtprob; D <- diag(Sigma[t,])
    
    for ( i in 1:numRegimes) {
      Q[[i]] <- (1 - theta1[i] - theta2[i])*Q.bar + 
        theta1[i]*t(rstd[(t-1),,drop=FALSE]) %*% rstd[(t-1),,drop=FALSE] +
        theta2[i]*Q[[i]]

      M <- diag(diag(Q[[i]])^(-0.5)); C <- M %*% Q[[i]] %*% M
      H <- D %*% C %*% D
      likelihood[i] <- dmvnorm(returns[t,], sigma = H)
      
      rownames(H) <- colnames(H) <- colnames(returns)
      Cov[[t]][[i]] <- H
    }
    
    mixLik <- as.numeric(
      matrix(rep(1, numRegimes), nrow = 1) %*% (filtprob*likelihood))
    
    if(mixLik == 0) loglik <- -1e6
    else {
      filtprob <- (filtprob*likelihood)/mixLik
      prob_t[t, ] <- filtprob
      loglik <- loglik + log(mixLik)
    }
  }
  
  list(loglik=loglik, Cov = Cov, filtProb = prob_t)
}

#' Calculate the log-liklihood
#' 
#' @details
#' This is an internal function that is used as an objective by the DEoptim
#' optimizer. Given a set of paramters it returns the negative of the 
#' log-liklihood.
#' 
#' @param params contains paramters for estimating the transition probability 
#'                matrix and parameters required for fitting a DCC model to each
#'                regime.
#' @param returns xts obect of returns
#' @param Sigma   fitted volatility for each asset return
#' @param numRegimes number of regimes to fit to the data
#' @param Q initial value of the correlation
#' @param Q.bar constant unconditional correlation matrix
#' 
#'  
.isdcc.loglik <- function(params, returns, Sigma, numRegimes, Q, Q.bar) {
  
  taus <- params[1:(numRegimes*(numRegimes-1))]; 
  
  theta1 <- params[(1 + numRegimes*(numRegimes-1)):(numRegimes*numRegimes)]; 
  theta2 <- params[(1+numRegimes*numRegimes):length(params)]
  
  if (any(colSums(rbind(theta1, theta2)) >=1)) return(1e6)
  
  results <- .helper.loglik(taus, theta1, theta2, returns, Sigma, 
                            numRegimes, Q, Q.bar)
  -results$loglik
}

#' Fit an Independent Regime Switching Model
#' 
#' @references 
#' Lee, H.-T. (2010). Regime switching correlation hedging. Journal of Banking &
#  Finance, 34, 2728-2741

#' @details
#' This method takes in returns data and the number of regimes and fits 
#' sepearate covariances to each regime using the Expectation Maximization
#' algorithm decribed in (Lee, 2010). IS-DCC model avoids the path dependency
#' problem observed in other regime switching models and makes the solution more
#' tractable by running a separate DCC process for each regime.
#' 
#' Fitting the IS-DCC model to data corresponds works in two steps. In the first
#' step a time varying univariate volatility process, GARCH(1,1) is fitted to 
#' each time series. In the second step DCC parameters for each state
#' are estimated along with the transition probabilities corresponding to the
#' Hidden Markov model. This is done by maximising the log-likelihood of 
#' observing the residuals
#' 
#' @importFrom fGarch garchFit coef
#' @importFrom DEoptim DEoptim
#' @importFrom mvtnorm dmvnorm
#' @importFrom lhs maximinLHS
#' 
#' @param  R xts object of asset returns
#' @param  numRegimes number of regimes to fit to the data
#' @param  transMatbounds bounds on the parameter tau as described in (Lee, 2010).
#'          Each paramter is defaulted to lie in the range (2,10)
#' @param  dccBounds  bounds on the paramter theta as described in (Lee, 2010).
#'          Each paramter is defaulted to lie in the range (0,1)
#' @param  w proportion of entries to consider in initializing correlation for 
#'          for each regime. It is defualted to split data equally across 
#'          all regimes
#' @param  ... addition control paramters that can be passed to the control object
#'        in DEoptim 
#'        
#' @examples 
#' \dontrun{
#'  data("largereturn")
#'  model <- isdccfit(largesymdata, numRegimes=2, maxiter=50, paralletType=0)  
#' }
#' 
#' 
#' @export
#' 
#' @author Rohit Arora
#' 
isdccfit <- function(R, numRegimes=NA, 
                     transMatbounds = c(2,10), dccBounds = c(0,1),
                     w = NA, ...) {
  
  if (!is.xts(R)) stop("Only xts object required")
  if (!(abs(numRegimes - round(numRegimes)) < .Machine$double.eps)) 
    stop("Regimes must be an integer")
  if (length(transMatbounds) != 2 || any(is.na(transMatbounds)))
    stop("Transition Matrix bounds are incorrect")
  if ( transMatbounds[2]  < transMatbounds[1] )
    stop("Transition Matrix bounds are invalid")
  if (length(dccBounds) != 2 || any(is.na(dccBounds)))
    stop("DCC paremeter bounds are incorrect")
  if ( dccBounds[2]  < dccBounds[1] )
    stop("DCC parameter bounds are invalid")
  if ( dccBounds[1]  < 0 )
    stop("DCC parameter bounds must be positive")
  
  if (any(is.na(w))) w <- rep(1/numRegimes, numRegimes)  
  
  
  garchFit <- apply(R, 2, function(data)
    garchFit(~ garch(1,1), data = data, trace = FALSE))  
  
  sigma_t <- lapply(garchFit, function(garch) sqrt(garch@h.t))
  sigma_t <- do.call(cbind,sigma_t)
  
  stdReturn_t <- R/sigma_t
  
  Q0 <- .initQ(stdReturn_t, numRegimes, w);  Q.bar <- cov(stdReturn_t)
  
  initParams <- .optim.params(numRegimes, tau = transMatbounds, 
                              theta = dccBounds, ...)
  
  isParallel <- (initParams$add.args$parallelType > 0)
  
  cl <- NULL
  if (isParallel) {
    cl <- makeCluster(detectCores())
    clusterExport(cl, varlist = c(".helper.loglik", ".TransitionProb", 
                               ".initFilterProb", "dmvnorm"))
    registerDoParallel(cl)  
  }

  fit <- DEoptim(fn = .isdcc.loglik, 
                 lower=initParams$lower, upper=initParams$upper,
                 returns = R, Sigma = sigma_t, numRegimes = numRegimes, 
                 Q=Q0, Q.bar = Q.bar, control = initParams$add.args)
  
  if (isParallel)  { 
    stopCluster(cl)
    registerDoSEQ()
  }
  
  params <- fit$optim$bestmem
  
  garchParams <- do.call(rbind, lapply(garchFit, coef))
  
  taus <- params[1:(numRegimes*(numRegimes-1))]
  names(taus) <- sapply(1:length(taus), function(i) paste("tau",i,sep=""))
  
  theta1 <- params[(1 + numRegimes*(numRegimes-1)):(numRegimes*numRegimes)];
  names(theta1) <- sapply(1:length(theta1), function(i) paste("theta1_",i,sep=""))
  
  theta2 <- params[(1+numRegimes*numRegimes):length(params)]
  names(theta2) <- sapply(1:length(theta2), function(i) paste("theta2_",i,sep=""))
  
  result <- 
    .helper.loglik(taus, theta1, theta2, R, sigma_t, 
                   numRegimes, Q0, Q.bar)
  
  dates <- as.character(index(R))
  prob <- result$filtProb; rownames(prob) <- dates
  cov <- result$Cov; names(cov) <- dates
  
  fit <- list(logLik = result$loglik, filtProb = prob, 
       param = list(garch = garchParams, ISDCCParams = c(taus, theta1, theta2)),
       cov = cov)
  
  class(fit) <- "isdcc"
  fit
}

#'Implied State plot
#' 
#' @details
#' Plot implied states using the fitted Independent Switching DCC model
#' 
#' @param x model of the type isdcc obtained by fitting an IS-DCC model to the data
#' @param  y type of plot. takes values 1/2. 1 = Implied States, 2 = Smoothed Proabability
#' @param ... additional arguments unused
#' @author Rohit Arora
#' @examples 
#' \dontrun{
#'  data("largereturn")
#'  model <- isdccfit(largesymdata, numRegimes=2)
#'  plot(model)
#' }
#' 
#' @method plot isdcc
#' @export
#' 
plot.isdcc <- function(x, y = c(1,2), ...){
  
  which.plot <- y[1]
  
  prob <- x$filtProb
  states <- apply(prob, 1, which.max)
  dates <- rownames(prob)
  
  df <- data.frame(dates = 1:length(states), states = states, 
                   prob = apply(prob, 1, max))

  year.dates <- format(as.Date(dates), "%Y")
  ind <- sapply(unique(year.dates), function(val) 
    which.max(year.dates == val))
  ind.ind <- seq.int(1, length(ind), length.out = min(10,length(ind)))
  ind <- ind[ind.ind]

  p <- ggplot(data = df, aes(x = dates)) + 
    scale_x_continuous(breaks = ind, labels=year.dates[ind])

  p <- 
    if (which.plot == 1) 
      p + geom_step(aes(y = states),direction = "hv") + ylab("States") + 
      scale_y_continuous(breaks = c(min(states),max(states))) + 
      ggtitle("Implied States")
  else if (which.plot == 2) 
    p + geom_line(aes(y = prob, color = as.factor(states))) + 
    ylab("Probability") + ggtitle("Smoothed Probability") + 
    scale_colour_discrete(name = "States")

  p <- p + theme(plot.title = element_text(size = 20, face = "bold", vjust = 1),
             axis.title=element_text(size=14,face="bold"))
  
  options(warn = -1)
  print(p)
  options(warn = 0)
  p
}