#' Estimate discrete choice model with random parameters
#' 
#' Estimation of discrete choice models such as Binary (logit and probit), 
#' Poisson and Ordered (logit and probit) model with random coefficients for cross-sectional and panel data using simulated maximum likelihood.
#' 
#' @name Rchoice
#' @param formula a symbolic description of the model to be estimated. The \code{formula} consists in two parts. The first one is reserved for standard variables with fixed and random parameters. The second one is reserved for variables that enter in the mean of the random parameters. See for example \code{\link{rFormula}},
#' @param new an updated formula for the update method,
#' @param data the data. It may be a \code{pdata.frame} object or an ordinary \code{data.frame}, 
#' @param subset an optional vector specifying a subset of observations,
#' @param weights an optional vector of weigths,
#' @param na.action a function wich indicated what should happen when the data
#' contains \code{NA}'s,
#' @param start a vector of starting values,
#' @param family the distribution to be used. It might be \code{family = binomial("probit")} for a Probit Model, \code{family = binomial("logit")} for a Logit model, \code{family = ordinal("probit")} for an Ordered Probit Model, \code{family = ordinal("logit")} for a Ordered Logit Model for an Ordered Logit Model, and  \code{family = "poisson"} for a Poisson Model, 
#' @param ranp a named vector whose names are the random parameters and values the distribution:
#' "\code{n}" for normal, "\code{ln}" for log-normal, "\code{cn}" for truncated normal, "\code{u}" for uniform, "\code{t}" for triangular, "\code{sb}" for Johnson Sb,
#' @param R the number of draws if \code{ranp} is not \code{NULL},
#' @param haltons only relevant if \code{ranp} is not \code{NULL}. If not \code{NULL}, halton sequence is used
#' instead of pseudo-random numbers. If \code{haltons=NA}, some default values are used for
#' the prime of the sequence and for the number of element dropped. Otherwise, \code{haltons} should
#' be a list with elements \code{prime} and \code{drop},
#' @param seed  the seed for the pseudo-random draws. This is only relevant if \code{haltons = NULL},
#' @param correlation only relevant if \code{ranp} is not \code{NULL}. If \code{TRUE}, the correlation between random parameters is taken into account,
#' @param panel if \code{TRUE} a panel data model is estimated,
#' @param index a string indicating the `id' for individuals in the data. This argument is not required if data is a \code{pdata.frame} object,
#' @param mvar only valid if \code{ranp} is not \code{NULL}. This is a named list, where the names correspond to the variables with random parameters, and the values correspond to the variables that enter in the mean of each random parameters,
#' @param print.init if \code{TRUE}, the initial values for the optimization procedure are printed,
#' @param init.ran initial values for standard deviation of random parameters. Default is 0.1,
#' @param gradient if \code{FALSE}, numerical gradients are used for the optimization procedure of models with random parameters,
#' @param digits number of digits,
#' @param width width,
#' @param x,object and object of class \code{Rchoice},
#' @param ... further arguments passed to \code{\link[maxLik]{maxLik}},
#' @export
#' @details
#' The models are estimated using the \code{maxLik} function from \code{\link[maxLik]{maxLik}} package.
#' 
#' 
#'  If \code{ranp} is not \code{NULL}, the random parameter  model is estimated.   
#'  A random parameter model or random coefficient models permits regression parameter to 
#'  vary across individuals according to some distribution. A fully parametric 
#'  random parameter model specifies the latent variable  \eqn{y^{*}} conditional on regressors
#'  \eqn{x} and given parameters \eqn{\beta_i} to have conditional density \eqn{f(y|x, \beta_i)} where
#'  \eqn{\beta_i} are iid with density \eqn{g(\beta_i|\theta_i)}. The density is assumed a priori by the user by the argument
#'  \code{ranp}. If the parameters are assumed to be normally distributed \eqn{\beta_i ~ N(\beta, \Sigma)}, then the random parameter are constructed as: \deqn{\beta_{ir}=\beta+L\omega_{ir}} where \eqn{LL'=\Sigma} and \eqn{\omega_{ir}} is the {r}-th draw from standard normal distribution for individual \eqn{i}. 
#'  
#'  
#'  Once the model is specified by the argument \code{family}, the model is estimated using 
#'  Simulated Maximum Likelihood (SML). The probabilities, given by \eqn{f(y|x, \beta_i)}, are simulated using \code{R} pseudo-draws if \code{halton=NULL} or \code{R} halton draws if \code{halton = NA}. The user can also specified the primes and the number of dropped elements for the halton draws. For example, if the model consists of two random parameters, the user can specify \code{haltons = list("prime" = c(2, 3), "drop" = c(11, 11))}. 
#'  
#'  
#'  A random parameter hierarchical model can be estimated by including heterogeneity in the mean of the 
#'  random parameters: \deqn{\beta_{ir}=\beta+\pi's_i+L\omega_{ir}} \pkg{Rchoice} manages the variables in the hierarchical model by the \code{formula} object: all the hierarchical variables (\eqn{s_i}) are included after the \code{|} symbol. The argument \code{mvar} indicate which variables enter in each random parameter. See examples below
#' 
#' @return An object of class ``\code{Rchoice}'', a list elements:
#' \item{coefficients}{the named vector of coefficients,}
#' \item{family}{type of model,}
#' \item{link}{distribution of the errors,}
#' \item{logLik}{a set of values of the maximum likelihood procedure,}   
#' \item{mf}{the model framed used,} 
#' \item{formula}{the formula (a Formula object),}
#' \item{time}{\code{proc.time()} minus the start time,}
#' \item{freq}{frequency of dependent variable,}
#' \item{draws}{type of draws used,}
#' \item{R.model}{\code{TRUE} if a random parameter model is fitted,}
#' \item{R}{number of draws used,}
#' \item{bi}{an array of dimension \eqn{N \times R \times K} with the individual parameters,}
#' \item{Qir}{matrix of dimension \eqn{N \times R} representing \eqn{P_{ir}/\sum_r P_{ir}},}
#' \item{ranp}{vector indicating the variables with random parameters and their distribution,}
#' \item{probabilities}{the fitted probabilities for each individuals,}
#' \item{residuals}{the residuals,}
#' \item{call}{the matched call.}
#' @aliases ordinal   
#' @author Mauricio Sarrias \email{msarrias86@@gmail.com}
#' @seealso \code{\link[Rchoice]{plot.Rchoice}}, \code{\link[Rchoice]{effect.Rchoice}}
#' @examples
#' ## Probit model
#' data("Workmroz")
#' probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc,  
#'                  data = Workmroz, family = binomial('probit'))
#' summary(probit)
#' 
#' ## Poisson model
#' data("Articles")
#' poisson <- Rchoice(art ~ fem + mar + kid5 + phd + ment, data = Articles, family = poisson)
#' summary(poisson)
#' 
#' ## Ordered probit model
#' data("Health")
#' oprobit <- Rchoice(newhsat ~ age + educ + hhinc + married + hhkids, 
#' data = Health, family = ordinal('probit'), subset = year == 1988)
#' summary(oprobit)
#' 
#' ## Poisson Model with Random Parameters
#' \dontrun{
#' poisson.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment, 
#'                        data = Articles,  family = poisson,
#'                        ranp = c(kid5 = "n", phd = "n", ment = "n"))
#' summary(poisson.ran)                      
#' 
#' ## Poisson Model with Correlated Random Parameters
#' poissonc.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment, 
#'                        data = Articles, 
#'                        ranp = c(kid5 = "n", phd = "n", ment = "n"), 
#'                        family = poisson, 
#'                        correlation =  TRUE)
#' summary(poissonc.ran)
#' 
#' ## Hierarchical Poisson Model
#' poissonH.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment | fem + phd,
#'                        data = Articles,
#'                        ranp = c(kid5 = "n", phd = "n", ment = "n"),
#'                        mvar = list(phd = c("fem"), ment = c("fem", "phd")),
#'                        family = poisson,
#'                        R = 10)
#' summary(poissonH.ran)
#' 
#' ## Probit Model with Random Effects and Random Parameters
#' data('Unions', package = 'pglm')
#' Unions$lwage <- log(Unions$wage)
#' union.ran <- Rchoice(union ~ age + exper + rural + lwage,
#'                      data = Unions[1:2000, ],
#'                      family = binomial('probit'),
#'                      ranp = c(constant = "n", lwage = "t"),
#'                      R = 10,
#'                      panel = TRUE,
#'                      index = "id",
#'                      print.init = TRUE)
#' summary(union.ran)
#' 
#' ## Ordered Probit Model with Random Effects and Random Parameters
#'oprobit.ran <- Rchoice(newhsat ~ age + educ + married + hhkids + linc,
#'                      data = Health[1:2000, ],
#'                      family = ordinal('probit'),
#'                      ranp = c(constant = "n", hhkids = "n", linc = "n"),
#'                      panel = TRUE,
#'                      index = "id",
#'                      R = 100,
#'                      print.init = TRUE)
#'summary(oprobit.ran)
#'} 
#'
#' 
#' @references
#' Greene, W. H. (2012). Econometric Analysis. 7 edition. Prentice Hall.
#' 
#' Train, K. (2009). Discrete Choice Methods with Simulation. Cambridge university press.
#' @import maxLik Formula
#' @importFrom plm pdata.frame
Rchoice <- function(formula, data, subset, weights, na.action, family,
                    start = NULL, ranp = NULL, R = 40, haltons = NA, 
                    seed = 10, correlation = FALSE, panel = FALSE, index = NULL,
                    mvar = NULL, print.init = FALSE, init.ran = 0.1, 
                    gradient = TRUE, ...){
  ####################
  # 1) Check arguments
  ####################
  start.time <- proc.time()
  callT <- match.call(expand.dots = TRUE)
  callF <- match.call(expand.dots = FALSE)
  
  ## version 2.0 uses rFormula
  formula <- callF$formula <- rFormula(formula)
  nframe  <- length(sys.calls())
  
  ## family
  if (is.character(family)){
    if (family %in% c("ordprobit", "ordlogit")){
      if (family == "ordprobit") family <- list(family = "ordinal", link = "probit")
      if (family == "ordlogit")  family <- list(family = "ordinal", link = "logit")
    }
    else  family <- get(family, mode = "function")
  }
  if (is.function(family)) family <- family()
  link   <- family$link
  family <- family$family
  
  ## check whether the model has random parameters and is hierarchical
  R.model <- !is.null(ranp)
  Hier <- is.hierarchical(formula)
  if (Hier){
    if (!R.model) stop('Hierarchical model needs ranp to be specified')
    if (is.null(mvar)) stop("mvar argument is NULL")
    if (!is.list(mvar)) stop("mvar is not a list")
    rvar <- names(mvar)[! (names(mvar) %in% names(ranp))]
    if (length(rvar) > 0){
      udstr <- paste("The following variables are not specified in the argument ranp:", 
                     paste(unique(rvar), collapse = ", "))
      stop(udstr)
    }
  }
  if (R.model){
    if (is.null(callT$method)) callT$method <- 'bfgs'
    if (is.null(callT$iterlim)) callT$iterlim <- 2000
  } else{
    if (is.null(callT$method) && (family == "ordinal")) callT$method <- 'bfgs'
    if (is.null(callT$method)) callT$method <- 'nr'
  }
  
  ## check whether is panel. If it is, check data is pdata.frame.
  ## if not create it.
  if (!is.null(index) && !panel){
    panel <- TRUE
    warning("panel was FALSE, while index was not NULL. Panel is now TRUE")
  }
  if (panel){
    if (inherits(data, "pdata.frame") && !is.null(index))
      warning("The index argument is ignored because data is a pdata.frame")
    if (!inherits(data, "pdata.frame") && is.null(index)) 
      stop("The data is not pdata.frame and index is needed")
    if (!inherits(data, "pdata.frame")) data <- pdata.frame(data, index)
  }
  
  ## for the Ordered Model, we need the constant
  if (family == "ordinal" && !has.intercept(formula, rhs = 1)){
    formula  <- callF$formula <- update(formula, ~. + 1)
    warning("Ordinal model need a constant ... updating formula")
  } 
  ####################
  # 2) Model Frame
  ####################
  mf <- callT
  m  <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf$data <- data
  mf <- eval(mf, parent.frame())
  
  ## check Panel
  if (panel){
    id <- attr(mf, "index")[[1]]
    if (is.null(id)) stop("No individual index")
  } else {
    id <- NULL
  }
  ##########################################
  # 3) Extract the elements of the model
  ##########################################
  if (Hier){
    S <- model.matrix(formula, mf , rhs = 2)
    for (i in 1:length(mvar)){
      rvar <- mvar[[i]][! (mvar[[i]] %in% colnames(S))]
      if(length(rvar) > 0){
        udstr <- paste("The following variables are not specified in the formula: ", 
                       paste(unique(rvar), collapse = ", "))
        stop(udstr)
      }
    }
  } else S <- NULL 
  y <- model.response(mf)
  X <- model.matrix(formula, mf, rhs = 1)
  if (has.intercept(formula, rhs = 1)){
    namesX <- colnames(X)
    namesX[1L] <- "constant"
    colnames(X) <- namesX
  } 
  freq <- table(y)
  
  ## Weights
  weights <- model.weights(mf)
  
  # Other warnings
  if (family == "ordinal") {
    y <- as.factor(y)
    J <- length(levels(y))
    if (J < 3) stop("The alternatives must be >=3 in y")
  }
  if (family == "binomial"){
    if (!all(y %in% c( 0, 1, TRUE, FALSE))){
      stop( "all dependent variables must be either 0, 1, TRUE, or FALSE")
    }
    if (!is.numeric(y)) y <- as.numeric(y)  
  }
  ########################################################
  #   4) Initial Values
  ########################################################
  
  ## Names of thresholds if ordered model
  names.kappa <- c()
  if (family == "ordinal")  names.kappa <- paste('kappa', 1:(J - 2) , sep = '.')
  
  ## Names for models with random parameters
  names.random <- c()
  if (R.model){
    ## Some warnings
    if (!correlation){
      ndist <- ranp[! (ranp %in% c("cn", "ln", "n", "u", "t", "sb"))]
      if (length(ndist) > 0){
        udstr <- paste("unknown distribution", paste(unique(ndist), collapse = ", "))
        stop(udstr)
      }
    } else {
      ndist <- ranp[! (ranp %in% c("cn", "ln", "n"))]
      if (length(ndist) > 0){
        udstr <- paste("Correlated parameters is suitable for distribution from normal, such as cn, ln or n")
        stop(udstr)
      }
    }
    novar <- names(ranp)[!((names(ranp) %in% colnames(X)))]
    if (length(novar) > 0 ){
      uvar <- paste("The following random variables are not in the data: ", paste(unique(novar), collapse = ", "))
      stop(uvar)
    }
    Vara <- sort(match(names(ranp), colnames(X))) 
    Varc <- (1:ncol(X))[- Vara]
    fixed <- !(length(Varc) == 0)
    Xa   <- X[ , Vara, drop = F]                        
    Xc   <- X[ , Varc, drop = F]
    allX <- if(fixed) cbind(Xc, Xa) else  Xa
    colnamesX <- colnames(X)
    names.f <- colnamesX[Varc] 
    names.b <- paste('mean', colnamesX[Vara], sep = '.')
    if (!correlation) {
      names.sd <- paste('sd', colnamesX[Vara], sep = '.')
    } else {
      names.sd <- c()
      Ka <- length(ranp)
      for (i in 1:Ka){
        names.sd <- c(names.sd,
                      paste('sd', names(ranp)[i], names(ranp)[i:Ka], sep = '.')
        )
      }
    }
    names.phi <- c()
    if (Hier) names.phi <- unlist(lapply(names(mvar), function(x) outer(x, mvar[[x]], FUN = paste, sep = ".")))
    names.random <- c(names.b, names.phi, names.sd)
  }  else {
    names.f <- colnames(X) 
  }
  
  all.names   <- c(names.kappa, names.f, names.random)
  if (is.null(start)){
    if (family == "ordinal"){
      if (R.model)  theta <- coef(lm(unclass(y) ~ allX - 1)) else theta <- coef(lm(unclass(y) ~ X - 1))
      z <- as.integer(table(unclass(y)))
      z <- (cumsum(z) / sum(z))[1:(J-1)]
      start.kappa <- switch(link,
                            "probit" = qnorm(z),
                            "logit"  = qlogis(z))
      theta       <- c(log(diff(start.kappa)), theta)
    } else {
      if (R.model)  theta <- coef(lm(y ~ allX - 1)) else theta <- coef(lm(y ~ X - 1))
    }
    names(theta) <- c(names.kappa, names.f)
  }
  
  ## Initial value if random parameters and start is null
  if (R.model && is.null(start)) {
    callst        <- call("maxLik")
    callst$start  <- theta
    callst$method <- switch(family,
                            "ordinal"  = 'bfgs',
                            "binomial" = 'nr',
                            "poisson"  = 'nr')
    callst$weights <- as.name("weights")
    callst$X <- as.name('allX') 
    callst$y <- as.name('y')
    callst$link  <- link
    if (family == "poisson")  callst$logLik <- as.name('lnpoisson')
    if (family == "binomial") callst$logLik <- as.name('lnbinary')
    if (family == "ordinal")  callst$logLik <- as.name('lnordered')
    start.fixed  <- coef(eval(callst, sys.frame(which = nframe)))
    if (is.null(start.fixed)) stop("attempt to find suitable starting values failed")
    start.random <- c(rep(0, length(names.phi)), rep(init.ran, length(names.sd)))
    theta        <- c(start.fixed, start.random)
    names(theta) <- all.names
  }
  
  ## If initial value start is not null
  if (!is.null(start)){
    theta <- start
    if (length(start) != length(all.names)) stop('Incorrect Number of Initial Parameters')
    names(theta) <- all.names
  }
  
  ## Check if some variable is log-normal or Johnson Sb
  if(R.model){
    if (any(ranp == "ln")){
      ln <- paste('mean', names(ranp[ranp == "ln"]), sep = '.')
      if (sum(theta[ln] < 0) >= 1)  stop("Some variables specified as ln have negative values in the non-random parameter model. Try using the negative of the variable.")
      theta[ln] <- log(theta[ln])
    }
    if (any(ranp == "sb")){
      sb <- paste('mean', names(ranp[ranp == "sb"]), sep = '.')
      theta[sb] <- exp(theta[sb]) / (1 + exp(theta[sb])) 
    }
  }
  if (print.init){
    cat("\nStarting Values:\n")
    print(theta)
  } 
 #######################################################################
 # 5) Estimate the model using maxLik and passing the correct arguments
 #######################################################################
 opt <- callT
 ## Maximization control arguments
 m <- match(c("print.level", "ftol", "tol", "reltol",
              "gradtol", "steptol", "lambdatol", "qrtol",
              "iterlim", "fixed", "activePar", "method"),
               names(opt), 0L)
 opt <- opt[c(1L, m)]
 opt$start <- theta
 ## Optimization code name
 opt[[1]] <- as.name('maxLik')
 
 #Variables
 opt[c('X', 'y')] <- list(as.name('X'), as.name('y'))
 
 ## Weights
 if (is.null(weights)) weights <- 1
 opt$weights <- as.name('weights')
 
 ## Link
 opt$link  <- link
 
 ## loglik for standard Models
 if (family == "poisson")  opt$logLik <- as.name('lnpoisson')
 if (family == "binomial") opt$logLik <- as.name('lnbinary')
 if (family == "ordinal")  opt$logLik <- as.name('lnordered')
 
 ## Arguments for random parameters
 if (R.model) {
   opt[c('R', 'seed', 'ranp', 'correlation','haltons', 'id', 'S', 'mvar')] <-
     list(as.name('R'), as.name('seed'), as.name('ranp'), as.name('correlation'), as.name('haltons'), as.name('id'),
          as.name('S'), as.name('mvar'))
   if (family == "binomial") opt$logLik <- as.name('lnlbinary.ran')
   if (family == "poisson")  opt$logLik <- as.name('lnlpoisson.ran')
   if (family == "ordinal")  opt$logLik <- as.name('lnordered.ran')
 }
 
  ## Optimizing the ML
  x <- eval(opt, sys.frame(which = nframe))
 ###################################
 # 6) Extract predicted probabilities, 
 # conditional means of ranp, etc.
 ##################################
 
 ## Get probability, and conditional beta
 opt$steptol <- opt$logLik <- opt$iterlim <- opt$method <- opt$print.level <- opt$tol <- opt$ftol <- NULL
 names(opt)[[2]] <- 'theta'
 betahat <- coef(x)
 if (!is.null(ranp)){
   opt$make.estb <- TRUE
   opt$gradient <- FALSE
   if (family == "binomial") opt[[1]]  <- as.name('lnlbinary.ran')
   if (family == "poisson")  opt[[1]]  <- as.name('lnlpoisson.ran')
   if (family == "ordinal")  opt[[1]]  <- as.name('lnordered.ran')
   if (!correlation) {
     diag.sd <- paste('sd',  names(ranp),  sep = '.')
     betahat <- ifelse(names(betahat) %in% diag.sd, abs(betahat), betahat)
     names(betahat) <- names(coef(x))
   }
   opt[[2]] <- betahat
   again <- eval(opt, sys.frame(which = nframe))
   probabilities <- attr(again, 'probabilities')
   bi      <- attr(again, 'bi')
   Qir     <- attr(again, 'Qir')
   x$estimate      <- betahat
 } else {
   if (family == "poisson")  opt[[1]] <- as.name('lnpoisson')
   if (family == "binomial") opt[[1]] <- as.name('lnbinary')
   if (family == "ordinal")  opt[[1]] <- as.name('lnordered')
   again <- eval(opt, sys.frame(which = nframe))
   probabilities <- drop(attr(again, 'probabilities'))
   bi <- Qir <- NULL
 }
 
 ## Ordered Model
 if (family == "ordinal"){
   J <- length(levels(y))
   attr(x$estimate, "alphas") <- x$estimate[1:(J - 2)]
   kappas <- cumsum(c(exp(x$estimate[1:(J - 2)])))
   names(kappas) <- names.kappa
   x$estimate[names.kappa]  <- kappas
   attr(x$estimate, "fixed") <- x$estimate[-c(1:(J - 2))]
 }
 resid <- drop(unclass(y) - probabilities) 
###########################
# 7) Put results in form
###########################
 logLik <- structure(list(
                    maximum     = logLik(x),
                    gradient    = x$gradient,
                    nobs        = nrow(X),
                    gradientObs = x$gradientObs,
                    hessian     = hessian(x),
                    iterations  = nIter(x),
                    type        = maximType(x),
                    code        = returnCode(x),
                    nparam      = length(x$estimate),
                    message     = returnMessage(x)),
                    class = "logLik"
                   )
  result <- structure(
                      list(
                        coefficients  = x$estimate,
                        family        = family,
                        link          = link,
                        logLik        = logLik,
                        mf            = mf,
                        formula       = formula,
                        time          = proc.time()-start.time,
                        freq          = freq,
                        draws         = haltons,
                        R.model       = R.model,
                        R             = R,
                        bi            = bi,
                        Qir           = Qir,
                        ranp          = ranp,
                        probabilities = probabilities,
                        residuals     = resid,
                        correlation   = correlation,
                        call          = callT),
                    class = 'Rchoice')
 result
}