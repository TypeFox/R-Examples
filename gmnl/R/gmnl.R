#' Estimate Multinomial Logit Models with Observed and Unobserved Individual Heterogeneity.
#' 
#' Estimate different types of multinomial logit models with observed and unobserved individual heterogneity, such as
#' MIXL, S-MNL, G-MNL, LC and MM-MNL models. These models are estimated using  Maximum Simulated Likelihood. It supports both cross-sectional and panel data.
#' 
#' 
#' @name gmnl
#' @param x,object and object of class \code{gmnl}.
#' @param formula a symbolic description of the model to be estimated. The formula is divided in five parts, each of them separated by the symbol \code{|}. The first part is reserved for alternative-specific variables with a generic coefficient. The second part corresponds to individual-specific variables with an alternative specific coefficients. The third part corresponds to alternative-specific variables with an alternative-specific coefficident. The fourth part is reserved for time-invariant variables that modify the mean of the random parameters. Finally, the fifth part is reserved for time-invariant variables that enter in the scale coefficient or in the probability assignment in models with latent classes. 
#' @param data the data of class \code{\link[mlogit]{mlogit.data}}.
#' @param subset an optional vector specifying a subset of observations.
#' @param weights an optional vector of weights. Default to 1.
#' @param na.action a function wich indicated what should happen when the data
#' contains \code{NA}'s.
#' @param model a string indicating which model is estimated. The options are "\code{mnl}" for the Multinomial Logit Model, "\code{mixl}" for the Mixed Logit Model, "\code{smnl}" for the Scaled Multinomial Logit Model, "\code{gmnl}" for the Generalized Multinomial Logit Model, "\code{lc}" for the Latent Class Multinomial Logit Model, and "\code{mm}" for the Mixed-Mixed Multinomial Logit Model.
#' @param start a vector of starting values.
#' @param ranp a named vector whose names are the random parameters and values the distribution:
#' "\code{n}" for normal, "\code{ln}" for log-normal, "\code{cn}" for truncated normal, "\code{u}" for uniform, "\code{t}" for triangular, "\code{sb}" for Sb Johnson.
#' @param R the number of draws of pseudo-random numbers if \code{ranp} is not \code{NULL}.
#' @param Q number of classes for LC or MM-MNL models.
#' @param haltons only relevant if \code{ranp} is not \code{NULL}. If \code{haltons = NULL}, pseudo-random numbers are used instead of Halton sequences. If \code{haltons=NA}, the first \eqn{K} primes are used to generates the Halton draws, where \eqn{K} is the number of random parameters, and 15 of the initial sequence of elements are dropped. Otherwise, \code{haltons} should be a list with elements \code{prime} and \code{drop}.
#' @param mvar only valid if \code{ranp} is not \code{NULL}. This is a named list, where the names correspond to the variables with random parameters, and the values correspond to the variables that enter in the mean of each random parameters.
#' @param seed seed for the random number generator. Default is \code{seed = 12345}.
#' @param correlation only relevant if \code{ranp} is not \code{NULL}. If true, the correlation across random parameters is taken into account.
#' @param bound.err only relevenat if model is \code{smnl} or \code{gmnl}. It indicates at which values the draws for the scale parameter are truncated. By default \code{bound.err = 2}, therefore a truncated normal distribution with truncation at -2 and +2 is used.
#' @param panel if \code{TRUE} a panel data model is estimated.
#' @param hgamma a string indicated how to estimate the parameter gamma. If "\code{direct}", then \eqn{\gamma} is estimated directly, if "\code{indirect}" then \eqn{\gamma ^*} is estimated, where \eqn{\gamma = \exp(\gamma^*)/(1 + \exp(\gamma^*))}.
#' @param reflevel the base alternative.
#' @param init.tau initial value for the \eqn{\tau} parameter.
#' @param init.gamma initial value  for \eqn{\gamma}.
#' @param notscale only relevant if model is \code{smnl} or \code{gmnl}. It is a vector indicating which variables should not be scaled.
#' @param print.init if \code{TRUE}, the initial values for the optimization procedure are printed.
#' @param gradient if \code{TRUE}, analytical gradients are used for the optimization procedure.
#' @param typeR if \code{TRUE}, truncated normal draws are used for the scale parameter, if \code{FALSE} the procedure suggested by Greene (2010) is used.
#' @param ... additional arguments to be passed to \code{\link[maxLik]{maxLik}}, which depend in the maximization routine.
#' @param new an updated formula for the \code{update} method.
#' @param digits the number of digits.
#' @param width width.
#' @param outcome if \code{TRUE}, then the \code{fitted} and \code{residuals} methods return a vector that corresponds to the chosen alternative, otherwise it returns a matrix where each column corresponds to each alternative.
#' @seealso \code{\link[mlogit]{mlogit}}, \code{\link[mlogit]{mlogit.data}},  \code{\link[maxLik]{maxLik}}, \code{Rchoice}  
#' 
#' 
#' @details Let the utility to person \eqn{i} from choosing alternative \eqn{j} on choice occasion \eqn{t} be: \deqn{U_{ijt} = \beta_{i}x_{ijt} + \epsilon_{ijt}} where \eqn{\epsilon_{ijt}} is i.i.d extreme value, and \eqn{\beta_i} vary across individuals. Each model estimated by \code{gmnl} depends on how \eqn{\beta_i} is specified. The options are the following:
#' \enumerate{
#' \item S-MNL if \eqn{\beta_i=\sigma_i\beta}, where the scale \eqn{\sigma_i} varies across individuals.
#' \item MIXL  if \eqn{\beta_i=\beta + s\eta_i}, where \eqn{\eta_i} is a draw from some distribution. For example, if \eqn{\beta_i\sim N(\beta, s^2)}, then \eqn{\eta_i\sim N(0, 1)}.
#' \item GMNL if \eqn{\beta_i=\sigma_i\beta + \gamma s\eta_i + \sigma_i(1-\gamma)s\eta_i}, where \eqn{\sigma_i} is the scale parameter, and \eqn{\gamma} is a parameter that controls how the variance of residual taste heterogeneity varies with scale.
#' \item LC if \eqn{\beta_i=\beta_q} with probability \eqn{w_{iq}} for \eqn{q = 1,...,Q}, where \eqn{Q} is the total number of classes.
#' \item MM-MIXL if  \eqn{\beta_i\sim f(\beta_q, \Sigma_q)} with probability \eqn{w_{iq}} for \eqn{q = 1,...,Q}, where \eqn{Q} is the total number of classes.
#' }
#' 
#' Observed heterogeneity can be also accommodated in the random parameters when the MIXL is estimated by including individual-specific covariates. Specifically, the vector of random coefficients is \deqn{\beta_i=\beta +\Pi z_i + L\eta_i} where \eqn{z_i} is a set of characteristics of individual \eqn{i} that influence the mean of the taste parameters; and \eqn{\Pi} is matrix of parameters. To estimate this model, the fourth part of the \code{formula} should be specified along with the \code{mvar} argument.
#' 
#' 
#' One can also allow the mean of the scale to differ across individuals by including individual-specific characteristics. Thus, the scale parameters can be written as \deqn{\exp(\bar{\sigma} + \delta h_i + \tau \upsilon_i)} where \eqn{h_i} is a vector of attributes of individual \eqn{i}. To estimate this model, the fifth part of the \code{formula} should include the variables that enter \eqn{h_i}.
#' 
#' For models with latent classes,  the class assignment is modeled as a semi-parametric multinomial logit format \deqn{w_{iq}= \frac{\exp(\gamma_q)}{\sum_{q=1}^Q\exp(\gamma_q)}} for \eqn{q = 1,...,Q, \gamma_1 = 0}. Latent class models (LC and MM-MIXL) requires at least that a constant should be specified in the fifth part of the \code{formula}. If the class assignment, \eqn{w_{iq}}, is also determined by socio-economic characteristics, these variables can be also included in the fifth part.  
#' 
#' 
#' Models that involve random parameters are estimated using Maximum Simulated Likelihood using the \code{maxLik} function of \code{\link[maxLik]{maxLik}} package.
#' 
#' @return An object of class ``\code{gmnl}'' with the following elements
#' \item{coefficients}{the named vector of coefficients,}
#' \item{logLik}{a set of values of the maximum likelihood procedure,}   
#' \item{mf}{the model framed used,} 
#' \item{formula}{the formula (a \code{gFormula} object),}
#' \item{time}{\code{proc.time()} minus the start time,}
#' \item{freq}{frequency of dependent variable,}
#' \item{draws}{type of draws used,}
#' \item{model}{the fitted model,}
#' \item{R}{number of draws used,}
#' \item{ranp}{vector indicating the variables with random parameters and their distribution,}
#' \item{residuals}{the residuals,}
#' \item{correlation}{whether the model is fitted assuming that the random parameters are correlated,}
#' \item{bi}{matrix of conditional expectation of random parameters,}
#' \item{Q}{number of classes,}
#' \item{call}{the matched call.}   
#' 
#' @references
#' \itemize{
#' \item Keane, M., & Wasi, N. (2013). Comparing alternative models of heterogeneity in consumer choice behavior. Journal of Applied Econometrics, 28(6), 1018-1045.
#' \item Fiebig, D. G., Keane, M. P., Louviere, J., & Wasi, N. (2010). The generalized multinomial logit model: accounting for scale and coefficient heterogeneity. Marketing Science, 29(3), 393-421.
#' \item Greene, W. H., & Hensher, D. A. (2010). Does scale heterogeneity across individuals matter? An empirical assessment of alternative logit models. Transportation, 37(3), 413-428.
#' \item Train, K. (2009). Discrete choice methods with simulation. Cambridge University Press.
#' }
#' @author Mauricio Sarrias \email{msarrias86@@gmail.com}
#' @import Formula maxLik truncnorm stats
#' @importFrom mlogit mlogit.data
#' @export
#' @examples
#' ## Examples using the Fishing data set from the AER package
#' data("TravelMode", package = "AER")
#' library(mlogit)
#' TM <- mlogit.data(TravelMode, choice = "choice", shape = "long", 
#'                  alt.levels = c("air", "train", "bus", "car"), chid.var = "individual")
#' \dontrun{
#' ## S-MNL model, ASCs not scaled
#' smnl <- gmnl(choice ~ wait + vcost + travel + gcost| 1, data = TM, 
#'              model = "smnl", R = 100, 
#'              notscale = c(1, 1, 1, rep(0, 4)))
#' summary(smnl)
#' 
#' ## MIXL model with observed heterogeneity
#' mixl.hier <- gmnl(choice ~ vcost + gcost + travel + wait | 1 | 0 | income + size - 1,
#'                  data = TM,
#'                  model = "mixl",
#'                  ranp = c(travel = "t", wait = "n"),
#'                  mvar = list(travel = c("income","size"), wait = c("income")),
#'                  R = 30,
#'                  haltons = list("primes"= c(2, 17), "drop" = rep(19, 2)))
#' summary(mixl.hier)
#' 
#' ## Examples using the Electricity data set from the mlogit package
#' data("Electricity", package = "mlogit")
#' Electr <- mlogit.data(Electricity, id.var = "id", choice = "choice",
#'                      varying = 3:26, shape = "wide", sep = "")
#'                      
#' ## Estimate a MIXL model with correlated random parameters
#' Elec.cor <- gmnl(choice ~ pf + cl + loc + wk + tod + seas| 0, data = Electr,
#'                  subset = 1:3000,
#'                  model = 'mixl',
#'                  R = 10,
#'                  panel = TRUE,
#'                  ranp = c(cl = "n", loc = "n", wk = "n", tod = "n", seas = "n"),
#'                  correlation = TRUE)
#' summary(Elec.cor)
#' cov.gmnl(Elec.cor)
#' se.cov.gmnl(Elec.cor)
#' se.cov.gmnl(Elec.cor, sd = TRUE)
#' cor.gmnl(Elec.cor)
#' 
#' ## Estimate a G-MNL model, where ASCs are also random
#' Electr$asc2 <- as.numeric(Electr$alt == 2)
#' Electr$asc3 <- as.numeric(Electr$alt == 3)
#' Electr$asc4 <- as.numeric(Electr$alt == 4)
#' 
#' Elec.gmnl <- gmnl(choice ~ pf + cl + loc + wk + tod + seas + asc2 + asc3 + asc4 | 0,
#'                  data = Electr,
#'                  subset = 1:3000,
#'                  model = 'gmnl',
#'                  R = 30,
#'                  panel = TRUE,
#'                  notscale = c(rep(0, 6), 1, 1, 1),
#'                  ranp = c(cl = "n", loc = "n", wk = "n", tod = "n", seas = "n",
#'                  asc2 = "n", asc3 = "n", asc4 = "n"))
#' summary(Elec.gmnl)
#' 
#' ## Estimate a LC model with 2 classes
#' Elec.lc <- gmnl(choice ~ pf + cl + loc + wk + tod + seas| 0 | 0 | 0 | 1,
#'                data = Electr,
#'                subset = 1:3000,
#'                model = 'lc',
#'                panel = TRUE,
#'                Q = 2)
#' summary(Elec.lc)
#' 
#' ## Estimate a MM-MIXL model
#' Elec.mm <- gmnl(choice ~ pf + cl + loc + wk + tod + seas| 0 | 0 | 0 | 1,
#'                  data = Electr,
#'                  subset = 1:3000,
#'                  model = 'mm',
#'                  R = 30,
#'                  panel = TRUE,
#'                  ranp = c(pf = "n", cl = "n", loc = "n", wk = "n", tod = "n",
#'                  seas = "n"),
#'                  Q = 2,
#'                  iterlim = 500)
#' summary(Elec.mm)
#' }
gmnl <- function(formula, data, subset, weights, na.action,
                 model = c("mnl", "mixl","smnl", "gmnl", "lc", "mm"),
                 start = NULL, ranp = NULL, R = 40, Q = 2, haltons = NA, mvar = NULL, 
                 seed = 12345, correlation = FALSE, bound.err = 2, panel = FALSE,
                 hgamma = c("direct", "indirect"), reflevel = NULL, 
                 init.tau = 0.1, init.gamma = 0.1, notscale = NULL, print.init = FALSE, 
                 gradient = TRUE, typeR = TRUE, ...){
  ####################
  # 1) Check arguments
  ####################
  start.time <- proc.time()
  callT <- match.call(expand.dots = TRUE)
  callF <- match.call(expand.dots = FALSE)
  
  formula <- callF$formula <- gFormula(formula)
  nframe  <- length(sys.calls())
  
  ## Check models 
  model <- match.arg(model)
  has.rand <- !is.null(ranp)
  if (model == "mixl" && !has.rand) stop("mixl model needs ranp to be specified")
  if (model == "gmnl" && !has.rand) stop("gmnl model needs ranp to be specified")
  if (model == "mm" && !has.rand)  stop("mn model needs ranp to be specified")
  if ((model == "lc" || model == "mn") && Q < 2) stop("Classes cannot be lower than 2")
  
  if (model == "mnl") {
    if (is.null(callT$method)) callT$method <- 'nr'
  } else {
    if (is.null(callT$method)) callT$method <- 'bfgs'
  }
  
  ########################
  # 2) Check Model Frame
  ########################
  if (!inherits(data, "mlogit.data")) stop("Data must be of class mlogit.data")
  mf <- callT
  m  <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## Change the reference level
  if (!is.null(reflevel)) attr(mf, "index")[["alt"]] <- relevel(attr(mf, "index")[["alt"]], reflevel)
  
  ## Indeces
  index <- attr(mf, "index")
  alt   <- index[["alt"]]
  chid  <- index[["chid"]]
  alt.lev <- levels(alt)
  J <- length(alt.lev)
  
  ## Check Panel
  if (panel) {
    if (model == "mnl") stop("Panel is not relevant for mnl model")
    id <- index[["id"]]
    if (is.null(id)) stop("No individual index")
    id <- split(index[["id"]], alt)[[1]]
  } else{
    #id <- NULL
    # This does not confuse the code when the data has id.var but the model is panel = FALSE
    id <- attr(mf, "index")[["id"]] <- NULL
  }
  
  ##########################################
  # 3) Extract the elements of the model
  ##########################################
  y <- model.response(mf)
  X <- model.matrix(formula, mf)
  
  ## Check variables for the mean and/or heterogeneity
  has.mvar <- has.othervar(formula, 4)
  if (has.mvar){
    if (model == "mnl" || model == "smnl") stop(paste("Variables for mean are not relevant for", paste(model, collapse = ": ")))
    if (is.null(mvar))  stop("mvar is null")
    if (!is.list(mvar)) stop("mvar is not a list")
    rvar <- names(mvar)[!(names(mvar) %in% names(ranp))]
    if (length(rvar) > 0) {
      udstr <- paste("The following variables are not specified in the argument ranp:", paste(unique(rvar), collapse = ", "))
      stop(udstr)
    }
    Z <- model.matrix(formula, mf, rhs = 4) 
    for (i in 1:length(mvar)) {
      rvar <- mvar[[i]][!(mvar[[i]] %in% colnames(Z))]
      if (length(rvar) > 0) {
        udstr <- paste("The following variables are not specified in the formula: ", paste(unique(rvar), collapse = ", "))
        stop(udstr)
      }
    }
  } else Z <- NULL
  
  has.het <- has.othervar(formula, 5)
  if ((model == "lc" || model == "mm") && !has.het) stop("lc needs variables for class probabilities")
  if (has.het){
    if (model == "mnl" || model == "mixl") stop(paste("variables for scale are not relevant for", paste(model, collapse = ": ")))
    if (model == "lc" || model == "mm") H <- model.matrix(formula, mf, rhs = 5, Q = Q) else H <- model.matrix(formula, mf, rhs = 5) 
  } else H <- NULL
  
  # Weights
  if (any(names(mf) == "(weights)")) {
    weights <- mf[["(weights)"]] <- mf[["(weights)"]]/mean(mf[["(weights)"]])
    weights <- split(weights, alt)[[1]]
  } else weights <- NULL
  
  freq <- table(alt[y])
  
  # List of Variables
  Xl <- vector(length = J, mode = "list")
  names(Xl) <- levels(alt)
  for (i in levels(alt)) {
    Xl[[i]] <- X[alt == i, , drop = FALSE]
  }
  yl <- split(y, alt)
  yl <- lapply(yl, function(x){x[is.na(x)] <- FALSE ; x})
  
  
  ########################################################
  #   4) Initial Values
  ########################################################
  
  # Names of means
  mean.names <- colnames(X) 
  
  ## Standard deviation of random parameters
  names.stds <- start.stds <- c()
  if(has.rand) {
    if (!correlation) {
      ndist <- ranp[!(ranp %in% c("cn", "ln", "n", "u", "t", "sb"))]
      if (length(ndist) > 0) {
        udstr <- paste("unknown distribution", paste(unique(ndist), collapse = ", "))
        stop(udstr)
      }
    } else {
      ndist <- ranp[!(ranp %in% c("cn", "ln", "n"))]
      if (length(ndist) > 0) {
        udstr <- paste("Correlated parameters is suitable for distribution from normal, such as cn, ln or n")
        stop(udstr)
      }
    }
    if (model == "gmnl") {
      ndist <- ranp[!(ranp %in% c("n", "u", "t", "ln"))]
      if (length(ndist) > 0) {
        udstr <- paste("Coefficients from G-MNL model can only be distributed as n, u, ln or t")
        stop(udstr)
      }
    }
    namesX <- colnames(X)
    novar <- names(ranp)[!((names(ranp) %in% namesX))]
    if (length(novar) > 0 ) {
      uvar <- paste("The following random variables are not in the data: ", paste(unique(novar), collapse = ", "))
      stop(uvar)
    }
    Vara <- sort(match(names(ranp), namesX)) 
    Varc <- (1:ncol(X))[-Vara]
    fixed <- !(length(Varc) == 0)
    Xa   <- lapply(Xl, function(x) x[, Vara, drop = F])                        
    Xc   <- lapply(Xl, function(x) x[, Varc, drop = F]) 
    allX <- if (fixed) mapply(cbind, Xc, Xa, SIMPLIFY = FALSE) else  Xa
    mean.names <- colnames(allX[[1]])
    nrap  <- length(Vara)
    if (!correlation) {
      names.stds <- paste("sd", namesX[Vara], sep = ".")
      start.stds <- rep(.1, nrap)
      names(start.stds) <- names.stds
    } else {
      names.stds <- c()
      Ka <- length(ranp)
      for (i in 1:Ka) {
        names.stds <- c(names.stds,
                        paste('sd', names(ranp)[i], names(ranp)[i:Ka], sep = '.'))
        start.stds <- rep(.1, .5 * nrap * (nrap + 1))
        names(start.stds) <- names.stds
      }
    }
  } else {
    mean.names <- colnames(X)
    allX <- Xl
  }
  
  ## Heterogeneity in scale or Class variables
  names.het <- start.het <- c()
  if (has.het) {
    if (model == "lc" | model == "mm") {
      start.het <- rep(0, ncol(H))
      names.het <- colnames(H)
      classes <- attr(H, "alt")
      Hl <- vector(length = Q, mode = "list")
      names(Hl) <- levels(classes)
      for (i in levels(classes)) {
        Hl[[i]] <- H[classes == i, , drop = FALSE]
      }
    } else {
      start.het <- rep(0, ncol(H))
      names.het <- paste("het", colnames(H), sep = ".")
      names(start.het) <- names.het 
    }
  }
  
  ## Heterogeneity in mean
  names.phi <- start.phi <- c()
  if (has.mvar) {
    names.phi <- unlist(lapply(names(mvar), function(x) outer(x, mvar[[x]], FUN = paste, sep = ".")))
    start.phi <- rep(0, length(names.phi))
    names(start.phi) <- names.phi
  }
  
  if (model == "lc") {
    lc.names <- c()
    for (i in 1:Q) {
      lc.names <- c(lc.names, paste('class', i, colnames(X), sep = '.'))
    }  
  }
  
  if (model == "mm") {
    lc.names <- c()
    ls.names <- c()
    for (i in 1:Q) {
      lc.names <- c(lc.names, paste('class', i, colnames(X), sep = '.'))
      ls.names <- c(ls.names, paste('class', i, names.stds, sep = '.'))
    }  
  }
  
  if (model == "smnl") allnames <- c(mean.names, "tau", names.het)
  if (model == "mixl") allnames <- c(mean.names, names.phi, names.stds)
  if (model == "gmnl") allnames <- c(mean.names, names.het, names.stds, "tau", "gamma")
  if (model == "lc")   allnames <- c(lc.names, names.het)
  if (model == "mm")   allnames <- c(lc.names, ls.names,  names.het)
  
  
  # If null start, estimate mean parameters from mnl
  if (!is.null(start)) {
    if (length(start) != length(allnames)) stop("Incorrect number of initial parameters")
    theta <- start
    names(theta) <- allnames
  } else {
    K <- ncol(allX[[1]])
    theta <- rep(0, K)
    names(theta) <- mean.names
    if (model != "mnl") {
      # Estimate mean parameters from mnl
      calls <- call("maxLik")
      calls$start  <- theta
      calls$method   <- 'nr'
      calls$weights <- as.name("weights")
      calls$X <- as.name('allX') 
      calls$y <- as.name('yl')
      calls$logLik <- as.name('ll.mlogit')
      mean <- coef(eval(calls, sys.frame(which = nframe)))
      if (model == "lc" || model == "mm") {
        lc.mean <- c()
        init.shift <- seq(-0.02, 0.02, length.out = Q)
        for (i in 1:Q) {
          lc.mean <- c(lc.mean, mean + init.shift[i])
        }
      }
      if (model  == "mm") {
        ls.mean <- c()
        init.shift <- seq(-0.02, 0.02, length.out = Q)
        for(i in 1:Q) {
          ls.mean <- c(ls.mean, start.stds + init.shift[i])
        }
      }
      theta <- switch(model,
                      "smnl" = c(mean, init.tau, start.het),
                      "mixl" = c(mean, start.phi, start.stds),
                      "gmnl" = c(mean, start.phi, start.het, start.stds, init.tau, init.gamma),
                      "lc"   = c(lc.mean, start.het),
                      "mm"   = c(lc.mean, ls.mean, start.het))
      names(theta) <- allnames
    }
  }
  
  if (model == "smnl" || model == "gmnl") {
    if (!is.null(notscale)) {
      if (length(notscale) != ncol(X)) stop("Not scaled variables vector is not the same length as initial parameters")
      names(notscale) <- mean.names
      wns <- names(notscale[notscale == 1])
      cat("\nThe following variables are not scaled:\n")
      print(wns)
    } else{
      notscale <- rep(0, ncol(X))
      names(notscale) <- mean.names
    }
  }
  if ((model == "mixl" || model == "gmnl" || model == "mm") & is.null(start)) {
    ln <- names(ranp[ranp == "ln"])
    sb <- names(ranp[ranp == "sb"])
    if (length(ln) != 0) {
      if (model == "mixl" || model == "gmnl") {
        if (sum(theta[ln] < 0) >= 1)  stop("Some variables specified as ln have negative values in the clogit")
        theta[ln] <- log(theta[ln])
      } else {
        log.names <- c()
        for (i in 1:Q) log.names <- c(log.names, paste('class', i, ln, sep = '.'))
        if (sum(theta[log.names] < 0) >= 1)  stop("Some variables specified as ln have negative values in the clogit")
        theta[log.names] <- log(theta[log.names])
      }
    }
    if (length(sb) != 0) {
      if (model == "mixl" || model == "gmnl") {
        if (sum(theta[sb] < 0) >= 1)  stop("Some variables specified as sb have negative values in the clogit")
        theta[sb] <- log(theta[sb])
      } else {
        sb.names <- c()
        for (i in 1:Q) sb.names <- c(sb.names, paste('class', i, sb, sep = '.'))
        if (sum(theta[sb.names] < 0) >= 1)  stop("Some variables specified as sb have negative values in the clogit")
        theta[sb.names] <- log(theta[sb.names])
      }
    }
  }  
  if (print.init) {
    cat("\nStarting Values:\n")
    print(theta)
  } 
  
  #######################################################################
  # 5) Estimate the model using maxLik and passing the correct arguments
  #######################################################################
  opt <- callT
  opt$start <- theta
  
  ## Maximization control arguments
  m <- match(c('method', 'print.level', 'iterlim',
               'start','tol', 'ftol', 'steptol', 'fixed', 'constraints'),
             names(opt), 0L)
  opt <- opt[c(1L, m)]
  
  ## Optimization code name
  opt[[1]] <- as.name('maxLik')
  
  #Variables
  opt[c('X', 'y')] <- list(as.name('Xl'), as.name('yl'))
  
  ## Loglikelihood function 
  if (model == "mnl") opt$logLik <- as.name('ll.mlogit') 
  
  ##Gradient
  opt$gradient <- as.name('gradient') 
  
  if (model == "smnl") {
    cat("Estimating SMNL model", "\n")
    opt$logLik <- as.name('ll.smlogit')
    opt[c('R', 'seed', 'bound.err', 'id', 'H', 'notscale', 'typeR')] <- list(as.name('R'), as.name('seed'), as.name('bound.err'), as.name('id'), as.name('H'), as.name('notscale'), as.name('typeR'))
  }
  if (model == "mixl") {
    cat("Estimating MIXL model", "\n")
    opt$logLik <- as.name('ll.mixlog')
    opt[c('R', 'seed', 'ranp', 'id', 'Z','correlation', 'haltons', 'mvar')] <- 
      list(as.name('R'), as.name('seed'), as.name('ranp'), as.name('id'), as.name('Z'), as.name('correlation'), as.name('haltons'),
           as.name('mvar'))
  }
  if (model == "gmnl") {
    cat("Estimating GMNL model", "\n")
    opt$logLik <- as.name('ll.gmlogit')
    hgamma <- match.arg(hgamma)
    opt[c('R', 'seed', 'ranp', 'id', 'H','correlation', 'haltons', 'bound.err', 'hgamma','notscale', 'typeR')] <- 
      list(as.name('R'), as.name('seed'), as.name('ranp'), as.name('id'), as.name('H'), as.name('correlation'), as.name('haltons')
           , as.name('bound.err'), as.name('hgamma'), as.name('notscale'), as.name('typeR'))
  }
  
  ## Weights
  if (is.null(weights)) weights <- 1
  opt$weights <- as.name('weights')
  
  if (model == "lc") {
    cat("Estimating LC model", "\n")
    opt$logLik <- as.name('ll.mlogitlc')
    opt[c('H', 'Q', 'id')] <- list(as.name('Hl'), as.name('Q'), as.name('id'))
  }
  
  if (model == "mm") { 
    cat("Estimating MM-MNL model", "\n")
    opt$logLik <- as.name('ll.mnlogit')
    opt[c('R', 'seed', 'ranp', 'id', 'H', 'correlation', 'haltons', 'Q')] <- 
      list(as.name('R'), as.name('seed'), as.name('ranp'), as.name('id'), as.name('Hl'), as.name('correlation'), as.name('haltons'),
           as.name('Q'))
  }
  x <- eval(opt, sys.frame(which = nframe))
  actpar <- activePar(x)
  ###################################
  # 6) Extract predicted probabilities, 
  # conditional means of ranp, etc.
  ##################################
  
  names(opt)[[2]] <- 'theta'
  betahat <- coef(x)
  opt$gradient <- FALSE
  opt$get.bi <- TRUE
  opt$fixed <- opt$steptol <- opt$iterlim <- opt$method <- opt$print.level <- opt$tol <- opt$ftol <- opt$logLik <- opt$start <- NULL
  opt$constraints <- NULL
  if (model != "mnl") {
    opt[[1]] <- switch(model,
                       "smnl" = as.name('ll.smlogit'),
                       "mixl" = as.name('ll.mixlog'),
                       "gmnl" = as.name('ll.gmlogit'),
                       "lc"   = as.name('ll.mlogitlc'),
                       "mm"   = as.name('ll.mnlogit'))
    if (has.rand) {
      if (!correlation) {
        if (model != "mm") {
          diag.sd <- paste("sd", namesX[Vara], sep = ".")
          betahat <- ifelse(names(betahat) %in% diag.sd, abs(betahat), betahat)
          names(betahat) <- names(coef(x))
        } else {
          betahat <- ifelse(names(betahat) %in% ls.names, abs(betahat), betahat)
          names(betahat) <- names(coef(x))
        }
      }
    }
    opt[[2]] <- betahat
    again <- eval(opt, sys.frame(which = nframe))
    bi  <- attr(again, 'bi')
    Qir <- attr(again, 'Qir')
    x$estimate      <- betahat
  } else {
    opt$hessian <- FALSE
    opt[[1]] <- as.name('ll.mlogit') 
    names(opt)[[2]] <- 'theta'
    opt[[2]] <- betahat
    again <- eval(opt, sys.frame(which = nframe))
    bi <- NULL
    Qir <- NULL
  }
  prob.alt <- attr(again, 'prob.alt')
  prob.ind <- attr(again, 'prob.ind')
  residuals <- Reduce("cbind", yl) - prob.alt
  colnames(prob.alt) <- colnames(residuals) <- levels(alt)
  
  ###########################
  # 7) Put results in form
  ###########################
  
  if (gradient) gradientObs <- x$gradientObs[, actpar] else gradientObs <- NULL
  
  logLik <- structure(list(
                          maximum     = logLik(x),
                          gradient    = x$gradient[actpar],
                          nobs        = nrow(X)/J,
                          gradientObs = gradientObs,
                          hessian     = hessian(x)[actpar, actpar],
                          iterations  = nIter(x),
                          type        = maximType(x),
                          code        = returnCode(x),
                          nparam      = length(x$estimate),
                          message     = returnMessage(x)),
                          class = "logLik")
     
     
  result <- structure(
                      list(
                        coefficients  = x$estimate[actpar],
                        logLik        = logLik,
                        mf            = mf,
                        formula       = formula,
                        time          = proc.time() - start.time,
                        freq          = freq,
                        draws         = haltons,
                        model         = model,
                        R             = R,
                        ranp          = ranp,
                        residuals     = residuals,
                        correlation   = correlation,
                        prob.alt      = prob.alt,
                        prob.ind      = prob.ind,
                        bi            = bi,
                        Qir           = Qir,
                        notscale      = notscale,
                        Q             = Q, 
                        call          = callT),
                      class = 'gmnl'
                    )
  result
}