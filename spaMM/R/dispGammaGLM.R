calc_CARdispGammaGLM <-function(data,family=spaMM_Gamma("inverse"),lev=NULL,
                                lambda.Fix, control) {
  if ( ! is.null(lev)) {
    resp <- data$resp/(1-lev) 
    resp[resp==0L] <- 1e-100
    resp[resp>1e150] <- 1e150 ## v1.2 fix for glm -> glm.fit -> .Call(C_Cdqrls, x[good, , drop = FALSE]*w problem
    data$resp <- resp
    weights <- (1-lev)/2 
    loclist <- list(data=data, family=family, weights=weights)
  } else loclist <- list(data=data, family=family)
  if (is.null(lambda.Fix)) {
    loclist$formula <- resp~1+adjd
    loclist$start <- c(1/mean(data$resp),0)
  } else {
    loclist$formula <- resp ~ adjd -1
    loclist$offset <- rep(1/lambda.Fix,nrow(data))
    loclist$start <- c(0)
  }
  loclist$control <- control
  suppressWarnings(do.call("spaMM_glm",loclist))
}


# a specialized version of glm, which further constructs the response value internally; the formula has no lhs 
calc_dispGammaGLM <- 
  function (formula, 
            dev.res, ## cf comments on data argument
            lev,
            data, ## provide only the predictor variables ! alternatively could proceed as calc_CARdispGammaGLM
            family = spaMM_Gamma(link=log), 
            offset=rep.int(0, NROW(dev.res)),
            na.action, start = NULL, etastart, mustart, 
            control, 
            try = FALSE, 
            ## more args from glm() def:
            subset, ##  not used but rethink.
            model = TRUE, ## whether to return the mf, cf end of code
            x = FALSE, ## whether to return the design matrix, cf end of code
            y = TRUE, ## whether to include y... idem 
            contrasts = NULL,
            method="glm.fit",
            ...) {
    Y <- as.numeric(dev.res/((1-lev))) 
  ## potential problem of any corrections is that v_h and lambda estimates may be inconsistent 
  Y[Y==0L] <- 1e-100 ## v_h ~0 occurs eg Poisson, y ~ 1+ ranef(lambda->very small)
  Y[Y>1e150] <- 1e150 ## v1.2 fix for glm -> glm.fit -> .Call(C_Cdqrls, x[good, , drop = FALSE]*w problem
  #
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset","X"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame) 
  if (missing(data)) data <- environment(formula)
  mf <- eval(mf, parent.frame()) ## formula argument is required to eval mf => for mt too
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf, contrasts) ## from glm() code
  #
  offset <- as.vector(model.offset(mf))
  intercept <- FALSE ## previously tried attr(mt, "intercept") > 0L but this is wrong
  #
  weights <- (1-lev)/2 # glm() code : weights <- as.vector(model.weights(mf))
  tryfit <- tryCatch.W.E(eval(call(method, 
                                   x = X, y = Y, weights = weights, offset = offset, family = family, 
                                   control = control, intercept = intercept)))
  fit <- tryfit$value
  warnmess <- fit$warning$message
  if (inherits(fit,"error")) {
    if (try) { ## was for the final glm_lambda, but deprecated.
      ## do nothing: allows the glm to have failed in this case 
    } else {
      ## we need a valid GLM: we fit again with a more controlled method
      fit <- eval(call("spaMM_glm.fit", 
                  x = X, y = Y, weights = weights, offset = offset, family = family, 
                  control = control, intercept = intercept))
      warnmess <- NULL
    }
  }
  if (model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x) fit$x <- X
  if (!y) fit$y <- NULL
  fit <- c(fit, list(call = match.call(), formula = formula, terms = mt, 
                     data = data, offset = offset, control = control, method= "glm.fit",
                     contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf),
                     warnmess=warnmess
                     ))
  class(fit) <- c(fit$class, c("glm", "lm"))
  fit
}