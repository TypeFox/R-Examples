#' Risk model for binary outcome
#'
#' @param model Formula specifying the risk model
#' @param D number of samples for the simulated annealing integration
#' @param risk Function for transforming a linear predictor into a probability.
#'   E.g., risk.logit for the logistic model, risk.probit for the probit model
#' @export


risk_binary <- function(model = Y ~ S.1 * Z, D = 5000, risk = risk.logit){

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    expanded <- expand_augdata(model, psdesign, D = D)

    trtmat <- expanded$noimp
    Y.trt <- expanded$noimp.Y

    untrt.expand <- expanded$imp
    Y.untrt <- expanded$imp.Y

    likelihood <- function(beta){

      trted <- risk(trtmat %*% beta)^Y.trt *
        (1 - risk(trtmat %*% beta))^(1 - Y.trt)

      # for each W in untrted, sample a bunch from cdf_sbarw and take the mean

      if(!is.null(untrt.expand) & !is.null(Y.untrt)){

        untrted <- matrix(risk(untrt.expand %*% beta)^Y.untrt *
          (1 - risk(untrt.expand %*% beta))^(1 - Y.untrt), nrow = D, byrow = TRUE)

      } else untrted <- matrix(1)

      -1 * (sum(log(trted)) + sum(log(colMeans(untrted))))

    }

    psdesign$risk.function <- function(data, beta){  ## P(D = 1 | S, Z)

      risk(as.vector(model.matrix(model[-2], data) %*% beta))

    }

    psdesign$likelihood <- likelihood
    psdesign$risk.model <- list(model = "binary", args = arglist)
    psdesign$nparam <- ncol(trtmat)
    psdesign$param.names <- colnames(trtmat)

    psdesign

  }
  ## return a likelihood closure

  class(rval) <- c("ps", "riskmodel")
  rval
}


#' Weibull risk model for time to event outcome
#'
#' @param model Formula specifying the risk model. The outcome should be a \link{Surv} object specifying right censoring
#' @param D number of samples for simulated annealing
#'
#' @export

risk_weibull <- function(model = Y ~ S.1 * Z, D = 5000 ){

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    expanded <- expand_augdata(model, psdesign, D = D)

    if(!inherits(expanded$noimp.Y, "Surv")) stop("Requires a survival outcome specified with Surv(time, event)")

    trtmat <- expanded$noimp
    Y.trt <- expanded$noimp.Y[, 1]
    delt.trt <- expanded$noimp.Y[, 2]

    untrt.expand <- expanded$imp
    Y.untrt <- expanded$imp.Y[, 1]
    delt.untrt <- expanded$imp.Y[, 2]

    likelihood <- function(beta){


      gamma0<-beta[1]
      beta0<-beta[-1]
      shapepram<-exp(gamma0)
      scalepram<-exp(trtmat %*% beta0)

      trtlike <- (log(shapepram) - log(scalepram)  +  (log(Y.trt) - log(scalepram)) * (shapepram-1)) * delt.trt - (Y.trt/scalepram)^shapepram


      if(!is.null(untrt.expand) & !is.null(Y.untrt)){

        scale.untrt <- exp(untrt.expand %*% beta0)
        untrted <- matrix((log(shapepram) - log(scale.untrt)  +
                             (log(Y.untrt) - log(scale.untrt)) * (shapepram-1)) * delt.untrt -
                            (Y.untrt/scale.untrt)^shapepram, nrow = D, byrow = TRUE)
      } else untrted <- matrix(1)

      -1 * (sum(trtlike) + sum(colMeans(untrted)))

    }

    psdesign$risk.function <- function(data, beta, t){  #P(Y < t | S, Z)

      mat <- model.matrix(model[-2], data)
      shape <- exp(beta[1])
      beta0 <- beta[-1]

      scale <- as.vector(exp(mat %*% beta0))

      1 - exp(-(t / scale)^shape)

    }

    psdesign$likelihood <- likelihood
    psdesign$risk.model <- list(model = "weibull", args = arglist )
    psdesign$nparam <- ncol(trtmat) + 1
    psdesign$param.names <- c("shape", colnames(trtmat))

    psdesign

  }
  ## return a likelihood closure

  class(rval) <- c("ps", "riskmodel")
  rval

}



#' Exponential risk model for time to event outcome
#'
#' @param model Formula specifying the risk model. The outcome should be a \link{Surv} object specifying right censoring
#' @param D number of samples for simulated annealing
#'
#' @export

risk_exponential <- function(model = Y ~ S.1 * Z, D = 5000 ){

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    expanded <- expand_augdata(model, psdesign, D = D)

    if(!inherits(expanded$noimp.Y, "Surv")) stop("Requires a survival outcome specified with Surv(time, event)")

    trtmat <- expanded$noimp
    Y.trt <- expanded$noimp.Y[, 1]
    delt.trt <- expanded$noimp.Y[, 2]

    untrt.expand <- expanded$imp
    Y.untrt <- expanded$imp.Y[, 1]
    delt.untrt <- expanded$imp.Y[, 2]

    likelihood <- function(beta){

      scalepram <- 1/exp(trtmat %*% beta)

      trtlike <- log(scalepram) * delt.trt - scalepram * Y.trt

      if(!is.null(untrt.expand) & !is.null(Y.untrt)){

        scale.untrt <- 1/exp(untrt.expand %*% beta)
        untrted <- matrix(log(scale.untrt) * delt.untrt - scale.untrt * Y.untrt, nrow = D, byrow = TRUE)
      } else untrted <- matrix(1)

      -1 * (sum(trtlike) + sum(colMeans(untrted)))

    }

    psdesign$risk.function <- function(data, beta, t){ # P(T < t | S, Z)

      scale <- as.vector(exp(model.matrix(model[-2], data) %*% beta))

      1 - exp(-scale * t)

    }

    psdesign$likelihood <- likelihood
    psdesign$risk.model <- list(model = "exponential", args = arglist )
    psdesign$nparam <- ncol(trtmat)
    psdesign$param.names <- colnames(trtmat)

    psdesign

  }
  ## return a likelihood closure

  class(rval) <- c("ps", "riskmodel")
  rval

}


#' Poisson risk model for count outcomes
#'
#' @param model Formula specifying the risk model. The outcome should be an integer of counts. This right side of the formula may contain an \link{offset} term.
#' @param D number of samples for simulated annealing
#'
#' @export

risk_poisson <- function(model = Y ~ S.1 * Z, D = 5000 ){

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    expanded <- expand_augdata(model, psdesign, D = D)

    trtmat <- expanded$noimp
    Y.trt <- expanded$noimp.Y

    untrt.expand <- expanded$imp
    Y.untrt <- expanded$imp.Y

    attr(trtmat, "hasoffset") <- attr(expanded, "hasoffset")

    likelihood <- function(beta){


      if(attr(trtmat, "hasoffset")){
        beta.o <- c(beta, 1)
      } else{
        beta.o <- beta
      }

      lambda <- 1/exp(trtmat %*% beta.o)

      trtlike <- dpois(Y.trt, lambda, log = TRUE)

      if(!is.null(untrt.expand) & !is.null(Y.untrt)){

        lambda.untrt <- 1/exp(untrt.expand %*% beta.o)
        untrted <- matrix(dpois(Y.untrt, lambda.untrt, log = TRUE), nrow = D, byrow = TRUE)
      } else untrted <- matrix(1)

      -1 * (sum(trtlike) + sum(colMeans(untrted)))

    }

    psdesign$risk.function <- function(data, beta, t = 0){ # P(Y <= t | S, Z)

      lambda <- as.vector(exp(model.matrix(model[-2], data) %*% beta))

      ppois(t, 1/lambda, lower.tail = FALSE)

    }

    psdesign$likelihood <- likelihood
    psdesign$risk.model <- list(model = "poisson", args = arglist )
    psdesign$nparam <- ncol(trtmat) + ifelse(attr(trtmat, "hasoffset"), -1, 0)
    psdesign$param.names <- colnames(model.matrix(model, psdesign$augdata))

    psdesign

  }
  ## return a likelihood closure

  class(rval) <- c("ps", "riskmodel")
  rval

}


#' Logit link function
#'
#' @param x A vector of linear predictors
#' @export
#'
#' @return A vector of probabilities

risk.logit <- function(x) {

  exp(x)/(1 + exp(x))

}

#' Probit link function
#'
#' @param x A vector of linear predictors
#' @export
#' @return A vector of probabilities
#'
risk.probit <- function(x) {

  pnorm(x)

}

#' Expand augmented data using the integration function
#'
#'
#' @param model Formula defining the risk model
#' @param psdesign An object of class \link{psdesign}, that contains at least 1 integration model
#' @param D The number of samples to take for the simulated annealing
#'
#' @keywords Internal

expand_augdata <- function(model, psdesign, D = 500){


  vars0 <- sapply(attr(terms(model), "variables"), deparse)[-c(1,2)]
  vars <- unlist(lapply(colnames(psdesign$augdata), function(x){
    if(length(grep(x, vars0, fixed = TRUE) > 0)){
      return(x)
    } else return(NULL)
  }))

  ## check for offset

  mf.c <- model.frame(model, psdesign$augdata)
  hasoffset <- ("offset" %in% names(attributes(terms(mf.c))))

  noimpdex <- rowSums(is.na(psdesign$augdata[, vars, drop = FALSE])) == 0

  trtmat <- model.matrix(model, psdesign$augdata[noimpdex, ])
  Y.trt <- psdesign$augdata[noimpdex, ]$Y

  if(hasoffset){
    vars0 <- sapply(attr(terms(model), "variables"), deparse)[-c(1,2)]
    vars <- unlist(lapply(colnames(psdesign$augdata), function(x){
      if(length(grep(x, vars0, fixed = TRUE) > 0)){
        return(x)
      } else return(NULL)
    }))

    noimpdex <- rowSums(is.na(psdesign$augdata[, vars, drop = FALSE])) == 0

    trtmat <- cbind(trtmat, model.offset(model.frame(model, psdesign$augdata[noimpdex, ])))


  }

  if(all(noimpdex)){

    rt <- list(noimp = trtmat, noimp.Y = Y.trt, imp = NULL, imp.Y = NULL)
    attr(rt, "hasoffset") <- hasoffset
    return(rt)

  } else {

    misvars <- vars[apply(is.na(psdesign$augdata[, vars, drop = FALSE]), MARGIN = 2, any)]
    impdex <- rowSums(is.na(psdesign$augdata[, misvars, drop = FALSE])) > 0
    dex <- (1:nrow(psdesign$augdata))[impdex]
    untrtobs <- psdesign$augdata[rep(dex, D), ]

    for(j in misvars){

      if(!j %in% names(psdesign$integration.models)) stop(paste("Missing values in", j, "but no integration model present."))

      untrtsamp <- c(psdesign$integration.models[[j]]$icdf_sbarw(runif(D)))
      untrtobs[, j] <- untrtsamp

    }

    untrt.expand <- model.matrix(model, untrtobs)
    Y.untrt <- untrtobs$Y

    if(hasoffset){
      untrt.expand <- cbind(untrt.expand, model.offset(model.frame(model, untrtobs)))

    }

    rt <- list(noimp = trtmat, noimp.Y = Y.trt, imp = untrt.expand, imp.Y = Y.untrt)
    attr(rt, "hasoffset") <- hasoffset
    return(rt)

  }

}

