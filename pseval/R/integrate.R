#' Parametric integration model for the missing S(1)
#'
#' @param formula Formula specifying the integration model for the surrogate
#'   under treatment. Generally the candidate surrogate will be on the left side
#'   in the formula, and the BIP or BIPs will be on the right side
#' @param family Assumed distribution for the integration model. Must be
#'   compatible with the \code{family} argument of \link{glm}. Currenly only
#'   Gaussian models are supported
#' @param ... Arguments passed to \link{glm}
#'
#' @export

integrate_parametric <- function(formula, family = gaussian, ...){

  stopifnot(identical(family, gaussian))

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    if(!"integration.models" %in% names(psdesign)) psdesign$integration.models <- NULL
    outname <- paste(formula[[2]])

    missdex <- !is.na(get(paste(formula[[2]]), psdesign$augdata))
    glmdat <- psdesign$augdata[missdex, ]

    if(is.factor(glmdat[, outname])){

      stop("Parametric integration models not valid for factors. Use integrate_nonparametric instead.")

    }

    cdfweights <- NULL # hack to get rid of note
    fit <- glm(formula, data = glmdat, weights = cdfweights, family = gaussian())

    psdesign$integration.models[[outname]]$model <- list(model = "parametric", args = arglist)

    mindelta <- subset(psdesign$augdata, !missdex)

    psdesign$integration.models[[outname]]$cdf_sbarw <-
      function(S.1){

        mu <- predict(fit, newdata = mindelta, type = "response")
        sd <- sd(fit$residuals)

        sapply(S.1, function(s) pnorm(s, mean = mu, sd = sd))

      }
    psdesign$integration.models[[outname]]$icdf_sbarw <-
      function(U.1){

        mu <- predict(fit, newdata = mindelta, type = "response")
        sd <- sd(fit$residuals)

        sapply(U.1, function(u) qnorm(u, mean = mu, sd = sd))

      }

    psdesign
  }

  class(rval) <- c("ps", "integration")
  rval

}

#' Nonparametric integration model for the missing S(1)
#'
#' Both S(1) and the BIP or set of BIPs must be categorical. This model integrates over the estimated distribution of S(1) | BIP
#'
#' @param formula Formula specifying the integration model for the surrogate
#'   under treatment. Generally the candidate surrogate will be on the left side
#'   in the formula, and the BIP or BIPs will be on the right side. In this case
#'   the BIP and the S(1) must be categorical.
#' @param ... Not currently used
#'
#' @export

integrate_nonparametric <- function(formula, ...){

  arglist <- as.list(match.call())
  rval <- function(psdesign){

    if(!"integration.models" %in% names(psdesign)) psdesign$integration.models <- NULL
    outname <- paste(formula[[2]])

    missdex <- !is.na(get(paste(formula[[2]]), psdesign$augdata))

    groups <- model.frame(formula, data = psdesign$augdata)
    mindelta <- subset(psdesign$augdata, !missdex)

    tab.nomiss <- table(model.frame(formula[-2], data = subset(psdesign$augdata, missdex)))
    tab.miss <- table(model.frame(formula[-2], data = mindelta))

    d.nomiss <- do.call("c", dimnames(tab.nomiss))
    d.miss <- do.call("c", dimnames(tab.miss))

    if(!identical(d.nomiss, d.miss)){

      missinglevels <- lapply(d.miss, function(d) {
        if(d %in% d.nomiss) return("") else return(d)
      })
      stop("Too many categories. ", paste("Cannot estimate probabilities for levels: ", paste(missinglevels, collapse = ", ")))

    }

    lookup <- prop.table(table(groups), margin = (1:length(dim(table(groups))))[-1]) # conditional probability
    impwith <- model.frame(formula[-2], mindelta)
    impwith <- lapply(impwith, as.character)

    if(is.factor(psdesign$augdata[[outname]])){
      outvect <- levels(psdesign$augdata[[outname]])
      } else {
      outvect <- sort(unique(psdesign$augdata[[outname]]))
    }
    psdesign$integration.models[[outname]]$model <- list(model = "nonparametric", args = arglist)

    psdesign$integration.models[[outname]]$cdf_sbarw <-
      function(S.1){

        sapply(S.1, function(s){

          dex <- c(list(lookup), list(s), as.list(impwith))
          res <- do.call("[", dex)

        })

      }
    psdesign$integration.models[[outname]]$icdf_sbarw <-
      function(U.1){

        dex <- c(list(lookup), list(TRUE), as.list(impwith))
        res <- do.call("[", dex)

        t(sapply(1:ncol(res), function(i){

          sample(outvect, size = length(U.1), prob = res[, i], replace = TRUE)

        }))

      }

    psdesign
  }

  class(rval) <- c("ps", "integration")
  rval

}




#' Bivariate normal integration models for the missing S(1)
#'
#' This model assumes that the pair [S(1), W] is bivariate normal, where W is
#' the BIP. The means, standard deviations, and correlation are estimated or
#' fixed before calling this function. Then the conditional normal formula is
#' applied in order to get the distribution of S(1) | W. That distribution is
#' used to integrate over the missing S(1) values. This method requires a BIP in the
#' design.
#'
#' @param x, expression identifying the variable to be integrated. Typically this is S.1 or S.0
#' @param mu, means of the pair of surrogates, missing one first
#' @param sd, standard deviations of the pair, missing one first
#' @param rho, the correlation between X1 and X2
#'
#' @export

integrate_bivnorm <- function(x = "S.1", mu = c(0, 0), sd = c(1, 1), rho = .2){

  outname <- x
  arglist <- as.list(match.call())

  rval <- function(psdesign){


    if(!"integration.models" %in% names(psdesign)) psdesign$integration.models <- NULL

    psdesign$integration.models[[outname]]$model <- list(model = "bivnorm", args = arglist)

    missdex <- !is.na(get(outname, psdesign$augdata))
    mindelta <- subset(psdesign$augdata, !missdex)

    if(is.factor(mindelta[, x])){

      stop("Parametric integration models not valid for factors. Use integrate_nonparametric instead.")

    }

    vmu <- mu[1] + sd[1]/sd[2] * rho * (mindelta$BIP - mu[2])
    vsd <- (1 - rho^2) * sd[1]^2




    psdesign$integration.models[[outname]]$cdf_sbarw <-
      function(s1){

        sapply(s1, function(s) pnorm(s, mean = vmu, sd = vsd))

      }
    psdesign$integration.models[[outname]]$icdf_sbarw <-
      function(U.1){

        sapply(U.1, function(s) qnorm(s, mean = vmu, sd = vsd))

      }

    psdesign
  }

  class(rval) <- c("ps", "integration")
  rval

}


#' Semiparametric integration model using the location-scale model
#'
#' @param formula.location Formula specifying the integration model for the location component of the surrogate
#'   under treatment. Generally the candidate surrogate will be on the left side
#'   in the formula, and the BIP or BIPs will be on the right side
#' @param formula.scale Formula specifying the integration model for the scale component of the surrogate
#'   under treatment. Generally the candidate surrogate will be on the left side
#'   in the formula, and the BIP or BIPs will be on the right side
#' @param ... Other parameters passed to \link{sp_locscale}
#'
#' @export

integrate_semiparametric <- function(formula.location, formula.scale, ...){
  arglist <- as.list(match.call())

  rval <- function(psdesign){

    if(!"integration.models" %in% names(psdesign)) psdesign$integration.models <- NULL

    outname <- paste(formula.location[[2]])

    missdex <- !is.na(get(outname, psdesign$augdata))

    if(is.factor(psdesign$augdata[, outname])){

      stop("Semiparametric integration models not valid for factors. Use integrate_nonparametric instead.")

    }

    fit <- sp_locscale(formula.location = formula.location, formula.scale = formula.scale,
                       data = psdesign$augdata[missdex, ], weights = psdesign$augdata$cdfweights[missdex], ...)

    stopifnot(fit$converge)

    psdesign$integration.models[[paste(formula.location[[2]])]]$model <- list(model = "semiparametric", args = arglist)

    mindelta <- subset(psdesign$augdata, !missdex)

    psdesign$integration.models[[paste(formula.location[[2]])]]$cdf_sbarw <-
      function(S.1){

        mu <- model.matrix(formula.location[-2], data = mindelta) %*% fit$beta.l
        sig <- exp(model.matrix(formula.scale[-2], data = mindelta) %*% fit$beta.s)

        sapply(S.1, function(s){
          sapply((s - mu)/sig, function(s.in) mean(fit$resid < s.in))
        })

      }
    psdesign$integration.models[[paste(formula.location[[2]])]]$icdf_sbarw <-
      function(U.1){

        mu <- model.matrix(formula.location[-2], data = mindelta) %*% fit$beta.l
        sig <- exp(model.matrix(formula.scale[-2], data =mindelta) %*% fit$beta.s)

        U.new <- sample(fit$resid, size = length(U.1), replace = TRUE)
        sapply(U.new, function(u) u * sig + mu)

      }

    psdesign
  }

  class(rval) <- c("ps", "integration")
  rval

}


#' Fit the semi-parametric location-scale model
#'
#' This estimates the location-scale model as described in Heagerty and Pepe (1999) using the Newton-Raphson method. The location and scale formulas must have the same outcome, but they may have different predictors.
#'
#' @param formula.location Formula specifying the model for the location
#' @param formula.scale Formula specifying the model for the scale
#' @param data Data used to estimate the model
#' @param weights Weights applied to the estimating equations
#' @param tol Convergence tolerance
#' @param maxit Maximum number of iterations
#'
#' @return A list containing the parameter estimates, the convergence indicator, and residuals
#'
#' @export

sp_locscale <- function(formula.location, formula.scale, data, weights, tol = 1e-6, maxit = 100){

  stopifnot(identical(formula.location[[2]], formula.scale[[2]]))


  mf.l <- model.frame(formula.location, data = data)
  mf.s <- model.frame(formula.scale, data = data)

  stopifnot(all.equal(model.response(mf.l), model.response(mf.s)))

  resp <- model.response(mf.l)
  if(missing(weights)) weights <- rep(1, length(resp))
  desmat.l <- model.matrix(formula.location, data = data)
  desmat.s <- model.matrix(formula.scale, data = data)

  nmu <- ncol(desmat.l)
  nsig <- ncol(desmat.s)
  nobs <- nrow(desmat.l)


  fitU <- glm(resp ~ desmat.l - 1, weights=weights)
  eta <- fitU$coef
  mu1 <- fitU$fitted.value

  res2 <- (resp - mu1)^2

  fitW <- glm(res2 ~ desmat.s - 1, weights = weights, family=quasi(variance="mu", link="log"))

  delta<-fitW$coef/2

  beta <- c(eta,delta)

  i <- 1
  for(i in 1:maxit){

    beta.l <- beta[1:nmu]
    beta.s <- beta[-c(1:nmu)]

    mu <- desmat.l %*% beta.l
    sig2 <- exp(desmat.s %*% beta.s)^2

    eeq.l <- matrix(rep(weights * (resp - mu) / sig2, nmu), byrow = FALSE, ncol = nmu)
    eeq.s <- matrix(rep(weights * ((resp - mu)^2 - sig2) / sig2, nsig), byrow = FALSE, ncol = nsig)

    eeq.all <- c(colMeans(desmat.l * eeq.l),
      colMeans(desmat.s * eeq.s))

    matrixUU <- desmat.l / matrix(rep(sig2, nmu), byrow=FALSE, nrow = nobs)
    H.UU<--t(desmat.l * weights) %*% matrixUU / nobs

    matrixUW <- desmat.s * matrix(rep((resp - mu) / sig2, nsig), byrow = FALSE, nrow = nobs)
    H.UW <- -2 * t(desmat.l * weights) %*% matrixUW / nobs

    matrixWW <- desmat.s * matrix(rep((resp - mu)^2 / sig2, nsig), byrow = FALSE, nrow = nobs)
    H.WW <- -2 * t(desmat.s * weights) %*% matrixWW / nobs

    H<-rbind(cbind(H.UU,H.UW),cbind(t(H.UW),H.WW))

    beta.0 <- beta - c(solve(H) %*% eeq.all)

    if(sum((beta.0 - beta)^2) < tol){
      break
    }
    beta <- beta.0

    i <- i + 1

  }

  beta <- beta.0
  converge <- i < maxit

  resid <- c((resp - desmat.l %*% beta[1:nmu])/exp(desmat.s %*% beta[-c(1:nmu)]))

  list(beta.l = beta[1:nmu], beta.s = beta[-c(1:nmu)], converge = converge, resid = resid)

}


