#' Fit cosinor model
#'
#' Given an outcome and time variable, fit the cosinor model with optional
#' covariate effects.
#'
#' @param formula Forumla specifying the model. Indicate the time variable with
#'   \code{time()} and covariate effects on the amplitude and acrophase with
#'   \code{amp.acro()}. See details for more information.
#' @param period Length of time for a complete period of the sine curve.
#' @param data Data frame where variable can be found
#' @param na.action What to do with missing data
#'
#' @details This defines special functions that are used in the formula to
#'   indicate the time variable and which covariates effect the amplitude. To
#'   indicate the time variable wrap the name of it in the function
#'   \code{time()}. To indicate a variable which affects the
#'   acrophase/amplitude, wrap the name in \code{amp.acro()}. This will then do
#'   all the tranformations for you. See examples for usage.
#'
#' @examples
#'
#' cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind)
#'
#' @references Tong, YL. Parameter Estimation in Studying Circadian Rhythms, Biometrics (1976). 32(1):85--94.
#'
#'
#' @export
#'



cosinor.lm <- function(formula, period = 12,
                       data, na.action = na.omit){

 # build time tranformations

  Terms <- terms(formula, specials = c("time", "amp.acro"))

  stopifnot(attr(Terms, "specials")$time != 1)
  varnames <- get_varnames(Terms)
  timevar <- varnames[attr(Terms, "specials")$time - 1]

  data$rrr <- cos(2 * pi * data[,timevar] / period)
  data$sss <- sin(2 * pi * data[,timevar] / period)

  spec_dex <- unlist(attr(Terms, "special")$amp.acro) - 1
  mainpart <- c(varnames[c(-spec_dex, - (attr(Terms, "special")$time - 1))], "rrr", "sss")
  acpart <- paste(sort(rep(varnames[spec_dex], 2)), rep(c("rrr", "sss"), length(spec_dex)), sep = ":")
  newformula <- as.formula(paste(rownames(attr(Terms, "factors"))[1],
                           paste(c(mainpart, acpart), collapse = " + "), sep = " ~ "))

  fit <- lm(newformula, data, na.action = na.action)

  mf <- fit

  r.coef <- c(FALSE, as.logical(attr(mf$terms, "factors")["rrr",]))
  s.coef <- c(FALSE, as.logical(attr(mf$terms, "factors")["sss",]))
  mu.coef <- c(TRUE, ! (as.logical(attr(mf$terms, "factors")["sss",]) |
                          as.logical(attr(mf$terms, "factors")["rrr",])))

  beta.s <- mf$coefficients[s.coef]
  beta.r <- mf$coefficients[r.coef]

  groups.r <- c(beta.r["rrr"], beta.r["rrr"] + beta.r[which(names(beta.r) != "rrr")])
  groups.s <- c(beta.s["sss"], beta.s["sss"] + beta.s[which(names(beta.s) != "sss")])

  amp <- sqrt(groups.r^2 + groups.s^2)
  names(amp) <- gsub("rrr", "amp", names(beta.r))

  acr <- atan(groups.s / groups.r)
  names(acr) <-  gsub("sss", "acr", names(beta.s))
  coef <- c(mf$coefficients[mu.coef], amp, acr)

  structure(list(fit = fit, Call = match.call(), Terms = Terms, coefficients = coef, period = period), class = "cosinor.lm")

}

#' Print cosinor model
#'
#' Given an outcome and time variable, fit the cosinor model with optional covariate effects.
#'
#' @param x cosinor.lm object
#' @param ... passed to summary
#'
#'
#' @export
#'


print.cosinor.lm <- function(x, ...){


  cat("Call: \n")
  print(x$Call)
  cat("\n Raw Coefficients: \n")
  print(x$fit$coefficients)
  cat("\n Transformed Coefficients: \n")
  t.x <- x$coefficients
  names(t.x) <- update_covnames(names(t.x))
  print(t.x)

}

#' Fit cosinor model
#'
#' Given an outcome and time variable, fit the cosinor model with optional covariate effects.
#'
#' @param formula Forumla specifying the model. Indicate the time variable with \code{time()} and covariate effects on the
#' amplitude and acrophase with \code{amp.acro()}. See details.
#' @param ... other arguments
#'
#' @details This defines special functions that are used in the formula to indicate the time variable
#' and which covariates effect the amplitude. To indicate the time variable wrap the name of it in the function
#' \code{time()}. To indicate a variable which affects the acrophase/amplitude, wrap the name in
#' \code{amp.acro()}. This will then do all the tranformations for you. See examples for usage.
#'
#' @examples
#'
#' cosinor.lm(Y ~ time(time) + X + amp.acro(X), data = vitamind)
#'
#' @export
#'


cosinor.lm.default <- function(formula, ...){

  UseMethod("cosinor.lm")

}

#' Extract variable names from terms object, handling specials
#'
#' @param Terms a terms object
#'
#' @keywords Internal
#'

get_varnames <- function(Terms){

  spec <- names(attr(Terms, "specials"))
  tname <- attr(Terms, "term.labels")

  dex <- unlist(sapply(spec, function(sp){

    attr(Terms, "specials")[[sp]] - 1

  }))

  tname2 <- tname
  for(jj in spec){

    gbl <- grep(paste0(jj, "("), tname2, fixed = TRUE)
    init <- length(gbl) > 0
    if( init ){
    jlack <- gsub(paste0(jj, "("), "", tname2, fixed = TRUE)
    tname2[gbl] <- substr(jlack[gbl], 1, nchar(jlack[gbl]) - 1)
    }

    }

   tname2

}

#' Replace covariate names with descriptive text
#'
#' @param names Coefficient names to update
#'
#' @export
#'

update_covnames <- function(names){

  covnames <- grep("(amp|acr|Intercept)", names, invert = TRUE, value = TRUE)

  lack <- names
  for(n in covnames){
  lack <- gsub(paste0(n, ":"), paste0("[", n, " = 1]:"), lack)
  lack <- gsub(paste0("^", n, "$"), paste0("[", n, " = 1]"), lack)
  }
  lack
}
