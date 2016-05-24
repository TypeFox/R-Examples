#' Combined List Estimator
#'
#' This function implements the combined list estimator described in Aronow, Coppock, Crawford, and Green (2015): Combining List Experiment and Direct Question Estimates of Sensitive Behavior Prevalence
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. Should be of the form Y ~ T + X1 + X2, where Y is the list response, T is the treatment indicator, and X1, X2, etc are pretreatment covariates. It is recommended that T be a numeric variable whose values are 0 for subjects in control and 1 for subjects in treatment.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which combined.list is called. It is good practice to include all variables used in the estimation (list response, treatment indicator, direct response, and optional pre-treatment covariates) in a dataframe rather than calling data from the global environent.
#' @param treat a character string giving the name of the treatment variable. Defaults to "treat".
#' @param direct a character string giving the name of the direct response variable. Defaults to "direct". The direct response variable itself must only contain the values 0 and 1, where 1 refers to subjects who answered "Yes" to the direct question.
#' @return a list containing conventional, direct, and combined prevalence estimates with associated standard errors as well as the results of two placebo tests.
#' @examples
#' # Load data from Aronow, Coppock, Crawford, and Green (2015)
#' data("combinedListExps")
#' 
#' # complete case analysis
#' combinedListExps <- na.omit(combinedListExps)
#' # Conduct estimation without covariate adjustment
#' out.1 <- combinedListDirect(list1N ~ list1treat, 
#'                             data = subset(combinedListExps, directsfirst==1), 
#'                             treat = "list1treat", direct = "direct1")
#' summary(out.1)
#'
#' # Conduct estimation with covariate adjustment
#' out.2 <- combinedListDirect(list1N ~ list1treat + gender + 
#'                             ideology + education + race, 
#'                             data = subset(combinedListExps, directsfirst==1), 
#'                             treat = "list1treat", direct = "direct1")
#' summary(out.2)
combinedListDirect <- function(formula, data = parent.frame(), treat="treat", direct="direct"){

  comblist.call <- match.call()

  # set up data frame, with support for standard and modified responses
  mf <- match.call(expand.dots = FALSE)

  # make all other call elements null in mf <- NULL in next line
  # you can remove all the options I have, but make sure to make NULL any option that remains like below
  mf$treat <- mf$direct <- NULL
  mf[[1]] <- as.name("model.frame")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)

  # define design, response data frames
  x.all <- model.matrix.default(attr(mf, "terms"), mf)
  #x.all <- x.all[,-1]
  x.all <- x.all[,setdiff(colnames(x.all), treat)]
  x.all <- as.matrix(x.all)
  y.all <- model.response(mf)
  t.all <- data[,paste(treat)]
  d.all <- data[,paste(direct)]

  # list-wise missing deletion
  na.x <- apply(is.na(x.all), 1, sum)
  na.y <- is.na(y.all)
  na.t <- is.na(t.all)
  na.d <- is.na(d.all)

  keepers <- (na.x==0 & na.y==0 & na.t==0 & na.d==0)

  x.nona <- x.all[keepers, , drop=FALSE]
  y.nona <- y.all[keepers]
  t.nona <- t.all[keepers]
  d.nona <- d.all[keepers]

  if(length(unique(t.nona))!=2){stop("The treatment variable must only contain two unique values.")}
  if(!all(d.nona %in% c(0,1))){stop("The direct response variable must be a numeric variable whose values are 0 or 1.")}

  if(class(t.nona) == "factor") {

    levels(t.nona) <- tolower(levels(t.nona))

    if (length(which(levels(t.nona) == "control")) == 1) {
      t.nona <- relevel(t.nona, ref = "control")
    } else {
      warning("Note: using the first level of the treatment variable as the control condition, but it is not labeled 'control'.")
    }

    t.nona <- as.numeric(t.nona) - 1

  } else {
    if(!all(t.nona %in% c(0,1))){stop("The treatment variable must be a numeric variable whose values are 0 or 1.")}
  }

  ## so that the output data has the same dimension as x.all and y.all
  data <- data[keepers, , drop = FALSE]

  # Direct Estimate
  d.bar <- mean(d.nona)
  direct.var.est <- var(d.nona)/length(d.nona)

  # Conventional List Estimate
  conv.fit <- lm(y.nona ~ t.nona + x.nona)
  conv.est <- coef(conv.fit)[2]
  conv.resids <- residuals(conv.fit)
  conv.var.t.1 <- var(conv.resids[t.nona==1])
  conv.var.t.0 <- var(conv.resids[t.nona==0])
  conv.var.est <- conv.var.t.1/sum(t.nona) + conv.var.t.0/sum(1-t.nona)

  # Combined List Estimate
  comb.fit <- lm(y.nona ~ t.nona + x.nona, subset=d.nona==0)
  intermediate.est <- coef(comb.fit)[2]
  mu.hat <- d.bar + (1-d.bar)*(intermediate.est)

  comb.resids <- residuals(comb.fit)

  comb.var.t.1 <- var(comb.resids[t.nona[d.nona==0]==1])
  comb.var.t.0 <- var(comb.resids[t.nona[d.nona==0]==0])

  gamma.hat <- mean(t.nona)
  comb.var.est <- ((((1-mu.hat)^2)/(1-d.bar))*d.bar +
                (1-d.bar)*(comb.var.t.1/gamma.hat + comb.var.t.0/(1-gamma.hat)))/length(t.nona)

  # Placebo Test I

  placebo.I.fit <- lm(y.nona ~ t.nona + x.nona, subset=d.nona==1)
  placebo.I.est <- as.numeric(coef(placebo.I.fit)[2])
  placebo.I.resids <- residuals(placebo.I.fit)
  placebo.I.var <- (var(placebo.I.resids[t.nona[d.nona==1]==1]))/sum(t.nona[d.nona==1]==1) +
    (var(placebo.I.resids[t.nona[d.nona==1]==0]))/sum(t.nona[d.nona==1]==0)
  placebo.I.se <- as.numeric(sqrt(placebo.I.var))
  placebo.I.p <- as.numeric(2*pnorm(- abs(placebo.I.est -1)/ placebo.I.se))
  placebo.I.n <- as.numeric(sum(d.nona==1))

  placebo.I.list <- list(estimate = placebo.I.est, se = placebo.I.se,
                         p = placebo.I.p, n = placebo.I.n)

  # Placebo Test II

  placebo.II.fit <- lm(d.nona ~ t.nona + x.nona)
  placebo.II.est <- as.numeric(coef(placebo.II.fit)[2])
  placebo.II.resids <- residuals(placebo.II.fit)
  placebo.II.var <- var(placebo.II.resids[t.nona==1])/(sum(t.nona==1)) +
                    var(placebo.II.resids[t.nona==0])/(sum(t.nona==0))
  placebo.II.se <- as.numeric(sqrt(placebo.II.var))
  placebo.II.p <- as.numeric(2*pnorm(- abs(placebo.II.est)/ placebo.II.se))
  placebo.II.n <- as.numeric(length(t.nona))

  placebo.II.list <- list(estimate = placebo.II.est, se = placebo.II.se,
                         p = placebo.II.p, n = placebo.II.n)


  return.object <- list(comb.est=mu.hat, comb.se=sqrt(comb.var.est),
                        direct.est=d.bar, direct.se=sqrt(direct.var.est),
                        conv.est=conv.est, conv.se=sqrt(conv.var.est),
                        placebo.I = placebo.I.list,
                        placebo.II = placebo.II.list,
                        call = match.call(),
                        data = data, x = x.nona, y = y.nona, treat = t.nona, direct=d.nona)

  class(return.object) <- "comblist"

  return(return.object)

}

print.comblist <- function(x, ...){

  cat("\n Combined List Estimates \n\nCall: ")

  dput(x$call)

  cat("\n Prevalence estimate\n")

  tb <- as.matrix(c(x$comb.est, x$comb.se))
  rownames(tb) <- c("Estimate", "Standard Error")
  colnames(tb) <- "Prevalence"
  print(tb)

}

summary.comblist <- function(object, ...) {
  structure(object, class = c("summary.comblist", class(object)))
}

print.summary.comblist <- function(x, ...){

  cat("\n Combined List Estimates \n\nCall: ")

  dput(x$call)

  cat("\n Prevalence estimates\n")

  estimates <- c(x$comb.est, x$direct.est, x$conv.est)
  ses <- c(x$comb.se, x$direct.se, x$conv.se)

  tb <- rbind(estimates, ses)
  colnames(tb) <- c("Combined", "Direct", "Conventional")
  rownames(tb) <- c("Estimate", "Standard Error")

  print(tb)

  cat("\n Placebo Test I
       Beta is the conventional list experiment estimate among those who answer 'Yes' to the direct question.
       Ho: beta = 1
       Ha: beta != 1
      \n ")
  pItb <- t(as.matrix(unlist(x$placebo.I)))
  colnames(pItb) <- c("Estimate", "SE", "p", "n")
  rownames(pItb) <- "beta"
  print(pItb)


  cat("\n Placebo Test II
       Delta is the average effect of the receiving the treatment list on the direct question response.
       Ho: delta = 0
       Ha: delta != 0
      \n ")
  pIItb <- t(as.matrix(unlist(x$placebo.II)))
  colnames(pIItb) <- c("Estimate", "SE", "p", "n")
  rownames(pIItb) <- "delta"
  print(pIItb)
}

#' Five List Experiments with Direct Questions
#'
#' A dataset containing the five list experiments in 
#' Aronow, Coppock, Crawford, and Green (2015)
#'
#'
#' @docType data
#' @keywords datasets
#' @name combinedListExps
#' @usage data(combinedListExps)
#' @format A data frame with 1023 observations and 23 variables
NULL

