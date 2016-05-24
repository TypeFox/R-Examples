# Here we have the environmental fit functions.

#' Environmental fit for a single dataset
#'
#' \code{all_environmental_fit} calculates the best expressions for colonization
#' and extinction rates given their dependency on environmental variables. \cr
#' \code{greedy_environmental_fit} calculates expressions for colonization and
#' extinction rates given their dependency on environmental variables using a
#' greedy algorithm. \cr \code{custom_environmental_fit} calculates the m.l.e.
#' of the parameters describing the relationship between colonization and
#' extinction rates and environmental variables.
#'
#' @param dataset A single dataset.
#' @param vector A vector indicating the columns with presence-absence data.
#' @param env The names of the environmental variables to be considered.
#' @param c Tentative colonization rate.
#' @param e Tentative extinction rate.
#' @param aic Tentative AIC to be improved by the optimizer.
#' @param params A vector with priors of the parameters in exp1 and exp2.
#' @param exp1 Expression for colonization.
#' @param exp2 Expression for extinction.
#'
#' @return A list with three components: a expression for colonization, a
#'   expression for extinction and the output of the optimization function, or
#'   the output of the optimization function in the custom environmental fit.
#'
#' @details \code{all_environmental_fit} calculates all the combinations of
#'   parameters, that increase exponentially with the number of parameters. We
#'   advise to keep low the number of parameters. \cr
#'   \code{greedy_environmental_fit} adds sequentially environmental variables
#'   to the expressions of colonization and extinction rates and fix one at a
#'   time until termination, when only adding one variable does not improve the
#'   AIC of the last accepted model.
#' @note AIC is recomended to be higher than the AIC of the most simple model
#'   (i.e. not including environmental variables).
#'
#' @seealso \code{\link{rates_calculator}}
#'
#' @examples \dontrun{
#' all_environmental_fit(idaho[[1]],3:23,c("idaho[[2]]$TOTAL.ppt",
#' "idaho[[2]]$ANNUAL.temp"),0.13,0.19,100000)
#' greedy_environmental_fit(idaho[[1]],3:23,c("idaho[[2]]$TOTAL.ppt",
#' "idaho[[2]]$ANNUAL.temp"),0.13,0.19,100000)
#' }
#' custom_environmental_fit(idaho[[1]], 3:23, c(-0.00497925, -0.01729602,
#' 0.19006501, 0.93486956), expression(params[1] * idaho[[2]]$TOTAL.ppt[i] +
#' params[3]), expression(params[2] * idaho[[2]]$ANNUAL.temp[i] + params[4]))
#'
#' @export
all_environmental_fit <- function(dataset, vector, env, c, e, aic) {
  co <- incounts(dataset, vector)
  t <- times(dataset, vector)
  env2 <- c(paste0(env, "[i]"), paste0(env, "[i]"))
  for (i in 1:length(env2)) {
    index <- (utils::combn(length(c(env, env)), i))
    for (j in 1:ncol(index)) {
      cexp <- c()
      eexp <- c()
      for (k in 1:nrow(index)) {
        if (index[k, j] > length(env)) {
          eexp <- c(eexp, index[k, j])} else {
            cexp <- c(cexp, index[k, j])
          }
      }
      n <- length(cexp)

      if (is.null(cexp)) {
        cexp <- paste0("params[", i + 1, "]")
      } else {

        cexp <- paste0("params[", 1:n, "] * ", env2[cexp], collapse = " + ")
        cexp <- paste0(cexp, " + params[", i + 1,"]")
      }

      if (is.null(eexp)) {
        eexp <- paste0("params[", i + 2,"]")
      } else {
        eexp <- paste0("params[", (n + 1):i,"] * ", env2[eexp], collapse =
                         " + ")
        eexp <- paste0(eexp, " + params[", i + 2,"]")
      }

      cexp <- parse(text = cexp)
      eexp <- parse(text = eexp)
      params <- c(rep(0, i), c, e)

      r <- stats::optim(params, dtlsolverenv, t = t, n00 = co[, 3], n10 = co[, 2],
                 n01 = co[, 1], n11 = co[, 4], exp1 = cexp, exp2 = eexp,
                 control = list(maxit = 5000000, fnscale = - 1))
           aic2 <- - 2 * r$value + 2 * (i + 2)
      if (aic2 < aic) {
        aic <- aic2
        print(list(cexp, eexp, r))
      }
    }
  }
}

#' @rdname all_environmental_fit
#' @export
custom_environmental_fit <- function(dataset, vector, params, exp1, exp2) {
  co <- incounts(dataset, vector)
  t <- times(dataset, vector)
  r <- stats::optim(params, dtlsolverenv, t = t, n00 = co[, 3], n10 = co[, 2],
             n01 = co[, 1], n11 = co[, 4], exp1 = exp1, exp2 = exp2,
             control = list(maxit = 5000000, fnscale = - 1))
  r
}

#' @rdname all_environmental_fit
#' @export
greedy_environmental_fit <- function(dataset, vector,env, c, e, aic) {
  co <- incounts(dataset, vector)
  t <- times(dataset, vector)
  env2 <- c(paste0(env, "[i]"), paste0(env, "[i]"))
  selected <- c()
  for (i in 1:length(env2)) {
    event <- F
    index <- 1:length(env2)
    if (i >= 3) {
      selection <- c()
      for (xx in 1:nrow(selected)) {
        selection <- c(selection, selected[xx, 1])
      }
    } else {
      selection <- selected
    }
    if (i >= 2) {
      index <- index[-selection]
      index <- rbind(index, matrix(rep(selected, length(index)),
                                   ncol = length(index)))
    } else {
      index <- matrix(index, nrow = 1)
    }
    jota <- NULL
    for (j in 1:ncol(index)) {
      cexp <- c()
      eexp <- c()
      for (k in 1:nrow(index)) {
        if (index[k, j] > length(env)) {
          eexp <- c(eexp, index[k, j])
        } else {
          cexp <- c(cexp, index[k, j])
        }
      }
      n <- length(cexp)
      if (is.null(cexp)) {
        cexp <- paste0("params[", i + 1, "]")
      } else {
        cexp <- paste0("params[", 1:n, "] * ", env2[cexp], collapse = " + ")
        cexp <- paste0(cexp, " + params[", i + 1, "]")
      }
      if (is.null(eexp)) {
        eexp <- paste0("params[", i + 2, "]")
      } else {
        eexp <- paste0("params[", (n + 1):i, "] * ", env2[eexp],
                       collapse = " + ")
        eexp <- paste0(eexp, " + params[", i + 2, "]")
      }
      cexp <- parse(text = cexp)
      eexp <- parse(text = eexp)
      params <- c(rep(0, i), c, e)
      r<-stats::optim(params, dtlsolverenv, t = t, n00 = co[, 3], n10 = co[, 2],
               n01 = co[, 1], n11 = co[, 4], exp1 = cexp, exp2 = eexp,
               control = list(maxit = 5000000, fnscale = -1), hessian = F)
      aic2 <- -2 * r$value + 2 * (i + 2)
      if (aic2 < aic) {
        jota <- index[1, j]
        aic <- aic2
        event <- T
        print(list(cexp, eexp, r))
      }
    }
    if (event) selected <- rbind(selected, jota) else {
      print("Last iteration found the solution for this algorithm.")
      break
    }
  }
}


#' Colonization and extinction rates calculator for expressions.
#'
#' \code{rates_calculator} Calculate colonization and extinction rates depending
#' of their expressions.
#'
#' @param params A vector with priors of the parameters in exp1 and exp2.
#' @param exp1 Expression for colonization.
#' @param exp2 Expression for extinction.
#' @param t Number of colonization and extinction pairs required.
#'
#' @return A matrix with the colonization and extinction rates.
#'
#' @examples rates_calculator(c(-0.00497925, -0.01729602, 0.19006501,
#' 0.93486956), expression(params[1] * idaho[[2]]$TOTAL.ppt[i] + params[3]),
#' expression(params[2] * idaho[[2]]$ANNUAL.temp[i] + params[4]), 21)
#'
#' @seealso \code{\link{all_environmental_fit}}
#'
#'
#' @export
rates_calculator <- function(params, exp1, exp2, t) {
  c <- vector(length = t)
  e <- vector(length = t)
  out <- data.frame()
  for (i in 1:t) {
    c[i] <- eval(exp1)
    e[i] <- eval(exp2)
  }
  out <- cbind(c, e)
  out
}

dtlsolverenv <- function(params, t, n00, n10, n01, n11, exp1, exp2){
  c <- exp1
  e <- exp2
  ll <- vector(length = length(t))

  for (i in 1:length(t)) {
    ceval <- eval(c)
    eeval <- eval(e)
    ll[i] <- n00[i] * log(1 - (ceval / (eeval + ceval)) * (1 - exp( - (eeval
                                                           + ceval) * t[i]))) +
      n01[i] * log((eeval / (eeval + ceval)) * (1 - exp( - (eeval + ceval) *
                                                                      t[i]))) +
      n10[i] * log((ceval / (eeval + ceval)) * (1 - exp( - (eeval + ceval) *
                                                                      t[i]))) +
      n11[i] * log(1 - (eeval / (eeval + ceval)) * (1 - exp( - (eeval +
                                                               ceval) * t[i])))

  }
  sum(ll)
}
