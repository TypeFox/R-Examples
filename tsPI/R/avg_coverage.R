#' Compute the average coverage of the prediction intervals computed by
#' naive plug-in method and \code{\link{arima_pi}}
#'
#' Computes expected coverage probabilities of the prediction intervals of
#' ARMA process by simulating time series from the known model.
#'
#'
#' @export
#' @seealso \code{\link{arima_pi}}.
#
#' @param phi vector containing the AR parameters
#' @param theta vector containing the MA parameters
#' @param d degree of differencing
#' @param n length of the time series
#' @param n_ahead length of the forecast horizon
#' @param nsim number of simulations used in importance sampling
#' @param nsim2 number of simulations used in computing the expected coverage
#' @param level desired coverage probability of the prediction intervals
#' @param prior prior to be used in importance sampling. Multiple choices are allowed.
#' @param return_all_coverages return raw results i.e. coverages for each simulations. When \code{FALSE} (default), summary statistics are returned.
#' @param ... additional arguments to \code{\link{arima_pi}}.
#' @return a list containing the coverage probabilities
#' @examples
#'
#'
#' \dontrun{
#' set.seed(123)
#' # takes a while, notice se, increase nsim2 to get more accurate results
#' avg_coverage_arima(phi = 0.9, n = 50, n_ahead = 10, nsim2 = 100)
#'
#' avg_coverage_arima(phi = 0.9, theta = -0.6, n = 50, n_ahead = 10, nsim2 = 100)
#' }
avg_coverage_arima <- function(phi = NULL, theta = NULL, d = 0, n, n_ahead = 1,
  nsim2, nsim = 100, level = 0.95, prior = "uniform", return_all_coverages = FALSE, ...){

  p <- length(phi)
  q <- length(theta)
  prior <- match.arg(prior,
    c("uniform", "approx_joint_jeffreys", "approx_marginal_jeffreys",
      "exact_joint_jeffreys", "exact_marginal_jeffreys", "custom"))

  if (p > 0 && !all(Mod(polyroot(c(1, -phi))) > 1))
    stop("True model is not stationary.")
 # if (q > 0 && !all(Mod(polyroot(c(1, -theta))) > 1))
 #   stop("True model is not invertible.")

  covprobs<-array(0, c(n_ahead, nsim2, 2))

  dimnames(covprobs)[[3]] <- c("plug-in", prior)
  count <- count2 <- 0
  for (i in 1:nsim2) {
    fit <- NULL
    fit$code <- 1
    x <- arima.sim(n = n - d, list(ar = phi, ma = theta, order = c(p, d, q)))

    fit <- try(suppressWarnings(arima(x, order = c(p,d,q),method="CSS-ML",optim.control=list(maxit=1000))),TRUE)
    # try again
    if (inherits(fit, "try-error") || fit$code!=0 || sum(diag(fit$var.coef)<1e-7)!=0 || ifelse(p>0,!all(Mod(polyroot(c(1, -fit$coef[1:p]))) > 1), FALSE)) {
      fit <- try(suppressWarnings(arima(x,order=c(p,d,q), method = "ML",optim.control=list(maxit=1000))), TRUE)
      if(inherits(fit, "try-error") || fit$code!=0 || sum(diag(fit$var.coef)<1e-7)!=0 || ifelse(p>0,!all(Mod(polyroot(c(1, -fit$coef[1:p]))) > 1),FALSE)){
        #still fails, give up
        count <- count + 1
        covprobs[, i, 2] <- covprobs[, i, 1] <- NA
        next
      }
    }

    model <- SSModel(x ~ SSMarima(ar = if (p > 0) fit$coef[1:p] else NULL,
      ma = if(q > 0) fit$coef[(p + 1):(p + q)] else NULL, d = d, Q = fit$sigma), H = 0)
    newdata <- SSModel(rep(NA, n_ahead) ~ SSMarima(ar = if(p > 0) fit$coef[1:p] else NULL,
      ma = if(q > 0) fit$coef[(p + 1):(p + q)] else NULL,d = d,Q = fit$sigma), H = 0)
    if(d == 0){
      model$P1inf[1,1] <- 0
      model$a1[1] <- fit$coef[p+q+1]
      newdata$P1inf[1,1]<-0
      newdata$a1[1]<-fit$coef[p+q+1]
    }
    pred <- predict(model, newdata = newdata, level = level, interval = "prediction")



    model<-SSModel(x~SSMarima(ar=phi,ma=theta,d=d,Q=1),H=0)
    newdata<-SSModel(rep(NA,n_ahead)~SSMarima(ar=phi,ma=theta,d=d,Q=1),H=0)
    if(d==0){
      model$P1inf[1,1]<-0
      newdata$P1inf[1,1]<-0
    }
    true_pred<-predict(model,newdata=newdata,se=TRUE)


    ipi <- try(arima_pi(x,order=c(p,d,q),nsim = nsim,level = level, n_ahead = n_ahead,
      prior = prior, median = FALSE, se_limits = FALSE, ...),TRUE)
    if(!inherits(ipi, "try-error")){
      covprobs[, i,2] <- pnorm(q=ipi[,"upr"],mean=true_pred[,"fit"], sd=true_pred[,"se.fit"]) -
        pnorm(q=ipi[,"lwr"],mean=true_pred[,"fit"], sd=true_pred[,"se.fit"])
    } else {
      count2 <- count2 + 1
    }
    covprobs[, i,1] <- pnorm(q=pred[,3],mean=true_pred[,"fit"], sd=true_pred[,"se.fit"]) - pnorm(q=pred[,2],mean=true_pred[,"fit"], sd=true_pred[,"se.fit"])

  }
  if(count > 0)
    warning(paste0("There were",count, " cases where the arima.sim generated a series which caused the estimation of the model to fail. These cases were set as NA when computing coverage probabilities. " ))
  if(count2 > 0)
      warning(paste0("There were ",count2, " cases where the the importance sampling method generated error. The coverage probability for these cases were set to zero when computing coverage probabilities. " ))

  if(!return_all_coverages){
    out <- vector("list", 2)
    names(out) <- c("plugin", prior)
    for(i in 1:2){
      out[[i]] <- matrix(0, n_ahead, 8)
      colnames(out[[i]]) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.","Standard error", "Failures")
      rownames(out[[i]]) <- paste("n_ahead =", 1:n_ahead)
      for(j in 1:n_ahead)
        out[[i]][j,] <- c(summary(na.exclude(covprobs[j,,i])), sd(covprobs[j,,i], na.rm = TRUE)/sqrt(nsim2-count), sum(is.na(covprobs[j,,i]))+(i==2)*count2)
    }
    return(out)
  } else return(covprobs)
}
