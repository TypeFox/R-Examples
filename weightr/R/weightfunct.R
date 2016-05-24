#' Estimate the Vevea and Hedges (1995) Weight-Function Model
#'
#' This function allows the user to estimate the Vevea and Hedges (1995) weight-function model for publication bias.
#' @param effect a vector of meta-analytic effect sizes.
#' @param v a vector of meta-analytic sampling variances; needs to match up with the vector of effects, such that the first element in the vector of effect sizes goes with the first element in the vector of sampling variances, and so on.
#' @param steps a vector of p-value cutpoints. The default only distinguishes between significant and non-significant effects (p < 0.05).
#' @param mods defaults to \code{NULL}. A formula specifying the linear model.
#' @param weights defaults to \code{FALSE}. A vector of prespecified weights for p-value cutpoints to estimate the Vevea and Woods (2005) model.
#' @param fe defaults to \code{FALSE}. Indicates whether to estimate a fixed-effects model.
#' @param table defaults to \code{FALSE}. Indicates whether to print a table of the p-value intervals specified and the number of effect sizes per interval.
#' @importFrom stats model.matrix optim pchisq pnorm qnorm
#' @keywords weightr
#' @details This function allows meta-analysts to estimate both the
#' weight-function model for publication bias that was originally published in
#' Vevea and Hedges (1995) and the modified version presented in Vevea and Woods
#' (2005). Users can estimate both of these models with and without predictors and
#' in random-effects or fixed-effects situations.
#'
#' The Vevea and Hedges (1995) weight-function model is a tool for modeling publication
#' bias using weighted distribution theory. The model first estimates an unadjusted
#' fixed-, random-, or mixed-effects model, where the observed effect sizes are
#' assumed to be normally distributed as a function of predictors. This unadjusted
#' model is no different from the traditional meta-analytic model. Next, the Vevea
#' and Hedges (1995) weight-function model estimates an adjusted model that includes
#' not only the original mean model, fixed-, random-, or mixed-effects, but a series
#' of weights for any pre-specified p-value intervals of interest. This produces mean,
#' variance component, and covariate estimates adjusted for publication bias, as well
#' as weights that reflect the likelihood of observing effect sizes in each specified
#' interval.
#'
#' It is important to remember that the weight for each
#' estimated p-value interval must be interpreted relative to the first interval,
#' the weight for which is fixed to 1 so that the model is identified. In other
#' words, a weight of 2 for an interval indicates that effect sizes in that p-value
#' interval are about twice as likely to be observed as those in the first interval.
#' Finally, it is also important to remember that the model uses p-value cutpoints
#' corresponding to one-tailed p-values. This allows flexibility in the selection
#' function, which does not have to be symmetric for effects in the opposite direction;
#' a two-tailed p-value of 0.05 can therefore be represented as p < .025 or p > .975.
#'
#' After both the unadjusted and adjusted meta-analytic models are estimated, a
#' likelihood-ratio test compares the two. The degrees of freedom for this test are
#' equal to the number of weights being estimated. If the likelihood-ratio test is
#' significant, this indicates that the adjusted model is a better fit for the data,
#' and that publication bias may be a concern.
#'
#' To estimate a large number of weights for p-value intervals, the Vevea and Hedges
#' (1995) model works best with large meta-analytic datasets. It may have trouble
#' converging and yield unreliable parameter estimates if researchers, for instance,
#' specify a p-value interval that contains no observed effect sizes. However,
#' meta-analysts with small datasets are still likely to be interested in assessing
#' publication bias, and need tools for doing so. Vevea and Woods (2005)
#' attempted to solve this problem by adapting the Vevea and Hedges (1995) model to
#' estimate fewer parameters. The meta-analyst can specify p-value cutpoints,
#' as before, and specify corresponding fixed weights for those cutpoints. Then the
#' model is estimated. For the adjusted model, only the variance component and mean
#' model parameters are estimated, and they are adjusted relative to the fixed weights.
#' For example, weights of 1 for each p-value interval specified describes a situation
#' where there is absolutely no publication bias, in which the adjusted estimates are
#' identical to the unadjusted estimates. By specifying weights that depart from 1 over various p-value intervals, meta-analysts can
#' examine how various one-tailed or two-tailed selection patterns would alter their
#' effect size estimates. If changing the pattern of weights drastically changes
#' the estimated mean, this is evidence that the data may be vulnerable to
#' publication bias.
#'
#' For more information, consult the papers listed in the References section here.
#' Also, feel free to email the maintainer of \code{weightr} at kcoburn@ucmerced.edu.
#' The authors are currently at work on a detailed package tutorial, which we
#' hope to publish soon.
#' @export
#' @return The function returns a list containing the following components: \code{unadj_est}, \code{unadj_se}, \code{adj_est}, \code{adj_se}, \code{z_unadj}, \code{z_adj}, \code{p_unadj}, \code{p_adj}, \code{ci.lb_unadj}, \code{ci.ub_unadj}, \code{ci.lb_adj}, and \code{ci.ub_adj}.
#'
#' For each component, the order of relevant parameters is as follows: variance component, mean or linear coefficients, and weights. (Note that if \code{weights}
#' are specified using the Vevea and Woods (2005) model, no standard errors, p-values, z-values, or confidence intervals
#' are provided for the adjusted model, as these are no longer meaningful. Also note that the variance component is not reported for fixed-effects models.)
#' \describe{
#'    \item{\code{unadj_est}}{the unadjusted model estimates}
#'    \item{\code{unadj_se}}{standard errors for unadjusted model estimates}
#'    \item{\code{adj_est}}{the adjusted model estimates}
#'    \item{\code{adj_se}}{standard errors for adjusted model estimates}
#'    \item{\code{z_unadj}}{z-statistics for the unadjusted model estimates}
#'    \item{\code{z_adj}}{z-statistics for the adjusted model estimates}
#'    \item{\code{p_unadj}}{p-values for the unadjusted model estimates}
#'    \item{\code{p_adj}}{p-values for the adjusted model estimates}
#'    \item{\code{ci.lb_unadj}}{lower bounds of the 95\% confidence intervals for the
#'    unadjusted model estimates}
#'    \item{\code{ci.ub_unadj}}{upper bounds of the 95\% confidence intervals for the
#'    unadjusted model estimates}
#'    \item{\code{ci.lb_adj}}{lower bounds of the 95\% confidence intervals for the
#'    adjusted model estimates}
#'    \item{\code{ci.ub_adj}}{upper bounds of the 95\% confidence intervals for the
#'    adjusted model estimates}
#'    }
#'
#' @references Coburn, K. M. & Vevea, J. L. (2015). Publication bias as a function
#' of study characteristics. Psychological Methods, 20(3), 310.
#'
#' Vevea, J. L. & Hedges, L. V. (1995). A general linear model for
#' estimating effect size in the presence of publication bias. Psychometrika, 60(3),
#' 419-435.
#'
#' Vevea, J. L. & Woods, C. M. (2005). Publication bias in research synthesis:
#' Sensitivity analysis using a priori weight functions. Psychological Methods, 10(4),
#' 428-443.
#' @examples
#' \dontrun{
#' # Uses the default p-value cutpoints of 0.05 and 1:
#'
#' weightfunct(effect, v)
#'
#' # Estimating a fixed-effects model, again with the default cutpoints:
#'
#' weightfunct(effect, v, fe=TRUE)
#'
#' # Specifying cutpoints:
#'
#' weightfunct(effect, v, steps=c(0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50, 1.00))
#'
#' # Including a linear model, where moderators are denoted as 'mod1' and mod2':
#'
#' weightfunct(effect, v, mods=~mod1+mod2)
#'
#' # Specifying cutpoints and weights to estimate Vevea and Woods (2005):
#'
#' weightfunct(effect, v, steps=c(0.01, 0.05, 0.50, 1.00), weights=c(1, .9, .7, .5))
#'
#' # Specifying cutpoints and weights while including a linear model:
#'
#' weightfunct(effect, v, mods=~mod1+mod2, steps=c(0.05, 0.10, 0.50, 1.00), weights=c(1, .9, .8, .5))
#' }

weightfunct <- function(effect, v, steps=c(0.05,1.00), mods=NULL,
                        weights=NULL, fe=FALSE, table=FALSE){

  neglike_unadj <- function(pars) {
    if(fe == FALSE){
      vc = pars[1]
      beta = pars[2:(npred+2)]

      mn = XX%*%beta
      eta = sqrt(v + vc)
      b = 1/2 * sum(((effect-mn)/eta)^2)
      c = sum(log(eta))
    }
    else{
      beta = pars[1:(npred+1)]
      mn = XX%*%beta
      eta = sqrt(v)
      b = 1/2 * sum(((effect-mn)/eta)^2)
      c = sum(log(eta))
    }
    return(b+c)
  }

  neglike_adj <- function(pars) {
    if(fe == FALSE){
      vc = pars[1]
      beta = pars[2:(npred+2)]
      if(is.null(weights) == FALSE){
        w = weights
      }
      else{
        w = c(1,pars[(npred+3):( (nsteps - 2) + (npred+3) )])
      }
      contrib = log(w[wt])

      mn = XX%*%beta
      a = sum(contrib)
      eta = sqrt(v + vc)
    }
    else{
      beta = pars[1:(npred+1)]
      if(is.null(weights) == FALSE){
        w = weights
      }
      else{
        w = c(1,pars[(npred+2):( (nsteps - 2) + (npred+2) )])
      }
      contrib = log(w[wt])

      mn = XX%*%beta
      a = sum(contrib)
      eta = sqrt(v)
    }
    b = 1/2 * sum(((effect-mn)/eta)^2)
    c = sum(log(eta))
    Bij <- matrix(rep(0,number*nsteps),nrow=number,ncol=nsteps)
    bi = -si * qnorm(steps[1])
    Bij[,1] = 1-pnorm((bi-mn)/eta)
    if(nsteps > 2){
      for(j in 2:(length(steps)-1)) {
        bi = -si * qnorm(steps[j])
        bilast = -si * qnorm(steps[j-1])
        Bij[,j] = pnorm((bilast-mn)/eta) - pnorm((bi-mn)/eta)
      }
    }
    bilast = -si * qnorm(steps[length(steps)-1])
    Bij[,length(steps)] = pnorm((bilast-mn)/eta)

    swbij = 0
    for(j in 1:length(steps)) swbij = swbij + w[j]*Bij[,j]
    d = sum(log(swbij))

    return(-a + b + c + d)
  }

  model.print <- function(results){
    digits <- 4
    cat("\n")
    cat("Unadjusted Model (k = ", length(effect), "):", sep="")
    cat("\n\n")
    if(fe == FALSE){
      cat("tau^2 (estimated amount of total heterogeneity): ", formatC(round(results[[1]]$par[1], digits = digits), digits = digits, format = "f"), " (SE = ", formatC(round(sqrt(diag(solve(results[[1]]$hessian)))[1], digits = digits),digits = digits, format = "f"), ")", sep="")
      cat("\n")
      cat("tau (square root of estimated tau^2 value): ", formatC(round(sqrt(results[[1]]$par[1]), digits = digits),digits = digits, format = "f"))
      cat("\n\n")
    }
    cat("Model Results:")
    cat("\n\n")
    if(fe == FALSE){
      unadj_est <- cbind(results[[1]]$par[2:(npred+2)])
      unadj_se <- cbind(sqrt(diag(solve(results[[1]]$hessian)))[2:(npred+2)])
      z_stat <- unadj_est/unadj_se
      p_val <- (2*pnorm(-abs(z_stat)))
      ci.lb <- unadj_est - qnorm(0.975) * unadj_se
      ci.ub <- unadj_est + qnorm(0.975) * unadj_se
    }
    if(fe == TRUE){
      unadj_est <- cbind(results[[1]]$par[1:(npred+1)])
      unadj_se <- cbind(sqrt(diag(solve(results[[1]]$hessian)))[1:(npred+1)])
      z_stat <- unadj_est/unadj_se
      p_val <- (2*pnorm(-abs(z_stat)))
      ci.lb <- unadj_est - qnorm(0.975) * unadj_se
      ci.ub <- unadj_est + qnorm(0.975) * unadj_se
    }
    res.table <- data.frame(matrix(c(unadj_est, unadj_se, z_stat, p_val, ci.lb, ci.ub), nrow=(npred+1), byrow=F),stringsAsFactors=FALSE)
    rowlabels <- rep(0, (npred+1))
    rowlabels[1] <- "Intercept"
    if(npred > 0){
      for(i in 2:(npred+1)){
        rowlabels[i] <- paste(c(colnames(XX)[i]))
      }
    }
    row.names(res.table) <- c(rowlabels)
    colnames(res.table) <- c("estimate","std.error","z-stat","p-val","ci.lb","ci.ub")
    res.table[,4] <- format.pval(res.table[,4])
    res.table[,c(1,2,3,5,6)] <- format(res.table[,c(1,2,3,5,6)], digits=4)
    print(res.table)

    cat("\n")
    cat("Adjusted Model (k = ", length(effect), "):", sep="")
    cat("\n\n")
    if(fe == FALSE){
      if(is.null(weights)){
        cat("tau^2 (estimated amount of total heterogeneity): ", formatC(round(results[[2]]$par[1], digits = digits),digits = digits, format = "f"), " (SE = ", formatC(round(sqrt(diag(solve(results[[2]]$hessian)))[1], digits = digits),digits = digits, format = "f"), ")", sep="")
      }
      if(is.null(weights) == FALSE){
        cat("tau^2 (estimated amount of total heterogeneity): ", formatC(round(results[[2]]$par[1], digits = digits),digits = digits, format = "f"), " (SE = ", "N/E", ")", sep="")
      }
      cat("\n")
      cat("tau (square root of estimated tau^2 value): ", formatC(round(sqrt(results[[2]]$par[1]), digits = digits),digits = digits, format = "f"))
      cat("\n\n")
    }
    cat("Model Results:")
    cat("\n\n")

    if(fe == FALSE){
      if(is.null(weights)){
        adj_int_est <- cbind(results[[2]]$par[2:( (nsteps - 1) + (npred+2) )])
        adj_int_se <- cbind(sqrt(diag(solve(results[[2]]$hessian)))[2:( (nsteps - 1) + (npred+2) )])
      }
      else{
        adj_int_est <- cbind(c(results[[2]]$par[2:( (npred+2) )], weights[2:length(weights)]))
        adj_int_se <- cbind(rep("N/E", length(results[[2]]$par[2:length(results[[2]]$par)])))
      }
    }
    if(fe == TRUE){
      if(is.null(weights)){
        adj_int_est <- cbind(results[[2]]$par[1:( (nsteps - 1) + (npred+1) )])
        adj_int_se <- cbind(sqrt(diag(solve(results[[2]]$hessian)))[1:( (nsteps - 1) + (npred+1) )])
      }
      else{
        adj_int_est <- cbind(c(results[[2]]$par[1:( (npred+1) )], weights[2:length(weights)]))
        adj_int_se <- cbind(rep("N/E", length(results[[2]]$par[1:length(results[[2]]$par)])))
      }
    }

    if(is.null(weights)){
      z_stat_int <- adj_int_est/adj_int_se
      p_val_int <- (2*pnorm(-abs(z_stat_int)))
      ci.lb_int <- adj_int_est - qnorm(0.975) * adj_int_se
      ci.ub_int <- adj_int_est + qnorm(0.975) * adj_int_se
    }
    else{
      z_stat_int <- cbind(rep("N/E", length(results[[2]]$par[1:length(results[[2]]$par)])))
      p_val_int <- cbind(rep("N/E", length(results[[2]]$par[1:length(results[[2]]$par)])))
      ci.lb_int <- cbind(rep("N/E", length(results[[2]]$par[1:length(results[[2]]$par)])))
      ci.ub_int <- cbind(rep("N/E", length(results[[2]]$par[1:length(results[[2]]$par)])))
    }
    res.table <- data.frame(matrix(c(adj_int_est, adj_int_se, z_stat_int, p_val_int, ci.lb_int, ci.ub_int), nrow=(npred+1+(nsteps-1)), byrow=F),stringsAsFactors=FALSE)

    rowlabels1 <- rep(0, (npred+1))
    rowlabels1[1] <- "Intercept"
    if(npred > 0){
      for(i in 2:length(rowlabels1)){
        rowlabels1[i] <- paste(c(colnames(XX)[i]))
      }
    }
    rowlabels2 <- rep(0, (nsteps-1))
    for(i in 1:(length(rowlabels2))){
      rowlabels2[i] <- paste(c(steps[i], "< p <", steps[i + 1]), collapse=" ")
    }
    row.names(res.table) <- c(rowlabels1,rowlabels2)
    colnames(res.table) <- c("estimate","std.error","z-stat","p-val","ci.lb","ci.ub")
    if(is.null(weights)){
      res.table[,"p-val"] <- format.pval(res.table[,"p-val"])
    }
    res.table[,c(1,2,3,5,6)] <- format(res.table[,c(1,2,3,5,6)], digits=4)
    print(res.table)

    if(is.null(weights)){
      cat("\n")
      cat("Likelihood Ratio Test:")
      cat("\n")
      df <- length(results[[2]]$par) - length(results[[1]]$par)
      lrchisq <- 2*(results[[1]]$value - results[[2]]$value)
      pvalue <- 1-pchisq(lrchisq,df)
      cat("X^2(df = ", df, ") = ", lrchisq, ", p-val = ", format.pval(pvalue), sep="")
    }

    if(table){
      cat("\n")
      cat("\n")
      cat("Number of Effect Sizes per Interval:")
      cat("\n")
      cat("\n")
      format(print(sampletable(p, pvalues, steps)))
    }
  }

  if(length(effect)!=length(v)){
    stop('Your vector of effect sizes and your vector of sampling variances are not the same length. Please check your data.')
  }

  if(identical(effect,v)){
    stop('Your vector of effect sizes is exactly the same as your vector of sampling variances. Please check your data.')
  }

  if(min(v) < 0){
    stop('Sampling variances cannot be negative. Please check your data.')
  }

  si <- sqrt(v)
  p <- 1-pnorm(effect/sqrt(v))
  if(max(steps)!=1){
    steps <- c(steps,1)
  }
  if(max(steps) > 1){
    stop('p-value cutpoints cannot be greater than 1.')
  }
  if(min(steps) < 0){
    stop('p-value cutpoints cannot be negative.')
  }
  if(length(unique(steps)) != length(steps)){
    stop('Two or more p-value cutpoints are identical.')
  }

  if(is.null(weights)){
    steps <- sort(steps)
  }
  if(is.null(weights) == FALSE){
    if(min(weights) < 0){
      stop('Weights for p-value intervals cannot be negative.')
    }
    if(length(weights)!=length(steps)){
      stop('The number of weights does not match the number of p-value intervals created.')
    }
    new <- cbind(steps, weights)
    steps <- new[order(steps),1]
    weights <- new[order(steps),2]
  }

  number <- length(effect)
  nsteps <- length(steps)


  wt <- rep(1,number)
  for(i in 1:number) {
    if(p[i] <= steps[1]) wt[i] = 1
    for(j in 2:nsteps) {
      if (steps[j-1] < p[i] && p[i] <= steps[j]) wt[i] = j
    }
    if(  p[i] > steps[nsteps-1] ) wt[i] = nsteps
  }

  intervaltally <- function(p, steps) {
    p1 <- cut(p, breaks=c(-Inf,steps), labels=steps)
    return(p1) }
  pvalues <- as.numeric(table(intervaltally(p, steps)))

  sampletable <- function(p, pvalues, steps){
    nsteps <- length(steps)
    results <- matrix(nrow=length(pvalues),ncol=1)
    results[,1] <- pvalues
    rowlabels <- c(0, length(results[,1]))
    rowlabels[1] <- paste(c("p-values <", steps[1]), collapse="")
    for(i in 2:nsteps){
      rowlabels[i] <- paste(c(steps[i - 1], "< p-values <", steps[i]), collapse=" ")
    }
    resultsb <- data.frame(results, row.names=c(rowlabels))
    colnames(resultsb) <- c("Frequency")
    return(resultsb)
  }

  if(sum(table(intervaltally(p,steps)) == 0) >= 1){
    warning('At least one of the p-value intervals contains no effect sizes, leading to estimation problems. Consider re-specifying the cutpoints.')
  }

  if(sum( table(intervaltally(p,steps)) > 0 & table(intervaltally(p, steps)) <= 3) >= 1){
    warning('At least one of the p-value intervals contains three or fewer effect sizes, which may lead to estimation problems. Consider re-specifying the cutpoints.')
  }

  if(is.null(mods)){

    if(fe == FALSE){

      XX <- cbind(rep(1,length(effect)))
      npred <- 0
      pars <- c(mean(v)/4, mean(effect), rep(1,(nsteps-1)))

      output_unadj <- optim(par=pars[1:2],fn=neglike_unadj,lower=c(0,-Inf),method="L-BFGS-B",hessian=TRUE)

      output_adj <- optim(par=pars,fn=neglike_adj,lower=c(0,-Inf, rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)

      results <- list(output_unadj,output_adj)

      model.print(results)

      if(is.null(weights)){
        results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                         adj_est=c(output_adj$par),adj_se=c(sqrt(diag(solve(output_adj$hessian)))),
                         z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         z_adj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         p_adj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.lb_adj=c(output_adj$par - qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))),
                         ci.ub_adj=c(output_adj$par + qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))))
      }
      if(is.null(weights) == FALSE){
        results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                         adj_est=c(output_adj$par),
                         z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))))
      }

    }

    if(fe == TRUE){

      XX <- cbind(rep(1,length(effect)))
      npred <- 0
      pars <- c(mean(effect), rep(1,(nsteps-1)))

      output_unadj <- optim(par=pars[1],fn=neglike_unadj,lower=c(-Inf),method="L-BFGS-B",hessian=TRUE)

      output_adj <- optim(par=pars,fn=neglike_adj,lower=c(-Inf, rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)

      results <- list(output_unadj,output_adj)

      model.print(results)

      if(is.null(weights)){
        results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                         adj_est=c(output_adj$par),adj_se=c(sqrt(diag(solve(output_adj$hessian)))),
                         z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         z_adj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         p_adj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.lb_adj=c(output_adj$par - qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))),
                         ci.ub_adj=c(output_adj$par + qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))))
      }
      if(is.null(weights) == FALSE){
        results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                         adj_est=c(output_adj$par),
                         z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))))
      }

    }


  }

  else{

    if(typeof(mods)=="language"){

      if(fe == FALSE){

        XX <- model.matrix(mods)
        npred <- dim(XX)[2]-1
        pars <- c(mean(v)/4, mean(effect), rep(0, npred), rep(1, (nsteps - 1)))

        output_unadj <- optim(par=pars[1:(npred+2)],fn=neglike_unadj,lower=c(0,rep(-Inf, (npred+1))),method="L-BFGS-B",hessian=TRUE)

        output_adj <- optim(par=pars,fn=neglike_adj,lower=c(0,rep(-Inf, (npred+1)),rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)

        results <- list(output_unadj,output_adj)

        model.print(results)

        if(is.null(weights)){
          results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                           adj_est=c(output_adj$par),adj_se=c(sqrt(diag(solve(output_adj$hessian)))),
                           z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           z_adj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           p_adj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.lb_adj=c(output_adj$par - qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))),
                           ci.ub_adj=c(output_adj$par + qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))))
        }
        if(is.null(weights) == FALSE){
          results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                           adj_est=c(output_adj$par),
                           z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))))
        }

      }

      if(fe == TRUE){

        XX <- model.matrix(mods)
        npred <- dim(XX)[2]-1
        pars <- c(mean(effect), rep(0, npred), rep(1, (nsteps - 1)))

        output_unadj <- optim(par=pars[1:(npred+1)],fn=neglike_unadj,lower=c(rep(-Inf, (npred+1))),method="L-BFGS-B",hessian=TRUE)

        output_adj <- optim(par=pars,fn=neglike_adj,lower=c(rep(-Inf, (npred+1)),rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)

        results <- list(output_unadj,output_adj)

        model.print(results)

        if(is.null(weights)){
          results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                           adj_est=c(output_adj$par),adj_se=c(sqrt(diag(solve(output_adj$hessian)))),
                           z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           z_adj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           p_adj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.lb_adj=c(output_adj$par - qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))),
                           ci.ub_adj=c(output_adj$par + qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))))
        }
        if(is.null(weights) == FALSE){
          results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                           adj_est=c(output_adj$par),
                           z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))))
        }

      }

    }

    else{

      stop('Moderators must be entered as "mods= ~ x1 + x2"')

    }

  }

  invisible(results2)

}
