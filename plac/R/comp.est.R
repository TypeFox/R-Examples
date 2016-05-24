#' Calculate the PLAC estimator when a time-dependent indicator presents
#'
#' Both a conditional approach Cox model and a pairwise likelihood augmented
#' estimator are fitted and the corresponding results are returned in a list.
#' @param ltrc.formula a formula of of the form \code{Surv(A, Y, D) ~ Z}, where
#'   \code{Z} only include the time-invariate covariates.
#' @param ltrc.data a data.frame of the LTRC dataset including the responses,
#'   time-invariate covariates and the jump times for the time-depnencent
#'   covariate.
#' @param id.var a name of the subject id in \code{data}.
#' @param td.var  a name of the time-dependent covariate in the output.
#' @param td.type the type of the time-dependent covariate. Either one of
#'   \code{c("none", "independent", "post-trunc", "pre-post-trunc")}. See
#'   Details.
#' @param t.jump a name of the jump time variable in \code{data}.
#' @param init.val a list of the initial values of the coefficients and the
#'   baseline hazard function for the PLAC estimator.
#' @param max.iter the maximal number of iteration for the PLAC estimator
#' @param print.result logical, if a brief summary of the regression coefficient
#'   estiamtes should be printed out.
#' @param ... other arguments
#' @useDynLib plac
#' @importFrom stats as.formula coef model.frame model.matrix model.response
#'   model.matrix pnorm
#' @importFrom Rcpp sourceCpp
#' @importFrom survival Surv tmerge coxph
#'
#' @details \code{ltrc.formula} should have the same form as used in \code{coxph()}; e.g., \code{Surv(A, Y, D) ~ Z1 + Z2}.
#'   where \code{A} is the truncation time (\code{tstart}), \code{Y} is the
#'   survival time (\code{tstop}) and \code{D} is the status indicator
#'   (\code{event}). \code{td.type} is used to determine which \code{C++}
#'   function will be invoked: either \code{PLAC_TI} (if \code{td.type =
#'   "none"}), \code{PLAC_TD} (if \code{td.type = "independent"}) or
#'   \code{PLAC_TDR}) (if \code{td.type \%in\% c("post-trunc",
#'   "pre-post-trunc")}). For \code{td.type = "post-trunc"}, the pre-truncation
#'   values for the time-dependent covariate will be set to be zero for all
#'   subjects.
#'
#' @return a list of model fitting results for both conditional approach and the
#'   PLAC estimators.
#'   \describe{
#'   \item{\code{Event.Time}}{Ordered distinct observed event times}
#'   \item{\code{b}}{Regression coefficients estiamtes}
#'   \item{\code{se.b}}{Model-based SEs of the regression coefficients
#'   estiamtes}
#'   \item{\code{H0}}{Estimated cumulative baseline hazard function}
#'   \item{\code{se.H0}}{Model-based SEs of the estimated cumulative baseline
#'   hazard function}
#'   \item{\code{sandwich}}{The sandwich estimator for (beta,
#'   lambda)}
#'   \item{\code{k}}{The number of iteration for used for the PLAC
#'   estimator}
#'   \item{\code{summ}}{A brief summary of the covariates effects}
#'   }
#' @references Wu, F. Kim, S. and Li, Y. "A Pairwise Likelihood Augmented
#'   Estimator for Left-Truncated Data with Time-Dependent Covariates."
#'   (\emph{in preparation})
#' @references Wu, F., Kim, S., Qin, J., Saran, R. and Li, Y. (2015) "A
#'   Pairwise-Likelihood Augmented Estimator for the Cox Model Under
#'   Left-Truncation." (Submitted to \emph{Journal of American Statistical
#'   Association}.) \url{http://biostats.bepress.com/umichbiostat/paper118/}
#' @examples
#' # When only time-invariant covariates are involved
#' dat1 = sim.ltrc(n = 50)$dat
#' PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z1 + Z2,
#'      ltrc.data = dat1, td.type = "none")
#' # When there is a time-dependent covariate that is independent of the truncation time
#' dat2 = sim.ltrc(n = 50, time.dep = TRUE,
#'                distr.A = "binomial", p.A = 0.8, Cmax = 5)$dat
#' PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z,
#'      ltrc.data = dat2, td.type = "independent",
#'      td.var = "Zv", t.jump = "zeta")
#' # When there is a time-dependent covariate that depends on the truncation time
#' dat3 = sim.ltrc(n = 50, time.dep = TRUE, Zv.depA = TRUE, Cmax = 5)$dat
#' PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z,
#'      ltrc.data = dat3, td.type = "post-trunc",
#'      td.var = "Zv", t.jump = "zeta")
#'
#' @export
PLAC = function(ltrc.formula, ltrc.data, id.var = "ID",
                td.var = NULL, td.type = "none", t.jump = NULL,
                init.val = NULL, max.iter = 100, print.result = TRUE, ...){
  if( !inherits(ltrc.formula, "formula") ) stop("A formula of the form 'Surv (A, Y, D) ~ Z' is required!")
  # grep the model.frame for later use
  mf = model.frame(formula = ltrc.formula, data = ltrc.data)
  # prepare the response and time-invariant covariate matrix
  X = as.matrix(model.response(mf))
  resp = gsub("Surv|\\(|\\)| ","",unlist(strsplit(names(mf)[1],",")))
  colnames(X) = resp
  n = nrow(X)
  W = unique(X[X[,3] == 1, 2])
  m = length(W)
  # at-risk processes
  Ind1 = SgInd(X, W)
  # truncation processes
  Ind2 = PwInd(X, W)
  # number of events at each w_k
  Dn = as.numeric(rle(X[X[,3] == 1, 2])$lengths)
  ZF = as.matrix(model.matrix(attr(mf, "terms"), data = mf)[ , -1], nrow = n)
  colnames(ZF) = names(mf)[-1]
  if( td.type == "none"){
    p = ncol(ZF)
    Z = t(ZF)
    cox.LTRC = survival::coxph(formula = ltrc.formula, data = ltrc.data, method="breslow", model = TRUE)
  }else{
    # for "post-trunc", all subjects have pre-trunc Zv = 0.
    if( td.type == "post-trunc" ){
      assign(td.var, rep(0, n))
      ZF0 = NULL
      eval(parse(text = paste0("ZF0 = cbind(ZF, ", td.var, ")")))
      ZFt = t(ZF0)
    }else{
      ZFt = t(ZF)
    }
    # the jump times of the time-dependent indicator (zeta)
    ZV = ltrc.data[[t.jump]]
    # need counting process expansion of the data
    eval(parse(text = paste0("data.count = tmerge(ltrc.data, ltrc.data, id = ", id.var,
                             ", death = event(", resp[2], ", ", resp[3], "), ",
                             td.var, " = tdc(", t.jump, "))")))
    data.count$id = NULL
    eval(parse(text = paste0("data.count = subset(data.count, !(",
                             t.jump, " <= ", resp[1],
                             "& tstart == 0))")))
    eval(parse(text = paste0("data.count$tstart[data.count$",
                             resp[1], "> data.count$tstart] = data.count$",
                             resp[1], "[data.count$",
                             resp[1], "> data.count$tstart]")))
    # covariate values at the observed survival times
    ZV_ = NULL
    eval(parse(text = paste0("ZV_ = subset(data.count, select = ", td.var, ", tstop == ", resp[2],")[[1]]")))
    ZFV_ = rbind(t(ZF), ZV_)
    p = nrow(ZFV_)
    IndZ = TvInd(ZV, W)
    Za = array(c(rep(ZF, each = m), IndZ), c(m, n, p))
    Z = matrix(aperm(Za, c(3, 1, 2)), m * p, n)
    cox.LTRC = survival::coxph(as.formula(paste0("Surv(tstart, tstop, death) ~ ",
                                                 paste(c(colnames(ZF), td.var),
                                                       collapse = "+"))),
                               data = data.count, model = TRUE)

  }
  b.cox = unname(coef(cox.LTRC))
  h.cox = survival::basehaz(cox.LTRC, centered = FALSE)$hazard
  # get h from H
  if(h.cox[1]!=0){
    H.cox = c(0, unique(h.cox))
    h.cox = diff(H.cox)
  }else{
    # if the first obs is censored..
    H.cox = unique(h.cox)
    h.cox = diff(H.cox)
  }
  newdata = model.matrix(cox.LTRC)[FALSE,]
  newdata = data.frame(rbind(newdata, rep(0, ncol(newdata))))
  summ.cox = summary(survival::survfit(cox.LTRC, newdata = newdata))
  se.H0.cox = c(0, summ.cox$std.err/summ.cox$surv)
  # set the initial values for b (coefficients) and h (baseline hazard)
  if( is.null(init.val) ){
    b_0 = b.cox
    h_0 = h.cox
  }else{
    b_0 = init.val$b_0
    h_0 = init.val$h_0
  }
  # Call the proper C++ wrapper function
  if( td.type == "none" ){
    cat("Calling PLAC_TI()...\n")
    plac.fit = PLAC_TI(Z, X, W, Ind1, Ind2, Dn, b_0, h_0, max.iter)
  }else if( td.type == "independent" ){
    cat("Calling PLAC_TD()...\n")
    plac.fit = PLAC_TD(Z, ZFV_, X, W, Ind1, Ind2, Dn, b_0, h_0, max.iter)
  }else if( td.type %in% c("post-trunc", "pre-post-trunc") ){
    cat("Calling PLAC_TDR()...\n")
    plac.fit = PLAC_TDR(ZFt, ZFV_, Z, X, W, Ind1, Ind2, Dn, b_0, h_0, max.iter)
  }
  # summarizing fitting results
  b = cbind(Cox = b.cox, PLAC = plac.fit$b.hat)
  se.b = cbind(Cox = sqrt(diag(cox.LTRC$var)), PLAC = plac.fit$se.b.hat)
  if( td.type == "none" ){
    rownames(b) = colnames(ZF)
  }else{
    rownames(b) = c(colnames(ZF), td.var)
  }
  rownames(se.b) = rownames(b)
  H0 = cbind(Cox = H.cox, PLAC = plac.fit$H0.hat)
  se.H0 = cbind(Cox = se.H0.cox, PLAC = plac.fit$se.H0.hat)

  rslt = cbind(est.Cox = b[ , 1], se.Cox = se.b[ , 1],
               p.Cox = pnorm(2 * abs(b[ , 1] / se.b[ , 1]), lower.tail = FALSE),
               est.PLAC = b[ , 2], se.PLAC = se.b[ , 2],
               p.PLAC = pnorm(2 * abs(b[ , 2] / se.b[ , 2]), lower.tail = FALSE))
  rslt = round(rslt, 3)
  if( print.result ){
    cat("Coefficient Estimates:\n")
    print(rslt)
  }

  return(invisible(list(Event.Time = summ.cox$time,
                        b = b, se.b = se.b,
                        H0 = H0, se.H0 = se.H0,
                        sandwich = plac.fit$swe,
                        k = plac.fit$iter,
                        summ = rslt)))
}

#' Calulate the Values of the cumulative Hazard at Fixed Times
#'
#' @param est an object of the class \code{plac.fit}.
#' @param t.eval time points at which the Lambda(t) is evaluated (for both conditional apporach and the PLAC estimator).
#' @importFrom stats stepfun
#' @return a list containing the estiamtes and SEs of Lambda(t) for both conditional apporach and the PLAC estimator.
#' @examples
#' dat1 = sim.ltrc(n = 100)$dat
#' est = PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z1 + Z2,
#'      ltrc.data = dat1, td.type = "none")
#' H = cum.haz(est, t.eval = seq(0.1, 0.9, 0.1))
#' H$L
#' H$se.L
#' @export
cum.haz = function(est, t.eval=c(0.25, 0.75)){

  # evaluation of Lambda_hat
  tim = est$Event.Time

  H0.cox = stepfun(tim, est$H0[,1])
  H0.cox.eval = H0.cox(t.eval)

  H0.plac = stepfun(tim,est$H0[,2])
  H0.plac.eval = H0.plac(t.eval)

  L = rbind(Cox=H0.cox.eval, PLAC=H0.plac.eval)

  se.H0.cox = stepfun(tim,est$se.H0[,1])
  se.H0.cox.eval = se.H0.cox(t.eval)

  se.H0.plac = stepfun(tim,est$se.H0[,2])
  se.H0.plac.eval = se.H0.plac(t.eval)

  se.L = rbind(Cox=se.H0.cox.eval, PLAC=se.H0.plac.eval)

  colnames(L) = t.eval
  colnames(se.L) = t.eval

  return(list(TrueL = t.eval,
              L=L, se.L=se.L,
              L.fun=list(Cox=H0.cox, PLAC=H0.plac),
              se.L.fun=list(Cox=se.H0.cox, PLAC=se.H0.plac)))

}
