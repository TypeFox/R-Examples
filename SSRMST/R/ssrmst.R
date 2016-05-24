#'  @name SSRMST-package
#'  @aliases  SSRMST-package
#'  @docType  package
#'  @title Sample Size Calculation using Restricted Mean Survival Time
#'  @description
#'  The differences in restricted mean survival times (RMST), a clinically interpretable model-free measure, can be one of the alternatives to the hazard ratio.
#'  The package calculates the study sample size and power in designing clinical trials using the differences in RMSTs.
#'  Two types of one-sided tests, non-inferiority and superiority tests, are prepared.
#'  @author Miki Horiguchi
#'  @details Please check the vignette for details: \code{browseVignettes(package = "SSRMST")}
#'  @references
#'  Uno H, Claggett B, Tian L, Inoue E, Gallo P, Miyata T, Schrag D, Takeuchi M, Uyama Y, Zhao L,
#'  Skali H, Solomon S, Jacobus S, Hughes M, Packer M, Wei LJ. Moving beyond the hazard ratio in
#'  quantifying the between-group difference in survival analysis. Journal of clinical Oncology 2014,
#'  32, 2380-2385.
#'
#'  Uno H, Wittes J, Fu H, Solomon SD, Claggett B, Tian L, Cai T, Pfeffer MA, Evans SR, Wei LJ.
#'  Alternatives to Hazard Ratios for Comparing the Efficacy or Safety of Therapies in non-inferiority Studies.
#'  Annals of Internal Medicine 2015, 163, 127-134.
#'  @seealso
#'  survival
#'  survRM2
#'  @examples
#'  #---Example data
#'
#'  ac_rate   = 15
#'  ac_period = 35
#'  tot_time  = 510
#'  tau       = 500
#'  scale0    = 8500
#'  scale1    = 8500
#'  margin    = 18
#'
#'  a = ssrmst(ac_rate, ac_period, tot_time, tau, scale0, scale1, margin=margin, ntest=20)
#'  print(a)
NULL

#' Sample Size Calculation using Restricted Mean Survival Time
#'
#' The package calculates the study sample size and power in designing clinical trials using the differences in restricted mean survival times (RMST).
#' Two types of one-sided tests, non-inferiority and superiority tests, are prepared.
#' Under certain conditions, 2,000 sets of realizations in default are generated for calculating confidence intervals of RMST differences.
#' Then the power is calculated, i.e., the chance that the lower bound of 2,000 confidence intervals of RMST differences falls above a margin.
#'
#'
#' @usage
#' ssrmst(ac_rate, ac_period, tot_time, tau, scale0, scale1, shape = 1, margin = 0,
#'        allocation1 = 0.5, one_sided_alpha = 0.025, seed=NULL, ntest=2000)
#'
#' @param ac_rate         Accrual rate: the number of patients per unit time.
#' @param ac_period       Accrual period: the time point at last accrual.
#' @param tot_time        Total study time: the time point at last follow-up.
#' @param tau             Truncation time point to calculate RMSTs.
#' @param scale0,scale1   Scale parameters for the Weibull distribution in both the control (arm0) and the treatment (arm1). Note that the value of the scale parameter in the treatment (arm1) needs to be larger than or equal to that in the control (arm0), because the difference of the RMSTs (arm1 minus arm0) is of interest.
#' @param shape           Shape parameter for the Weibull distribution in both arms. As we assume the proportional hazards, the shape parameters for the Weibull distribution in both arms are the same. When default (\code{shape = 1}), the Weibull distribution reduces to the Exponential distribution.
#' @param margin          Non-inferiority margin: a clinically acceptable difference in RMSTs. A value of minus \code{margin} is used to evaluate the power. When default (\code{margin = 0}), a superiority test is selected.
#' @param allocation1     Proportion of patients allocated to the treatment (arm1). Default value is 0.5.
#' @param one_sided_alpha Nominal type I error level as one-sided. When default (\code{one_sided_alpha = 0.025}), 0.95 confidence intervals of the difference in RMSTs are estimated to calculate the power.
#' @param seed            Random seed used for the sampling. Default is \code{NULL}.
#' @param ntest           Number of simulations. When default (\code{ntest = 2000}), 2,000 sets of realizations are generated for calculating confidence intervals of RMST differences.
#'
#' @return A list with components:
#' @return \item{result}{Total study population and expected number of events.}
#' @return \item{power}{Chance that the lower bound of 2,000 confidence intervals of differences in RMSTs falls above a value of minus margin in a non-inferiority test (or above 0 in a superiority test).}
#' @return \item{accrual rate}{Accrual rate used in the analyses.}
#' @return \item{accrual period}{Accrual period used in the analyses.}
#' @return \item{total study time}{Total study time used in the analyses.}
#' @return \item{margin}{Margin used in the analyses.}
#' @return \item{tau}{Truncation time point used in the analyses.}
#' @return \item{note}{Note regarding the truncation time, tau.}
#' @details For more details, please refer to the vignette: \code{browseVignettes(package = "SSRMST")}
#' @seealso
#'  survival
#'  survRM2
#' @references
#'  Uno H, Wittes J, Fu H, Solomon SD, Claggett B, Tian L, Cai T, Pfeffer MA, Evans SR, Wei LJ.
#'  Alternatives to Hazard Ratios for Comparing the Efficacy or Safety of Therapies in non-inferiority Studies.
#'  Annals of Internal Medicine 2015, 163, 127-134.
#' @name   ssrmst
#' @aliases ssrmst
#' @import survival survRM2
#' @importFrom stats rweibull runif
#' @examples
#'  #---Example data
#'
#'  #--Non-inferiority test
#'  ac_rate   = 15
#'  ac_period = 35
#'  tot_time  = 510
#'  tau       = 500
#'  scale0    = 8500
#'  scale1    = 8500
#'  margin    = 18
#'
#'  a = ssrmst(ac_rate, ac_period, tot_time, tau, scale0, scale1, margin=margin, ntest=20)
#'  print(a)
#'
#'
#'  #--Superiority test
#'  ac_rate   = 15
#'  ac_period = 35
#'  tot_time  = 510
#'  tau       = 500
#'  scale0    = 4000
#'  scale1    = 8500
#'  b = ssrmst(ac_rate, ac_period, tot_time, tau, scale0, scale1, ntest=20)
#'  print(b)
#' @export
ssrmst <-
  function(ac_rate, ac_period, tot_time, tau, scale0, scale1, shape=1,
           margin=0, allocation1=0.5, one_sided_alpha=0.025, seed=NULL, ntest=2000){

    ###--- initial check ---
    if (tau >= tot_time){
      stop(paste("The truncation time, tau, needs to be shorter than total study time, tot_time."))
    }
    if (allocation1<=0 | allocation1>=1){
      stop(paste("The proportion of patients allocated to the treatment (arm1), allocation1, needs to be between 0 and 1."))
    }
    if (ac_period > tot_time){
      n0 = round(ac_rate*tot_time*(1-allocation1))
      n1 = round(ac_rate*tot_time*allocation1)
    }
    if (ac_period <= tot_time){
      n0 = round(ac_rate*ac_period*(1-allocation1))
      n1 = round(ac_rate*ac_period*allocation1)
    }
    if (margin<0) {
      stop(paste("The margin needs to be larger than or equal to 0."))
    }
    if (shape<=0) {
      stop(paste("The value of the shape parameter for the Weibull distribution in both arms needs to be larger than 0."))
    }
    if (scale0 > scale1){
      stop(paste("The value of the scale parameter for the Weibull distribution in the treatment (arm1) needs to be larger than or equal to that in the control (arm0)."))
    }


    ###--- test (main part) ---
    answer     = NULL
    check      = NULL
    event_arm0 = NULL
    event_arm1 = NULL

    if (!is.null(seed)){
      set.seed(seed)
    }
    for (w in 1:ntest){

      ##-- data frame --
      E             = rweibull(n0, shape, scale0)
      C             = tot_time - runif(n0, 0, ac_period)
      time          = pmin(E,C)
      status        = as.numeric(E<=C)
      arm           = rep(0,n0)
      data0         = data.frame(time, status, arm)
      ind0          = data0$status==1
      event_arm0[w] = sum(data0$time[ind0]<=tot_time)

      E             = rweibull(n1, shape, scale1)
      C             = tot_time - runif(n1, 0, ac_period)
      time          = pmin(E,C)
      status        = as.numeric(E<=C)
      arm           = rep(1,n1)
      data1         = data.frame(time, status, arm)
      ind1          = data1$status==1
      event_arm1[w] = sum(data1$time[ind1]<=tot_time)

      data   = rbind(data0, data1)
      data   = data[data$time>0, ]

      ##-- tau check --
      idx = data$arm==0; tt = data$time[idx]; event = data$status[idx]; tau0 = max(tt[event==1]); tau0max = max(tt)
      idx = data$arm==1; tt = data$time[idx]; event = data$status[idx]; tau1 = max(tt[event==1]); tau1max = max(tt)

      tau_default = min(tau0,    tau1)
      tau_max     = min(tau0max, tau1max)

      if (tau >= tau_default & tau <= tau_max){
        check[w] = 1
      } else{
        check[w] = 0
      }
      if (tau >= tau_max){
        stop(paste("The truncation time, tau, needs to be shorter than or equal \n to the minimum of the largest observed time on each of the two arms: ", round(tau_max, digits=3)))
      }

      ##-- RMST calculation --
      ans = rmst2(data$time, data$status, data$arm, tau=tau, alpha=one_sided_alpha*2)

      lower = ans$unadjusted.result[1,2]
      if (lower > -margin){
        answer[w] = 1
      } else {
        answer[w] = 0
      }
    }


    ###--- expected number of events ---
    n0_event = round(sum(event_arm0)/ntest)
    n1_event = round(sum(event_arm1)/ntest)


    ###--- power ---
    power = sum(answer)/ntest


    ###--- final check ---
    if (sum(check) >= 1){
      NOTE = paste("The truncation time: tau =", tau, " was specified, but there are no observed events after tau=,", tau, "on either or both groups. Make sure that the size of riskset at tau=,", tau, "is large enough in each arm.")
    } else{
      NOTE = paste("The truncation time: tau =", tau, " was specified.")
    }


    ###--- output ---
    out = matrix(0,2,3)

    out[1,1] = n0+n1
    out[1,2] = n0
    out[1,3] = n1
    out[2,1] = n0_event + n1_event
    out[2,2] = n0_event
    out[2,3] = n1_event

    rownames(out) = c("Sample size", "Expected number of events")
    colnames(out) = c("Total", "arm0", "arm1")

    Z = list()

    Z$result      = out
    Z$power       = power
    Z$ac_rate     = ac_rate
    Z$ac_period   = ac_period
    Z$tot_time    = tot_time
    Z$margin      = margin
    Z$tau         = tau
    Z$note        = NOTE

    class(Z) = "ssrmst"

    Z
  }
NULL
#' @name print.ssrmst
#' @aliases print.ssrmst
#' @title print.ssrmst
#' @description S3 method for class 'ssrmst'
#' @param x Object to be printed.
#' @param ... Further arguments ignored in this function.
#' @export
print.ssrmst <- function(x, ...){

  ###--- superiority test ---
  if (x$margin==0){
    cat ("Superiority test \n")
    cat ("\n")
    print(x$result)
    cat ("\n")
    names(x$power) = "Power"
    print(x$power)
  }

  ###--- non-inferiority test ---
  if (x$margin>0){
    cat ("Non-inferiority test \n")
    cat ("\n")
    print(x$result)
    cat ("\n")
    names(x$power) = "Power"
    print(x$power)
  }

  invisible(x)
}


