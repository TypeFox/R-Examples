#' Fit Generalized Competing Event Model Based on Proportional Hazards Regression
#'
#' @description Fit a generalized competing event model by using Cox proportational hazards regression model
#' with \code{coxph} function in \code{survival} package.
#' @param formula1 a formula object for event(s) of interest, with a survival response returned by \code{Surv}
#' function on the left, and the covariate terms on the right.
#' @param formula2 a formula object for competing event(s), with a survival response returned by \code{Surv}
#' function on the left, and the covariate terms on the right.
#' @param formula3 a formula object for the composite set of all events, with a survival response returned by \code{Surv}
#' function on the left, and the covariate terms on the right.
#' @param surv1 a formula object for event(s) of interest, with a survival response returned by \code{Surv}
#' function on the left, and 1 on the right.
#' @param surv2 a formula object for competing event(s), with a survival response returned by \code{Surv}
#' function on the left, and 1 on the right.
#' @param data a data frame containing variables named in formula.
#' @param N the number of bootstrap replicates
#' @param M the number of bins for \eqn{\omega} or \eqn{\omega+} plots.
#' @param t survival time point for \eqn{\omega} or \eqn{\omega+} plots.
#' @details The \strong{gcerisk} package is designed to help investigators optimize risk-stratification methods for competing risks data, such as described in
#' Carmona R, Gulaya S, Murphy JD, Rose BS, Wu J, Noticewala S, McHale MT, Yashar CM, Vaida F, Mell LK. Validated competing event model for the stage I-II endometrial cancer population.
#' Int J Radiat Oncol Biol Phys. 2014;89:888-98. Standard risk models typically estimate the effects of one or more covariates on either
#' a single event of interest (such as overall mortality, or disease recurrence), or a composite set of events (e.g., disease-free survival, which combines events of interest with death from any cause).
#' This method is inefficient in stratifying patients who may be simultaneously at high risk for the event of interest but low risk for competing events, and who thus stand to gain the most from strategies to modulate the event of interest.
#' Compared to standard risk models, GCE models better stratify patients at higher (lower) risk for an event of interest and lower (higher) risk of competing events. GCE models focus on differentiating subjects based on
#' the ratio of the cumulative hazard (or cumulative hazard of the subdistribution) for the event of interest to the cumulative hazard (or cumulative hazard of the subdistribution) for all events (\eqn{\omega}),
#' and the ratio of the cumulative hazard (or cumulative hazard of the subdistribution) for the event of interest to the cumulative hazard (or cumulative hazard of the subdistribution) for competing events (\eqn{\omega+}).
#'
#' The \code{gcecox} function produces model estimates and confidence intervals from a generalized competing event model based on the Cox PH model for cause-specific hazards. The model assumes proportional hazards for the composite set of events.
#'
#' The function returns \eqn{\omega} and \eqn{\omega+} ratio estimates for the chosen covariates, with 95\% confidence intervals, and plots \eqn{\omega} and \eqn{\omega+} at time t within M ordered subsets of subjects as a function of increasing risk (based on the linear predictor, i.e. the inner product of a subject's data vector and the coefficient vector).
#' @import survival cmprsk ggplot2 stats
#' @examples
#' # sample data to test
#' data(Sample)
#' test <- Sample
#' rm(list=setdiff(ls(), "test"))
#' test <- transform(test, LRF_OR_DF_FLAG = as.numeric(test$LRFFLAG | test$DFFLAG))
#' test <- transform(test, LRF_OR_DF_MO = pmin(test$LRFMO, test$DFMO))
#' test <- transform(test, CMFLAG = as.numeric(test$OSFLAG & !test$LRFFLAG & !test$DFFLAG))
#' test <- transform(test, ACMFLAG = as.numeric(test$LRF_OR_DF_FLAG | test$CMFLAG))
#' test <- transform(test, ACM_MO = pmin(test$LRF_OR_DF_MO, test$OSMO))
#'
#' formula1 <- Surv(LRF_OR_DF_MO, LRF_OR_DF_FLAG) ~ age + gender + smoke20 +
#' etohheavy + higrade + BMI + black
#' formula2 <- Surv(OSMO, CMFLAG) ~ age + gender + smoke20 + etohheavy + higrade + BMI + black
#' formula3 <- Surv(ACM_MO, ACMFLAG) ~ age + gender + smoke20 + etohheavy + higrade + BMI + black
#' surv1 <- Surv(LRF_OR_DF_MO, LRF_OR_DF_FLAG) ~ 1
#' surv2 <- Surv(OSMO, CMFLAG) ~ 1
#' N <- 100
#' M <- 5
#' t <- 60
#'
#' fitgce.cox <- gcecox(formula1, formula2, formula3, surv1, surv2, test, N, M, t)

#' @author Hanjie Shen, Ruben Carmona, Loren Mell
#' @references
#' \itemize{
#' \item Carmona R, Gulaya S, Murphy JD, Rose BS, Wu J, Noticewala S, McHale MT, Yashar CM, Vaida F, Mell LK. (2014) Validated competing event model for the stage I-II endometrial cancer population. Int J Radiat Oncol Biol Phys.89:888-98.
#' \item Carmona R, Green GB, Zakeri K, Gulaya S, Xu B, Verma R, Williamson C, Rose BS, Murphy JD, Vaida F, Mell LK. (2015) Novel method to stratify elderly patients with head and neck cancer. J Clin Oncol 33 (suppl; abstr 9534).
#' \item Carmona R, Zakeri K, Green GB, Triplett DP, Murphy JD, Mell LK. (2015) Novel method to stratify elderly patients with prostate cancer. J Clin Oncol 33 (suppl; abstr 9532).
#' }
#' @return
#' \item{$coef1}{generalized competing event model coefficients (log (\eqn{\omega} ratio))}
#' \item{$coef2}{generalized competing event model coefficients (log (\eqn{\omega+} ratio))}
#' \item{$result1}{result table for generalized competing event model containing exponential of coefficients (\eqn{\omega} ratio) and 95\% confidence intervals}
#' \item{$result2}{result table for generalized competing event model containing exponential of coefficients (\eqn{\omega+} ratio) and 95\% confidence intervals}
#' \item{$omegaplot1}{\eqn{\omega} plot for generalized  competing evet model}
#' \item{$omegaplot2}{\eqn{\omega+} plot for generalized  competing evet model}
#' \item{$omegaplot3}{plot of \eqn{\omega} vs time}
#' @export


#### gce function by cox
gcecox <- function(formula1, formula2, formula3, surv1, surv2, data, N, M, t)
{
  set.seed(seed = 2015)

  # coefficients
  covnames <- attr(terms(formula1), "term.labels")

  CA_CPH <- coxph(formula1, data)
  CM_CPH <- coxph(formula2, data)
  All_CPH <- coxph(formula3, data)
  Beta1 <- CA_CPH$coef
  Beta2 <- CM_CPH$coef
  Beta <- All_CPH$coef

  Beta12 <- matrix(nrow = length(covnames), ncol = 1, data = 0)
  Betanew <- matrix(nrow = length(covnames), ncol = 1, data = 0)
  rownames(Beta12) <- covnames
  colnames(Beta12) <- " "
  rownames(Betanew) <- covnames
  colnames(Betanew) <- " "
  for (i in 1:length(covnames)){
    Beta12[i,] <- coef(CA_CPH)[i] - coef(CM_CPH)[i]
  }
  for (i in 1:length(covnames)){
    Betanew[i,] <- coef(CA_CPH)[i] - coef(All_CPH)[i]
  }

  Beta12 <- round(Beta12,5)
  Betanew <- round(Betanew,5)

  # variance
  var1 <- diag(CA_CPH$var)
  var2 <- diag(CM_CPH$var)
  var3 <- diag(All_CPH$var)

  coefvar1.boot <- matrix(nrow = N,ncol = length(covnames), data = NA)
  colnames(coefvar1.boot) <- covnames
  for (j in 1:N){
    for (i in 1:length(covnames)){
      Beta1dist <- sample(rnorm(1000, Beta1[i],sqrt(var1)[i]), 1000, replace = T)
      Betadist <- sample(rnorm(1000, Beta[i],sqrt(var3)[i]), 1000, replace = T)
      Betanewdist <- Beta1dist - Betadist
      coefvar1.boot[j,i] <- var(Betanewdist)
    }
  }
  coefvar1 <- matrix(nrow = length(covnames),ncol = 1, data = NA)
  rownames(coefvar1) <- covnames
  colnames(coefvar1) <- " "
  for (i in 1:length(covnames)){
    coefvar1[i,] <- mean(coefvar1.boot[,i])
  }


  coefvar2.boot <- matrix(nrow = N,ncol = length(covnames), data = NA)
  colnames(coefvar2.boot) <- covnames
  for (j in 1:N){
    for (i in 1:length(covnames)){
      Beta1dist <- sample(rnorm(1000, Beta1[i],sqrt(var1)[i]), 1000, replace = T)
      Beta2dist <- sample(rnorm(1000, Beta2[i],sqrt(var2)[i]), 1000, replace = T)
      Beta12dist <- Beta1dist - Beta2dist
      coefvar2.boot[j,i] <- var(Beta12dist)
    }
  }
  coefvar2 <- matrix(nrow = length(covnames),ncol = 1, data = NA)
  rownames(coefvar2) <- covnames
  colnames(coefvar2) <- " "
  for (i in 1:length(covnames)){
    coefvar2[i,] <- mean(coefvar2.boot[,i])
  }


  # 95% CI
  ci1 <- matrix(NA,length(covnames),2)
  colnames(ci1) <- c("2.5%","97.5%")
  rownames(ci1) <- covnames
  df <- dim(data)[1]
  for (i in 1:length(covnames)) {
    ci1[i,1] <- (Betanew[i]+c(-1,1)*qt(0.975, df = df)*sqrt(coefvar1[i]))[1]
    ci1[i,2] <- (Betanew[i]+c(-1,1)*qt(0.975, df = df)*sqrt(coefvar1[i]))[2]
  }

  ci2 <- matrix(NA,length(covnames),2)
  colnames(ci2) <- c("2.5%","97.5%")
  rownames(ci2) <- covnames
  df <- dim(data)[1]
  for (i in 1:length(covnames)) {
    ci2[i,1] <- (Beta12[i]+c(-1,1)*qt(0.975, df = df)*sqrt(coefvar2[i]))[1]
    ci2[i,2] <- (Beta12[i]+c(-1,1)*qt(0.975, df = df)*sqrt(coefvar2[i]))[2]
  }


  # Omega plot for B1-B
  # Competing event risk score
  riskscorenew <- numeric(dim(data)[1])
  for (i in 1:length(covnames)){
    riskscorenew <- riskscorenew + Betanew[i]*with(data,get(covnames[i]))
  }

  # Omega values and plots for competing event risk score
  data$normCER <- (riskscorenew - mean(riskscorenew))/sd(riskscorenew)
  l <- quantile(data$normCER, prob=seq(0,1,1/M))
  data$normCERomega <- M
  for (i in 1:M){
    data$normCERomega[data$normCER >= l[i] & data$normCER < l[i+1]] <- i
  }

  omegas <- seq(0,0,len=M)

  for (i in 1:M){
    # cum hazard for cancer death
    om.fit1 <- summary(survfit(surv1, data[data$normCERomega == i,]))
    H.hat.cancer.death <- -log(om.fit1$surv)
    yca <- H.hat.cancer.death[is.finite(H.hat.cancer.death)]
    xca <- om.fit1$time[is.finite(H.hat.cancer.death)]
    fitlm <- lm(yca~xca)
    point <- data.frame(xca = t)
    H.hat.catime <- predict(fitlm, point, interval ="prediction")[1]


    # cum hazard for competing event
    om.fit2 <- summary(survfit(surv2,data[data$normCERomega == i,]))
    H.hat.competing.event <- -log(om.fit2$surv)
    ycm <- H.hat.competing.event[is.finite(H.hat.competing.event)]
    xcm <- om.fit2$time[is.finite(H.hat.competing.event)]
    fitlm2 <- lm(ycm~xcm )
    point2 <- data.frame(xcm = t)
    H.hat.cmtime <- predict(fitlm2, point2, interval ="prediction")[1]


    omegas[i] <- H.hat.catime/(H.hat.catime + H.hat.cmtime)
  }


  y1 <- omegas
  x1 <- seq(min(data$normCER), max(data$normCER), len = M)
  z1 <- qplot(x1,y1, xlab = "Risk Score", ylab = expression(omega))



  # Omega plot for B1-B2
  # Competing event risk score
  riskscore12 <- numeric(dim(data)[1])
  for (i in 1:length(covnames)){
    riskscore12 <- riskscore12 + Beta12[i]*with(data,get(covnames[i]))
  }

  # Omega values and plots for competing event risk score
  data$normCER <- (riskscore12 - mean(riskscore12))/sd(riskscore12)
  l <- quantile(data$normCER, prob=seq(0,1,1/M))
  data$normCERomega <- M
  for (i in 1:M){
    data$normCERomega[data$normCER >= l[i] & data$normCER < l[i+1]] <- i
  }

  omegas <- seq(0,0,len=M)

  for (i in 1:M){
    # cum hazard for cancer death
    om.fit1 <- summary(survfit(surv1, data[data$normCERomega == i,]))
    H.hat.cancer.death <- -log(om.fit1$surv)
    yca <- H.hat.cancer.death[is.finite(H.hat.cancer.death)]
    xca <- om.fit1$time[is.finite(H.hat.cancer.death)]
    fitlm <- lm(yca~xca)
    point <- data.frame(xca = t)
    H.hat.catime <- predict(fitlm, point, interval ="prediction")[1]


    # cum hazard for competing event
    om.fit2 <- summary(survfit(surv2,data[data$normCERomega == i,]))
    H.hat.competing.event <- -log(om.fit2$surv)
    ycm <- H.hat.competing.event[is.finite(H.hat.competing.event)]
    xcm <- om.fit2$time[is.finite(H.hat.competing.event)]
    fitlm2 <- lm(ycm~xcm )
    point2 <- data.frame(xcm = t)
    H.hat.cmtime <- predict(fitlm2, point2, interval ="prediction")[1]


    omegas[i] <- H.hat.catime/H.hat.cmtime
  }

  y2 <- omegas
  x2 <- seq(min(data$normCER), max(data$normCER), len = M)
  z2 <- qplot(x2,y2, xlab = "Risk Score", ylab = expression(omega("+")))

  # Omega vs Time plot
  omegas <- numeric(t)

  for (i in 1:t){
    om.fit1 <- summary(survfit(surv1, data))
    H.hat.cancer.death <- -log(om.fit1$surv)
    yca <- H.hat.cancer.death[is.finite(H.hat.cancer.death)]
    xca <- om.fit1$time[is.finite(H.hat.cancer.death)]
    fitlm <- lm(yca~xca)
    point <- data.frame(xca = i)
    H.hat.catime <- predict(fitlm, point, interval ="prediction")[1]

    om.fit2 <- summary(survfit(surv2,data))
    H.hat.competing.event <- -log(om.fit2$surv)
    ycm <- H.hat.competing.event[is.finite(H.hat.competing.event)]
    xcm <- om.fit2$time[is.finite(H.hat.competing.event)]
    fitlm2 <- lm(ycm~xcm )
    point2 <- data.frame(xcm = i)
    H.hat.cmtime <- predict(fitlm2, point2, interval ="prediction")[1]


    omegas[i] <- H.hat.catime/(H.hat.catime + H.hat.cmtime)
  }

  y3 <- omegas
  x3 <- 1:t
  z3 <- qplot(x3,y3, ylim = c(0,2*max(y3)), xlab = "Time", ylab = expression(omega))

  # result tables
  table1 <- matrix(0,length(covnames),3)
  rownames(table1) <- covnames
  colnames(table1) <- c("exp(coef)", "lower .95", "upper .95")
  table1[,1] <- round(exp(Betanew),5)
  table1[,2] <- round(exp(ci1),5)[,1]
  table1[,3] <- round(exp(ci1),5)[,2]

  table2 <- matrix(0,length(covnames),3)
  rownames(table2) <- covnames
  colnames(table2) <- c("exp(coef)", "lower .95", "upper .95")
  table2[,1] <- round(exp(Beta12),5)
  table2[,2] <- round(exp(ci2),5)[,1]
  table2[,3] <- round(exp(ci2),5)[,2]

  return(list(coef1 = Betanew, coef2 = Beta12, result1 = table1, result2 = table2, omegaplot1 = z1, omegaplot2 = z2, omegaplot3 = z3))

}
