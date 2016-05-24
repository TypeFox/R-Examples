#' Fit Generalized Competing Event Model Based on Fine Gray Regression
#'
#' @description Fit a generalized competing event model by using Fine Gray
#' regression model with \code{crr} function in \code{cmprsk} package.
#' @param ostime1 vector of times for event(s) of interest.
#' @param ostime2 vector of times for competing event(s).
#' @param ostime3 vector of times for the composute set of all events.
#' @param cod1 vector with 1 for event(s) pf interest, 2 for competing event(s) and 0 for censored observations.
#' @param cod2 vector with 1 for the composite set of all events and 0 for censored observations.
#' @param data a data frame containing all covariates.
#' @param covnames vector of names for all covariates in data.
#' @param N the number of bootstrap replicates
#' @param M the number of bins for \eqn{\omega} or \eqn{\omega+} plots.
#' @param t survival time point for \eqn{\omega} or \eqn{\omega+} plots.
#' @import survival cmprsk ggplot2 stats
#' @details The \strong{gcerisk} package is designed to help investigators optimize risk-stratification methods for competing risks data, such as described in
#' Carmona R, Gulaya S, Murphy JD, Rose BS, Wu J, Noticewala S, McHale MT, Yashar CM, Vaida F, Mell LK. Validated competing event model for the stage I-II endometrial cancer population.
#' Int J Radiat Oncol Biol Phys. 2014;89:888-98. Standard risk models typically estimate the effects of one or more covariates on either
#' a single event of interest (such as overall mortality, or disease recurrence), or a composite set of events (e.g., disease-free survival, which combines events of interest with death from any cause).
#' This method is inefficient in stratifying patients who may be simultaneously at high risk for the event of interest but low risk for competing events, and who thus stand to gain the most from strategies to modulate the event of interest.
#' Compared to standard risk models, GCE models better stratify patients at higher (lower) risk for an event of interest and lower (higher) risk of competing events. GCE models focus on differentiating subjects based on
#' the ratio of the cumulative hazard (or cumulative hazard of the subdistribution) for the event of interest to the cumulative hazard (or cumulative hazard of the subdistribution) for all events (\eqn{\omega}),
#' and the ratio of the cumulative hazard (or cumulative hazard of the subdistribution) for the event of interest to the cumulative hazard (or cumulative hazard of the subdistribution) for competing events (\eqn{\omega+}).
#'
#' The \code{gcefg} function produces model estimates and confidence intervals from a generalized competing event model based on the Fine-Gray model for subdistribution hazards. In the subdistribution hazards model, the function H(t)= -log(1-F(t)) represents the cumulative hazard of the subdistribution for the cumulative distribution function F(t).
#' The model assumes proportional subdistribution hazards for the composite set of events.
#'
#' The function returns \eqn{\omega} and \eqn{\omega+} ratio estimates for the chosen covariates, with 95\% confidence intervals, and plots \eqn{\omega} and \eqn{\omega+} at time t within M ordered subsets of subjects as a function of increasing risk (based on the linear predictor, i.e. the inner product of a subject's data vector and the coefficient vector).
#' @examples
#' # sample data to test
#' data(Sample)
#' test <- Sample
#' d <- trunc(dim(test)[1]*0.1)
#' set.seed(seed=2015)
#' s <- sample(dim(test)[1],d,replace = FALSE)
#' test <- test[s,]
#' rm(list=setdiff(ls(), "test"))
#' test <- transform(test, LRF_OR_DF_FLAG = as.numeric(test$LRFFLAG | test$DFFLAG))
#' test <- transform(test, LRF_OR_DF_MO = pmin(test$LRFMO, test$DFMO))
#' test <- transform(test, CMFLAG = as.numeric(test$OSFLAG & !test$LRFFLAG & !test$DFFLAG))
#' test <- transform(test, ACMFLAG = as.numeric(test$LRF_OR_DF_FLAG | test$CMFLAG))
#' test <- transform(test, ACM_MO = pmin(test$LRF_OR_DF_MO, test$OSMO))
#'
#' cod1 <- test$ACMFLAG
#' cod1[test$LRF_OR_DF_FLAG == 1] <- 1
#' cod1[test$CMFLAG == 1] <- 2
#' cod2 <- test$ACMFLAG
#' ostime1 <- test$LRF_OR_DF_MO
#' ostime2 <- test$OSMO
#' ostime3 <- test$ACM_MO
#' covnames <- c("age", "gender", "smoke20", "etohheavy", "higrade", "BMI", "black")
#'
#' N <- 50
#' M <- 5
#' t <- 60
#'
#' fitgce.fg <- gcefg(ostime1, ostime2, ostime3, cod1, cod2, test, covnames, N, M, t)

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


#### gce function by fine gray
gcefg <- function(ostime1, ostime2, ostime3, cod1, cod2, data, covnames, N, M, t)
{
  set.seed(seed = 2015)

  # coefficients
  covdata <- data[,covnames[1:length(covnames)]]
  fit1 <- crr(ostime1,cod1,covdata,failcode=1)
  fit2 <- crr(ostime2,cod1,covdata,failcode=2)
  fit3 <- crr(ostime3,cod2,covdata,failcode=1)
  Beta1 <- fit1$coef
  Beta2 <- fit2$coef
  Beta <- fit3$coef

  Beta12 <- matrix(nrow = length(covnames),ncol = 1, data = NA)
  rownames(Beta12) <- covnames
  colnames(Beta12) <- " "
  for (i in 1:length(covnames)){
    Beta12[i,] <- fit1$coef[i] - fit2$coef[i]
  }
  Betanew <- matrix(nrow = length(covnames),ncol = 1, data = NA)
  rownames(Betanew) <- covnames
  colnames(Betanew) <- " "
  for (i in 1:length(covnames)){
    Betanew[i,] <- fit1$coef[i] - fit3$coef[i]
  }

  Beta12 <- round(Beta12,5)
  Betanew <- round(Betanew,5)

  # variance
  var1 <- diag(fit1$var)
  var2 <- diag(fit2$var)
  var3 <- diag(fit3$var)

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


  ## Omega plot for B1-B
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

  omegas <- seq(0,0,len = M)

  for (i in 1:M){

    l <- dim(data[data$normCERomega == i,])[1]

    # cum hazard for cancer death
    H.hat.catime <- rep(0,l)
    for (j in 1:l){
      yca <- predict(fit1, as.matrix(covdata[data$normCERomega == i,])[j,])[,2]
      yca <- -log(1-yca)
      xca <- predict(fit1, as.matrix(covdata[data$normCERomega == i,])[j,])[,1]
      fitlm <- lm(yca~xca)
      point <- data.frame(xca = t)
      H.hat.catime[j] <- predict(fitlm, point, interval ="prediction")[1]
    }

    # cum hazard for all events
    H.hat.alltime <- rep(0,l)
    for (j in 1:l){
      yall <- predict(fit3, as.matrix(covdata[data$normCERomega == i,])[j,])[,2]
      yall <- -log(1-yall)
      xall <- predict(fit3, as.matrix(covdata[data$normCERomega == i,])[j,])[,1]
      fitlm <- lm(yall~xall)
      point <- data.frame(xall = t)
      H.hat.alltime[j] <- predict(fitlm, point, interval ="prediction")[1]
    }

    omegas[i] <- mean(H.hat.catime)/mean(H.hat.alltime)
  }


  N1 <- table(data$normCERomega)
  y1 <- omegas
  x1 <- seq(min(data$normCER), max(data$normCER), len = M)
  z1 <- qplot(x1, y1, xlab = "Risk Score", ylab = expression(omega))

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


  omegas <- seq(0,0,len = M)

  for (i in 1:M){

    l <- dim(data[data$normCERomega == i,])[1]

    # cum hazard for cancer death
    H.hat.catime <- rep(0,l)
    for (j in 1:l){
      yca <- predict(fit1, as.matrix(covdata[data$normCERomega == i,])[j,])[,2]
      yca <- -log(1-yca)
      xca <- predict(fit1, as.matrix(covdata[data$normCERomega == i,])[j,])[,1]
      fitlm <- lm(yca~xca)
      point <- data.frame(xca = t)
      H.hat.catime[j] <- predict(fitlm, point, interval ="prediction")[1]
    }


    # cum hazard for competing event
    H.hat.cmtime <- rep(0,l)
    for (j in 1:l){
      ycm <- predict(fit2, as.matrix(covdata[data$normCERomega == i,])[j,])[,2]
      ycm <- -log(1-ycm)
      xcm <- predict(fit2, as.matrix(covdata[data$normCERomega == i,])[j,])[,1]
      fitlm <- lm(ycm~xcm)
      point <- data.frame(xcm = t)
      H.hat.cmtime[j] <- predict(fitlm, point, interval ="prediction")[1]
    }

    omegas[i] <- mean(H.hat.catime)/mean(H.hat.cmtime)
  }


  N2 <- table(data$normCERomega)
  y2 <- omegas
  x2 <- seq(min(data$normCER), max(data$normCER), len = M)
  z2 <- qplot(x2,y2, xlab = "Risk Score", ylab = expression(omega("+")))

  # Omega vs Time plot
  omegas <- numeric(t)

  for (i in 1:t){
    H.hat.catime <- rep(0,l)
    for (j in 1:l){
      yca <- predict(fit1, as.matrix(covdata)[j,])[,2]
      yca <- -log(1-yca)
      xca <- predict(fit1, as.matrix(covdata)[j,])[,1]
      fitlm <- lm(yca~xca)
      point <- data.frame(xca = i)
      H.hat.catime[j] <- predict(fitlm, point, interval ="prediction")[1]
    }

    H.hat.alltime <- rep(0,l)
    for (j in 1:l){
      yall <- predict(fit3, as.matrix(covdata)[j,])[,2]
      yall <- -log(1-yall)
      xall <- predict(fit3, as.matrix(covdata)[j,])[,1]
      fitlm2 <- lm(yall~xall)
      point2 <- data.frame(xall = i)
      H.hat.alltime[j] <- predict(fitlm2, point2, interval ="prediction")[1]
    }

    omegas[i] <- mean(H.hat.catime)/mean(H.hat.alltime)
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



