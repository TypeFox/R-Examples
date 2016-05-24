
# Currently, the lowest dose group (whatever its value is) is the basis for iLOGEL and BMD estimates.

#' @title Dose-response Modeling with Smoothing Splines
#' @usage
#' dichot_smooth(dosecolumn = "", 
#'                targetcolumn = "", 
#'                k = 4, 
#'                return_predict = FALSE, 
#'                write_predict = TRUE,
#'                STD_bias = TRUE,
#'                data = NA)
#' @aliases dichot_smooth
#' @description
#' Generates a spline model given dose and target response columns.
#' @details
#' This function generates a spline model with the input dose and target response
#' columns, plots the spline-estimated dose-response function with its upper and lower
#' 95 percent confidence bounds in green and red respectively along with the actual data, and returns
#' key metrics related to the dose-response function.  Note that the confidence bounds depicted on the
#' plot are for the dose-response function itself, and not for the raw data.
#' 
#' The parameter 'k', defaulted to 4, defines the number of dimensions the spline function will use in
#' estimating the response relation.  With 2 reserved for each end of the smooth, the default
#' allows for 2 bends in the smooth.  In the case that this appears to overfit the data,
#' the user may choose to override the default to 3, which would allow only one bend.
#' @param dosecolumn   Name of dose column of interest in dataframe.
#' @param targetcolumn  Name of response column of interest in dataframe.
#' @param k  Dimension of the basis used to represent the smooth term; see Details.
#' @param return_predict  If TRUE (default FALSE), returns dataframe of predicted values.
#' @param write_predict  If TRUE (the default), writes the dataframe of predicted values to a .csv file in the working directory.
#' @param STD_bias  If TRUE (the default), calculates the slope transition dose, a bootstrapped and resource-intensive computation.
#' @param data   Input dataframe.
#' @return A plot of the spline-estimated dose-response function along with the actual data.
#' @importFrom stats binomial terms
#' @keywords internal

dichot_smooth <- function (dosecolumn="", 
                           targetcolumn="", 
                           k = 4, 
                           return_predict = FALSE,
                           write_predict = TRUE,
                           STD_bias = TRUE,
                           data=NA) {
  
  x <- data
  targetvariable <- x[,targetcolumn]
  dose <- x[,dosecolumn]
  alpha <- .10
  
  # v1.9.0 allow user to define dimension of basis to smooth term, k
  # mgcv::s does not evaluate variables
  if (k == 4) {
    spline <- mgcv::gam(targetvariable~s(dose, k=4), family=binomial, data=x)
  } else {
    if (k == 3) {
      spline <- mgcv::gam(targetvariable~s(dose, k=3), family=binomial, data=x)
    } else {
      stop("Dimension k of smooth term basis must be 3 or 4.")}
  }
  
  min <- min(dose)
  max <- max(dose)
  
  # Was (max-min)/1000; this new variant will produce steps of .001
  # when the dose range is >1 & <=10, .01 for range >10 & <=100, .0001 for range >.1 & <=1, etc.
  step <- .001*(10^(floor(log10(max-min))))  
  
  
  # !obs! as.ordered could change the order from that entered?
  dose_fac <- as.ordered(sort(dose))    # Needed to get the counts for the standard errors below.
  ns <- as.vector(summary(dose_fac))    # Gets a vector of counts for each dose group.
  n <- sum(ns)
  avgn <- n/nlevels(dose_fac)
  
  # This adds the unique values of dose to a new predictdata dataframe, gets rid of duplicates, and sorts it.
  predictdata <- data.frame(cbind(dose = sort(unique(c(seq(min, max, by = step),sort(unique(dose)))))))  
  
  CIfactor <- stats::qt(c(alpha/2), df = spline$df.residual, lower.tail=FALSE)
  predictnew <- stats::predict(spline, predictdata, type="response", se.fit=TRUE)
  predictdata$fit <- as.vector(predictnew$fit)
  predictdata$se <- as.vector(predictnew$se)  # Important: This is the standard error for a single case. 
  predictdata$lcl90 <- predictdata$fit-(CIfactor*predictdata$se)
  predictdata$ucl90 <- predictdata$fit+(CIfactor*predictdata$se)
  
  f <- get("firstDeriv", envir = environment(drsmooth))
  d1_spline <- f(spline, n = length(predictdata$dose)) 
  
  # Sets an eps at the same value as in the first.deriv function for the finite differencing that follows. 
  eps <- 1e-7   
  # The d1 values are in logit space.
  d1_logit <- as.vector(d1_spline$dose$deriv)  
  # This generates a vector, fit_logit, that puts the predicted values into logit space.
  # Packages boot & car now conflict (v.1.9.0) because specific package calls were not 
  # either specified or designated with @importFrom or using correct namespace definition. 
  fit_logit <- boot::logit(predictdata$fit)  
  fit_logit_at_eps_using_d1 <- fit_logit + (d1_logit*eps)
  # This gives the amount of change at that dose in response space for the finite difference increment, 
  # divided by the finite difference.
  # By adding the package call to inv.logit, we attempt to avoid breaking mid-stream when the wrong logit is found.  
  predictdata$d1 <- (boot::inv.logit(fit_logit_at_eps_using_d1) - boot::inv.logit(fit_logit))/eps  
  
  d1se_logit <- as.vector(d1_spline$dose$se.deriv)
  d1_lcl90_logit <- as.vector(d1_logit-(CIfactor*d1se_logit))
  d1_ucl90_logit <- as.vector(d1_logit+(CIfactor*d1se_logit))
  
  # The above generate the lcl and ucl values of the first derivative in logit space.
  # Now repeat the process above using those lcl and ucl slope values to get the bounds of the first derivative in response space.
  
  fit_logit_at_eps_using_d1_lcl90 <- fit_logit+(d1_lcl90_logit*eps)
  predictdata$d1_lcl90 <- (boot::inv.logit(fit_logit_at_eps_using_d1_lcl90) - boot::inv.logit(fit_logit))/eps
  
  fit_logit_at_eps_using_d1_ucl90 <- fit_logit+(d1_ucl90_logit*eps)
  predictdata$d1_ucl90 <- (boot::inv.logit(fit_logit_at_eps_using_d1_ucl90) - boot::inv.logit(fit_logit))/eps
  
  # Oddly, this is the most straightforward way of getting the se of d1 in response space.
  predictdata$d1se <- (predictdata$d1_ucl90-predictdata$d1)/CIfactor   
  
  predict <- predict(spline, type="response", se.fit=TRUE)
  x$fit <- as.vector(predict$fit)
  x$se <- as.vector(predict$se)
  x$lcl90 <- x$fit-(CIfactor*x$se)
  x$ucl90 <- x$fit+(CIfactor*x$se)
  
  zerodose90ucl <- predictdata$fit[1] + (CIfactor*((predictdata$fit[1]*(1-predictdata$fit[1]))^.5)/((n)^.5))  
  # Note that this is the standard error of the mean for the zeroes/ones raw data -- 
  # i.e., the standard deviation of the mean (the mean is the proportion itself).
  # Using the overall n for the study is appropriate here, rather than n only for zero-dose group.
  # Consider the situation where there aren't dose groups, but rather doses are simply scattered 
  # through the range (as would an IV in a survey situation) to intuitively understand why.
  # (The n at any particular dose isn't relevant -- the error is estimated for the whole dataset -- 
  # just as the residuals in a regular regression are estimated based on the total n.)  
  
  # AB v. 1.9.0 The following min calls often generate and print warnings.
  # Suppressing them here, following conversations with CT 03/17/15 prior to SOT demos 
  # and to serve as placeholders for future solutions that may include estimation 
  # outside the observed dose range.
  
  loaellogical <- which(predictdata$fit>zerodose90ucl)
  suppressWarnings(loael <- round(predictdata$dose[min(loaellogical)], digits = 4))
  
  # The following does BMDs for 1sd and for 10% extra risk.
  
  # See EPA definition in document in email from Chad Re: TiO2 data sent at 7/11/2014 at 9:30am.
  BMD1perc_response_target <- predictdata$fit[1] + ((1-predictdata$fit[1])*.01)   
  BMD1perclogical <- which(predictdata$fit>BMD1perc_response_target)
  suppressWarnings(BMD1perc <- round(predictdata$dose[min(BMD1perclogical)], digits = 4))
  # The dose at which the upper bound of the dose-response function intersects the BMR and slope>0.
  BMDL1perclogical <- which(predictdata$ucl90>BMD1perc_response_target&predictdata$d1_lcl90>0)  
  #	See http://onlinelibrary.wiley.com/doi/10.1002/jat.1298/pdf (Figure 2, especially)
  suppressWarnings(BMDL1perc <- round(predictdata$dose[min(BMDL1perclogical)], digits = 4))
  
  # AB adding 5% BMD/L
  # See EPA definition in document in email from Chad Re: TiO2 data sent at 7/11/2014 at 9:30am.
  BMD5perc_response_target <- predictdata$fit[1] + ((1-predictdata$fit[1])*.05)   
  BMD5perclogical <- which(predictdata$fit>BMD5perc_response_target)
  suppressWarnings(BMD5perc <- round(predictdata$dose[min(BMD5perclogical)], digits = 4))
  # The dose at which the upper bound of the dose-response function intersects the BMR and slope>0.
  BMDL5perclogical <- which(predictdata$ucl90>BMD5perc_response_target&predictdata$d1_lcl90>0)  
  #  See http://onlinelibrary.wiley.com/doi/10.1002/jat.1298/pdf (Figure 2, especially)
  suppressWarnings(BMDL5perc <- round(predictdata$dose[min(BMDL5perclogical)], digits = 4))
  
  # See EPA definition in document in email from Chad Re: TiO2 data sent at 7/11/2014 at 9:30am.
  BMD10perc_response_target <- predictdata$fit[1] + ((1-predictdata$fit[1])*.10)   
  BMD10perclogical <- which(predictdata$fit>BMD10perc_response_target)
  suppressWarnings(BMD10perc <- round(predictdata$dose[min(BMD10perclogical)], digits = 4))
  # The dose at which the upper bound of the dose-response function intersects the BMR and slope>0.
  BMDL10perclogical <- which(predictdata$ucl90>BMD10perc_response_target&predictdata$d1_lcl90>0)  
  #	See http://onlinelibrary.wiley.com/doi/10.1002/jat.1298/pdf (Figure 2, especially)
  suppressWarnings(BMDL10perc <- round(predictdata$dose[min(BMDL10perclogical)], digits = 4))
  
  # AB v. 1.9.0 fix plot errors:  "dose" and graphmax.
  # AB 4/24/15 graphmax appears to intend to test max of ucl90 > .9; corrected.
  # AB 6/5/15 all plotting now corrected and encapsulated in spline.plot
  # However, an issue remains in 1:1 dichotmous "human" data, i.e. one dose per response.
  # graphmax must be 1 to accomodate proportion 1.0 for plotting raw data in that case.
  # if (max(predictdata$ucl90<.9)) {graphmax <- (max(predictdata$ucl90)*.1) + max(predictdata$ucl90)} else {graphmax <- 1}
  # Note: adding this targetvariable check is only appropriate for dichotomous
  # if (max(predictdata$ucl90) < .9 & max(targetvariable) < .9) {graphmax <- (max(predictdata$ucl90)*.1) + max(predictdata$ucl90)} else {graphmax <- 1}
  
  if (is.na(loael)) {
    warning('No ST dose can be identified in these data.')
    # AB 3/17/15 Appears to assume input file dose column is called "dose" and therefore matching the
    # predictdata column called "dose" and failing.  Corrected to use column "dose" as defined in predictdata.
    # matplot(predictdata[,dosecolumn], predictdata[, c("fit","lcl90","ucl90")], type="l", pch=19, lty=1, ylim=c(min(x[,targetcolumn]),graphmax),
    #     matplot(predictdata[,"dose"], predictdata[, c("fit","lcl90","ucl90")], 
    #             type="l", pch=19, lty=1, ylim=c(min(x[,targetcolumn]),graphmax),
    #             xlab="Dose", ylab="Predictions, CIs, and Raw Data")
    #     matpoints(x[,dosecolumn], x[,targetcolumn], type="p", pch=19)
    drsmooth::spline.plot(dosecolumn, targetcolumn, data_type = "dichotomous", k=k, data=data)
  } else {
    
    lclloaellogical <- which(predictdata$ucl90>zerodose90ucl)
    suppressWarnings(lcl90loael <- round(predictdata$dose[min(lclloaellogical)], digits = 4))
    
    uclloaellogical <- which(predictdata$lcl90>zerodose90ucl)
    suppressWarnings(ucl90loael <- round(predictdata$dose[min(uclloaellogical)], digits = 4))
    
    tdlogical <- which(predictdata$d1_lcl90>0)
    suppressWarnings(td <- round(predictdata$dose[min(tdlogical)], digits = 4))
    
    n <- nrow(x)
    predictdata$se_of_d1se <- predictdata$d1se/(n^.5)
    predictdata$d1_lcl90_u <- predictdata$d1_lcl90+(1.644854*predictdata$se_of_d1se)
    predictdata$d1_lcl90_l <- predictdata$d1_lcl90-(1.644854*predictdata$se_of_d1se)
    
    lcl90tdlogical <- which(predictdata$d1_lcl90_u>0)
    suppressWarnings(lcl90td <- round(predictdata$dose[min(lcl90tdlogical)], digits = 4))
    
    ucl90tdlogical <- which(predictdata$d1_lcl90_l>0)
    suppressWarnings(ucl90td <- round(predictdata$dose[min(ucl90tdlogical)], digits = 4))
    
    # The following line initializes a vector for storing the td values from the bootstrap process to follow.
    # AB 6/5/15 v1.9. Make this calculation driven by parameter STD_bias.
    if (STD_bias == TRUE) {
      tdsampleresults <- as.vector(rep(NA, 1000))
      
      # And now the loop to perform the bootstrap.
      for (i in 1:1000) {
        
        sampledata <- x[sample(nrow(x),replace=TRUE),]
        sampletargetvariable <- sampledata[,targetcolumn]
        sampledosevariable <- sampledata[,dosecolumn]
        
        # v1.9.0 allow user to define dimension of basis to smooth term, k
        # mgcv::s does not evaluate variables
        if (k == 4) {
          samplespline <- mgcv::gam(targetvariable~s(dose, k=4), family=binomial, data=x)
        } else {
          if (k == 3) {
            samplespline <- mgcv::gam(targetvariable~s(dose, k=3), family=binomial, data=x)
          } 
        } 
          
        samplemin <- min(sampledata[,dosecolumn])
        samplemax <- max(sampledata[,dosecolumn])
        samplestep <- .001*(10^(floor(log10(samplemax-samplemin)))) 
        predictsampledata <- data.frame(cbind(dose = seq(samplemin, samplemax, by = samplestep)))
        CIfactor <- stats::qt(c(alpha/2), df = samplespline$df.residual, lower.tail=FALSE)
        predictsamplenew <- predict(samplespline, predictsampledata, type="response", se.fit=TRUE)
        predictsampledata$fit <- as.vector(predictsamplenew$fit)
        predictsampledata$se <- as.vector(predictsamplenew$se)
        
        predictsampledata$lcl90 <- predictsampledata$fit-(CIfactor*predictsampledata$se)
        predictsampledata$ucl90 <- predictsampledata$fit+(CIfactor*predictsampledata$se)
        
        f <- get("firstDeriv", envir = environment(drsmooth))
        d1_samplespline <- f(samplespline, n = length(predictsampledata$dose))
        
        eps <- 1e-7   # Sets an eps at the same value as in the first.deriv function for a finite differencing to follow. 
        d1_logit <- as.vector(d1_samplespline$dose$deriv)  # The d1 values are in logit space.
        # Suppress warnings for SOT demo
        # This generates a vector, fit_logit, that puts the predicted values into logit space.
        suppressWarnings(fit_logit <- boot::logit(predictsampledata$fit)) 
        fit_logit_at_eps_using_d1 <- fit_logit + (d1_logit*eps)
        # This gives the amount of change at that dose in response space for the finite difference increment, divided by the finite difference.
        predictsampledata$d1 <- (boot::inv.logit(fit_logit_at_eps_using_d1) - boot::inv.logit(fit_logit))/eps  
        d1se_logit <- as.vector(d1_samplespline$dose$se.deriv)
        d1_lcl90_logit <- as.vector(d1_logit-(CIfactor*d1se_logit))
        d1_ucl90_logit <- as.vector(d1_logit+(CIfactor*d1se_logit))
        
        # The above generate the lcl and ucl values of the first derivative in logit space.
        # Now repeat the process above using those lcl and ucl slope values to get the bounds of the first derivative in response space.
        
        fit_logit_at_eps_using_d1_lcl90 <- fit_logit+(d1_lcl90_logit*eps)
        predictsampledata$d1_lcl90 <- (boot::inv.logit(fit_logit_at_eps_using_d1_lcl90) - boot::inv.logit(fit_logit))/eps
        
        fit_logit_at_eps_using_d1_ucl90 <- fit_logit+(d1_ucl90_logit*eps)
        predictsampledata$d1_ucl90 <- (boot::inv.logit(fit_logit_at_eps_using_d1_ucl90) - boot::inv.logit(fit_logit))/eps
        
        # Oddly, this is the most straightforward way of getting the se of d1 in response space.
        predictsampledata$d1se <- (predictsampledata$d1_ucl90-predictsampledata$d1)/CIfactor   
        
        tdsamplelogical <- which(predictsampledata$d1_lcl90>0)
        suppressWarnings(tdsample <- predictsampledata$dose[min(tdsamplelogical)])
        
        tdsampleresults[i] <- tdsample
      }
      
      tdsampleresults_s <- sort(tdsampleresults)
      
      tdbias <- round(td-mean(tdsampleresults_s), digits = 3)
    } else {
      tdbias = NA
    }
    # AB 3/17/15 Appears to assume input file dose column is called "dose" and therefore matching the
    # predictdata column called "dose" and failing.  Corrected to use column "dose" as defined in predictdata.
    # AB 3/24/15 graphmax already corrected above.
    # AB 6/5/15 all plotting now corrected and encapsulated in spline.plot
    #matplot(predictdata[,dosecolumn], predictdata[, c("fit","lcl90","ucl90")], type="l", pch=19, lty=1, ylim=c(min(x[,targetcolumn]),graphmax),
    #     matplot(predictdata[,"dose"], predictdata[, c("fit","lcl90","ucl90")], type="l", pch=19, lty=1, ylim=c(min(x[,targetcolumn]),graphmax),
    #             xlab="Dose", ylab="Predictions, CIs, and Raw Data")
    #     
    #     doses <- as.vector(sort(unique(dose)))
    #     responseproportions <- tapply(targetvariable, as.factor(sort(dose)),mean)
    #     
    #     matpoints(doses, responseproportions, type="p", pch=19)
    
    drsmooth::spline.plot(dosecolumn, targetcolumn, data_type = "dichotomous", k=k, data=data)
    
    neg2ll <- -2 * stats::logLik(spline)
    # R's AIC function counts the intercept as a parameter, so to be consistent we need to add 1 here,
    # because summary(spline)$edf doesn't include it.
    df <- summary(spline)$edf + 1   
    aic <- stats::AIC(spline)
    
    output <- c("STD"=td, "STD_l"=lcl90td, "STD_u"=ucl90td, "STD_bias"=tdbias, "iLOGEL"=loael, "iLOGEL_l"=lcl90loael, "iLOGEL_u"=ucl90loael)
    output2 <- c("BMD1perc"=BMD1perc, "BMDL1perc"=BMDL1perc, "BMD5perc"=BMD5perc, "BMDL5perc"=BMDL5perc, "BMD10perc"=BMD10perc, "BMDL10perc"=BMDL10perc)
    output3 <- c("-2LL"=neg2ll, "df"=df, "AIC"=aic)
    
    print(output)
    cat("\n")
    print(output2)
    cat("\n")
    print(output3)
    cat("\n")
    cat("NOTE: The AIC is -2LL+(2*df), where df (the model degrees of freedom, 
            or total number of model parameters) includes the intercept term.")
    # v.1.9.0 requested unique and meaningful output file name
    # write.csv(predictdata, "predictions_etc.csv", row.names=FALSE)
    if (write_predict == TRUE) {
      now <- Sys.time()
      file_prefix <- paste(deparse(substitute(data)), "predictions", lubridate::now(), sep = " ")
      file_name <- paste(file_prefix, ".csv", sep = "")
      utils::write.csv(predictdata, file_name, row.names=FALSE)}
    if (return_predict == TRUE) {return(predictdata)}
  }
}

