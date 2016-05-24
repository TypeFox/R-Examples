#' @title Dose-response Modeling with Smoothing Splines
#' @usage
#' drsmooth(dosecolumn = "", 
#'          targetcolumn = "", 
#'          data_type = "", 
#'          k = 4,
#'          return_predict = FALSE, 
#'          write_predict = TRUE, 
#'          STD_bias = TRUE,
#'          data = NA)
#' @aliases drsmooth
#' @description
#' Generates a spline model given dose and target response columns.
#' @details
#' This function generates a spline model with the input dose and target response
#' columns, plots the spline-estimated dose-response function with its upper and lower
#' 95 percent confidence bounds in green and red respectively along with the actual data, and returns
#' key metrics related to the dose-response function.  Currently, the program will use 
#' the lowest dose group (whatever its value is) as the basis for iLOGEL and BMD estimates
#' Normally we would assume that would be zero, but it is not forced to be.
#' 
#' Note that the confidence bounds depicted on the plot are for the dose-response function itself, 
#' and not for the raw data.
#' 
#' The parameter 'k', defaulted to 4, defines the dimensions allowed the smooth term in
#' estimating the response relation.  With 2 reserved for each end of the smooth, the default
#' allows for 2 bends in this univariate smooth.  In the case that this appears to overfit the data,
#' the user may choose to override the default to 3, which would allow only one bend.
#' @param dosecolumn   Name of dose column of interest in dataframe.
#' @param targetcolumn  Name of response column of interest in dataframe.
#' @param data_type  Allowed values "continuous" or "dichotomous".
#' @param k  Dimension of the basis used to represent the smooth term; see Details.
#' @param return_predict  If TRUE (default FALSE), returns dataframe of predicted values.
#' @param write_predict  If TRUE (the default), writes the dataframe of predicted values to a .csv file in the working directory.
#' @param STD_bias  If TRUE (the default), calculates the bias of the slope transition dose, a bootstrapped and resource-intensive computation.
#' @param data   Input dataframe.
#' @return A plot of the spline-estimated dose-response function along with the actual data.
#' Also, several key metrics are reported:
#'
#' STD (slope transition dose): The lowest dose at which the slope of the dose-response function
#' is significantly (90% two-sided) positive.
#'
#' STD_l and STD_u: The 90 percent lower and upper confidence bounds on the STD.
#'
#' STD_bias (experimental): An estimate of the bias associated with the STD:  the difference between the point estimate and
#' the mean of 1000 bootstrapped STDs.
#' 
#' iLOGEL (experimental: interpolated lowest observed effect level) The lowest dose at which the predicted response
#' exceeds the 90 percent upper confidence bound of the response at zero dose.  This value
#' can be anywhere within the dose range -- hence "interpolated."
#'
#' iLOGEL_l and iLOGEL_u: The 90 percent lower and upper confidence bounds on the iLOGEL.
#'
#' For data_type = "continuous":
#' BMD1sd and BMD10: Benchmark doses corresponding to a 1sd and 10% increase in control response, respectively.
#' BMDL1sd and BMDL10: 90 percent (two-sided) lower bounds on the indicated BMDs.
#' 
#' For data_type = "dichotomous":
#' BMD1perc, BMD5perc, BMD10perc: Benchmark doses corresponding to a 1%, 5%, and 10% increase in control response, respectively.
#' BMDL1perc, BMDL5perc, BMDL10perc: 90 percent (two-sided) lower bounds on the indicated BMDs.
#'
#' -2LL, the number of parameters associated with the spline model, and the AIC.
#'
#' @examples
#' # Produces and plots spline model with confidence bounds, and prints key metrics.
#' # For the plot only, see spline.plot
#' # The STD_bias is defaulted here to FALSE to speed run time.
#' # For continuous outcomes
#' data(DRdata)
#' drsmooth("dose", "MF_Log", data_type = "continuous", k = 4, return_predict = FALSE, 
#' write_predict = FALSE, STD_bias = FALSE, data=DRdata)
#' 
#' # For dichotomous outcomes
#' data(DIdata)
#' # If necessary, convert summarized dataframe into 1 row per case dataframe (see drsmooth::expand)
#' DIdata_expanded <- expand(dosecolumn = "Dose", targetcolumn = "Tumor", ncolumn = "n", data = DIdata)
#' 
#' drsmooth("Dose", "Tumor", data_type = "dichotomous", return_predict = FALSE, 
#' write_predict = FALSE, STD_bias = FALSE, data=DIdata_expanded)
#' @export

drsmooth <- function (dosecolumn="", 
                      targetcolumn="", 
                      data_type = "continuous", 
                      k = 4,
                      return_predict = FALSE, 
                      write_predict = TRUE, 
                      STD_bias = TRUE,
                      data=NA) {
  
  if (data_type == "dichotomous") {
    dichot_smooth(dosecolumn = dosecolumn, 
                  targetcolumn = targetcolumn, 
                  k = k,
                  return_predict = return_predict,
                  write_predict = write_predict,
                  STD_bias = STD_bias,
                  data = data)} 
  else {
    if (data_type == "continuous") {
      x <- data
      targetvariable <- x[,targetcolumn]
      dose <- x[,dosecolumn]
      alpha <- .10
      
      # v1.9.0 allow user to define dimension of basis to smooth term, k
      # mgcv::s does not evaluate variables
      if (k == 4) {
        spline <- mgcv::gam(targetvariable~s(dose, k=4), data=x) } else
        {if (k == 3) {
          spline <- mgcv::gam(targetvariable~s(dose, k=3), data=x)
        } else {
          stop("Dimension k of smooth term basis must be 3 or 4.")}
        }
      
      min <- min(dose)
      max <- max(dose)
      # Was (max-min)/1000; this new variant will produce steps of .001 when the dose range
      #  is >1 & <=10, .01 for range >10 & <=100, .0001 for range >.1 & <=1, etc. 
      step <- .001*(10^(floor(log10(max-min))))   
      
      # This adds the unique values of dose to the predictdata dataframe, gets rid of duplicates, and sorts it.
      predictdata <- data.frame(cbind(dose = sort(unique(c(seq(min, max, by = step),sort(unique(dose)))))))  
      CIfactor <- stats::qt(c(alpha/2), df = spline$df.residual, lower.tail=FALSE)
      predictnew <- predict(spline, predictdata, type="response", se.fit=TRUE)
      predictdata$fit <- as.vector(predictnew$fit)
      predictdata$se <- as.vector(predictnew$se)
      predictdata$lcl90 <- predictdata$fit-(CIfactor*predictdata$se)
      predictdata$ucl90 <- predictdata$fit+(CIfactor*predictdata$se)
      
      # Given the newly variable length of the dose vector in predictdata, this needed to be changed.
      f <- get("firstDeriv", envir = environment(drsmooth))
      d1_spline <- f(spline, n = length(predictdata$dose))   
      
      predictdata$d1 <- as.vector(d1_spline$dose$deriv)
      predictdata$d1se <- as.vector(d1_spline$dose$se.deriv)
      predictdata$d1_lcl90 <- as.vector(predictdata$d1-(CIfactor*predictdata$d1se))
      predictdata$d1_ucl90 <- as.vector(predictdata$d1+(CIfactor*predictdata$d1se))
      
      predict <- predict(spline, type="response", se.fit=TRUE)
      x$fit <- as.vector(predict$fit)
      x$se <- as.vector(predict$se)
      x$lcl90 <- x$fit-(CIfactor*x$se)
      x$ucl90 <- x$fit+(CIfactor*x$se)
      
      # AB v1.9.0 Some of the following min calls often generate and print warnings.
      # Suppressing them here, to serve as placeholders for future solutions that may include estimation 
      # outside the observed dose range or returning informative messages, i.e.
      # "...falls (below/above) the observed dose range ..."
      
      # AB SOT testing.  The original code (uncommented) appears to assume that the first fitted
      # value is the lowest.  This is not always the case.
      # AB proposed fix, currently commented out, for initial negative d-r does not assume first fitted value is the lowest.
      # zerodose90ucl <- min(predictdata$fit) + (CIfactor*sd(spline$residuals)) 
      
      # v1.9 commented out and changed to line following. GH bug fix of iLogel CI.
      # zerodose90ucl <- predictdata$ucl90[1]
      zerodose90ucl <- predictdata$fit[1] + (CIfactor*stats::sd(spline$residuals))
      loaellogical <- which(predictdata$fit>zerodose90ucl)
      suppressWarnings(loael <- round(predictdata$dose[min(loaellogical)], digits = 4))
      
      # The following does BMDs for 1sd and for 10% higher response.
      # Note that for dichotomous outcome data, the 10% higher risk (which on the surface would 
      # seem to be analogous to what's below) is not done this way.
      # See the email from Chad Re: TiO2 data sent at 7/11/2014 at 9:30am, 
      # where the attachment gives the formulas for dichotomous data from EPA documentation.
      
      BMD1sd_response_target <- predictdata$fit[1] + stats::sd(spline$residuals)
      BMD1sdlogical <- which(predictdata$fit>BMD1sd_response_target)
      suppressWarnings(BMD1sd <- round(predictdata$dose[min(BMD1sdlogical)], digits = 4))
      BMDL1sdlogical <- which(predictdata$ucl90>BMD1sd_response_target)
      suppressWarnings(BMDL1sd <- round(predictdata$dose[min(BMDL1sdlogical)], digits = 4))
      
      BMD10perc_response_target <- predictdata$fit[1] + (0.1*predictdata$fit[1])
      BMD10perclogical <- which(predictdata$fit>BMD10perc_response_target)
      suppressWarnings(BMD10perc <- round(predictdata$dose[min(BMD10perclogical)], digits = 4))
      BMDL10perclogical <- which(predictdata$ucl90>BMD10perc_response_target)
      suppressWarnings(BMDL10perc <- round(predictdata$dose[min(BMDL10perclogical)], digits = 4))
      
      
      if (is.na(loael)) {
        warning('No ST dose can be identified in these data.')
        # AB 3/17/15 Use column "dose" as defined in predictdata, not assumed to be the name of input "dosecolumn" 
        # AB 6/5/15 all plotting now corrected and encapsulated in spline.plot
        drsmooth::spline.plot(dosecolumn, targetcolumn, data_type, k=k, data=data)
      } else {
        
        lclloaellogical <- which(predictdata$ucl90>zerodose90ucl)
        lcl90loael <- round(predictdata$dose[min(lclloaellogical)], digits = 4)
        
        uclloaellogical <- which(predictdata$lcl90>zerodose90ucl)
        suppressWarnings(ucl90loael <- round(predictdata$dose[min(uclloaellogical)], digits = 4))
        
        tdlogical <- which(predictdata$d1_lcl90>0)
        td <- round(predictdata$dose[min(tdlogical)], digits = 4)
        
        n <- nrow(x)
        predictdata$se_of_d1se <- predictdata$d1se/(n^.5)
        predictdata$d1_lcl90_u <- predictdata$d1_lcl90+(1.644854*predictdata$se_of_d1se)
        predictdata$d1_lcl90_l <- predictdata$d1_lcl90-(1.644854*predictdata$se_of_d1se)
        
        lcl90tdlogical <- which(predictdata$d1_lcl90_u>0)
        lcl90td <- round(predictdata$dose[min(lcl90tdlogical)], digits = 4)
        
        ucl90tdlogical <- which(predictdata$d1_lcl90_l>0)
        suppressWarnings(ucl90td <- round(predictdata$dose[min(ucl90tdlogical)], digits = 4))
        
        # And now the loop to perform the bootstrap.
        # AB 6/5/15 v1.9. Make this calculation driven by parameter STD_bias.
        if (STD_bias == TRUE) {
          
        # The following line initializes a vector for storing the td values from the bootstrap process to follow.
        tdsampleresults <- as.vector(rep(NA, 1000))

          for (i in 1:1000) {
            
            sampledata <- x[sample(nrow(x),replace=TRUE),]
            sampletargetvariable <- sampledata[,targetcolumn]
            sampledosevariable <- sampledata[,dosecolumn]
            
            # v1.9.0 allow user to define dimension of basis to smooth term, k
            # mgcv::s does not evaluate variables
            if (k == 4) {
              samplespline <- mgcv::gam(sampletargetvariable~s(dose, k=4), data=sampledata) }
            else {
              if (k == 3) {
                samplespline <- mgcv::gam(sampletargetvariable~s(dose, k=3), data=sampledata) }
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
            
            predictsampledata$d1 <- as.vector(d1_samplespline$dose$deriv)
            predictsampledata$d1se <- as.vector(d1_samplespline$dose$se.deriv)
            predictsampledata$d1_lcl90 <- as.vector(predictsampledata$d1-(CIfactor*predictsampledata$d1se))
            predictsampledata$d1_ucl90 <- as.vector(predictsampledata$d1+(CIfactor*predictsampledata$d1se))
            
            tdsamplelogical <- which(predictsampledata$d1_lcl90>0)
            suppressWarnings(tdsample <- predictsampledata$dose[min(tdsamplelogical)])
            
            tdsampleresults[i] <- tdsample
          }
          
          tdsampleresults_s <- sort(tdsampleresults)
          
          tdbias <- round(td-mean(tdsampleresults_s), digits = 3)
        } else {
          tdbias = NA
        }
        # AB 3/17/15 Use column "dose" as defined in predictdata, not assumed to be the name of input "dosecolumn" 
        # AB 6/5/15 all plotting now corrected and encapsulated in spline.plot
        drsmooth::spline.plot(dosecolumn, targetcolumn, data_type, k=k, data=data)
        
        neg2ll <- -2 * stats::logLik(spline)
        # R's AIC function counts the intercept as a parameter, 
        # so to be consistent we need to add 1 here because summary(spline)$edf doesn't include it.
        df <- summary(spline)$edf + 1   
        aic <- stats::AIC(spline)
        
        output <- c("STD"=td, "STD_l"=lcl90td, "STD_u"=ucl90td, "STD_bias"=tdbias, "iLOGEL"=loael, "iLOGEL_l"=lcl90loael, "iLOGEL_u"=ucl90loael)
        output2 <- c("BMD1sd"=BMD1sd, "BMDL1sd"=BMDL1sd, "BMD10perc"=BMD10perc, "BMDL10perc"=BMDL10perc)
        output3 <- c("-2LL"=neg2ll, "df"=df, "AIC"=aic)
        
        print(output)
        cat("\n")
        print(output2)
        cat("\n")
        print(output3)
        cat("\n")
        cat("NOTE: The AIC is -2LL+(2*df), where df (the model degrees of freedom, 
            or total number of model parameters) includes the intercept term.")
        
        if (write_predict == TRUE) {
          now <- Sys.time()
          file_prefix <- paste(deparse(substitute(data)), "predictions", lubridate::now(), sep = " ")
          file_name <- paste(file_prefix, ".csv", sep = "")
          utils::write.csv(predictdata, file_name, row.names=FALSE)}
        if (return_predict == TRUE) {return(predictdata)}
      }
      
    }
    else {
      stop("Please define an allowed data_type")
    }
  }
}
