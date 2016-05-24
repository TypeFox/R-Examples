################################################################################
# TODO LIST
# TODO: Calculate per dye channel.

# NOTE: Column names used for calculations with data.table is declared
# in globals.R to avoid NOTES in R CMD CHECK.

################################################################################
# CHANGE LOG (last 20 changes)
# 30.11.2015: Added 'NB!' in the description.
# 30.11.2015: Remove rows with NA. Added 'what' parameter to 'addData'.
# 24.11.2015: Added message to print data to be used in regression.
# 28.08.2015: Added importFrom
# 20.08.2015: Now use 'sigma' instead of 'standard error of the intercept'.
# 26.06.2015: Changed to a one-sided critical t-value (alpha/2 -> alpha).
# 24.06.2015: Added some debug information.
# 18.06.2015: Flipped the signs when calculating 'Lower' and 'AT6' and added 'lower.tail = FALSE'.
#   Will give identical result but is easier to follow.
# 30.05.2015: First version.

#' @title Calculate Analytical Threshold
#'
#' @description
#' Calculate analytical thresholds estimate using linear regression.
#'
#' @details
#' Calculate the analytical threshold (AT) according to method 6 as
#' outlined in the reference. In short serial dilutions are analysed
#' and the average peak height is calculated. Linear regression or
#' Weighted linear regression with amount of DNA as the predictor for
#' the peak height is performed.
#' Method 6: A simplified version of the upper limit approach.
#' AT6 = y-intercept + t-statistic * standard error of the regression.
#' Assumes the y-intercept is not different from the mean blank signal.
#' The mean blank signal should be included in the confidence range
#' ('Lower' to 'AT6' in the resulting data frame).
#' NB! This is an indirect method to estimate AT and should be verified
#' by other methods. From the reference: A way to determine the validity
#' of this approach is based on whether the y-intercept +- (1-a)100%
#' contains the mean blank signal. If the mean blank signal is included
#' in the y-intercept band, the following relationship [i.e. AT6] can be used
#' to determine the AT. However, it should be noted that the ATs derived in
#' this manner need to be calculated for each color and for all preparations
#' (i.e., different injections, sample preparation volumes, post-PCR cleanup,
#' etc.).
#' 
#' @param data data.frame containing at least columns 'Sample.Name', 'Marker',
#'  'Allele', and 'Height'.
#' @param ref data.frame containing at least columns 'Sample.Name', 'Marker',
#' and 'Allele'.
#' @param amount data.frame containing at least columns 'Sample.Name' and 'Amount'.
#' If NULL 'data' must contain a column 'Amount'.
#' @param weighted logical to calculate weighted linear regression (weight=1/se^2).
#' @param alpha numeric [0,1] significance level for the t-statistic.
#' @param ignore.case logical to indicate if sample matching should ignore case.
#' @param debug logical to indicate if debug information should be printed.
#' 
#' @return data.frame with columns 'Amount', 'Height', 'Sd', 'Weight', 'N',
#'  'Alpha', 'Lower', 'Intercept', and 'AT6'.
#' 
#' @export
#' 
#' @importFrom stats sd qt lm
#' @importFrom utils str head tail
#' 
#' @seealso \code{\link{calculateAT6_gui}}, \code{\link{calculateAT}},
#'  \code{\link{calculateAT_gui}}, \code{\link{lm}}
#' 
#' @references
#'  J. Bregu et.al.,
#'   Analytical thresholds and sensitivity: establishing RFU thresholds for
#'   forensic DNA analysis, J. Forensic Sci. 58 (1) (2013) 120-129,
#'   ISSN 1556-4029, DOI: 10.1111/1556-4029.12008.
#' \url{http://onlinelibrary.wiley.com/doi/10.1111/1556-4029.12008/abstract}
#' 
#' 

calculateAT6 <- function(data, ref, amount=NULL, weighted=TRUE, alpha=0.05,
                        ignore.case=TRUE, debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Parameters:")
    print("data")
    print(str(data))
    print("ref")
    print(str(ref))
    print("amount")
    print(str(amount))
    print("weighted")
    print(weighted)
    print("alpha")
    print(alpha)
    print("ignore.case")
    print(ignore.case)
  }
  
  # Check data ----------------------------------------------------------------
  
  if(!"Sample.Name" %in% names(data)){
    stop("'data' must contain a column 'Sample.Name'")
  }

  if(!"Marker" %in% names(data)){
    stop("'data' must contain a column 'Marker'")
  }
  
  if(!"Allele" %in% names(data)){
    stop("'data' must contain a column 'Allele'")
  }
  
  if(!"Height" %in% names(data)){
    stop("'data' must contain a column 'Height'")
  }

  if(is.null(amount)){
    if(!"Amount" %in% names(data)){
      stop("'data' must contain a column 'Amount'")
    }
  } else {
    if(!"Sample.Name" %in% names(amount)){
      stop("'amount' must contain a column 'Sample.Name'")
    }
    if(!"Amount" %in% names(amount)){
      stop("'amount' must contain a column 'Amount'")
    }
  }
  
  # Check if slim format.  
  if(sum(grepl("Allele", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }
  
  if(sum(grepl("Height", names(data))) > 1){
    stop("'data' must be in 'slim' format",
         call. = TRUE)
  }

  if(!"Sample.Name" %in% names(ref)){
    stop("'ref' must contain a column 'Sample.Name'")
  }
  
  if(!"Marker" %in% names(ref)){
    stop("'ref' must contain a column 'Marker'")
  }
  
  if(!"Allele" %in% names(ref)){
    stop("'ref' must contain a column 'Allele'")
  }
  
  # Check if slim format.  
  if(sum(grepl("Allele", names(ref))) > 1){
    stop("'ref' must be in 'slim' format",
         call. = TRUE)
  }
    
  # Check parameters.  
  if(!is.logical(weighted)){
    stop("'weighted' must be logical",
         call. = TRUE)
  }

  if(!is.numeric(alpha) | alpha < 0 | alpha > 1){
    stop("'alpha' must be numeric [0,1]",
         call. = TRUE)
  }
  
  if(!is.logical(ignore.case)){
    stop("'ignore.case' must be logical",
         call. = TRUE)
  }
  
  # Prepare -------------------------------------------------------------------

  # Filter profiles.
  data <- filterProfile(data=data, ref=ref)

  # Add heterozygous indicator to ref.
  if(!"Heterozygous" %in% names(ref)){
    ref <- calculateHeterozygous(data=ref, debug=debug)
    message("Heterozygous indicator added to 'ref'.")
  }
  
  # Add heterozygous indicator to data.
  if(!"Heterozygous" %in% names(data)){
    data <- addData(data=data, new.data=ref, by.col="Sample.Name",
                    then.by.col="Marker", exact=FALSE, ignore.case=ignore.case,
                    what="Heterozygous", debug=debug)
    message("Heterozygous indicator added to 'data'.")
  }

  # Calculate the average peak height.
  dfHeight <- calculateHeight(data=data, na=0, add=FALSE, exclude=NULL, debug=debug)
  
  # Add amount to data.
  if(is.null(amount)){
    dfHeight <- addData(data=dfHeight, new.data=data, by.col="Sample.Name",
                        then.by.col=NULL, exact=TRUE, ignore.case=ignore.case,
                        what="Amount", debug=debug)
    message("Amount added to 'data'.")
  } else {
    dfHeight <- addData(data=dfHeight, new.data=amount, by.col="Sample.Name",
                        then.by.col=NULL, exact=TRUE, ignore.case=ignore.case,
                        what="Amount", debug=debug)
    message("Amount added to 'data'.")
  }
  
  message("Processed data to be used in regression:")
  print(dfHeight)
  
  # Convert -------------------------------------------------------------------

  # Convert to numeric.
  if(!is.numeric(dfHeight$H)){
    dfHeight$H <-  as.numeric(dfHeight$H)
    message("'H' must be numeric. 'data' converted!")
  }
  
  # Convert to numeric.
  if(!is.numeric(dfHeight$Amount)){
    dfHeight$Amount <-  as.numeric(dfHeight$Amount)
    message("'Amount' must be numeric. 'data' converted!")
  }

  # Remove NA rows.  
  if(any(is.na(dfHeight$Amount))){
    dfHeight <- dfHeight[!is.na(dfHeight$Amount), ]
    message("Removed rows where Amount=NA.")
  }

  # Analyse -------------------------------------------------------------------
  
  # Convert to data.table.
  dt <- data.table::data.table(dfHeight)
  
  # Calculate the standard deviation per amount.
  dfSd <- dt[, list(H=mean(H), Sd=sd(H)), by=Amount]
  
  # Get number of data points.
  n <- nrow(dfSd)

  if(weighted){
    
    # Calculate weights.
    weight <- 1/dfSd$Sd^2
    
    # Perform weighted linear regression.
    fit <- lm(dfSd$H ~ dfSd$Amount, weights=weight)
    
    # Extract estimates and standard error of the regression.
    coeff <- c(summary(fit)$coef[1:2], summary(fit)$sigma)
    
    if(debug){
      print("coeff:")
      print(coeff)
      print("Intercept:")
      print(coeff[1])
      print("Slope:")
      print(coeff[2])
      print("Std.Error:")
      print(coeff[3])
    }
    
    # Create result data frame.
    res <- data.frame(Amount=dfSd$Amount, Height=dfSd$H, Sd=dfSd$Sd, Weight=weight, N=n, Alpha=alpha,
                     Lower=as.numeric(coeff[1] - qt(alpha, n-1, lower.tail=FALSE) * coeff[3]),
                     Intercept=as.numeric(coeff[1]),
                     AT6=as.numeric(coeff[1] + qt(alpha, n-1, lower.tail=FALSE) * coeff[3]),
                     Std.Error=coeff[3], Slope=coeff[2])
    
  } else {
    
    # Perform linear regression.
    fit <- lm(dfSd$H ~ dfSd$Amount)

    # Extract estimates and standard error of the regression.
    coeff <- c(summary(fit)$coef[1:2], summary(fit)$sigma)
    
    if(debug){
      print("coeff:")
      print(coeff)
      print("Intercept:")
      print(coeff[1])
      print("Slope:")
      print(coeff[2])
      print("Std.Error:")
      print(coeff[3])
    }
    
    # Create result data frame.
    res <- data.frame(Amount=dfSd$Amount, Height=dfSd$H, Sd=dfSd$Sd, Weight=NA, N=n, Alpha=alpha,
                     Lower=as.numeric(coeff[1] - qt(alpha, n-1, lower.tail=FALSE) * coeff[3]),
                     Intercept=as.numeric(coeff[1]),
                     AT6=as.numeric(coeff[1] + qt(alpha, n-1, lower.tail=FALSE) * coeff[3]),
                     Std.Error=coeff[3], Slope=coeff[2])
    
    
  }

  if(debug){
    print("str(res)")
    print(str(res))
    print("head(res)")
    print(head(res))
    print("tail(res)")
    print(tail(res))
  }

  if(debug){
    print(paste("EXIT:", match.call()[[1]]))
  }
  
  # Return result.
  return(res)
  
}