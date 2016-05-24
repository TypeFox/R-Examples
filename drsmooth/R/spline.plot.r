#' @name spline.plot
#' @title Plot Spline
#' @usage
#' spline.plot(dosecolumn = "", targetcolumn = "", k = 4, data_type = "", data = NA)
#' @description
#' This function generates a spline model with the input dose and target response
#' columns, and plots the spline-estimated dose-response function with its upper and lower
#' 95 percent confidence bounds in green and red respectively along with the actual data.
#' Note that the confidence bounds depicted on the plot are for the dose-response function
#' itself, and not for the raw data.
#' @param dosecolumn   Name of dose column.
#' @param targetcolumn   Name of response column.
#' @param k  Dimension of the basis used to represent the smooth term.
#' @param data_type  Allowed values "continuous" or "dichotomous".
#' @param data   Input dataframe.
#' @return A plot of the spline-estimated dose-response function along with the actual data.
#' @examples
#' # Produces and plots the spline model with confidence bounds.
#' # For the same plot with key metrics, see drsmooth().
#' # For continuous outcomes:
#' data(DRdata)
#' spline.plot("dose", "MF_Log",  k = 4, data_type = "continuous", data=DRdata) 
#' 
#' # For dichotomous outcomes:
#' data(DIdata)
#' # If necessary, convert summarized dataframe into 1 row per case dataframe (see drsmooth::expand)
#' DIdata_expanded <- expand(dosecolumn = "Dose", targetcolumn = "Tumor", ncolumn = "n", data = DIdata)
#' spline.plot("Dose", "Tumor", k = 4, data_type = "dichotomous", data=DIdata_expanded)
#' @importFrom stats binomial terms
#' @export

spline.plot <- function (dosecolumn="", targetcolumn="", k = 4, data_type = "", data=NA) {
  
  targetvariable <- data[,targetcolumn]
  dose <- data[,dosecolumn]
  
  # mgcv::s does not evaluate variables
  if (k == 4) {
    if (data_type == "dichotomous") {
      spline <- mgcv::gam(targetvariable~s(dose, k=4), family=binomial, data=data)} else
      {
        spline <- mgcv::gam(targetvariable~s(dose, k=4), data=data)
      }
  } else {
    if (k == 3) {
      if (data_type == "dichotomous") {
        spline <- mgcv::gam(targetvariable~s(dose, k=3), family=binomial, data=data)} else
        {
          spline <- mgcv::gam(targetvariable~s(dose, k=3), data=data)
        }
    } else {
      stop("Dimension k of smooth term basis must be 3 or 4.")
    }
  }
  min <- min(dose)
  max <- max(dose)
  step <- (max-min)/1000
  predictdata <- data.frame(cbind(dose = seq(min, max, by = step)))
  
  alpha <- .10
  CIfactor <- stats::qt(c(alpha/2), df = spline$df.residual, lower.tail=FALSE)
  predictnew <- stats::predict(spline, predictdata, type="response", se.fit=TRUE)
  predictdata$fit <- as.vector(predictnew$fit)
  predictdata$se <- as.vector(predictnew$se)
  predictdata$lcl90 <- predictdata$fit-(CIfactor*predictdata$se)
  predictdata$ucl90 <- predictdata$fit+(CIfactor*predictdata$se)
  
  if (data_type == "dichotomous") {
    responseproportions <- tapply(targetvariable, as.factor(sort(dose)),mean)
   # if  ((max(predictdata$ucl90) < .9) & (max(targetvariable) < .9)) {
    if  ((max(predictdata$ucl90) < .9) & (max(responseproportions) < .9)) {
      graphmax <- (max(predictdata$ucl90)*.1) + max(predictdata$ucl90)}
    else {graphmax <- 1}
    
    graphics::matplot(predictdata$dose, predictdata[, c("fit","lcl90","ucl90")], 
            type="l", pch=19, lty=1, ylim=c(min(data[,targetcolumn]), graphmax),
            xlab="Dose", ylab="Predictions, CIs, and Raw Data") 
    doses <- as.vector(sort(unique(dose)))

    graphics::matpoints(doses, responseproportions, type="p", pch=19)
    
  } else {
    graphmax <- max(data[,targetcolumn])
    
    graphics::matplot(predictdata$dose, predictdata[, c("fit","lcl90","ucl90")], 
            type="l", pch=19, lty=1, ylim=c(min(data[,targetcolumn]), graphmax),
            xlab="Dose", ylab="Predictions, CIs, and Raw Data")
    graphics::matpoints(data[,dosecolumn], data[,targetcolumn], type="p", pch=19)
  }
}