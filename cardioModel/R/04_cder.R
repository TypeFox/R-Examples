cder.regression <- function(myData, my.nlme, delay.confidence) {
  regression_parameter <- data.frame()

  myResiduals <- data.frame(RES=residuals(my.nlme, type = c("pearson")))
  myFitted  <- data.frame(FIT=fitted(my.nlme))

  #Since i used "na.omit" for NLME have to subset only non-NA EXPOSURE
  myData <- myData[!is.na(myData$EXPOSURE),]
  myOutput <- cbind(myData, myFitted, myResiduals)

  #Calculate CDER
  original_names = names(myOutput)
  for (name in setdiff(c("COHORT", "PERIOD", "ID", "SEQUENCE", "NTAFD"), original_names)){
    myOutput[[name]] <- 0
  }
  myOutput <- myOutput[with(myOutput, order(COHORT, PERIOD, ID, SEQUENCE, NTAFD)),]
  myOutput <- myOutput[, original_names]

  myOutput$CDER <- NA
  myOutput$C2 <- NA
  myOutput$C1 <- myOutput$EXPOSURE
  myOutput$T2 <- NA
  myOutput$T1 <- myOutput$NTAFD

  my.nrow <- nrow(myOutput)
  my.t2 <- data.frame(T2=myOutput$T1[2:my.nrow])
  my.na <- data.frame(T2=0)
  my.t2 <- rbind(my.t2, my.na)
  myOutput$T2 <- my.t2$T2

  my.c2 <- data.frame(C2=myOutput$C1[2:my.nrow])
  my.na <- data.frame(C2=0)
  my.c2 <- rbind(my.c2, my.na)
  myOutput$C2 <- my.c2$C2

  for (i in 1:nrow(myOutput)){
    if (myOutput$T2[i] > myOutput$T1[i]){
      myOutput$CDER[i] <-
        (myOutput$C2[i]-myOutput$C1[i])/(myOutput$T2[i]-myOutput$T1[i])
    }
    else {
      myOutput$CDER[i] <- NA
    }
  }
  myOutput.cder <- myOutput[!is.na(myOutput$CDER),]

  #All values
  fit <- lm(myOutput.cder$RES ~ myOutput.cder$CDER)
  regression_parameter[1,"SLOPE.CDER"] <- coef(fit)["myOutput.cder$CDER"]
  regression_parameter[1,"VAR.SLOPE.CDER"] <- vcov(fit)["myOutput.cder$CDER","myOutput.cder$CDER"]

  #' upper bound of a 2-side 98% ci is used as the 99 percent upper confidence bound for the slope
  my.upper.CI <- function (parameter, variance, delay.confidence) {
    parameter+qnorm(delay.confidence)*(variance^0.5)
  }

  regression_parameter$upper.bound.99 <-
    my.upper.CI(regression_parameter$SLOPE.CDER, regression_parameter$VAR.SLOPE.CDER, delay.confidence)

  if (regression_parameter$upper.bound.99<0) {
    regression_parameter$Drug.Effect.Delay <- "yes"
  } else {
    regression_parameter$Drug.Effect.Delay <- "no"
  }
  return(regression_parameter$Drug.Effect.Delay)
}
