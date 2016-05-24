prepareData <- function(data){
    ## Dismiss values of normalized RMSE above 1
    data <- data[data$nrmse <= 1,]
    
    ## RMSEc and difSD according to Joliff et al.
    data$nrmsec <- with(data, sqrt(nrmse^2 - nmbe^2))
    data$difSD <- with(data, sdm - sdo)
    data
}

## RMSE Circles
makeCircles <- function(data, type, cuts){

    radius <- switch(type,
                     quantile = {
                         quantile(data$nrmse, probs = cuts,
                                  na.rm = TRUE)
                     },
                     at = cuts)

    circle <- expand.grid(theta = seq(0, 2*pi,length = 100),
                          r = radius)
    circle$x <- with(circle, r * sin(theta))
    circle$y <- with(circle, r * cos(theta))
    circle
}
