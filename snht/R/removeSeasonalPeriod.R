##' Remove Seasonal Period
##' 
##' This function estimates the seasonal period of a time series via a GAM
##' model.  The data series with the seasonal period removed is then returned.
##' 
##' @usage removeSeasonalPeriod(x, period, time = 1:length(x))
##' @param x The time series to be analyzed.
##' @param period The period of the seasonality of the data.
##' @param time If not provided, then the observations are assumed to occur at
##' integer times 1, 2, ..., length(x).  Otherwise, the time vector may specify
##' when these observations occur.
##' 
##' @return Returns a vector of data with the seasonality component removed.
##' 
##' @author Josh Browning (jbrownin@@mines.edu)
##' 
##' @importFrom mgcv gam
##' @importFrom stats glm predict
##'

removeSeasonalPeriod = function(x, period, time = 1:length(x)){
    timeOfPeriod = (time - time[1]) %% period
    validObservations = sum(!is.na(x))
    if(validObservations < 3){
        stop("Can't remove seasonal period with so few observations!")
    } else if(validObservations < 15){
        warning("Very few observations available; seasonal trend removed via ",
            "glm model.")
        mod = glm(x ~ timeOfPeriod)
    } else {
        mod = mgcv::gam(x ~ s(timeOfPeriod))
    }
    x = x - predict(mod, newdata=data.frame(timeOfPeriod=timeOfPeriod))
    return(x)
}