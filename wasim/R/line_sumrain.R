line_sumrain <- function(xdata, cum_sum_rain, theMax, ...){
    while(any(cum_sum_rain>=0)){
        cum_sum_rain[cum_sum_rain>=0] <- cum_sum_rain[cum_sum_rain>=0] - theMax
    }
    cum_sum_rain <- cum_sum_rain+theMax
    cum_sum_rain[c(FALSE,diff(cum_sum_rain)<0)] <- NA
    lines(1:length(xdata),cum_sum_rain,  type="l", ...)
}
