#' Running mean of a vector
#' 
#' Function to calculate the running mean of a numeric vector
#' 
#' 
#' @param vec numeric vector
#' @param runn_mean number of vector elements to use for calculating the
#' running mean
#' @return numeric vector containing the running mean
#' @author Eike Luedeling
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' 
#' plot(runn_mean(rnorm(1000),150))
#' 
#' @export runn_mean
runn_mean<-function(vec,runn_mean)
{ww <- vec
rr <- vec
for (dd in 1:length(ww)) {
  if (dd < ceiling(runn_mean/2)) {
    rr[dd] <- mean(ww[1:(dd + floor(runn_mean/2))])
  }
  if ((dd >= ceiling(runn_mean/2)) & (dd <= length(ww) - 
                                      ceiling(runn_mean/2))) {
    rr[dd] <- mean(ww[(dd - floor(runn_mean/2)):(dd + 
                                                   floor(runn_mean/2))])
  }
  if (dd > (length(ww) - ceiling(runn_mean/2))) {
    rr[dd] <- mean(ww[(dd - floor(runn_mean/2)):length(ww)])
  }
}
return(rr)}
