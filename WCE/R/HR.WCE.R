HR.WCE <- function(x, vecnum, vecdenom, allres = FALSE){
cutoff <- ncol(x$WCEmat)
if (length(vecnum) != cutoff | length(vecdenom) != cutoff) stop("At least one of the vector provided as the numerator or denominator is not of proper length.")
if (allres == FALSE){
best <- which.min(x$info.criterion)
hr <- exp(x$WCEmat[best,]%*%vecnum)/exp(x$WCEmat[best,]%*%vecdenom)
hr} else {
hr <- exp(x$WCEmat%*%vecnum)/exp(x$WCEmat%*%vecdenom)
colnames(hr) <- 'HR'
hr}}