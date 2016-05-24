runmean <- function(x, y, fr = 1, tr = 0.2){
  #
  # running mean using interval method
  #
  # fr controls amount of smoothing
  # tr is the amount of trimming
  #
  # Missing values are automatically removed.
  #
  eout <- FALSE
  xout <- FALSE
  pyhat <- TRUE
  #if(eout && xout) xout<-F
  temp <- cbind(x,y)
  temp <- elimna(temp) # Eliminate any rows with missing values
  #if(eout){
  #  flag<-outfun(temp,plotit=FALSE)$keep
  #  temp<-temp[flag,]
  #}
  #if(xout){
  #  flag<-outfun(x,plotit=FALSE)$keep
  #  temp<-temp[flag,]
  #}
  x <- temp[,1]
  y <- temp[,2]
  #pyhat <- as.logical(pyhat)
  rmd <- c(1:length(x))
  for(i in 1:length(x)) rmd[i] <- mean(y[near(x,x[i],fr)],tr)
  return(rmd)
}