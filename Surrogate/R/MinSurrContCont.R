MinSurrContCont <- function(T0T0, T1T1, Delta, T0T1=seq(from=0, to=1, by=.01)){

sigma_delta_T <- T0T0 + T1T1 - (2 * sqrt(T0T0*T1T1) * T0T1)
rho2_min <- 1 - (Delta/sigma_delta_T)   
T0T1 <- T0T1[rho2_min >=-1 & rho2_min <=1]
sigma_delta_T <- sigma_delta_T[rho2_min >=-1 & rho2_min <=1] 
rho2_min <- rho2_min[rho2_min >=-1 & rho2_min <=1] 

fit <- 
  list(T0T1=T0T1, Sigma.Delta.T=sigma_delta_T, Rho2.Min = rho2_min, Call=match.call())   

class(fit) <- "MinSurrContCont"
fit
}


plot.MinSurrContCont <- function(x, main, col, Type="Percent", Labels=FALSE, Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), ...) {

  
  Object <- x 
  if (missing(main)) {main = " "}
  if (missing(col)) {col=8}
  
  dat <- data.frame(cbind(Object$T0T1, Object$Rho2.Min))
  colnames(dat) <- c("T0T1", "Rho2.Min")
  dev.new()
  par=Par
  
  if (Type=="Freq"){
  
    h <- hist(dat$Rho2.Min, ...)
    h$density <- h$counts/sum(h$counts)
    cumulMidPoint <- ecdf(x=dat$Rho2.Min)(h$mids)
    labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
    
    if (Labels==FALSE){
      plot(h, freq=T, xlab=expression(rho[min]^2), ylab="Frequency", col=col, main=main)}
    if (Labels==TRUE){
      plot(h, freq=T, xlab=expression(rho[min]^2), ylab="Frequency", col=col, main=main, labels=labs)}
    }
  
   if (Type=="Percent"){
    
    h <- hist(dat$Rho2.Min, ...)
    h$density <- h$counts/sum(h$counts)
    cumulMidPoint <- ecdf(x=dat$Rho2.Min)(h$mids)
    labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
    
    if (Labels==FALSE){
      plot(h, freq=F, xlab=expression(rho[min]^2), ylab="Percentage", col=col, main=main)}
    if (Labels==TRUE){
      plot(h, freq=F, xlab=expression(rho[min]^2), ylab="Percentage", col=col, main=main, labels=labs)}
      }
  
  if (Type=="CumPerc"){
    h <- hist(dat$Rho2.Min, breaks=length(dat$Rho2.Min), ...)
    h$density <- h$counts/sum(h$counts)
    cumulative <- cumsum(h$density)
    plot(x=h$mids, y=cumulative, xlab=expression(rho[min]^2), ylab="Cumulative percentage", col=0, main=main)
    lines(x=h$mids, y=cumulative)
  }


}


summary.MinSurrContCont <- function(object, ..., Object){

  if (missing(Object)){Object <- object} 
  cat("\nFunction call:\n\n")
  print(Object$Call)
  cat("\n\n\n# Rho2.Min results summary (Inf values are excluded)")
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")  
  cat("Mean (SD) Rho^2_min: ", format(round(mean(Object$Rho2.Min[which(is.finite(Object$Rho2.Min))]), 4), nsmall = 4), " (", format(round(sd(Object$Rho2.Min[which(is.finite(Object$Rho2.Min))]), 4), nsmall = 4), ")", 
      "  [min: ", format(round(min(Object$Rho2.Min[which(is.finite(Object$Rho2.Min))]), 4), nsmall = 4), "; max: ",  format(round(max(Object$Rho2.Min[which(is.finite(Object$Rho2.Min))]), 4), nsmall = 4), "]", sep="")
  cat("\n\nQuantiles of the Rho2.Min distribution: \n\n")
  quant <- quantile(Object$Rho2.Min, probs = c(.05, .10, .20, .50, .80, .90, .95))
  print(quant)
}

