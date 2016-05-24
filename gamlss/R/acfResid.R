acfResid <- function(obj=NULL, resid=NULL)
{
  resid <- if (is.null(obj))  resid
           else resid(obj)
 # if (!is.gamlss(obj)) stop("this is not a gamlss object")
  r <- resid-mean(resid)
  r2 <- r^2
  r3 <- r^3#/sd(r)^2
  r4 <- r^4#/sd(r)^4
  #nf<-layout(matrix(c(1,3,5,7,2,4,6,8),4,2))
  nf <- layout(matrix(c(0,11,12,13,14,9,1,3,5,7,10,2,4,6,8),5,3), widths = c(1.5,10,10), 
               heights = c(2, 5,5,5,5), respect = TRUE)
  #layout.show(nf) 
  op <- par(mar = c(1, 1, 1, 1), lab = c(5, 5, 7)) 
  on.exit({par(op); layout(1, widths = 1, heights = 1)})
  #op <-  par(mfrow=c(4,2))
  # mu
  acf.new <- acf(r, plot = FALSE)
  # 1
  plot(acf.new, xlim = c(2, length(acf.new$acf)), ylim = range(acf.new$acf[-1]), main=NULL)
  # 2
  pacf(r, main="residuals")
  # sigma
  acf.new2 <- acf(r2, plot = FALSE)
  # 3
  plot(acf.new2, xlim = c(2, length(acf.new2$acf)), ylim = range(acf.new2$acf[-1]), main="residuals^2")
  # 4
  pacf(r2, main="residuals^2")  
  # nu 
  acf.new3 <- acf(r3, plot = FALSE)
  # 5
  plot(acf.new3, xlim = c(2, length(acf.new3$acf)), ylim = range(acf.new3$acf[-1]), main="residuals^3")
  # 6
  pacf(r3, main="residuals^3") 
  # tau
  acf.new4 <- acf(r4, plot = FALSE)
  # 7
  plot(acf.new4, xlim = c(2, length(acf.new4$acf)), ylim = range(acf.new4$acf[-1]), main="residuals^4")
  # 8
  pacf(r4, main="residuals^4") 
  # 9
  #plot(resid(obj), type="n")
  plot.new()
  text(.5,0, labels="ACF", pos=3, cex=1.5)
  # 10
  plot.new()
  text(0.5,0, labels="PACF", pos=3, cex=1.5)
  # 11
  plot.new()
  text(.5,0.5, labels="r", pos=3, cex=1.5)
  # 12
  plot.new()
  text(.5,0.5, labels="r^2", pos=3, cex=1.5)
  # 13
  plot.new()
  text(.5,0.5, labels="r^3", pos=3, cex=1.5)
  # 14
  plot.new()
  text(.5,0.5, labels="r^4", pos=3, cex=1.5)
}
