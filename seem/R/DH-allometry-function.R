# Parameter estimation height vs dbh allometry
# Linear estimation for various values of b4 to perform regression on b3 for each
# The b3,b4 pair with smallest residual sum-of-squares error is selected
# This is followed by nonlinear regression where the starting values are the pair selected above
# arguments data frame with D.H values and array with Hamx and Dmax values

 D.H.est <- function(D.H,D.Hmax,species.name,pdfout=F,linregout=F) {

  # Assign variables from arguments and constant breast height.
  D <- D.H$D; H <- D.H$H; Dmax<- D.Hmax[1]; Hmax <- D.Hmax[2]
  H0 <- 1.37
  # define set of b4 and allocate
  b4 <- seq(0.2,1.8,0.1); n <- length(b4)
  b3 <- array(); RSS <- array()
  # Loop through different b4 values and perform linear regression on log transformed data for each.
  for (i in 1:n){
   # Transforming variables; use negative of log to see fitted line as positive.
   Hlog <- - log(1 - ((H-H0)/(Hmax-H0))^(1/b4[i]))
   # Starting linear regression
   allom.lm <- lm(Hlog ~ D -1)
   b3[i] <- allom.lm$coeff
   # Predicted height. Note: negative sign used in exp since -log was used above.
   H.est <-DH.allom(D,param=c(Hmax,-b3[i],b4[i]))
   RSS[i] <- sum((H.est-H)^2)

   if(linregout==T){
   # Graph scatter plot and fitted line
    if(pdfout==F) win.graph()
    mat<- matrix(1:2,2,1,byrow=T)
    layout(mat,c(7,7),c(3.5,3.5),respect=TRUE)
    par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")
    plot(D,H, xlab="Diameter (DBH, cm)", ylab="Height (m)"); lines(D,H.est)
    title(paste("Linear Estimation,",species.name,",","Dmax=",Dmax,",",
         "Hmax=",Hmax,",","b3=",round(b3[i],4),",","b4=",round(b4[i],4),",","RSS=",round(RSS[i],4)),cex.main=0.7)
    plot(fitted(allom.lm), resid(allom.lm),xlab="Fitted Values", ylab="Residuals")
    abline(h=0)
   }
  }

 # Find the linear regression that resulted in the minimum residual sum-of-squares error.  
 # Select this b3, b4 pair as starting values of non linear regression.
  for(i in 1:n){
   if(identical(RSS[i],min(RSS))){startb3 <- b3[i]; startb4 <- b4[i]}
  }
 # Nonlinear Regression.
  nl.est <- nls(H ~  DH.allom(D,param=c(Hmax,-b3,b4)), data= D.H, start= list(b3=startb3, b4=startb4))
 # Graph scatter plot with fitted line and the residuals vs fitted vales.
  b3.est <- coef(nl.est)[1];  b4.est <- coef(nl.est)[2]
  H.est <- DH.allom(D,param=c(Hmax,-b3.est,b4.est))
  RSS.est <- round(sum((H.est-H)^2),4)
  b3.est <- round(b3.est,4);  b4.est <- round(b4.est,4)
  if(pdfout==F) win.graph()
   mat<- matrix(1:2,2,1,byrow=T)
   layout(mat,c(7,7),c(3.5,3.5),respect=TRUE)
   par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")
  plot(D, H, xlab="Diameter [cm]", ylab="Height [m]")
  lines(D, H.est)
  title(paste("Non Linear Estimation,", species.name,",","Dmax=",Dmax,",", 
        "Hmax=", Hmax,",","b4=", b4.est,",","b3=", b3.est,",","RSS=",RSS.est), cex.main=0.7)
  plot(H.est, resid(nl.est), xlab="Fitted Values", ylab="Residuals")
  abline(h=0)
  H.est <- round(H.est,2)
  return(list(b3=b3.est,b4=b4.est,H.est=H.est,RSS.est=RSS.est))
}

