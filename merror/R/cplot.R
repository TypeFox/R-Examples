"cplot" <- function(df,i,j,leg.loc="topleft",regress=FALSE,lw=1,t.size=1,alpha.beta.sigma=NULL)
{
#  print(head(df))
  
  del <- dim(df)[2] + 1
#  print(del)

  if(is.null(alpha.beta.sigma)) {
    a <- ncb.od(df)$sigma.table$alpha.ncb[-del]
    b <- ncb.od(df)$sigma.table$beta[-del]
    s <- ncb.od(df)$sigma.table$sigma[-del]
  }
  else {
    a <- alpha.beta.sigma[1,]
    b <- alpha.beta.sigma[2,]
    s <- alpha.beta.sigma[3,]
  }  
#  print(a)
#  print(b)
#  print(s)

  A21 <- a[i]-a[j]*b[i]/b[j]
  B21 <- b[i]/b[j]
#  print(A21)
#  print(B21)
  
  A12 <- a[j]-a[i]*b[j]/b[i]
  B12 <- b[j]/b[i]
#  print(A12)
#  print(B12)

  plot(df[,i]~df[,j],data=df,xlab=names(df)[j],ylab=names(df)[i],las=1)
  abline(0,1,lty=2,col="green",lwd=lw)
  abline(A21,B21,col="red",lwd=lw)


  title(main=paste("Calibration Curve: ",
  names(df)[i]," = ",format(round(A21,3))," + ",format(round(B21,3))," ",names(df)[j],
  " and ",
  names(df)[j]," = ",format(round(A12,3))," + ",format(round(B12,3))," ",names(df)[i],
  "\nScale Adjusted Imprecision SDs - ",names(df)[i],": ",format(round(s[i]/b[i],3)),
  " - ",names(df)[j],": ",format(round(s[j]/b[j],3)),
  sep=""),cex.main=t.size)

  if(regress) {
    lm.21 <- lm(df[,i]~df[,j])
    print(summary(lm.21))

    abline(lm.21,col="blue",lty=2,lwd=lw)

    lm.12 <- lm(df[,j]~df[,i])
    print(summary(lm.12))
    cf <- coef(lm.12)

    abline(-cf[1]/cf[2],1/cf[2],col="black",lty=2,lwd=lw)
    
    legend(leg.loc,c("No Bias Model","Calibration Curve",
      paste(names(df)[i]," regressed on ",names(df)[j],sep=""),
      paste(names(df)[j]," regressed on ",names(df)[i],sep="")),
      lty=c(2,1,2,2),col=c("green","red","blue","black"),lwd=lw)
      
    title(sub="Regression lines are not calibration lines and are for comparison only.",cex.sub=0.75)
  }
  
  else {
    legend(leg.loc,c("No Bias Model","Calibration Curve"),
      lty=c(2,1),col=c("green","red"),lwd=lw)
  }
   
}
