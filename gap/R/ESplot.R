ESplot <- function(ESdat,SE=TRUE,logscale=TRUE,alpha=0.05,xlim=c(-2,8),v=1,...)
{
   n <- dim(ESdat)[1]
   models <- ESdat[,1]
   ES <- ESdat[,2]
   if (SE)
   {
      z <- abs(qnorm(alpha/2))
      if (logscale)
      {
         SElogES <- ESdat[,3]
         LCL <- exp(log(ES)-z*SElogES)
         UCL <- exp(log(ES)+z*SElogES)
      } else {
         SEplain <- ESdat[,3]
         LCL <- ES-z*SEplain
         UCL <- ES+z*SEplain
      }
   } else {
      LCL <- ESdat[,3]
      UCL <- ESdat[,4]
   }
   x <- xlim[1]:xlim[2]
   plot(x,x,type="n",xlab="",ylab="",axes=FALSE)
   for(i in 1:n)
   {
       points(ES[i],i,pch=22,bg="black")
       segments(LCL[i],i,UCL[i],i)
   }
   axis(1)
   par(las=1)
   abline(v=v)
   axis(2,labels=models,at=1:n,lty="blank",hadj=0)
}
