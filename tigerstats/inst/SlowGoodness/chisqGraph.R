# Modified from pchisqGC.R in package tigerstats

chisqGraph <- function(bound,region="above",df=NA,xlab="chi_square_statistic",graph=FALSE) {
  if (!is.numeric(bound)) stop("Specify a numerical boundary")
  if (bound < 0)  stop("The chi-square statistic must be at least 0")
  if (is.na(df)) stop("Specify the degrees of freedom using the argument df")
  if (!(region %in% c("below","above"))) stop("Specify either \"region=\"below\" or
\"region=\"above\"")
  if (df < 0) stop("Degrees of freedom must be positive")

  if (region=="below")  {
    area <- pchisq(bound,df=df)
    if (graph && df==1) warning("No graph produced for region below when df=1")
    if (graph) {
    bound <- round(bound,2)
    upper <- max(qchisq(.9999,df=df),bound+1)
    lower <- 0
    curve(dchisq(x,df=df),from=lower,to=upper,ylab="density",axes=FALSE,n=501,xlab=xlab,
          main=paste("Chi-Square Curve, df = ",df,"\nShaded Area = ",round(area,4)))
    axis(1,at=c(lower,bound,upper),labels=c(as.character(0),as.character(bound),""))
    axis(2)
   x.coords <- c(lower,seq(lower,bound,length.out=301),bound)
   y.coords <- c(0,dchisq(seq(lower,bound,length.out=301),df=df),0)
   polygon(x.coords,y.coords,col="lightblue",cex=2)
    }
  }

  if (region=="above")  {
    area <- pchisq(bound,df=df,lower.tail=FALSE)
    if (graph) {
    bound <- round(bound,2)
    upper <- max(qchisq(.9999,df=df),bound+1)
    lower <- 0
    curve(dchisq(x,df=df),from=lower,to=upper,ylab="density",axes=FALSE,n=501,xlab=xlab,
          main=paste("Chi-Square Curve, df = ",df,"\nShaded Area = ",round(area,4)))
    axis(1,at=c(lower,bound,upper),labels=c(as.character(0),as.character(bound),""))
    axis(2)
    x.coords <- c(bound,seq(bound,upper,length.out=301),upper)
    y.coords <- c(0,dchisq(seq(bound,upper,length.out=301),df=df),0)
    polygon(x.coords,y.coords,col="lightblue",cex=2)
    }
  }
  
 

}#end of chisqGraph