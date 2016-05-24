plot.zipath=function(x, xvar=c("norm","lambda"),type=c("count", "zero"), label=FALSE,...){
  xvar=match.arg(xvar)
  type=match.arg(type)
 if(type=="count"){
 b <- x$coefficients$count[-1,]
 df <- apply(abs(b)>0, 2, sum)
 plotCoef(b,lambda=x$lambda.count,df=df,label=label,xvar=xvar,...)
} else{
 b <- x$coefficients$zero[-1,]
 df <- apply(abs(b)>0, 2, sum)
 plotCoef(b,lambda=x$lambda.zero,df=df,label=label,xvar=xvar,...)
}
}
