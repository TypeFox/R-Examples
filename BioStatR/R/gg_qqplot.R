gg_qqplot <- function(df,var,qdist=qnorm,params=list(),qq.line=TRUE,color="red",alpha=.5) 
{
  force(params)
  y <- quantile((df[var])[!is.na(df[var])], c(0.25, 0.75))
  mf <- names(formals(qdist))
  m <- match(names(formals(qdist)), names(params), 0L)
  uparams <- params[m]
  x <- do.call("qdist",c(list(p=c(0.25, 0.75)),uparams))
  if(qq.line){
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
    }
  p <- ggplot(df, aes_string(sample=var)) + stat_qq(alpha = alpha,distribution=qdist,dparams=params)
  if(qq.line){
    p <- p + geom_abline(slope = slope, intercept = int, color=color)  
  cat(paste(c("1st quartile : ",x[1],"\n3rd quartile : ",x[2],"\nIntercept : ",int,"\nSlope : ",slope,"\n"),sep=""))
    }
  return(p)
}
