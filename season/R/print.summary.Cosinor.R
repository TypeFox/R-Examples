# Oct 2011
print.summary.Cosinor <- function(x, ...){
  ## report results
  if (class(x)!="summary.Cosinor"){
    stop("Object must be of class 'summary.Cosinor'")
  } 
# fix the digits, October 2011
 if(x$text!=TRUE){
    x$phase=round(x$phase,x$digits)
    x$lphase=round(x$lphase,x$digits)
 }

  cat('Cosinor test\n')
  cat('Number of observations =',x$n,'\n')
  cat('Amplitude =',round(x$amp,x$digits),x$amp.scale,'\n')
  cat('Phase:',x$phase,'\n')
  cat('Low point:',x$lphase,'\n')
  if(x$type=='hourly'){
    cat('Significant circadian pattern based on adjusted significance level of',
        x$alpha/2,' = ',x$significant,'\n', ...)
  }
  if(x$type!='hourly'){
    cat('Significant seasonality based on adjusted significance level of',
        x$alpha/2,' = ',x$significant,'\n', ...)
  }
}
