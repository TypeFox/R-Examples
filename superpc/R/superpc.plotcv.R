superpc.plotcv=
function (object, cv.type=c("full","preval"),smooth = TRUE, smooth.df = 10, call.win.metafile=FALSE, ...)
{

cv.type=match.arg(cv.type)

if(cv.type=="full"){
 scor=object$scor
 smooth=FALSE
 }
else{
    scor = object$scor.preval
}
k=nrow(scor)
    if (smooth) {
        for (j in 1:nrow(scor)) {
         if(is.null(smooth.df)){
               om=!is.na(scor[j, ])
               junk=smooth.spline(object$th[om], scor[j,om ])
               scor[j,om] = predict(junk,object$th[om])$y
              }
            if(!is.null(smooth.df)){
               om=!is.na(scor[j, ])
             junk=smooth.spline(object$th[om], scor[j,om ], df=smooth.df)
            scor[j,om] =predict(junk,object$th[om])$y
             }
        }
    }

if(object$type=="survival"){
if(cv.type=="full"){    ymax = max(object$scor.upper[!is.na(object$scor.upper)], qchisq(0.95, nrow(scor)))}
if(cv.type=="preval"){    ymax = max(scor[!is.na(scor)], qchisq(0.95, nrow(scor)))}
}

if(object$type=="regression"){
  # df of denom for f is average sample size in validation fold

n.mean=0
for(i in 1:object$n.fold){
   n.mean=n.mean+length(object$folds[[i]])/object$n.fold
}
   denom.df=n.mean-1-nrow(scor)
if(cv.type=="full"){    ymax = max(object$scor.upper[!is.na(object$scor.upper)], qf(0.95, nrow(scor), denom.df))}
if(cv.type=="preval"){    ymax = max(scor[!is.na(scor)], qf(0.95, nrow(scor), denom.df))}
}

if(call.win.metafile){win.metafile()}

# if(object$type=="survival"){ ylab="Likelihood ratio test statistic"}
#if(object$type=="regression"){ ylab="F statistic"}

ylab="Likelihood ratio test statistic"

    matplot(object$th, t(scor), xlab = "Threshold", ylab = ylab, ylim = c(0, ymax), lty=rep(1,k))
    matlines(object$th, t(scor), lty=rep(1,k), ...)


    for (j in 1:k) {
      if(object$type=="survival"){  abline(h = qchisq(0.95, j), lty = 2, col = j)}
       if(object$type=="regression"){
         # df of denom for f is average sample size in validation fold
         abline(h = qf(0.95, j, denom.df), lty = 2, col = j)
       }
if(cv.type=="full"){
delta=((-1)^j)*diff(object$th)[1]/4
error.bars(object$th+delta*(j>1),t(object$scor.lower[j,]),
 t(object$scor.upper[j,]),lty=2, col=j)
}
    }


if(call.win.metafile){dev.off()}
return(TRUE)

}

error.bars <-function(x, upper, lower, width = 0.005, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}
jitter<-
function(x)
{
        x + 0.03 * abs(x) * sign(rnorm(length(x)))
}


