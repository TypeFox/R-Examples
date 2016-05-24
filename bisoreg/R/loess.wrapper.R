## `loess.wrapper` <- function(x,y,span.vals=seq(0.25,1.00,by=0.05),x0=x,folds=5){
`loess.wrapper` <- function(x,y,span.vals=seq(0.25,1.00,by=0.05),folds=5){
  #require(bootstrap)
  mae <- numeric(length(span.vals))
  theta.fit <- function(x,y,span) loess(y~x,span=span)
  theta.predict <- function(fit,x0) predict(fit,newdata=x0)
  ii=0
  for(span in span.vals){
    ii <- ii+1
    y.cv <- crossval(x,y,theta.fit,theta.predict,span=span,ngroup=folds)$cv.fit
    fltr <- !is.na(y.cv)
    mae[ii] <- mean(abs(y[fltr]-y.cv[fltr]))
    }
  span <- span.vals[which.min(mae)]
  out <- loess(y~x,span=span)
  return(out)
  }

