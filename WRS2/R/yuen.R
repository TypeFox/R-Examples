yuen <- function(formula, data, tr = 0.2){

  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  xy <- split(model.extract(mf, "response"), mf[,2])
  faclevels <- names(xy)
  x <- xy[[1]]
  y <- xy[[2]]
  
  if (tr==0.5) warning("Comparing medians should not be done with this function!")
  

  alpha <- 0.05
  if(is.null(y)){
    if(is.matrix(x) || is.data.frame(x)){
      y=x[,2]
      x=x[,1]
    }
    if(is.list(x)){
      y=x[[2]]
      x=x[[1]]
    }
  }
  #if(tr==.5)stop("Using tr=.5 is not allowed; use a method designed for medians")
  if(tr>.25)print("Warning: with tr>.25 type I error control might be poor")
  x<-x[!is.na(x)]  # Remove any missing values in x
  y<-y[!is.na(y)]  # Remove any missing values in y
  h1<-length(x)-2*floor(tr*length(x))
  h2<-length(y)-2*floor(tr*length(y))
  q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
  q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
  df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
  crit<-qt(1-alpha/2,df)
  dif<-mean(x,tr)-mean(y,tr)
  low<-dif-crit*sqrt(q1+q2)
  up<-dif+crit*sqrt(q1+q2)
  test<-abs(dif/sqrt(q1+q2))
  yuen<-2*(1-pt(test,df))
  
  result <- list(test = test, conf.int = c(low, up), p.value = yuen, df = df, diff = dif, call = cl)
  class(result) <- "yuen"
  result
}
