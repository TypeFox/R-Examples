yuenbt <- function(formula, data, tr = 0.2, nboot = 599){
#

  side <- TRUE
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
  alpha<-.05
  nullval<-0
  pr<-TRUE
  plotit<-FALSE
  op<-1
  side<-as.logical(side)
  p.value<-NA
  yuenbt<-vector(mode="numeric",length=2)
  #if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  x<-x[!is.na(x)]  # Remove missing values in x
  y<-y[!is.na(y)]  # Remove missing values in y
  xcen<-x-mean(x,tr)
  ycen<-y-mean(y,tr)
  if(!side){
    if(pr)warning("p-value computed only when side = TRUE")
  }
  test<-(mean(x,tr)-mean(y,tr))/sqrt(trimse(x,tr=tr)^2+trimse(y,tr=tr)^2)
  datax<-matrix(sample(xcen,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  datay<-matrix(sample(ycen,size=length(y)*nboot,replace=TRUE),nrow=nboot)
  top<-apply(datax,1,mean,tr)-apply(datay,1,mean,tr)
  botx<-apply(datax,1,trimse,tr)
  boty<-apply(datay,1,trimse,tr)
  tval<-top/sqrt(botx^2+boty^2)
 
  if(side)tval<-abs(tval)
  tval<-sort(tval)
  icrit<-floor((1-alpha)*nboot+.5)
  ibot<-floor(alpha*nboot/2+.5)
  itop<-floor((1-alpha/2)*nboot+.5)
  se<-sqrt((trimse(x,tr))^2+(trimse(y,tr))^2)
  yuenbt[1]<-mean(x,tr)-mean(y,tr)-tval[itop]*se
  yuenbt[2]<-mean(x,tr)-mean(y,tr)-tval[ibot]*se
  if(side){
    yuenbt[1]<-mean(x,tr)-mean(y,tr)-tval[icrit]*se
    yuenbt[2]<-mean(x,tr)-mean(y,tr)+tval[icrit]*se
    p.value<-(sum(abs(test)<=abs(tval)))/nboot
  }
  mdiff <- mean(x,tr)-mean(y,tr)
  
  result <- list(test = test, conf.int = yuenbt, p.value = p.value, df = NA, diff = mdiff, call = cl)
  class(result) <- "yuen"
  result
  
}
