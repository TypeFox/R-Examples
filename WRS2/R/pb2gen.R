pb2gen <- function(formula, data, est = "mom", nboot = 599){

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
  est <- get(est)  
  
  alpha=.05  
  pr=TRUE
  x<-x[!is.na(x)] # Remove any missing values in x
  y<-y[!is.na(y)] # Remove any missing values in y
  #if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  #if(pr)print("Taking bootstrap samples. Please wait.")
  datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
  bvecx<-apply(datax,1,est)
  bvecy<-apply(datay,1,est)
  bvec<-sort(bvecx-bvecy)
  low<-round((alpha/2)*nboot)+1
  up<-nboot-low
  temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
  sig.level<-2*(min(temp,1-temp))
  se<-var(bvec)
  
  result <- list(test = est(x)-est(y), conf.int = c(bvec[low],bvec[up]), p.value = sig.level, call = cl)
  class(result) <- "pb2"
  result
  
}
