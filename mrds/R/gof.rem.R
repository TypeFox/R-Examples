# Compute chi-square gof test for rem models
#
# model ddf model object
# breaks distance cut points
# nc number of distance classes
#
# return list with chi-square value, df and p-value
# documented in ?ddf.gof
gof.rem <- function(model,breaks=NULL,nc=NULL){
  width <- model$meta.data$width
  left <- model$meta.data$left
  xmat <- model$mr$mr$data
  n <- nrow(xmat)

  # Set up omega index
  #   1 - detected by primary/secondary
  #   2 - detected by secondary only
  xmat$omega <- rep(2,dim(xmat)[1])
  xmat$omega[xmat$detected==0] <- 1

  # If number of classes for histogram intervals was not set
  # compute a reasonable default
  if(is.null(nc)){
    nc<-round( sqrt(min(length(xmat$distance[xmat$detected==1]),
                 length(xmat$distance[xmat$detected==0]) )),0)
  }

  # Set up default break points - need to allow user-defined values
  if(is.null(breaks)){
    breaks <- left + ((width-left)/nc)*(0:nc)
  }else{
    nc <- length(breaks)-1
  }

  # Get predicted values for mr component
  predict.list <- predict(model$mr)
  p1 <- predict.list$p1
  p2 <- predict.list$p2
  p.omega <- data.frame(object=rep(1:n,2),
                        omega=c(rep(1,n),rep(2,n)),
                        distance=rep(xmat$distance,2),
                        prob=rep(0,2*n))
  p.omega$prob[p.omega$omega==2] <- p1*1/(p1+p2-p1*p2)
  p.omega$prob[p.omega$omega==1] <- p2*(1-p1)/(p1+p2-p1*p2)
  expected.2 <- by(p.omega$prob,list(as.factor(p.omega$omega),
                                     cut(p.omega$distance,breaks,
                                         include.lowest=TRUE)),
                   sum,na.rm=TRUE)

  # Get predicted values for ds component
  expected.1 <- rep(0,nc)
  for(j in 1:nc){
    expected.1[j] <- sum(predict(model,compute=TRUE,int.range=matrix(c(breaks[j],breaks[j+1]),nrow=1))$fitted/model$fitted,na.rm=TRUE)
  }

  # Compute observed values of distance bins
  observed.count.1 <- table(cut(xmat$distance,breaks,include.lowest=TRUE))
  observed.count.2 <- table(as.factor(xmat$omega),
                            cut(xmat$distance,breaks,include.lowest=TRUE))
  chisq.1 <- sum((observed.count.1-expected.1)^2/expected.1,na.rm=TRUE)
  chisq.2 <- sum((observed.count.2-expected.2)^2/expected.2,na.rm=TRUE)
  df.1 <- nc-1-length(model$ds$ds$par)
  if(df.1<=0){
    df.1 <- NA
    p.1 <- NA
  }else{
    p.1 <- 1-pchisq(chisq.1,df.1)
  }

  df.2 <- 2*nc-length(model$mr$par)
  if(df.2<=0){
    df.2 <- NA
    p.2 <- NA
  }else{
    p.2 <- 1-pchisq(chisq.2,df.2)
  }

  return(list(chi1=list(observed=observed.count.1,
                        expected=expected.1,
                        chisq=chisq.1,
                        p=p.1,
                        df=df.1),
         chi2=list(observed=observed.count.2,
                   expected=expected.2[1:2,],
                   chisq=chisq.2,
                   p=p.2,
                   df=df.2),
         pooled.chi=list(chisq=chisq.1+chisq.2,
                         df=2*nc-length(model$par)-1,
                         p=1-pchisq(chisq.1+chisq.2,
                                    2*nc-length(model$par)-1))))
}
