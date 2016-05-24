gof.trial <-
function(model,breaks=NULL,nc=NULL)
#
# gof.trial - Computes chi-square gof test for trial models
#
# Arguments: 
#
# model  - ddf model object
# breaks - distance cut points
# nc     - number of distance classes
#
# Value:
# 
# result - lists with observed,expected, chi-square value, df and p-value
#
#  Functions used: predict(predict.trial)
#
{
    width <- model$meta.data$width 
    left <- model$meta.data$left
    xmat=model$mr$mr$data
    data=eval(model$data)
    data=data[data$observer==1&data$object %in% as.numeric(names(model$fitted)),]
    n=dim(xmat)[1]
#
#   Set up omega index; 1 - detected by secondary only, 2 - detected by both
#
    xmat$omega=rep(1,dim(xmat)[1])
    xmat$omega[xmat$timesdetected==2]=2
#
#   If number of classes for histogram intervals was not set compute a reasonable default
#
    if(is.null(nc))
    nc<-round( sqrt(min(length(xmat$distance[xmat$observer==1&xmat$detected==1]),
                   length(xmat$distance[xmat$observer==1&xmat$timesdetected==2]) )),0)
#
#   Set up default break points unless some are specified
#
    if(is.null(breaks))
       breaks <- left + ((width-left)/nc)*(0:nc)
    else
       nc=length(breaks)-1
#
#   Get predicted values for mr component
#
    xmat$detected=1
    p1=predict(model$mr,xmat,compute=TRUE,integrate=FALSE)$fitted
    p.omega=data.frame(object=rep(1:n,2),omega=c(rep(1,n),rep(2,n)),distance=rep(xmat$distance,2),prob=rep(0,2*n))
    p.omega$prob[p.omega$omega==1]=(1-p1)
    p.omega$prob[p.omega$omega==2]=p1
    expected.2=by(p.omega$prob,list(as.factor(p.omega$omega),cut(p.omega$distance,breaks,include.lowest=TRUE)),sum,na.rm=TRUE)
#
#   Get predicted values for ds component
#
    expected.1=rep(0,nc)
    for (j in 1:nc)
       expected.1[j]=sum(predict(model,compute=TRUE,int.range=matrix(c(breaks[j],breaks[j+1]),nrow=1))$fitted/model$fitted,na.rm=TRUE)   
#
#   Compute observed values of distance bins
#
    observed.count.1=table(cut(data$distance,breaks,include.lowest=TRUE))
    observed.count.2=table(as.factor(xmat$omega),cut(xmat$distance,breaks,include.lowest=TRUE))
    chisq.1=sum((observed.count.1-expected.1)^2/expected.1,na.rm=TRUE)
    chisq.2=sum((observed.count.2-expected.2)^2/expected.2,na.rm=TRUE)
    df.1=nc-1-length(model$ds$ds$par)
    if(df.1<=0)
    {   
       df.1=NA
       p.1=NA
    }
    else
       p.1=1-pchisq(chisq.1,df.1)
    df.2=nc-length(model$mr$par)
    if(df.2<=0)
    {   
       df.2=NA
       p.2=NA
    }
    else
       p.2=1-pchisq(chisq.2,df.2)
    return(list(chi1=list(observed=observed.count.1,expected=expected.1,chisq=chisq.1,p=p.1,df=df.1),
           chi2=list(observed=observed.count.2,expected=expected.2[1:2,],chisq=chisq.2,p=p.2,df=df.2),
           pooled.chi=list(chisq=chisq.1+chisq.2,df=2*nc-length(model$par)-1,p=1-pchisq(chisq.1+chisq.2,2*nc-length(model$par)-1))))
}
