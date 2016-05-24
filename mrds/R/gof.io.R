gof.io <-
function(model,breaks=NULL,nc=NULL)
#
# gof.io - Computes chi-square gof test for io models
#
# Arguments: 
#
# model  - ddf model object
# breaks - distance cut points
# nc     - number of distance classes
#
# Value:
# 
# result - list with chi-square value, df and p-value
#
#  Functions used: predict(predict.io)
#
{
    width <- model$meta.data$width
    left <- model$meta.data$left
    xmat=model$mr$mr$data
    n=dim(xmat)[1]/2
#
#   Set up omega index; 1 - detected by primary only, 2 - detected by secondary only, 3 - detected by both
#
    xmat=xmat[xmat$observer==1,]
    xmat$omega=rep(1,dim(xmat)[1])
    xmat$omega[xmat$timesdetected==2]=3
    xmat$omega[xmat$timesdetected==1&xmat$detected==0]=2
#
#   If number of classes for histogram intervals was not set compute a reasonable default
#
    if(is.null(nc))
    nc<-round( sqrt(min(length(xmat$distance[xmat$observer==1&xmat$detected==1]),
                   length(xmat$distance[xmat$observer==1&xmat$timesdetected==2]) )),0)
#
#   Set up default break points - need to allow user-defined values
#
    if(is.null(breaks))
      breaks <- left + ((width-left)/nc)*(0:nc)
    else
       nc=length(breaks)-1
#
#   Get predicted values for mr component
#
    predict.list=predict(model$mr)
    p1=predict.list$p1
    p2=predict.list$p2
    p.omega=data.frame(object=rep(1:n,3),omega=c(rep(1,n),rep(2,n),rep(3,n)),distance=rep(xmat$distance,3),prob=rep(0,3*n))
    p.omega$prob[p.omega$omega==1]=p1*(1-p2)/(p1+p2-p1*p2)
    p.omega$prob[p.omega$omega==2]=p2*(1-p1)/(p1+p2-p1*p2)
    p.omega$prob[p.omega$omega==3]=p1*p2/(p1+p2-p1*p2)
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
    observed.count.1=table(cut(xmat$distance,breaks,include.lowest=TRUE))
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
    df.2=2*nc-length(model$mr$par)
    if(df.2<=0)
    {   
       df.2=NA
       p.2=NA
    }
    else
       p.2=1-pchisq(chisq.2,df.2)
    return(list(chi1=list(observed=observed.count.1,expected=expected.1,chisq=chisq.1,p=p.1,df=df.1),
           chi2=list(observed=observed.count.2,expected=expected.2[1:3,],chisq=chisq.2,p=p.2,df=df.2),
           pooled.chi=list(chisq=chisq.1+chisq.2,df=3*nc-length(model$par)-1,p=1-pchisq(chisq.1+chisq.2,3*nc-length(model$par)-1))))
}
