gains <- function(actual, predicted, groups=10, ties.method=c("max","min","first","average","random"),
      conf=c("none","normal","t","boot"), boot.reps=1000, conf.level=0.95, optimal=FALSE,
      percents=FALSE)
{

    if (length(actual) != length(predicted)) {     #Check that the vectors are of equal length
        stop("Error: The actual and predicted vectors are not of equal length")
    }

    if (!is.numeric(actual)) {
        stop("Error: The vector of actuals is not numeric")
    }

    if (!is.numeric(predicted)) {
        stop("Error: The vector of predicted values is not numeric")
    }

    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if (!is.numeric(groups) | !is.wholenumber(groups)) {
        stop("Error: The number of groups must be an integer")
    }

    ties.method <- match.arg(ties.method)
    conf <- match.arg(conf)

    if (!is.wholenumber(boot.reps)) {
        stop("Error: The number of bootstrap replications must be an integer")
    }

    if (!is.numeric(conf.level) | conf.level < 0.001 | conf.level > 0.999) {
        stop("Error: conf.level must be a number between 0.001 and 0.999")
    }

   ci.gains <- function(x, conf.type, conf.level, boot.reps)
   {
      if (conf.type=="normal") {
         n<-length(x)
         xbar<-mean(x)
         zstar<-qnorm((1-conf.level)/2,lower.tail=F)
         s<-sd(x)
         sqrt.n<-sqrt(n)
         obj<-c(xbar - zstar*(s/sqrt.n), xbar + zstar*(s/sqrt.n))
      } else if (conf.type=="t") {
         n<-length(x)
         xbar<-mean(x)
         tstar<-qt((1-conf.level)/2,df=n-1,lower.tail=F)
         s<-sd(x)
         sqrt.n<-sqrt(n)
         obj<-c(xbar - tstar*(s/sqrt.n), xbar + tstar*(s/sqrt.n))
      } else if (conf.type=="boot") {
         resamples<-lapply(1:boot.reps, function(i) sample(x,replace=T))
         boot.means<-sapply(resamples,mean)
         obj<-quantile(boot.means, probs=c((1-conf.level)/2, 1 - ((1-conf.level)/2)),names=FALSE)
      }
      return(obj)
    }

    total.n<-length(predicted)
    num.groups<-length(unique(predicted))

    if (num.groups < groups) {  #Just return all distinct scores
       warning("Warning: Fewer distinct predicted values than groups requested")
       pred.group<-rank(-1*predicted,ties.method="min")
    } else if (num.groups == groups)  {  #Just return all distinct scores
       pred.group<-rank(-1*predicted,ties.method="min")
    } else {  
       num.groups<-groups
       pred.rank<-total.n+1-rank(predicted,ties.method=ties.method)
       pred.group<-ceiling(pred.rank/(total.n/groups))    }

    total.resp<-sum(actual)
    total.resp.rate<-total.resp/total.n

    n<-as.vector(tapply(predicted,pred.group,FUN=length))
    cume.n<-cumsum(n)
    depth<-100*round(cume.n/total.n,2)
    resp<-as.vector(tapply(actual,pred.group,FUN=sum))
    cume.resp<-cumsum(resp)    
    resp.rate<-resp/n
    cume.resp.rate<-cume.resp/cume.n
    cume.pct.of.total<-cume.resp/total.resp
    lift<-round(100*(resp.rate/total.resp.rate),0)
    cume.lift<-round(100*(cume.resp.rate/total.resp.rate),0)
    mean.prediction<-as.vector(tapply(predicted,pred.group,FUN=mean))
    min.prediction<-as.vector(tapply(predicted,pred.group,FUN=min))
    max.prediction<-as.vector(tapply(predicted,pred.group,FUN=max))

    if (conf != "none") {
       conf.int<-tapply(actual,pred.group,
         FUN="ci.gains",conf.type=conf,conf.level=conf.level,boot.reps=boot.reps)
    }

    if (optimal==TRUE) {
       opt.rank<-total.n+1-rank(actual,ties.method="random")
       opt.group<-cut(opt.rank,c(0,cume.n),labels=F)
       opt.resp<-as.vector(tapply(actual,opt.group,FUN=sum))
       opt.cume.resp<-cumsum(opt.resp)
       opt.resp.rate<-opt.resp/n
       opt.cume.resp.rate<-opt.cume.resp/cume.n
       opt.lift<-round(100*(opt.resp.rate/total.resp.rate),0)
       opt.cume.lift<-round(100*(opt.cume.resp.rate/total.resp.rate),0)
    }

    if (conf != "none" & optimal==TRUE) {
    obj<-list(depth=depth,obs=n,cume.obs=cume.n,mean.resp=resp.rate,
              cume.mean.resp=cume.resp.rate, cume.pct.of.total=cume.pct.of.total,
              lift=lift,cume.lift=cume.lift,mean.prediction=mean.prediction,
              min.prediction=min.prediction, max.prediction=max.prediction,
              conf=conf,optimal=optimal,
              num.groups=num.groups, percents=percents,
              conf.lower=unlist(conf.int,use.names=F)[seq(from=1,to=2*num.groups,by=2)],
              conf.upper=unlist(conf.int,use.names=F)[seq(from=2,to=2*num.groups,by=2)],
              opt.lift=opt.lift, opt.cume.lift=opt.cume.lift)
    } else if (conf != "none" & optimal!=TRUE) {
    obj<-list(depth=depth,obs=n,cume.obs=cume.n,mean.resp=resp.rate,
              cume.mean.resp=cume.resp.rate, cume.pct.of.total=cume.pct.of.total,
              lift=lift,cume.lift=cume.lift,mean.prediction=mean.prediction,
              min.prediction=min.prediction, max.prediction=max.prediction,
              conf=conf,optimal=optimal,
              num.groups=num.groups, percents=percents,
              conf.lower=unlist(conf.int,use.names=F)[seq(from=1,to=2*num.groups,by=2)],
              conf.upper=unlist(conf.int,use.names=F)[seq(from=2,to=2*num.groups,by=2)])
    } else if (conf == "none" & optimal==TRUE) {
    obj<-list(depth=depth,obs=n,cume.obs=cume.n,mean.resp=resp.rate,
              cume.mean.resp=cume.resp.rate, cume.pct.of.total=cume.pct.of.total,
              lift=lift,cume.lift=cume.lift,mean.prediction=mean.prediction,
              min.prediction=min.prediction, max.prediction=max.prediction,
              conf=conf,optimal=optimal,
              num.groups=num.groups, percents=percents,
              opt.lift=opt.lift, opt.cume.lift=opt.cume.lift)
    } else if (conf == "none" & optimal!=TRUE) {
    obj<-list(depth=depth,obs=n,cume.obs=cume.n,mean.resp=resp.rate,
              cume.mean.resp=cume.resp.rate, cume.pct.of.total=cume.pct.of.total,
              lift=lift,cume.lift=cume.lift,mean.prediction=mean.prediction,
              min.prediction=min.prediction, max.prediction=max.prediction,
              conf=conf,optimal=optimal,
              num.groups=num.groups, percents=percents)
    }

    class(obj) <- "gains"
    return(obj)
}

