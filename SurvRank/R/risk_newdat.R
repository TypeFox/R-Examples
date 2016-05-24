#' Main function of SurvRank.
#'
#' Main input function for SurvRank.
#' @param dat_new a new data set that is not used for the model building but only for prediction
#' @param sel_names the variables that were selected (from riskscore_fct) (see \code{\link{CVrankSurv_fct}})
#' @param dat_old the data set used to fit the survival model
#' @param cv.out number of cross-validation folds for the prediction
#' @param c.time as defined in UnoC{survAUC} time; a positive number restricting the upper limit of the time range under consideration
#' @param detail TRUE do the prediction and Uno's C-Statistic computation for the models using 1:\code{sel_names} variables FALSE only save the statistics for the different cross validation folds
#' @param plot TRUE do a plot of the survival curves FALSE no plot
#' @param surv.tab Defaults to c(0.5). Calculates for selected features survival curves. \code{surv.tab} determines quantiles of predictions.
#' @param mcox TRUE a cox model is fitted FALSE a Cox model with ridge penalty using \code{cv.out} cross-validation folds is fitted
#' @keywords SurvRank
#' @export
#' @details details to follow
#' @return Output of the \code{risk_newdat}, basically a list containing the following elements
#' \item{\code{unocv}}{Matrix of censoring-adjusted C-statistic by Uno et al. for the different cross-validation folds and if \code{detail=T} as well for different number of variables }
#' \item{\code{unoi}}{if \code{detail=T} Vector of censoring-adjusted C-statistic by Uno et al. for the different number of variables, if \code{detail=FALSE} it correspons to \code{uno_new} }
#' \item{\code{rs}}{model prediction for the new data set}
#' \item{\code{sfit.tab}}{survfit object according to \code{surv.tab} seperation}
#' \item{\code{sfit.diff}}{surfdiff: Tests if there is a difference between two or more survival curves using the G-rho family of tests, or for a single curve against a known alternative}
#' \item{\code{model}}{model output for \code{dat_old} and using the variables given by \code{sel_names}}
#' \item{\code{uno_new}}{the censoring-adjusted C-statistic by Uno et al. using the prediction for \code{dat_new}}
#' Additionally if \code{plot} is \code{T}, the survival curves given by \code{sfit.tab} are plotted
#' @examples
#' ## Simulating a survival data set
#' N=100; p=10; n=3
#' x=data.frame(matrix(rnorm(N*p),nrow=N,p))
#' beta=rnorm(n)
#' mx=matrix(rnorm(N*n),N,n)
#' fx=mx[,seq(n)]%*%beta/3
#' hx=exp(fx)
#' ty=rexp(N,hx)
#' tcens=1-rbinom(n=N,prob=.3,size=1)
#' y=Surv(ty,tcens)
#' data=list()
#' data$x<-x; data$y<-y
#' ## CV object
#' out<-CVrankSurv_fct(data,2,3,3,fs.method="cox.rank")
#' ## The variables selected from the \code{\link{riskscore_fct}}
#' selected<-riskscore_fct(out,data,list.t="weighted")$selnames
#' ## Applying the risk_newdat function
#' x=data.frame(matrix(rnorm(N*p),nrow=N,p))
#' beta=rnorm(n)
#' mx=matrix(rnorm(N*n),N,n)
#' fx=mx[,seq(n)]%*%beta/3
#' hx=exp(fx)
#' ty=rexp(N,hx)
#' tcens=1-rbinom(n=N,prob=.3,size=1)
#' y=Surv(ty,tcens)
#' data_new=list()
#' data_new$x<-x; data_new$y<-y
#' risk<-risk_newdat(data_new,selected,data)

risk_newdat = function(dat_new,sel_names,dat_old, cv.out=10, c.time=NA, detail=NA, plot=F, surv.tab = c(0.5),mcox=T){
  res=list()
  if(is.na(c.time)){
    c.time=max(dat_new$y[,1])                  #set the time to the maximum survival time in the data
  }
  if (mcox==T){                                #fit a cox regression model using the old names and the selected features
    beta_hat = survival::coxph(dat_old$y~.,data = data.frame(subset(dat_old$x,select=sel_names)))
    pred_new = predict(beta_hat,newdata=data.frame(dat_new$x),type="lp")
  }
  if (mcox==F){                               #fit a GLM with ridge penalty using cv.out-fold cross-validation
    beta_hat = glmnet::cv.glmnet(x=as.matrix(data.frame(subset(dat_old$x,select=sel_names))),y=dat_old$y,alpha=0,family = "cox",nfolds = cv.out)
    pred_new = predict(beta_hat,newx=as.matrix(data.frame(subset(dat_new$x,select=sel_names))),type="link",s="lambda.min")
  }
  uno_new = survAUC::UnoC(dat_old$y,dat_new$y,pred_new,time = c.time)      #Compute the censoring-adjusted C-statistic by Uno et al.
  unoi= numeric(length(sel_names))
  fold.new = crossvalFolds(y = dat_new$y[,2],k=cv.out)          #creating cv.out cross-validation folds using the new data
  if(!is.na(detail)){             #if detail=T, fit a cox model with 1: sel_names variables an predict the new data
    if(detail=="full"){
      unocv = matrix(NA,nrow=length(sel_names),ncol=cv.out)
      rownames(unocv) = sel_names;colnames(unocv) = paste("kcv",1:cv.out,sep="")
      for (i in 1:length(sel_names)){
        sn = sel_names[1:i]
        betai = survival::coxph(dat_old$y ~ ., data = data.frame(subset(dat_old$x,select=sn)))
        predi = predict(betai,newdata=data.frame(dat_new$x),type="lp")
        unoi[i] = survAUC::UnoC(dat_old$y,dat_new$y,predi,c.time)          #compute the Uno statistic
        for (j in 1:cv.out){                                        #do the prediction an uno statistic computation for the cross-validation folds
          predcv = predict(betai,newdata=data.frame(dat_new$x[fold.new!=j,]),type="lp")
          unocv[i,j] = survAUC::UnoC(dat_old$y,dat_new$y[fold.new!=j,],predcv,time=c.time)
        }
      }
    }
    if(detail=="last"){
      unocv = matrix(NA,nrow=1,ncol=cv.out)
      colnames(unocv) = paste("kcv",1:cv.out,sep="")
      for (i in length(sel_names)){
        sn = sel_names[1:i]
        betai = survival::coxph(dat_old$y ~ ., data = data.frame(subset(dat_old$x,select=sn)))
        predi = predict(betai,newdata=data.frame(dat_new$x),type="lp")
        unoi = survAUC::UnoC(dat_old$y,dat_new$y,predi,c.time)
        for (j in 1:cv.out){
          predcv = predict(betai,newdata=data.frame(dat_new$x[fold.new!=j,]),type="lp")
          unocv[1,j] = survAUC::UnoC(dat_old$y,dat_new$y[fold.new!=j,],predcv,time=c.time)
        }
      }
    }
    res$unocv = unocv
    res$unoi = unoi
  }
  res$rs = pred_new
  a = pred_new
  a = pred_new
  a1 = as.integer(cut(a, breaks = c(min(a), quantile(a, probs = surv.tab),  #cut the prediction into intervals and code them according the interval they fall in
                                    max(a)), include.lowest = T, labels = c(0:length(surv.tab))))
  a2 = data.frame(y=dat_new$y[which(a1==min(a1) | a1==max(a1)),],a2=a1[which(a1==min(a1) | a1==max(a1))])
  res$sfit.tab = survival::survfit(dat_new$y ~ a1, conf.type = "log-log")           #creates survival curves
  res$sfit.ex = survival::survfit(a2$y ~ a2$a2, conf.type = "log-log")
  if(plot==T){
    plot(res$sfit.tab, conf.int = F, col = 1:length(unique(a1)),
         lwd = 2, cex.axis = 1.7, cex.lab = 1.7, xlab = "time",
         ylab = "SP",xlim=c(0,10))
    plot(res$sfit.ex, conf.int = T, col = 1:length(unique(a2)),
         lwd = 2, cex.axis = 1.7, cex.lab = 1.7, xlab = "time",
         ylab = "SP",xlim=c(0,10))
  }
  res$sfit.diff = survival::survdiff(dat_new$y ~ as.factor(a1))         #is there a difference between the survival curves?
  res$sfit.ex.diff = survival::survdiff(a2$y ~ a2$a2)
  res$model = beta_hat
  res$uno_new = uno_new
  return(res)
}

