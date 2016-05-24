ICBayes.formula <-
function (formula, data, ...) 
  {
    # 'data' must be data.frame
    mf<-model.frame(formula=formula,data=data)  
    y<-model.response(mf)
    LR<-matrix(SurvtoLR(y),ncol=2)
    x<-model.matrix(attr(mf,'terms'),data=mf)  # get design matrix
    xcov<-x[,-1]  # delete intercept column
    est<-ICBayes.default(L=LR[,1],R=LR[,2],xcov=xcov,...)
    #est<-do.call('ICBayes',c(list(L=LR[,1],R=LR[,2],xcov=xcov),list(...)))
    est$call <- match.call()
    est$formula <- formula
    est
}
