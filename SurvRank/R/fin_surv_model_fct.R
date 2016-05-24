#' Building the final survRank model
#'
#' This function builds the final survRank model
#' @param f ranking approach function. One of \code{fsSurvRankConc}, \code{fsSurvRankGlmnet}, \code{fsSurvRankRf}, \code{fsSurvRankBoost}, \code{fsSurvRankCox}, \code{fsSurvRankRandCox}, \code{fsSurvRankRpart}, \code{fsSurvRankWang}, \code{fsSurvRankitBMA} or NA, no calculation
#' @param data same list used as input in \code{\link{CVrankSurv_fct}}
#' @param cv.out number of folds in outer cross validation loop (for estimation of the predictive accuracy)
#' @param nr.var Number of variables up to which stepwise selection should be carried out. Has to be smaller than n number of observations.
#' @param c.time as defined in UnoC{survAUC} time; a positive number restricting the upper limit of the time range under consideration
#' @keywords SurvRank
#' @export
#' @details details to follow
#' @return Output of the \code{riskscore_fct}, basically a list containing the following elements
#' \item{\code{ranking}}{ranking of the variables in the data set using the ranking apporach function}
#' \item{\code{used.rank}}{variables used in the survRank model according to the number of parameters to be used}
#' \item{\code{model}}{cox regression model for the selected features}
#' \item{\code{sum.model}}{summary of the fitted cox model}
fin_surv_model_fct = function(f,data,cv.out,nr.var=10,c.time=NA){
  # finales survRank model
  res=list()
  out.s=matrix(NA,nrow=cv.out,ncol=nr.var)
  fold = crossvalFolds(data$y[,2],k=cv.out)
  for(k in 1:cv.out){
    data.in = list()
    data.in$x = data$x[fold!=k,]
    data.in$y = data$y[fold!=k,]
    ranking.i = tryCatch(f(x=data.in$x,y=data.in$y), error = function(e) NA)    #ranking in the considered cross-validation loop
    if(is.na(c.time)){
      c.time = max(data.in$y[,1])
    }
    if(is.na(ranking.i)){
      ranking.i = list()
      ranking.i$rank = colnames(data.in$x)
      ranking.i$nrcoef = ncol(data.in$x)
    }
    for (j in 1:min(nr.var,ranking.i$nrcoef,na.rm=T)){                  #fit the model going for growing number of covariates
      glm.in = tryCatch(survival::coxph(y~.,data = data.frame(y=data.in$y,subset(data.in$x,select=ranking.i$rank[1:j]))), error = function(e) NA)
      pred.in = tryCatch(predict(glm.in,newdata = data.frame(data$x[fold==k,]),type="lp"),error = function(e) NA)
      out.s[k,j] = tryCatch(survAUC::UnoC(Surv.rsp = data.in$y,Surv.rsp.new=data$y[fold==k,],lpnew = pred.in,time = c.time),error = function(e) NA)
    }
  }
  out.score=apply(out.s,2,mean,na.rm=T)          #determine the number of covariates where Uno's C-statistic takes the maximum
  nr.par.in = which.max(out.score)              #set this as the number of parameters to use
  res$ranking = f(x=data$x,y=data$y)            #rank the data
  res$used.rank = res$ranking$rank[1:nr.par.in]   #take the n.par.in first variables
  res$model = survival::coxph(y~.,data=data.frame(y=data$y,subset(data$x,select=res$used.rank)))
  res$sum.model = summary(res$model)
  return(res)
}

