#' Wrapper function using f as ranking method
#'
#' This wrapper function passes the ranking method to further functions and computes the outer survival AUCs. Setting up for parallel computing of different folds via foreach.
#' @param input see details in \code{\link{riskscore_fct}}
#' @param ... other arguments, not used now
#' @keywords internal
#' @export

fixRank_fct = function(data,f,fix.var,cv.out,t.times,ncl=2,c.time=10,...){
  pb = utils::txtProgressBar(style=3,min=1,max=t.times*cv.out)
  res = ranking = data.out  = list()
  auc.out = matrix(NA,nrow=cv.out,ncol=t.times)
  i=NULL
  # bei LOO
  if (cv.out==length(data$y[,2])){
    fold.outer = c(1:length(data$y[,2]))
    t.times=1
  }
  # determine prediction accuracy
  for (t in 1:t.times){
    utils::setTxtProgressBar(pb,t*cv.out)
    # nicht outer LOO
    if(cv.out!=length(data$y[,2])){
      fold.outer = crossvalFolds(data$y[,2],k=cv.out)
    }
    cl=parallel::makeCluster(ncl,type="PSOCK")
    doParallel::registerDoParallel(cl)
    output = foreach::foreach(i = 1:cv.out, .packages=c('SurvRank'),.errorhandling = "pass",.verbose = F) %dopar% {
      #for(i in 1:cv.out){
      l=(t-1)*cv.out+i
      data.out$x = data$x[fold.outer!=i,]
      data.out$y = data$y[fold.outer!=i,]
      ranking[[l]] = tryCatch(f(x=data.out$x,y=data.out$y), error = function(e) NA)
      glm.in = tryCatch(survival::coxph(y~.,data = data.frame(y=data.out$y,subset(data.out$x,select=ranking[[l]]$rank[1:fix.var]))), error = function(e) NA)
      pred.in = tryCatch(predict(glm.in,newdata = data.frame(data$x[fold.outer==i,]),type="lp",reference="sample"),error = function(e) NA)
      #auc.out[i,t] = tryCatch(UnoC(Surv.rsp = data.out$y,Surv.rsp.new=data$y[fold.outer==i,],lpnew = pred.in,time = c.time),error = function(e) NA)
      auc.out = tryCatch(survAUC::UnoC(Surv.rsp = data.out$y,Surv.rsp.new=data$y[fold.outer==i,],lpnew = pred.in,time = c.time),error = function(e) NA)
      list(auc.out)
    }
    parallel::stopCluster(cl)
    auc.out[,t] = unlist(output)
  }
  #res$used.ranks=ranking
  res$auc.mat = auc.out
  return(res)
}
