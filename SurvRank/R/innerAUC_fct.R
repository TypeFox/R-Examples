#' Calculates inner and outer survival AUCs of training and testset
#'
#' work-horse function for all ranking methods with inner and outer CV loops
#' @param f ranking function
#' @param data data input
#' @param t current t.times argument
#' @param cv.out number of outer CV
#' @param cv.in number of inner CV
#' @param i current counter index
#' @param fold current fold
#' @param data.out training data of outer CV
#' @param nr.var maximum number of variables into model
#' @param out.s predefined output matrix for inner tAUC evaluations
#' @param sd1 factor to which sparser solutions should be chosen. Not maximum Survival AUC in inner loop is used in stepwise selection, instead \code{max(survAUC)*sd1} leading to sparser solutions
#' @param c.time as defined in package \code{survAUC} time; a positive number restricting the upper limit of the time range under consideration.
#' @param ranking predefined ranking list
#' @param used.rank predefined list
#' @param used.rank1se predefined list
#' @param pred.in predefined list for inner predictions
#' @param pred.out predefined list for outer predictions
#' @param pred.out1 predefined list for outer predictions with \code{sd1} factor
#' @param auc.out vector of survival AUCs
#' @param auc.out1 vector of survival AUCs with \code{sd1} factor
#' @keywords internal
#' @export

innerAUC_fct = function(f,data,t,cv.out,cv.in,i,fold,data.out,nr.var,out.s,sd1,c.time,ranking,used.rank,used.rank1se,pred.in,pred.out,pred.out1,auc.out,auc.out1){
  l=(t-1)*cv.out+i                      # counter
  data.out$x = data$x[fold!=i,]     #create the training data set for this loop
  data.out$y = data$y[fold!=i,]
  fold.in = crossvalFolds(data.out$y[,2],k=cv.in)  #create stratified cross validation folds for the inner cv-loop
  for(k in 1:cv.in){                     #start the inner cross validation loop
    data.in = list()
    data.in$x = data.out$x[fold.in!=k,]       #create the inner training data set
    data.in$y = data.out$y[fold.in!=k,]
    ranking.i = tryCatch(f(x=data.in$x,y=data.in$y), error = function(e) NA)   #do the ranking of the variables using the predefined function
    if(is.na(ranking.i)){
      ranking.i = list()
      ranking.i$rank = colnames(data.in$x)  #if the ranking function produces an error, just take the variables as they are in the data set
      ranking.i$nrcoef = ncol(data.in$x)
    }
    for (j in 1:min(nr.var,ranking.i$nrcoef,na.rm=T)){    #start a new loop going from 1 to the number of ranked variables
      glm.in = tryCatch(survival::coxph(y~.,data = data.frame(y=data.in$y,subset(data.in$x,select=ranking.i$rank[1:j]))), error = function(e) NA) #fit a cox model for the inner training set using the ranked variables from 1 to j
      pred.in = tryCatch(predict(glm.in,newdata = data.frame(data.out$x[fold.in==k,]),type="lp",reference="sample"),error = function(e) NA) #predict the inner test set
      out.s[k,j] = tryCatch(survAUC::UnoC(Surv.rsp = data.in$y,Surv.rsp.new=data.out$y[fold.in==k,],lpnew = pred.in,time = c.time),error = function(e) NA) #calculate the estimated C-Statistic according to Uno for the given data and model
    }
  }
  out.score=apply(out.s,2,mean,na.rm=T)  #take the mean of the C-statistics for the different number of parameters
  nr.par.in = which.max(out.score)      # number of parameters where the C-statistic takes the largest value
  if(length(nr.par.in)==0){           # if the maximum is at zero, set the number of parameters to the number of variables
    nr.par.in=nr.var
  }
  if(nr.par.in==1){nr.par.1se=1}
  else{nr.par.1se=tryCatch(ceiling(approx(x=c(0.5,out.score[1:nr.par.in]),y=0:nr.par.in,xout=max(out.score,na.rm=T)*sd1)$y),error=function(e) NA)}
  if(is.nan(nr.par.1se) | is.na(nr.par.1se)){nr.par.1se=nr.var}
  ranking[[l]] = f(x=data.out$x,y=data.out$y)       #apply the ranking using the whole training set
  used.rank[[l]] = ranking[[l]]$rank[1:nr.par.in]   # number of parameters that are used
  used.rank1se[[l]]=ranking[[l]]$rank[1:nr.par.1se]
  glm.out = tryCatch(survival::coxph(y~.,data = data.frame(y=data.out$y,subset(data.out$x,select=ranking[[l]]$rank[1:nr.par.in]))),error = function(e) NA) #fit a cox model for the training set using the ranked variables from 1 to the previously determined number of parameters
  glm.out1 = tryCatch(survival::coxph(y~.,data = data.frame(y=data.out$y,subset(data.out$x,select=ranking[[l]]$rank[1:nr.par.1se]))),error = function(e) NA) #fit a cox model for the training set using the ranked variables from 1 to the more sparsely chosen number of parameters
  pred.out[[l]] = tryCatch(predict(glm.out,data.frame(data$x[fold==i,]),type="lp",reference="sample"),error = function(e) NA) #predict the test set
  pred.out1[[l]] = tryCatch(predict(glm.out1,data.frame(data$x[fold==i,]),type="lp",reference="sample"),error = function(e) NA) #predict the test set
  auc.out[l] = tryCatch(survAUC::UnoC(Surv.rsp=data.out$y,Surv.rsp.new=data$y[fold==i,],lpnew=pred.out[[l]],time=c.time),error = function(e) NA) #calculate the AUC according to Uno (C-Statistic)
  auc.out1[l] = tryCatch(survAUC::UnoC(Surv.rsp=data.out$y,Surv.rsp.new=data$y[fold==i,],lpnew=pred.out1[[l]],time=c.time),error = function(e) NA) #calculate the AUC according to Uno (C-Statistic)
  return(list(ranking = ranking[[l]], pred.in = out.score, pred.out=pred.out[[l]],used.rank = used.rank[[l]], auc.out = auc.out[l],
              used.rank1se = used.rank1se[[l]],auc.out1 = auc.out1[l])) #return a list with the ranking, the AUC for the different number of parameters, the final prediction, the AUC of the final model, the used rank, as well as the used rank and the AUC for the sparser models

}
