
superpc.fit.to.outcome<- function(fit, data.test,score, competing.predictors=NULL,  print=TRUE, iter.max=5){


type=fit$type

if(type=="survival"){temp.list=makelist(data.test$y, data.test$censoring.status, score)}
if(type=="regression"){temp.list=makelist(data.test$y,NULL,  score)}

if(!is.null(competing.predictors)){
 temp.list=c(temp.list,competing.predictors)
}


 if(type=="survival"){
   require(survival)
   results<-coxph(Surv(y, censoring.status)~., data=temp.list, control=coxph.control(iter.max=iter.max))
}

 else{
   results<-lm(data.test$y~.,  data=temp.list)
}


if(print){print(summary(results))}


ss=summary(results)
if(type=="survival"){ test.stat=ss$logtest[1]
                      df=ss$logtest[2]
                      pvalue=ss$logtest[3]
                     }
if(type=="regression"){ test.stat=ss$fstat[1]
                      df=ss$fstat[2:3]
                      pvalue=1-pf(test.stat,df[1],df[2])
                     }

teststat.table=matrix(c(test.stat, df, pvalue), nrow=1)
if(length(df)==1){dflabel="df"}
if(length(df)==2){dflabel=c("df1", "df2")}

dimnames(teststat.table)=list(NULL,c("test statistic",dflabel,"p-value"))


return(list(results=results, teststat.table=teststat.table,  coeftable=ss$coef))
}

makelist=function (y, censoring.status, predictors)
{
    val = list(y = y)
    if (!is.null(censoring.status)) {
        val$censoring.status = censoring.status
    }
    if (!is.matrix(predictors)) {
        val$score.1 = predictors
    }

    if (is.matrix(predictors)) {
        if (ncol(predictors) > 3) {
            stop("Can't have > 3 principal components")
        }
predictor.type=dimnames(predictors)[[2]]

if(is.null(dimnames(predictors)[[2]])){
  predictor.type=rep("continuous",ncol(predictors))
 }
        score1 = predictors[, 1]
        if(predictor.type[1]=="factor") {
            score1 = as.factor(score1)
        }
        val$score.1 = score1
        if (ncol(predictors) > 1) {
            score2 = predictors[, 2]
 if(predictor.type[2]=="factor") {
                score2 = as.factor(score2)
            }
            val$score.2 = score2
        }
        if (ncol(predictors) > 2) {
            score3 = predictors[, 3]
 if(predictor.type[3]=="factor") {
                score3 = as.factor(score3)
            }
            val$score.3 = score3
        }
    }
    return(val)
}
