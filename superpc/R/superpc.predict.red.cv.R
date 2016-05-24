superpc.predict.red.cv <- function(fitred, fitcv, data, threshold,  sign.wt="both"){

 # try reduced predictor on cv folds, via full cross-validation

                           
  this.call=match.call()

  type=fitred$type

  n.components=fitred$n.components


  n.fold<-length(fitcv$folds)

  shrinkages<- fitred$shrinkages
  n.shrinkages<- length(shrinkages)
  cur.vall<- array(NA,c(n.shrinkages,ncol(data$x),n.components))

import.cv=array(NA,c(nrow(data$x), n.fold,n.components))

 lrtest.reduced<-array(NA,c(n.fold,n.shrinkages, n.components))

  for(j in 1:n.fold){
    cat(j,fill=TRUE)
    fit.temp<-list(feature.scores=fitcv$featurescores.fold[,j], type=type)
    ii<-fitcv$folds[[j]]
    
    data1<-list(x=data$x[,-ii],y=data$y[-ii],censoring.status=data$censoring.status[-ii])
    data2<-list(x=data$x[,ii],y=data$y[ii],censoring.status=data$censoring.status[ii])
    junk<- superpc.predict.red(fit.temp, data1,data2, threshold, shrinkages=shrinkages, n.components=n.components,  compute.lrtest=TRUE, sign.wt=sign.wt)
 lrtest.reduced[j,,]=junk$lrtest.reduced
import.cv[,j,]=junk$import
  }

 mean.na <- function(x) {
            mean(x[!is.na(x)])
        }
        se.na <- function(x) {
            val = NA
            if (sum(!is.na(x)) > 0) {
                val = sqrt(var(x[!is.na(x)])/sum(!is.na(x)))
            }
            return(val)
        }

   llr= apply(log(lrtest.reduced), c(2,3), mean.na)
        se.llr = apply(log(lrtest.reduced), c(2,3), se.na)
        lrtest.reduced.lower = exp(llr - se.llr)
        lrtest.reduced.upper = exp(llr + se.llr)
        lrtest.reduced <- exp(llr)


return(list(shrinkages=shrinkages, lrtest.reduced=lrtest.reduced,  lrtest.reduced.lower= lrtest.reduced.lower, lrtest.reduced.upper=lrtest.reduced.upper,  n.components=n.components, num.features=fitred$num.features,  sign.wt=sign.wt,import.cv=import.cv, type=type,call=this.call))
}
