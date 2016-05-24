featureselection.meta <-
function (gnExpMat, survivaltime, censor){
        zstat = NULL

        for (i in 1:ncol(gnExpMat)){
                 cox.t = coxph(Surv (survivaltime, censor)~.,data = as.data.frame (gnExpMat[,i]))
                zstat = c(zstat, summary(cox.t)[[6]][4])
        }
        return (zstat)
}

