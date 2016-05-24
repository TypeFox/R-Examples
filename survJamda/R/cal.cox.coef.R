cal.cox.coef <-
function (gnExpMat, survivaltime, censor){
	cox.coef = NULL

        max.col = ifelse (is.matrix(gnExpMat), ncol(gnExpMat), 1)

        for (i in 1:max.col){
                if(is.matrix(gnExpMat))
                        var = gnExpMat[,i]
                else
                        var = gnExpMat
                cox.t = coxph(Surv (survivaltime, censor)~var)          
                cox.coef = c(cox.coef, cox.t$coef)
        }
  return (cox.coef)
}

