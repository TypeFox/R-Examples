`predictProb2.CoxBoost` <-
function(object, newdata, times, complexity, ...){
    predict(object, type="risk", newdata=as.matrix(newdata[,object$xnames,drop=FALSE]),
       times=times, at.step=complexity)
}

