selectCox <- function(formula,data,rule="aic"){
    ## require(rms)
    ## require(prodlim)
    fit <- rms::cph(formula, data, surv=TRUE)
    bwfit <- rms::fastbw(fit,rule=rule)
    if (length(bwfit$names.kept)==0){
        newform <- update(formula,".~1")
        newfit <- prodlim(newform,data=data)
    }
    else{
        newform <- update(formula,paste(".~",paste(bwfit$names.kept,collapse="+")))
        ## newform <- reformulate(bwfit$names.kept, formula[[2]])
        newfit <- cph(newform,data, surv=TRUE)
    }
    out <- list(fit=newfit,In=bwfit$names.kept)
    out$call <- match.call()
    class(out) <- "selectCox"
    out
}

##' @export
predictSurvProb.selectCox <- function(object,newdata,times,...){
   predictSurvProb(object[[1]],newdata=newdata,times=times,...)
 }
