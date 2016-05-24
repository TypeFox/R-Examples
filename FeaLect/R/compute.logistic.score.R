compute.logistic.score <- function(F_,L_,considered.features,training.samples,validating.samples,linear.scores,report.fitting.failure=TRUE) {	
    ## This function, computes the total logistic score.logistic scorelogistic scorelogistic scorelogistic
    ## OUTPUT = 	fits a logistic regression model and computes: 1/(1+exp(aX+b))
    ## INPUT: 	feature.matrix, each row is a feature
    ## 			linear.cofs, a vector containing linear for each sample
    ##			logistic.cof, a vector containg 

    ##plot(linear.scores);
    feature.matrix <- F_[training.samples, considered.features]
    ## rms library is needed for logistic regression to derive probablities from the scores computed by linear models. Used to be package Design.
    lrm.result <- try(lrm(L_ ~ linear.scores),silent=!report.fitting.failure)
    if(inherits(lrm.result, "try-error"))
        stop("lrm() failed to fit a logistic regession model.")
    logistic.cofs <- coef(lrm.result)
    ## logistic model is fitted using the validating samples.
    logistic.scores <- 1/(1+exp(-(logistic.cofs[1]+ linear.scores * logistic.cofs[2])))
    ## 1/(1+exp(-(a+bX)))
    ## The output of lrm() is aimed at getting 0 and 1 instead of -1 and 1.	
    ##message("lrm()-end"); plot(linear.scores,ylim=c(-2,2),col="green"); points(L_); points(logistic.scores,col="red"); a();	
    ## CHECK POINT for testing						
    return(list(logistic.scores=logistic.scores, logistic.cofs=logistic.cofs))
}##End compute.logistic.score <- function.	

