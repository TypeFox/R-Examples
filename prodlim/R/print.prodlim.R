#' Print objects in the prodlim library
#' 
#' Pretty printing of objects created with the functionality of the `prodlim'
#' library.
#' 
#' 
#' @aliases print.prodlim print.neighborhood print.Hist
#' @param x Object of class \code{prodlim}, \code{Hist} and
#' \code{neighborhood}.
#' @param \dots Not used.
#' @author Thomas Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{summary.prodlim}}, \code{\link{predict.prodlim}}
#' @keywords survival
#' @export 
"print.prodlim" <- function(x,...) {
    cat("\n")
    cat("Call: ")
    print(x$call)
    cat("\n")
    model <- x$model
    ##   message("Estimation method:")
    if (!(model %in% c("survival","competing.risks"))) stop("Under construction")
    if (model=="survival")
        if (x$cens.type=="intervalCensored"){
            message(switch(x$covariate.type,"NPMLE",
                           "Stratified NPMLE estimator",
                           "Stratified NPMLE estimator",
                           "Stratified NPMLE estimator")," for the",ifelse(x$covariate.type==1," "," conditional "),ifelse(x$reverse==FALSE,"event time ","censoring time "),"survival function")
            message(paste("\nIteration steps:",x$n.iter,"\n"))
            ##   summary(x)
            cat("\n")
        }
        else{
            message(switch(x$covariate.type,"Kaplan-Meier estimator",
                           "Stratified Kaplan-Meier estimator",
                           "Stone-Beran estimator",
                           "Stratified Stone-Beran estimator")," for the",ifelse(x$covariate.type==1," "," conditional "),ifelse(x$reverse==FALSE,"event time ","censoring time "),"survival function")
        }
    cat("\n")
    ##   discrete.predictors <- extract.name.from.special(grep("strata.",names(x$X),value=TRUE),pattern="strata\\.")
    ##   continuous.predictors <- extract.name.from.special(grep("NN.",names(x$X),value=TRUE),pattern="NN\\.")
    discrete.predictors <- x$discrete.predictors
    continuous.predictors <- x$continuous.predictors
    if (!is.null(x$cluster))
        message("\nCluster-correlated data:\n\n cluster variable: ",x$cluster,"\n")
    format.disc <- function(name){
        paste(name," (",
              paste(x$xlevels[[name]],collapse=", ",sep=""),")",
              collapse=", ",sep="")
    }
    message(#"Predictor space:\n\n",
            switch(x$covariate.type,
                   "No covariates",{
                       if (length(discrete.predictors)==1){
                           c("Discrete predictor variable: ", format.disc(discrete.predictors))
                       }else{
                           c("Discrete predictor variables:\n", sapply(discrete.predictors,function(x)paste("\n - ",format.disc(x))))
                       }},
                   c("Continuous predictors: ",continuous.predictors),
                   c("  Discrete predictor variables: ",
                     paste(discrete.predictors,collapse=", "),
                     "\nContinuous predictor variables: ",
                     continuous.predictors)))
    summary(x$model.response,verbose=TRUE)
    if (!is.null(x$na.action)){
        cat("\n",
            length(x$na.action),
            ifelse(length(x$na.action)==1,
                   " observation",
                   " observations")," deleted due to missing values.\n",sep="")
    }
}
