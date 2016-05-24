#' Predict individual survival probabilities
#' 
#' Function to extract the predicted probabilities at the individual event
#' times that have been used for fitting a prodlim object.
#' 
#' 
#' @param object A fitted object of class "prodlim".
#' @param lag Integer. `0' means predictions at the individual times, 1 means
#' just before the individual times, etc.
#' @return A vector of survival probabilities.
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{predict.prodlim}},\code{\link{predictSurv}},
#' @keywords survival
#' @examples
#' 
#'  SurvFrame <- data.frame(time=1:10,status=rbinom(10,1,.5))
#'  x <- prodlim(formula=Hist(time=time,status!=0)~1,data=SurvFrame)
#'  predictSurvIndividual(x,lag=1)
#' 
#' @export
predictSurvIndividual <- function(object,
                                  lag=1){
    obs.times <- as.numeric(object$model.response[,1])
    if (object$covariate.type==1){
        locOBS <- match(obs.times,object$time,nomatch=FALSE)
        if (any(locOBS==FALSE)) stop("Can't locate all individual observation times" )
        psurv <- c(rep(1,lag),object$surv)[locOBS]}
    else{
        N <- length(obs.times)
        if (is.null(object$model.matrix)) stop("Cannot find model.matrix, need to set x=TRUE in call of prodlim.")
        findex <- row.match(object$model.matrix,object$X)
        ## if (any(is.na(findex)))
        ## stop("Cannot identify all rows in object$model.matrix in ")
        psurv <- .C("predict_individual_survival",pred=double(N),as.double(object$surv),as.double(object$time),as.double(obs.times),as.integer(object$first.strata[findex]),as.integer(object$size.strata[findex]),as.integer(N),as.integer(lag),NAOK=FALSE,PACKAGE="prodlim")$pred}
    psurv
}
