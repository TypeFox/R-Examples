#####################
## Shewhart Charts ##
#####################

#' @include main.R model.R
NULL

#' Shewhart charts.
#'
#' @examples
#' X<-rnorm(100);
#' ##calibrate to ARL 100
#' chartShew <- new("SPCShew",model=SPCModelNormal(),twosided=TRUE)
#' \dontrun{
#' SPCproperty(data=X,nrep=500,
#'             property="calARL",chart=chartShew,params=list(target=100),
#'             covprob=c(0.7,0.9))
#'
#' chartShewOneSided <- new("SPCShew",model=SPCModelNormal(),twosided=FALSE)
#' SPCproperty(data=X,nrep=500,
#'             property="calARL",chart=chartShewOneSided,
#'             params=list(target=100),covprob=c(0.7,0.9))
#'
#' ##calibrate to a hitting probability of 0.01 in 100 steps
#' SPCproperty(data=X,nrep=500,
#'             property="calhitprob",
#'             chart=chartShew,params=list(target=0.01,nsteps=100))
#' SPCproperty(data=X,nrep=500,
#'             property="calhitprob",chart=chartShewOneSided,params=list(target=0.01,nsteps=100))
#'
#' ## work out  for ARL for a fixed threshold of 4
#' SPCproperty(data=X,nrep=500,
#'             property="ARL",chart=chartShew,params=list(threshold=4))
#' SPCproperty(data=X,nrep=500,
#'             property="ARL",chart=chartShewOneSided,
#'             params=list(threshold=4))
#'
#' SPCproperty(data=X,nrep=500,
#'             property="hitprob",chart=chartShew,params=list(nsteps=100,threshold=4))
#'
#' SPCproperty(data=X,nrep=500,
#'             property="hitprob",chart=chartShewOneSided,params=list(nsteps=100,threshold=4))
#' }
#' 
#' X<-rnorm(100)
#' chartShew <- new("SPCShew",model=SPCModelNormal())
#' \dontrun{
#' SPCproperty(data=X,nrep=500,
#'             property="calARL", chart=chartShew,
#'             params=list(target=1000))
#' SPCproperty(data=X,nrep=500,
#'             property="calhitprob",chart=chartShew,
#'             params=list(target=0.01,nsteps=100))
#' SPCproperty(data=X,nrep=10,chart=chartShew,
#'             property="ARL",params=list(threshold=3))
#' SPCproperty(data=X,nrep=500,
#'             property="hitprob",
#'             chart=chartShew,params=list(nsteps=100,threshold=4))
#' }
#' @slot twosided  TRUE if a two-sided chart should be used. Default FALSE.
#' @export
setClass("SPCShew",contains=c("SPCchart"),
         slots=list(twosided="logical"),
         prototype=list(twosided=FALSE))
#' @describeIn runchart Simply computes the updates appropriate for
#' the Shewhart chart and returns them.
#' @export
setMethod("runchart", signature="SPCShew", function(chart,newdata,xi){
    updates(chart,xi=xi, data=newdata)
})
#' @describeIn updates Computes the updates taking into account if the
#' chart is one-sided or two-sided.
#' @export
setMethod("updates", signature="SPCShew", function(chart,xi,data){
    if (!chart@twosided)
        chart@model$updates(xi,data)
    else
        abs(chart@model$updates(xi,data))
})

#' @describeIn getcdfupdates Computes the CDF of the updates taking
#' into account if the chart is one-sided or two-sided.
#' @export
setMethod("getcdfupdates", signature="SPCShew", function(chart,P, xi){
    if (!chart@twosided)
        chart@model$getcdfupdates(P=P,xi=xi)
    else{
        f <- chart@model$getcdfupdates(P=P,xi=xi)
        fcaglad <- chart@model$getcdfupdates(P=P,xi=xi,cadlag=FALSE)
        function(x) ifelse(x>=0,f(abs(x))-fcaglad(-abs(x)),0)
}})



#' @describeIn getq Implements the properties \code{ARL},
#' \code{calARL}, \code{hitprob} and \code{calhitprob}.
#' @export
setMethod("getq", signature="SPCShew",function(chart,property,params){
    switch(property,
           ARL=list(
               q=function(P,xi)
                   -log(1-getcdfupdates(chart, xi=xi, P=P)(params$threshold)),
               trafo=function(x) exp(x),
               lowerconf=FALSE,
               format=function(res)
                   paste("A threshold of  ", params$threshold,
                         " gives an in-control ARL of at least ",
                         format(res,digits=4), ".", sep="",collapse="")
           ),
           calARL=list(
               q=function(P,xi){
                   pobs <- getcdfupdates(chart, xi=xi, P=P)
                   uniroot(function(x) params$target-(1/(1-pobs(x))),
                           lower=-1,upper=1,extendInt="downX")$root
               },
               trafo=function(x) x,
               lowerconf=TRUE,
               format=function(res)
                   paste("A threshold of ", format(res,digits=4),
                         " gives an in-control ARL of at least ",
                         params$target, ".", sep="",collapse="")
           ),
           hitprob=list(
               q=function(P,xi){
                   survprob1step <- getcdfupdates(chart, xi=xi, P=P)(params$threshold)
                   log(1-survprob1step^params$nsteps)-params$nsteps*log(survprob1step)
               },
               trafo=function(x) exp(x)/(1+exp(x)),
               lowerconf=TRUE,
               format=function(res) paste("A threshold of  ", params$threshold,
                   " gives an in-control false alarm probability of at most ",
                   format(res,digits=4),
                   " within ",params$nsteps," steps.",
                   sep="",collapse="")
           ),
           calhitprob=list(
               q= function(P,xi){
                   pobs <- getcdfupdates(chart, xi=xi, P=P)
                   uniroot(function(x) params$target-(1-pobs(x)^params$nsteps),
                           lower=-1,upper=1,extendInt="upX")$root
               },
               trafo=function(x) x,
               lowerconf=TRUE,
               format=function(res) paste("A threshold of ",
                   format(res,digits=4),
                   " gives an in-control false alarm probability of at most ",
                   params$target, " within ",params$nsteps, " steps.",
                   sep="",collapse="")
           ),
           stop(paste("Property",property,"not implemented."))
           )
})



