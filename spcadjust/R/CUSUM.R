########################
## Basic CUSUM Charts ##
########################

#' @include main.R model.R CUSUMlib.R
NULL

#' CUSUM Charts
#'
#' Class extending SPCChart with a basic  CUSUM charts implementation.
#'
#' The only slot this class contains is the data model. This data
#' model should already incorporate the negative mean for in-control
#' updates that is typical for CUSUM charts.
#'
#' Let \eqn{U_t, t=1,2,\dots} be the updates from the data model. Then
#' the CUSUM chart is given by \eqn{S_0=0} and
#' \deqn{S_t=max(S_{t-1}+U_t,0)}
#' 
#' @examples
#' X <-  rnorm(1000)
#' chart <- new("SPCCUSUM",model=SPCModelNormal(Delta=1))
#' \dontrun{
#' SPCproperty(data=X,nrep=10,chart=chart,
#'             property="calARL",params=list(target=100))
#' SPCproperty(data=X,nrep=10,chart=chart,
#'             property="calhitprob",params=list(target=0.05,nsteps=1e3))
#' SPCproperty(data=X,nrep=10,chart=chart,
#'             property="ARL",params=list(threshold=3))
#' }
#' SPCproperty(data=X,nrep=10,chart=chart,
#'             property="hitprob",params=list(threshold=3,nsteps=1e3))
#' #increase the number of repetitions nrep for real applications.
#'
#' @export
setClass("SPCCUSUM",contains="SPCchart")

#' @describeIn runchart Generic function for running CUSUM
#' charts. Relies on \code{\link{updates}} being implemented for the
#' chart.
#' @export
setMethod("runchart", signature="SPCCUSUM", function(chart,newdata,xi){
    R <- cumsum(chart@model$updates(xi=xi, data=newdata))
    R - cummin(R)
})

#' @describeIn getq Implements the properties \code{ARL},
#' \code{calARL}, \code{hitprob} and \code{calhitprob}.
#'
#' @export
setMethod("getq", signature="SPCCUSUM",function(chart,property,params){
    if (is.null(params[["gridpoints"]])) params$gridpoints=75;
    if (grepl("cal",property)&&is.null(params$target))
        stop("Argument params contains no element target (needed for calibration).")
    if (grepl("hitprob",property)){
        if (is.null(params$nsteps))
            stop("Argument params does not contain an element nsteps (needed for hitting probabilities).")
        else if (params$nsteps<1|round(params$nsteps)!=params$nsteps)
            stop("nsteps has to be a positive integer.")
        
    }
    if (is.element(property,c("ARL","hitprob"))){
        if (is.null(params$threshold))
            stop("Argument params does not contain an element threshold.")
        else
            if (params$threshold<0) stop("Negative threshold.")
    }
    switch(property,
           "calARL"= 
               list(q= function(P,xi)
                   log(calibrateARL_Markovapprox(pobs=getcdfupdates(chart,xi=xi, P=P),
                                                 ARL=params$target,
                                                 gridpoints=params$gridpoints)),
                    trafo=function(x) exp(x),
                    lowerconf=TRUE,
                    format=function(res)
                        paste("A threshold of ", format(res,digits=4),
                              " gives an in-control ARL of at least ",
                              params$target, ".", sep="",collapse="")
                    ),
           "ARL"=
           list(q= function(P,xi){
               as.double(log(ARL_CUSUM_Markovapprox(c=params$threshold,pobs=getcdfupdates(chart,xi=xi, P=P),gridpoints=params$gridpoints)))
           },
                trafo=function(x) exp(x),
                lowerconf=FALSE,
                format=function(res)
                paste("A threshold of  ", params$threshold,
                      " gives an in-control ARL of at least ",
                      format(res,digits=4), ".", sep="",collapse="")
                ),
           "hitprob"=
               list(q=function(P,xi){
                            res <- hitprob_CUSUM_Markovapprox(c=params$threshold,pobs=getcdfupdates(chart,xi=xi, P=P),n=params$nsteps,gridpoints=params$gridpoints);
                            as.double(log(res/(1-res)))
                        },
                    trafo=function(x) exp(x)/(1+exp(x)),
                    lowerconf=TRUE,
                    format=function(res) paste("A threshold of  ", params$threshold,
                        " gives an in-control false alarm probability of at most ",
                        format(res,digits=4),
                        " within ",params$nsteps," steps.",
                        sep="",collapse="")
                    ),
           "calhitprob"=
               list(q=function(P,xi)
                   log(calibratehitprob_Markovapprox(pobs=getcdfupdates(chart,xi=xi, P=P),
                                                     hprob=params$target,
                                                     n=params$nsteps,
                                                     gridpoints=params$gridpoints)),
                    trafo=function(x) exp(x),
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

