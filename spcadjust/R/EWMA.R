########################
## Basic EWMA Charts ##
########################

#' @include main.R model.R EWMAlib.R
NULL

#' EWMA Charts
#'
#' Class extending SPCChart with a basic  EWMA charts implementation.
#'
#' Let \eqn{Y_t, t=1,2,\dots} be the updates from the data model. Then
#' the EWMA chart is given by \eqn{Q_0=0} and
#' \deqn{Q_t=lambda Y_t+(1-lambda) Q_{t-1}}
#'
#' @examples
#' X <-  rnorm(1000)
#' chart <- new("SPCEWMA",model=SPCModelNormal(Delta=0),lambda=0.8)
#' \dontrun{
#' SPCproperty(data=X,nrep=10,chart=chart,
#'             property="calARL",params=list(target=100))
#' SPCproperty(data=X,nrep=10,chart=chart,
#'             property="calhitprob",params=list(target=0.05,nsteps=1e3))
#' }
#' SPCproperty(data=X,nrep=10,chart=chart,
#'             property="ARL",params=list(threshold=3))
#' SPCproperty(data=X,nrep=10,chart=chart,
#'             property="hitprob",params=list(threshold=3,nsteps=1e3))
#' #increase the number of repetitions nrep for real applications.
#'
#' @slot model  The data model. The  data  model should center the in-control
#' updates such that they have mean 0.
#' @slot lambda The smoothing constant, \eqn{0<lambda<=1}.
#' @export
setClass("SPCEWMA",contains="SPCchart",slots=list(lambda="numeric"))

#' @describeIn runchart Generic function for running EWMA
#' charts. Relies on \code{\link{updates}} being implemented for the
#' chart.
#' @export
setMethod("runchart", signature="SPCEWMA", function(chart,newdata,xi){
 updates <- chart@model$updates(xi=xi, data=newdata)
 Zi <- 0
 vapply(updates, function(y) Zi <<-  chart@lambda*y+(1-chart@lambda)*Zi, 0) 
})


#' @describeIn getq  Implements the properties \code{ARL},
#' \code{calARL}, \code{hitprob} and \code{calhitprob}.
#'
#' @export
setMethod("getq", signature="SPCEWMA",function(chart,property,params){
    if (is.null(params[["gridpoints"]])) params$gridpoints=100;
    switch(property,
           "calARL"=
               list(q= function(P,xi)
                   log(calibrateARL_Markovapprox(f=ARL_EWMA_Markovapprox,
                                                 pobs=getcdfupdates(chart,xi=xi, P=P),
                                                 ARL=params$target,
                                                 gridpoints=params$gridpoints,
                                                 lambda=chart@lambda)),
                    trafo=function(x) exp(x),
                    lowerconf=TRUE,
                    format=function(res)
                        paste("A threshold of +/-  ", format(res,digits=4),
                              " gives an in-control ARL of at least ",
                              params$target, ".", sep="",collapse="")
                    ),
           "ARL"=
               list(q= function(P,xi){
                   as.double(log(ARL_EWMA_Markovapprox(c=params$threshold,pobs=getcdfupdates(chart,xi=xi, P=P),gridpoints=params$gridpoints,lambda=chart@lambda)))
               },
                    trafo=function(x) exp(x),
                    lowerconf=FALSE,
                    format=function(res)
                        paste("A threshold of +/-  ", params$threshold,
                              " gives an in-control ARL of at least ",
                              format(res,digits=4), ".", sep="",collapse="")
                    ),
           "hitprob"=
               list(q=function(P,xi){
                            res <- hitprob_EWMA_Markovapprox(c=params$threshold,pobs=getcdfupdates(chart,xi=xi, P=P),n=params$nsteps,gridpoints=params$gridpoints,lambda=chart@lambda);
                            as.double(log(res/(1-res)))
                        },
                    trafo=function(x) exp(x)/(1+exp(x)),
                    lowerconf=TRUE,
                    format=function(res) paste("A threshold of +/-  ", params$threshold,
                        " gives an in-control false alarm probability of at most ",
                        format(res,digits=4),
                        " within ",params$nsteps," steps.",
                        sep="",collapse="")
                    ),
           "calhitprob"=
               list(q=function(P,xi)
                   log(calibratehitprob_Markovapprox(f=hitprob_EWMA_Markovapprox,
                                                     pobs=getcdfupdates(chart,xi=xi, P=P),
                                                     hprob=params$target,
                                                     n=params$nsteps,
                                                     gridpoints=params$gridpoints,
                                                     lambda=chart@lambda)),
                    trafo=function(x) exp(x),
                    lowerconf=TRUE,
                    format=function(res) paste("A threshold of +/-  ",
                        format(res,digits=4),
                        " gives an in-control false alarm probability of at most ",
                        params$target, " within ",params$nsteps, " steps.",
                        sep="",collapse="")
                    ),
           stop(paste("Property",property,"not implemented."))
           )
})

