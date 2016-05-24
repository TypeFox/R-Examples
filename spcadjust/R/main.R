############################################
## Base Classes and Computational Methods ##
############################################

## for NAMESPACE file
#' @import methods stats graphics utils
#' @include model.R
NULL

#' Virtual Base Class for Control Charts
#'
#' Virtual S4 base class for Control Charts.
#'
#' @slot model The data model to be used in the chart. Must be of type
#' \code{\link{SPCDataModel}}.
#'
#' @export
setClass("SPCchart",representation=list(model="SPCDataModel"))

#' Estimate Chart Parameters
#'
#' Estimates the parameters used to run a control chart from given
#' data.
#'
#' @param chart the control chart to be used.
#' @param data the data to be used.
#'
#' @return The parameter values for running the control chart \code{chart}.
#' @export
setGeneric("xiofdata",def=function(chart,data){standardGeneric("xiofdata")})

#' @describeIn xiofdata Standard method which simply first applies
#' PofData to get a model and then xiofP to get the parameters.
#' @export
setMethod("xiofdata", signature="SPCchart", function(chart,data){
    chart@model$xiofP(chart@model$Pofdata(data))
})



#' Updates of a Control Chart
#'
#' Computes updates of a control chart using the given parameters and the given data.
#'
#' @param chart the control chart.
#' @param xi the parameters used for running the chart.
#' @param data the observed data.
#'
#' @return A vector of the same length as data.
#' @export
setGeneric("updates", def= function(chart,xi,data)  standardGeneric("updates"))

#' @describeIn updates Standard method which simply first applies
#' getresiduals from the data model.
#' @export
setMethod("updates", signature="SPCchart", function(chart,xi,data){
    chart@model$updates(xi,data)
})


#' CDF of Updates of a Control Chart
#'
#' Consider running a control chart with given parameters with data
#' coming from a given probability model. This function computes the
#' cumulative distribution function (CDF) of the updates of the
#' control charts as they would be computed by the method \code{\link{updates}}.
#'
#' @param chart the chart to be used.
#' @param P the probability model from which data is generated.
#' @param xi the parameters of the control chart.
#'
#' @return A function mapping one-dimensional numerical values into
#' the interval [0,1], having all properties of a cumulative
#' distribution function.
#'
#' @export
setGeneric("getcdfupdates", def= function(chart, P, xi) standardGeneric("getcdfupdates"))


#' @describeIn getcdfupdates Standard method which simply first applies
#' getcdfresiduals from the data model.
#' @export
setMethod("getcdfupdates", signature="SPCchart", function(chart,P, xi){
    chart@model$getcdfupdates(P=P,xi=xi)
})


#' Runs a chart
#'
#' Generic method for running a chart on new data using given
#' parameters \code{xi}.
#' 
#' @param chart the chart to be used.
#' @param newdata the new observed data.
#' @param xi the parameters to be used in running the chart.
#'
#' @return The path of the chart over time.
#'
#' @export
setGeneric("runchart",def=function(chart,newdata,xi){standardGeneric("runchart")})

#' Returns a List to Compute Properties of a chart
#'
#' Returns functions to compute desired properties of a given
#' control chart.
#' 
#' @param chart the chart to be used.
#' @param property the name of the property.
#' @param params additional parameters needed for the computations.
#'
#' @return A list with the elements \code{q}, \code{trafo}, \code{lowerconf}, \code{format}.
#' \itemize{
#' \item{\code{q(P,xi)}: The transformed property of interest. To improve the bootstrap a log transform is used for \code{calARL},\code{calhitprob} and \code{ARL}, and a logit transform for \code{hitprob}. This function depends on the distribution of updates \code{P} and the chart parameters \code{xi}.} 
#' \item{\code{trafo(x)}: The inverse of the transformation of the property used in the bootstrap. Needed to back-transform the result to the correct scale. }
#' \item{\code{lowerconf}: Logical value. TRUE if a lower confidence interval should be reported, FALSE otherwise. Default is TRUE for properties \code{calARL},  \code{calhitprob} and \code{hitprob} and FALSE for \code{ARL}. }
#' \item{\code{format(res)}: Output summary given as a text string. }
#' }
#' @export
setGeneric("getq",def=function(chart,property,params){standardGeneric("getq")})




#################
#Computational routines

#' Computes bootstrap adjusted properties for control charts
#'
#' Computes bootstrap adjusted properties for control charts.
#'
#' @param data The observed data.
#' @param nrep The number of bootstrap repetitions. Default 500.
#' @param chart The chart to be used.
#' @param property The property to be computed. A string. Must be implemented by the chart.
#' @param params Additional parameters for computing the property.
#' @param covprob The coverage probability of the adjustment. Default 0.9.
#' @param quiet Logical value indicating if progress bar should be suppressed. Default FALSE.
#' @param reportdistr Logical value indicating if the ecdf of the bootstrap distribution should be plotted. Default FALSE.
#' @param parallel Number of cores to use for parallel computations
#' (using mclapply from the package parallel). Defaults to 1. If set
#' to \code{Inf} then the number of cores is automatically detected
#' and all but one are used.
#' 
#'
#' @return An object of type SPCpropertyres.
#'
#' @seealso  \code{\link{SPC2sidedconfint}}
#'
#' @examples
#'  # calibrate CUSUM chart to an in-control ARL of 100.
#'  # run with a larger number of replications in real examples!
#'  X <-  rnorm(100) #observed data
#'  chart <- new("SPCCUSUM",model=SPCModelNormal(Delta=1)) # CUSUM chart with normal observations
#' SPCproperty(data=X,nrep=15,chart=chart,property="calARL", params=list(target=100))
#'
#' @export
SPCproperty <- function(data, nrep=500,chart, property,params,covprob=0.9,quiet=FALSE,reportdistr=FALSE,parallel=1){
    model <- chart@model
    Phat <-  model$Pofdata(data)
    qfunc <- getq(chart,property,params)
    if(parallel>1){
        if (!requireNamespace("parallel", quietly = TRUE)) {
            stop("parallel needed for this function to work. Please install it.",
                 call. = FALSE)
        }
        if (parallel==Inf)
            cores <- parallel::detectCores()-1
        else
            cores <- parallel;
        qdifffunc <- function(...){
            qdiff <- replicate(ceiling(nrep/cores),{
                Phatstar<-model$Pofdata(model$resample(Phat))
                xihatstar <- model$xiofP(Phatstar)
                qfunc$q(xi=xihatstar,P=Phatstar)-qfunc$q(xi=xihatstar,P=Phat)
            })
            qdiff
        }
        qdiff <- do.call(c,parallel::mclapply(seq_len(cores),qdifffunc,mc.cores=cores))
    }else{
        if (!quiet){
            pb <- txtProgressBar(min=0,max=nrep+1,initial=0,style=3)
            pbpos <- 0;
        }
      qdiff <- replicate(nrep,{
        if (!quiet){
         setTxtProgressBar(pb, pbpos)
          pbpos <<- pbpos+1;
        }
        Phatstar<-model$Pofdata(model$resample(Phat))
        xihatstar <- model$xiofP(Phatstar)
        qfunc$q(xi=xihatstar,P=Phatstar)-qfunc$q(xi=xihatstar,P=Phat)
      })
        if (!quiet) setTxtProgressBar(pb, pbpos)
    }
    raw <- qfunc$q(xi=model$xiofP(Phat),P=Phat)
    if (parallel==1&&!quiet) setTxtProgressBar(pb, nrep+1)
    if (qfunc$lowerconf)
        res <- qfunc$trafo(raw+quantile(-qdiff,covprob))
    else
        res <- qfunc$trafo(raw-quantile(qdiff,covprob))
    if (parallel==1&&!quiet) close(pb)
    if(reportdistr)
      plot(ecdf(qdiff), main="ECDF of bootstrap distribution")
    new("SPCpropertyres", res=res, raw=qfunc$trafo(raw),
        covprob=covprob,chart=chart,property=property,
        nrep=nrep,params=params,restext=sapply(res,qfunc$format))
}

#' Results of SPCproperty.
#'
#' @slot nrep number of repetitions used in the simulation.
#' @slot chart the chart used.
#' @slot property the property of interest, \code{ARL}, \code{calARL}, \code{calhitprob} or \code{hitprob}. 
#' @slot covprob the probability of the guaranteed conditional performance.
#' @slot res the guaranteed conditional performance.
#' @slot raw the unadjusted result.
#' @slot params additional parameters used for computing this property.
#' @slot restext a readable version of the result.
#'
#' @export
setClass("SPCpropertyres", representation=list(nrep="numeric", chart="SPCchart",property="character", covprob="numeric",res="numeric",raw="numeric",params="list",restext="character"))
#' @describeIn SPCpropertyres Prints the object nicely.
#' @param object the result to be shown.
#' @export
setMethod("show", "SPCpropertyres",  function(object){
    for (i in 1:length(object@covprob)){
        cat(paste(strwrap(exdent=2,paste(object@covprob[i]*100,"% CI: ",object@restext[i])),collapse="\n"),"\n")
    }
    cat("Unadjusted result: ", format(object@raw,digits=4),"\n")
    cat("Based on ",object@nrep, "bootstrap repetitions.\n")
})


#' Computes a two-sided confidence interval for properties of a control chart.
#'
#' @param ... Parameters to be passed to SPCproperty
#' @param covprob The coverage probability of the adjustment.
#'
#' @return The desired confidence interval, a vector of length 2.
#'
#' @seealso \code{\link{SPCproperty}}
#'
#' @examples
#' # Compute 2-sided CI for the ARL of a CUSUM control chart assuming normality.
#'  X <-  rnorm(100) #observed data
#'  chart <- new("SPCCUSUM",model=SPCModelNormal(Delta=1)) # CUSUM chart with normal observations
#'  SPC2sidedconfint(data=X,nrep=100,covprob=0.95,
#'             property="ARL",chart=chart,params=list(threshold=4))
#'
#' @export
SPC2sidedconfint <- function(covprob=0.9,...){
    sort(SPCproperty(covprob=c((1-covprob)/2,1-(1-covprob)/2),...)@res)
}

