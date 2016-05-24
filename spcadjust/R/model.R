#' Data Model for SPC charts
#'
#' This is the basic structure for converting observations (data) into
#' updates for control charts. Classes of this type also have the
#' ability to generate new data sets (resampling).
#' 
#' Every element of this class has to consist of a list of the
#' following functions: updates, getcdfupdates, Pofdata,
#' resample, xiofP, which have to be of a specific form.  The arguments generally have the
#' following meaning: \code{xi} denotes the parameter vector needed to create
#' updates for running the chart from observed data, \code{data} is
#' observed data, \code{P} is a data model.
#'
#' \itemize{
#' \item \code{updates(xi,data)}:  Returns updates for the chart using the parameter \code{xi} and the observed data \code{data}.
#' \item \code{Pofdata(data)}:  Estimates a probability model from the data.
#' \item xiofP(P): Computes the parameter \code{xi} needed to compute updates from an (estimated) probability model \code{P}.
#' \item \code{resample(P)}: Generates a new data set from the probability model \code{P}.
#' \item \code{getcdfupdates(P,xi,cadlag=TRUE)}: Returns the cumulative distribution function (CDF) of updates of data generated from the probability model \code{P} and computed using the parameter \code{xi}. The CDF has to be a function of one argument that also accepts vectors. If cadlag is  TRUE then the CDF is right-continuous (i.e. \eqn{F(x)=P(X\le x)}). If cadlag is FALSE then the CDF is left-continuous (i.e. \eqn{F(x) = P(X<x)}).
#' }
#'
#' @name SPCDataModel-class
#' @aliases SPCDataModel
#' @family SPCDataModel
#'
#' @seealso \code{\link{SPCModelNormal}},
#' \code{\link{SPCModelNonpar}},
#' \code{\link{SPCModelNonparCenterScale}}
#' 
#' @exportClass SPCDataModel
setOldClass("SPCDataModel")

#' Parametric Model for Normally Distributed Data with Centering and Scaling
#'
#' Returns a data model for univariate observations using normality
#' assumptions with updates that center and scale the observations and
#' potentially subtract half of a constant \code{Delta}. Subtracting
#' \eqn{Delta/2} is useful for CUSUM charts.
#' 
#' The parameters to the function have the following meaning.
#' \itemize{
#' \item data: a numeric vector.
#' \item xi: a list with two elements:
#'       \itemize{
#'       \item mu: the mean.
#'       \item sd: the standard deviation.
#'       }
#' \item P: a list with three elements:
#'       \itemize{
#'       \item mu: the mean.
#'       \item sd: the standard deviation.
#'       \item m: the number of data points to resample.
#'       }
#' }
#'
#' The main operations are defined as follows:
#' \itemize{
#' \item updates(xi,data): returns the centered and scale version of the data from which \eqn{Delta/2} has been subtracted, i.e.
#'  \deqn{\frac{data-mu-Delta/2}{sd}}{(data-mu-Delta/2)/sd}.
#' \item resample(P): resamples \code{m} new data points from a normal distribution if mean \code{mu} and standard deviation \code{sd}.
#' }
#' 
#' @param Delta Half of this constant is subtracted for updates (before centering and scaling).
#' @return An object of class \code{\link{SPCDataModel}}.
#' @export
SPCModelNormal <- function(Delta=0){
    structure(
        list(
            updates=function(xi,data) (data-xi$mu-Delta/2)/xi$sd,
            Pofdata=function(data){
                list(mu= mean(data), sd= sd(data), m=length(data))
            },
            xiofP=function(P)  P,
            resample=function(P) P$sd*rnorm(P$m)+P$mu,
            getcdfupdates=function(P, xi,cadlag=TRUE) {
                muupd <- (P$mu-xi$mu-Delta/2)/xi$sd
                sdsqupd <- P$sd/xi$sd
                function(x) pnorm(x, mean=muupd, sd=sdsqupd)
            }
            ),
        class="SPCDataModel"
        )
}

#' Generic Model for Nonparametric Resampling
#'
#' Generic model that allows nonparametric resampling (with
#' replacement) of the data. The transformation of data into updates
#' needs to be defined via the arguments.
#' 
#' The parameters to the functions being returned have the following meaning.
#' \itemize{
#' \item data: a numeric vector or a matrix where the rows contain the observations.
#' \item xi: depends on the parameter xiofP.
#' \item P: The \code{data} with no modifications (thus either a numeric vector or a matrix).
#' }
#'
#' The main operations are defined as follows:
#' \itemize{
#' \item resample(P): generates a new data set of the same size by either resampling the values (if the data is a vector) or by resampling the rows of the data matrix (both use resampling with replacement).
#' }
#' 
#' @param updates function that computes updates.
#' @param xiofP function that computes xi given P.
#'
#' @return An object of class SPCDataModel.
#' @export
SPCModelNonpar <- function(updates,xiofP){
    structure(
        list(updates=updates,
             getcdfupdates=function(P, xi,cadlag=TRUE){
                 if (cadlag){
                     ecdf(updates(xi=xi,data=P))
                 }else{ #P(X<x)
                     x <- sort(updates(xi=xi,data=P))
                     vals <- unique(x)
                     count <- cumsum(tabulate(match(x, vals))/length(x))
                     approxfun(vals, c(0,count[-length(count)]),
                               method = "constant", yleft = 0, yright = 1, f = 1, ties = "ordered")
                 }
             },
             Pofdata=function(data) data,
             resample=function(P){
                 if (is.vector(P))
                     sample(P,replace=TRUE)
                 else
                     P[sample.int(dim(P)[1],replace=TRUE),]
             },
             xiofP=xiofP             
             ),
        class="SPCDataModel")
}
#' Nonparametric Resampling with Centering and Scaling
#'
#' Nonparametric resampling of univariate observations. Updates are
#' centered and scaled transformations of the data (with a constant
#' potentially being subtracted).
#'
#' Calls \code{\link{SPCModelNonpar}} to generate the data object, so it only needs to specify the meaning of \code{xi}, the parameter needed to compute updates and the definition of the updates.
#' \itemize{
#' \item xi: a list with two elements:
#'       \itemize{
#'       \item mu: the mean.
#'       \item sd: the standard deviation.
#'       }
#' \item updates(xi,data): returns the centered and scale version of the data from which \eqn{Delta/2} has been subtracted, i.e.
#'  \deqn{\frac{data-mu-Delta/2}{sd}}{(data-mu-Delta/2)/sd}
#' }
#'
#' @param Delta how much to subtract before scaling.
#' 
#' @return An object of class SPCDataModel.
#' 
#' @examples
#' X <- rnorm(1000)
#'
#' #CUSUM chart
#' model <- SPCModelNonparCenterScale(1)
#' chart <- new("SPCCUSUM",model=model)
#' SPCproperty(data=X,nrep=10,property="calARL",
#'             chart=chart,params=list(target=100))
#'
#' #Shewhart chart
#' model <- SPCModelNonparCenterScale(0)
#' chart <- new("SPCCUSUM",model=model)
#' SPCproperty(data=X,nrep=10,property="calARL",
#'             chart=chart,params=list(target=100))
#' @seealso \code{\link{SPCModelNonpar}}
#' @return An object of class SPCDataModel.
#' @export
SPCModelNonparCenterScale <- function(Delta=0){
    SPCModelNonpar(updates=function(xi,data)(data-xi$mu-Delta/2)/xi$sd,
                  xiofP=function(P) list(mu=mean(P),sd=sd(P)))
}
