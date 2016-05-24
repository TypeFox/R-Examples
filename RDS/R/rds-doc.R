
#' This package provides functionality for carrying out estimation
#'     with data collected using Respondent-Driven Sampling. This includes
#'     Heckathorn's RDS-I and RDS-II estimators as well as Gile's Sequential
#'     Sampler estimator.
#' @import ggplot2 reshape2 gridExtra methods
#' @importFrom stats approx coef fft median optim pnorm printCoefmat pt qnorm qt quantile resid residuals rnorm runif sd spline symnum uniroot var vcov xtabs
#' @importFrom graphics abline axis hist legend lines par plot points segments strwidth symbols
#' @importFrom scales hue_pal
#' @docType package
#' @name RDS
#' @useDynLib RDS
NULL


#' A Simulated RDS Data Set
#' @description This is a faux set used to demonstrate RDS functions and analysis.
#' It is used is some simple examples and has categorical variables "X", "Y" and "Z".
#' @docType data
#' @keywords datasets
#' @format An rds.data.frame object
#' @references Gile, Krista J., Handcock, Mark S., 2010 \emph{Respondent-driven Sampling: An Assessment of Current Methodology},  \emph{Sociological Methodology}, 40, 285-327. 
#' @seealso \code{\link{fauxsycamore}}, \code{\link{fauxmadrona}}
#' @examples 
#' data(faux)
#' RDS.I.estimates(rds.data=faux,outcome.variable='X')
#' @name faux
NULL


#' A Simulated RDS Data Set with no seed dependency
#' 
#' This is a faux set used to illustrate how the estimators perform under
#' different populations and RDS schemes.
#' 
#' The population had N=1000 nodes.  In this case, the sample size is 500 so
#' that there is a relatively small sample fraction (50\%). There is homophily
#' on disease status (R=5) and there is differential activity by disease status
#' whereby the infected nodes have mean degree twice that of the uninfected
#' (w=1.8).
#' 
#' In the sampling, the seeds are chosen randomly from the full population, so
#' there is no dependency induced by seed selection.
#' 
#' Each sample member is given 2 uniquely identified coupons to distribute to
#' other members of the target population in their acquaintance.  Further each
#' respondent distributes their coupons completely at random from among those
#' they are connected to.
#' 
#' 
#' 
#' Here are the results for this data set and the sister \code{fauxsycamore}
#' data set:
#' 
#' \tabular{rlllllll}{ \bold{Name} \tab \bold{City} \tab \bold{Type} \tab
#' \bold{Mean} \tab \bold{RDS I (SH)} \tab \bold{RDS II (VH)} \tab \bold{SS}
#' \cr fauxsycamore \tab Oxford\tab seed dependency, 70\% \tab
#' 0.2408 \tab 0.1087 \tab 0.1372 \tab 0.1814\cr fauxmadrona \tab
#' Seattle\tab no seed dependency, 50\% \tab 0.2592 \tab 0.1592 \tab 0.1644
#' \tab 0.1941}
#' 
#' Even with only 50\% sample, the VH is substantially biased , and the SS
#' does much better.
#' 
#' @name fauxmadrona
#' @aliases fauxmadrona fauxmadrona.network
#' @docType data
#' @format An \code{rds.data.frame}
#' @seealso \code{\link{fauxsycamore}}, \code{\link{faux}}
#' @references Gile, Krista J., Handcock, Mark S., 2010 \emph{Respondent-driven
#' Sampling: An Assessment of Current Methodology}, \emph{Sociological
#' Methodology}, 40, 285-327.
#' @source The original network is included as
#' \code{fauxmadrona.network} as a \code{network} object.  \cr The data set
#' also includes the \code{data.frame} of the RDS data set as
#' \code{fauxmadrona}.  \cr Use \code{data(package="RDS")} to get a full list
#' of datasets.
#' 
#' @keywords datasets
NULL




#' A Simulated RDS Data Set with extreme seed dependency
#' 
#' This is a faux set used to demonstrate RDS functions and analysis.  The
#' population had N=715 nodes.  In this case, the sample size is 500 so that
#' there is a relatively large sample fraction (70\%). There is homophily on
#' disease status (R=5) and there is differential activity by disease status
#' whereby the infected nodes have mean degree twice that of the uninfected
#' (w=1.8).
#' 
#' In the sampling the seeds are chosen randomly from the infected population,
#' so there is extreme dependency induced by seed selection.
#' 
#' Each sample member is given 2 uniquely identified coupons to distribute to
#' other members of the target population in their acquaintance.  Further each
#' respondent distributes their coupons completely at random from among those
#' they are connected to.
#' 
#' With 70\% sample, the VH is substantially biased, so the SS (and presumably
#' MA) do much better.  We expect the MA to perform a bit better than the SS.
#' 
#' It is network 702 and its sample from YesYes on mosix. Look for
#' "extract702.R" \cr The original network is included as
#' \code{fauxsycamore.network} as a \code{network} object.  \cr The data set
#' also includes the \code{data.frame} of the RDS data set as
#' \code{fauxsycamore}.  \cr Use \code{data(package="RDS")} to get a full list
#' of datasets.
#' 
#' @name fauxsycamore
#' @aliases fauxsycamore fauxsycamore.network
#' @docType data
#' @format An rds.data.frame plus the original network as a network object
#' @seealso \code{\link{faux}}, \code{\link{fauxmadrona}}
#' @references Gile, Krista J., Handcock, Mark S., 2009.
#' \emph{Respondent-driven Sampling: An Assessment of Current Methodology},
#' \emph{Sociological Methodology}, 40, 285-327.
#' @keywords datasets
NULL




