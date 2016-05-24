#' Fit a density surface model to segment-specific estimates of abundance
#' or density.
#'
#' Fits a density surface model (DSM) to detection adjusted counts from a
#' spatially-referenced distance sampling analysis. \code{\link{dsm}} takes
#' observations of animals, allocates them to segments of line (or strip
#' transects) and optionally adjusts the counts based on detectability using a
#' supplied detection function model. A generalized additive model, generalized
#' mixed model or generalised is then used to model these adjusted counts based
#' on a formula involving environmental covariates.
#'
#' The response (LHS of `formula`) can be one of the following:
#' \tabular{ll}{
#'   \code{n}, \code{count}, \code{N}, \code{abundance} \tab count in each segment\cr
#'   \code{Nhat}, \code{abundance.est} \tab estimated abundance per segment, estimation is via a Horvitz-Thompson estimator. This should be used when there are covariates in the detection function.\cr
#'   \code{presence} \tab interpret the data as presence/absence (remember to change the \code{family} argument to \code{binomial()})\cr
#'   \code{D}, \code{density}, \code{Dhat}, \code{density.est} \tab density per segment\cr
#'  }
#'
#' The offset used in the model is dependent on the response:
#' \tabular{ll}{
#'   count \tab area of segment multiplied by average probability of detection in the segment\cr
#'   estimated count \tab area of the segment\cr
#'   presence \tab zero\cr
#'   density \tab zero\cr
#'  }
#'
#' In the latter two cases (density and presence estimation) observations can be weighted by segment areas via the \code{weights=} argument. By default (\code{weights=NULL}), when density or presence are estimated the weights are set to the segment areas (using \code{segment.area} or by calculating \code{2*}(strip width)\code{*Effort}) Alternatively \code{weights=1} will set the weights to all be equal.  A third alternative is to pass in a vector of length equal to the number of segments, containing appropriate weights.
#'
#' @section Units:
#'
#' It is often the case that distances are collected in metres and segment lengths are recorded in kilometres. \code{dsm} allows you to provide a conversation factor (\code{convert.units}) to multiply the areas by. For example: if distances are in metres and segment lengths are in kilometres setting \code{convert.units=1000} will lead to the analysis being in metres. Setting \code{convert.units=1/1000} will lead to the analysis being in kilometres. The conversion factor will be applied to `segment.area` if that is specified.
#'
#' @section Large models:
#'
#' For large models, \code{engine="bam"} with \code{method="fREML"} may be useful. Models specified for \code{bam} should be as \code{gam}. READ \code{\link{bam}} before using this option; this option is considered EXPERIMENTAL at the moment. In particular note that the default basis choice (thin plate regression splines) will be slow and that in general fitting is less stable than when using \code{gam}. For negative binomial response, theta must be specified when using \code{bam}.
#'
#'
#' @param formula formula for the surface. This should be a valid \code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}} formula. See "Details", below, for how to define the response.
#' @param ddf.obj result from call to \code{\link{ddf}} or \code{\link[Distance]{ds}}. If \code{ddf.obj} is \code{NULL} then strip transects are assumed.
#' @param segment.data segment data, see \code{\link{dsm-data}}.
#' @param observation.data observation data, see \code{\link{dsm-data}}.
#' @param engine which fitting engine should be used for the DSM (\code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}}/\code{\link{bam}}).
#' @param convert.units conversion factor to multiply the area of the segments by. See 'Units' below.
#' @param family response distribution (popular choices include \code{\link{quasipoisson}}, \code{\link{Tweedie}} and \code{\link{negbin}}). Defaults to \code{quasipossion}.
#' @param group if \code{TRUE} the abundance of groups will be calculated rather than the abundance of individuals. Setting this option to \code{TRUE} is equivalent to setting the size of each group to be 1.
#' @param control the usual \code{control} argument for a \code{gam}; \code{keepData} must be \code{TRUE} for variance estimation to work.
#' @param availability an availability bias used to scale the counts/estimated  counts by. If we have \code{N} animals in a segment, then \code{N/availability} will be entered into the model. Uncertainty in the availability is not handled at present.
#' @param gamma parameter to \code{gam()} set to a value of 1.4 (from advice in Wood (2006)) such that the \code{gam()} is inclined to not 'overfit' when GCV is used to select the smoothing parameter (ignored for REML, see \code{link{gam}} for further details).
#' @param strip.width if \code{ddf.obj}, above, is \code{NULL}, then this is where the strip width is specified (i.e. for a strip transect survey). This is sometimes (and more correctly) referred to as the half-width, i.e. right truncation minus left truncation.
#' @param segment.area if `NULL` (default) segment areas will be calculated by multiplying the `Effort` column in `segment.data` by the (right minus left) truncation distance for the `ddf.obj` or by `strip.width`. Alternatively a vector of segment areas can be provided (which must be the same length as the number of rows in `segment.data`) or a character string giving the name of a column in `segment.data` which contains the areas. If \code{segment.area} is specified it takes precident.
#' @param weights weights for each observation used in model fitting. The default, \code{weights=NULL}, weights each observation by its area (see Details). Setting a scalar value (e.g. \code{weights=1}) all observations are equally weighted.
#' @param \dots anything else to be passed straight to \code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}}/\code{\link{bam}}.
#' @return a \code{\link{glm}}/\code{\link{gam}}/\code{\link{gamm}} object, with an additional element, \code{ddf} which holds the detection function object.
#' @author David L. Miller
# @seealso
#' @references Hedley, S. and S. T. Buckland. 2004. Spatial models for line transect sampling. JABES 9:181-199.
#'
#' Miller, D. L., Burt, M. L., Rexstad, E. A., Thomas, L. (2013), Spatial models for distance sampling data: recent developments and future directions. Methods in Ecology and Evolution, 4: 1001-1010. doi: 10.1111/2041-210X.12105 (Open Access, available at http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12105/abstract)
#'
#' Wood, S.N. 2006. Generalized Additive Models: An Introduction with R. CRC/Chapman & Hall.
#' @export
#' @importFrom stats quasipoisson
#' @importFrom utils packageVersion
#'
#' @examples
#' library(Distance)
#' library(dsm)
#'
#' # load the Gulf of Mexico dolphin data (see ?mexdolphins)
#' data(mexdolphins)
#'
#' # fit a detection function and look at the summary
#' hr.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance),
#'                key = "hr", adjustment = NULL)
#' summary(hr.model)
#'
#' # fit a simple smooth of x and y
#' mod1 <- dsm(N~s(x,y), hr.model, mexdolphins$segdata, mexdolphins$obsdata)
#' summary(mod1)
#'
#' # create an offset (in metres)
#' # each prediction cell is 444km2
#' off.set <- 444*1000*1000
#'
#' # predict over a grid
#' mod1.pred <- predict(mod1, mexdolphins$preddata, off.set)
#'
#' # calculate the predicted abundance over the grid
#' sum(mod1.pred)
#'
#' # plot the smooth
#' plot(mod1)
dsm <- function(formula, ddf.obj, segment.data, observation.data,
                engine="gam", convert.units=1,
                family=quasipoisson(link="log"), group=FALSE, gamma=1.4,
                control=list(keepData=TRUE), availability=1, strip.width=NULL,
                segment.area=NULL, weights=NULL, ...){

  # if we have a model fitted using Distance, then just pull out the
  # ddf component
  if(!is.null(ddf.obj)){
    if(all(class(ddf.obj)=="dsmodel")){
      ddf.obj <- ddf.obj$ddf
    }
  }

  ## check the formula
  response <- as.character(formula)[2]
  possible.responses <- c("D", "density", "Dhat", "density.est",
                          "N", "abundance", "count", "n",
                          "Nhat", "abundance.est",
                          "presence")
  if(!(response %in% possible.responses)){
    stop(paste("Model must be one of:",
               paste(possible.responses,collapse=", ")))
  }

  ## check that the necessary columns exist in the data
  # NB this doesn't return anything just throws an error if something
  #    bad happens
  check.cols(ddf.obj, segment.data, observation.data, strip.width,segment.area)

  ## build the data
  dat <- make.data(response, ddf.obj, segment.data, observation.data,
                   group, convert.units, availability, strip.width,
                   segment.area)

  ## if we are not modelling density/presence, then add in the offset
  ##  to the formula
  if(!(response %in% c("D","density","Dhat","density.est","presence"))){
    formula <- as.formula(paste(c(as.character(formula)[c(2,1,3)],
                                "+ offset(off.set)"), collapse=""))
  }else{
    # set the weights if we are doing density or presence estimation
    if(is.null(weights)){
      weights <- dat$segment.area
    }else if(length(weights)==1){
      weights <- rep(1, nrow(dat))
    }
  }

  # if we're using a gamm engine
  if(engine == "gamm"){

    # warn if using an old version of mgcv
    mgcv.version <- as.numeric(strsplit(as.character(packageVersion("mgcv")),
                                        "\\.")[[1]])
    if(mgcv.version[1]<1 | (mgcv.version[2]<7 |
                            (mgcv.version[2]==7 & mgcv.version[3]<24))){
      message("You are using mgcv version < 1.7-24, please update to at least 1.7-24 to avoid fitting problems.")
    }

    # unsupported
    control$keepData <- NULL
  }

  ## run the engine
  if(engine %in% c("gam","bam","glm","gamm")){
    args <- list(formula = formula,
                 family  = family,
                 data    = dat,
                 gamma   = gamma,
                 weights = weights,
                 control = control, ...)

    fit <- withCallingHandlers(do.call(engine, args),
                               warning=matrixnotposdef.handler)
  }else{
    stop("engine must be one of 'gam', 'gamm', 'bam' or 'glm'")
  }


  ## save knots
  if("knots" %in% names(match.call())){
    fit$knots <- get(as.character(match.call()$knots))
  }

  ## add the detection function object into the gam/gamm/glm object
  if(engine == "gamm"){
    fit$gam$ddf <- ddf.obj
    fit$gam$data <- dat
    fit$gam$gamma <- gamma
    # yucky way to get dsm.var/dsm.var.prop to work because gamm()
    #  doesn't store the call()
    fit$gam$gamm.call.list <- list(formula = formula,
                                   family  = family,
                                   data    = dat,
                                   gamma   = gamma,
                                   control = control)
  }else{
    fit$ddf <- ddf.obj
    fit$gamma <- gamma
  }

  class(fit) <- c("dsm",class(fit))

  ## return the model
  return(fit)

}
