#' Fit detection functions and calculate abundance from line or point transect data
#'
#' This function fits detection functions to line or point transect data and then (provided that survey information is supplied) calculates abundance and density estimates. The examples below illustrate some basic types of analysis using \code{ds()}.
#'
#' @param data a \code{data.frame} containing at least a column called
#'        \code{distance}. NOTE! If there is a column called \code{size} in
#'        the data then it will be interpreted as group/cluster size, see the
#'        section "Clusters/groups", below. One can supply data as a "flat file"
#'        and not supply \code{region.table}, \code{sample.table} and
#'        \code{obs.table}, see "Data format", below and \code{\link{flatfile}}.
#' @param truncation either truncation distance (numeric, e.g. 5) or percentage (as a string, e.g. "15\%"). Can be supplied as a \code{list} with elements \code{left} and \code{right} if left truncation is required (e.g. \code{list(left=1,right=20)} or \code{list(left="1\%",right="15\%")} or even \code{list(left="1",right="15\%")}).
#' By default for exact distances the maximum observed distance is used as the right truncation. When the data is binned, the right truncation is the largest bin end point. Default left truncation is set to zero.
#' @param transect indicates transect type "line" (default) or "point".
#' @param formula formula for the scale parameter. For a CDS analysis leave this as its default \code{~1}.
#' @param key key function to use; "hn" gives half-normal (default), "hr" gives hazard-rate and "unif" gives uniform. Note that if uniform key is used, covariates cannot be included in the model.
#' @param adjustment adjustment terms to use; "cos" gives cosine (default),
#'        "herm" gives Hermite polynomial and "poly" gives simple polynomial.
#'        "cos" is recommended. A value of \code{NULL} indicates that no
#'        adjustments are to be fitted.
#' @param order orders of the adjustment terms to fit (as a vector/scalar), the
#'        default value (\code{NULL}) will select via AIC up to order 5. If a single number is given, that number is expanded to be \code{seq(term_min, order, by=1)} where \code{term.min} is the appropriate minimum order for this type of adjustment. For cosine
#'        adjustments, valid orders are integers greater than 2 (except when a
#'        uniform key is used, when the minimum order is 1). For Hermite
#'        polynomials, even integers equal or greater than 2 are allowed and for
#'        simple polynomials even integers equal or greater than 2 are allowedi (though note these will be multiplied by 2, see Buckland et al, 2001 for details on their specification). By default, AIC selection will try up to 5 adjustments, beyond that you must specify these manually, e.g. \code{order=2:6} and perform your own AIC selection.
#' @param scale the scale by which the distances in the adjustment terms are
#'        divided. Defaults to "width", scaling by the truncation
#'        distance. If the key is uniform only "width" will be used. The other
#'        option is "scale": the scale parameter of the detection
#' @param cutpoints if the data are binned, this vector gives the cutpoints of
#'        the bins. Ensure that the first element is 0 (or the left truncation
#'        distance) and the last is the distance to the end of the furthest bin.
#'        (Default \code{NULL}, no binning.)
#'        Note that if \code{data} has columns \code{distbegin} and
#'        \code{distend} then these will be used as bins if \code{cutpoints}
#'        is not specified. If both are specified, \code{cutpoints} has
#'        precedence.
#' @param monotonicity should the detection function be constrained for monotonicity weakly (\code{"weak"}), strictly (\code{"strict"}) or not at all (\code{"none"} or \code{FALSE}). See Montonicity, below. (Default \code{"strict"}). By default it is on for models without covariates in the detection function, off when covariates are present.
#' @param dht.group should density abundance estimates consider all groups to be
#'        size 1 (abundance of groups) \code{dht.group=TRUE} or should the
#'        abundance of individuals (group size is taken into account),
#'        \code{dht.group=FALSE}. Default is \code{FALSE} (abundance of
#'        individuals is calculated).
#' @param region.table \code{data.frame} with two columns:
#'        \tabular{ll}{ \code{Region.Label} \tab label for the region\cr
#'                     \code{Area} \tab area of the region\cr}
#'        \code{region.table} has one row for each stratum. If there is no
#'        stratification then \code{region.table} has one entry with \code{Area}
#'        corresponding to the total survey area.
#' @param sample.table \code{data.frame} mapping the regions to the samples (
#'        i.e. transects). There are three columns:
#'        \tabular{ll}{\code{Sample.Label} \tab label for the sample\cr
#'                     \code{Region.Label} \tab label for the region that the
#'                          sample belongs to.\cr
#'                     \code{Effort} \tab the effort expended in that sample
#'                          (e.g. transect length).\cr}
#' @param obs.table \code{data.frame} mapping the individual observations
#'        (objects) to regions and samples. There should be three columns:
#'        \tabular{ll}{\code{object} \tab \cr
#'                     \code{Region.Label} \tab label for the region that the
#'                          sample belongs to.\cr
#'                     \code{Sample.Label} \tab label for the sample\cr}
#' @param convert.units conversion between units for abundance estimation,
#'        see "Units", below. (Defaults to 1, implying all of the units are
#'        "correct" already.)
#'
#' @param method optimization method to use (any method usable by
#'        \code{\link{optim}} or \pkg{optimx}). Defaults to
#'        "nlminb".
#'
#' @param debug.level print debugging output. 0=none, 1-3 increasing level of
#'        debugging output.
#'
#' @param quiet surpress non-essential messages (useful for bootstraps etc).
#'              Default value FALSE.
#'
#' @param initial.values a \code{list} of named starting values, see
#'        \code{\link{mrds-opt}}. Only allowed when AIC term selection is not used.
#'
#' @return a list with elements:
#'        \tabular{ll}{\code{ddf} \tab a detection function model object.\cr
#'                     \code{dht} \tab abundance/density information (if survey
#'                      region data was supplied, else \code{NULL}).}
#'
#' @section Details:
#'
#' If abundance estimates are required then the \code{data.frame}s \code{region.table} and \code{sample.table} must be supplied. If \code{data} does not contain the columns \code{Region.Label} and \code{Sample.Label} thenthe \code{data.frame} \code{obs.table} must also be supplied. Note that stratification only applies to abundance estimates and not at the detection function level.
#'
#' @section Clusters/groups:
#'  Note that if the data contains a column named \code{size} and \code{region.table}, \code{sample.table} and \code{obs.table} are supplied, cluster size will be estimated and density/abundance will be based on a clustered analsis of the data. Setting this column to be \code{NULL} will perform a non-clustred analysis (for example if "\code{size}" means something else in your dataset).
#'
#' @section Truncation:
#' The right truncation point is by default set to be largest observed distance or bin end point. This is a default will not be appropriate for all data and can often be the cause of model convergence failures. It is recommended that one plots a histogram of the observed distances prior to model fitting so as to get a feel for an appropriate truncation distance. (Similar arguments go for left truncation, if appropriate). Buckland et al (2001) provide guidelines on truncation.
#'
#' When specified as a percentage, the largest \code{right} and smallest \code{left} percent distances are discarded. Percentages cannot be supplied when using binned data.
#'
#'  @section Binning: Note that binning is performed such that bin 1 is all distances greater or equal to cutpoint 1 (>=0 or left truncation distance) and less than cutpoint 2. Bin 2 is then distances greater or equal to cutpoint 2 and less than cutpoint 3 and so on.
#'
#' @section Monotonicity: When adjustment terms are used, it is possible for the detection function to not always decrease with increasing distance. This is unrealistic and can lead to bias. To avoid this, the detection function can be constrained for monotonicity (and is by default for detection functions without covariates).
#'
#'  Monotonicity constraints are supported in a similar way to that described in Buckland et al (2001). 20 equally spaced points over the range of the detection function (left to right truncation) are evaluated at each round of the optimisation and the function is constrained to be either always less than it's value at zero (\code{"weak"}) or such that each value is less than or equal to the previous point (monotonically decreasing; \code{"strict"}). See also \code{\link{check.mono}} in \code{mrds}.
#'
#' Even with no monotonicity constraints, checks are still made that the detection function is monotonic, see \code{\link{check.mono}}.
#'
# THIS IS STOLEN FROM mrds, sorry Jeff!
#' @section Units:
#'  In extrapolating to the entire survey region it is important that
#'  the unit measurements be consistent or converted for consistency.
#'  A conversion factor can be specified with the \code{convert.units}
#'  variable.  The values of \code{Area} in \code{region.table}, must be made
#'  consistent with the units for \code{Effort} in \code{sample.table} and the
#'  units of \code{distance} in the \code{data.frame} that was analyzed.  It is
#'  easiest if the units of \code{Area} are the square of the units of
#'  \code{Effort} and then it is only necessary to convert the units of
#'  \code{distance} to the units of \code{Effort}. For example, if \code{Effort}
#'   was entered in kilometers and \code{Area} in square kilometers and
#'  \code{distance} in meters then using \code{convert.units=0.001} would
#'  convert meters to kilometers, density would be expressed in square
#'  kilometers which would then be consistent with units for \code{Area}.
#'  However, they can all be in different units as long as the appropriate
#'  composite value for \code{convert.units} is chosen.  Abundance for a survey
#'  region can be expressed as: \code{A*N/a} where \code{A} is \code{Area} for
#'  the survey region, \code{N} is the abundance in the covered (sampled)
#'  region, and \code{a} is the area of the sampled region and is in units of
#'  \code{Effort * distance}.  The sampled region \code{a} is multiplied by
#'   \code{convert.units}, so it should be chosen such that the result is in
#'  the same units as \code{Area}.  For example, if \code{Effort} was entered
#'  in kilometers, \code{Area} in hectares (100m x 100m) and \code{distance}
#'  in meters, then using \code{convert.units=10} will convert \code{a} to
#'  units of hectares (100 to convert meters to 100 meters for distance and
#'  .1 to convert km to 100m units).
#'
#'  @section Data format: One can supply \code{data} only to simply fit a detection function. However, if abundance/density estimates are necessary further information is required. Either the \code{region.table}, \code{sample.table} and \code{obs.table} \code{data.frame}s can be supplied or all data can be supplied as a "flat file" in the \code{data} argument. In this format each row in data has additional information that would ordinarily be in the other tables. This usually means that there are additional columns named: \code{Sample.Label}, \code{Region.Label}, \code{Effort} and \code{Area} for each observation. See \code{\link{flatfile}} for an example.
#'
#' @author David L. Miller
#' @seealso \code{\link{flatfile}}
#' @export
#'
#' @importFrom stats quantile as.formula
#' @references
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L., and Thomas, L. (2001). Distance Sampling. Oxford University Press. Oxford, UK.
#'
#' Buckland, S.T., Anderson, D.R., Burnham, K.P., Laake, J.L., Borchers, D.L., and Thomas, L. (2004). Advanced Distance Sampling. Oxford University Press. Oxford, UK.
#'
#' @examples
#'
#' # An example from mrds, the golf tee data.
#' library(Distance)
#' data(book.tee.data)
#' tee.data<-book.tee.data$book.tee.dataframe[book.tee.data$book.tee.dataframe$observer==1,]
#' ds.model <- ds(tee.data,4)
#' summary(ds.model)
#' plot(ds.model)
#'
#' \dontrun{
#' # same model, but calculating abundance
#' # need to supply the region, sample and observation tables
#' region <- book.tee.data$book.tee.region
#' samples <- book.tee.data$book.tee.samples
#' obs <- book.tee.data$book.tee.obs
#'
#' ds.dht.model <- ds(tee.data,4,region.table=region,
#'              sample.table=samples,obs.table=obs)
#' summary(ds.dht.model)
#'
#' # specify order 2 cosine adjustments
#' ds.model.cos2 <- ds(tee.data,4,adjustment="cos",order=2)
#' summary(ds.model.cos2)
#'
#' # specify order 2 and 3 cosine adjustments, turning monotonicity
#' # constraints off
#' ds.model.cos23 <- ds(tee.data,4,adjustment="cos",order=c(2,3),
#'                    monotonicity=FALSE)
#' # check for non-monotonicity -- actually no problems
#' check.mono(ds.model.cos23$ddf,plot=TRUE,n.pts=100)
#'
#' # include both a covariate and adjustment terms in the model
#' ds.model.cos2.sex <- ds(tee.data,4,adjustment="cos",order=2,
#'                         monotonicity=FALSE, formula=~as.factor(sex))
#' # check for non-monotonicity -- actually no problems
#' check.mono(ds.model.cos2.sex$ddf,plot=TRUE,n.pts=100)
#'
#' # truncate the largest 10% of the data and fit only a hazard-rate
#' # detection function
#' ds.model.hr.trunc <- ds(tee.data,truncation="10%",key="hr",adjustment=NULL)
#' summary(ds.model.hr.trunc)
#'}
#'
ds <- function(data, truncation=ifelse(is.null(cutpoints),
                                     ifelse(is.null(data$distend),
                                            max(data$distance),
                                            max(data$distend)),
                                     max(cutpoints)),
             transect=c("line","point"),
             formula=~1, key=c("hn","hr","unif"),
             adjustment=c("cos","herm","poly"),
             order=NULL, scale=c("width","scale"),
             cutpoints=NULL, dht.group=FALSE,
             monotonicity=ifelse(formula==~1,
                                 "strict",
                                 "none"),
             region.table=NULL, sample.table=NULL, obs.table=NULL,
             convert.units=1, method="nlminb", quiet=FALSE, debug.level=0,
             initial.values=NULL){

  # this routine just creates a call to mrds, it's not very exciting
  # or fancy, it does do a lot of error checking though


  # truncation
  if(is.null(truncation)){
    stop("Please supply truncation distance or percentage.")
  }else if(any(unlist(lapply(truncation,is.character))) &
           (!is.null(cutpoints) |
            any(c("distbegin","distend") %in% colnames(data))
           )){
    stop("Truncation cannot be supplied as a percentage with binned data")
  }else{
    # if we have left truncation too...
    if(is.list(truncation)){
      if((any(names(truncation)=="left") &
          any(names(truncation)=="right")) &
         length(truncation)==2){

        # check for each of left and right that we have % or distance...
        # left
        if(is.double(truncation$left) & length(truncation$left)==1){
          left <- truncation$left
        }else if(is.character(truncation$left) & length(truncation$left)==1){
          # % string to number
          truncation$left <- as.numeric(sub("%","",truncation$left))
          left <- quantile(data$distance,probs=(truncation$left/100))
        }else{
          stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
        }
        # right
        if(is.double(truncation$right) & length(truncation$right)==1){
          width <- truncation$right
        }else if(is.character(truncation$right) & length(truncation$right)==1){
          # % string to number
          truncation$right <- as.numeric(sub("%","",truncation$right))
          width <- quantile(data$distance,probs=1-(truncation$right/100))
        }else{
          stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
        }
      }else{
        stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
      }

    # just right truncation
    }else if(is.double(truncation) & length(truncation)==1){
      width <- truncation
      left <- NULL
    }else if(is.character(truncation) & length(truncation)==1){
      # % string to number
      truncation <- as.numeric(sub("%","",truncation))
      width <- quantile(data$distance,probs=1-(truncation/100))
      left <- NULL
    }else{
      stop("Truncation must be supplied as a single number/string or a list with elements \"left\" and \"right\".")
    }
  }

  # check the data, format into the correct tables if we have a flat file
  data <- checkdata(data, region.table, sample.table, obs.table, formula)
  region.table <- data$region.table
  sample.table <- data$sample.table
  obs.table    <- data$obs.table
  data         <- data$data

  ### binning
  if(is.null(cutpoints)){
    if(any(names(data)=="distend") & any(names(data)=="distbegin")){
      message("Columns \"distbegin\" and \"distend\" in data: performing a binned analysis...")
      binned <- TRUE
      breaks <- sort(unique(c(data$distend,data$distbegin)))
    }else{
      binned <- FALSE
      breaks <- NULL
    }
  }else{
    # make sure that the first bin starts 0 or left
    if(!is.null(left)){
      if(cutpoints[1]!=left){
        stop("The first cutpoint must be 0 or the left truncation distance!")
      }
    }else if(cutpoints[1]!=0){
      stop("The first cutpoint must be 0 or the left truncation distance!")
    }

    # remove distbegin and distend if they already exist
    if(any(names(data)=="distend") & any(names(data)=="distbegin")){
      message("data already has distend and distbegin columns, removing them and appling binning as specified by cutpoints.")
      data$distend <- NULL
      data$distbegin <- NULL
    }
    # send off to create.bins to make the correct columns in data
    data <- create.bins(data,cutpoints)
    binned <- TRUE
    breaks <- cutpoints
  }

  # transect type
  point <- switch(match.arg(transect),
                  "line"=FALSE,
                  "point"=TRUE,
                  stop("Only \"point\" or \"line\" transects may be supplied.")
                 )

  # key and adjustments
  key <- match.arg(key)
  # keep the name for the key function
  key.name <- switch(key,
                     hn="half-normal",
                     hr="hazard-rate",
                     unif="uniform"
                    )

  # no uniform key with no adjustments
  if(is.null(adjustment) & key=="unif"){
    stop("Can't use uniform key with no adjustments.")
  }
  # no covariates with uniform
  if((as.formula(formula)!=~1) & key=="unif"){
    stop("Can't use uniform key with covariates.")
  }
  # uniform key must use width scaling
  scale <- match.arg(scale)
  if(key=="unif"){
    scale <- "width"
  }

  # check we have an allowed adjustment
  if(!is.null(adjustment)){
    adjustment <- match.arg(adjustment)
  }

  # if the user supplied order=0, that's equivalent to adjustment=NULL
  if(!is.null(order) & all(order==0)){
    adjustment <- NULL
  }

  if(!is.null(adjustment)){

    if(!is.null(order)){
      aic.search <- FALSE
      if(any(order != ceiling(order))){
          stop("Adjustment orders must be integers.")
      }

      #if(formula != ~1){
      #  stop("Cannot use both adjustments and covariates, choose one!")
      #}

      # check for each adjustment type
      order <- sort(order)
      if(adjustment=="poly"){
        if(any(order/2 != ceiling(order/2))){
          stop("Adjustment orders must be even for Hermite and simple polynomials.")
        }
      }
      if((adjustment=="herm" | adjustment=="cos") & key!="unif"){
        if(any(order==1)){
          stop("Adjustment orders for Hermite polynomials and cosines must start at 2.")
        }
      }

      # if a single number is provided do adjmin:order
      if(length(order)==1){
        # this is according to p. 47 of IDS.
        if(adjustment=="poly"){
          order <- 1:order
        }else{
          order <- 2:order
        }

        # for Fourier...
        if(key=="unif" & adjustment=="cos"){
          order <- 1:order
        }

        if(adjustment=="herm" | adjustment=="poly"){
          order <- 2*order
          order <- order[order<=2*order]
        }
      }


    }else{

      # if there are covariates then don't do the AIC search
      if(formula != ~1){
        aic.search <- FALSE
        message("Cannot perform AIC adjustment term selection when covariates are used.")
      }else{
      # otherwise go ahead and set up the candidate adjustment orders
        aic.search <- TRUE
        max.order <- 5

        # this is according to p. 47 of IDS.
        if(adjustment=="poly"){
          order <- seq(1,max.order)
        }else{
          order <- seq(2,max.order)
        }

        # for Fourier...
        if(key=="unif" & adjustment=="cos"){
          order <- c(1,order)
        }

        if(adjustment=="herm" | adjustment=="poly"){
          order <- 2*order
          order <- order[order<=2*max.order]
        }

      }
    }

    # keep the name for the adjustments
    adj.name <- switch(adjustment,
                       cos="cosine",
                       herm="Hermite",
                       poly="simple polynomial"
                      )

  }else{
    aic.search <- FALSE
  }


  # monotonicity
  if(is.logical(monotonicity)){
    if(!monotonicity){
      mono <- FALSE
      mono.strict <- FALSE
    }
  }else if(monotonicity=="none"){
    mono <- FALSE
    mono.strict <- FALSE
  }else if(monotonicity=="weak"){
    mono <- TRUE
    mono.strict <- FALSE
  }else if(monotonicity=="strict"){
    mono <- TRUE
    mono.strict <- TRUE
  }else{
    stop("monotonicity must be one of \"none\", FALSE, \"weak\" or \"strict\".")
  }

  # can't do monotonicity and covariates, fail!
  if(mono & formula!=as.formula("~1")){
    stop("Monotonicity cannot be enforced with covariates.")
  }

  # set up the control options
  control <- list(optimx.method=method, showit=debug.level)

  # if initial values were supplied, pass them on
  if(!is.null(initial.values) & !aic.search){
    control$initial <- initial.values
  }else if(!is.null(initial.values) & aic.search){
    stop("Cannot supply initial values when using AIC term selection")
  }

  ### Actually fit some models here

  # construct the meta data object...
  meta.data <- list(width = width,point = point,binned = binned,
                    mono=mono, mono.strict=mono.strict)
  if(!is.null(left)){
    meta.data$left <- left
  }
  if(binned){
    meta.data$breaks <- breaks
  }

  # if we are doing an AIC-based search then, create the indices for the
  # for loop to work along, else just give the length of the order object
  if(aic.search){
    if(key!="unif"){
      for.ind <- c(0,seq(along=order))
    }else{
      for.ind <- seq(along=order)
    }
    message("Starting AIC adjustment term selection.")
  }else if(!is.null(adjustment)){
    for.ind <- length(order)
  }else{
    for.ind <- 1
  }

  # dummy last model
  last.model<-list(criterion=Inf)

  # loop over the orders of adjustments
  for(i in for.ind){
    # construct model formulae
    # CDS model
    if(formula==as.formula("~1")){
      model.formula <- paste("~cds(key =\"", key,"\", formula = ~1",sep="")
    # MCDS model
    }else{
      model.formula <- paste("~mcds(key = \"",key,"\",",
                                  "formula =~",as.character(formula)[2],sep="")
    }

    # build a message to let the user know what is being fitted
    this.message <- paste("Fitting ",key.name," key function",sep="")

    # adjustments?
    # this handles the case when we have adjustments but are doing AIC search
    # so want to fit a key function alone to begin with.
    if(!is.null(adjustment) & i!=0){
      if(length(order[1:i])==1){
        order.str <- order[1:i]
      }else{
        order.str <- paste("c(",paste(order[1:i],collapse=","),")",sep="")
      }

      model.formula <- paste(model.formula,",",
                           "adj.series=\"",adjustment,
                           "\",adj.order=",order.str,",",
                           "adj.scale=\"",scale,"\"",sep="")

      this.message <- paste(this.message,
                            " with ", adj.name,"(",
                            paste(order[1:i],collapse=","),
                            ") adjustments", sep="")
    }

    model.formula <- paste(model.formula,")",sep="")

    # tell the user what is being fitted
    message(this.message)

    # actually fit a model
    # wrap everything around this so we don't print out a lot of useless
    # stuff...
    model <- suppressPackageStartupMessages(
               suppressWarnings(try(
                                  ddf(dsmodel = as.formula(model.formula),
                                      data = data, method = "ds",
                                      control=control,
                                      meta.data = meta.data),silent=TRUE)))


    # if that worked
    if(any(class(model)!="try-error")){
      if(model$ds$converge==0){

        model$name.message <- sub("^Fitting ","",this.message)

        # need this to get plotting to work!
        model$call$dsmodel <- as.formula(model.formula)

        message(paste("AIC=",round(model$criterion,3)))

        if(aic.search){
          # if this models AIC is worse (bigger) than the last
          # return the last model and stop looking.
          if(model$criterion>last.model$criterion){
            model <- last.model
            message(paste0("\n\n",model$name.message," selected!"))
            break
          }else{
            # otherwise keep this, best model
            last.model <- model
          }
        }
      }else{
        message("  Model failed to converge.")
      }
    }else{
      if(last.model$criterion == Inf & length(last.model)==1){
        message("\n\nAll models failed to fit!\n")
        model <- NULL
      }else{
        message(paste0("\n\nError in model fitting, returning: ",
                       sub("^Fitting ","",last.model$name.message)))
        message(paste0("\n  Error: ",model[1],"\n"))
        model <- NULL
        model <- last.model
      }
      break
    }
  }

  if(is.null(model)){
    stop("No models could be fitted.")
  }

  ## Now calculate abundance/density using dht()
  if(!is.null(region.table) & !is.null(sample.table)){

    # if obs.table is not supplied, then data must have the Region.Label and
    # Sample.Label columns
    if(is.null(obs.table)){
      if(c("Region.Label","Sample.Label") %in% names(data)){
        message("No obs.table supplied but data does not have Region.Label or Sample.Label columns, only estimating detection function.\n")
      }
    }

    # from ?dht:
    # For animals observed in tight clusters, that estimator gives the
    # abundance of groups (group=TRUE in options) and the abundance of
    # individuals is estimated as s_1/p_1 + s_2/p_2 + ... + s_n/p_n, where
    # s_i is the size (e.g., number of animals in the group) of each
    # observation(group=FALSE in options).

    dht.res <- dht(model,region.table,sample.table,obs.table,
                 options=list(#varflag=0,group=TRUE,
                              group=dht.group,
                              convert.units=convert.units),se=TRUE)
  }else{
    # if no information on the survey area was supplied just return
    # the detection function stuff
    dht.res <- NULL

    if(!quiet){
      message("No survey area information supplied, only estimating detection function.\n")
    }
  }

  # construct return object
  ret.obj <- list(ddf = model,
                dht = dht.res)

  # give it some class
  class(ret.obj) <- "dsmodel"

  return(ret.obj)

}
