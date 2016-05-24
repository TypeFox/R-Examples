#' Create detection function object
#'
#' Creates and populates a specific list structure to define a detection
#' function object and its data. The \code{ddfobj} is used throughout the
#' package as a calling argument to various functions.
#'
#' @param model model list with key function and possibly adjustment functions,
#'   scale formula, and shape formula
#' @param xmat model data frame
#' @param meta.data list of options describing data like width, etc
#' @param initial vector of initial values for parameters of the detection
#'   function
#' @return Distance sampling function object list with elements that all can be
#'   null except type: \item{type}{type of detection function
#'   hn,hr,gamma,unif,logistic} \item{xmat}{model data frame}
#'   \item{intercept.only}{TRUE if scale = ~1 and any shape formula =~1}
#'   \item{scale}{sublist with elements (can be NULL i.e., unif key):formula,
#'   parameters, design matrix (dm)} \item{shape}{sublist with elements (power
#'   of hazard rate or gamma) (can be NULL i.e., unif or hn key):formula,
#'   parameters, design matrix (dm)} \item{adjustment}{sublist with elements
#'   (is NULL if no adjustments used):series,order,scale,parameters}
#'   \item{g0}{sublist with elements (not used at present):formula,parameters,
#'   design matrix(dm), link}
#' @note Internal function not meant to be called by user
#' @author Jeff Laake
#' @seealso \code{\link{detfct}}, \code{\link{ddf}}
create.ddfobj <- function(model,xmat,meta.data,initial){

  # Create empty object and get values from cds or mcds function
  ddfobj <- vector("list")
  point <- meta.data$point
  modpaste <- paste(model)
  modelvalues <- eval(parse(text=modpaste[2:length(modpaste)]))

  # Specify key function type
  ddfobj$type <- modelvalues$key

  if(ddfobj$type=="logistic"){
    stop("Logistic detection function has been temporarily disabled")
  }

  if(!ddfobj$type%in%c("gamma","hn","hr","unif","th1","th2")){
    stop("Invalid value for detection key function =",ddfobj$type,
         "  Only hn, hr, gamma, unif, th1, th2 allowed")
  }

  # Set adjustment function options
  if(!is.null(modelvalues$adj.series)){

    if(is.null(modelvalues$adj.order)){
      stop("You must specify adjustment order(s) via adj.order")
    }

    if(is.null(modelvalues$adj.exp)) modelvalues$adj.exp <- FALSE

    ddfobj$adjustment <- list(series = modelvalues$adj.series,
                              order  = modelvalues$adj.order,
                              scale  = modelvalues$adj.scale,
                              exp    = modelvalues$adj.exp)
  }else{
    ddfobj$adjustment <- NULL
  }

  # Assign scale and shape (if any) formulae
  if(ddfobj$type!="unif"){
    if(is.null(modelvalues$formula)){
      ddfobj$scale <- list(formula="~1")
    }else{
      ddfobj$scale <- list(formula=paste(as.character(modelvalues$formula),
                                      collapse=""))
    }
  }

  if(!is.null(modelvalues$shape.formula)){
    ddfobj$shape <- list(formula=paste(as.character(modelvalues$shape.formula),
                                       collapse=""))
  }else if(ddfobj$type%in%c("hr","gamma","th1","th2")){
    ddfobj$shape <- list(formula=~1)
  }else{
    ddfobj$shape <- NULL
  }

  # Create model data frame and design matrices
  ddfobj$xmat <- create.model.frame(xmat,as.formula(ddfobj$scale$formula),
                                    meta.data,as.formula(ddfobj$shape$formula))
  if(ddfobj$type !="unif"){
    ddfobj$scale$dm <- setcov(ddfobj$xmat,ddfobj$scale$formula)
    ddfobj$scale$parameters <- rep(0,ncol(ddfobj$scale$dm))
    # Next determine if scale covariate model is intercept only.
    ddfobj$intercept.only <- FALSE
    if(ddfobj$scale$formula == "~1" &
        (is.null(ddfobj$shape) || ddfobj$shape$formula== "~1")){
      ddfobj$intercept.only <- TRUE
    }
  }else{
    ddfobj$intercept.only <- TRUE
  }

  if(!is.null(ddfobj$shape)){
    ddfobj$shape$dm <- setcov(ddfobj$xmat,ddfobj$shape$formula)
    ddfobj$shape$parameters <- rep(0,ncol(ddfobj$shape$dm))
  }

  # Compute initialvalues unless uniform
  initialvalues <- setinitial.ds(ddfobj,width=meta.data$width,initial,point)

  # Delete columns of dm that end up as NA from initialvalues
  if(!is.null(ddfobj$scale)){
    if(!ddfobj$intercept.only){
      if(any(is.na(initialvalues$scale))){
        warning("Model is not full rank - not all parameters are estimable.")
      }
      ddfobj$scale$dm[,!is.na(initialvalues$scale)]
    }
    ddfobj$scale$parameters <- initialvalues$scale[!is.na(initialvalues$scale)]
  }

  if(!is.null(ddfobj$shape)){
    ddfobj$shape$parameters <- initialvalues$shape
  }

  if(!is.null(ddfobj$adjustment)){
    ddfobj$adjustment$parameters <- initialvalues$adjustment
  }

  # If the model includes both covariates (in scale/shape) and adjustments
  if(!is.null(ddfobj$adjustment) & ddfobj$type!="unif"){
    if(ddfobj$scale$formula != ~1){
      # commented-out line below prevents models with adjustments and
      # covariates being fitted
      #stop("\nAdjustment functions cannot be used with scale covariates")
    }else{
      # shape covariates
      if(!is.null(ddfobj$shape)){
        if(ddfobj$shape$formula != ~1){
          stop("\nAdjustment functions cannot be used with shape covariates")
        }
      }
    }
  }

  ## check that we didn't specify a model that doesn't make sense

  # can't constrain monotonicity with covariates
  # don't use monotonicity unless we have adjustment terms
  if(meta.data$mono){
    if(ddfobj$type !="unif"){
      if(ddfobj$scale$formula!="~1"){
        stop("Covariate models cannot be constrained for monotonicity.")
      }
    }

    if(is.null(ddfobj$adjustment)){
      message("Key only models do not require monotonicity contraints. Not constraining model for monotonicity.")
      meta.data$mono <- FALSE
    }
  }


  return(ddfobj)
}
