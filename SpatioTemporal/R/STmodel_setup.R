######################################
## FUNCTIONS THAT HANDLE COVARIATES ##
######################################
##Functions in this file:
## processLUR        EX:ok
## processST         EX:ok
## processLocation   EX:ok

#############################################################
## Helper functions that process LUR and ST specifications ##
#############################################################
##' Function that create covariate specifications for \code{\link{createSTmodel}},
##' and compare the covariates requested (both geographic and spatio-temporal) with
##' those available in \code{STdata}.
##'
##' Several options exist for \code{LUR.in}
##' \describe{
##'   \item{\code{LUR.in=NULL}}{Only an intercept for all beta-fields.}
##'   \item{\code{LUR.in="all"}}{Use all elements in \code{STdata$covars},
##'     \emph{NOT} recommended.}
##'   \item{\code{LUR.in=list(...)}}{Use different covariates for each, specified
##'     by the different components of the list.}
##'   \item{\code{LUR.in=vector}}{Use the same covariates for all beta-field.}
##' }
##' For the two last options the vector/list-elements can contain either:
##' \describe{
##'   \item{integer}{This will be used as \code{names(STdata$covars)[int]} to
##'     extract a character vector (see below) of covariates.}
##'   \item{character}{The character vector will be used to create a
##'     \code{\link[stats:formula]{formula}} (see below), through: \cr
##'     \code{as.formula(paste("~", paste(unique(chars), collapse="+")),
##'           env=.GlobalEnv)}
##'   }
##'   \item{\code{\link[stats:formula]{formula}}}{The formula will be used as
##'     \code{ \link[stats:model.matrix]{model.matrix}(formula, STdata$covars) }
##'     to create a covariate matrix.}
##' }
##' Setting any element(s) of the list to \code{NULL} implies \emph{only an
##' intercept} for the corresponding temporal trend(s).
##'
##' \code{ST.in} should be a vector specifying the spatio-temporal covariates to
##' use; the vector either give names or layers in \code{STdata$SpatioTemporal}
##' to use, compare \code{character} and \code{integer} options for \code{LUR.in}
##' above.
##' 
##' If covariates are specified using names these should match \cr
##' \code{dimnames(STdata$SpatioTemporal)[[3]]}, unmatched elements are
##' dropped with a warning.
##' 
##' @title Internal Function that do Covariate Selection
##' @param STdata \code{STdata} object with observations, covariates, trends, etc;
##'   see \code{\link{mesa.data.raw}}.
##' @param LUR.in A vector or list indicating which geographic covariates to use.
##' 
##' @return A list of LUR specifications, as \code{\link[stats:formula]{formula}};
##'   or a ST specification as a character vector.
##' 
##' @example Rd_examples/Ex_processLUR.R
##' 
##' @author Johan Lindström
##' @family STmodel functions
##' @export
processLUR <- function(STdata, LUR.in){
  ##check class belonging
  stCheckClass(STdata, "STdata", name="STdata")

  ##find the number of temporal trends (including intercept)
  if( is.null(STdata$trend) ){
    nt <- 1
  }else{
    nt <- dim(STdata$trend)[2]
  }
  
  ##There are a few different options for LUR.in
  if( is.null(LUR.in) ){
    ##only constant
    LUR <- vector("list", nt)
  }else if( is.list(LUR.in) ){
    ##already a list, make sure the length is right
    if( length(LUR.in)!=nt ){
      stop("Length of 'LUR.in' does not match number of temporal trends.")
    }
    LUR <- LUR.in
  }else{
    if( length(LUR.in)==1 && LUR.in=="all" ){
      ##use all covariates
      LUR <- names(STdata$covars)
    }
    ##LUR a single vector, replicate for all trends.
    LUR <- vector("list", nt)
    for(i in 1:nt){
      LUR[[i]] <- LUR.in
    }
  }
  ##convert each element in the list to a formulae
  for(i in 1:nt){
    ##already a formula object, next iteration please
    if( class(LUR[[i]])=="formula" ){
      next
    }
    ##only intercept, ~1 and next iteration please
    if( length(LUR[[i]])==0 ){
      LUR[[i]] <- as.formula("~1", env=.GlobalEnv)
      next
    }
    if( is.numeric(LUR[[i]]) ){
      ##convert numeric to names of covariates
      if( any(LUR[[i]] > dim(STdata$covars)[2]) ){
        warning( sprintf("Unable to find LUR with indecies larger than %d",
                         dim(STdata$covars)[2]) )
        LUR[[i]] <- LUR[[i]][ LUR[[i]] <= dim(STdata$covars)[2] ]
        ##safe guard if we drop all the covariates, revert to constant
        if( length(LUR[[i]])==0 ){
          LUR[[i]] <- as.formula("~1", env=.GlobalEnv)
          next
        }
      }
      LUR[[i]] <- names(STdata$covars)[ LUR[[i]] ]
    }
    if( is.character(LUR[[i]]) ){
      ##convert character string to formula
      LUR[[i]] <- as.formula(paste("~", paste(unique(LUR[[i]]), collapse="+")),
                             env=.GlobalEnv)
    }else{
      stop(paste("Unknown LUR specification:", LUR[[i]]))
    }
  }##for(i in 1:nt)
  ##validate that we only request variables that are available
  for(i in 1:nt){
    stCheckFields(STdata$covars, all.vars(LUR[[i]]), name="STdata$covars")
  }##for(i in 1:nt)
  ##add names for the components
  names(LUR) <- c("const",names(STdata$trend)[-which(names(STdata$trend)=="date")])
  ##return LUR specification
  return(LUR)
}##function processLUR

##' @param ST.in A vector indicating which spatio-temporal covariates to use.
##' @rdname processLUR
##' @export
processST <- function(STdata, ST.in){
  ##check class belonging
  stCheckClass(STdata, "STdata", name="STdata")

  ##no ST-covariates, return empty vector.
  if( length(ST.in)==0 ){
    return( character(0) )
  }
  if( is.null(STdata$SpatioTemporal) || dim(STdata$SpatioTemporal)[3]==0 ){
    warning("No spatio-temporal covariate exists in STdata.")
    return( character(0) )
  }

  if( is.numeric(ST.in) ){
    ##pick by number
    if( any(ST.in > dim(STdata$SpatioTemporal)[3]) ){
      warning( sprintf("Unable to find SpatioTemporal with indecies larger than %d.",
                       dim(STdata$SpatioTemporal)[3]) )
      ST.in <- ST.in[ ST.in <= dim(STdata$SpatioTemporal)[3] ]
    }
    ST.in <- dimnames( STdata$SpatioTemporal )[[3]][ST.in]
  }else if( is.character(ST.in) ){
    ##pick by name
    if( length(ST.in)==1 && ST.in=="all" ){
      ##use all covariates
      ST.in <- dimnames( STdata$SpatioTemporal )[[3]]
    }
  }else{
    stop( paste("Unknown ST specification:", ST.in) )
  }

  ##retain only matched names
  if( !all(ST.in %in% dimnames( STdata$SpatioTemporal )[[3]] )){
    tmp <- ST.in[ !(ST.in %in% dimnames( STdata$SpatioTemporal )[[3]]) ]
    warning( paste("Unmatched ST-covariates:", paste(tmp, collapse=", ")) )
  }
  ST.in <- ST.in[ ST.in %in% dimnames(STdata$SpatioTemporal)[[3]] ]
  ##ensure no duplicates
  return( unique(ST.in) )
}##function processST

##########################################################
## Helper function that process location specifications ##
##########################################################

##' Function that creates a data.frame of locations (and auxillirary information)
##' from \code{STdata$covars}, used by \code{\link{createSTmodel}}.
##'
##' The \code{locations} list specifies what should go in the locations data.frame,
##' in addition to thing listed below \code{STdata$covars$ID} is always added.
##' Each of the fields below should contain names (as character) of columns in
##' \code{STdata$covars}
##' \describe{
##'   \item{coords}{The x,y-coordinates for monitors}
##'   \item{coords.beta,coords.nu}{Alternative x,y-coordinates for monitors,
##'     used when computing distance-matrices for the beta- and nu-fields. Allows
##'     the use of non-stationary covariance structures thourgh the deformation
##'     method of Damian (2003), given a precomputed deformation.}
##'   \item{long.lat}{The long,lat-coordinates for monitors}
##'   \item{others}{Additional fields in \code{STdata$covars} that should be added
##'     to the location data.frame}
##' }
##' 
##' @title Internal Function that Extracts Locations
##' @param STdata \code{STdata} object with observations, covariates, trends, etc;
##'   see \code{\link{mesa.data.raw}}.
##' @param locations A list specifying which fields in \code{STdata$covars} that
##'   should be used for what in the location data.frame, see details.
##' 
##' @return A data.frame with location information for all the sites.
##'
##' @references
##' D. Damian, P. D. Sampson, P. Guttorp. (2003) Variance modeling for
##'   nonstationary processes with temporal replications. J. Geophys. Res.:
##'   D24(108)
##' 
##' @example Rd_examples/Ex_processLocation.R
##' 
##' @author Johan Lindström
##' @family STmodel functions
##' @export
processLocation <- function(STdata, locations){
  ##check class belonging
  stCheckClass(STdata, "STdata", name="STdata")
  ##check for fields
  stCheckFields(locations, "coords", name="locations")
  
  ##default beta and nu coordinates to the over all coordinates.
  if( is.null(locations$coords.beta) ){
    locations$coords.beta <- locations$coords
  }
  if( is.null(locations$coords.nu) ){
    locations$coords.nu <- locations$coords
  }
  ##create an index vector of all the location specific things
  ##to extract from covars
  I <- c("ID", locations$coords, locations$coords.beta,
         locations$coords.nu, locations$long.lat, unique(locations$others))
  ##names of these
  names.loc <- c("ID", "x", "y", "x.beta", "y.beta", "x.nu", "y.nu")
  if( !is.null(locations$long.lat) ){
    names.loc <- c(names.loc, "long", "lat")
  }
  ##drop locations$others that match existing fields.
  locations$others <- unique(locations$others)
  if( any(locations$others %in% names.loc) ){
    I <- locations$others[locations$others %in% names.loc]
    warning( paste("Dropping locations$others already in names.loc:",
                   paste(locations$others[I], collapse=", ")) )
    locations$others <- locations$others[!I]
  }
  names.loc <- c(names.loc, locations$others)
  ##check for size consistency
  if( length(names.loc)!=length(I) ){
    stop("Length missmatch; most likely some of 'locations$coords' are not length=2.")
  }
  ##test that these actually exist in STdata
  stCheckFields(STdata$covars, I, name="STdata$covars")
  ##create location matrix.
  locations <- STdata$covars[,I,drop=FALSE]
  names(locations) <- names.loc

  ##return the location data.frame
  return( locations )
}##function processLocation

