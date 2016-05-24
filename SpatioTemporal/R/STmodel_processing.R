###############################################
## FUNCTIONS THAT MANIPULATE STmodel OBJECTS ##
###############################################
##Functions in this file:
## updateCovf    EX:ok
## createLUR     EX:INTERNAL
## createST      EX:INTERNAL

#####################################################################
## Helper function that updates covariance function specifications ##
#####################################################################

##' Updates/sets the covariance functions for \code{STmodel} objects. Used by
##' \code{\link{createSTmodel}}.
##'
##' The covariance function is specified using lists for \code{cov.beta} and
##' \code{cov.nu}. The lists should contain the following elements:
##' \describe{
##'   \item{covf}{The type of covariance function(s), see
##'     \code{\link{namesCovFuns}}.}
##'   \item{nugget}{For the beta-fields: a vector of \code{TRUE}/\code{FALSE}
##'     indicating if each beta-field should contain a nugget. \cr
##'     For the nu-field: Either \code{TRUE}/\code{FALSE} for constant nugget/no
##'     nugget; a formula; or length=1 character vector. For the latter two the
##'     nugget is allowed to vary as \code{exp(B*theta)} where: \cr
##'  \code{nugget = as.formula(paste("~", paste(cov.nu$nugget, collapse="+")))} \cr
##'  \code{ covars = model.frame(nugget, covars, drop.unused.levels=TRUE)} \cr
##'  \code{ B=model.matrix(nugget, covars)} \cr
##'  \code{ B=as.matrix(B) } \cr
##'     The resulting regression matrix is stored as \code{STmodel$cov.nu$nugget.matrix}
##'     giving nugget for the observed locations. Unobserved locations are assumed
##'     to have a zero nugget.}
##'   \item{random.effect}{Only used for \code{cov.nu}, \code{TRUE}/\code{FALSE}
##'     indicating if a random.effect for the mean value should be included, see
##'     \code{\link{makeSigmaNu}}.}
##' }
##' 
##' @title Update Covariance Functions in \code{STmodel} Objects
##' @param STmodel \code{STmodel} object with observations, covariates, trends, etc;
##'   see \code{\link{mesa.model}}.
##' @param cov.beta,cov.nu Covariance specification for the beta- and nu-fields
##'   should contain fields \code{covf}, \code{nugget}, and \code{random.effect}
##'   (for the nu field); see details for description of these fields.
##'   For \code{cov.beta} the fields should contain one element for each
##'   smooth-temporal trend / beta-field if the fields have only one element,
##'   these elements are repeated implying the same covariance for all beta-fields.
##' @return updated version of \code{STmodel} with new covariance specifications.
##' 
##' @example Rd_examples/Ex_updateCovf.R
##' 
##' @author Johan Lindström
##' @family STmodel functions
##' @family covariance functions
##' @export
updateCovf <- function(STmodel, cov.beta=STmodel$cov.beta,
                       cov.nu=STmodel$cov.nu){
  ##check class belonging
  stCheckClass(STmodel, "STmodel", name="STmodel")
  ##check for fields
  stCheckFields(cov.beta, c("covf","nugget"), name="cov.beta")
  stCheckFields(cov.nu, c("covf","nugget","random.effect"), name="cov.nu")

  ##find the number of temporal trends (including intercept)
  nt <- dim(STmodel$trend)[2]
  
  ##expand elements in cov.beta
  for(i in 1:length(cov.beta) ){
    if( length(cov.beta[[i]])==1 ){
      cov.beta[[i]] <- rep(cov.beta[[i]], nt)
    }else if( length(cov.beta[[i]])!=nt ){
      stop( sprintf("Covariance specification, %s, for beta has %d!=%d elements.",
                    names(cov.beta)[i], length(cov.beta[[i]]), nt) )
    }
  }##for(i in 1:length(cov.beta) )
  ##check that covariance functions are valid
  tmp <- unlist(lapply( parsCovFuns( c(cov.beta$covf, cov.nu$covf) ), is.null))
  if( any(tmp) ){
    stop( paste("Unknown covariance specification(s):",
                paste(names(tmp)[tmp], collapse=", ")) )
  }
  ##check random effect
  if( !is.logical(cov.nu$random.effect) ){
    stop("'random.effect' specification must be 'logical'.")
  }
  ##check nugget - beta
  if( !is.logical(cov.beta$nugget) ){
    stop("'nugget' specification for beta-fields must be 'logical'.")
  }
  ##check iid and nugget - beta
  if( any(cov.beta$covf=="iid" & !cov.beta$nugget) ){
    stop("beta-fields: Covariance model iid requires a nugget.")
  }
  ##check nugget - nu
  ##... but first drop unobserved locations from the covars
  covars <- STmodel$locations
  Ind <- STmodel$locations$ID %in% unique(STmodel$obs$ID)
  covars <- covars[Ind,,drop=FALSE]
  
  if( is.logical(cov.nu$nugget) ){
    ##either constant or nothing.
    if( cov.nu$nugget[1]==TRUE ){
      cov.nu$nugget <- as.formula("~1", env=.GlobalEnv)
    }else{
      cov.nu$nugget <- FALSE
      cov.nu$nugget.matrix <- matrix(NULL, dim(covars)[1], 0)
    }
  }
  if( class(cov.nu$nugget)=="formula" || is.character(cov.nu$nugget) ){
    if( is.character(cov.nu$nugget) ){
      cov.nu$nugget <- as.formula(paste("~", paste(cov.nu$nugget, collapse="+")),
                                  env=.GlobalEnv)
    }
    covars.tmp <- model.frame(cov.nu$nugget, covars, drop.unused.levels=TRUE)
    ##A formula has been specified, use this
    tmp <- tryCatch(model.matrix(cov.nu$nugget, covars.tmp),
                    error = function(e){
                      tmp <- model.matrix(cov.nu$nugget, covars);
                      ##find constant columns
                      test <- apply(tmp, 2, function(x){length(unique(x))})
                      test <- test!=1 | names(test)=="(Intercept)"
                      warning(
  paste("Some covariates are constant or factors with only one level",
        "in the specification of the nugget. This will result in",
        "unidentifiable parameters. To avoid this the following elemets",
        "of the formula have been dropped:",
        paste(names(test[!test]),collapse=", "), sep="\n") )
                      ##drop constant columns
                      tmp <- tmp[,test, drop=FALSE];
                      return(tmp)})
    ##store nugget matrix.
    cov.nu$nugget.matrix <- as.matrix( tmp )
  }else if( !is.logical(cov.nu$nugget) ){
    stop( paste("Unknown specification of nugget for nu-field:", cov.nu$nugget) )
  }
  ##check iid and nugget - beta
  if( dim(cov.nu$nugget.matrix)[2]==0 && cov.nu$covf=="iid" ){
    stop("nu-field: Covariance model iid requires a nugget.")
  }
  ##add rownames
  rownames(cov.nu$nugget.matrix) <- covars$ID

  ##return
  STmodel$cov.beta <- cov.beta
  STmodel$cov.nu <- cov.nu
  return( STmodel )
}##function updateCovf


#####################################################
## Helper functions that add LUR and ST covariates ##
#####################################################

##' Extracts the requested geographic and spatio-temporal covariates from
##' a \code{STmodel} object and formats them into suitable matrices.
##' For \strong{INTERNAL} use by \code{\link{createSTmodel}}
##' 
##' @title Add Covariate Fields to \code{STdata} Object.
##' @param STdata \code{STdata} object with observations, covariates, trends,
##'   etc; see \code{\link{mesa.data.raw}}.
##' @param LUR.list specification of covariates; e.g. output from
##'   \code{\link{processLUR}}.
##' 
##' @return \code{STdata} with added fields \code{LUR}, or \code{ST} and
##'   \code{ST.all}.
##' 
##' @author Johan Lindström
createLUR <- function(STdata, LUR.list){
  ##check class belonging
  stCheckClass(STdata, "STdata", name="STdata")

  ##extract important fields
  covars <- STdata$covars
  scale.covars <- STdata$scale.covars

  ##should we scale the covariates
  if( !is.null(scale.covars) ){
    ##check that we have all the covariates
    stCheckFields(scale.covars$mean, names(covars), name="scale.covars$mean")
    stCheckFields(scale.covars$sd, names(covars), name="scale.covars$sd")

    ##ensure right order
    scale.covars$mean <- scale.covars$mean[ names(covars) ]
    scale.covars$sd <- scale.covars$sd[ names(covars) ]
    
    ##scale all non NA:s
    I <- !is.na(scale.covars$mean)
    n1 <- dim(covars)[1]
    n2 <- sum(I)
    covars[,I] <- scale(covars[,I,drop=FALSE], center=scale.covars$mean[I],
                        scale=scale.covars$sd[I])
  }
  ##create a list of LUR:s, one for each temporal trend
  STdata$LUR.all <- vector("list", length(LUR.list))
  names(STdata$LUR.all) <- names(LUR.list)
  for( i in 1:length(LUR.list) ){
    STdata$LUR.all[[i]] <- model.matrix(LUR.list[[i]], covars)
    rownames(STdata$LUR.all[[i]]) <- covars$ID
  }
  ##also create a reduced LUR for observations only
  if( length(STdata$obs$obs)!=0 ){
    STdata$LUR <- STdata$LUR.all
    I.idx <- 1:max(match(STdata$obs$ID,STdata$covars$ID))
    for( i in 1:length(STdata$LUR) ){
      STdata$LUR[[i]] <- STdata$LUR[[i]][I.idx,,drop=FALSE]
    }
  }
  
  return(STdata)
}##function createLUR


##' @param ST.list specification of spatio-temporal covariates; e.g. output from
##'   \code{\link{processST}}.
##' @rdname createLUR
createST <- function(STdata, ST.list){
  ##check class belonging
  stCheckClass(STdata, "STdata", name="STdata")

  ##extract important fields
  SpatioTemporal <- STdata$SpatioTemporal
  obs <- STdata$obs

  if( length(ST.list)!=0 ){
    ##ST-covariate for all locations
    STdata$ST.all <- SpatioTemporal[,,ST.list,drop=FALSE]

    ##and now for the observations
    STdata$ST <- matrix(0, dim(obs)[1], length(ST.list))

    ##dates and ID:s for the spatio-temporal covaraites
    ST.dates <- convertCharToDate( rownames(SpatioTemporal) )
    ST.ID <- colnames(SpatioTemporal)
    ##correspinding index
    I <- (match(obs$ID,ST.ID)-1)*length(ST.dates) + match(obs$date,ST.dates)
    ##extract these points from the Spatio temporal covariate
    for(i in 1:length(ST.list)){
      tmp <- SpatioTemporal[,,ST.list[i],drop=FALSE]
      STdata$ST[,i] <- tmp[I]
    }
    colnames(STdata$ST) <- ST.list
    STdata$ST <- data.matrix(STdata$ST)
  }
  return( STdata )
}##function createST
