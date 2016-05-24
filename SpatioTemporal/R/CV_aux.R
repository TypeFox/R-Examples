###################################
## Functions for crossvalidation ##
###################################
##Functions in this file:
## createCV           EX:ok
## dropObservations   EX:ok
## stCheckInternalCV  EX:INTERNAL
## predictNaive       EX:ok
## computeLTA         EX:ok

##' Creates a matrix that specifies cross-validation schemes.
##' 
##' The number of observations left out of each group can be rather uneven;
##' the main goal of \code{createCV} is to create CV-groups such that the groups
##' contain roughly the same \emph{number of locations} ignoring the number of
##' observations at each location. If there are large differences in the number
##' of observations at differnt locations one could use the \code{subset} option
##' to create different CV-groupings for different types of locations. If
##' \code{Icv.vector=FALSE}, the groups can then be combined as \cr
##' \code{I.final = I.1 | I.2 | I.3}.
##' 
##' The \code{option} input determines which sites to include in the
##' cross-validation. Possible options are \code{"all"}, \code{"fixed"},
##' \code{"comco"}, \code{"snapshot"} and \code{"home"}.
##' \describe{
##'   \item{all}{Uses all available sites, possibly subset according to
##'              \code{subset}. The sites will be grouped with sites
##'              seperated by less than \code{min.dist} being put in the
##'              same CV-group.}
##'   \item{fixed}{Uses only sites that have \cr
##'      \code{STmodel$locations$type \%in\% c("AQS","FIXED")}.
##'      Given the subsettting the sites will be grouped as for \code{all}.}
##'   \item{home}{Uses only sites that have \cr
##'      \code{STmodel$locations$type \%in\% c("HOME")}.
##'      Given the subsettting the sites will be grouped as for \code{all}.}
##'   \item{comco}{Uses only sites that have \cr
##'      \code{STmodel$locations$type \%in\% c("COMCO")}.
##'      The sites will be grouped together if they are from the same road
##'      gradient. The road gradients are grouped by studying the name of the
##'      sites. With "?" denoting one or more letters and "#" denoting one or
##'      more digits the names are expected to follow "?-?#?#", for random sites,
##'      and "?-?#?#?" for the gradients (with all but the last letter being the
##'      same for the entire gradient).}
##' }
##' 
##' @title Define Cross-Validation Groups
##'
##' @param STmodel Model object for which to determine cross-validation.
##' @param groups Number of cross-validation groups, zero gives leave-one-out
##'   cross-validation.
##' @param min.dist Minimum distance between locations for them to end up in
##'   separate groups. Points closer than \code{min.dist} will be forced into the
##'   same group. A high value for \code{min.dist} can result in fewer
##'   cross-validation groups than specified in \code{groups}.
##' @param random If \code{FALSE} repeated calls to the function will return the
##'   same grouping, if \code{TRUE} repeated calls will give different
##'   CV-groupings. Ensures that simulation studies are reproducable.
##' @param subset A subset of locations for which to define the cross-validation
##'   setup. Only sites listed in \code{subset} are dropped from one of the
##'   cross-validation groups; in other words sites \emph{not in} \code{subset}
##'   are used for estimation and preidiction of \emph{all} cross-validation
##'   groups. This option is \emph{ignored} if \code{option!="all"}.
##' @param option For internal MESA Air usage, see Details below.
##' @param Icv.vector Attempt to return a vector instead of a matrix. If the same
##'   observation is in several groups a matrix will still be returned.
##' 
##' @return Return a vector, with each element giving the CV-group (as an
##'   integer) of each observation; Or a (number or observations) - by -
##'   (groups) logical matrix; each column defines a cross-validation set with
##'   the \code{TRUE} values marking the observations to be left out.
##' 
##' @example Rd_examples/Ex_createCV.R
##' 
##' @author Johan Lindström
##' @family STmodel functions
##' @family cross-validation functions
##' @export
createCV <- function(STmodel, groups=10, min.dist=.1, random=FALSE,
                     subset=NA,
                     option=c("all","fixed","comco","snapshot","home"),
                     Icv.vector=TRUE){
  ##check class belonging
  stCheckClass(STmodel, "STmodel", name="STmodel")
  
  ##options - for internal MESA usage:
  ## all      (just divide over all sites)
  ## fixed    (only AQS and FIXED sites)
  ## comco    (only comco, keeping gradients together, ignores min.dist)
  ## snapshot (same as comco)
  ## home     (only home)
  option <- match.arg( option )
  
  if(option=="all"){
    Ind <- 1:dim(STmodel$locations)[1]
    if( all(!is.na(subset)) ){
      ##subset the data - first cast to character
      subset <- as.character(subset)
      ##check that all requested names exist in the colnames
      if( any(!(subset %in% STmodel$locations$ID)) )
        warning( sprintf("Subset names not found: %s",
                         paste(subset[!(subset %in% STmodel$locations$ID)],
                               collapse=", ")) )
      Ind <- which(STmodel$locations$ID %in% subset)
    }
  }else if(option=="fixed"){
    Ind <- which(toupper(STmodel$locations$type) %in% c("AQS","FIXED"))
  }else if(option=="home"){
    Ind <- which(toupper(STmodel$locations$type) %in% c("HOME"))
  }else{ ##remaining options are comco/snapshot
    Ind <- which(toupper(STmodel$locations$type) %in% c("COMCO"))
  }
  if(length(Ind)==0)
    stop( paste("No sites found for option:",option) )
  
  if(!random){
    ##set the seed to ensure that repeated calls returns the same grouping
    if( exists(".Random.seed") )
      old.seed <- .Random.seed ##save current seed so it can be re-set
    else
      old.seed <- NULL
    set.seed(0,kind="Mersenne-Twister")
  }

  if( option %in% c("all","fixed","home") ){
    ##extract distance matrix and find stations closer than min.dist
    D <- STmodel$D.beta[Ind,Ind]
    D[lower.tri(D,diag=TRUE)] <- NA
    Ind.coloc <- apply(D, 1, function(x){which(x<min.dist)} )
    Ind.G <- list()
    if( any(sapply(Ind.coloc,length)>0) ){
      ##find Ind.coloc's that contain colocated sites
      Ind.coloc.I <- unname(which(sapply(Ind.coloc,length)>0))
      i <- 1
      while( length(Ind.coloc.I)>0 ){
        ##extract the i:th colocated sites and add it to the grouping vector
        Ind.G[[i]] <- unique(c(Ind.coloc.I[1],Ind.coloc[[Ind.coloc.I[1]]]))
        ##drop the site from Ind.coloc.I
        Ind.coloc.I <- Ind.coloc.I[-1]
        ##see if this leads to additional matches
        ##find all Ind.coloc elemets that contain things from Ind.G[[i]]
        ##This is done both looking for matches
        Ind.tmp <- which(sapply(Ind.coloc,function(x)
                                any(x %in% Ind.G[[i]])))
        ##and adding the Ind.G[[i]] elements since the "diagonal"
        ##isn't found thtough matches
        Ind.tmp <- unique(c(Ind.tmp,Ind.G[[i]]))
        ##Finnaly let's only consider things not yet used.
        Ind.tmp <- Ind.tmp[Ind.tmp %in% Ind.coloc.I]
        ##Add and repeat until no more new elements match
        while( length(Ind.tmp)>0 ){
          Ind.G[[i]] <- unique(c(Ind.G[[i]], Ind.tmp,
                                 unlist(Ind.coloc[Ind.tmp])))
          Ind.coloc.I <- Ind.coloc.I[!(Ind.coloc.I %in% Ind.tmp)]
          Ind.tmp <- which(sapply(Ind.coloc,function(x)
                                  any(x %in% Ind.G[[i]])))
          Ind.tmp <- unique(c(Ind.tmp,Ind.G[[i]]))
          Ind.tmp <- Ind.tmp[Ind.tmp %in% Ind.coloc.I]
        }
        i <- i+1
      }
      for(i in 1:length(Ind.G))
        Ind.G[[i]] <- sort(Ind[Ind.G[[i]]])
    }
    ##extract points that are not colocated
    Ind.R <- Ind[!(Ind %in% unlist(Ind.G))]
  }else{ ##the snapshot sites
    ##drop starting non-numerical and then numerical characters
    ##from site ID, if the first character of whats left=="R" then
    ##this is a random site. -> seprate indecies into random and gradient
    ##random sites are on the format [A-Z]-[A-Z][0-9]R[0-9]
    ##so remove the leading characters and see if we find an R.
    Ind.R <- Ind[substr(sub("^[A-Z]+-[A-Z]+[0-9]+","",
                            STmodel$locations$ID[Ind]),1,1)=="R"]
    ##non random sites are the rest.
    Ind <- Ind[!(Ind %in% Ind.R)]
    ##now lets group the gradient sites
    grad.id <- unique(sub("[A-Z]+$","",
                          as.character(STmodel$locations$ID[Ind])))
    Ind.G <- list()
    for(i in 1:length(grad.id))
      Ind.G[[i]] <- Ind[grep(grad.id[i],STmodel$locations$ID[Ind])]
  }##if( option %in% c("all","fixed","home") ){...}else{...}

  ##we have no created an Ind.G and a Ind.R variable, lets permute and
  ##create the CV-indicator matrix
  if(groups==0){ ##leave one out CV, place each group seperately
    Ind.fin <- Ind.G
    if( length(Ind.R)!=0 )
      for(i in 1:length(Ind.R))
        Ind.fin[[i+length(Ind.G)]] <- Ind.R[i]
    groups <- length(Ind.fin)
  }else{    
    groups <- min(c(groups,length(Ind.R)+length(Ind.G)))
    ##randomly permute the sites
    Ind.R <- sample(Ind.R,length(Ind.R),replace=FALSE)
    if(length(Ind.G)!=0){
      Ind.G <- Ind.G[sample(length(Ind.G),length(Ind.G),replace=FALSE)]
    }
    ##collect sites to be in each group, first random sites
    Ind.R <- matrix(c(Ind.R,rep(NA,ceiling(length(Ind.R)/groups)*groups -
                                length(Ind.R))), nrow=groups)
    ##and then collect the grouped sites
    if(length(Ind.G)==0){
      Ind <- matrix(NA,groups,1)
    }else{
      Ind <- 1:length(Ind.G)
      Ind <- matrix(c(Ind,rep(NA,ceiling(length(Ind)/groups)*groups -
                              length(Ind))), nrow=groups)
      ##flips the matrix so that we have more groups at the bottom, since
      ##there are more random sites at the top.
      Ind <- Ind[seq(dim(Ind)[1],1,-1),,drop=FALSE]
    }
    ##compose a list, each element containing one cv-group
    Ind.fin <- list()
    for(i in 1:groups){
      Ind.fin[[i]] <- c(Ind.R[i,],unlist(Ind.G[Ind[i,]]))
      Ind.fin[[i]] <- Ind.fin[[i]][!is.na(Ind.fin[[i]])]
    }
  }##if(groups==0) else
  if(!random){
    ##re-set the random number generator
    if( !is.null(old.seed) )
      .Random.seed <- old.seed
  }
  if( any(sapply(Ind.fin,length)==0) )
    stop("Some CV-groups contain NO sites, try reducing the number of requested groups.")
  Ind.cv <- matrix(FALSE,dim(STmodel$obs)[1],groups)
  for(i in 1:groups)
    Ind.cv[,i] <- STmodel$obs$idx %in% Ind.fin[[i]]

  ##return just a vector?
  Ind.cv <- stCheckInternalCV(Ind.cv, Icv.vector)
  
  ##create a suitable cross-validation matrix, each colum
  ##marks the elements to drop for a specific CV.
  return(Ind.cv)
}##function createCV

##' Drops observations from \code{STmodel}, removing marked observations
##' along with the corresponding locations and recomputes a number of relevant
##' elements.
##' 
##' @title Drop Observations from a STmodel
##' 
##' @param STmodel Model object from which to drop observations.
##' @param Ind.cv A logical vector with one element per observation in \cr
##'   \code{STmodel$obs}. Observations marked with the \code{TRUE}
##'   will be dropped from the data structure. Use \code{\link{createCV}} to
##'   create the logical vector.
##' 
##' @return Returns the \code{STmodel} without the observations marked by
##'   \code{Ind.cv}. Only observed locations are retained.
##' 
##' @example Rd_examples/Ex_dropObservations.R
##' 
##' @author Johan Lindström
##' @family STmodel functions
##' @family cross-validation functions
##' @export
dropObservations <- function(STmodel, Ind.cv){
  ##check class belonging
  stCheckClass(STmodel, "STmodel", name="STmodel")
  
  ##check that the size of the indicator is consistent.
  N <- dim(STmodel$obs)[1]
  if(( is.vector(Ind.cv) && length(Ind.cv)!=N ) ||
     ( !is.vector(Ind.cv) && any(dim(Ind.cv) != c(N,1)) ) ){
    stop("Length of Ind.cv must match dim(STmodel$obs)[1]")
  }
  
  ##drop observations, temporal trends and spatio-temporal covars.
  STmodel$obs <- STmodel$obs[!Ind.cv,,drop=FALSE]
  STmodel$F <- STmodel$F[!Ind.cv,,drop=FALSE]
  if( length(STmodel$ST)!=0 ){
    STmodel$ST <- STmodel$ST[!Ind.cv,,drop=FALSE]
  }
  
  ##find which locations that remain
  IND <- STmodel$locations$ID %in% unique(STmodel$obs$ID)
  ##drop locations
  STmodel$locations <- STmodel$locations[IND,,drop=FALSE]
  ##and spatio-temporal covariates
  if( length(STmodel$ST.all)!=0 ){
    STmodel$ST.all <- STmodel$ST.all[,IND,,drop=FALSE]
  }
  ##find missing parts in LUR (repeat due to possible missing locations)
  for(i in 1:length(STmodel$LUR)){
    IND <- rownames(STmodel$LUR[[i]]) %in% unique(STmodel$obs$ID)
    STmodel$LUR[[i]] <- STmodel$LUR[[i]][IND,,drop=FALSE]
  }
  STmodel$LUR.all <- STmodel$LUR
  ##and recompute the indecies
  STmodel$obs$idx <- match(STmodel$obs$ID, STmodel$locations$ID)
  ##...and some other things
  STmodel <- createSTmodelInternalDistance(STmodel)
  ##...and the covariance
  STmodel <- updateCovf(STmodel)

  ##return stripped structure
  return(STmodel)
}##function dropObservations

stCheckInternalCV <- function(Ind.cv, Icv.vector=TRUE){
  if( is.null(Ind.cv) ){
    stop("Ind.cv must be defined (is NULL)")
  }else if( is.vector(Ind.cv) ){
    ##ensure integers
    Ind.cv <- as.integer(Ind.cv)
    uInd.cv <- sort(unique(Ind.cv))
    ##find missing/empty groups
    I.rep <- match(1:max(uInd.cv),uInd.cv)
    I.drop <- which(is.na(I.rep))
    if( length(I.drop!=0) ){
      ##renumber CV groups
      I.rep[ !is.na(I.rep) ] <- 1:sum(!is.na(I.rep))
      I.rep <- c(I.rep,0)
      Ind.cv[Ind.cv==0] <- length(I.rep)
      Ind.cv <- I.rep[Ind.cv]
    }
  }else{
    I.drop <- which( colSums(Ind.cv)==0 )
    ##renumber CV groups
    if( length(I.drop!=0) ){
      Ind.cv <- Ind.cv[, -I.drop, drop=FALSE]
    }
    if( Icv.vector && max(rowSums(Ind.cv))==1 ){
      Ind.cv <- apply(Ind.cv, 1, function(x){ x=which(x); if(length(x)!=1){x=0};
                                              return(x)})
    }
  }##if( is.vector(Ind.cv) ){...}else{...}
  if( length(I.drop!=0) ){
    warning( paste("Empty cross validation group(s) dropped:",
                   paste(I.drop, collapse=", ")) )
  }
  return(Ind.cv)
}##function stCheckInternalCV

##' Computes naive predictions that are based on a few sites. These predictions
##' can then be used, e.g. in \code{\link{summary.predCVSTmodel}}, to evaluate
##' how much better the spatial-temporal model performs compared to simple
##' (temporal) predictions.
##' \cr
##' The function requires one of \code{location} and \code{type} to be
##' specified, if both are given \code{location} \emph{will be used over}
##' \code{type}. If \code{type} is given locations such that \cr
##' \code{as.character(STmodel$locations$type) %in% type} \cr will be
##' used.
##' 
##' Given a set of locations the function computes 4 sets of naive prediction
##' for the observations in \code{STmodel}:
##' \describe{
##'   \item{smooth.fixed}{The smooth trend in \code{STmodel$trend} is fit to
##'                       \emph{all} observations at the sites in
##'                       \code{locations} using a linear regression. The
##'                       resulting smooth is used as a naive prediction for
##'                       all locations.}
##'   \item{avg.fixed}{The temporal average over sites in \code{locations} is
##'                    used as a naive prediction.}
##'   \item{smooth.closest.fixed}{This fits the smooth trend in
##'                               \code{STmodel$trend} to each site in
##'                               \code{locations}; using the smooth at the
##'                               closest fixed site as a naive prediction.}
##'   \item{closest.fixed}{This uses the observations at the closest site in
##'                        \code{locations} to predict observations at all other
##'                        sites.}
##' }
##' 
##' @title Naive Temporal Predictions
##' 
##' @param STmodel \code{STmodel} object for which to compute simple
##'   predictions. 
##' @param locations Locations on which to base the naive predictions.
##' @param type The type of sites to base the predictions on, uses the
##'   (optional) field \code{STmodel$locations$type}.
##' @return A list with items:
##'   \item{pred}{A (number of observations) - by - (6) data.frame containing
##'               the four naive predictions described under \code{details},
##'               along with dates and IDs.}
##'   \item{locations}{The locations used for the naive predictions.}
##' 
##' @example Rd_examples/Ex_predictNaive.R
##' 
##' @author Johan Lindström
##' @family STmodel functions
##' @family cross-validation functions
##' @export
predictNaive <- function(STmodel, locations=NULL, type=NULL){
  ##check class belonging
  stCheckClass(STmodel, "STmodel", name="STmodel")

  ##matrix holding the predictions
  pred <- data.frame(smooth.fixed=double( length(STmodel$obs$obs) ),
                     avg.fixed=double( length(STmodel$obs$obs) ),
                     smooth.closest.fixed=double( length(STmodel$obs$obs) ),
                     closest.fixed=double( length(STmodel$obs$obs) ))
  pred$ID <- STmodel$obs$ID
  pred$date <- STmodel$obs$date
  ##select sites to base the predictions on
  if( !is.null(locations) ){
    if( is.character(locations) ){
      IND.loc <- locations
    }else{
      IND.loc <- STmodel$locations$ID[as.double(locations)]
    }
  }else if( !is.null(type) ){
    IND.loc <- STmodel$locations$ID[which(STmodel$locations$type %in% type)]
  }else{
    stop("Both locations and type are NULL, you must specify one.")
  }
  ##check for specified locations
  I <- IND.loc %in% unique(STmodel$obs$ID)
  if( any(!I) ){
    warning( paste("unobserved locations", paste(IND.loc[!I], collapse=", ")) )
    IND.loc <- IND.loc[I]
  }
  if( length(IND.loc)==0 ){
    stop("No sites to use for naive prediction found")
  }
  ##indicator picking out the sites to use for naive prediction.
  IND <- (STmodel$obs$ID %in% IND.loc)
  ##indicator for which is the closest fixed site
  IND.closest <- apply(STmodel$D.beta[IND.loc,,drop=FALSE],2,order)
  ##special case for only one site
  if( length(IND.loc)==1 || !is.matrix(IND.closest)){
    IND.closest <- matrix(IND.closest, nrow=1)
    colnames(IND.closest) <- colnames(STmodel$D.beta)
  }
  
  ##1) fit temporal smooths to all the fixed sites
  pred[,1] <- c(STmodel$F %*% lm(STmodel$obs$obs[IND] ~
                                 STmodel$F[IND,]-1)$coefficients)
  
  ##2) use average of the fixed sites to predict
  date.tmp <- as.double(STmodel$obs$date)
  for(t in unique(date.tmp)){
    I <- date.tmp==t
    pred[I,2] <- mean(STmodel$obs$obs[IND & I])
  }
  
  ##3) fit temporal smooth to the closest fixed site to predict the data
  ##start by fitting the temporal to each of the fixed sites
  pred.tmp <- matrix(NA,length(STmodel$obs$obs),length(IND.loc))
  for(i in 1:length(IND.loc)){
    IND.tmp <- STmodel$obs$ID==IND.loc[i]
    pred.tmp[,i] <- STmodel$F %*%
      lm(STmodel$obs$obs[IND.tmp] ~
         STmodel$F[IND.tmp,]-1)$coefficients
  }
  for(i in unique(STmodel$obs$ID)){
    pred[pred$ID==i,3] <- pred.tmp[STmodel$obs$ID==i,IND.closest[1,i]]
  }

  ##4) use observations from the closest available site
  pred.tmp <- matrix(NA,length(STmodel$obs$obs),length(IND.loc))
  for(i in 1:length(IND.loc)){
    IND.tmp <- STmodel$obs$ID==IND.loc[i]
    obs.tmp <- STmodel$obs$obs[IND.tmp]
    pred.tmp[,i] <- obs.tmp[match(STmodel$obs$date,
                                  STmodel$obs$date[IND.tmp])]
  }
  for(i in unique(STmodel$obs$ID)){
    pred[pred$ID==i,4] <-
      apply(pred.tmp[pred$ID==i, IND.closest[,i], drop=FALSE], 1,
            function(x) switch(all(is.na(x))+1,x[min(which(!is.na(x)))],NA))
  }
  
  return( list(pred=pred, locations=IND.loc) )
}##predictNaive

##' Computes the long term average of observations and cross-validated
##' predictions for each of the sites in \code{object}. The long term averages
##' are computed using \emph{only} timepoints that have observations, this
##' applies to both the observed and predicted. Also the function allows for a
##' transformation: if requested the transformation is applied \emph{before} the
##' averaging.
##' 
##' @title Computes the Long Term Average for Each Sites.
##' 
##' @param object A \code{predCVSTmodel} object, the result of
##'   \code{\link{predictCV.STmodel}}. 
##' @param transform Transform observations (\emph{without} bias correction) and
##'   predictions \emph{before} 
##'   computing averages; e.g. \code{transform=exp} gives the long term averages
##'   as  \code{mean( exp(obs) )} and \code{mean( exp(pred) )}.
##' @return Returns a (number of locations) - by - 4 matrix with the observed
##'   and predicted value (using the three different model parts) for each
##'   location. 
##'
##' @example Rd_examples/Ex_computeLTA.R
##' 
##' @author Johan Lindström
##' @family predCVSTmodel functions
##' @family cross-validation functions
##' @export
computeLTA <- function(object, transform=function(x){return(x)}){
  ##check class belonging
  stCheckClass(object, "predCVSTmodel", name="object")

  ##is data already transformed - LTA is a bad idea
  if( !is.null(object$opts$transform) && !missing(transform)){
    stop("Data already transformed by predictCV")
  }
  
  ##check transformation
  if( !is.function(transform) ){
    stop("'transform' should be a function")
  }

  ##apply transformation
  I <- names(object$pred.obs) %in% c("obs","EX.mu","EX.mu.beta","EX","EX.pred")
  out <- transform( object$pred.obs[,I,drop=FALSE] )
  ##split by site
  out <- split(out, object$pred.obs$ID)
  ##drop NA:s
  for(i in 1:length(out)){
    I <- !apply( is.na(out[[i]]), 1, any)
    out[[i]] <- out[[i]][I,,drop=FALSE]
  }
  ##compute averages for each site, and convert to data.frame
  out <- as.data.frame( t(sapply(out, colMeans)) )
  
  return(out)
}##function computeLTA
