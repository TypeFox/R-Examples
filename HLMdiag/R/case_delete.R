#' @export
case_delete <- function(model, ...){
  UseMethod("case_delete", model)
}

#' @export
#' @rdname case_delete.mer
#' @method case_delete default
#' @S3method case_delete default
case_delete.default <- function(model, ...){
  stop(paste("there is no case_delete() method for objects of class",
             paste(class(model), collapse=", ")))
}

#' Case Deletion for \code{mer}/\code{lmerMod} objects
#'
#'This function is used to iteratively delete groups corresponding to the
#'levels of a hierarchical linear model. It uses \code{lmer()} to fit
#'the models for each deleted case (i.e. uses brute force). To investigate
#'numerous levels of the model, the function will need to be called multiple
#'times, specifying the group (level) of interest each time.
#'
#' @export
#' @method case_delete mer
#' @S3method case_delete mer
#' @aliases case_delete
#'@param model the original hierarchical model fit using \code{lmer()}
#'@param group a variable used to define the group for which cases will be
#'deleted.  If this is left \code{NULL} (default), then the function will delete
#'individual observations.
#'@param type the part of the model for which you are obtaining deletion
#'diagnostics: the fixed effects (\code{"fixef"}), variance components
#'(\code{"varcomp"}), or \code{"both"} (default).
#'@param delete index of individual cases to be deleted.  For higher level
#'units specified in this manner, the \code{group} parameter must also be
#'specified.  If \code{delete = NULL} then all cases are iteratively deleted.
#' @param ... do not use
#'@return a list with the following compontents:
#' \describe{
#'   \item{\code{fixef.original}}{the original fixed effects estimates}
#'   \item{\code{ranef.original}}{the original predicted random effects}
#'   \item{\code{vcov.original}}{the original variance-covariance matrix for the fixed effects}
#'   \item{\code{varcomp.original}}{the original estimated variance components}
#'   \item{\code{fixef.delete}}{a list of the fixed effects estimated after case deletion}
#'   \item{\code{ranef.delete}}{a list of the random effects predicted after case deletion}
#'   \item{\code{vcov.delete}}{a list of the variance-covariance matrices for the fixed 
#'      effects obtained after case deletion}
#'   \item{\code{fitted.delete}}{a list of the fitted values obtained after case
#'      deletion}
#' \item{\code{varcomp.delete}}{a list of the estimated variance components obtained after
#'      case deletion}
#' }
#'@author Adam Loy \email{loyad01@@gmail.com}
#' @keywords models regression
#'@references Christensen, R., Pearson, L.M., and Johnson, W. (1992)
#'Case-Deletion Diagnostics for Mixed Models, \emph{Technometrics}, \bold{34}, 38
#'-- 45.
#'
#'Schabenberger, O. (2004) Mixed Model Influence Diagnostics, in
#'\emph{Proceedings of the Twenty-Ninth SAS Users Group International
#'Conference}, SAS Users Group International.
#'@examples
#'
#'data(sleepstudy, package = 'lme4')
#'fm <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
#'
#' # Deleting every Subject
#' fmDel <- case_delete(model = fm, group = "Subject", type = "both")
#'
#' # Deleting only subject 308
#' del308 <- case_delete(model = fm, group = "Subject", type = "both", delete = 308)
#' 
#' # Deleting a subset of subjects
#' delSubset <- case_delete(model = fm, group = "Subject", type = "both", delete = 308:310)
#'
case_delete.mer <- function(model, group = NULL, type = c("both", "fixef", "varcomp"), 
                        delete = NULL, ...){
  if(!is(model, "mer")) stop("model must be of class 'mer'")
  if(!model@dims["LMM"]){
    stop("case_delete is currently only implemented for mixed/hierarchical models.")
  }
  
  flist <- model@flist
  
  fixef.delete   <- NULL
  vcov.delete    <- NULL
  varcomp.delete <- NULL
  ranef.delete   <- NULL
  fitted.delete  <- NULL
  
  type <- match.arg(type) #default is "both"
  if( is.null(group) ){ # SINGLE CASE DELETION DIAGNOSTICS
    n <- model@dims[["n"]]
    modframe <- model@frame

    if( is.null(delete) ) {
      for(i in 1:n){
        model.delete <- lme4::lmer(formula = formula(model), data = model@frame[-i,])
        
        if(type %in% c("both", "varcomp")){
          if(length(lme4::getME(model.delete, "flist")) == 1) {
            ranef.delete[[i]] <- data.frame(deleted = i, 
                                            id = rownames(lme4::ranef(model.delete)[[1]]), 
                                            lme4::ranef(model.delete)[[1]])
          }
          else{
            ranef.delete[[i]] <- lme4::ranef(model.delete)
            ranef.delete[[i]] <- lapply(ranef.delete[[i]], function(x){
              x$id <- rownames(x)
              x$deleted <- i
              return(x)
            })
          }
          
          varcomp.delete[[i]] <- varcomp.mer(model.delete)
        }
        
        if(type %in% c("both", "fixef")){
          fixef.delete[[i]] <- c(deleted = i, lme4::fixef(model.delete))
          vcov.delete[[i]]  <- as.matrix(vcov(model.delete))
        }
        
        fitted.delete[[i]] <- data.frame(deleted = i, model.delete@frame, fitted(model.delete))
        
      }
    }
    else {
      model.delete   <- lme4::lmer(formula = formula(model), data = model@frame[-delete,])
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- lme4::fixef(model.delete)
        vcov.delete    <- as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.mer(model.delete)
        ranef.delete   <- lme4::ranef(model.delete)
        if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]]
      }
      fitted.delete  <- fitted(model.delete)
    }
  }

  else{ # MULTIPLE CASE DELETION DIAGNOSTICS
    
    if(!group %in% names(flist)) {
      stop(paste(group, "is not a valid grouping factor for this model."))
    }
    
    
    if( is.null(delete) ){
      data.delete <- split(model@frame, model@frame[, group])
      data.delete <- lapply(data.delete, function(df){
        index <- unique( df[, group ] )
        if(class(index) != "character") index <- as.character(index)
        data.delete[[ index ]] <- NULL
        do.call('rbind', data.delete)
      })
      
      model.delete <- lapply(data.delete, lme4::lmer, formula = formula(model))
      
      
      if(length(flist) == 1) {
        ranef.delete <- lapply(model.delete, function(x){
          data.frame(deleted = setdiff(model@frame[, group], x@frame[, group]),
                     id = rownames(lme4::ranef(x)[[1]]), lme4::ranef(x)[[1]])
        })
      }
      else{
        ranef.delete  <- lapply(model.delete, lme4::ranef)
        deleted.group <- rownames(lme4::ranef(model)[[group]])
        
        ranef.delete <- lapply(1:length(ranef.delete), function(x){
          ranef.list <- ranef.delete[[x]]
          lapply(ranef.list, function(y) {
            y$id <- rownames(y)
            y$deleted <- deleted.group[x]
            return(y)
          })
        })
      }
      
      varcomp.delete <- lapply(model.delete, varcomp.mer)
      
      if(type %in% c("both", "fixef")){
        fixef.delete <- lapply(model.delete, lme4::fixef)
        
        vcov.delete <- lapply(model.delete, vcov)
        vcov.delete <- lapply(vcov.delete, as.matrix)
      }
      
      
      fitted.delete <- lapply(model.delete, function(x){
        data.frame(deleted = setdiff(model@frame[, group], x@frame[, group]),
                   x@frame, fitted(x))
      })
    }
    else{
      index <- !model@frame[,group] %in% delete
      model.delete   <- lme4::lmer(formula = formula(model), data = model@frame[index,])
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- lme4::fixef(model.delete)
        vcov.delete    <-  as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.mer(model.delete)
        ranef.delete   <- lme4::ranef(model.delete)
        if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]]
      }
      fitted.delete  <- fitted(model.delete)
    }
    
    
  }
  
  # Organizing results
  if(is.null(delete)) {
    if(type %in% c("both", "fixef")){
      fitted.delete <- do.call('rbind', fitted.delete)
      #if(model@dims[["p"]] > 1) 
      fixef.delete  <- do.call('rbind', fixef.delete)
    }
    
    
    if(type %in% c("both", "varcomp")){
      if(length(lme4::getME(model, "flist")) == 1) {
        ranef.delete <- do.call('rbind', ranef.delete)
      }
      else {
        flist <- names(model@flist)
        temp  <- NULL
        for(i in 1:length(flist)) {
          temp[[i]] <- ldply(ranef.delete, function(x) x[[i]])
        }
        ranef.delete <- temp
        names(ranef.delete) <- names(lme4::ranef(model))
      }
    }
  }
  
  fixef.original <- model@fixef
  ranef.original <- lme4::ranef(model)
  if(length(ranef.original) == 1) ranef.original <- ranef.original[[1]]

  vcov.original <- as.matrix(vcov(model))
  varcomp.original <- varcomp.mer(model)

  val <- list(fixef.original = fixef.original, ranef.original = ranef.original,
              vcov.original = vcov.original, varcomp.original = varcomp.original,
              fixef.delete = fixef.delete, ranef.delete = ranef.delete,
              vcov.delete = vcov.delete, fitted.delete = fitted.delete,
              varcomp.delete = varcomp.delete)

  attr(val, "type") <- type
  class(val) <- "case_delete"
  return(val)
}



#' @export
#' @rdname case_delete.mer
#' @method case_delete lmerMod
#' @S3method case_delete lmerMod
case_delete.lmerMod <- function(model, group = NULL, type = c("both", "fixef", "varcomp"), 
                            delete = NULL, ...){
  if(!isNestedModel(model)){
    stop("case_delete is currently only implemented for mixed/hierarchical models.")
  }
  
  flist <- model@flist
  
  fixef.delete   <- NULL
  vcov.delete    <- NULL
  varcomp.delete <- NULL
  ranef.delete   <- NULL
  fitted.delete  <- NULL
  
  type <- match.arg(type) #default is "both"
  if( is.null(group) ){ # SINGLE CASE DELETION DIAGNOSTICS
    n <- lme4::getME(model, "n")
    modframe <- model@frame
    
    if( is.null(delete) ) {
      for(i in 1:n){
        model.delete <- lme4::lmer(formula = formula(model), data = model@frame[-i,])
        
        if(type %in% c("both", "varcomp")){
          if(length(lme4::getME(model.delete, "flist")) == 1) {
            ranef.delete[[i]] <- data.frame(deleted = i, 
                                            id = rownames(lme4::ranef(model.delete)[[1]]), 
                                            lme4::ranef(model.delete)[[1]])
          }
          else{
            ranef.delete[[i]] <- lme4::ranef(model.delete)
            ranef.delete[[i]] <- lapply(ranef.delete[[i]], function(x){
              x$id <- rownames(x)
              x$deleted <- i
              return(x)
            })
          }
          
          varcomp.delete[[i]] <- varcomp.mer(model.delete)
        }
        
        if(type %in% c("both", "fixef")){
          fixef.delete[[i]] <- c(deleted = i, lme4::fixef(model.delete))
          vcov.delete[[i]]  <- as.matrix(vcov(model.delete))
        }
        
        fitted.delete[[i]] <- data.frame(deleted = i, model.delete@frame, fitted(model.delete))
        
      }
    }
    else {
      model.delete   <- lme4::lmer(formula = formula(model), data = model@frame[-delete,])
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- lme4::fixef(model.delete)
        vcov.delete    <- as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.mer(model.delete)
        ranef.delete   <- lme4::ranef(model.delete)
        if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]]
      }
      fitted.delete  <- fitted(model.delete)
    }
  }
  
  else{ # MULTIPLE CASE DELETION DIAGNOSTICS
    
    if(!group %in% names(flist)) {
      stop(paste(group, "is not a valid grouping factor for this model."))
    }
    
    
    if( is.null(delete) ){
      data.delete <- split(model@frame, model@frame[, group])
      data.delete <- lapply(data.delete, function(df){
        index <- unique( df[, group ] )
        if(class(index) != "character") index <- as.character(index)
        data.delete[[ index ]] <- NULL
        do.call('rbind', data.delete)
      })
      
      model.delete <- lapply(data.delete, lme4::lmer, formula = formula(model))
      
      
      if(length(flist) == 1) {
        ranef.delete <- lapply(model.delete, function(x){
          data.frame(deleted = setdiff(model@frame[, group], x@frame[, group]),
                     id = rownames(lme4::ranef(x)[[1]]), lme4::ranef(x)[[1]])
        })
      }
      else{
        ranef.delete  <- lapply(model.delete, lme4::ranef)
        deleted.group <- rownames(lme4::ranef(model)[[group]])
        
        ranef.delete <- lapply(1:length(ranef.delete), function(x){
          ranef.list <- ranef.delete[[x]]
          lapply(ranef.list, function(y) {
            y$id <- rownames(y)
            y$deleted <- deleted.group[x]
            return(y)
          })
        })
      }
      
      varcomp.delete <- lapply(model.delete, varcomp.mer)
      
      if(type %in% c("both", "fixef")){
        fixef.delete <- lapply(model.delete, lme4::fixef)
        
        vcov.delete <- lapply(model.delete, vcov)
        vcov.delete <- lapply(vcov.delete, as.matrix)
      }
      
      
      fitted.delete <- lapply(model.delete, function(x){
        data.frame(deleted = setdiff(model@frame[, group], x@frame[, group]),
                   x@frame, fitted(x))
      })
    }
    else{
      index <- !model@frame[,group] %in% delete
      model.delete   <- lme4::lmer(formula = formula(model), data = model@frame[index,])
      
      if(type %in% c("both", "fixef")) {
        fixef.delete   <- lme4::fixef(model.delete)
        vcov.delete    <-  as.matrix(vcov(model.delete))
      }
      
      if(type %in% c("both", "varcomp")) {
        varcomp.delete <- varcomp.mer(model.delete)
        ranef.delete   <- lme4::ranef(model.delete)
        if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]]
      }
      fitted.delete  <- fitted(model.delete)
    }
    
    
  }
  
  # Organizing results
  if(is.null(delete)) {
    if(type %in% c("both", "fixef")){
      fitted.delete <- do.call('rbind', fitted.delete)
      #if(model@dims[["p"]] > 1) 
      fixef.delete  <- do.call('rbind', fixef.delete)
    }
    
    
    if(type %in% c("both", "varcomp")){
      if(length(lme4::getME(model, "flist")) == 1) {
        ranef.delete <- do.call('rbind', ranef.delete)
      }
      else {
        flist <- names(model@flist)
        temp  <- NULL
        for(i in 1:length(flist)) {
          temp[[i]] <- ldply(ranef.delete, function(x) x[[i]])
        }
        ranef.delete <- temp
        names(ranef.delete) <- names(lme4::ranef(model))
      }
    }
  }
  
  fixef.original <- lme4::fixef(model)
  ranef.original <- lme4::ranef(model)
  if(length(ranef.original) == 1) ranef.original <- ranef.original[[1]]
  
  vcov.original <- as.matrix(vcov(model))
  varcomp.original <- varcomp.mer(model)
  
  val <- list(fixef.original = fixef.original, ranef.original = ranef.original,
              vcov.original = vcov.original, varcomp.original = varcomp.original,
              fixef.delete = fixef.delete, ranef.delete = ranef.delete,
              vcov.delete = vcov.delete, fitted.delete = fitted.delete,
              varcomp.delete = varcomp.delete)
  
  attr(val, "type") <- type
  class(val) <- "case_delete"
  return(val)
}


# #' @export
# #' @rdname case_delete
# #' @method case_delete lme
# #' @S3method case_delete lme
# case_delete.lme <- function(model, group = NULL, type = c("both", "fixef", "varcomp"), 
                                # delete = NULL, ...){
  # if(!isNestedModel(model)){
    # stop("case_delete is currently only implemented for mixed/hierarchical models.")
  # }
  
  # flist <- model$groups
  
  # fixef.delete   <- NULL
  # vcov.delete    <- NULL
  # varcomp.delete <- NULL
  # ranef.delete   <- NULL
  # fitted.delete  <- NULL
  
  # type <- match.arg(type) #default is "both"
  # if( is.null(group) ){ # SINGLE CASE DELETION DIAGNOSTICS
    # n <- model$dims$N
    # modframe <- model$data
    # dataformula <- formula(modframe)
    # randcall <- model$call$random    
    
    # if( is.null(delete) ) {
      # for(i in 1:n){
        # if(is.null(randcall)) {
          # model.delete <- lme(formula(model), data = groupedData(dataformula, data = modframe[-i,]))
        # } else{
          # model.delete <- lme(formula(model), data = groupedData(dataformula, data = modframe[-i,]), random = randcall)
        # }
        
        # if(type %in% c("both", "varcomp")){
          # if(length(flist) == 1) {
            # ranef.delete[[i]] <- data.frame(deleted = i, 
                                            # id = rownames(ranef(model.delete)), 
                                            # ranef(model.delete))
          # }
          # else{
            # ranef.delete[[i]] <- ranef(model.delete)
            
            # if(is.list(ranef.delete[[i]])) { 
              # ranef.delete[[i]] <- lapply(ranef.delete[[i]], function(x){
                # x$id <- rownames(x)
                # x$deleted <- i
                # return(x)
              # })
            # } else{
              # ranef.delete[[i]]$id <- rownames(ranef.delete[[i]])
              # ranef.delete[[i]]$deleted <- i
            # }
            
          # }
          
          # varcomp.delete[[i]] <- varcomp(model.delete)
        # }
        
        # if(type %in% c("both", "fixef")){
          # fixef.delete[[i]] <- c(deleted = i, fixef(model.delete))
          # vcov.delete[[i]]  <- as.matrix(vcov(model.delete))
        # }
        
        # fitted.delete[[i]] <- data.frame(deleted = i, model.delete$data, fitted(model.delete))
        
      # }
    # }
    # else {
      # if(is.null(randcall)) {
        # model.delete <- lme(formula(model), data = modframe[-delete,])
      # } else{
        # model.delete <- lme(formula(model), data = modframe[-delete,], random = randcall)
      # }
      
      # if(type %in% c("both", "fixef")) {
        # fixef.delete   <- fixef(model.delete)
        # vcov.delete    <- as.matrix(vcov(model.delete))
      # }
      
      # if(type %in% c("both", "varcomp")) {
        # varcomp.delete <- varcomp(model.delete)
        # ranef.delete   <- ranef(model.delete)
        # if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]]
      # }
      # fitted.delete  <- fitted(model.delete)
    # }
  # }
  
  # else{ # MULTIPLE CASE DELETION DIAGNOSTICS
    # modframe <- model$data
    # dataformula <- formula(modframe)
    # modframe <- as.data.frame(modframe)
    # randcall <- model$call$random
    
    # if(!group %in% names(flist)) {
      # stop(paste(group, "is not a valid grouping factor for this model."))
    # }
    
    
    # if( is.null(delete) ){
      # data.delete <- split(modframe, modframe[, group])
      # data.delete <- lapply(data.delete, function(df, dataformula){
        # data.delete[[ unique( df[, group ] ) ]] <- NULL
        # temp <- do.call('rbind', data.delete)
        # groupedData(dataformula, temp)
      # }, dataformula = dataformula)
      
      # if(is.null(randcall)) {
        # model.delete <- lapply(data.delete, lme, fixed = formula(model))
      # } else{
        # model.delete <- lapply(data.delete, lme, fixed = formula(model), random = randcall)
      # }
      
      
      # if(length(flist) == 1) {
        # ranef.delete <- lapply(model.delete, function(x){
          # data.frame(deleted = setdiff(modframe[, group], x$data[, group]),
                     # id = rownames(ranef(x)), ranef(x))
        # })
      # }
      # else{
        # ranef.delete  <- lapply(model.delete, ranef)
        # deleted.group <- rownames(ranef(model))
        
        # ranef.delete <- lapply(1:length(ranef.delete), function(x){
          # ranef.list <- ranef.delete[[x]]
          # lapply(ranef.list, function(y) {
            # y$id <- rownames(y)
            # y$deleted <- deleted.group[x]
            # return(y)
          # })
        # })
      # }
      
      # varcomp.delete <- lapply(model.delete, varcomp)
      
      # if(type %in% c("both", "fixef")){
        # fixef.delete <- lapply(model.delete, fixef)
        
        # vcov.delete <- lapply(model.delete, vcov)
        # vcov.delete <- lapply(vcov.delete, as.matrix)
      # }
      
      
      # fitted.delete <- lapply(model.delete, function(x){
        # data.frame(deleted = setdiff(modframe[, group], x$data[, group]),
                   # x$data, fitted(x))
      # })
    # }
    # else{
      # index <- !modframe[,group] %in% delete
      # if(is.null(randcall)) {
        # model.delete   <- lme(formula(model), data = groupedData(dataformula, data = modframe[index,]))
      # } else{
        # model.delete   <- lme(formula(model), data = groupedData(dataformula, data = modframe[index,]), random = randcall)
      # }
      
      # if(type %in% c("both", "fixef")) {
        # fixef.delete   <- fixef(model.delete)
        # vcov.delete    <-  as.matrix(vcov(model.delete))
      # }
      
      # if(type %in% c("both", "varcomp")) {
        # varcomp.delete <- varcomp(model.delete)
        # ranef.delete   <- ranef(model.delete)
        # if( length(flist) == 1 ) ranef.delete <- ranef.delete[[1]]
      # }
      # fitted.delete  <- fitted(model.delete)
    # }
    
    
  # }
  
  # # Organizing results
  # if(is.null(delete)) {
    # if(type %in% c("both", "fixef")){
      # fitted.delete <- do.call('rbind', fitted.delete)
      # #if(model@dims[["p"]] > 1) 
      # fixef.delete  <- do.call('rbind', fixef.delete)
    # }
    
    
    # if(type %in% c("both", "varcomp")){
      # if(ncol(flist) == 1) {
        # ranef.delete <- do.call('rbind', ranef.delete)
      # }
      # else {
        # flist <- colnames(flist)
        # temp  <- NULL
        # for(i in 1:ncol(flist)) {
          # temp[[i]] <- ldply(ranef.delete, function(x) x[[i]])
        # }
        # ranef.delete <- temp
        # names(ranef.delete) <- names(ranef(model))
      # }
    # }
  # }
  
  # fixef.original <- fixef(model)
  # ranef.original <- ranef(model)
  # if(length(ranef.original) == 1) ranef.original <- ranef.original[[1]]
  
  # vcov.original <- as.matrix(vcov(model))
  # varcomp.original <- varcomp(model)
  
  # val <- list(fixef.original = fixef.original, ranef.original = ranef.original,
              # vcov.original = vcov.original, varcomp.original = varcomp.original,
              # fixef.delete = fixef.delete, ranef.delete = ranef.delete,
              # vcov.delete = vcov.delete, fitted.delete = fitted.delete,
              # varcomp.delete = varcomp.delete)
  
  # attr(val, "type") <- type
  # class(val) <- "case_delete"
  # return(val)
# }