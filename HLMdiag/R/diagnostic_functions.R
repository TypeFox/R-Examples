#' Calculating influence diagnostics for HLMs.
#'
#' @description
#' This group of functions is used to compute deletion diagnostics for a
#' hierarchical linear model based on the building blocks returned by
#' \code{case_delete}. 
#'
#' @details
#' The primary function is \code{diagnostics} which returns either a
#' list or data frame of influence measures depending on whether
#' \code{type = "both"} (\code{list}) or if only one aspect of the model 
#' is selected (\code{data.frame}).
#' If \code{type = "both"}, then a list with Cook's distance, MDFFITS,
#' COVTRACE, and COVRATIO are returned for the fixed effects and
#' relative variance change (RVC) is returned for the variance components.
#'
#' The methods \code{cooks.distance}, \code{mdffits}, \code{covtrace},
#' \code{covratio}, and \code{rvc} can be used for direct computation
#' of the corresponding diagnostic quantities from an object of class
#' \code{case_delete}.
#' 
#' @note
#' The results provided by this function will give exact values of the 
#' diagnostics; however, these are computationally very slow. Approximate
#' versions of \code{cooks.distance}, \code{mdffits}, \code{covtrace},
#' \code{covratio} are implemented in HLMdiag, and can be called directly on
#' the \code{mer} object.
#'
#' @aliases cooks.distance.case_delete mdffits.case_delete covtrace.case_delete covratio.case_delete rvc.case_delete
#' @param object an object containing the output returned by \code{case_delete()}
#' @author Adam Loy \email{loyad01@@gmail.com}
#' @keywords models regression
#' @references Christensen, R., Pearson, L.M., and Johnson, W. (1992)
#' ``Case-Deletion Diagnostics for Mixed Models, \emph{Technometrics},
#' \bold{34}, 38 -- 45.
#'
#' Dillane, D. (2005), Deletion Diagnostics for the Linear Mixed Model,''
#' Ph.D. thesis, Trinity College Dublin.
#'
#' Schabenberger, O. (2004) Mixed Model Influence Diagnostics,
#' in \emph{Proceedings of the Twenty-Ninth SAS Users Group International Conference},
#' SAS Users Group International.
#'
#' @export
#' @seealso \code{\link{leverage.mer}}, 
#' \code{\link{cooks.distance.mer}}, \code{\link{mdffits.mer}},
#' \code{\link{covratio.mer}}, \code{\link{covtrace.mer}}
#' @examples
#' data(sleepstudy, package = 'lme4')
#' fm <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' 
#' # Subject level deletion and diagnostics
#' subject.del  <- case_delete(model = fm, group = "Subject", type = "both")
#' subject.diag <- diagnostics(subject.del)
diagnostics <- function(object){
  type <- attributes(object)$type
  if(type %in% c("fixef", "both")){
    ids <- as.vector(rownames(object$fixef.delete, do.NULL = FALSE, prefix = ""))
  }
  else{
    ids <- as.vector(names(object$varcomp.delete))
    if(is.null(ids)) ids <- 1:length(object$varcomp.delete)
  }
  if(type  %in% c("fixef", "both")){
    if(!is(object$fixef.delete, "matrix")) {
      res1 <- data.frame(COOKSD = cooks.distance(object),
                         MDFFITS = mdffits(object),
                         COVTRACE = covtrace(object),
                         COVRATIO = covratio(object))
    }
    else {
      res1 <- data.frame(IDS = ids, COOKSD = cooks.distance(object),
                         MDFFITS = mdffits(object),
                         COVTRACE = covtrace(object),
                         COVRATIO = covratio(object))
    }
    if(type == "fixef") return(res1)
  }
  
  if(type %in% c("varcomp", "both")){
    if(!is(object$varcomp.delete, "list")) { res2 <- data.frame(rvc(object)) }
    else res2 <- data.frame(IDS = ids, rvc(object))
    
    if(type == "varcomp") return(res2)
  }
  
  if(type == "both"){
    res <- list(fixef_diag = res1, varcomp_diag = res2)
    return(res)
  } 
}

#' @export
#' @rdname diagnostics
#' @param model an object containing the output returned by \code{case_delete()}.
#' This is only named differently to agree with the generic.
#' @S3method cooks.distance case_delete
#' @method cooks.distance case_delete
cooks.distance.case_delete <- function(model, ...){
  p <- length(model$fixef.original)

  if(is(model$fixef.delete, "matrix")) {
    groups <- rownames(model$fixef.delete, do.NULL = FALSE, prefix = "")
    cook <- NULL
    for(i in 1:length(groups)){
      change.fixef <- as.matrix(model$fixef.original - model$fixef.delete[i,])
      cook <- c(cook, t(change.fixef) %*% solve( as.matrix( model$vcov.original ) ) %*% change.fixef / p)
    }
  }
  else{
    change.fixef <- as.matrix(model$fixef.original - model$fixef.delete)
    cook <- t(change.fixef) %*% solve( as.matrix( model$vcov.original ) ) %*% change.fixef / p
  }

  return(cook)
}


#' @export
#' @rdname diagnostics
#' @S3method mdffits case_delete
#' @method mdffits case_delete
mdffits.case_delete <- function(object, ...){
  p <- length(object$fixef.original)

  if(is(object$fixef.delete, "matrix")) {
    groups <- rownames(object$fixef.delete, do.NULL = FALSE, prefix = "")
    MDFFITS <- NULL
    for(i in 1:length(groups)){
      change.fixef <- as.matrix(object$fixef.original - object$fixef.delete[i,])
      MDFFITS <- c(MDFFITS, t(change.fixef) %*% solve( as.matrix(object$vcov.delete[[i]]) ) %*% change.fixef / p)
    }
  }
  else{
    change.fixef <- as.matrix(object$fixef.original - object$fixef.delete)
    MDFFITS <- t(change.fixef) %*% solve( as.matrix(object$vcov.delete) ) %*% change.fixef / p
  }

  return(MDFFITS)
}


#' @export
#' @rdname diagnostics
#' @S3method covtrace case_delete
#' @method covtrace case_delete
covtrace.case_delete <- function(object, ...){
  p <- length(object$fixef.original)

  if(is(object$vcov.delete, "list")) {
    groups <- rownames(object$fixef.delete, do.NULL = FALSE, prefix = "")
    COVTRACE <- NULL
    for(i in 1:length(groups)){
      V.original <- as.matrix(object$vcov.original)
      V.delete <- as.matrix(object$vcov.delete[[i]])
      COVTRACE <- c(COVTRACE, abs(sum(diag(solve(V.original) %*% V.delete)) - p))
    }
  }
  else{
    V.original <- as.matrix(object$vcov.original)
    V.delete <- as.matrix(object$vcov.delete)
    COVTRACE <- abs(sum(diag(solve(V.original) %*% V.delete)) - p)
  }

  return(COVTRACE)	
}



#' @export
#' @rdname diagnostics
#' @S3method covratio case_delete
#' @method covratio case_delete
covratio.case_delete <- function(object, ...){
  if(is(object$vcov.delete, "list")) {
    groups <- rownames(object$fixef.delete, do.NULL = FALSE, prefix = "")
    COVRATIO <- NULL
    for(i in 1:length(groups)){
      V.original <- as.matrix(object$vcov.original)
      V.delete <- as.matrix(object$vcov.delete[[i]])
      COVRATIO <- c(COVRATIO, det(V.delete) / det(V.original))
    }
  }
  else{
    V.original <- as.matrix(object$vcov.original)
    V.delete <- as.matrix(object$vcov.delete)
    COVRATIO <- det(V.delete) / det(V.original)
  }
  
  return(COVRATIO)
}



#' @export
#' @rdname diagnostics
#' @param ... do not use
#' @S3method rvc case_delete
#' @method rvc case_delete
rvc.case_delete <- function(object, ...){
	if(class(object$varcomp.delete) == "list") {
	  res <- do.call('rbind', lapply(object$varcomp.delete, function(x){ (x / object$varcomp.original) - 1}))
	}
  else{
    res <- (object$varcomp.delete / object$varcomp.original) - 1
  }
	return(res)
}
