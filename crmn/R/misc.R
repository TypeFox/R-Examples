##' Show method for nFit
##'
##' @title Show nfit
##' @param object the \code{nFit} object
##' @return prints some basic information
##' @aliases show_nfit
##' @author Henning Redestig
show_nfit <- function(object) {
  cat(method(object), "normalization model\n")
  cat("========================\n")
  if(method(object) %in% c("crmn")) {
    cat("Effect of experiment design on standards:\n")
    cat("-----------------------------------------\n")
    print(anova(sFit(object)$fit$fit))
    cat("\nCaptured Tz:\n")
    cat("--------------\n")
    print(summary(sFit(object)$fit$pc))
    cat("\n")
  }
  if(method(object) %in% c("crmn", "nomis")) {
    cat("R2 from Tz to analytes:\n")
    cat("-----------------------\n")
    print(1 - sum(residuals(model(object)$fit)^2) /
          (sum(fitted(model(object)$fit)^2) + sum(residuals(model(object)$fit)^2)))
  }
}

##' Simple plot function for a CRMN normalization model.
##' 
##' Shows Tz and the optimization (if computed) of the PCA model. The
##' number of components used for normalization should not exceed the
##' maximum indicated by Q2. The structure shown in the Tz plot
##' indicate the analytical variance which is exactly independent of
##' the experimental design. The corresponding loading plot shows how
##' this structure is capture by the used ISs.
##' @title Plot a statistics for CRMN normalization model
##' @param x an \code{nFit} object
##' @param y not used
##' @param ... passed on to the scatter plot calls
##' @return nothing
##' @seealso \code{slplot}
##' @examples
##' data(mix)
##' nfit <- normFit(mix, "crmn", factors="type", ncomp=2)
##' plot(nfit)
##' @export
##' @method plot nFit
##' @aliases plot.nFit plot,nFit-method
##' @author Henning Redestig
plot.nFit <- function(x, y=NULL,...) {
  if(!method(x) == "crmn")
    stop("Can only plot CRMN normalization models")
  pcMod <- sFit(x)$fit$pc
  if(!is.null(sFit(x)$q2)) {
    par(mfrow=c(1,3))
    mindim <- min(length(drop(sFit(x)$q2)), nPcs(sFit(x)$fit$pc))
    xx <- rbind(drop(sFit(x)$q2)[1:mindim],
                drop(cumsum(sFit(x)$r2))[1:mindim])
    barplot(xx, beside=TRUE,
               sub=expression(T[z] ~~ optimization), ylim=c(0,1.1))
    legend(x="topleft", fill=c("white", "grey"),
           legend=c(expression(R^2),
             expression(Q^2)))
  }
  else
    par(mfrow=c(1,2))    
  plot(scores(pcMod)[,1:2], xlab=expression(T[z1]), ylab=expression(T[z2]),...)
  plot(loadings(pcMod)[,1:2], xlab=expression(P[z1]), ylab=expression(P[z2]),...)
}
setMethod("plot", "nFit", plot.nFit)

##' Subset an data set to only contain the labeled internal standards.
##'
##' @title Accessor for the Internal Standards
##' @param object an \code{ExpressionSet}
##' @param where Column index or name in fData which equals
##' \code{what} for the ISs
##' @param what What the column \code{where} equals for ISs
##' @param ... not used
##' @return subsetted dataset
##' @aliases standards_eset standards,ExpressionSet,missing-method 
##' @examples
##' data(mix)
##' standards(mix)
##' fData(mix)$test <- fData(mix)$tag
##' standards(mix, where="test")
##' @author Henning Redestig
standards_eset <- function(object, where="tag", what="IS", ...) {
  if(is.null(fData(object)[,where]))
    stop(paste("No column named", where, "in the feature data"))
  if(all(!fData(object)[,where] %in% what))
    stop(paste("No rows tagged as", what))
  object[fData(object)[,where] %in% what,]
}

##' Subset an data set to only contain the labeled internal standards.
##'
##' @title Accessor for the Internal Standards
##' @param object an \code{matrix} or \code{data.frame}
##' @param standards a logical vector indicating which rows are internal standards
##' @param ... not used
##' @return subsetted dataset
##' @aliases standards_other standards,matrix,logical-method
##' standards,data.frame,logical-method 
##' @examples
##' data(mix)
##' standards(exprs(mix), fData(mix)$tag == 'IS')
##' @author Henning Redestig
standards_other <- function(object, standards, ...) {
  if(all(!standards))
    stop("No standards")
  object[standards,,drop=FALSE]
}

##' Subset an expression set to remove the internal standards
##'
##' @title Accessor for the analytes
##' @param object an \code{ExpressionSet}
##' @param where Column index or name of fData which equals
##' \code{what} for the ISs (and something else for the analytes)
##' @param what What the column \code{where} does not equal for
##' analytes. Can be vector values too.
##' @param ... not used
##' @aliases analytes_eset analytes,ExpressionSet,missing-method 
##' @return \code{ExpressionSet}
##' @examples
##' data(mix)
##' analytes(mix)
##' fData(mix)$test <- fData(mix)$tag
##' analytes(mix, where="test")
##' @author Henning Redestig 
analytes_eset<- function(object, where="tag", what="IS", ...) {
  if(is.null(fData(object)[,where]))
    stop(paste("No column named", where, "in the feature data"))
  chosen <- !fData(object)[,where] %in% what
  if(all(!chosen))
    stop(paste("All rows are tagged as", what))
  object[chosen,]
}

##' Subset an expression set to remove the internal standards
##'
##' @title Accessor for the analytes
##' @param object an \code{ExpressionSet}
##' @param standards a logical vector indicating which rows are
##' internal standards
##' @param ... not used
##' @aliases analytes_other analytes,matrix,logical-method
##' analytes,data.frame,logical-method 
##' @return \code{ExpressionSet}
##' @examples
##' data(mix)
##' analytes(exprs(mix), fData(mix)$tag == 'IS')
##' @author Henning Redestig
analytes_other<- function(object, standards, ...)
  object[!standards,,drop=FALSE]

##' Construct a design matrix 
##'
##' @title Make X
##' @param object an \code{ExpressionSet}
##' @param factors column names from the pheno data of \code{object}
##' @param ... not used
##' @return a design matrix
##' @aliases makeX_eset makeX,ExpressionSet,character-method
##' @examples
##' data(mix)
##' makeX(mix, "runorder")
##' @author Henning Redestig
##' @noRd
makeX_eset <- function(object, factors, ...) {
  x <- pData(object)[,factors,drop=FALSE]
  ## construct a Y matrix (the experiment related information)
  nm <- sapply(x, class) %in% c("numeric", "integer")
  fac <- x[,!nm,drop=FALSE]
  mod <- NULL
  for(i in 1:ncol(fac))
    if(length(levels(fac[,i])) > 1)
      mod <- cbind(mod, model.matrix(~-1+fac[,i]))
  num <- x[,nm,drop=FALSE]
  mod <- cbind(mod, num)
  as.matrix(mod)
}  

##' Construct a design matrix 
##'
##' Convenience function that just return the given design matrix. 
##' @title Make X
##' @param object, not used
##' @param factors a design matrix
##' @param ... not used
##' @return the same design matrix
##' @aliases makeX_other makeX,ANY,matrix-method
##' @examples
##' data(mix)
##' makeX(mix, model.matrix(~pData(mix)[,"runorder"]))
##' @author Henning Redestig
##' @noRd
makeX_other <- function(object, factors, ...)
  factors

##' Drop unused factor levels in a data frame.
##'
##' @title Drop unused levels
##' @param x the data frame
##' @author Henning Redestig
##' @examples
##' iris[1:10,]$Species
##' dropunusedlevels(iris[1:10,])$Species
##' @export
dropunusedlevels <- function (x)  {
    if (!is.data.frame(x)) 
        stop("only data frames")
    for (i in 1:length(x)) if (is.factor(x[, i])) 
        x[, i] <- factor(x[, i])
    x
}
