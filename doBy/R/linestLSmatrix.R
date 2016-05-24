## #######################################################################
##
## Construction of LSmatrix and computation of LSmeans
##
## Author: Søren Højsgaard
##
## #######################################################################


LSmatrix <- function(object, effect=NULL, at=NULL){
  UseMethod("LSmatrix")
}

## FIXME: LSmatrix.default: Should be a check of what 'object' is
LSmatrix.default <- function(object, effect=NULL, at=NULL){
    res <- .get_linest_list( object, effect, at )
    res <- .finalize_linest_list ( res )
    class(res) <- c("LSmatrix", "matrix")
    res
}

.finalize_linest_list <- function (aa){
    res               <- lapply( aa, function( mm ) apply( mm, 2, mean ) )
    res               <- do.call(rbind, res)
    attr(res, "at")   <- attr(aa, "at")
    attr(res, "grid") <- attr(aa, "grid")
    attr(res, "offset") <- attr(aa, "offset")
    res
}

print.LSmatrix <- function(x,...){
  prmatrix(x)
  ## atr <- attributes(x)[c("at","grid")]
  ## aa <- !unlist(lapply(atr, is.null))
  ## str(atr[aa])
  invisible(x)
}

.get_linest_list <- function(object, effect=NULL, at=NULL){
    ##cat(".get_linest_list\n")
    trms     <- delete.response( terms(object) )
    fact.lev <- get_xlevels( object )            ## factor levels
    ##cat("fact.lev:\n"); print(fact.lev)
    cov.ave  <- .get_covariate_ave( object, at )  ## average of covariates (except those mentioned in 'at')
    ##cat("cov.ave:\n"); print(cov.ave)
    vartype  <- get_vartypes( object )           ## which are factors and which are numerics
    ##cat("vartype:\n"); print(vartype)
    at.factor.name <- intersect( vartype$factor, names(at) )
    cov.ave.name   <- names( cov.ave )
    effect         <- setdiff( effect, at.factor.name )

    #' tmp <- list(fact.lev=fact.lev, cov.ave=cov.ave, vartype=vartype, at.factor.name=at.factor.name, cov.ave.name=cov.ave.name, effect=effect, at=at)
    #' print(tmp)

    if (is.null(effect)){
        if (length( at.factor.name ) > 0){
            new.fact.lev <- at[ at.factor.name ]
        } else {
            new.fact.lev <- NULL
        }
    } else {
        new.fact.lev  <- set_xlevels( fact.lev, at=at )
        new.fact.lev  <- new.fact.lev[c(effect, at.factor.name)]#
    }
    if (is.null(new.fact.lev)){
        ##cat("No 'effect' and no 'at'-factors; hence just a global average... \n")
        ## print(fact.lev)
        ## print(cov.ave.name)

        if ( length(fact.lev) > 0 ){
            ##cat("yes there are factors\n")
            newdata <- expand.grid( fact.lev )
            if (length( cov.ave.name ) > 0){
                ##cat("yes there are covariates\n")
                newdata[, cov.ave.name] <- cov.ave
            }
        } else {
            if (length( cov.ave.name ) > 0){
                ##cat("yes there are covariates\n")
                newdata <- matrix(unlist(cov.ave), nrow=1L)
                colnames(newdata) <- cov.ave.name
                newdata <- as.data.frame( newdata )
            } else {
                ##cat("there are no factors or covariates\n")
                newdata <- data.frame(1)
            }
        }

        XXlist <- list(get_X(object, newdata))
        ## cat("XXlist:\n"); print(XXlist)
        attr(XXlist,"at")   <- at[intersect(vartype$numeric, names(at))]
        attr(XXlist,"grid") <- NULL
    } else {
        ##cat("The general case; there are 'effect' factors or 'at' factors...\n")
        grid.data <- expand.grid(new.fact.lev)
        grid.data <- as.data.frame(lapply(grid.data, as.character), stringsAsFactors=FALSE)
        XXlist    <- list()
        for (ii in 1:nrow(grid.data)){
            config    <- grid.data[ ii, ,drop=FALSE ]
            fact.lev2 <- set_xlevels(fact.lev,  at=config)

            newdata   <- expand.grid( fact.lev2 )
            newdata[, cov.ave.name]  <- cov.ave
            XX             <- get_X(object, newdata, at)
            XXlist[[ ii ]] <- XX
        }

        grid.data[, names(cov.ave) ] <- cov.ave
        attr(XXlist,"at") <- at
        attr(XXlist,"grid") <- grid.data
        attr(XXlist,"offset") <- attr(XX, "offset")
    }
    class(XXlist) <- "linestList"
    XXlist
}

## --------------------------------------------------------------------

LSmeans <- function(object, effect=NULL, at=NULL, level=0.95,...){
    UseMethod("LSmeans")
}

LSmeans.default <- function(object, effect=NULL, at=NULL, level=0.95,...){
    K   <- LSmatrix(object, effect=effect, at=at)
    out <- linest(object, K, level=level, ...)
    out
}

LSmeans.lmerMod <- function(object, effect=NULL, at=NULL, level=0.95, adjust.df=TRUE, ...){
    K   <- LSmatrix(object, effect=effect, at=at)
    out <- linest(object, K, level=level, adjust.df=adjust.df, ...)
    out
}

popMeans         <- LSmeans
popMeans.default <- LSmeans.default
popMeans.lmerMod <- LSmeans.lmerMod

## --------------------------------------------------------------------


setOldClass("LSmatrix")

setAs("LSmatrix","matrix",
      function(from){
          attr(from,"at")<- attr(from,"grid")<-NULL
          class(from)<-"matrix"
          from
      })

setAs("LSmatrix","Matrix",
      function(from){
          attr(from,"at")<- attr(from,"grid")<-NULL
          class(from)<-"matrix"
          from
          as(from,"Matrix")
      })

