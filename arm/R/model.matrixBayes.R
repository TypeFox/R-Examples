#setMethod("model.matrix.bayes", signature(object = "bayesglm"),
model.matrixBayes <- function(object, data = environment(object),
        contrasts.arg = NULL, xlev = NULL, keep.order=FALSE, drop.baseline=FALSE,...)
{
    #class(object) <- c("terms", "formula")
    t <- if( missing( data ) ) {
          terms( object )
         }else{
            terms.formula(object, data = data, keep.order=keep.order)
         }
    attr(t, "intercept") <- attr(object, "intercept")
    if (is.null(attr(data, "terms"))){
      data <- model.frame(object, data, xlev=xlev)
    }else {
        reorder <- match(sapply(attr(t,"variables"), deparse, width.cutoff=500)[-1], names(data))
        if (anyNA(reorder)) {
            stop( "model frame and formula mismatch in model.matrix()" )
        }
        if(!identical(reorder, seq_len(ncol(data)))) {
            data <- data[,reorder, drop = FALSE]
        }
    }
    int <- attr(t, "response")
    if(length(data)) {      # otherwise no rhs terms, so skip all this

        if (drop.baseline){
          contr.funs <- as.character(getOption("contrasts"))
        }else{
          contr.funs <- as.character(list("contr.bayes.unordered", "contr.bayes.ordered"))
        }

        namD <- names(data)
        ## turn any character columns into factors
        for(i in namD)
            if(is.character( data[[i]] ) ) {
                data[[i]] <- factor(data[[i]])
                warning( gettextf( "variable '%s' converted to a factor", i ), domain = NA)
            }
        isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)
        isF[int] <- FALSE
        isOF <- vapply(data, is.ordered, NA)
        for( nn in namD[isF] )            # drop response
            if( is.null( attr( data[[nn]], "contrasts" ) ) ) {
                contrasts( data[[nn]] ) <- contr.funs[1 + isOF[nn]]
            }
        ## it might be safer to have numerical contrasts:
        ##    get(contr.funs[1 + isOF[nn]])(nlevels(data[[nn]]))
        if ( !is.null( contrasts.arg ) && is.list( contrasts.arg ) ) {
            if ( is.null( namC <- names( contrasts.arg ) ) ) {
                stop( "invalid 'contrasts.arg' argument" )
            }
            for (nn in namC) {
                if ( is.na( ni <- match( nn, namD ) ) ) {
                    warning( gettextf( "variable '%s' is absent, its contrast will be ignored", nn ), domain = NA )
                }
                else {
                    ca <- contrasts.arg[[nn]]
                    if( is.matrix( ca ) ) {
                        contrasts( data[[ni]], ncol( ca ) ) <- ca
                    }
                    else {
                        contrasts( data[[ni]] ) <- contrasts.arg[[nn]]
                    }
                }
            }
        }
    } else {               # internal model.matrix needs some variable
        isF  <-  FALSE
        data <- data.frame(x=rep(0, nrow(data)))
    }
    #ans  <- .Internal( model.matrix( t, data ) )
    ans  <- model.matrix.default(object=t, data=data)
    cons <- if(any(isF)){
              lapply( data[isF], function(x) attr( x,  "contrasts") )
            }else { NULL }
    attr(ans, "contrasts" ) <- cons
    ans
}
#)

#setMethod("model.matrix.bayes", signature(object = "bayesglm.h"),
#model.matrix.bayes.h <- function (object, data = environment(object),
#            contrasts.arg = NULL,
#            xlev = NULL, keep.order = FALSE, batch = NULL, ...)
#{
#    class(object) <- c("formula")
#    t <- if (missing(data)) {
#        terms(object)
#    }
#    else {
#        terms(object, data = data, keep.order = keep.order)
#    }
#    attr(t, "intercept") <- attr(object, "intercept")
#    if (is.null(attr(data, "terms"))) {
#        data <- model.frame(object, data, xlev = xlev)
#    }
#    else {
#        reorder <- match(sapply(attr(t, "variables"), deparse,
#            width.cutoff = 500)[-1], names(data))
#        if (any(is.na(reorder))) {
#            stop("model frame and formula mismatch in model.matrix()")
#        }
#        if (!identical(reorder, seq_len(ncol(data)))) {
#            data <- data[, reorder, drop = FALSE]
#        }
#    }
#    int <- attr(t, "response")
#    if (length(data)) {
#        contr.funs <- as.character(getOption("contrasts"))
#        contr.bayes.funs <- as.character(list("contr.bayes.unordered",
#            "contr.bayes.ordered"))
#        namD <- names(data)
#        for (i in namD) if (is.character(data[[i]])) {
#            data[[i]] <- factor(data[[i]])
#            warning(gettextf("variable '%s' converted to a factor", i), domain = NA)
#        }
#        isF <- sapply(data, function(x) is.factor(x) || is.logical(x))
#        isF[int] <- FALSE
#        isOF <- sapply(data, is.ordered)
#        if (length(batch) > 1) {
#            ba <- batch[isF[-1]]
#        }
#        else if (length(batch) == 1) {
#            ba <- rep(batch, length(isF[-1]))
#        }
#        else {
#            ba <- rep(0, length(isF[-1]))
#        }
#        iin <- 1
#        for (nn in namD[isF]) if (is.null(attr(data[[nn]], "contrasts"))) {
#            if (ba[[iin]] > 0) {
#                contrasts(data[[nn]]) <- contr.bayes.funs
#            }
#            else {
#                contrasts(data[[nn]]) <- contr.funs
#            }
#            iin <- iin + 1
#        }
#        if (!is.null(contrasts.arg) && is.list(contrasts.arg)) {
#            if (is.null(namC <- names(contrasts.arg))) {
#                stop("invalid 'contrasts.arg' argument")
#            }
#            for (nn in namC) {
#                if (is.na(ni <- match(nn, namD))) {
#                  warning(gettextf("variable '%s' is absent, its contrast will be ignored",
#                    nn), domain = NA)
#                }
#                else {
#                  ca <- contrasts.arg[[nn]]
#                  if (is.matrix(ca)) {
#                    contrasts(data[[ni]], ncol(ca)) <- ca
#                  }
#                  else {
#                    contrasts(data[[ni]]) <- contrasts.arg[[nn]]
#                  }
#                }
#            }
#        }
#    }
#    else {
#        isF <- FALSE
#        data <- list(x = rep(0, nrow(data)))
#    }
#    ans <- .Internal(model.matrix(t, data))
#    cons <- if (any(isF)) {
#        lapply(data[isF], function(x) attr(x, "contrasts"))
#    }
#    else {
#        NULL
#    }
#    attr(ans, "contrasts") <- cons
#    ans
#}
##)
