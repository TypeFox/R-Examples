
lm <- function(formula, ...){
UseMethod("lm")
}

lm.default <- stats::lm

lm.design <- function (formula, ..., response = NULL, degree = NULL, FUN = mean, 
    use.center=NULL, use.star=NULL, use.dummies=FALSE){
    if (!"design" %in% class(formula)) 
        stop("lm.design works on class design objects only")
    
    di <- design.info(formula)
    ## capture designs from conf.design
    if (is.null(di)) stop("lm.design does not work for class design from package conf.design")
    if (is.null(use.center)) if (di$type=="ccd") use.center <- TRUE else use.center <- FALSE
    if (is.null(use.star)) if (di$type=="ccd") use.star <- TRUE else use.star <- FALSE

    fo <- formula(formula, ..., response = response, degree = degree, 
        FUN = deparse(substitute(FUN)), use.center=use.center, use.star=use.star, use.dummies=use.dummies)
    if (di$repeat.only | (length(grep("param", di$type)) > 0 & 
        length(grep("wide", di$type)) == 0) | (length(grep("center", di$type)) > 0 & !use.center & !use.star)) 
        aus <- lm(fo, data = model.frame(fo, data = NULL), ...)
    else aus <- lm(fo, data = model.frame(fo, data = formula), ...)
    class(aus) <- c("lm.design",class(aus))
    
        if (length(grep("splitplot", di$type)) > 0) {
        mm <- model.matrix(aus)
        coefs <- coef(aus)[-1]
        nms <- names(coefs)
        hilf <- apply(mm[, 1 + (1:di$nfac.WP), drop = FALSE], 
            1, paste, collapse = "")
        pchs <- rep("*", length(coefs))
        pchs[1:di$nfac.WP] <- "o"
        for (j in setdiff(1:(length(coefs)), 1:di$nfac.WP)) {
            if (!length(table(paste(hilf, mm[, nms[j]], sep = ""))) > 
                di$nWPs) 
                pchs[j] <- "o"
        }
        aus$WholePlotEffects <- names(coefs[pchs=="o"])
    }
    aus
}

summary.lm.design <- function(object, ...){
   aus <- summary.lm(object)
   class(aus) <- c("summary.lm.design", class(aus))
   ## cater for split plot design warning
   aus$WholePlotEffects <- object$WholePlotEffects
   aus
}

print.summary.lm.design <- function(x, ...){
    cat("Number of observations used:", sum(x$df[1:2]),"\n")
    cat("Formula:\n")
    fop <- formula(x)
    attributes(fop) <- NULL 
    print(fop)
    
    ## printing x with the summary.lm method
    NextMethod(x, ...)
        
    if (!is.null(x$WholePlotEffects)){
       cat("WARNING: This is a split plot design, whole plot effects may have larger variance!\n")
       cat("         p-values for whole plot effects may be misleadingly low!\n")
       cat("The whole plot effects are:\n")
       print(x$WholePlotEffects, quote=FALSE)
       }
}

print.lm.design <-function(x, ...){
    cat("Number of observations used:", nrow(model.frame(x)),"\n")
    cat("Formula:\n")
    fop <- formula(x)
    attributes(fop) <- NULL 
    print(fop)
    
    ## printing x with the lm method
    NextMethod(x, ...)
    
    if (!is.null(x$WholePlotEffects)){
       cat("WARNING: This is a split plot design, whole plot effects may have larger variance!\n")
       cat("         p-values for whole plot effects may be misleadingly low!\n")
       cat("The whole plot effects are:\n")
       print(x$WholePlotEffects, quote=FALSE)
       }
}

coef.lm.design <- function(object, ...){
   class(object) <- c("aov", class(object))
   coef(object, ...)
   }