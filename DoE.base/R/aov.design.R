aov <- function(formula, ...){
     UseMethod("aov")
}
aov.default <- stats::aov
aov.design <- function (formula, ..., response = NULL, degree = NULL, FUN = mean, 
    use.center = FALSE) 
{
    if (!"design" %in% class(formula)) 
        stop("aov.design works on class design objects only")
    di <- design.info(formula)
    ## capture designs from conf.design
    if (is.null(di)) 
    if (is.null(di)) stop("aov.design does not work for class design from package conf.design")

    fo <- formula(formula, ..., response = response, degree = degree, 
        FUN = deparse(substitute(FUN)), use.center = use.center)
    if (di$repeat.only | (length(grep("param", di$type)) > 0 & 
        length(grep("wide", di$type)) == 0) | (length(grep("center", 
        di$type)) > 0 & !use.center)) 
        aus <- aov(fo, data = model.frame(fo, data = NULL), ...)
    else aus <- aov(fo, data = model.frame(fo, data = formula), 
        ...)
    class(aus) <- c("aov.design", class(aus))

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

summary.aov.design <- function(object, ...){
   aus <- summary.aov(object)
   class(aus) <- c("summary.aov.design", class(aus))
   fop <- formula(terms(object))
   attributes(fop) <- NULL
   attr(aus, "formula") <- fop
   attr(aus, "WholePlotEffects") <- object$WholePlotEffects
   aus
}

print.summary.aov.design <- function(x, ...){
    cat("Number of observations used:", sum(x[[1]]$Df) + 1,"\n")
    cat("Formula:\n")
    print(attr(x, "formula"))
    
    ## printing x with the summary.aov method
    NextMethod(x, ...)   

    if (!is.null(attr(x, "WholePlotEffects"))){
       cat("WARNING: This is a split plot design, whole plot effects may have larger variance!\n")
       cat("         p-values for whole plot effects may be misleadingly low!\n")
       cat("The whole plot effects are:\n")
       print(attr(x, "WholePlotEffects"), quote=FALSE)
       }}

print.aov.design <-function(x, ...){
    cat("Number of observations used:", nrow(model.frame(x)),"\n")
    cat("Formula:\n")
    fop <- formula(x)
    attributes(fop) <- NULL 
    print(fop)
    ## printing x with the aov method
    NextMethod(x, ...)   

    if (!is.null(x$WholePlotEffects)){
       cat("WARNING: This is a split plot design, whole plot effects may have larger variance!\n")
       cat("         p-values for whole plot effects may be misleadingly low!\n")
       cat("The whole plot effects are:\n")
       print(x$WholePlotEffects, quote=FALSE)
       }}
