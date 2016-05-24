boot.relimp.formula <- function(formula, data, weights, na.action, ..., subset=NULL){
    if (missing(formula)) stop("formula missing")
    if (missing(na.action)) 
        na.action <- getOption("na.action")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())

    terms <- attr(mf,"terms")
    resp <- attr(terms,"response")
    if (resp != 1 ) stop("incorrect formula") 
    if (max(attr(terms,"order")) > 2) stop ("formula contains terms of order higher than 2")
    if (attr(terms,"intercept") != 1) stop ("model must contain intercept")

    if (!is.null(attr(mf,"na.action"))) warning(naprint(attr(mf,"na.action")))
    if (!is.null(dim(model.response(mf)))){
          if (ncol(model.response(mf))>1) stop("too many response variables")
          }
    ## selection of columns from model needed because of e.g. lm(y~x1+x2+x3-x2)
    ## selection of columns from model based on formula below 
    ##               works even in case of multi-column terms such as poly(x2,3)

    #DATA <- as.matrix(mf[,c(resp,which(rowSums(attr(terms,"factors"))>0))])
    #weights <- mf$"(weights)"
    
    bt <- do.call("boot.relimp", list(lm(mf), ...))
    bt
}
