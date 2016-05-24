mvPCA <- function(X, score = "identity", estimate= "outer", na.action=na.fail, ...)
    {
    DNAME <- deparse(substitute(X))
    score <- match.arg(score,c("identity","sign","symmsign","rank"))
    estimate <- match.arg(estimate,c("outer","inner"))
    
    X<-na.action(X)
    
    
    
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    n<-dim(X)[1]
    p<-dim(X)[2]
    #if (is.null(colnames(X))) colnames(X) <- paste(rep("X",p),1:p,sep=".")
    #c.names <- colnames(X)
    
    covM <- switch(score,
              "identity" = {
                      METHOD <- "regular covariance matrix"  
                      (n-1) /n * cov(X) 
                      }
              ,
              "sign" = {
                       switch(estimate,
                                "outer" =
                                {
                                METHOD <- "SCOV"
                                SCov(X,spatial.median(X))
                                }
                                ,
                                "inner" =
                                {
                                METHOD <- "Tyler's shape matrix"
                                HR.Mest(X,...)$scatter
                                }
                                )
                       } 
              ,
              "symmsign" = {
                           switch(estimate,
                                "outer" =
                                {
                                METHOD <- "TCOV"
                                SSCov(X)
                                }
                                ,
                                "inner" = 
                                {
                                METHOD <- "Duembgen's shape matrix"
                                duembgen.shape(X,...)
                                }
                                )
                           } 
              ,
              "rank" = {
                           switch(estimate,
                                "outer" = 
                                {
                                METHOD <- "RCOV"
                                RCov(X)
                                }
                                ,
                                "inner" = 
                                {
                                METHOD <- "rank shape matrix"
                                rank.shape(X,...)
                                }
                                )
                           } 
                      )
    
    covM <- covM / sum(diag(covM)) * p
    
    edc <- eigen(covM, symmetric = TRUE)
    ev <- edc$values
    if (any(neg <- ev < 0)) {
        if (any(ev[neg] < -9 * .Machine$double.eps * ev[1])) 
            stop("covariance matrix is not non-negative definite")
        else ev[neg] <- 0}
    cn <- paste("Comp.", 1:p, sep = "")
    names(ev) <- cn
    loadings <- edc$vectors
    colnames(loadings) <- cn
    if (!is.null(colnames(X))) 
        {
        c.names <- colnames(X)
        rownames(loadings) <- c.names
        }
    scores <- X %*% edc$vectors
    colnames(scores) <- cn
    
    res <- list(EigenV = ev, loadings=loadings, scores=scores, dname=DNAME, method=METHOD, n.obs = n, p = p)
    class(res) <- "mvPCA"
    res
    }


print.mvPCA <- function(x,...)
    {
    cat("PCA for",x$dname,"based on", x$method, "\n")
    cat("\nStandardized eigenvalues:\n")
    print(x$EigenV,...)
    cat("\n", x$p, "variables and", x$n.obs, "observations.\n")
    invisible(x)
    }
    
summary.mvPCA <-
function (object, loadings = FALSE, cutoff = 0.1, ...) 
{
    object$cutoff <- cutoff
    object$print.loadings <- loadings
    class(object) <- "summary.mvPCA"
    object
}

print.summary.mvPCA <-
function (x, digits = 3, loadings = x$print.loadings, cutoff = x$cutoff, ...) 
{
    vars <- x$EigenV / x$p
    cat("Importance of components:\n")
    print(rbind("Proportion of Variation" = vars, 
        "Cumulative Proportion" = cumsum(vars)))
    if (loadings) {
        cat("\nLoadings:\n")
        cx <- format(round(x$loadings, digits = digits))
        cx[abs(x$loadings) < cutoff] <- paste(rep(" ", nchar(cx[1, 
            1], type = "w")), collapse = "")
        print(cx, quote = FALSE, ...)
    }
    invisible(x)
}

screeplot.mvPCA <-
function (x, npcs = min(10, length(x$EigenV)), type = c("barplot", 
    "lines"), main = deparse(substitute(x)), ...) 
{
    main
    type <- match.arg(type)
    pcs <- x$EigenV / x$p
    xp <- seq_len(npcs)
    if (type == "barplot") 
        barplot(pcs[xp], names.arg = names(pcs[xp]), main = main, 
            ylab = "Proportion of variation", ...)
    else {
        plot(xp, pcs[xp], type = "b", axes = FALSE, main = main, 
            xlab = "", ylab = "Variances", ...)
        axis(2)
        axis(1, at = xp, labels = names(pcs[xp]))
    }
    invisible()
}

plot.mvPCA <-
function (x, main = deparse(substitute(x)), ...) 
screeplot.mvPCA(x, main = main, ...)

predict.mvPCA <-
function (object, newdata, ...) 
{
    if (missing(newdata)) 
        return(object$scores)
    if (length(dim(newdata)) != 2) 
        stop("'newdata' must be a matrix or data frame")
    p <- object$p
    nm <- rownames(object$loadings)
    if (!is.null(nm)) {
        if (!all(nm %in% colnames(newdata))) 
            stop("'newdata' does not have named columns matching one or more of the original columns")
        newdata <- newdata[, nm]
    }
    else {
        if (NCOL(newdata) != p) 
            stop("'newdata' does not have the correct number of columns")
    }
    newdata %*% object$loadings
}

#biplot.mvPCA <-
#function (x, choices = 1:2, ...) 
#{
#    if (length(choices) != 2) 
#        stop("length of choices must be 2")
#   
#    stats:::biplot.default(x$scores[, choices], x$loadings[, choices], ...)
#    invisible()
#}
