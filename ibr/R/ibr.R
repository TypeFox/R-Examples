ibr <- function(formula,data,subset,criterion="gcv",df=1.5,Kmin=1,Kmax=1e+06,smoother="k",kernel="g",rank=NULL,control.par=list(),cv.options=list()) {
    cl <- match.call() 
    mf <- match.call()
    m <- match(c("formula", "data", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- "na.fail"
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    attr(mt,"intercept") <- 0
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) stop("offset are not allowed")
    if (is.empty.model(mt)) stop("no model ?")
    if (any(apply(mf,2,is.factor)))  stop("no factor are allowed")
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf)
    attributes(x) <- attributes(x)[c("dim","dimnames")]
    res <- ibr.fit(x,y,criterion,df,Kmin,Kmax,smoother,kernel,rank,control.par,cv.options)
    res$call <- cl
    res$terms <- mt
    class(res) <- c("ibr", "list")
   return(res)
}
