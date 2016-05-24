### partial fits                            
gamplot <- function(object, newdata = NULL) {

     if (is.null(newdata)) {
         x <- object$data$input
         pr <- function(obj) fitted(obj)
     } else {
         x <- newdata
         pr <- function(obj) predict(obj, newdata = x)
     }                                                      
     lp <- matrix(0, ncol = length(object$data$input),      
                  nrow = NROW(x[[1]]))
     ens <- object$ensemble
     ensss <- object$ensembless
     nu <- object$control$nu
     mstop <- nrow(ens)
     for (m in 1:mstop) {
         xselect <- ens[m,"xselect"]
         lp[,xselect] <- lp[,xselect] + nu * pr(ensss[[m]])
     }
     colnames(lp) <- colnames(object$data$input)
     lp
}

#output partial function estimation and the associated covariates
#modified from mboost library gamplot function
gamplot.out <- function(x, which = NULL, ask = FALSE && dev.interactive(),
    type = "b", ylab = expression(f[partial]), add_rug = TRUE, ...) {

    lp <- gamplot(x)
    input <- x$data$input
    ### <FIXME>: y ~ bbs(x) means that we only have access to x via
    ### the environment of its dpp function
    tmp <- lapply(input, function(x)
        eval(expression(x), envir = environment(attr(x, "dpp"))))
#note: this "tmp" is the original input "x" for the function bbs. 
#But for bspatial function, there are arguments x and y, then again
# "tmp" is only x, not y  
    input <- as.data.frame(tmp)
    names(input) <- names(tmp)
    ### </FIXME>
    if (is.null(which)) which <- (1:ncol(input))[tabulate(x$ensemble,
                                                         nbins = ncol(input)) > 0]
    if (is.numeric(which)) which <- names(input)[which]

    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }

   xp <- sapply(which, function(w) input[[w]])
   yp <- sapply(which, function(w) lp[,w])
   list(xp=xp,yp=yp)
   # out <- sapply(which, function(w) {
   #     xp <- input[[w]]
   #     yp <- lp[,w]
   #    ox <- order(xp)
    #    plot(xp[ox], yp[ox], xlab = w, type = type,
    #         ylab = ylab, ylim = range(lp[,which]), ...)
    #    abline(h = 0, lty = 3)
    #    if (add_rug) rug(input[[w]])
    #})
    #rm(out)
}
