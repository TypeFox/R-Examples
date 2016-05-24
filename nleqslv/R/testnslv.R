
testnslv <- function(x, fn, jac=NULL, ...,
                            method=c("Newton", "Broyden"),
                            global=c("cline", "qline", "gline", "pwldog", "dbldog", "hook", "none"),
                            Nrep=0L, title=NULL
                          )
{
    # utility functions
    catmsg <- function(m,g,res) {
        cat(sprintf("Error (method=%s global=%s): %s\n",m,g,attr(res,"condition")$message))
    }

    makeerrlist <- function(m,g,cpusecs=NULL) {
        if(is.null(cpusecs)) {
            z <- list(Method=m, Global=g, termcd=NA, Fcnt=NA, Jcnt=NA, Iter=NA, Message="ERROR",Fnorm=NA)
        } else {
            z <- list(Method=m, Global=g, termcd=NA, Fcnt=NA, Jcnt=NA, Iter=NA, Message="ERROR",Fnorm=NA,
                             cpusecs=cpusecs)
        }
        z
    }

    makereslist <- function(m,g,res,cpusecs=NULL) {
        fnorm <- sum(res$fvec^2)/2
        # necessary to test for termcd<0 and >6 otherwise R errors later in output
        # see R-help about switch
        if(res$termcd < 0 ) stop("User supplied jacobian most likely incorrect: cannot continue") else
        if(res$termcd > 7 ) message <- "BADCD" else
            message <- switch(res$termcd, "Fcrit", "Xcrit", "Stalled", "Maxiter", "Illcond", "Singular", "BadJac")

        if(is.null(cpusecs)) {
           z <- list(Method=m, Global=g, termcd=res$termcd, Fcnt=res$nfcnt, Jcnt=res$njcnt,
                            Iter=res$iter, Message=message, Fnorm=fnorm)
        } else {
           z <- list(Method=m, Global=g, termcd=res$termcd, Fcnt=res$nfcnt, Jcnt=res$njcnt,
                            Iter=res$iter, Message=message, Fnorm=fnorm, cpusecs=cpusecs)
        }
        z
    }

    methods <- match.arg(method, c("Newton", "Broyden"), several.ok=TRUE)
    globals <- match.arg(global, c("cline", "qline", "gline", "pwldog", "dbldog", "hook", "none"), several.ok=TRUE)

    my.call <- match.call()
    reslist <- vector("list", length(methods)*length(globals))

    # use try to avoid process stopping for Jacobian with non-finite values
    # if that happens, go to next method/global
    # avoidable fatal user errors will also lead to useless next method/global
    idx <- 1
    for(m in methods)
        for(g in globals) {
            if( Nrep >= 1) {
                mytime <- system.time( for(k in seq_len(Nrep)) {
                                res <- try(nleqslv(x, fn, jac, ..., method=m, global=g), silent=TRUE)
                                if(inherits(res,"try-error")) break
                            }, gcFirst = FALSE)
                cpus <- mytime[3]
            } else {
                res <- try(nleqslv(x, fn, jac, ..., method=m, global=g),silent=TRUE)
                cpus <- NULL
            }
            if(inherits(res,"try-error")) {
                catmsg(m,g,res)
                z <- makeerrlist(m,g,cpus)
            } else {
                z <- makereslist(m,g,res,cpus)
            }
            reslist[[idx]] <- z
            idx <- idx+1
        }

# from http://stackoverflow.com/questions/4512465/what-is-the-most-efficient-way-to-cast-a-list-as-a-data-frame?rq=1

    ## @Martin Morgan's Map() sapply() solution:
    f <- function(x) function(i) sapply(x, `[[`, i)
    z <- as.data.frame(Map(f(reslist), names(reslist[[1]])), stringsAsFactors=FALSE)

    res <- list()
    res$out <- z
    res$call <- my.call
    res$title <- title
    class(res) <- "test.nleqslv"
    res
}

print.test.nleqslv <- function(x, digits=4, width.cutoff=45L, ...) {
    if(!inherits(x, "test.nleqslv"))
        stop("method is only for test.nleqslv objects")

    # calculate total number of function evaluations if numeric jacobian used

    cat("Call:\n",paste0(deparse(x$call, width.cutoff=width.cutoff), collapse = "\n"), "\n\n", sep = "")
    if(is.null(x$title)) cat("Results:\n") else cat("Results: ",x$title,"\n", sep="")
    print(x$out, digits=digits,...)
    invisible(x)
}
