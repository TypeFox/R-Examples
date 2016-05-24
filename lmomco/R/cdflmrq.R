"cdflmrq" <-
function(x,para, paracheck=FALSE) {
    if(! are.parlmrq.valid(para)) return()

    SMALL  <- .Machine$double.eps
    SMALLo <- 1 - SMALL
    hi <- qualmrq(SMALLo, para, paracheck=FALSE)

    "afunc" <- function(f, ax=NULL) return(qualmrq(f, para, paracheck=FALSE) - ax)
    f <- sapply(1:length(x), function(i) {
                  ax <- x[i]
                  if(ax == 0) return(0)
                  if(ax <  0) return(NA)
                  if(! is.finite(ax) | ax > hi) return(1)
                  rt <- NULL
                  try(rt <- uniroot(afunc, lower=SMALL, upper=SMALLo,
                                    ax=ax), silent=FALSE)
                  ifelse(is.null(rt), return(NA), return(rt$root))  })
    names(f) <- NULL
    return(f)
}



