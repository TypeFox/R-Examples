"gvlma.lm" <-
function(lmobj, alphalevel = 0.05, timeseq)
  {
    if (class(lmobj) != "lm")
      stop("First argument to gvlma.lm must be an lm object.")
    cl <- match.call()
    n <- length(fitted(lmobj))
    # check timeseq
    if (missing(timeseq)) v <- 1:n else v <- timeseq
    if (length(v) != n) stop("Argument timeseq of incorrect length.")
    #
    z <- computegvlma(lmobj, alphalevel, v)
    #
    z$GlobalTest$call <- cl
    z
  }

