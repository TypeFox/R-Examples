"gvlma.form" <-
function(formula, data, alphalevel = 0.05, timeseq = 1:nrow(data), ...)
  {
    cl <- match.call()
    lmc <- cl
    lmc$alphalevel <- NULL
    lmc$timeseq <- NULL
    lmc[[1]] <- as.name("lm")
    lmobj <- eval(lmc, parent.frame())
    # determine whether "subset" is provided
    extras <- match.call(expand.dots = FALSE)$...
    subsetgiven <- FALSE
    if (length(extras)>0)
      {
        subsetindex <- match("subset", names(extras))
        if (subsetindex > 0) subsetgiven <- TRUE
      }

    n <- length(fitted(lmobj))
    # get reasonable timeseq
    goodlength <- length(timeseq) == n
    if (missing(timeseq))
      v <- 1:n
    else
      {
        if (!subsetgiven)
          {
            if (!goodlength)
              stop("Argument timeseq of incorrect length.")
            else
              v <- timeseq
          }
        else
          if (goodlength) v <- timeseq
          else v <- timeseq[eval(extras[[subsetindex]])]
      }

    ##
    z <- computegvlma(lmobj, alphalevel, v)
    z$call <- lmobj$call
    ##
    z$GlobalTest$call <- cl
    z
  }

