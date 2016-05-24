"summaryNAltraj" <-
function(x)
  {
    if (!inherits(x,"ltraj"))
      stop("x should be of class 'ltraj'")
    nna <- unlist(lapply(x, function(i) length(i$x[is.na(i$x)])))
    n <- unlist(lapply(x, function(i) length(i$x)))
    ani <- unlist(lapply(x, function(i) attr(i, "id")))
    burst <- unlist(lapply(x, function(i) attr(i, "burst")))
    so <- data.frame(id=ani, burst=burst,
                     missing.values=nna, N.reloc=n, percent = 100*nna/n)
    return(so)
  }

