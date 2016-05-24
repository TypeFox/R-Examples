"quagep" <-
function(f, para, paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
       if(! are.pargep.valid(para)) return()
    }
    attributes(para$para) <- NULL
    B <- para$para[1]
    K <- para$para[2]
    H <- para$para[3]

    ix <- seq(1:length(f))
    ops <- options(warn=-1)
    x <- -B * log(1 + (1/H) * log(1 - f^(1/K) * (1-exp(-H)) ) )
    for(i in ix[is.nan(x)]) {
       warning("The ",i,"(th) value of 'f' results in NaN (assuming then f == 1), ",
               "decrementing from the Machine's small to an f that just hits non NaN")
       j <- 0
       while(1) {
          j <- j + 1
          aF <- 1 - .Machine$double.eps^(1/j)
          aX <- -B * log(1 + (1/H) * log(1 - aF^(1/K) * (1-exp(-H)) ) )
          if(! is.nan(aX)) {
             x[i] <- aX
             break
          }
       }
    }
    options(ops)
    return(x)
}

