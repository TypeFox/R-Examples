"cdfgld" <-
function(x,para,paracheck=TRUE) {

    # Check that the parameters are valid one time
    # then use the paracheck switch on quagld for
    # an extreme speed up on this algorithm.
    if(paracheck == TRUE) {
      if(! are.pargld.valid(para)) return()
    }

    sup <- supdist(par, paracheck=FALSE)
    lo  <- sup$support[1]
    hi  <- sup$support[2]
    flo <- sup$nonexceeds[1]
    fhi <- sup$nonexceeds[2]

    f <- sapply(1:length(x), function(i) {
            QUAx <- x[i]
            if(QUAx <= lo) return(flo)
            if(QUAx >= hi) return(fhi)

            fn <- function(F) {
                      #print(F) # clearly 0,1 are first two but then a negative!
                      # These two traps seemingly are not needed, except uniroot
                      # if it sees [0,1] interval but at least for the upper
                      # limit is Inf then a negative nonexceedance gets fed in
                      # WHA   March 4, 2015
                      if(F < 0) F <- 0  # WHY?
                      if(F > 1) F <- 1  # WHY?
                      qua <- quagld(F,para,paracheck=FALSE)
                      val <- QUAx - qua
                return(val)
            }

            root <- uniroot(fn,c(flo,fhi), trace=10)
            return(root$root) })
    names(f) <- NULL
    return(f)
}

