"cdfkap" <-
function(x,para) {
    if(! are.parkap.valid(para)) return()

    #  SMALL IS A SMALL NUMBER, USED TO TEST WHETHER X IS
    #  EFFECTIVELY AT AN ENDPOINT OF THE DISTRIBUTION
    SMALL <- 1e-15

    U <- para$para[1]
    A <- para$para[2]
    G <- para$para[3]
    H <- para$para[4]

    f <- sapply(1:length(x), function(i) {
                                Y <- (x[i]-U)/A
                                if(G == 0) {
                                   Y <- exp(-Y)
                                } else {
                                   ARG <- 1-G*Y
                                   if(ARG > SMALL) {
                                      Y <- exp(-1*(-log(ARG)/G))
                                   } else {
                                      if(G < 0) return(0)
                                      if(G > 0) return(1)
                                      stop("should not be here in execution")
                                   }
                                }
                                if(H == 0) return(exp(-Y))
                                ARG <- 1-H*Y
                                if(ARG > SMALL) return(exp(-1*(-log(ARG)/H)))
                                return(0) })
  names(f) <- NULL
  return(f)
}

