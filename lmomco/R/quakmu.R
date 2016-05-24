"quakmu" <-
function(f, para, paracheck=TRUE, getmed=FALSE, qualo=NA, quahi=NA, verbose=FALSE,
                    marcumQ=TRUE, marcumQmethod=c("chisq", "delta", "integral")) {

   marcumQmethod <- match.arg(marcumQmethod)
   if(! check.fs(f)) return()
   if(paracheck == TRUE & ! are.parkmu.valid(para)) return()

   LOWER <- ifelse(is.na(qualo), .Machine$double.eps, qualo)

   if(is.na(quahi)) {
      order <- .Machine$double.exponent
      first <- TRUE
      while(order > -order) {
         tmp.upper <- 10^order
         if(verbose) message("Order=", order, ", upper=",round(tmp.upper, digits=4),
                             " returns CDF=", appendLF=FALSE)
         F <- cdfkmu(tmp.upper, para, getmed=FALSE,
                     marcumQ=marcumQ, marcumQmethod=marcumQmethod)
         order <- order - 0.5;
         if(verbose) message(F)
         #if(! first    & is.na(F)) break
         if(! is.na(F) &   F != 1 & ! first) break
         if(! is.na(F) &   F == 1)    first <- FALSE
      }
      if(order == -.Machine$double.exponent) {
         warning("Could not sweep down to an enhanced estimate of upper, resetting")
         UPPER <- 10^.Machine$double.exponent
      } else {
         order <- order + 0.5
         UPPER <- 10^order
      }
   } else {
      UPPER <- quahi
   }
   if(verbose) message("   'LOWER'=",LOWER)
   if(verbose) message("   'UPPER'=",UPPER)

   "afunc" <- function(X, Fx=NA) {
         theF <- cdfkmu(X, para, getmed=FALSE, marcumQ=marcumQ,
                                               marcumQmethod=marcumQmethod)
         if(is.na(theF)) theF <- 1
         #message("theX=", X,"  theF=",theF)
         #message("Fx=", Fx)
         return(Fx - theF)
   }

   x <- sapply(1:length(f), function(i) {
                 Fx <- f[i]
                 if(Fx < LOWER) return(0)
                 if(Fx == Inf) # is this ok?
                 tmpkmurt <- NULL # It seems(?) that uniroot for years
                 # if erroring within a try() would leave tmpkmurt as is but
                 # testing in Feb 2014 (R:3.1.2) shows that an rm(tmpkmurt) seems
                 # to occur!  So have to ask the search path if the variable is present
                 # because is.null() will not work, so intercept and then ifelse.
                 # Pick an arbitrary means of populating the last error message
                 try(besselK(x=-1), silent=TRUE) # get to reset the 'last error message'
                 try( tmpkmurt <- uniroot(afunc, c(LOWER,UPPER), Fx=Fx)$root, silent=TRUE)
                 err <- grep("opposite sign", geterrmessage()) # does the error have this
                 if(length(err) != 0) { # this is a nasty unnatural hack to trap an error
                    warning("non-opposite sign trap for F=",Fx,", returning NA")
                    return(NA)
                 } # see the commented out code for a "micro" version of this
                 # yet, it seems to work propertly, as WHA has long expected and
                 # designed logic around
                 ifelse(is.null(tmpkmurt), return(NA), return(tmpkmurt)) })
   return(x)
}

#"afunc" <- function(x) return(abs(x))
#results <- sapply(1:2, function(i) {
#   rtmp <- NULL
#   try(rtmp <- uniroot(afunc, c(1,2))$root, silent=FALSE)
#   print(rtmp) })
