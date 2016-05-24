"theoLmoms.max.ostat" <-
function(para=NULL, cdf=NULL, pdf=NULL, qua=NULL, nmom=4,
         switch2minostat=FALSE, showterms=FALSE, ...) {

   if(is.null(para)) stop("parameter list not specified and it needs to be in lmomco style")

   if(is.null(qua)) {
      if(is.null(cdf)) stop("cdf function using lmomco parameter style not specified")
      if(is.null(pdf)) stop("pdf function using lmomco parameter style not specified")
   }

   enn <- vector(mode="numeric", length=nmom)
   lms <- lmr <- enn
   for(r in 1:nmom) {
     mo <- NA
     if(is.null(qua)) {
        mo <- ifelse(switch2minostat,
                     expect.min.ostat(r, para=para, cdf=cdf, pdf=pdf, ...),
                     expect.max.ostat(r, para=para, cdf=cdf, pdf=pdf, ...))
     } else {
        mo <- ifelse(switch2minostat,
                     expect.min.ostat(r, para=para, qua=qua, ...),
                     expect.max.ostat(r, para=para, qua=qua, ...))
     }
     enn[r] <- mo
     if(is.na(mo)) next
     series <- 0
     the.terms <- vector(mode="numeric", length=r)
     for(k in r:1) {
       term <- (-1)^(r-k) * k^(-1) * exp( lchoose(r-1,k-1) + lchoose(r+k-2,k-1) )
       if(switch2minostat) term <- (-1)^(r-1) * term
       the.terms[k] <- term
       series[k] <- term*enn[k]
     }
     if(showterms) {
        message("  Order (r)=", r, ": terms (m=r:1) on E[X_m:(m|1)]s --> ",
                                paste(rev(the.terms), collapse="  "))
        #print(rev(series)) # only for debugging
     }
     # Note that the (-1)^(r-1) term in typeset mathematics would be pulled out from
     # the summation, but it is inserted here so that the multipliers on the Exx have
     # the correct sign in showterms is TRUE
     lms[r] <- sum(series)
   }
   if(nmom >= 1) lmr[1] <- NA
   if(nmom >= 2) {
      lmr[2] <- lms[2]/lms[1]
      lmr[! is.finite(lmr[2])]  <- NA
      if(is.nan(lmr[2])) lmr[2] <- NA
   }
   if(nmom >= 3) {
      lmr[3:nmom] <- lms[3:nmom]/lms[2]
      lmr[! is.finite(lmr[3:nmom])] <- NA
      lmr[is.nan(lmr[3:nmom])] <- NA
   }
   src <- ifelse(switch2minostat, "theoLmoms.min.ostat",
                                  "theoLmoms.max.ostat")
   z <- list(lambdas=lms,
             ratios=lmr,
             trim=NULL,
             lefttrim=NULL,
             righttrim=NULL,
             source=src)
   return(z)
}

"theoLmoms.min.ostat" <-
function(...) {
   z <- theoLmoms.max.ostat(switch2minostat=TRUE, ...)
   return(z)
}

