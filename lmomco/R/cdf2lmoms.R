"cdf2lmoms" <-
function(para, nmom=6, fdepth=0, silent=TRUE, lambegr=1, ...) {
   if(lambegr > nmom) {
       warning("lmambegr can not be greater than number of L-moments requested by nmom")
       lambegr <- nmom
   }
   lams <- rats <- vector(mode="numeric", length=nmom)
   lams[1:nmom] <- rats[1:nmom] <- NA
   for(r in lambegr:nmom) lams[r] <- cdf2lmom(r, para, fdepth=fdepth, silent=silent, ...)
   rats[1] <- NA; rats[2] <- lams[2]/lams[1]
   ratbegr <- 3
   if(lambegr > ratbegr) ratbegr <- lambegr
   rats[ratbegr:nmom] <- sapply(ratbegr:nmom, function(r) { # notice the indexing!!
                          ifelse(is.na(lams[2]), return(NA), return(lams[r]/lams[2])) })
   z <- list(lambdas=lams,
             ratios=rats,
             trim=NULL, leftrim=NULL, rightrim=NULL,
             source="cdf2lmoms")
   return(z)
}
