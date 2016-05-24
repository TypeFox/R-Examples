"lmomsRCmark" <-
function(x, rcmark=NULL, nmom=5, flip=NA, flipfactor=1.1) {
    n <- length(x);
    if(nmom > n)
        stop("More L-moments requested by parameter 'nmom' than data points available in 'x'");
    if(length(unique(x)) == 1)
        stop("all values are equal--Lmoments can not be computed");
    if(is.null(rcmark)) rcmark <- rep(0,n);
    if(n != length(rcmark))
        stop("sample size != right-censoring marker (rcmark)");
    rcmark <- as.numeric(rcmark);

    if(! is.na(flip)) {
       if(is.logical(flip)) {
         if(flip) {
            if(flipfactor < 1) {
               warning("flipfactor < 1, setting to unity");
               flipfactor <- 1;
            }
            flip <- flipfactor*max(x);
            x <- flip - x;
         }
       } else {
         x <- as.numeric(flip) - x;
       }
    }

    
    ix <- sort(x, index.return=TRUE)$ix;
     x <- x[ix]; rcmark <- rcmark[ix];

    L <- R <- rep(NA, nmom);
    for(i in 1:nmom) {
       L[i] <- lmomRCmark(x, rcmark=rcmark, r=i, sort=FALSE);
       if(i == 2) R[i] <- L[2]/L[1];
       if(i >= 3) R[i] <- L[i]/L[2];
    }
    z <- list(lambdas=L,
              ratios=R,
              trim=0,
              lefttrim=NULL,
              rightrim=NULL,
              n=n,
              n.cen=sum(rcmark),
              flip=flip,
              source="lmomsRCmark");
    return(z);
}



