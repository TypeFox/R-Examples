"fliplmoms" <-
function(lmom, flip=NULL, checklmom=TRUE) {
    if(length(lmom$lambdas) == 0) {
    	warning("need vectorized L-moments, see lmorph()?");
    	return();
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return();
    }
    if(length(lmom$flip) == 0) {
      if(is.null(flip)) {
        warning("The flip is not provided in the lmom list or as argument");
        return();
      } else {
        lmom$flip <- flip
      }
    }
   odd <- seq(3,length(lmom$lambda), by=2);
   lmom$lambdas[odd] <- -1*lmom$lambdas[odd];
   lmom$ratios[odd]  <- -1*lmom$ratios[odd];
   lmom$lambdas[1]   <- lmom$flip    - lmom$lambdas[1];
   lmom$ratios[2]    <- lmom$lambda[2]/lmom$lambdas[1];
   return(lmom);
}
