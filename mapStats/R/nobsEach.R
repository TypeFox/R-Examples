  nobsEach <- function(divs) {
      #number of observations in each class
      nobs <- table( findInterval(x=divs$var, vec=divs$brks, all.inside=TRUE))
      ngroups <- length(divs$brks) -1
      #if empty groups
      if (ngroups > length(nobs)) {
         empty_groups <- (1:ngroups)[ ! (1:ngroups) %in% as.numeric(names(nobs)) ]
         te <- as.table(rep(0, length(empty_groups)))
         names(te) <- empty_groups
         #merge
         nobs <- c(nobs, te)
	 #reorder
         nobs <- nobs[ order(as.numeric(names(nobs))) ]
        }
     as.vector(nobs)
   }
       