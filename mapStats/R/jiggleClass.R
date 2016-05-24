jiggleClass <- function(x) {
  
   ngroups <- length(x$brks) -1

   if (ngroups >=2) {
   
   #eliminate any empty groups:
   nobs <- nobsEach(x)
   
   if (any(nobs==0)) { 
       warning("Some class intervals are empty. Number of groups will decrease.  Consider decreasing ngroups.") 
    
       x$brks <- x$brks[ c(1, which(nobs !=0) +1 )]
       ngroups <- length(x$brks) -1
       
       nobs <- nobsEach(x) 
   }

   #jiggle the end points
   r <- range(x$var, na.rm=TRUE)
   jiggle_end <- 0.01*(abs(r[2]-r[1]))

   x$brks[1] <- x$brks[1] - jiggle_end
   x$brks[ ngroups + 1] <- x$brks[ ngroups + 1 ] + jiggle_end


   #jiggle middle points by a bit for rounding error 
   s <- sort(x$var)
   
   ind <- cumsum(nobs)[ -ngroups ]
   if (attributes(x)$intervalClosure =="right") {
     
     #increase the breaks a little bit
     
     d <- 0.001*(s[ind +1] -x$brks[2:ngroups])
     x$brks[ 2:ngroups ] <- x$brks[ 2:ngroups ] + d
    }
   else {
     s <- c(x$brks[1], s)
     d <- 0.001*(s[ind + 1] -x$brks[2:ngroups])
     x$brks[ 2:ngroups ] <- x$brks[ 2:ngroups ] - d     
    }

   }     
  
   x

}
