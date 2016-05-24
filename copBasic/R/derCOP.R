"derCOP" <-
function(cop=NULL, u, v,
         delu=.Machine$double.eps^0.50,
         derdir=c("left", "right", "center"), ...) {

    derdir <- match.arg(derdir)
    if(length(u) == 1) {
       if(u - delu < 0) derdir <- "left"
       if(u + delu > 1) derdir <- "right"
       #str(cop)
       if(derdir == "left") {
          return((cop(u+delu,v,...) - cop(u,v, ...))/delu)
       } else if(derdir == "right") {
          return((cop(u,v,...)      - cop(u-delu,v, ...))/delu)
       } else {
          return((cop(u+delu,v,...) - cop(u-delu,v, ...))/(2*delu))
       }
    } else {
      if(length(u) != length(v)) {
        #warning("length of u and v are not equal, so using only first element of v")
        v <- rep(v[1], length(u))
      }
      der <- vector(mode="numeric", length(u))
      for(i in 1:length(u)) {
         tmpdir <- derdir; au <- u[i]; av <- v[i]
         if(au - delu < 0) tmpdir <- "left"
         if(au + delu > 1) tmpdir <- "right"
         if(tmpdir == "left") {
           der[i] <- (cop(au+delu,av,...) - cop(au,av, ...))/delu
         } else if(tmpdir == "right") {
           der[i] <- (cop(au,av,...)      - cop(au-delu,av, ...))/delu
         } else {
           der[i] <- (cop(au+delu,av,...) - cop(au-delu,av, ...))/(2*delu)
         }
      }
      return(der)
    }
}


