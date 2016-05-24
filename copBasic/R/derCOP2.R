"derCOP2" <-
function(cop=NULL, u, v,
         delv=.Machine$double.eps^0.50,
         derdir=c("left", "right", "center"), ...) {

    derdir <- match.arg(derdir)
    if(length(v) == 1) {
       if(v - delv < 0) derdir <- "left"
       if(v + delv > 1) derdir <- "right"
       #str(cop)
       if(derdir == "left") {
         return((cop(u,v+delv,...) - cop(u,v, ...))/delv)
       } else if(derdir == "right") {
         return((cop(u,v,...)      - cop(u,v-delv, ...))/delv)
       } else {
         return((cop(u,v+delv,...) - cop(u,v-delv, ...))/(2*delv))
       }
    } else {
       if(length(v) != length(u)) {
          #warning("length of v and u are not equal, so using only first element of u")
          u <- rep(u[1], length(v))
       }
       der <- vector(mode="numeric", length(u))
       for(i in 1:length(v)) {
          tmpdir <- derdir; au <- u[i]; av <- v[i]
          if(av - delv < 0) tmpdir <- "left"
          if(av + delv > 1) tmpdir <- "right"
          if(tmpdir == "left") {
             der[i] <- (cop(au,av+delv,...) - cop(au,av, ...))/delv
          } else if(tmpdir == "right") {
             der[i] <- (cop(au,av,...)      - cop(au,av-delv, ...))/delv
          } else {
             der[i] <- (cop(au,av+delv,...) - cop(au,av-delv, ...))/(2*delv)
          }
       }
       return(der)
    }
}

