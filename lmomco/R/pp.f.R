"pp.f" <-
function(f, x) {
   if(f > 1 | f < 0) {
       warning("Error with probability argument not in [0,1]");
       return(NA);
   }
   r <- rank(x, ties.method="first");
   return(qbeta(f,r,length(r)+1-r));
}

"pp.median" <-
function(x) {
   return(pp.f(0.5,x));
}

