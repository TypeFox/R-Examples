MandelIterate <-
function(z_0) {
#
# Return number of iterations before
# the value escapes to infinity
#   
#
   i <- 1 
   z <- 0
   ins <- 1 
   while(i < 51) {
     i <- i+1
     z <- z*z + z_0
     if(abs(z) > 2 ) { ins <- i; break; }
    } 
 return(ins);
}

