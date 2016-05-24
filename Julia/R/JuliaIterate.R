JuliaIterate <-
function(z,C) {
#
# Return number of iterations before
# the value escapes to infinity
#
   i <- 1 
   ins <- 1 
   while(i < 51) {
     i <- i+1
     z <- z*z + C
     if(abs(z) > 2 ) { ins <- i; break; }
    } 
 return(ins);
}

