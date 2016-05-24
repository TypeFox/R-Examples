return.ibd.val <-
function(col.dat){
   #expects a vector of zeros with at most one 1
   #returns the index of where the 1 is, or zero if 1
   # is not present in col.dat    
return(which.max(col.dat)*max(col.dat))
}

