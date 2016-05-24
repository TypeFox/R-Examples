polarization <-
function (V, old=FALSE) {
   # Calculates polarization (SOM Project)
   # Arguments: as Agreement function (V = frequnecy vector)
   p <- (1 - agreement(V,old))/2  # pass arguments to agreement function, rescale
   return(p)
   }
