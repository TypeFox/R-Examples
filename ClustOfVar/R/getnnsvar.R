getnnsvar <-
function(diss, flag) {
# Inputs.  diss: full distance matrix.
#          flag: "live" rows indicated by 1 are to be processed.
# Returns. List of: nn, nndiss.
#          nn:   list of nearest neighbor of each row.
#          nndiss: nearest neigbbor distance of each row.

   nn <- rep(0, nrow(diss))
   nndiss <- rep(0.0, nrow(diss))
   MAXVAL <- 1.0e12
   if (nrow(diss) != ncol(diss)) stop("Invalid input first parameter.")
   if (nrow(diss) != length(flag)) stop("Invalid inputs 1st/2nd parameters.")
 # if (nrow(diss) != length(nn)) stop("Invalid inputs 1st/3rd parameters.")
 # if (nrow(diss) != length(nndiss)) stop("Invalid inputs 1st/4th parameters.")

   for (i1 in 1:nrow(diss)) {
       if (flag[i1] == 1) {
          minobs <- -1
          mindis <- MAXVAL
          for (i2 in 1:ncol(diss)) {
              if ( (diss[i1,i2] < mindis) && (i1 != i2) ) {
                 mindis <- diss[i1,i2]
                 minobs <- i2
              }
          }
          nn[i1] <- minobs
          nndiss[i1] <- mindis
       }
   }
   list(nn = nn, nndiss = nndiss)
}

