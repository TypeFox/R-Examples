# ++++++++++++++++++++++++++++++++++++++++++++
# Generate the indexes of the neighbors polygones on a regular grid
# of polygones: the neighbors of J are the polygones
# located in the square composed by "nvois" polygones around J.
#NOTE:  The polygones are numbered from the left to the right,
# from bottom to top
# ++++++++++++++++++++++++++++++++++++++++++++
generVois <- function (np=10, nvois=1) {
  npoly <- np*np
  retour <- matrix(NA, ncol=2)
  for (j in 1:npoly) {
    # The extreme polygons of the neighborhood at the bottom,
    # on the left and on the right
    # Verify that these polys are in the frame
    re <- j%%np
    # number of free compartments on the left and on the right of j
    if (re==0)
      {
        nlibresg <- np-1
        nlibresd <- 0
      }    else {
        nlibresg <- re-1
        nlibresd <- np-re
      }
    nvg <- min(nvois, nlibresg)
    nvd <- min(nvois, nlibresd)
    # The extreme polygons of the neighborhood at the bottom,
    # on the left and on the right, are then:
     bas <- c(j-(np*nvois)-nvg, j-(np*nvois) + nvd)

  # Generate all  the indexes of the squares from those:
    i<-1
    while( i <= ((nvois*2)+1)) {
      if (all(bas>npoly)) break
      while ( all(bas<=0)) {
        bas <- bas+np
        i<-i+1
        next
      }
        
      for (iv in bas[1]:bas[2])
        retour <- rbind(retour, c(j, iv))
      bas <- bas+np
    i<-i+1
    } # fin i
  } # fin j

    # Remove the polys outside the frame
      retour <- retour[retour[,2]>0,, drop=FALSE]
      retour <- retour[retour[,2]<=npoly, , drop=FALSE]
  return(retour[2:nrow(retour),,drop=FALSE])
} #end generVois
  
      
      
                                    
