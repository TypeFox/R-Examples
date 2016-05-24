distance <-
function( locs1, locs2, geodesic=FALSE ) {##### treats 1st column as lat
    #dyn.load("~/Files/Creations/C/distance.so")
n1 <- nrow(locs1)
n2 <- nrow(locs2)
#D.out <- matrix(0, n1, n2)
d.out <- rep(0, n1 * n2)
if( geodesic ) { 
D.Mx <- .C("distance_geodesic_AB", as.double(locs1[ ,1]*pi/180), as.double(locs1[ ,2]*pi/180),  
   as.double(locs2[ ,1]*pi/180), as.double(locs2[ ,2]*pi/180), as.double(d.out), as.integer(n1), as.integer(n2) )[[5]]
} else { 
D.Mx <- .C("distance_AB", as.double(locs1[ ,1]), as.double(locs1[ ,2]),  
   as.double(locs2[ ,1]), as.double(locs2[ ,2]), as.double(d.out), as.integer(n1), as.integer(n2) )[[5]]
}
D.out <- matrix(D.Mx, n1, n2)
return(D.out)
}
