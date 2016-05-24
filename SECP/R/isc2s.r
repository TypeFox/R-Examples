# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# isc2s() and isc3s() functions calculates the boundary coordinates 
# for the isotropic set cover on the 2D & 3D square lattice 
# with a fixed point in the lattice center. 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# k, x - maximal set cover size and vector of lattice sizes;
# o, r - set cover center and radius.
# Value: 
# - list of boundary coordinates and sizes 
#   for the isotropic set cover on a 2D & 3D square lattice 
#   with a fixed point in the lattice center.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
isc2s <- function(k=12, x=rep(95, times=2), 
                  o=(x+1)/2, r=min(o-2)^(seq(k)/k)) {
  k <- length(r <- unique(round(r)))
  return(rbind(r=2*r+1, 
               x1=o[1]-r, x2=o[1]+r, 
               y1=o[2]-r, y2=o[2]+r))
}
isc3s <- function(k=12, x=rep(95, times=3), 
                  o=(x+1)/2, r=min(o-2)^(seq(k)/k)) {
  k <- length(r <- unique(round(r)))
  return(rbind(r=2*r+1,
               x1=o[1]-r, x2=o[1]+r, 
               y1=o[2]-r, y2=o[2]+r, 
               z1=o[3]-r, z2=o[3]+r))
}