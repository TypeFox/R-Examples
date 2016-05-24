# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function: 
# asc2s() and asc3s() functions calculates the boundary coordinates 
# for the anisotropic set cover on a 2D & 3D square lattice 
# with a fixed edge & face along the lattice boundary. 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# k, x - maximal set cover size and vector of lattice sizes;
# dir - index of the variable component: x) dir=1; y) dir=2; z) dir=3; 
# r -  variable size of set cover elements.
# Value: 
# - list of boundary coordinates and sizes 
#   for the anisotropic set cover on a 2D & 3D square lattice 
#   with a fixed edge & face along the lattice boundary.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
asc2s <- function(k=12, x=rep(95, times=2), dir=2,
                  r=(x[dir]-3)^(seq(k)/k)) {
  k <- length(r <- unique(round(r)))
  return(switch(dir, 
                rbind(r=r+1, 
                      x1=2, x2=r+2, 
                      y1=2, y2=x[2]-1),
                rbind(r=r+1,
                      x1=2, x2=x[1]-1,
                      y1=2, y2=r+2)))
}
asc3s <- function(k=12, x=rep(95, times=3), dir=3,
                  r=(x[dir]-3)^(seq(k)/k)) {
  k <- length(r <- unique(round(r)))
  return(switch(dir, 
                rbind(r=r+1, 
                      x1=2, x2=r+2, 
                      y1=2, y2=x[2]-1,
                      z1=2, z2=x[3]-1),
                rbind(r=r+1,
                      x1=2, x2=x[1]-1,
                      y1=2, y2=r+2,
                      z1=2, z2=x[3]-1),
                rbind(r=r+1,
                      x1=2, x2=x[1]-1,
                      y1=2, y2=x[2]-1,
                      z1=2, z2=r+2)))
}