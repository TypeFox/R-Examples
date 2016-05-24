## 
## Functions for the creation of special clusters
## 


##' equal_sizes
##'
##' generate a matrix of equal particle sizes
##' @title equal_sizes
##' @param a a
##' @param b b
##' @param c c
##' @param N N
##' @return matrix Nx3
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
equal_sizes <- function(a, b, c, N)
  data.frame(a = rep(a, N), b = rep(b, N), c = rep(c, N))


##' equal_angles
##'
##' generate a matrix of equal angles
##' @title equal_angles
##' @param phi phi
##' @param theta theta
##' @param psi psi
##' @param N N
##' @return matrix Nx3
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
equal_angles <- function(phi=0, theta=0, psi=0, N)
cbind(phi=rep(phi, N), theta=rep(theta, N), psi=rep(psi, N))

##' helix curve
##'
##' add particles on an helix
##' @title helix
##' @param R0 radius in nm
##' @param pitch pitch in nm
##' @param N number of particles
##' @param delta twist angle between two particles
##' @param delta0 angle shift at origin
##' @param n.smooth  number of interpolation points (for plotting purposes)
##' @param right logical, handedness
##' @return list with r,  sizes,  invalpha,  angles, R0 and smooth interp. points
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
##' @examples 
##' cl <- helix(500, 1000, 36, delta=pi/6, n.smooth=1e3) ; str(cl)
##' \dontrun{
##' require(rgl)
##' open3d()
##' spheres3d(cl$smooth, radius=1,col=2)
##' ## ellipsoids are oriented following the helix
##' sizes <- equal_sizes(40, 20, 20,NROW(cl$positions))
##' rgl.ellipsoids(cl$positions, sizes, cl$angles, col="gold") 
##' }

helix <- function(R0=500, pitch=600, N=5, 
                  delta=pi/8, delta0=pi/2, n.smooth=100*N, right=TRUE){

  handedness <- (-1)^(!right)
  phase <- seq(from=delta0, by=delta, length=N)
  phase2 <- seq(from=delta0, max(phase), length=n.smooth)
  
  xyz <- function(ph){
    x = R0*cos(ph)
    y = R0*sin(ph)
    z <- handedness * ph * pitch / (2*pi)
    
    data.frame(x=x, y=y, z=z)
}
  positions <- xyz(phase)
  centering <- mean(positions$z)
  positions$z <- positions$z - centering
  positions2 <- xyz(phase2)
  positions2$z <- positions2$z - centering
  xp <-  - positions$y 
  yp <-   positions$x 
  zp <-  handedness * pitch / (2*pi)
  n <- sqrt(xp^2+yp^2+zp^2) 
  
  phi <-  atan2(yp, xp)
  psi <-  asin(zp/n)

  list(positions=as.matrix(positions),
       angles = cbind(phi, pi/2, psi),
       radius=R0, smooth=as.matrix(positions2))
}

##' cluster_dimer
##'
##' cluster with two nanorods
##' first rod along x at (0, 0, -d/2)
##' second rod at (0, 0, d/2)
##' @title cluster_dimer
##' @param d center-to-center distance
##' @param dihedral dihedral angle
##' @param alpha1 angle first rod
##' @param alpha2 angle second rod
##' @param a semi axis
##' @param b semi axis
##' @return list with r,  sizes,  angles
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
cluster_dimer <- function(d=a, 
                             dihedral=0, alpha1=0, alpha2=0,
                             a=35e-3, b=12e-3){
 
  r <- cbind(c(0,0), c(0, 0), c(-d/2, d/2))
  sizes <- equal_sizes(a=a, b=b, c=b, N=2)  
  angles <- cbind(c(dihedral, 0), c(pi/2, pi/2), c(alpha1, alpha2))
  list(r=r, sizes=sizes, angles=angles)
  
}

##' cluster_dimer_end
##'
##' cluster with two nanorods
##' first rod along x at (0, 0, -d/2)
##' second rod at (0, 0, d/2)
##' @title cluster_dimer_end
##' @param d end-to-end distance
##' @param dihedral dihedral angle
##' @param a semi axis
##' @param b semi axis
##' @param rescale logical, rescale the z coordinates so that d is the center-to-center distance
##' @return list with r,  sizes,  angles
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
cluster_dimer_end <- function(d=a, 
                          dihedral=0, 
                          a=35e-3, b=12e-3, rescale=TRUE){
   if(rescale)
     d <- sqrt(d^2 - 2*a^2*(1 - cos(dihedral)))
  
  r <- rbind(c(a,0,0), 
             c(a*cos(dihedral),a*sin(dihedral),d))
  sizes <- equal_sizes(a=a, b=b, c=b, N=2)  
  angles <- rbind(c(0, pi/2, 0), 
                  c(dihedral, pi/2, 0))
  list(r=r, sizes=sizes, angles=angles)
  
}

##' cluster_chain
##'
##' linear chain of parallel nanorods
##' @title cluster_chain
##' @param N number of rods
##' @param pitch pitch
##' @param a semi axis
##' @param b semi axis
##' @param c semi axis
##' @return list with r,  sizes, angles
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
cluster_chain <- function(N, pitch=500, a=50, b=30,c=b){

  r <- as.matrix(expand.grid(x=0, y=seq(1,N) * pitch, z=0))
  sizes <- equal_sizes(a=a, b=b, c=c, N=N)
  angles <- equal_angles(N=N)
  list(r=r, sizes=sizes, angles=angles)
  
}

##' cluster_helix
##'
##' helical cluster of ellipsoids
##' @title cluster_helix
##' @param N number of particles
##' @param R0 radius of helix
##' @param pitch pitch of helix
##' @param delta angle between particles
##' @param delta0 initial angle
##' @param right logical, helicity
##' @param a ellipsoid semi-axis 
##' @param b ellipsoid semi-axis 
##' @param c ellipsoid semi-axis 
##' @param angles type of angular orientation
##' @param seed random seed for reproducibility
##' @param ... extra arguments (ignored)
##' @return list 
##' @author baptiste Auguie
##' @export
##' @family user_level cluster
cluster_helix <- function(N=5, R0=12, pitch=15, 
                             delta=pi/2, delta0=0, right=TRUE,
                             a=5, b=a, c=a,
                             angles=c("helix", "random", "fixed"),
                             seed=123, ...){

 hel <- helix(R0=R0, pitch=pitch, N=N, delta=delta, delta0=delta0, right=right)
 nr <- NROW(hel$angles)
 r <- hel$positions
 set.seed(seed) # always have same first particles
 angles <- switch(match.arg(angles),
                  "helix" = hel$angles,
                  "fixed" = cbind(rep(0, nr),
                    rep(0, nr),
                    rep(0, nr)),
                  "random" = matrix(cbind(runif(nr, 0, 2*pi),
                    runif(nr, 0, 2*pi),
                    runif(nr, 0, 2*pi)), ncol=3, byrow=T))
 
 sizes <- equal_sizes(a, b, c, N)
 list(r=r, sizes=sizes, angles=angles)
}

