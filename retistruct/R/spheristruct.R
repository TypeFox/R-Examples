## Return next index in path
path.next <- function(i, g, h) {
  return(ifelse(h[i]==i, g[i], h[i]))
}

## Return sequence of indicies in path between i and j, governed by
## pointer vector p
path <- function(i, j, g, h) {
  if (i == j) {
    return(i)
  } else {
    return(c(i, path(path.next(i, g, h), j, g, h)))
  }
}

## Return sequence of indicies in path between i and j, governed by
## pointer vector p
path.length <- function(i, j, g, h, P) {
  if (any(is.na(c(i, j)))) {
    stop("i or j contains NA")
  }
  if (i == j) {
    return(0)
  } else {
    if (h[i] == i) {
      return(sqrt(sum((P[i,] - P[g[i],])^2)) + path.length(g[i], j, g, h, P))
    } else {
      return(path.length(h[i], j, g, h, P))
    }
  }
}

## order.Rset(Rset, gf, hf)
##
## It is nice to create Rset as an ordered set
order.Rset <- function(Rset, gf, hf) {
  ## To to this, join the path from the first two members of the set.
  R12 <- path(Rset[1], Rset[2], gf, hf)
  R21 <- path(Rset[2], Rset[1], gf, hf)
  Rset <- c(R12[-1], R21[-1])
  return(Rset)
}

## Convert a matrix containing on each line the indicies of the points
## forming a segment, and convert this to two sets of ordered pointers
segments2pointers <- function(S) {
  g <- c()
  j <- 1                                # Row of S
  k <- 1                                # Column of S
  while(nrow(S) > 0) {
    i <- S[j,3-k]                       # i is index of the next point
    g[S[j,k]] <- i                      # Set the pointer to i
    S <- S[-j,,drop=FALSE]              # We have used this row of S
    if (nrow(S) == 0) {
      return(g)
    }
    j <- which(S[,1] == i)            # Is i in the first column of S?
    if (length(j) > 1) {
      stop("The segment list is not valid as it contains an element more than twice.")
    }
    if (length(j)) {              # If so, set the current column to 1
      k <- 1
    } else {
      j <- which(S[,2] == i) # Otherwise, look for i in the second column
      k <- 2
      if (!length(j)) {
        stop(paste("No matching index for point", i, "in S."))
        return(NULL)
      }
    }
  }
  return(g)
}

## Convert a set of ordered pointers to a matrix containing on each
## line the indicies of the points forming a segment
pointers2segments <- function(g) {
  S1 <- which(!is.na(g))
  S2 <- g[S1]
  return(cbind(S1, S2))
}

##' This function creates merged and transformed versions (all
##' suffixed with \code{t}) of a number of existing variables, as well
##' as a matrix \code{Bt}, which maps a binary vector representation
##' of edge indicies onto a binary vector representation of the
##' indicies of the points linked by the edge.
##' @title  Merge stitched points and edges 
##' @param t A \code{StitchedOutline} object in which points that have
##' been added by stitching have been triangulated
##' @return Adds following fields to input
##' \item{\code{Pt}}{Transformed point locations}
##' \item{\code{Tt}}{Transformed triangulation}
##' \item{\code{Ct}}{Transformed connection set}
##' \item{\code{Cut}}{Transformed symmetric connection set}
##' \item{\code{Bt}}{Transformed binary vector representation
##' of edge indicies onto a binary vector representation of the
##' indicies of the points linked by the edge}
##' \item{\code{Lt}}{Transformed edge lengths}
##' \item{\code{ht}}{Transformed correspondences}
##' \item{\code{u}}{Indicies of unique points in untransformed space}
##' \item{\code{U}}{Transformed indicies of unique points in untransformed space}
##' \item{\code{Rset}}{The set of points on the rim (which has been reoorded)}
##' \item{\code{Rsett}}{Transformed set of points on rim}
##' \item{\code{i0t}}{Transformed index of the landmark}
##' \item{H}{mapping from edges onto corresponding edges}
##' \item{Ht}{Transformed mapping from edges onto corresponding edges}
##' @author David Sterratt
##' @export
mergePointsEdges <- function(t) {
  h <- t$h
  T <- t$T
  Cu <- t$Cu
  L <- t$L
  P <- t$P
  gf <- t$gf
  
  ## Form the mapping from a new set of consecutive indicies
  ## the existing indicies onto the existing indicies
  u <- unique(h)

  ## Transform the point set into the new indicies
  Pt  <- P[u,]

  ## Transform the point correspondance mapping to the new index space  
  ht <- c()
  for (i in 1:length(h)) {
    ht[i] <- which(u == h[i])
  }

  ## DOESN'T WORK
  ## Form the inverse mapping from the existing indicies to the new
  ## set of consecutive indicies
  ## uinv <- c()
  ## uinv[u] <- 1:length(u)
  ## ht <- uinv[h[u]]

  ## Transform the triangulation to the new index space
  Tt  <- matrix(ht[T], ncol=3)

  ## Tansform the forward pointer into the new indicies
  gft <- ht[gf]

  ## Determine H, the mapping from edges onto corresponding edges
  Cut <- matrix(ht[Cu], ncol=2)
  Cut <- t(apply(Cut, 1, sort))
  M <- nrow(Cut)
  H <- rep(0, M)
  for (i in 1:M) {
    if (!H[i]) {
      H[i] <- i
      for (j in i:M) {
        if (identical(Cut[i,], Cut[j,])) {
          H[j] <- i
        }
      }
    }
  }

  ## Form the mapping from a new set of consecutive edge indicies
  ## onto the existing edge indicies
  U <- unique(H)

  ## Transform the edge set into the new indicies
  Cut <- Cut[U,]

  ## Transform the edge correspondance mapping to the new index space  
  Ht <- c()
  for (i in 1:length(H)) {
    Ht[i] <- which(U == H[i])
  }

  ## Create the lengths of the merged edges by averaging
  Lt <- c()
  for (k in 1:length(U)) {
    is <- which(Ht == k)
    ## if (length(is)>1) {
    ##   print(L[is])
    ## }
    Lt[k] <- mean(L[is])
  }

  ## Transform the rim set
  Rset <- order.Rset(t$Rset, t$gf, t$hf)
  Rsett <- unique(ht[Rset])
  i0t <- ht[t$i0]

  ## Create the symmetric connection set
  Ct <- rbind(Cut, Cut[,2:1])

  ## Matrix to map line segments onto the points they link
  ## Bt <- Matrix(0, nrow(Pt), nrow(Ct), sparse=TRUE)
  Bt <- matrix(0, nrow(Pt), nrow(Ct))
  for (i in 1:nrow(Ct)) {
    Bt[Ct[i,1],i] <- 1
  }

  m <- merge(list(Pt=Pt, Tt=Tt, Ct=Ct, Cut=Cut, Bt=Bt, Lt=Lt, ht=ht, u=u, U=U,
                  Rset=Rset, Rsett=Rsett, i0t=i0t, P=P, H=H, Ht=Ht), t)
  class(m) <- class(t)
  return(m)
}

##' Stretch the mesh in the flat retina to a circular outline
##'
##' @title Stretch mesh
##' @param Cu Edge matrix
##' @param L Lengths in flat outline
##' @param i.fix Indicies of fixed points
##' @param P.fix Coordinates of fixed points
##' @return New matrix of 2D point locations
##' @author David Sterratt
stretchMesh <- function(Cu, L, i.fix, P.fix) {
  N <- max(Cu)
  M <- length(L)
  C <- matrix(0, 2*N, 2*N)
  for (i in 1:M) {
    C[2*Cu[i,1]-1:0,2*Cu[i,2]-1:0] <- diag(2) / L[i]
    C[2*Cu[i,2]-1:0,2*Cu[i,1]-1:0] <- diag(2) / L[i]
  }

  dupC <- duplicated(C)
  if (any(dupC)) {
    i <- which(dupC)
    Ci <- C[i,]
    dups <- which(apply(C, 1, function(x) {identical(x, Ci)}))
    message(paste("dups", dups))
    message(paste("Ci", Ci))
    for (d in dups) {
      message(paste("d", d, ":", which(C[d,]==1)))
    }
  }
  
  ind <- as.vector(rbind(2*i.fix-1, 2*i.fix))
  A <- C[-ind, -ind]
  B <- C[-ind,  ind]
  P <- matrix(t(P.fix), ncol(B), 1)
  D <- diag(apply(cbind(A, 2*B), 1, sum))
  if (is.infinite(det(D))) stop ("det(D) is infinite")
  Q <- 2 * solve(D - A) %*% B %*% P
  Q <- matrix(Q, nrow(Q)/2, 2, byrow=TRUE)
  R <- matrix(0, nrow(Q) + length(i.fix), 2)
  R[i.fix,] <- P.fix
  R[-i.fix,] <- Q
  return(R)
}

##' This takes the mesh points from the flat outline and maps them to
##' the curtailed sphere. It uses the area of the flat outline and
##' \code{phi0} to determine the radius \code{R} of the sphere. It
##' tries to get a good first approximation by using the function
##' \code{\link{stretchMesh}}.
##'
##' @title Project mesh points in the flat outline onto a sphere
##' @param r \code{Outline} object to which the following information
##' has been added with \code{\link{mergePointsEdges}}:
##' \describe{
##' \item{\code{Pt}}{The mesh point coordinates.}
##' \item{\code{Rsett}}{The set of points on the rim.}
##' \item{\code{A.tot}}{The area of the flat outline.}}
##' @return \code{reconstructedOutline} object containing the
##' following extra information
##' \item{\code{phi}}{Lattitude of mesh points.}
##' \item{\code{lmabda}}{Longitude of mesh points.}
##' \item{\code{R}}{Radius of sphere.}
##' @author David Sterratt
##' @export
projectToSphere <- function(r) {
  Pt <- r$Pt
  Rsett <- r$Rsett
  i0t <- r$i0t
  A.tot <- r$A.tot
  Cut <- r$Cut
  Lt <- r$Lt
  phi0 <- r$phi0
  lambda0 <- r$lambda0
    
  Nt <- nrow(Pt)
  Nphi <- Nt - length(Rsett)

  ## From this we can infer what the radius should be from the formula
  ## for the area of a sphere which is cut off at a lattitude of phi0
  ## area = 2 * PI * R^2 * (sin(phi0)+1)
  R <- sqrt(A.tot/(2*pi*(sin(phi0)+1)))

  ## Find lengths between successive points on rim
  C <- matrix(NA, nrow(Pt), nrow(Pt))
  for (i in 1:nrow(Cut)) {
    C[Cut[i,1],Cut[i,2]] <- Lt[i]
    C[Cut[i,2],Cut[i,1]] <- Lt[i]
  }
  L.Rsett <- rep(NA, length(Rsett))
  for (i in 1:length(Rsett)) {
    L.Rsett[i] <- C[Rsett[i],Rsett[mod1(i+1, length(Rsett))]]
  }
  ## Check that this length matches the length computed from the AnnotatedOutline
  ## FIXME - this doesn't work for one retina - need to check why
  ## if (sum(L.Rsett) != getFlatRimLength(r)) {
  ##  stop("Internal error: Mismatch in rim lengths")
  ## }
  ## Stretch mesh points to circle
  Ps <- stretchMesh(Cut, Lt, Rsett, circle(L=L.Rsett))
  x <- Ps[,1]
  y <- Ps[,2]
  phi <- -pi/2 + sqrt(x^2 + y^2)*(phi0+pi/2)
  phi[Rsett] <- phi0
  lambda <- atan2(y, x)
  lambda <- lambda - lambda[i0t] + lambda0

  p <- merge(list(phi=phi, lambda=lambda, R=R,
                  phi0=phi0, lambda0=lambda0, Ps=Ps),
             r)
  class(p) <- addClass("reconstructedOutline", r)
  return(p)
}

##
## Energy/error functions
## 

## Calculate lengths of connections on sphere
compute.lengths <- function(phi, lambda, Cu, R) {
  ## Use the upper triagular part of the connectivity matrix Cu
  phi1    <- phi[Cu[,1]]
  lambda1 <- lambda[Cu[,1]]
  phi2    <- phi[Cu[,2]]
  lambda2 <- lambda[Cu[,2]]
  l <- R*central.angle(phi1, lambda1, phi2, lambda2)

  return(l)
}

## Calculate lengths of connections on sphere
compute.areas <- function(phi, lambda, T, R) {
  P <- R * cbind(cos(phi)*cos(lambda),
                 cos(phi)*sin(lambda),
                 sin(phi))

  ## Find areas of all triangles
  areas <- -0.5/R * dot(P[T[,1],], extprod3d(P[T[,2],], P[T[,3],]), 2)
  
  return(areas)
}

##' Piecewise, smooth function that increases linearly with negative arguments. 
##' \deqn{    f(x) = \left\{
##'        \begin{array}{ll}
##'          -(x - x_0/2) & x < 0 \\
##'          \frac{1}{2x_0}(x - x_0)^2 & 0 < x <x_0 \\
##'          0 & x \ge x_0
##'          \end{array} \right.
##' }
##'
##' @title Piecewise smooth function used in area penalty
##' @param x Main argument
##' @param x0 The cutoff parameter. Above this value the function is zero.
##' @return The value of the function.
##' @author David Sterratt
f <- function(x, x0) {
  y <- x

  c1 <- x <= 0
  c2 <- (0 < x) & (x < x0)
  c3 <- x0 <= x

  y[c1] <- -(x[c1] - x0/2)
  y[c2] <- 1/2/x0*(x0 - x[c2])^2
  y[c3] <- 0

  return(y)
}

##' Derivative of \code{\link{f}}
##'
##' @title Piecewise smooth function used in area penalty
##' @param x Main argument
##' @param x0 The cutoff parameter. Above this value the function is zero.
##' @return The value of the function.
##' @author David Sterratt
fp <- function(x, x0) {
  y <- x

  c1 <- x <= 0
  c2 <- (0 < x) & (x < x0)
  c3 <- x0 <= x

  y[c1] <- -1
  y[c2] <- -1/x0*(x0 - x[c2])
  y[c3] <- 0

  return(y)
}

##' The function that computes the energy (or error) of the
##' deformation of the mesh from the flat outline to the sphere. This
##' depends on the locations of the points given in spherical
##' coordinates. The function is designed to take these as a vector
##' that is received from the \code{optim} function.
##'
##' @title The deformation energy function
##' @param p Parameter vector of \code{phi} and \code{lambda}
##' @param Cu The upper part of the connectivity matrix
##' @param C The connectivity matrix
##' @param L Length of each edge in the flattened outline
##' @param B Connectivity matrix
##' @param T Triangulation in the flattened outline
##' @param A Area of each triangle in the flattened outline
##' @param R Radius of the sphere
##' @param Rset Indicies of points on the rim
##' @param i0 Index of fixed point on rim
##' @param phi0 Lattitude at which sphere curtailed
##' @param lambda0 Longitude of fixed points
##' @param Nphi Number of free values of \code{phi}
##' @param N Number of points in sphere
##' @param alpha Area scaling coefficient
##' @param x0 Area cutoff coefficient
##' @param nu Power to which to raise area
##' @param verbose How much information to report
##' @return A single value, representing the energy of this particular
##' configuration
##' @author David Sterratt
E <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi, N,
              alpha=1, x0,  nu=1, verbose=FALSE) {
  ## Extract phis and lambdas from parameter vector
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  ## Find cartesian coordinates of points
  P <- R * cbind(cos(phi)*cos(lambda),
                 cos(phi)*sin(lambda),
                 sin(phi))

  ## Compute elastic energy
  return(Ecart(P, Cu, L, T, A, R,
               alpha, x0, nu, verbose))

}

##' The function that computes the gradient of the  energy (or error)
##' of the deformation of the mesh from the flat outline to the
##' sphere. This depends on the locations of the points given in
##' spherical coordinates. The function is designed to take these as a
##' vector that is received from the \code{optim} function.
##'
##' @title The deformation energy gradient function
##' @param p Parameter vector of \code{phi} and \code{lambda}
##' @param Cu The upper part of the connectivity matrix
##' @param C The connectivity matrix
##' @param L Length of each edge in the flattened outline
##' @param B Connectivity matrix
##' @param T Triangulation in the flattened outline
##' @param A Area of each triangle in the flattened outline
##' @param R Radius of the sphere
##' @param Rset Indicies of points on the rim
##' @param i0 Index of fixed point on rim
##' @param phi0 Lattitude at which sphere curtailed
##' @param lambda0 Longitude of fixed points
##' @param Nphi Number of free values of \code{phi}
##' @param N Number of points in sphere
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cutoff coefficient
##' @param nu Power to which to raise area
##' @param verbose How much information to report
##' @return A vector representing the derivative of the energy of this
##' particular configuration with respect to the parameter vector
##' @author David Sterratt
dE <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi, N,
               alpha=1, x0, nu=1, verbose=FALSE) {
  ## Extract phis and lambdas from parameter vector
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  cosp <- cos(phi)
  cosl <- cos(lambda)
  sinl <- sin(lambda)
  sinp <- sin(phi)
  ## Find cartesian coordinates of points
  P <- R * cbind(cosp*cosl,
                 cosp*sinl,
                 sinp)

  ## Compute force in Cartesian coordinates
  dE.dp <- -Fcart(P, C, L, T, A, R,
                  alpha, x0, nu, verbose)

  ## Convert to Spherical coordinates
  dp.dphi <- R * cbind(-sinp * cosl,
                       -sinp * sinl,
                       cosp)
  dp.dlambda <- R * cbind(-cosp * sinl,
                          cosp * cosl,
                          0)

  dE.dphi    <- rowSums(dE.dp * dp.dphi)
  dE.dlambda <- rowSums(dE.dp * dp.dlambda)

  ## Return, omitting uncessary indicies
  return(c(dE.dphi[-Rset], dE.dlambda[-i0]))
}

##' The function that computes the energy (or error) of the
##' deformation of the mesh from the flat outline to the sphere. This
##' depends on the locations of the points given in spherical
##' coordinates. The function is designed to take these as a vector
##' that is received from the \code{optim} function.
##'
##' @title The deformation energy function
##' @param P N-by-3 matrix of point coordinates
##' @param Cu The upper part of the connectivity matrix
##' @param L Length of each edge in the flattened outline
##' @param T Triangulation in the flattened outline
##' @param A Area of each triangle in the flattened outline
##' @param R Radius of sphere
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cutoff coefficient
##' @param nu Power to which to raise area
##' @param verbose How much information to report
##' @return A single value, representing the energy of this particular
##' configuration
##' @author David Sterratt
Ecart <- function(P, Cu, L, T, A, R,
                  alpha=1, x0, nu=1, verbose=FALSE) {
  ## Compute elastic energy

  ## using the upper triagular part of the
  ## connectivity matrix Cu to extract coordinates of end-points of
  ## each edge
  ## Compute lengths of edges
  ## l <- vecnorm(P2 - P1)
  l <- 2*R*asin(vecnorm(P[Cu[,2],] - P[Cu[,1],])/2/R)
  if (verbose==2) { print(l) }

  ## Compute spring energy
  E.E <- 0.5/sum(L)*sum((l - L)^2/L)
  if (verbose>=1) { print(E.E) }

  ## Compute areal penalty term if alpha is nonzero
  E.A <- 0
  if (alpha) {
    ## Find signed areas of all triangles
    a <- -0.5/R * dot(P[T[,1],], extprod3d(P[T[,2],], P[T[,3],]), 2)

    ## Now compute area energy
    E.A <- sum((A/mean(A))^nu*f(a/A, x0=x0))
    ## E.A <- sum(f(a/A, x0=x0))
  }
  return(E.E + alpha*E.A)
}

##' The function that computes the gradient of the  energy (or error)
##' of the deformation of the mesh from the flat outline to the
##' sphere. This depends on the locations of the points given in
##' spherical coordinates. The function is designed to take these as a
##' vector that is received from the \code{optim} function.
##'
##' @title The deformation energy gradient function
##' @param P N-by-3 matrix of point coordinates
##' @param C The connectivity matrix
##' @param L Length of each edge in the flattened outline
##' @param T Triangulation in the flattened outline
##' @param A Area of each triangle in the flattened outline
##' @param R Radius of sphere
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cutoff coefficient
##' @param nu Power to which to raise area
##' @param verbose How much information to report
##' @return A vector representing the derivative of the energy of this
##' particular configuration with respect to the parameter vector
##' @author David Sterratt
##' @useDynLib retistruct
Fcart <- function(P, C, L, T, A, R,
                  alpha=1, x0, nu=1, verbose=FALSE) {
  ## Compute derivative of elastic energy
  
  ## Lengths of springs
  dP <- P[C[,2],] - P[C[,1],]
  l <- 2*R*asin(vecnorm(dP)/2/R)
  if (verbose==2) { print(l) }

  ## Compute general scaling factor
  fac <- 1/sum(L)*(l - c(L, L))/c(L, L)/c(L, L) #sqrt(1-(d/2/R)^2)/d

  ## Now compute the derivatives
  ## This method is slower than using dense matrix multiplication
  ## dF <- fac * dP

  ## F.E <- matrix(0, nrow(P), 3)
  ## for (i in 1:nrow(C)) {
  ##   F.E[C[i,1],] <- F.E[C[i,1],] + dF[i,]
  ## }

  ## Now compute the derivatives
  ## F.E <- as.matrix(B %*% (fac * dP))

  F.E <- matrix(0, nrow(P), 3)
  F.E <- .Call("sum_force_components", fac*dP, C, F.E, PACKAGE="retistruct")

  ## Compute the derivative of the area component if alpha is nonzero
  if (alpha) {
    dEdpi <- matrix(0, nrow(P), 3)
    ## Here follows computation of the derivative - it's a bit
    ## complicated!
    
    ## Expand triangulation so that every point is the first point
    ## once. The number of points is effectively tripled.
    T <- rbind(T, T[,c(2,3,1)], T[,c(3,1,2)])
    A <- c(A, A, A)

    ## Compute the derivative of area with respect to the first points
    dAdPt1 <- -0.5/R * extprod3d(P[T[,2],], P[T[,3],])
    
    ## Find areas of all triangles
    a <- dot(P[T[,1],], dAdPt1, 2)
    
    ## Now convert area derivative to energy derivative
    dEdPt1 <- -(A/mean(A))^nu*fp(a/A, x0=x0)/A*dAdPt1
    ## dEdPt1 <- -fp(a/A, x0=x0)/A*dAdPt1

    ## Map contribution of first vertex of each triangle back onto the
    ## points
    ## for(m in 1:nrow(T)) {
    ##   dEdpi[T[m,1],] <- dEdpi[T[m,1],] - dEdPt1[m,]
    ## }
    dEdpi <- .Call("sum_force_components", -dEdPt1, T, dEdpi , PACKAGE="retistruct")
    
    F.E <- F.E - alpha*dEdpi
  }
  return(F.E)
}

##' Restore points to spherical manifold after an update of the
##' Lagrange integration rule
##'
##' @title Restore points to spherical manifold
##' @param P Point positions as N-by-3 matrix
##' @param R Radius of sphere
##' @param Rset Indicies of points on rim
##' @param i0 Index of fixed point
##' @param phi0 Cutoff of curtailed sphere in radians
##' @param lambda0 Longitude of fixed point on rim
##' @return Points projected back onto sphere
##' @author David Sterratt
Rcart <- function(P, R, Rset, i0, phi0, lambda0) {
  
  ## Now ensure that Lagrange constraint is obeyed

  ## Points on rim
  P[Rset,1:2] <- R*cos(phi0)*P[Rset,1:2]/vecnorm(P[Rset,1:2])
  P[Rset,3]   <- R*sin(phi0)
  
  ## All points lie on sphere
  P[-Rset,] <- R*P[-Rset,]/vecnorm(P[-Rset,])

  ## Fixed point is set
  P[i0,] <- R*c(cos(phi0)*cos(lambda0),
                cos(phi0)*sin(lambda0),
                sin(phi0))

  return(P)
}

##' In the projection of points onto the sphere, some triangles maybe
##' flipped, i.e. in the wrong orientation.  This functions determines
##' which triangles are flipped by computing the vector pointing to
##' the centre of each triangle and comparing this direction to vector
##' product of two sides of the triangle.
##'
##' @title Determine indicies of triangles that are flipped
##' @param P Points in Cartesian coordinates
##' @param Tt Triangulation of points
##' @param R Radius of sphere
##' @return List containing:
##' \item{\code{flipped}}{Indicies of in rows of \code{Tt} of flipped triangles.}
##' \item{\code{cents}}{Vectors of centres.}
##' \item{\code{areas}}{Areas of triangles.}
##' @author David Sterratt
flipped.triangles.cart <- function(P, Tt, R) {
  ## Plot any flipped triangles
  ## First find verticies and find centres and normals of the triangles
  P1 <- P[Tt[,1],]
  P2 <- P[Tt[,2],]
  P3 <- P[Tt[,3],]
  cents <- (P1 + P2 + P3)/3
  normals <- 0.5 * extprod3d(P2 - P1, P3 - P2)

  areas <- -0.5/R * dot(P1, extprod3d(P2, P3), 2)
  
  flipped <- (-dot(cents, normals, 2) < 0)
  return(list(flipped=flipped, cents=cents, areas=areas))
}

##' In the projection of points onto the sphere, some triangles maybe
##' flipped, i.e. in the wrong orientation.  This functions determines
##' which triangles are flipped by computing the vector pointing to
##' the centre of each triangle and comparing this direction to vector
##' product of two sides of the triangle.
##'
##' @title Determine indicies of triangles that are flipped
##' @param phi Vector of lattitudes of points
##' @param lambda Vector of longitudes of points
##' @param Tt Triangulation of points
##' @param R Radius of sphere
##' @return List containing:
##' \item{\code{flipped}}{Indicies of in rows of \code{Tt} of flipped triangles.}
##' \item{\code{cents}}{Vectors of centres.}
##' \item{\code{areas}}{Areas of triangles.}
##' @author David Sterratt
flipped.triangles <- function(phi, lambda, Tt, R) {
  return(flipped.triangles.cart(sphere.spherical.to.sphere.cart(phi, lambda, R), Tt, R))
}


##' Optimise the mapping from the flat outline to the sphere
##'
##' @title Optimise mapping
##' @param r reconstructedOutline object
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cut-off coefficient
##' @param nu Power to which to raise area
##' @param method Method to pass to \code{optim}
##' @param plot.3d If \code{TRUE} make a 3D plot in an RGL window
##' @param dev.flat Device handle for plotting flatplot updates to. If
##' \code{NA} don't make any flat plots
##' @param dev.polar Device handle for plotting polar plot updates
##' to. If \code{NA} don't make any polar plots.
##' @param control Control argument to pass to \code{optim}
##' @return reconstructedOutline object
##' @author David Sterratt
##' @export
optimiseMapping <- function(r, alpha=4, x0=0.5, nu=1, method="BFGS",
                             plot.3d=FALSE, dev.flat=NA, dev.polar=NA,
                             control=list()) {
  phi <- r$phi
  lambda <- r$lambda
  R <- r$R
  phi0 <- r$phi0
  lambda0 <- r$lambda0
  Tt <- r$Tt
  A <- r$A
  Cut <- r$Cut
  Ct <- r$Ct
  Pt <- r$Pt
  Lt <- r$Lt
  Bt <- r$Bt
  Rsett <- r$Rsett
  i0t <- r$i0t
  Nt <- nrow(Pt)  
  Nphi <- Nt - length(Rsett)
  
  ## Optimisation and plotting 
  opt <- list()
  opt$p <- c(phi[-Rsett], lambda[-i0t])
  opt$conv <- 1
  count <- 0
  while (opt$conv) {
    ## Optimise
    opt <- optim(opt$p, E, gr=dE,
                  method=method,
                  T=Tt, A=A, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R,
                  alpha=alpha,  N=Nt, x0=x0, nu=nu,
                  Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi,
                  verbose=FALSE, control=control)
    
    ## Report
    E.tot <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
               alpha=alpha,  N=Nt, x0=x0, nu=nu,
               Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
    E.l <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
               alpha=0,  N=Nt, x0=x0, nu=nu,
               Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)

    ft <- flipped.triangles(phi, lambda, Tt, R)
    nflip <- sum(ft$flipped)
    message(sprintf("E = %8.5f | E_L = %8.5f | E_A = %8.5f | %3d flippped triangles", E.tot, E.l, E.tot - E.l,  nflip))
    if (nflip) {
      print(data.frame(rbind(id=which(ft$flipped),
                             A=A[ft$flipped],
                             a=ft$areas[ft$flipped])))
    }
    
    ## Decode p vector
    phi          <- rep(phi0, Nt)
    phi[-Rsett]  <- opt$p[1:Nphi]
    lambda       <- rep(lambda0, Nt)
    lambda[-i0t] <- opt$p[Nphi+1:(Nt-1)]

    ## Plot
    if (plot.3d) {
      sphericalplot(list(phi=phi, lambda=lambda, R=R,
                          Tt=Tt, Rsett=Rsett, gb=r$gb, ht=r$ht),
                     datapoints=FALSE)
    }

    if (!is.na(dev.flat)) {
      dev.set(dev.flat)
      flatplot(r, grid=TRUE, strain=TRUE,
               datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE)
    }

    if (!is.na(dev.polar)) {
      dev.set(dev.polar)
      r$phi <- phi
      r$lambda <- lambda
      projection(r)
    }
  }

  o <- merge(list(phi=phi, lambda=lambda, opt=opt, nflip=sum(ft$flipped),
                  E.tot=E.tot, E.l=E.l),
             r)
  o$mean.strain    <- mean(abs(getStrains(o)$spherical$strain))
  o$mean.logstrain <- mean(abs(getStrains(o)$spherical$logstrain))
  class(o) <- class(r)
  return(o)
}

##' Optimise the mapping from the flat outline to the sphere
##'
##' @title Optimise mapping
##' @param r reconstructedOutline object
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area penalty cutoff coefficient
##' @param nu Power to which to raise area
##' @param method Method to pass to \code{optim}
##' @param plot.3d If \code{TRUE} make a 3D plot in an RGL window
##' @param dev.flat Device handle for plotting grid to
##' @param dev.polar Device handle for plotting ploar plot to
##' @param ... Extra arguments to pass to \code{\link{fire}}
##' @return reconstructedOutline object
##' @author David Sterratt
##' @export
solveMappingCart <- function(r, alpha=4, x0=0.5, nu=1, method="BFGS",
                               plot.3d=FALSE, dev.flat=NA, dev.polar=NA, ...) {
  phi <- r$phi
  lambda <- r$lambda
  R <- r$R
  phi0 <- r$phi0
  lambda0 <- r$lambda0
  Tt <- r$Tt
  A <- r$A
  Cut <- r$Cut
  Ct <- r$Ct
  Pt <- r$Pt
  Lt <- r$Lt
  Bt <- r$Bt
  Rsett <- r$Rsett
  i0t <- r$i0t
  Nt <- nrow(Pt)  
  Nphi <- Nt - length(Rsett)
  
  ## Optimisation and plotting 
  opt <- list()
  opt$x <- sphere.spherical.to.sphere.cart(phi, lambda, R)
  opt$conv <- 1

  ## Compute "mass" for each node
  minL <- rep(Inf, nrow(Pt))
  for (i in 1:nrow(Cut)) {
    minL[Cut[i,1]] <- min(minL[Cut[i,1]], Lt[i])
    minL[Cut[i,2]] <- min(minL[Cut[i,2]], Lt[i])
  }
  m <- 1/minL
  m <- m/mean(m)
  count <- 50
  
  while (opt$conv && count) {
    ## Optimise
    opt <- fire(opt$x,
                force=function(x) {Fcart(x, Ct, Lt, Tt, A, R, alpha, x0, nu)},
                restraint=function(x) {Rcart(x, R, Rsett, i0t, phi0, lambda0)},
                dt=1,
                nstep=200,
                m=m, verbose=TRUE, ...) 
    count <- count - 1
    ## Report
    E.tot <- Ecart(opt$x, Cu=Cut, L=Lt, R=R, T=Tt, A=A,
                   alpha=alpha, x0=x0, nu=nu)
    E.l <- Ecart(opt$x, Cu=Cut, L=Lt, R=R, T=Tt, A=A,
                 alpha=0, x0=x0, nu=0)

    s <- sphere.cart.to.sphere.spherical(opt$x, R)
    phi <-    s[,"phi"]
    lambda <- s[,"lambda"]
    ft <- flipped.triangles(phi, lambda, Tt, R)
    nflip <- sum(ft$flipped)
    message(sprintf("E = %8.5f | E_L = %8.5f | E_A = %8.5f | %3d flippped triangles", E.tot, E.l, E.tot - E.l,  nflip))
    if (nflip) {
      print(data.frame(rbind(id=which(ft$flipped),
                             A=A[ft$flipped],
                             a=ft$areas[ft$flipped])))
    }

    ## Plot
    if (plot.3d) {
      sphericalplot(list(phi=phi, lambda=lambda, R=R,
                          Tt=Tt, Rsett=Rsett, gb=r$gb, ht=r$ht),
                     datapoints=FALSE)
    }

    if (!is.na(dev.flat)) {
      dev.set(dev.flat)
      flatplot(r, grid=TRUE, strain=TRUE,
                datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE)
    }

    if (!is.na(dev.polar)) {
      dev.set(dev.polar)
      r$phi <- phi
      r$lambda <- lambda
      projection(r)
    }
  }

  o <- merge(list(phi=phi, lambda=lambda, opt=opt, nflip=sum(ft$flipped),
                  E.tot=E.tot, E.l=E.l),
             r)
  o$mean.strain    <- mean(abs(getStrains(o)$spherical$strain))
  o$mean.logstrain <- mean(abs(getStrains(o)$spherical$logstrain))
  class(o) <- class(r)
  return(o)
}



