#############################################################################
## quaternion2rotation()
#############################################################################
#' @title Convert Quaternion into a Rotation Matrix
#' 
#' @description The affine/rotation matrix \eqn{R} is calculated from the quaternion
#' parameters.
#' 
#' @details The quaternion representation is chosen for its compactness in representing
#' rotations.  The orientation of the \eqn{(x,y,z)} axes relative to the
#' \eqn{(i,j,k)} axes in 3D space is specified using a unit quaternion
#' \eqn{[a,b,c,d]}, where \eqn{a^2+b^2+c^2+d^2=1}{a*a+b*b+c*c+d*d=1}.  The
#' \eqn{(b,c,d)} values are all that is needed, since we require that
#' \eqn{a=[1-(b^2+c^2+d^2)]^{1/2}}{a=sqrt(1.0-(b*b+c*c+d*d))} be non-negative.
#' The \eqn{(b,c,d)} values are stored in the (\code{quatern_b},
#' \code{quatern_c}, \code{quatern_d}) fields.
#' 
#' @aliases quaternion2rotation quaternion2mat44
#' @param nim is an object of class \code{nifti}.
#' @param tol is a very small value used to judge if a number is essentially
#' zero.
#' @param b is the quaternion \eqn{b} parameter.
#' @param c is the quaternion \eqn{c} parameter.
#' @param d is the quaternion \eqn{d} parameter.
#' @return The (proper) \eqn{3{\times}3}{3x3} rotation matrix or
#' \eqn{4{\times}4}{4x4} affine matrix.
#' @author Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @references NIfTI-1\cr \url{http://nifti.nimh.nih.gov/}
#' @examples
#' 
#' ## This R matrix is represented by quaternion [a,b,c,d] = [0,1,0,0]
#' ## (which encodes a 180 degree rotation about the x-axis).
#' (R <- quaternion2rotation(1, 0, 0))
#' @name quaternion2rotation
#' @rdname quaternion2rotation
#' @export
quaternion2rotation <- function(b, c, d, tol=1e-7) {
  ## compute a parameter from b,c,d
  a <- 1 - (b*b + c*c + d*d)
  if (a < tol) {                      # special case
    a <- 1 / sqrt(b*b + c*c +d*d)
    b <- a * b
    c <- a * c
    d <- a * d                        # normalize (b,c,d) vector
    a <- 0                            # a = 0 ==> 180 degree rotation
  } else {
    a <- sqrt(a)                     # angle = 2*arccos(a)
  } # a <- sqrt(1 - (b*b+c*c+d*d))
  R <- matrix(c(a*a+b*b-c*c-d*d, 2*b*c+2*a*d, 2*b*d-2*a*c,  # column 1
                2*b*c-2*a*d, a*a+c*c-b*b-d*d, 2*c*d+2*a*b,  # column 2
                2*b*d+2*a*c, 2*c*d-2*a*b, a*a+d*d-c*c-b*b), # column 3
              3, 3)
  return(R)
}

#############################################################################
## quaternion2mat44()
#############################################################################
#' @rdname quaternion2rotation
#' @export
quaternion2mat44 <- function(nim, tol=1e-7) {
  qb <- nim@"quatern_b"
  qc <- nim@"quatern_c"
  qd <- nim@"quatern_d"
  qx <- nim@"qoffset_x"
  qy <- nim@"qoffset_y"
  qz <- nim@"qoffset_z"
  dx <- pixdim(nim)[2]
  dy <- pixdim(nim)[3]
  dz <- pixdim(nim)[4]
  qfac <- pixdim(nim)[1]
  R <- matrix(0, nrow=4, ncol=4)
  b <- qb
  c <- qc
  d <- qd
  ## last row is always [ 0 0 0 1 ]
  R[4,1] <- R[4,2] <- R[4,3] <- 0.0
  R[4,4] <- 1.0
  ## compute a parameter from b,c,d
  a <- 1 - (b*b + c*c + d*d)
  if (a < tol) {                      # special case
    a <- 1 / sqrt(b*b + c*c +d*d)
    b <- a * b
    c <- a * c
    d <- a * d                        # normalize (b,c,d) vector
    a <- 0                            # a = 0 ==> 180 degree rotation
  } else {
    a <- sqrt(a)                     # angle = 2*arccos(a)
  }
  ## load rotation matrix, including scaling factors for voxel sizes
  xd <- ifelse(dx > 0, dx, 1)         # make sure are positive
  yd <- ifelse(dy > 0, dy, 1)
  zd <- ifelse(dz > 0, dz, 1)
  if (qfac < 0) {
    zd <- -zd                         # left handedness?
  }
  R[1,1] <- (a*a + b*b - c*c - d*d) * xd
  R[1,2] <- 2 * (b*c - a*d) * yd
  R[1,3] <- 2 * (b*d + a*c) * zd
  R[2,1] <- 2 * (b*c + a*d) * xd
  R[2,2] <- (a*a + c*c - b*b - d*d) * yd
  R[2,3] <- 2 * (c*d - a*b) * zd
  R[3,1] <- 2 * (b*d - a*c) * xd
  R[3,2] <- 2 * (c*d + a*b) * yd
  R[3,3] <- (a*a + d*d - c*c - b*b) * zd
  ## load offsets
  R[1,4] <- qx
  R[2,4] <- qy
  R[3,4] <- qz
  return(R)
}
