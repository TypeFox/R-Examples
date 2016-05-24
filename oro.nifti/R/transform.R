##
##
## Copyright (c) 2009-2011 Brandon Whitcher and Volker Schmid
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
## $Id: transform.R 332 2010-01-29 16:54:07Z bjw34032 $
##

############################################################################
## performPermutation
############################################################################
#' @name performPermutation
#' @title Transform array with orthogonal permutation matrix
#'
#' @description Given an orthogonal permutation matrix \eqn{T}, an array of
#' dimensions and a one-dimensional representation of data.  It will return a
#' transformed array with the transformed dimensions.
#' 
#' @details This function is mainly used by the \code{\link{reorient}} function
#' to transform nifti data into neuroradiological convention.
#' 
#' @param T is an orthogonal matrix.
#' @param real.dimensions is a one-dimensional array, representing the length
#' of dimensions in data.
#' @param data is a one-dimensional representation of the data to be
#' transformed.
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @author Andrew Thornton \email{zeripath@@users.sourceforge.net}
#' @seealso \code{\link{reorient}},\code{\link{inverseReorient}}
#' @export
#' @rdname performPermutation
performPermutation <- function(T, real.dimensions, data, verbose=FALSE) {
  workingdims <- (
    function(r) {
      lr <- length(r)
      if (lr <= 5) {
        c(r, rep(1, 5 - lr))
      } else {
        stop("array has dim > 5")
      }
    }
  )(real.dimensions) # An anonymous function
  
  if (verbose) {
    cat(" ## performPermutation", fill=TRUE)
    cat("  trans =", fill=TRUE)
    print(T)
    cat("  sum(T != 0) =", sum(T != 0), fill=TRUE)
    cat("       det(T) =", det(T), fill=TRUE)
    cat("  sum(T != 0) == 3 && det(T) != 0 is",
        (sum(T != 0) == 3 && det(T) != 0), fill=TRUE)
  }
  if (sum(T != 0) == 3 && det(T) != 0) {
    ## Now ensure T is descaled and work out the permutation
    trans <- sign(T)
    perms <- abs(trans %*% 1:3)
    if (length(perms) != length(workingdims)) {
      perms <- (c(perms, (length(perms)+1):length(workingdims)))
    }
    reverselist <- c(trans %*% rep(1,3) < 0, 
                     rep(FALSE, length(workingdims) - 3))
    if (any(reverselist[2:length(reverselist)]) || 
          any(perms != 1:length(perms))) {
      if (verbose) {
        cat(" ## if(any(reverselist[2...", fill=TRUE)
      }
      if (verbose) {
        cat("need to permute", fill=TRUE)
      }
      translatedData <- array(data, workingdims)
      ## Now if we have to do a permutation or reverse more than the first
      ## index we will be going slow anyway, so...
      prs <- (
        function(reverse, dims) {
          function(x) { 
            if (reverse[x]) {
              rev(1:dims[x])
            } else {
              1:dims[x]
            }
          }
        }
      )(reverselist, workingdims) # An anonymous function
      translatedData <- translatedData[prs(1), prs(2), prs(3), prs(4), prs(5),
                                       drop=FALSE]
      out <- array(aperm(translatedData, perms), real.dimensions)
    } else {
      if (reverselist[1]) {
        if (verbose) {
          cat(" ## if(reverselist[1]) {", fill=TRUE)
        }
        ## We just need to reverse the first index.
        out <- array(array(data, workingdims)[workingdims[1]:1,,,,],
                     real.dimensions)
      } else {
        if (verbose) {
          cat(" ## else {", fill=TRUE)
          cat(" ## real.dimensions =", real.dimensions, fill=TRUE)
        }
        out <- array(data, real.dimensions)
        if (verbose) {
          cat(" ## dim(out) =", dim(out), fill=TRUE)
        }
      }
    }
  } else {
    stop("Transformation is not simple, cannot reorient!")
  }
  return(out)
}

############################################################################
## reorient
############################################################################
#' @name reorient
#' @title Reorient Image using NIfTI header
#' 
#' @description Transforms in the NIfTI header are parsed and normalized
#' versions of these transforms are applied.
#' 
#' @details This function utilizes the \code{performPermutation} function
#' internally.
#' 
#' @aliases reorient inverseReorient
#' @param nim is an object of class \code{nifti}.
#' @param data is an array associated with \code{nim}.
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @param invert stores the inverse transform.
#' @param tol is a very small value used to judge if a number is essentially
#' zero.
#' @author Andrew Thornton \email{zeripath@@users.sourceforge.net},\cr
#' Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{performPermutation}}
#' @rdname reorient
#' @export
reorient <- function(nim, data, verbose=FALSE, invert=FALSE, tol=1e-7) {
  ## From nifti1.h there are three different methods of orienting the
  ## i,j,k data into x,y,z space.

  ## Method 1 is the default for ANALYZE 7.5 files, and will be dealt
  ## with last
  ##
  ## This function will try to reorient the data into an i, j, k space
  ## where increasing (+) (i,j,k) is correlated with
  ## (LEFT,ANTERIOR,SUPERIOR)
  real.dimensions <- nim@"dim_"[2:(1+nim@"dim_"[1])]
  if (nim@"qform_code" > 0) {
    if (verbose) {
      cat(" ## reorient = Method 2", fill=TRUE)
    }
    ## This method determines [x] coordinates by the pixdim[] scales
    ## (the values of which we don't care about just the signs):
    S <- diag(pixdim(nim)[2:4])
    ## A scaling factor, (because quaternions have to have det(R)=1):
    qfac <- pixdim(nim)[1]
    ## which is either = 1 or -1
    if (abs(qfac) != 1) { 
      if (verbose) {
        cat("ScalingFactor pixdim(nim)[1] =", qfac,
            "!= -1 or 1. Defaulting to 1", fill=TRUE)
      }
      qfac <- 1 
    }
    ## and applied to S[3,3]
    S[3,3] <- S[3,3] * qfac
    ## a quaternion rotation matrix, R:
    R <- quaternion2rotation(nim@"quatern_b", nim@"quatern_c", nim@"quatern_d")
    ## And a shift, which we'll have to update if we do any reversals
    shift <- c(nim@"qoffset_x", nim@"qoffset_y", nim@"qoffset_z")
    ## X is determined by X <- R %*% S %*% (I - 1) + shift
    ## Now we can reorient the data only if the X axes are aligned with
    ## the I axes i.e. there are only 3 non-zero values in the matrix RS 
    RS <- R %*% S
    ## Now descale RS and work out the permutation
    if (verbose) {
      cat("  RS =", fill=TRUE)
      print(RS)
      cat("  sum(RS != 0) =", sum(RS != 0), fill=TRUE)
      cat("       det(RS) =", det(RS), fill=TRUE)
      cat("  sum(RS != 0) == 3 && det(RS) != 0 is",
          (sum(RS != 0) == 3 && det(RS) != 0), fill=TRUE)
    }
    trans <- sign(ifelse(abs(RS) < tol, 0, RS)) # sign(ceiling(zapsmall(RS)))
    ## We will need to do something with the trans later...
    trans[1,1] <- -1 * trans[1,1]
  } else {
    if (nim@"sform_code" > 0) {
      if (verbose) {
        cat(" # reorient = Method 3", fill=TRUE)
      }
      ## [x] is given by a general affine transformation from [i]
      ##
      ## [x] <- A %*% (i-1) + shift
      ## 
      ## where A <- nim@srow_[x,y,z][1:3] and shift <- srow_x,y,z[4]
      S <- array(dim=c(3,4))
      S[1,] <- nim@"srow_x"
      S[2,] <- nim@"srow_y"
      S[3,] <- nim@"srow_z"
      shift <- S[,4]
      A <- S[,1:3]
      if (verbose) {
        cat("  A =", fill=TRUE)
        print(A)
        cat("  sum(A != 0) =", sum(A != 0), fill=TRUE)
        cat("       det(A) =", det(A), fill=TRUE)
        cat("  sum(A != 0) == 3 && det(A) != 0 is",
            (sum(A != 0) == 3 && det(A) != 0), fill=TRUE)
      }
      trans <- sign(ifelse(abs(A) < tol, 0, A))
      trans[1,1] <- -1 * trans[1,1]
    } else {
      if (verbose) {
        cat(" # reorient = Method 1", fill=TRUE)
      }
      scaling <- diag(pixdim(nim)[2:4])
      trans <- sign(ifelse(abs(scaling) < tol, 0, scaling))
      ## Method 1 by default has +x going LEFT so no sign-change for trans[1,1]
    }
  }
  if (invert) {
    trans <- qr.solve(trans)
  }
  pp <- performPermutation(trans, real.dimensions, data, verbose)
  return(pp)
}

############################################################################
## inverseReorient
############################################################################
#' @rdname reorient
#' @export
inverseReorient <- function(nim, verbose=FALSE) {
  reorient(nim, nim@.Data, verbose=verbose, invert=TRUE)
}

############################################################################
## integerTranslation
############################################################################
#' @title integerTranslation
#' 
#' @description ...
#' 
#' @details ...
#' 
#' @aliases integerTranslation invertIntegerTranslation
#' @param nim is an object of class \code{nifti}.
#' @param data is ...
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @return ...
#' @author Andrew Thornton \email{zeripath@@users.sourceforge.net}
#' @export integerTranslation
#' @rdname integerTranslation
#' @name integerTranslation
integerTranslation <- function(nim, data, verbose=FALSE) {
  ## 3D IMAGE (VOLUME) ORIENTATION AND LOCATION IN SPACE:
  ## There are 3 different methods by which continuous coordinates can
  ## attached to voxels.  The discussion below emphasizes 3D volumes,
  ## and the continuous coordinates are referred to as (x,y,z).  The
  ## voxel index coordinates (i.e., the array indexes) are referred to
  ## as (i,j,k), with valid ranges:
  ##   i = 0 .. dim[1]-1
  ##   j = 0 .. dim[2]-1  (if dim[0] >= 2)
  ##   k = 0 .. dim[3]-1  (if dim[0] >= 3)
  ## The (x,y,z) coordinates refer to the CENTER of a voxel.  In
  ## methods 2 and 3, the (x,y,z) axes refer to a subject-based
  ## coordinate system, with
  ##   +x = Right  +y = Anterior  +z = Superior.
  ## This is a right-handed coordinate system.  However, the exact
  ## direction these axes point with respect to the subject depends on
  ## qform_code (Method 2) and sform_code (Method 3).
  dims <- 2:(1+nim@"dim_"[1])
  if (nim@"qform_code" <= 0 && nim@"sform_code" <= 0 ) {
    if (verbose) {
      cat("  dims =", nim@"dim_"[dims], fill=TRUE)
    }
    return(array(data, nim@"dim_"[dims]))
  } else {
    i <- 0:(nim@"dim_"[2]-1)
    j <- 0:(nim@"dim_"[3]-1)
    k <- 0:(nim@"dim_"[4]-1)
    ijk <- cbind(rep(i, nim@"dim_"[3] * nim@"dim_"[4]),
                 rep(rep(j, each=nim@"dim_"[2]), nim@"dim_"[4]),
                 rep(k, each=nim@"dim_"[2] * nim@"dim_"[3]))
    index.ijk <- (ijk[,1] +
                  ijk[,2] * nim@"dim_"[2] +
                  ijk[,3] * nim@"dim_"[2] * nim@"dim_"[3])
    ## check for qform codes
    if (nim@"qform_code" > 0) {
      if (verbose) {
        cat("  NIfTI-1: qform_code > 0", fill=TRUE)
      }
      qfac <- pixdim(nim)[1]
      R <- quaternion2rotation(nim@"quatern_b", nim@"quatern_c", nim@"quatern_d")
      ## HACK!!! To ensure matrix is integer-valued
      R <- ceiling(R)
      ##
      qoffset <- c(nim@"qoffset_x", nim@"qoffset_y", nim@"qoffset_z")
      if (qfac < 0) {
        R[3,3] <- -R[3,3]
      }
      if (all(abs(R) == diag(3))) {
        ## HACK!!! Multiply x-dimension for proper orientation in R
        R[1,] <- -R[1,]
        xyz <-
          t(sweep(R %*% t(sweep(ijk, 2, as.array(pixdim(nim)[2:4]), "*")),
                  1, as.array(qoffset), "+"))
        index.xyz <- (xyz[,1] +
                      xyz[,2] * nim@"dim_"[2] +
                      xyz[,3] * nim@"dim_"[2] * nim@"dim_"[3])      
        if (verbose) {
          cat("  dims =", nim@"dim_"[dims], fill=TRUE)
        }
        return(array(data[order(index.xyz)], nim@"dim_"[dims]))
      } else {
        stop("-- rotation matrix is NOT (approximately) diagonal with +/- 1s --")
      }
    }
    ## check for sform codes
    if (nim@"sform_code" > 0) {
      if (verbose) cat("  NIfTI-1: sform_code > 0", fill=TRUE)
      xyz <- matrix(0, length(data), 3)
      xyz[,1] <- (nim@"srow_x"[1] * ijk[,1] + nim@"srow_x"[2] * ijk[,2] +
                  nim@"srow_x"[3] * ijk[,3] + nim@"srow_x"[4])
      ## HACK!!! Multiply x-dimension for proper orientation in R
      xyz[,1] <- -xyz[,1]
      xyz[,2] <- (nim@"srow_y"[1] * ijk[,1] + nim@"srow_y"[2] * ijk[,2] +
                  nim@"srow_y"[3] * ijk[,3] + nim@"srow_y"[4])
      xyz[,3] <- (nim@"srow_z"[1] * ijk[,1] + nim@"srow_z"[2] * ijk[,2] +
                  nim@"srow_z"[3] * ijk[,3] + nim@"srow_z"[4])
      index.xyz <- (xyz[,1] +
                    xyz[,2] * nim@"dim_"[2] +
                    xyz[,3] * nim@"dim_"[2] * nim@"dim_"[3])
      if (verbose) {
        cat("  dims =", nim@"dim_"[dims], fill=TRUE)
      }
      return(array(data[order(index.xyz)], nim@"dim_"[dims]))
    }
  }
  invisible()
}

############################################################################
## invertIntegerTranslation
############################################################################
#' @export
#' @rdname integerTranslation
invertIntegerTranslation <- function(nim, verbose=FALSE) {
  dims <- 2:(1+nim@"dim_"[1])
  if (nim@"qform_code" <= 0 && nim@"sform_code" <= 0) {
    if (verbose) {
      cat("  dims =", nim@"dim_"[dims], fill=TRUE)
    }
    return(nim@.Data)
  } else {
    i <- 0:(nim@"dim_"[2]-1)
    j <- 0:(nim@"dim_"[3]-1)
    k <- 0:(nim@"dim_"[4]-1)
    ijk <- cbind(rep(i, nim@"dim_"[3] * nim@"dim_"[4]),
                 rep(rep(j, each=nim@"dim_"[2]), nim@"dim_"[4]),
                 rep(k, each=nim@"dim_"[2] * nim@"dim_"[3]))
    index.ijk <- (ijk[,1] +
                    ijk[,2] * nim@"dim_"[2] +
                    ijk[,3] * nim@"dim_"[2] * nim@"dim_"[3])
    ## check for qform codes
    if (nim@"qform_code" > 0) {
      if (verbose) {
        cat("  NIfTI-1: qform_code > 0", fill=TRUE)
        cat("  dims =", nim@"dim_"[dims], fill=TRUE)
      }
      qfac <- pixdim(nim)[1]
      R <- quaternion2rotation(nim@"quatern_b", nim@"quatern_c", nim@"quatern_d")
      qoffset <- c(nim@"qoffset_x", nim@"qoffset_y", nim@"qoffset_z")
      if (qfac < 0) {
        R[3,3] <- -R[3,3]
      }
      if (all(abs(R) == diag(3))) {
        ## HACK!!! Multiply x-dimension for proper orientation in R
        R[1,] <- -R[1,]
        xyz <-
          t(sweep(R %*% t(sweep(ijk, 2, as.array(pixdim(nim)[2:4]), "*")),
                  1, as.array(qoffset), "+"))
        index.xyz <- (xyz[,1] +
                        xyz[,2] * nim@"dim_"[2] +
                        xyz[,3] * nim@"dim_"[2] * nim@"dim_"[3])      
        if (verbose) {
          cat("  dims =", nim@"dim_"[dims], fill=TRUE)
        }
        return(nim@.Data[order(index.xyz)])
      } else {
        stop("-- rotation matrix is NOT diagonal with +/- 1s --")
      }
      ## stop("-- qform_code > 0 not implemented --")
    }
    ## check for sform codes
    if (nim@"sform_code" > 0) {
      if (verbose) {
        cat("  NIfTI-1: sform_code > 0", fill=TRUE)
        cat("  dims =", nim@"dim_"[dims], fill=TRUE)
      }
      xyz <- matrix(0, length(nim@.Data), 3)
      xyz[,1] <- (nim@"srow_x"[1] * ijk[,1] + nim@"srow_x"[2] * ijk[,2] +
                    nim@"srow_x"[3] * ijk[,3] + nim@"srow_x"[4])
      ## HACK!!! Multiply x-dimension for proper orientation in R
      xyz[,1] <- -xyz[,1]
      xyz[,2] <- (nim@"srow_y"[1] * ijk[,1] + nim@"srow_y"[2] * ijk[,2] +
                    nim@"srow_y"[3] * ijk[,3] + nim@"srow_y"[4])
      xyz[,3] <- (nim@"srow_z"[1] * ijk[,1] + nim@"srow_z"[2] * ijk[,2] +
                    nim@"srow_z"[3] * ijk[,3] + nim@"srow_z"[4])
      index.xyz <- (xyz[,1] +
                      xyz[,2] * nim@"dim_"[2] +
                      xyz[,3] * nim@"dim_"[2] * nim@"dim_"[3])
      return(nim@.Data[order(index.xyz)])
    }
  }
  invisible()
}

############################################################################
## translateCoordinate
############################################################################
#' @name translateCoordinate
#' @title Translate Voxel Coordinates
#' 
#' @description Translates a voxel index into the continuous coordinate space
#' defined by the NIfTI qform and sform information.
#' 
#' @details This function takes as input a \code{nifti} object and an index
#' vector in the voxel space of the object and translates that voxel index
#' into the continuous coordinate space defined by the object's qform and
#' sform.
#' 
#' Please note:
#' \enumerate{ 
#' \item By default the index \code{i} varies most rapidly, etc.  
#' \item The ANALYZE 7.5 coordinate system is \tabular{ccl}{ +x \tab = \tab 
#' Left\cr +y \tab = \tab Anterior\cr +z \tab = \tab Superior } (A 
#' left-handed co-ordinate system).
#' \item The three methods below give the locations of the voxel centres in 
#' the x,y,z system.  In many cases programs will want to display the data 
#' on other grids.  In which case the program will be required to convert 
#' the desired (x,y,z) values in to voxel values using the inverse 
#' transformation.  
#' \item Method 2 uses a factor \code{qfac} which is either -1 or 1.  
#' \code{qfac} is stored in \code{pixdim[0]}.  If \code{pixdim[0]} != 1 or 
#' -1, which should not occur, we assume 1.  
#' \item The units of the \code{xyzt} are set in \code{xyzt_units} field.  
#' }
#' @param i An index vector in \code{nim}.
#' @param nim An object of class \code{nifti}.
#' @param verbose Provide detailed output to the user.
#' @return A \code{nifti}-class object with translated coordinates.
#' @author Andrew Thornton \email{zeripath@@users.sourceforge.net}
#' @examples
#' 
#' ffd <- readNIfTI(file.path(system.file("nifti", package="oro.nifti"),
#'                            "filtered_func_data"))
#' xyz <- c(1,1,1)
#' translateCoordinate(xyz, ffd, verbose=TRUE)
#' xyz <- trunc(dim(ffd)[1:3]/2)
#' translateCoordinate(xyz, ffd, verbose=TRUE)
#' @rdname transformCoordinate
#' @export
translateCoordinate <- function(i, nim, verbose=FALSE) {
  ## 3D Image orientation and location in space (as per nifti1.h)
  ##
  ## There are three different methods by which xyzt(i,j,k) may be determined
  ## I will henceforth write x for realspace co-ord of i a voxel in the data
  ## array with indices i,j,k i.e. nim@.Data[i,j,k].
  ## 
  ## NB 1: nifti refers to voxel co-ords as:
  ## i = 0:dim[1] etc however, R indices are 1 based
  ##
  ## NB 2: x(i) refers to the centre of the voxel.
  ## 
  ## NB 3: Methods 2 & 3 have subject based co-ords with increasing x,y,z going
  ## Right, Anteriorly, Superiorly respectively
  ## 
  ## This is a right-handed coordinate system. However, the exact direction
  ## these axes point with respect to the subject depends on qform_code
  ## (Method 2) and sform_code (Method 3).
  ##
  ## More NBs:
  ## 1. By default the i index varies most rapidly, etc.
  ## 2. The ANALYZE 7.5 coordinate system is
  ##    +x = Left  +y = Anterior  +z = Superior
  ## (A left-handed co-ordinate system)
  ## 3. The three methods below give the locations of the voxel centres in the 
  ##    x,y,z system. In many cases programs will want to display the data on
  ##    other grids. In which case the program will be required to convert the
  ##    desired (x,y,z) values in to voxel values using the inverse
  ##    transformation.
  ## 4. Method 2 uses a factor qfac which is either -1 or 1. qfac is stored in
  ##    pixdim[0]. if pixdim[0]!= 1 or -1, which should not occur, we assume 1.
  ## 5. The units of the xyzt are set in xyzt_units field

  if (length(i) == 1) {
    i <- rep(i, 3)
  }
  if (! is.matrix(i)) {
    i <- as.matrix(i)
  }
  if (verbose) {
    cat("Input voxel coordinate:", fill=TRUE)
    print(i)
  }
  ## Method 2. when qform_code > 0, which should be the "normal" case
  if (nim@"qform_code" > 0) { 
    if (verbose) {
      cat("QForm_code =", nim@"qform_code", ": Orientation by Method 2.", 
          fill=TRUE)
    }
    ## The [x] coordinates are given by the pixdim[] scales:
    scaling <- diag(pixdim(nim)[2:4])
    ##  a quaternion rotation matrix, R:
    R <- quaternion2rotation(nim@"quatern_b", nim@"quatern_c", nim@"quatern_d")
    ## A scaling factor, (because quaternions have to have det(R)=1):
    qfac <- pixdim(nim)[1]
    ## and a shift 
    shift <- c(nim@"qoffset_x", nim@"qoffset_y", nim@"qoffset_z")
    ## This method is intended to represent "scanner-anatomical"
    ## coordinates, which are often embedded in the image header, e.g. DICOM
    ## fields (0020,0032), (0020,0037), (0028,0030), and (0018,0050). These
    ## represent the nominal orientation and location of the data. This method
    ## can also be used to represent "aligned" coordinates, which would
    ## typically result from post-acquisition alignment of the volume to a
    ## "standard" orientation e.g.the same subject on another day, or a rigid
    ## rotation to true anatomical orientation from the tilted position of the
    ## subject in the scanner.
    ## 
    ## [x] =  _R_%*%scaling%*%[i]+shift
    ##
    ## Where scaling[3] = scaling[3]*scalingFactor
    ##
    ## first enforce qfac = 1 or -1
    if (qfac != 1 && qfac != -1) { 
      if (verbose) {
        cat("ScalingFactor pixdim(nim)[1]) =", qfac,
            "!= -1 or 1. Defaulting to 1", fill=TRUE)
      }
      qfac <- 1 
    }
    scaling[3,3] <- scaling[3,3] * qfac
    ## Now we can reorient the data only if the X axes are aligned with the I
    ## axes i.e. there are only 3 non-zero values in the matrix RS 
    RS <- R %*% scaling
    return(RS %*% (i - 1) + shift)
  }
  ## Method 3. when sform_code > 0
  if (nim@"sform_code" > 0) {
    if (verbose) {
      cat("SForm_code =", nim@"sform_code", ": Orientation by Method 3.",
          fill=TRUE)
    }
    ## [x] is given by a general affine transformation from [i]
    ##
    ## [x] <- A%*%(i-1) + shift
    ## 
    ## where A <- nim@srow_[x,y,z][1:3] and shift <- srow_x,y,z[4]
    S <- array(dim=c(3,4))
    S[1,] <- nim@"srow_x"
    S[2,] <- nim@"srow_y"
    S[3,] <- nim@"srow_z"
    shift <- S[,4]
    A <- S[,1:3]
    return (A %*% (i - 1) + shift)
  }
  ## Method 1. The `old' way used only if "qform_code" is 0
  ## The co-ord mapping from [i] to [x] is the ANALYZE 7.5 way.
  ## A simple scaling relationship applies
  ##
  ## x <- pixdim[1:3] * i[1:3]
  if (verbose) {
    cat("QForm_code and SForm_code unset: Orientation by Method 1.", fill=TRUE)
  }
  ## nifti1.h pixdim[1] <- ourpixdim[2]
  scaling <- diag(pixdim(nim)[2:4])
  ## Remember i. <- i - 1
  return(scaling %*% (i - 1))
}
