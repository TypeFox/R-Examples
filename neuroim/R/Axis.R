#' @include AllGeneric.R
NULL
#' @include AllClass.R
NULL


None <- new("NamedAxis", axis="None")

NullAxis <- new("AxisSet", ndim=as.integer(0))

## add permutation vector for each NamedAxis

LEFT_RIGHT <- new("NamedAxis", axis="Left-to-Right", direction=c(1,0,0))
RIGHT_LEFT <- new("NamedAxis", axis="Right-to-Left", direction=c(-1,0,0))
ANT_POST   <- new("NamedAxis", axis="Anterior-to-Posterior", direction=c(0,-1,0))
POST_ANT   <- new("NamedAxis", axis="Posterior-to-Anterior", direction=c(0,1,0))
INF_SUP    <- new("NamedAxis", axis="Inferior-to-Superior", direction=c(0,0,1))
SUP_INF    <- new("NamedAxis", axis="Superior-to-Inferior", direction=c(0,0,-1))


matchAxis <- function(firstAxis) {
  switch(toupper(firstAxis),
         "LEFT"=LEFT_RIGHT,
         "RIGHT"=RIGHT_LEFT,
         "ANTERIOR"=ANT_POST,
         "POSTERIOR"=POST_ANT,
         "INFERIOR"=INF_SUP,
         "SUPERIOR"=SUP_INF)
         
}

TIME <- new("NamedAxis", axis="Time")

TimeAxis <- new("AxisSet1D", ndim=as.integer(1), i=TIME)

AxisSet1D <- function(i) {
  new("AxisSet1D", ndim=as.integer(1), i=i)	
}

AxisSet2D <- function(i, j) {
	new("AxisSet2D", ndim=as.integer(2), i=i, j=j)	
}

AxisSet3D <- function(i, j, k) {
	new("AxisSet3D", ndim=as.integer(3), i=i, j=j, k=k)	
}


#' permMat
#' @export
#' @rdname permMat-methods
setMethod(f="permMat", signature=signature(x = "AxisSet2D"),
          def=function(x, ...) { 
            cbind(x@i@direction, x@j@direction)
          })



#' @export
#' @rdname dropDim-methods
setMethod(f="dropDim", signature=signature(x = "AxisSet2D", dimnum="numeric"),
          def=function(x, dimnum) {  	
            stopifnot(length(dimnum) == 1)
            if (dimnum == 1) {
              AxisSet1D(x@j) 
            } else if (dimnum == 2) {
              AxisSet1D(x@i)
            } else {
              stop(paste("illegal dimnum: ",dimnum, "for axis with 2 dimensions"))
            }
          })


#' @export
#' @rdname dropDim-methods
setMethod(f="dropDim", signature=signature(x = "AxisSet2D", dimnum="missing"),
          def=function(x, dimnum) {
            AxisSet1D(x@i)
          })


#' @export
#' @rdname dropDim-methods
setMethod(f="dropDim", signature=signature(x = "AxisSet3D", dimnum="numeric"),
          def=function(x, dimnum) {    
            stopifnot(length(dimnum) == 1)
            if (dimnum == 1) {
              AxisSet2D(x@j, x@k) 
            } else if (dimnum == 2) {
              AxisSet2D(x@i, x@k)
            } else if (dimnum == 3) {
              AxisSet2D(x@i, x@j)
            } else {
              stop(paste("illegal dimnum: ", dimnum, " for axis with 2 dimensions"))
            }
          })


#' @export
#' @rdname dropDim-methods
setMethod(f="dropDim", signature=signature(x = "AxisSet3D", dimnum="missing"),
           def=function(x, dimnum) {
             AxisSet2D(x@i, x@j)
           })
           
       
        
#' @export
#' @rdname ndim-methods
setMethod(f="ndim",signature=signature(x= "AxisSet"), 
          def=function(x, ...) { 
            x@ndim 
          })


#' show an \code{NamedAxis}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("NamedAxis"), 
		def=function(object) {
			cat(print(object@axis))
		})

#' print a \code{NamedAxis}
#' @param x the object
#' @param ... extra arguments
#' @export
setMethod(f="print", signature=signature("NamedAxis"), 
		def=function(x, ...) {
			x@axis
		})

#' show an \code{AxisSet1D}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("AxisSet1D"), 
		def=function(object) {
			cat("instance of class:", class(object), "\n\n")
			cat("Axis 1:", print(object@i@axis), "\n")
		})

#' print a \code{AxisSet2D} instance
#' @param x the object
#' @param ... extra args
#' @export
setMethod(f="print", signature=signature("AxisSet2D"), 
		def=function(x, ...) {
			paste(x@i@axis, "-", x@j@axis)
		})

#' show an \code{AxisSet2D}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("AxisSet2D"), 
		def=function(object) {
			cat("instance of class:", class(object), "\n\n")
			cat("Axis 1:", object@i@axis, "\n")
			cat("Axis 2:", object@j@axis, "\n")
		})

#' print a \code{AxisSet3D} instance
#' @param x the object
#' @param ... extra args
#' @export
setMethod(f="print", signature=signature("AxisSet3D"), 
		def=function(x, ...) {
			paste(x@i@axis, " -- ", x@j@axis, " -- ", x@k@axis)
		})


#' show an \code{AxisSet3D}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("AxisSet3D"), 
		def=function(object) {
			cat("instance of class:", class(object), "\n\n")
			cat("Axis 1:", object@i@axis, "\n")
			cat("Axis 2:", object@j@axis, "\n")
			cat("Axis 3:", object@k@axis, "\n")
		})

#' show an \code{AxisSet4D}
#' @param object the object
#' @export
setMethod(f="show", signature=signature("AxisSet4D"), 
		def=function(object) {
			cat("instance of class:", class(object), "\n\n")
			cat("Axis 1:", print(object@i), "\n")
			cat("Axis 2:", print(object@j), "\n")
			cat("Axis 3:", print(object@k), "\n")
			cat("Axis 4:", print(object@l), "\n")
			
		})



OrientationList2D <- list(
	SAGITTAL_AI = AxisSet2D(ANT_POST, INF_SUP), 
	SAGITTAL_PI = AxisSet2D(POST_ANT, INF_SUP),
	SAGITTAL_PS = AxisSet2D(POST_ANT, SUP_INF),
	SAGITTAL_AS = AxisSet2D(ANT_POST, SUP_INF),
	SAGITTAL_IA = AxisSet2D(INF_SUP, ANT_POST),
	SAGITTAL_IP = AxisSet2D(INF_SUP, POST_ANT),
	SAGITTAL_SP = AxisSet2D(SUP_INF, POST_ANT),
	SAGITTAL_SA = AxisSet2D(SUP_INF, ANT_POST),

	CORONAL_LI = AxisSet2D(LEFT_RIGHT, INF_SUP),
	CORONAL_RI = AxisSet2D(RIGHT_LEFT, INF_SUP),
	CORONAL_RS = AxisSet2D(RIGHT_LEFT, SUP_INF),
	CORONAL_LS = AxisSet2D(LEFT_RIGHT, SUP_INF),
	CORONAL_IL = AxisSet2D(INF_SUP, LEFT_RIGHT),
	CORONAL_IR = AxisSet2D(INF_SUP, RIGHT_LEFT),
	CORONAL_SR = AxisSet2D(SUP_INF, RIGHT_LEFT),
	CORONAL_SL = AxisSet2D(SUP_INF, LEFT_RIGHT),


	AXIAL_LA = AxisSet2D(LEFT_RIGHT, ANT_POST),
	AXIAL_RA = AxisSet2D(RIGHT_LEFT, ANT_POST),
	AXIAL_RP = AxisSet2D(RIGHT_LEFT, POST_ANT),
	AXIAL_LP = AxisSet2D(LEFT_RIGHT, POST_ANT),
	AXIAL_AL = AxisSet2D(ANT_POST, LEFT_RIGHT),
	AXIAL_AR = AxisSet2D(ANT_POST, RIGHT_LEFT),
	AXIAL_PL = AxisSet2D(POST_ANT, LEFT_RIGHT),
	AXIAL_PR = AxisSet2D(POST_ANT, RIGHT_LEFT))

OrientationList3D <- list(
	SAGITTAL_AIL = AxisSet3D(ANT_POST, INF_SUP, LEFT_RIGHT), 
	SAGITTAL_PIL = AxisSet3D(POST_ANT, INF_SUP, LEFT_RIGHT),
	SAGITTAL_PSL = AxisSet3D(POST_ANT, SUP_INF, LEFT_RIGHT),
	SAGITTAL_ASL = AxisSet3D(ANT_POST, SUP_INF, LEFT_RIGHT),
	SAGITTAL_IAL = AxisSet3D(INF_SUP,  ANT_POST, LEFT_RIGHT),
	SAGITTAL_IPL = AxisSet3D(INF_SUP,  POST_ANT, LEFT_RIGHT),
	SAGITTAL_SPL = AxisSet3D(SUP_INF,  POST_ANT, LEFT_RIGHT),
	SAGITTAL_SAL = AxisSet3D(SUP_INF,  ANT_POST, LEFT_RIGHT),

	SAGITTAL_AIR = AxisSet3D(ANT_POST, INF_SUP, RIGHT_LEFT), 
	SAGITTAL_PIR = AxisSet3D(POST_ANT, INF_SUP, RIGHT_LEFT),
	SAGITTAL_PSR = AxisSet3D(POST_ANT, SUP_INF, RIGHT_LEFT),
	SAGITTAL_ASR = AxisSet3D(ANT_POST, SUP_INF, RIGHT_LEFT),
	SAGITTAL_IAR = AxisSet3D(INF_SUP,  ANT_POST, RIGHT_LEFT),
	SAGITTAL_IPR = AxisSet3D(INF_SUP,  POST_ANT, RIGHT_LEFT),
	SAGITTAL_SPR = AxisSet3D(SUP_INF,  POST_ANT, RIGHT_LEFT),
	SAGITTAL_SAR = AxisSet3D(SUP_INF,  ANT_POST, RIGHT_LEFT),

	CORONAL_LIA = AxisSet3D(LEFT_RIGHT, INF_SUP, ANT_POST),
	CORONAL_RIA = AxisSet3D(RIGHT_LEFT, INF_SUP, ANT_POST),
	CORONAL_RSA = AxisSet3D(RIGHT_LEFT, SUP_INF, ANT_POST),
	CORONAL_LSA = AxisSet3D(LEFT_RIGHT, SUP_INF, ANT_POST),
	CORONAL_ILA = AxisSet3D(INF_SUP,    LEFT_RIGHT, ANT_POST),
	CORONAL_IRA = AxisSet3D(INF_SUP,    RIGHT_LEFT, ANT_POST),
	CORONAL_SRA = AxisSet3D(SUP_INF,    RIGHT_LEFT, ANT_POST),
	CORONAL_SLA = AxisSet3D(SUP_INF,    LEFT_RIGHT, ANT_POST),

	CORONAL_LIP = AxisSet3D(LEFT_RIGHT, INF_SUP, POST_ANT),
	CORONAL_RIP = AxisSet3D(RIGHT_LEFT, INF_SUP, POST_ANT),
	CORONAL_RSP = AxisSet3D(RIGHT_LEFT, SUP_INF, POST_ANT),
	CORONAL_LSP = AxisSet3D(LEFT_RIGHT, SUP_INF, POST_ANT),
	CORONAL_ILP = AxisSet3D(INF_SUP,    LEFT_RIGHT, POST_ANT),
	CORONAL_IRP = AxisSet3D(INF_SUP,    RIGHT_LEFT, POST_ANT),
	CORONAL_SRP = AxisSet3D(SUP_INF,    RIGHT_LEFT, POST_ANT),
	CORONAL_SLP = AxisSet3D(SUP_INF,    LEFT_RIGHT, POST_ANT),


	AXIAL_LAI = AxisSet3D(LEFT_RIGHT, ANT_POST, INF_SUP),
	AXIAL_RAI = AxisSet3D(RIGHT_LEFT, ANT_POST, INF_SUP),
	AXIAL_RPI = AxisSet3D(RIGHT_LEFT, POST_ANT, INF_SUP),
	AXIAL_LPI = AxisSet3D(LEFT_RIGHT, POST_ANT, INF_SUP),
	AXIAL_ALI = AxisSet3D(ANT_POST,   LEFT_RIGHT, INF_SUP),
	AXIAL_ARI = AxisSet3D(ANT_POST,   RIGHT_LEFT, INF_SUP),
	AXIAL_PLI = AxisSet3D(POST_ANT,   LEFT_RIGHT, INF_SUP),
	AXIAL_PRI = AxisSet3D(POST_ANT,   RIGHT_LEFT, INF_SUP),

	AXIAL_LAS = AxisSet3D(LEFT_RIGHT, ANT_POST, SUP_INF),
	AXIAL_RAS = AxisSet3D(RIGHT_LEFT, ANT_POST, SUP_INF),
	AXIAL_RPS = AxisSet3D(RIGHT_LEFT, POST_ANT, SUP_INF),
	AXIAL_LPS = AxisSet3D(LEFT_RIGHT, POST_ANT, SUP_INF),
	AXIAL_ALS = AxisSet3D(ANT_POST,   LEFT_RIGHT, SUP_INF),
	AXIAL_ARS = AxisSet3D(ANT_POST,   RIGHT_LEFT, SUP_INF),
	AXIAL_PLS = AxisSet3D(POST_ANT,   LEFT_RIGHT, SUP_INF),
	AXIAL_PRS = AxisSet3D(POST_ANT,   RIGHT_LEFT, SUP_INF))
	
	
	
#' given three named axes return AxisSet3D singleton
#' @param axis1 the first axis
#' @param axis2 the second axis
#' @param axis3 the third axis	
#' @export 
matchAnatomy3D <- function(axis1, axis2, axis3) {
	for (orient in OrientationList3D) {
		if (identical(orient@i,axis1) && identical(orient@j,axis2) && identical(orient@k, axis3)) {
			return(orient)
		}
	}
	
	stop("did not find legal matching anatomical orientation for axes: ",
			axis1, axis2, axis3)
		
}

#' given two named axes return AxisSet2D singleton
#' @param axis1 the first axis
#' @param axis2 the second axis
#' @export  
matchAnatomy2D <- function(axis1, axis2) {
	for (orient in OrientationList2D) {
		if (identical(orient@i,axis1) && identical(orient@j,axis2)) {
			return(orient)
		}
	}
	
	stop("did not find legal matching anatomical orientation for axes: ",
			axis1, axis2)
	
}


.nearestAnatomy <- function(mat44) {
  mat33 <- mat44[1:3, 1:3]
  #mat33 <- sweep(mat33, 2, sqrt(apply(mat33 * mat33, 2, sum)), "/")
  icol <- mat33[,1]
  jcol <- mat33[,2]
  kcol <- mat33[,3]
  
  ## normalize icol
  icol <- icol / sqrt(sum(icol^2))
  
  ## normalize jcol
  jcol <- jcol / sqrt(sum(jcol^2))
  
  orthogonalize <- function(col1, col2) {
    dotp <- sum(col1*col2)
    if (abs(dotp > 1.e-4)) {
      col2 <- col2 - (dotp * col1)
      norm <- sqrt(sum(col2^2))
      col2 <- col2 / norm
    }
    col2
  }
  
  jcol <- orthogonalize(icol, jcol)
   
  knorm <- sqrt(sum(kcol^2))
  if (knorm == 0.0) {
    kcol[1] <- icol[2] * jcol[3] - icol[3] * jcol[2]
    kcol[2] <- icol[3] * jcol[1] - jcol[3] * icol[1]
    kcol[3] <- icol[1] * jcol[2] - icol[2] * jcol[1]     
  } else {
    kcol <- kcol / knorm
  }
  
  ## orthogonalize k to i
  kcol <- orthogonalize(icol, kcol)
  kcol <- orthogonalize(jcol, kcol)
  
  Q <- cbind(icol, jcol, kcol)
  P <- matrix(0,3,3)
  detQ <- det(Q)
  if (detQ == 0.0) {
    stop("invalid matrix input, determinant is 0")
  }
  
  vbest = -666
  ibest = 1
  pbest = 1
  qbest = 1
  rbest = 1
  
  jbest = 2
  kbest = 3
  for (i in 1:3) {     
    for (j in 1:3) {
      if (i == j) next
      for (k in 1:3) {
        if (i == k || j == k) next
          P <- matrix(0,3,3)
          for (p in c(-1,1)) {
            for (q in c(-1,1)) {
              for (r in c(-1,1)) {
                P[1, i] <- p
                P[2, j] <- q
                P[3, k] <- r
                detP <- det(P)
                if (detP * detQ <= 0.0) next
                M <- P %*% Q
                crit <- sum(diag(M))
                if (crit > vbest) {
                  vbest = crit
                  ibest = i
                  jbest = j
                  kbest = k
                  pbest = p
                  qbest = q
                  rbest = r
                }                         
              }
            }
          }
      }
    }
  }
  
  .getAxis <- function(num) {
    switch(as.character(as.integer(num)),
                "1"=LEFT_RIGHT,
                "-1"=RIGHT_LEFT,
                "2"=POST_ANT,
                "-2"=ANT_POST,
                "3"=INF_SUP,
                "-3"=SUP_INF)
  }
  
  ax1 <- .getAxis(ibest*pbest)
  ax2 <- .getAxis(jbest*qbest)
  ax3 <- .getAxis(kbest*rbest)
  
  AxisSet3D(ax1,ax2,ax3)
  
  
                
  
}

