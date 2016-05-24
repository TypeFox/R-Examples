#'Generate a plant with randomly distributed leaves
#'
#'@description Generates a plant crown with specified total leaf area or number of leaves,
#'in one of six crown shapes. These include cone, paraboloid, cylinder, and
#'others. Leaves are randomly distributed in the crown shapes.
#'
#'Six crown shapes are implemented, use one of the following abbreviations for
#'the argument \code{crownshape}: 
#'\describe{ 
#'\item{BOX}{A box shape; radius in X and Y directions assumed to be the same.} 
#'\item{CONE}{A cone shape.}
#'\item{HELIP}{A half ellipsoid (the top half of a full ellipsoid).}
#'\item{PARA}{A paraboloid.} 
#'\item{ELIP}{A full ellipsoid.} 
#'\item{CYL}{Cylinder shape.} }
#'
#'The leaf angle distribution (\code{LAD}) must be specified using the
#'\code{\link{angledist}} function in the \code{LeafAngle} package. To assign a
#'constant leaf angle of 45 degrees, use this command: \preformatted{ LAD =
#'angledist("const", 45) } To assign a spherical LAD, use this command:
#'\preformatted{ LAD = angledist("spherical") }
#'
#'@param nleaves Number of leaves, can be left unspecified if leaf area (LA) is
#'given.
#'@param leaflen The leaf length (constant value for all leaves).
#'@param LA Total leaf area, optional (m2)
#'@param radius Crown radius (mm).
#'@param height Crown height (mm).
#'@param crownbase Height to crownbase (mm).
#'@param crownshape One of six pre-specified crown shapes. See Details.
#'@param LAD Leaf angle distribution.
#'@param lfile Name of the leaf file (see \code{\link{readl}}), if left
#'unspecified a triangle is used.
#'@param writefile Logical. If TRUE, writes a Q file to disk.
#'@param quiet Logical. If FALSE, no messages are printed to the console.
#'@return Produces a plant in the Q file format, see
#'\code{\link{constructplant}} for details.  If \code{writefile=TRUE}, the Q
#'file is written to the current working directory. Otherwise, returns a
#'dataframe that can be used directly in \code{constructplant}, see Details.
#'@note Currently only random distributions are supported. I have at this time
#'no idea how to produce non-random (clumped or uniform) distributions in 3D!
#'@author Remko Duursma. Uses code from Belinda Medlyn (SURFACE subroutine in
#'MAESTRA).
#'@seealso \code{\link{constructplant}}
#'@keywords misc
#'@examples
#'
#'\dontrun{
#'# An ellipsoid shape crown, write a Q file to disk.
#'# Specify approximate total leaf area.
#'randomplant(radius=400, height=2000, shape="ELIP",
#'leaflen=40, LA=1.5, writefile=TRUE)
#'
#'# Embed the function in a call to 'constructplant', giving a plant in the 
#'# YplantQMC format.
#'# Because no leaf file is specified, it uses a built-in triangle-shaped leaf.
#' coneplant <- constructplant(randomplant(radius=250, height=500, shape="CONE",
#' leaflen=25, LA=1))
#'
#'}
#'
#'@export
#'@importFrom LeafAngle drawsample
randomplant <- function(nleaves=500, leaflen=30, LA=NULL, 
                        radius=250, height=500,crownbase=0,
						            crownshape=c("BOX","CONE","HELIP","PARA","ELIP","CYL"),
						            LAD=LeafAngle::angledist("spherical"),
                        lfile=NULL,
						            writefile=FALSE, quiet=FALSE
						){
						
					
	shape <- match.arg(crownshape)

	# Get leaf shape from lfile, otherwise it will use triangleleaf (shape=0.5).
	if(!is.null(lfile)){
	  l <- readl(lfile)
	  if(length(l)>1 && !quiet)
	    message("Using first leaf type in leaf file only.")
	  l <- l[[1]]
	  leafshape <- l$leafshape
	} else {
	  leafshape <- 0.5
	}
  
	if(!is.null(LA)){
		laleaf <- leaflen^2 * leafshape

		nleaves <- 10^6 * LA / laleaf
		if(!quiet)message("Aiming for leaf area = ", signif(LA,4), " m2")
	}
	

	w <- radius
	
	# Make better guess of number of leaves:
	# (crown volume relative to box).
	boxvol <- height * (radius*2)^2
	crownvol <- switch(shape,
			BOX = boxvol,
			CONE = 1/3 * pi * height * radius^2,   
			HELIP = 2/3 * pi * height * radius^2,
			PARA = 1/2 * pi * height * radius^2,
			ELIP = 4/3 * pi * (height/2) * radius^2,
			CYL = pi * radius^2 * height)
	relcrownvol <- crownvol / boxvol
	
	n <- floor(nleaves / relcrownvol)  # approximate!
	
	xyz <- data.frame(X=runif(n, -w,w), 
                  Y = runif(n, -w,w),
				  Z = crownbase + runif(n, 0,height),
				  ang = drawsample(LAD,n,degrees=TRUE),
				  AZ = runif(n,0,360),
				  OR = runif(n,0,360),
				  len = leaflen)

	if(shape != "BOX"){

		# SURFACE from Maestra:
		crownrh <- function(H,JSHAPE=c("CONE","HELIP","PARA","ELIP","CYL")){
			JSHAPE <- match.arg(JSHAPE)
			
			relr <- switch(JSHAPE,
				CONE = 1-H,
				HELIP = sqrt(1-H^2),
				PARA = sqrt(1-H),
				ELIP = sqrt(1 - ((H-1/2)^2)/((1/2)^2)),
				CYL = 1)
			return(relr)
		}
		
		# Select only leaves that fall within the shaped hull.
		maxh <- max(xyz$Z)
		relr <- sqrt(xyz$X^2 + xyz$Y^2) / w
		crownr <- crownrh(xyz$Z / maxh, shape)
		xyz <- xyz[relr < crownr,]
		
	}

	#
	LAactual <- 10^-6 * sum(leafshape * xyz$len^2)
	if(!quiet){
		message("Actual leaf area = ", signif(LAactual,4), " m2")
		message("Number of leaves = ", nrow(xyz))
	}		  
	
	# Q file format (see constructplant()).	
	qfile <- xyz

	if(!writefile)
		return(qfile)
	else {
		filen <- paste(shape,"_",n,".Q")
		write.table(qfile, filen, row.names=FALSE, col.names=TRUE)
		return(invisible(qfile))
	}

}
