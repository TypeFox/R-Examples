#' @title Create a movie of plotRSA plots, with changing surface and/or rotation
#'
#' @description
#' Create a movie of plotRSA plots, with changing surface and/or rotation
#'
#' @details
#' \code{frames} is a list of the first, intermediate, and the final parameters of the surface. Each scalar parameter defined in \code{frames} is interpolated between steps in order to create a smooth sequence of plots. Logical and character parameters are inherited from the first frame. Plots are saved as individual still pictures in a subfolder called \code{name} and finally glued together using ffmpeg. Hence, a ffmpeg installation is needed to create the movie (the still pictures can be produced without ffmpeg).
#'
#' @export
#' @param name Name for the subfolder containing all still pictures, and for the final movie file.
#' @param frames A list of lists: Each list contains parameters which are passed to the plotRSA function. See \code{\link{plotRSA}} for details.
#' @param dur Duration of the movie in miliseconds
#' @param fps Frame per second (defaults to 30)
#' @param width Width of the final movie in pixels
#' @param height Height of the final movie in pixels
#' @param mirror If \code{TRUE}, the frame sequence is mirrored at the end so that the movie ends at frame 1.
#' @param savetodisk If \code{TRUE} the files are saved to the disk. If \code{FALSE}, the movie is only shown on the screen
#' @param clean Should the still images be deleted?
#'
#' @seealso \code{\link{plotRSA}}
#'
#' @examples
#' \dontrun{
#'movieRSA(name="SD0",
#' 		frames <- list(
#'		 	step1 = list(b0=0, xy=-.40, x2=.20, y2=.20, 
#'				rotation=list(x=-63, y=32, z=15),
#'				legend=FALSE, zlim=c(0, 4), param=FALSE),
#'		 	step2 = list(b0=0, xy=-.10, x2=.05, y2=.05, 
#'				rotation=list(x=-54, y=39, z=25)),
#'			step3 = list(b0=0, xy=-.40, x2=.20, y2=.20, 
#'				rotation=list(x=-45, y=45, z=35))
#'		 ),
#'		 mirror=TRUE, fps=30, dur=5000)
#' }


movieRSA <- function(name, frames, dur=2000, fps=30, width=800, height=600, mirror=TRUE, savetodisk=TRUE, clean=TRUE) {
	
	if (dir.create(name) == FALSE) 
		stop("Could not create directory for movie pictures!")
	
	# mirror keyframes if requested
	if (mirror==TRUE) {
		L <- length(frames)
		for (i in 1:(length(frames)-1)) {
			frames[[L+i]] <- frames[[L-i]]
		}
	}
	
	# compute number of still pictures, depending on duration and frame rate. Adjust to make whole integers for scenes
	max.n <- round(round((dur/1000)*fps)/(length(frames)-1))*(length(frames)-1)
	#max.n <- 10
	VAR <- names(frames[[1]])	# v keeps the variable names used in the frames
	
	# Produce transitions/ interpolations
	trans <- list()
	s <- 1
	for (i in 1:max.n) {
		trans[[i]] <- frames[[1]]
		scene.n <- max.n/(length(frames)-1)	# get number of frames within that transition
		
		# get the frame number from which to interpolate
		startFrame <- cut(1:max.n, length(frames)-1, labels=FALSE)[i]
		endFrame <- startFrame + 1
		
		# initialize secondary counter
		if (s > scene.n) {s <- 1} # s is the # of the current transition ("scene")

		for (v in VAR) {
			# scalar in frames: make interpolation
			if (is.numeric(frames[[startFrame]][[v]]) & length(frames[[startFrame]][[v]]) == 1) {
				trans[[i]][[v]] <- interpol(frames[[startFrame]][[v]], frames[[endFrame]][[v]], s, scene.n)
			} else if (is.list(frames[[startFrame]][[v]]) & (v %in% c("rotation", "label.rotation"))) {
				V2 <- names(frames[[startFrame]][[v]])
				for (v2 in V2) {
					trans[[i]][[v]][v2] <- interpol(frames[[startFrame]][[v]][[v2]], frames[[endFrame]][[v]][[v2]], s, scene.n)
				}
			} else {
				# everything else: keep parameters from first frame...
				trans[[i]][[v]] <- frames[[1]][[v]]
			}
		} # of v
			
		s <- s + 1
	}
	
	
	p <- list()				# p is a list of the produced plots
	for (i in 1:max.n) {
		print(paste0(i,"/",max.n))
		p[[i]] <- do.call(plotRSA, trans[[i]])
		
		if (savetodisk==TRUE) {
			picName <- paste0(name, "/", name, "_", sprintf("%04.0f", i), ".jpg")
			jpeg(picName, width=width, height=height, quality=95)			
		}
		print(p[[i]])
		if (savetodisk==TRUE) {dev.off()}
		print(p[[i]])
	}
	
	# finally: build movie
	if (savetodisk==TRUE) {
		# delete any existing movie file
		unlink(paste(name,".avi",sep=""))
	
		# point system to R's working directory
		system(paste0("cd ", gsub(" ", "\\ ", getwd(), fixed=TRUE)))
	
		# show & execute the command line expression for ffmpeg to glue the pictures together
		 print(paste(paste0("ffmpeg -r ", fps, " -i ", name, "/", name, "_%04d.jpg -sameq -r ", fps," -b:v 5000K ",  paste0(name, "/", name, ".avi"))))
		system(paste(paste0("ffmpeg -r ", fps, " -i ", name, "/", name, "_%04d.jpg -sameq -r ", fps," -b:v 5000K ",  paste0(name, "/", name, ".avi"))))
		
		if (clean==TRUE) {
			file.remove(Sys.glob(paste0(name, "/", "*.jpg")))
		}
	}
	
 }
 
 
 
# helper function interpolation function
interpol <- function(x1, x2, step, max.n) {
	return(x1 + seq(0, x2-x1, length.out=max.n)[step])
}