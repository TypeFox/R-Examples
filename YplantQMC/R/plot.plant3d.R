#'Plots a plant in 3D
#'
#'@description Based on a plant constructed with \code{\link{constructplant}}, plots the
#'plant in 3D using the rgl package. Optionally adds the convex hull to the
#'plot.
#'
#'See the function \code{\link{crownhull}} for details of the convex hull
#'algorithm.
#'
#'@details Some of the hard work, figuring out the X,Y,Z coordinates of the leaf tip and
#'leaf sides, is directly taken from the original Yplant Delphi code (in the
#'non-visible function \code{madeleafdirection}).
#'
#'
#'@aliases plot.plant3d plot.plant3dlist plot.stand3d
#'@param x Object of class \code{plant3d} or \code{plant3dlist}.
#'@param noleaves Logical. If TRUE, leaves are not plotted.
#'@param squarewidth Width of the square to be plotted at the ground surface.
#'@param addcrownhull Logical. If TRUE, a convex hull is added to the plot.
#'@param hullalpha Opacity of the convex hull.
#'@param leaffill Logical. If TRUE, the default, fills in the leaves.
#'@param leafoutline Logical. If TRUE, plots leaf outlines.
#'@param stems Logical. If TRUE, plots the stems (black) and branches .
#'@param cylinderstems Logical. If TRUE, plots the stems as cylinder sections
#'given the diameter in the plant file.
#'@param leafcolor Color of the filled leaves. Can only specify one color.
#'@param markleaves,markcolor A vector of leaf numbers to 'mark', using
#''markcolor'.
#'@param stemcol,branchcol,petiolecol Color of the stems, branches, petioles
#'@param add Logical. If TRUE (defaults to FALSE), adds this plant to an
#'existing scene.
#'@param shiftxyz Vector of x,y,z coordinates to shift the plant (see
#'Examples).
#'@param whichrows Which of the plants in the \code{plant3dlist} object to
#'plot. Defaults to all.
#'@param png If TRUE, saves a PNG snapshot of every plant plotted.
#'@param pngsuff A suffix to add to the filename (defaults to pfile.PNG).
#'@param keepopen Whether to keep all windows open. Max number of windows is
#'10.
#'@param windowrect Size of the \code{rgl} window (see \code{\link{open3d}}).
#'Affects PNG quality.
#'@param skipexisting If PNG=TRUE, whether to skip plants that are already
#'saved to disk.
#'@param addfilename If TRUE, places a label of the pfile on the plot
#'(currently poorly placed).
#'@param \dots Further arguments passed to \code{\link{plot3d}} in the
#'\code{rgl} package.
#'@author Remko Duursma. Robert Pearcy for the Y-plant code.
#'@seealso \code{\link{crownhull}}
#'@keywords misc
#'@examples
#'
#'\dontrun{
#'# Construct a plant like this:
#'# myplant <- constructplant("somepfile.p","aleaf.l")
#'
#'# Standard plot (using built-in sugar maple):
#'plot(sugarmaple)
#'
#'# Add a convex hull.
#'plot(sugarmaple, addcrownhull=TRUE)
#'
#'# Plot two plants in one scene (second plant is moved 750mm in the x-direction)
#'plot(sugarmaple)
#'plot(sugarmaple, shiftxyz=c(750,0,0), add=TRUE)
#'
#'# Mark the first 10 leaves in red (i.e. first 10 leaves in P file):
#'plot(sugarmaple, markleaves=1:10)
#'
#'# Mark all leaves on the plant that have a leaf angle > 45.
#'plot(toona, markleaves=which(toona$leafdata$ang > 45))
#'
#'# Plot the stems (and branches) only:
#'plot(pilularis, noleaves=TRUE)
#'
#'# Plot many plants.
#'# First organize a 'leafplantkey', a comma-separated file that links each pfile to an lfile.
#'# See the online manual for an example, or the help file for 'constructplant'.
#'myplants <- readplantlist(lpk="leafplantkey.txt")
#'plot(myplants, png=TRUE, addfilename=TRUE)
#'
#'}
#'


#'@rdname plot.plant3d
#'
#'@method plot plant3d
#'@S3method plot plant3d
#'@importFrom rgl lines3d
#'@importFrom rgl segments3d
#'@importFrom rgl polygon3d
plot.plant3d <- function(x,
					  noleaves=FALSE,
					  squarewidth=250,
					  addcrownhull=FALSE,
					  hullalpha=0.4,
					  leaffill=TRUE,
					  leafoutline=TRUE,
					  stems=TRUE,
					  cylinderstems=stems,
					  leafcolor="forestgreen",
					  markleaves=NULL,
					  markcolor="red",
					  stemcol="black",
					  branchcol="black",
					  petiolecol="brown",
					  add=FALSE,shiftxyz=c(0,0,0),
					  ...
					  ){  

  
	plant <- x
  
  if(!all(shiftxyz==0)){
    plant <- shiftplant(plant, shiftxyz[1],shiftxyz[2],shiftxyz[3])
  }
  
	if(noleaves){
		leaffill <- FALSE
		leafoutline <- FALSE
	}
	
	inputformat <- plant$inputformat
	if(stems && inputformat == "Q")stems <- FALSE  # Q format does not support stems.

	
	# Plot grey square on the ground:
	if(!add)open3d()
	s <- squarewidth
	M <- matrix(c(-s,-s,0,
				  s,-s,0,
				  s,s,0,
				  -s,s,0,
				  -s,-s,0), ncol=3, byrow=TRUE)
	if(s > 0)lines3d(M, col="darkgrey", add=add)

	# Add crownhull, maybe
	if(addcrownhull){
  			ch <- crownhull(plant, alpha=hullalpha)
  }
	
	# Plot stems, branches, and petioles.
	if(plant$inputformat == "P"){
	
	Nnodes <- nrow(plant$pdata)
	
	# Stem segments.
	if(stems){
		Ms <- list()
		for(i in 1:Nnodes){
			Ms[[i]] <- rbind(rbind(plant$stems[[i]]$xyz$from, plant$stems[[i]]$xyz$to))
		}
		Ms <- do.call("rbind", Ms)
		segments3d(Ms, col=stemcol)
			
		# Branches.
		Ms <- list()
		for(i in 1:Nnodes){
			Ms[[i]] <- rbind(rbind(plant$branches[[i]]$xyz$from, plant$branches[[i]]$xyz$to))
		}
		Ms <- do.call("rbind", Ms)
		segments3d(Ms, col=branchcol)
	}
	# Add cylinder sections.
	if(stems & cylinderstems)plotstemsections(plant)

	
	# Petioles.
	Ms <- list()
	for(i in 1:Nnodes){
		Ms[[i]] <- rbind(rbind(plant$petioles[[i]]$xyz$from, plant$petioles[[i]]$xyz$to))
	}
	Ms <- do.call("rbind", Ms)
	segments3d(Ms, col=petiolecol)
	}
	
	if(leafoutline){
		LM <- list()
		np <- nrow(plant$leaves[[1]]$XYZ)
		for(i in 1:plant$nleaves){
			LM[[i]] <- plant$leaves[[i]]$XYZ
			LM[[i]] <- rbind(LM[[i]], LM[[i]][1,])  # duplicate first point to complete polygon.
			if(np %% 2 > 0)LM[[i]] <- rbind(LM[[i]], LM[[i]][np,])
			nr <- nrow(LM[[i]])
			LM[[i]] <- LM[[i]][rep(1:nr,each=2),]
			LM[[i]] <- LM[[i]][-c(1,nr*2),]
		}
		LM <- do.call("rbind", LM)
		segments3d(LM, col="black")
	}
	

	if(leaffill & !noleaves){
    
    if(is.null(markleaves)){
      markleaves <- 1:plant$nleaves
    }
    
    for(i in markleaves){      
      tr <- try(polygon3d(plant$leaves[[i]]$XYZ, col=leafcolor),silent=TRUE)
      
      # If failed, try one more time (see Note in rgl::triangulate)
      if(inherits(tr, "try-error")){
        tr <- try(polygon3d(plant$leaves[[i]]$XYZ, col=leafcolor),silent=TRUE)
      }
    }
    
	}

options(warn=0)
  return(invisible(LM))
}



#'@rdname plot.plant3d
#'
#'@method plot plant3dlist
#'@S3method plot plant3dlist
#'@importFrom rgl open3d
#'@importFrom rgl rgl.close
#'@importFrom rgl snapshot3d
#'@importFrom rgl title3d
plot.plant3dlist <- function(x, whichrows=NA, png=FALSE, pngsuff="", keepopen=!png,
                             windowrect=c(50,50,800,800), squarewidth=25, skipexisting=FALSE, addfilename=FALSE, ...){
  
  plants <- x
  pfiles <- attr(plants, "pfiles")
  lfiles <- attr(plants, "lfiles")
  nplants <- attr(plants, "nplants")
  
  if(whichrows == "all" || is.na(whichrows))whichrows <- 1:nplants
  
  # nplants is in attributes of plants
  
  if(keepopen && length(whichrows)>10)stop("You don't want this many windows open ... set png=TRUE or reduce number of plots.")
  
  for(k in whichrows){
    
    pfilename <- pfiles[k]
    lfilename <- lfiles[k]
    
    outputname <- paste(pfilename,pngsuff,".png",sep="")
    
    # skipexisting=TRUE skips PNGs that are already produced.
    if(skipexisting && file.exists(outputname))next
    
    open3d(windowRect=windowrect)
    p <- try(plot(plants[[k]], add=TRUE, squarewidth=squarewidth, ...))
    if(inherits(p,"try-error")){
      warning("Failed to plot",pfilename)
      rgl.close()
      next
    }
    if(addfilename)title3d(sub=pfiles[k], line=-1, col="black")
    if(png)snapshot3d(outputname)
    if(!keepopen)rgl.close()
  }
}


#'@rdname plot.plant3d
#'
#'@method plot stand3d
#'@S3method plot stand3d
plot.stand3d <- function(x,...){
  
  n <- length(x$plants)
  
  # open rgl canvas
  open3d()
  
  # Draw plot boundary
  pb <- x$plotbox
  M <- matrix(c(pb[1],pb[2],0,
                pb[3],pb[2],0,
                pb[3],pb[4],0,
                pb[1],pb[4],0,
                pb[1],pb[2],0), ncol=3, byrow=TRUE)
  lines3d(M, col="darkgrey")
  
  # Add plants
  for(i in 1:n)plot(x$plants[[i]], 
                    add=TRUE, 
                    squarewidth=0, 
                    ...)
  
}
