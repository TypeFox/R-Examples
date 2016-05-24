plot.fmrisegment <- function(x, 
		                      anatomic = NULL, 
							  slice = 1, 
							  view = c( "axial", "coronal", "sagittal"),
							  zlim.u = NULL, 
							  zlim.o = NULL, 
							  col.o = c( rainbow( 64, start = 2/6, end = 4/6), rainbow( 64, start = 0, end = 1/6)), 
		                      col.u = grey(0:127/127),
							  verbose = FALSE,
							  ...) {
						
   if (verbose) cat( "plot.fmrisegment: entering function\n")
   view <- match.arg(view)
   
   if (is.null(anatomic)) anatomic <- array( 0, dim = x$dim[1:3])
  
   if ("fmridata" %in% class(anatomic)) {

	   if (verbose) cat( "plot.fmrisegment: calculate exact overlay\n")
	   img <- show.segmentslice(x, anatomic, slice =  slice, view = view, col.u = col.u, col.o = col.o, zlim.u, zlim.o)
	   
   } else {
	   if ( any(dim(anatomic) != dim(x$cbeta))) {
		   stop( "dimension of anatomic does not match overlay data!")
	   } else {
		   
		   if (verbose) cat( "plot.fmrisegment: perform simple overlay\n") 
		   ## select correct overlay slice according to view
		   if (view == "axial") {
			   imgdata.o <- x$cbeta[ , , slice]
			   imgdata.o[ x$segm[ , , slice] == 0] <- NA ## no significant voxel
			   imgdata.u <- anatomic[ , , slice]
			   aspect <- x$delta[2]/x$delta[1]
		   } else if (view == "coronal") {
			   imgdata.o <- x$cbeta[ , slice, ]
			   imgdata.o[ x$segm[ , slice, ] == 0] <- NA ## no significant voxel
			   imgdata.u <- anatomic[ , slice, ]
			   aspect <- x$delta[3]/x$delta[1]
		   } else if (view == "sagittal") {
			   imgdata.o <- x$cbeta[ slice, , ]
			   imgdata.o[ x$segm[ slice, , ] == 0] <- NA ## no significant voxel
			   imgdata.u <- anatomic[ slice, , ]
			   aspect <- x$delta[3]/x$delta[2]
		   } 
		   
		   ## user defined data limits to scale the image contrast
		   ## not sure whether this is what the user wants
		   if (any(!is.na(imgdata.o))) {
			   if (is.null(zlim.o)) {
				   zlim.o <- range( abs(imgdata.o), na.rm = TRUE)
			   } else {
				   if (length(zlim.o) != 2) stop("zlim.o not length 2")
				   if (zlim.o[2] < zlim.o[1]) stop("zlim.o[2] < zlim.o[1]")
				   imgdata.o[imgdata.o > zlim.o[2]] <- zlim.o[2]
				   imgdata.o[imgdata.o < zlim.o[1]] <- zlim.o[1]
				   imgdata.o[imgdata.o < -zlim.o[2]] <- -zlim.o[2]
				   imgdata.o[imgdata.o > -zlim.o[1]] <- -zlim.o[1]
			   }
		   }
		   if (is.null(zlim.u)) {
			   zlim.u <- range(imgdata.u, na.rm = TRUE)
		   } else {
			   if (length(zlim.u) != 2) stop("zlim.u not length 2")
			   if (zlim.u[2] < zlim.u[1]) stop("zlim.u[2] < zlim.u[1]")
			   imgdata.u[imgdata.u > zlim.u[2]] <- zlim.u[2]
			   imgdata.u[imgdata.u < zlim.u[1]] <- zlim.u[1]
		   }
		   
		   ## create the break points for the color scale
		   if (any(!is.na(imgdata.o))) {
			   zlim.o <- quantile( abs(imgdata.o), c( 0, 0.9, 1), na.rm = TRUE)
			   breaks.o <- c( -zlim.o[3], seq( -zlim.o[2], -zlim.o[1], length = 63), 0, seq( zlim.o[1], zlim.o[2], length = 63), zlim.o[3])
		   }
		   breaks.u <- seq( zlim.u[1], zlim.u[2], length = 129)
		   
		   ## plot the image 
		   graphics::image(1:dim(imgdata.u)[1], 1:dim(imgdata.u)[2], imgdata.u, col = col.u, asp = aspect, axes = FALSE, xlab = "",	ylab = "")
		   if (any(!is.na(imgdata.o))) {
			   graphics::image(1:dim(imgdata.o)[1], 1:dim(imgdata.o)[2], imgdata.o, asp = aspect, col = col.o, breaks = breaks.o, add = TRUE)
		   }
		   
		   ## finally create img for adimpro
		   img <- array(0, dim = c( dim(imgdata.u), 3))
		   for (i in 1:dim(imgdata.u)[1]) {
			   for (j in 1:dim(imgdata.u)[2]) {
				   if (!is.na(imgdata.o[ i, j])) { # use overlay
					   ind <- (0:128)[imgdata.o[ i, j] < breaks.o]
					   level <- ifelse(length(ind) == 0, 128, min(ind))
					   img[ i, j, ] <- as.integer( col2rgb( col.o[level])) * 256
				   } else { # use underlay
					   ind <- (0:128)[imgdata.u[ i, j] < breaks.u]
					   level <- ifelse(length(ind) == 0, 128, min(ind))
					   img[ i, j, ] <- as.integer( col2rgb( col.u[level])) * 256
				   }
			   }
		   }
		   
	   }	  
   }
		
   if (verbose) cat( "plot.fmrisegment: exiting function\n")
   return(invisible(img))
   
}

# can be called from fmri.gui or the command line (by plot)
# calls fmri.view2d or fmri.view3d with the fitting data
plot.fmridata <- function(x, anatomic = NULL , maxpvalue = 0.05, spm = TRUE,
                          pos = c(-1,-1,-1), type="slice", slice =  1, view = "axial",
                          zlim.u = NULL, zlim.o = NULL, col.o = heat.colors(256), 
                          col.u = grey(0:255/255), cutOff = c(0,1),...) {
  require(tcltk)
  if (!require(tkrplot))
    stop("required package tkrplot not found. Please install from cran.r-project.org")  
  cutOff[cutOff<0] <- 0
  cutOff[cutOff>1] <- 1
  if (cutOff[1] > cutOff[2]) cutOff[2] <- 1     
  inputStuff <- list(anatomic,maxpvalue,cutOff)

  if ("fmripvalue" %in% class(x)) {
    if ("fmridata" %in% class(anatomic)) {
      img <- show.slice(x, anatomic, maxpvalue = maxpvalue, slice =  slice, 
                        view = view, col.u, col.o, zlim.u, zlim.o)
      displayImage(img, ...)
      return(invisible(img))
    } else {
      signal <- x$pvalue
      cat("maxpvalue",maxpvalue,"\n")
      cat(sum(signal<maxpvalue),sum(signal<0.05),"mean signal",mean(signal),"\n")
      signal[signal > maxpvalue] <- 1
      cat(sum(signal<maxpvalue),"mean signal",mean(signal),"\n")
      signal[signal < 1e-10] <- 1e-10
      signal <- -log(signal)
      anatomic <- scaleAnatomic(anatomic,cutOff,dim(x$pvalue))
      # re-scale signal to 0.5 ... 1
      scale <- c(-log(maxpvalue),max(signal[is.finite(signal)]))
#      cat("scale",scale,"\n")
#      if (diff(scale)!=0) {
#        signal <- 0.5 + 0.5 * (signal - scale[1]) / diff(scale)
#      } else {
#        signal <- 0
#      }
#      # create an overlay
#      anatomic[signal >= 0.5] <- signal[signal >= 0.5]
#      anatomic[is.na(anatomic)] <- 0
#      anatomic[is.infinite(anatomic)] <- 0
      anatomic <- signalOverlay(signal,anatomic,scale)
      if (type == "3d" || type == "3D") {
        tt <- fmri.view3d(anatomic,col=mri.colors(255,255)$all,weights=x$weights,
                          scale=scale,scalecol=mri.colors(255,255)$col,
                          type= "pvalue",maxpvalue=maxpvalue,pos=pos)
      } else {
        posNew <- position(anatomic)
        fmri.view2d(anatomic,col=mri.colors(255,255)$all,weights=x$weights,scale=scale,
                    scalecol=mri.colors(255,255)$col,type= "pvalue",maxpvalue=maxpvalue,
                    posNew=position(anatomic),localx=x,inputStuff=inputStuff)
      }
    }
  } else if ("fmrispm" %in% class(x)) {
    signal <- if (spm) x$cbeta/sqrt(x$var) else x$cbeta
    # re-scale signal to 0 ... 1
    signal <- scaleSignal(signal,cutOff)
    if (type == "3d" || type == "3D") {
      if (spm) {
        tt <- fmri.view3d(signal$signal,col=mri.colors(255,255)$gray,weights=x$weights,
               scale=signal$scale,scalecol=mri.colors(255,255)$gray,type="spm",pos=pos)
      } else {
        quant <- qt(1-maxpvalue,
                 if(!is.null(x$df)) x$df else abs(diff(dim(attr(x, "design"))))) 
        tt <- fmri.view3d(signal$signal,sigma=sqrt(x$var),col=mri.colors(255,255)$gray,
                 weights=x$weights,scale=signal$scale,scalecol=mri.colors(255,255)$gray, 
                 type="spm",hrf=x$hrf, quant = quant,pos=pos)
      } 
    } else {
      fmri.view2d(signal$signal,col=mri.colors(255,255)$gray,weights=x$weights,
              scale=signal$scale,scalecol=mri.colors(255,255)$gray,type= "spm",
              maxpvalue=maxpvalue,posNew=position(x),localx=x,inputStuff=inputStuff)
    }
  } else if ("fmridata" %in% class(x)) {
    signal <- extract.data(x)
    signal <- scaleSignal(signal,cutOff)
    # re-scale signal to 0 ... 1
    if (type == "3d" || type == "3D") {
      tt <- fmri.view3d(signal$signal,col=mri.colors(255,255)$gray,weights=x$weights, 
              scale=signal$scale,scalecol=mri.colors(255,255)$gray, type="data",pos=pos)
    } else {
      fmri.view2d(signal$signal,col=mri.colors(255,255)$gray, weights=x$weights,     
                  scale=signal$scale,scalecol=mri.colors(255,255)$gray,type= "data",
                  maxpvalue=maxpvalue,posNew=position(x),localx=x,inputStuff=inputStuff)
    }
  } else {
    cat("Sorry, no plot for this class implemented\n
         Falling back to generic function, but this may fail!")
    plot(x)
  }
  if (exists("tt")) invisible(tt)
}

scaleSignal <- function(signal,cutOff){
    scale <- range(signal,finite=TRUE)
    if (diff(scale) != 0) {
      signal <-  (signal - scale[1]) / diff(scale)
    } else {
      signal <- 0
    }
    signal[is.na(signal)] <- 0
    signal[is.infinite(signal)] <- 0
    if (diff(scale) != 0){  
      signal[signal<cutOff[1]] <- cutOff[1]
      signal[signal>cutOff[2]] <- cutOff[2]    
      signal <- signal - cutOff[1]
      signal <- signal/(cutOff[2]-cutOff[1])  
    }    
    scale <- scale*(cutOff[2]-cutOff[1])
    list(signal=signal,scale=scale)
}
signalOverlay <- function(signal,anatomic,scale){
      if (diff(scale) != 0) {
        signal <- 0.5 + 0.5 * (signal - scale[1]) / diff(scale)
      } else if (scale[1] == 0) {
        signal <- 1
      } else {
        signal <- 1
      }
      # create an overlay
      anatomic[signal > 0.5] <- signal[signal > 0.5]
      anatomic[is.na(anatomic)] <- 0
      anatomic[is.infinite(anatomic)] <- 0
anatomic
}
scaleAnatomic <- function(anatomic,cutOff,ddim){
      if (is.null(anatomic)) anatomic <- array(0,dim=ddim)  
      # re-scale anatomic to 0 ... 0.5
      if (diff(range(anatomic)) !=0) {
        anatomic <- 0.5 * (anatomic - range(anatomic,finite=TRUE)[1]) /
                                  diff(range(anatomic,finite=TRUE))
      }
      anatomic[anatomic>cutOff[2]*0.5] <- cutOff[2]*0.5
      anatomic[anatomic<cutOff[1]*0.5] <- cutOff[1]*0.5
      anatomic <- anatomic - cutOff[1]*0.5
      anatomic <- anatomic/(cutOff[2]-cutOff[1])          
}
position <- function(obj){
  dt <- dim(obj)
  if(is.null(dt)) dt <- obj$dim
  c(1:dt[1],1:dt[2],1:dt[3])
}
displayImage <- function(img, ...){
      hex <- c(0:9, LETTERS[1:6])
      hex <- paste(hex[(0:255)%/%16+1],hex[(0:255)%%16+1],sep="")
      color <- paste("#",hex[img[,,1]%/%256+1],hex[img[,,2]%/%256+1],
                         hex[img[,,3]%/%256+1],sep="")
      xxx <- seq(1,dim(img)[1],length=dim(img)[1])
      yyy <- seq(1,dim(img)[2],length=dim(img)[2])
      zzz <- matrix(1:prod(dim(img)[1:2]),nrow = dim(img)[1],ncol = dim(img)[2])
      # display the image
      image(xxx, yyy, zzz, col = color, asp = 1, xlab="",ylab="", ...)
}

mri.colors <- function (n1, n2, factor=n1/(n1+n2), from=0, to=.2) {
    colors1 <- gray((0:n1)/(n1+n2))
    colors2 <- if(n2>0) hsv(h = seq(from,to,length=n2+1),
                   s = seq(from = n2/(n2+factor*n1) - 1/(2 * (n2+factor*n1)), to =
                     1/(2 * (n2+factor*n1)), length = n2+1),
                   v = 1) else NULL
    list(all=c(colors1,colors2),gray=colors1,col=colors2)
}
