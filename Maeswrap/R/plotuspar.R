#' Plot the understorey PAR points
#' 
#' Reads the 'uspar.dat' file in the current working directory, and plots the
#' incident or absorbed (or diffuse, or direct) PAR at the understorey points.
#' Either produces a plot, or makes a pdf with a plot for each hour of the
#' selected day.
#' 
#' If addNarrow is TRUE, attempts to read the trees.dat file in the current
#' working directory as well. Prints a warning when this file cannot be opened.
#' 
#' @param what Either 'diff', 'apar', 'ipar', or 'beam' (the default).
#' @param dataset If left alone, reads the uspar dataset.
#' @param day Which day to use, if left alone uses the first day only.
#' @param hour Which hour to plot. If left alone, makes a plot for each hour.
#' @param xlim,ylim X- and Y-axis limits.
#' @param makepdf Logical. If TRUE, produces a pdf in the working directory.
#' @param outputfile Name of the pdf file.
#' @param scaleeach Logical. Rescale grey scale for each plot, or same for all
#' hours?
#' @param addNarrow Logical. Add an arrow pointing North.
#' @return A lattice device, or a pdf.
#' @author Remko Duursma
#' @export
#' 
#' @examples
#' \dontrun{
#' 
#' # Plot one hour of the first day, showing incident PAR on understorey:
#' plotuspar("ipar", day=1,hour=12,makepdf=FALSE)
#' 
#' # Make pdf of the whole day, plotting beam radiation:
#' plotuspar("beam", day=1, outputfile="beam uspar")
#' 
#' }
#' @importFrom tools file_ext
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @rdname plotuspar
plotuspar <- function(what=c("PARbeam","PARtotal","PARdiffuse","APAR"), dataset=NULL, 
                      day=1, hour=NA, xlim=NULL, ylim=NULL,
					  makepdf=FALSE, outputfile = "aparunderstorey.pdf", scaleeach=TRUE,
					  addNarrow = TRUE
){

	what <- match.arg(what)

	if(addNarrow)Bearing <- try(readPAR("trees.dat", "bearing", "plot"),silent=TRUE)
	if(addNarrow && inherits(Bearing, "try-error")){
		warning("Could not read bearing\n")
		addNarrow <- FALSE
	}
	
	x0 <- try(readPAR("ustorey.dat", "X0", "control"))
	y0 <- try(readPAR("ustorey.dat", "Y0", "control"))
	xmax <- try(readPAR("ustorey.dat", "XMAX", "control"))
	ymax <- try(readPAR("ustorey.dat", "YMAX", "control"))
  
	if(any(sapply(c(x0,y0,xmax,ymax),inherits,"try-error"))){
		boxdraw <- FALSE
		warning("Could not read x0,y0,xmax and/or ymax.")
	} else boxdraw <- TRUE
	
	if(is.null(dataset)){
        dataset <- readuspar()
  }
  
  dataset$Z <- dataset[,what]
	nhour <- max(dataset$hour)
	
  if(!is.na(day))
	  dataset <- dataset[dataset$day==day,]
	
	if(!is.na(hour))
		dataset <- dataset[dataset$hour==hour,]
		
	nhours <- length(unique(dataset$hour))
	
  z <- list();k <- 1
  
  for(ihour in unique(dataset$hour)){
      
    DATA <- dataset[dataset$hour == ihour,]
    
    if(nhour == 24){
      TIME <- format(strptime(as.character(ihour), format="%H"), format="%H:%M")
    }
    if(nhour == 48){
       HOUR <- as.character(ihour %/% 2)
       MIN <- as.character( 30 * ihour %% 2)
       hourmin <- paste(HOUR,MIN, sep=" ")
       TIME <- format(strptime(hourmin, format="%H %M"), format="%H:%M")
    }
      
	  if(is.null(xlim))xlim <- c(0,max(DATA$X))
	  if(is.null(ylim))ylim <- c(0,max(DATA$Y))
	  
	  X0 <- xlim[1] + 0.1*(xlim[2] - xlim[1])
	  Y0 <- ylim[1] + 0.1*(ylim[2] - ylim[1])
	  LEN <- 0.1 * (xlim[2] - xlim[1])
	  
	  if(scaleeach)
		  AT <- seq(0,max(DATA$Z),length=50)
	  else
		  AT <- seq(0,max(dataset$Z),length=50)
		
	  if(all(AT == 0))next
	  	  
    z[[k]] <- with(DATA, levelplot(Z ~ X*Y, 
        main=TIME, xlim=xlim,ylim=ylim,
        at=AT, panel=function(x,y,...){
			                  panel.levelplot(x,y,...)
            			if(addNarrow)addarrow(X0,Y0,len=LEN,bearing=Bearing,addto="lattice")
            			if(boxdraw)lrect(x0,y0,xmax,ymax,border="darkgrey")
            		},
        col.regions=grey(seq(0,1,length=50))))
      k <- k + 1
    }
	
	if(makepdf && file_ext(outputfile) != "pdf"){
		outputfile <- paste(outputfile,".pdf",sep="")
	}
  
  if(makepdf){
    pdf(outputfile,onefile=TRUE)
    on.exit(dev.off())
  }
    for(i in seq_along(z)){
      print(z[[i]])
    }

}
