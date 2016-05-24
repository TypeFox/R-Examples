#'@rdname ModifyPfiles
#'@title Modify Yplant input files.
#'
#'@description Reads a plant file (a Yplant input file with extension \code{.p}, known as a
#'\code{pfile}), and writes a new pfile, by modifying angles, internode
#'lengths, or segment diameters.
#'
#'The function \code{modifypfile} is a general function that can be used to modify
#'any variable in a \code{pfile}. The following variables can be changed in a
#'\code{pfile} :
#'
#'\describe{ 
#'\item{Az, An}{ Azimuth and angle of stem sections} 
#'\item{Or}{Orientation angle of the midrib of the leaf} 
#'\item{Az.1, An.1}{ Branch sections} 
#'\item{Az.2, An.2}{ Petioles } 
#'\item{Az.3, An.3}{ Leaves } }
#'
#'The function \code{replaceangles} is a specialized version for replacing leaf
#'angles (and do some specific error checking, or sampling from leaf angle
#'distributions). The function \code{changeinternodes} is a special function
#'for lengths of woody segments (technically not quite the same as
#'internodes!).
#'
#'For \code{replaceangles}, an object of class \code{angledist} can be
#'constructed with the \code{\link{angledist}} function in package
#'\code{LeafAngle} (See Example below) or from fitting to a sample with
#'\code{\link{fitdistribution}}.  If angles are provided but the vector is too
#'short, the vector will be sampled with replacement.
#'
#'@aliases replaceangles changeinternodes modifypfile
#'@param whichangle By default the leaf angle, or one of the other angles (see
#'Details).
#'@param pfile Name of the .p file.
#'@param outputfile Name of the new file.
#'@param distobj An object of class \code{\link{angledist}}, that is a leaf
#'angle distribution (See examples).
#'@param newangles A vector of angles to be used in the new p file (optional,
#'in stead of distobj).
#'@param whichvar Name of variable in p file to modify (see \code{\link{readp}}
#'for list of variables.
#'@param newvalues Vector of new values for the variable that is to be
#'replaced.
#'@param method for \code{changeinternodes}, a percentage change ("perc") or a
#'new constant value ("constant").
#'@param changeperc If method="perc", change the internodes by this percentage
#'of their orginal value.
#'@param consvalue If method="constant", change all internodes to this constant
#'value.
#'@return A new \code{pfile} is created, by default with the name "tmp.p",
#'unless the argument \code{outputfile} is set.
#'@author Remko Duursma
#'@seealso \code{\link{readp}}, \code{\link{fitdistribution}},
#'\code{\link{readl}}, \code{\link{constructplant}}
#'@keywords misc
#'@examples
#'
#'
#'\dontrun{
#'
#'# Replace angles by sampling from an ellipsoidal distribution:
#'mydist <- angledist("ellipsoid", distpars=0.7)
#'replaceangles(pfile="someplant.p", distobj=mydist)
#'
#'# Make constant angles:
#'replaceangles(pfile="someplant.p", newangles=45)
#'
#'# Change new file name:
#'replaceangles(pfile="someplant.p", outputfile="someplant 45degrees.p", newangles=45)
#'
#'# Change petiole orientation, choose pfile with dialog box:
#'replaceangles("Az.2", newangles=runif(300, 0, 360))
#'
#'# Modify various variables in a pfile, until we end up with an artificial plant,
#'# useful for testing.
#'# Order of changes, in this case (although it does not matter!):
#'# Leaf azimuth, leaf orientation, leaf angle, petiole length, petiole angle.
#'modifypfile("originalplant.p", whichvar="Az.3", outputfile="testplant.p", newvalues=45)
#'modifypfile("testplant.p", whichvar="Or", outputfile="testplant.p", newvalues=45)
#'modifypfile("testplant.p", whichvar="An.3", outputfile="testplant.p", newvalues=-45)
#'modifypfile("testplant.p", whichvar="L.2", outputfile="testplant.p", newvalues=10)
#'modifypfile("testplant.p", whichvar="An.2", outputfile="testplant.p", newvalues=45)
#'
#'
#'}
#'
#'
replaceangles <- function(whichangle="An.3",
						  pfile=NA,
						  outputfile = NA,
					    newangles = NULL,
						  distobj = NULL
							){

	if(is.na(pfile)){
		if(.Platform$OS.type != "windows" || !interactive())
			stop("Please provide a plant (.P) file")
		pfile <- file.choose()
	}						  
	
	if(is.na(outputfile)){
		filer <- gsub("\\.P$","",pfile,ignore.case=TRUE)
		outputfile <- paste0(filer, "-Modified.p")
	}
	
	prestuff <- readLines(pfile, warn=FALSE)[1:5]

	# Read angles
	pdata <- readp(pfile)
	angles <- getangles(pfile, whichangle)
	N <- length(angles)
	
	# Sample from distribution object, or use provided angles.
	if(!is.null(distobj)){
		newangles <- drawsample(distobj, N, degrees=TRUE)
	}
	if(!is.null(newangles)){
		if(length(newangles) < N){
			if(length(newangles)==1)newangles <- rep(newangles, N)
			if(length(newangles)>1)newangles <- sample(newangles, N, replace=TRUE)
		} else
			newangles <- newangles[1:N]
	}
	if(is.null(newangles) & is.null(distobj))stop("No new angles specified\n")
	
	# Replace angles with new angles.
	pdata[pdata$Lt >= 1, whichangle] <- round(newangles,3)
	
	# Write the data.
	writeLines(prestuff, outputfile)
	
	write.table(pdata, outputfile, append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)

}

#' @rdname ModifyPfiles
modifypfile <- function(  pfile=NA,
                          whichvar=NA,
                          outputfile = "tmp.p",
                          newvalues = NULL
){
  
  
  if(is.na(pfile)){
    if(.Platform$OS.type != "windows" || !interactive())
      stop("Please provide a plant (.P) file")
    pfile <- file.choose()
  }					
  
  if(is.na(whichvar))stop("Need a variable name (set 'whichvar' argument).")
  
  prestuff <- readLines(pfile, warn=FALSE)[1:5]
  
  # Read old values
  pdata <- readp(pfile)
  oldvals <- pdata[,whichvar]
  N <- length(oldvals)
  
  # Sample from distribution object, or use provided angles.
  if(!is.null(newvalues)){
    if(length(newvalues) < N)
      if(length(newvalues)==1)newvalues <- rep(newvalues, N)
    if(length(newvalues)>1)newvalues <- sample(newvalues, N, replace=TRUE)
    else
      newvalues <- newvalues[1:N]
  }
  
  # Replace values.
  pdata[, whichvar] <- newvalues
  
  # Write the data.
  writeLines(prestuff, outputfile)
  
  write.table(pdata, outputfile, append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
  
}

#' @rdname ModifyPfiles
changeinternodes <- function(pfile=NA,
                             outputfile = "tmp.p",
                             method=c("perc","constant"),
                             changeperc=50,
                             consvalue = NA
){
  
  if(is.na(pfile)){
    if(.Platform$OS.type != "windows" || !interactive())
      stop("Please provide a plant (.P) file")
    pfile <- file.choose()
  }						 
  
  method <- match.arg(method)
  
  prestuff <- readLines(pfile, warn=FALSE)[1:5]
  pdata <- readp(pfile)
  
  # Internode lengths are both branch and stem lengths:
  branchlens <- pdata$L.1[pdata$L.1 > 0]
  stemlens <- pdata$L[pdata$L > 0.1]  # 0.1 is used for dummy stem segments.
  
  if(method == "perc"){
    branchlens <- branchlens * changeperc/100
    stemlens <- stemlens * changeperc/100
  }
  if(method == "constant"){
    branchlens <- consvalue
    stemlens <- consvalue
  }
  
  # replace.
  pdata$L.1[pdata$L.1 > 0] <- round(branchlens,4)
  pdata$L[pdata$L > 0.1] <- round(stemlens,4)
  
  # Write the data.
  writeLines(prestuff, outputfile)
  
  write.table(pdata, outputfile, append=TRUE, sep=" ", row.names=FALSE, col.names=FALSE)
  
}
