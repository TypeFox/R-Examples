#' Plot the stand in 3D
#' 
#' @description Reads the MAESTRA trees file, and plots the stand in 3D. Supports all
#' MAESTRA crown shapes except the box shape.  Looks for the 'trees.dat' file
#' in the current working directory, unless specified (see Examples).  The XY
#' coordinates *must be present* in the 'trees.dat' file. Users will typically
#' only use the 'Plotstand' function.
#' 
#' Optionally reads the crown shape from the 'str.dat' file, and plots the
#' correct crown shape for each species in the stand by reading the
#' multi-species namelists in 'confile.dat' and 'trees.dat'.
#' 
#' The target trees are colored red (unless specified otherwise, see Details),
#' if the 'itargets' is specified in the confile.
#' 
#' Attempts to read indivradx, indivrady, indivhtcrown, indivdiam, and
#' indivhttrunk namelists from the 'trees.dat' file. If any of these fail, the
#' 'all' versions are tried ('allradx', etc.).  Although MAESTRA runs fine when
#' no XY coordinates are provided, this plot function crashes. A future
#' implementation will calculate XY coordinates in the same way as MAESTRA.
#' 
#' If the 'strfiles' parameter is set in 'confile.dat' (one str.dat file for
#' each species in the stand), these files are opened and used to set the crown
#' shape by species. Alternatively, you may specify crownshape as a parameter,
#' and override reading of str.dat by setting readstrfiles=FALSE.
#' 
#' The 'nz' and 'nalpha' arguments specify the 'smoothness' of the crowns:
#' higher values provide more detailed triangulation of the crowns, at the
#' expense of speed.
#' 
#' @details For large stands, the plot takes quite a while to complete. This code is far 
#' from optimized for speed, because I am patient. Also, minimize the rgl window to greatly speed it up.
#' 
#' @param treesfile By default, the 'trees.dat' file in the current dir.
#' @param strfile Not used, yet.
#' @param crownshape Character,
#' "cone","elipsoid","ellipsoid","halfellipsoid","paraboloid","cylinder", or
#' abbreviation.
#' @param addNarrow Logical. Add arrow pointing North?
#' @param xyaxes Logical. Add annotated X and Y axes?
#' @param labcex Relative size of X and Y axis labels.
#' @param axiscex Relative size of X and Y axis annotation.
#' @param verbose If TRUE, writes more info to the screen while plotting.
#' @param CL Crown length (m).
#' @param CW Crown width (m).
#' @param readstrfiles Read the 'str.dat' file(s) to find out crown shape?
#' @param targethighlight Plot the target trees in red?
#' @param crowncolor The color of the tree crowns. Default, obviously,
#' 'forestgreen'.
#' @param stemcolor The color of the tree stems. Default 'brown'.
#' @param HCB Height of crown base (m).
#' @param X,Y X- and Y-coordinates of tree stem base (m).
#' @param dbh Stem diameter (m). Converted to m if appears to be in cm.
#' @param nz Number of z divisions (increase number to get smoother crowns).
#' @param nalpha Number of angular divisions (increase number to get smoother
#' crowns).
#' @param idate If multiple dates are provided for tree size variables, which
#' one to display.
#' @param path The folder where the input files are stored.
#' @param \dots Further parameters passed (to plottree, or triangles3d).
#' @return An rgl device is opened.
#' @author Remko Duursma
#' @examples
#' 
#' 
#' \dontrun{
#' 
#' # Plot the 'trees.dat' file in the current working directory:
#' Plotstand()
#' 
#' # Open a dialog box to select a trees.dat file:
#' Plotstand(file.choose())
#' 
#' # Save a snapshot to a .png file.
#' # Note: make sure to move the 3D plot into view 
#' # (so that other windows are not blocking it!)
#' snapshot3d('myforest.png')
#' 
#' # For publication-quality graphs:
#' Plotstand(nz=50, nalpha=50)
#' 
#' 
#' }
#' 
#' @export
#' @rdname Plotstand
#' @importFrom rgl par3d
#' @importFrom rgl title3d
#' @importFrom rgl axes3d
#' @importFrom rgl text3d
Plotstand <- function(treesfile="trees.dat", 
          				    strfile="str.dat",
          					  crownshape=c("cone","ellipsoid","round","halfellipsoid","paraboloid","cylinder"), 
          					  readstrfiles=TRUE,
          					  targethighlight=TRUE,
          					  addNarrow=TRUE, 
                      xyaxes=TRUE,
                      labcex=1,
                      axiscex=1,
                      verbose=FALSE,
                      idate=1,
                      path="",
                      ...){

  o <- getwd()
  on.exit(setwd(o))
  if(path != "")setwd(path)
  
	notrees <- readPAR(treesfile, "notrees", "plot")
	
	crownshapes <- rep(NA,notrees)
	
	haveconfile <- file.exists("confile.dat")
	
	if(targethighlight & !haveconfile){
		warning("No confile.dat found - target trees not highlighted")
		targethighlight <- FALSE
		crowncolors <- rep("forestgreen",notrees)
	}
	
	if(targethighlight){
		crowncolors <- rep("forestgreen",notrees)
		itargets <- readPAR("confile.dat", "itargets", "treescon", fail=FALSE)
		if(all(is.na(itargets))){
			warning("itargets not read in confile.dat")
		} else crowncolors[itargets] <- "red"
	}
	
	if(!haveconfile){
		strfiles <- strfile
	} else {
		strfiles <- readPAR("confile.dat","strfiles","species",fail=FALSE)
		if(all(is.na(strfiles)))
			strfiles <- strfile
		else {
			strfiles <- gsub("'","",strfiles)
		}
		if(!all(file.exists(strfiles)))
			stop("Not all strfiles are in the current working directory.")
	}
	
	if(readstrfiles){
		species <- readPAR("trees.dat","ispecies","speclist",fail=FALSE)
		if(all(is.na(species)))
      species <- rep(1,notrees)
    
		for(i in 1:notrees){
			crownshapes[i] <- tolower(readPAR(strfiles[species[i]],"cshape","canopy"))
	  }
		crownshapes <- gsub("'","",crownshapes)
	
	} else crownshapes[] <- match.arg(crownshape)  
	

  xycoor <- readPAR(treesfile, "xycoords", "xy", fail=FALSE)
  if(all(is.na(xycoor)))stop("XY coordinates must be present in trees.dat file")
  xycoor <- matrix(xycoor,ncol=2,byrow=TRUE)
  X <- xycoor[,1]
  Y <- xycoor[,2]
  
  Bearing <- readPAR(treesfile, "bearing", "plot")
  
  # varname w/o 'indiv' or 'all'!
  readVar <- function(varname, idate=1){
    
    indvar <- paste0("indiv",varname)
    allvar <- paste0("all",varname)
    
    vals <- readPAR(treesfile, "values", indvar,fail=FALSE)
    whichvar <- indvar
    
    if(all(is.na(vals))){
      vals <- rep(readPAR(treesfile, "values", allvar)[1], notrees)
      whichvar <- allvar
    }
      
    # dates?
    ndates <- readPAR(treesfile, "nodates", whichvar)
    if(all(is.na(ndates)) || ndates==1)return(vals)
    
    # else...
    vals <- matrix(vals, ncol=ndates, byrow=TRUE)
    return(vals[,idate])
  }
  
  radx <- readVar("radx", idate) 
  CW <- 2*radx
     
  CL <- readVar("htcrown", idate)
	DBH <- readVar("diam", idate)
	if(max(DBH) > 3)DBH <- 0.01 * DBH
  
  # Avoid small DBH, rgl does not like that.
  DBH <- pmin(0.05, DBH)
  
  HCB <- readVar("httrunk", idate)
  
  if(any(is.na(c(CL,CW,X,Y,HCB,DBH))))stop("Missing values somewhere!")
  
  Openstand(treesfile)

  for(i in 1:notrees){
	if(verbose)message("Plotting tree number : ", i)
      plottree(crownshape=crownshapes[i], CL=CL[i], CW=CW[i], 
          HCB=HCB[i], X=X[i], Y=Y[i], dbh=DBH[i], crowncolor=crowncolors[i])
  }
     
	if(addNarrow){
		X0 <- readPAR(treesfile,"x0","plot", fail=FALSE)
		if(is.na(X0))X0 <- 0
		Y0 <- readPAR(treesfile,"y0","plot",fail=FALSE)
		if(is.na(Y0))Y0 <- 0
		
		Xmax <- readPAR(treesfile,"xmax","plot")
		Ymax <- readPAR(treesfile,"ymax","plot")
		addarrow(x0=X0 + 0.1*Xmax,y0=Y0 + 0.1*Ymax,
			 len=0.1*(Ymax-X0),bearing=Bearing)
  }
	if(xyaxes){
		par3d(cex=axiscex)
	  axes3d(c('x-','y-'))
	  par3d(cex=labcex)
	  title3d(xlab="X",ylab="Y")
	}
   
}
