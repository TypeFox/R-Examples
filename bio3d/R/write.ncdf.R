`write.ncdf` <-
function(x, trjfile="R.ncdf", cell=NULL){

  ##- Load ncdf4 package
  oops <- requireNamespace("ncdf4", quietly = TRUE)
  if(!oops)
    stop("Please install the ncdf4 package from CRAN")

  ## Error checking
  if(!is.matrix(x))
    stop("input x should be a natom by nframe numeric matrix of coordinates")
  
  if(!is.null(cell)) {
     if(!is.matrix(cell))
       stop("input cell should be a 6 by nframe numeric matrix")
  }

  nframe <- nrow(x)
  natom <- ncol(x)/3

  ## Define dimensions
  frame <- ncdf4::ncdim_def(name="frame", units="", vals=c(1:nframe),
                        unlim=TRUE, create_dimvar=FALSE)

  spatial <- ncdf4::ncdim_def(name="spatial", units="", vals=1:3, #c(1:3),#"xyz",
                          unlim=FALSE, create_dimvar=TRUE)

  atom <- ncdf4::ncdim_def(name="atom", units="", vals=c(1:natom),
                       unlim=FALSE, create_dimvar=FALSE)

  if(!is.null(cell)) {
     label <- ncdf4::ncdim_def(name="label", units="", vals=1:5, 
             unlim=FALSE, create_dimvar=FALSE)
     cell_spatial <- ncdf4::ncdim_def(name="cell_spatial", units="", 
            vals=1:3, unlim=FALSE, create_dimvar=TRUE)
     cell_angular <- ncdf4::ncdim_def(name="cell_angular", units="", 
            vals=1:3, unlim=FALSE, create_dimvar=TRUE)
  }
     

##  label <- dim.def.ncdf(name="label", units="", vals=1:5, ##???
##                        unlim=FALSE, create_dimvar=FALSE)

##  cells <- dim.def.ncdf(name="cell_spatial", units="", vals=1:3,# "abc",
##                        unlim=FALSE, create_dimvar=TRUE)

##  cella <- dim.def.ncdf(name="cell_angular", units="", vals=1:3,
##                        #vals=c("alpha", "beta", "gamma"),
##                        unlim=FALSE, create_dimvar=TRUE)

  ## Define variables
  time <- ncdf4::ncvar_def(name="time", units="picosecond", dim=frame,
                       missval=1e+30, prec="float") #"single" float
  coor <- ncdf4::ncvar_def(name="coordinates", units="angstrom", missval=1e+30,
                       dim=list(spatial,atom,frame), prec="float")#"single" float
  if(!is.null(cell)) {
    cell_lengths <- ncdf4::ncvar_def(name="cell_lengths", units="angstrom", 
                 missval=1e+30, dim=list(cell_spatial, frame), prec="double")
    cell_angles <- ncdf4::ncvar_def(name="cell_angles", units="degree",
                 missval=1e+30, dim=list(cell_angular, frame), prec="double")
  }
##  cell.len <- var.def.ncdf(name="cell_lengths", units="angstrom", missval=1e+30,
##                           dim=list(cells,frame),  prec="double")
##  cell.ang <- var.def.ncdf(name="cell_angles", units="degrees", missval=1e+30,
##                           dim=list(cella,frame),  prec="double")
  
  ## Create the file
  if(!is.null(cell)) {
     ncw <- ncdf4::nc_create( trjfile, list(time, coor, cell_lengths, cell_angles))
  } else {
     ncw <- ncdf4::nc_create( trjfile, list(time, coor))#, cell.len, cell.ang) )
  }

  ## Write data to file
  if(is.null(rownames(x)))
     ncdf4::ncvar_put(ncw, time, c(1:nframe), start=1, count=nframe)
  else
     ncdf4::ncvar_put(ncw, time, as.numeric(rownames(x)), start=1, count=nframe)
  ncdf4::ncvar_put( ncw, coor, t(x), start=c(1,1,1), count=c(3,natom,nframe))
  if(!is.null(cell)) {
     ncdf4::ncvar_put( ncw, cell_lengths, t(cell[,1:3]), start=c(1,1), count=c(3,nframe))
     ncdf4::ncvar_put( ncw, cell_angles, t(cell[,4:6]), start=c(1,1), count=c(3,nframe))
  }

  ## Define Required Attributes
  ncdf4::ncatt_put(ncw, varid=0, attname="Conventions", attval="AMBER")
  ncdf4::ncatt_put(ncw, varid=0, attname="ConventionVersion", attval="1.0")
  ncdf4::ncatt_put(ncw, varid=0, attname="program",attval="bio3d")
  ncdf4::ncatt_put(ncw, varid=0, attname="programVersion", attval="1.2")
  null <- ncdf4::nc_close(ncw)
}

