#' @export
#' @rdname Plotstand
plottree <- function(crownshape=c("cone","elipsoid","ellipsoid","round","halfellipsoid","paraboloid","cylinder"), 
	CL=1, CW=1, HCB=1, X=0, Y=0, dbh=0.3, crowncolor="forestgreen", stemcolor="brown",
    nz=25, nalpha=25,
    ...){

    shape <- match.arg(crownshape)
	if(shape == "elipsoid")shape <- "ellipsoid"
	if(shape == "round")shape <- "ellipsoid"
    H <- HCB + CL
    dbase <- dbh * (H / (H - 1.3))
	if(!is.finite(dbase))dbase <- dbh
	
	# tree crown
	m1 <- coord3dshape(shape,CW=CW,CL=CL,z0=HCB,x0=X,y0=Y,nz=nz,nalpha=nalpha)
	
	# stem
	m2 <- coord3dshape("cone",CW=dbase,CL=H,z0=0,x0=X,y0=Y,nz=nz,nalpha=nalpha)
	
	plot3dtriangles(m1, col=crowncolor,...)
	plot3dtriangles(m2, col=stemcolor,...)

}