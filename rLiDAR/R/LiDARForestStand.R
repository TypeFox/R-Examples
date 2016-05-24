#'3D stand visualization of LiDAR-derived individual trees
#'
#'@description Draws a 3D scatterplot for individual trees detected from LiDAR data.   
#'
#'@usage LiDARForestStand(crownshape = c("cone", "ellipsoid", "halfellipsoid",
#'                 "paraboloid", "cylinder"), CL = 4, CW = 8, HCB = 10, 
#'                   X = 0, Y = 0, dbh = 0.3, crowncolor = "forestgreen", 
#'                  stemcolor = "chocolate4", resolution="high", mesh=TRUE) 
#'
#'@param crownshape shape of individual tree crown: "cone", "ellipsoid","halfellipsoid", "paraboloid" or "cylinder". Default is "halfellipsoid".
#'@param CL crown length.
#'@param CW crown diameter.
#'@param HCB height at canopy base.
#'@param X x-coordinate.
#'@param Y y-coordinate.
#'@param dbh diameter at breast height (1.73 m).
#'@param crowncolor crown color.
#'@param stemcolor stem color.
#'@param resolution crown resolution: "low", "medium" and "high".
#'@param mesh Logical, if TRUE (default) returns a tree crown mesh model, and if FALSE returns a tree crown line mode.  
#'@return Returns a 3-D scatterplot of the individual trees as identified automatically from the LiDAR.   
#'@author Carlos Alberto Silva and Remko Duursma. Uses code by Remko Duursma (\emph{Maeswrap} package,see "Plotstand").
#'@references \url{http://maespa.github.io/}
#'@examples
#'
#'\dontrun{
#'#=======================================================================#
#'# EXAMPLE 01: Plotting single trees
#'#=======================================================================#
#'
#'# cone crown shape
#'library(rgl)
#'open3d() 
#'LiDARForestStand(crownshape = "cone", CL = 10, CW =7, 
#'            HCB = 5, X =0, Y = 0, dbh = 0.4, crowncolor = "forestgreen", 
#'            stemcolor = "chocolate4", resolution="high", mesh=TRUE) 
#'                        
#'# ellipsoid crown shape 
#'open3d()
#'LiDARForestStand(crownshape = "ellipsoid", CL = 10, CW =7, 
#'            HCB = 5, X =0, Y = 0, dbh = 0.4, crowncolor = "forestgreen", 
#'            stemcolor = "chocolate4", resolution="high", mesh=TRUE) 
#'                        
#'# halfellipsoid crown shape 
#'open3d()
#'LiDARForestStand(crownshape = "halfellipsoid", CL = 10, CW =7, 
#'            HCB = 5, X =0, Y = 0, dbh = 0.4, crowncolor = "forestgreen", 
#'            stemcolor = "chocolate4", resolution="high", mesh=TRUE) 
#'                        
#'# paraboloid crown shape 
#'open3d()
#'LiDARForestStand(crownshape = "paraboloid", CL = 10, CW =7, 
#'            HCB = 5, X =0, Y = 0, dbh = 0.4, crowncolor = "forestgreen", 
#'            stemcolor = "chocolate4", resolution="high", mesh=TRUE)
#'
#'# cylinder crown shape 
#'open3d()
#'LiDARForestStand(crownshape = "cylinder", CL = 10, CW =7, 
#'            HCB = 5, X =0, Y = 0, dbh = 0.4, crowncolor = "forestgreen", 
#'            stemcolor = "chocolate4", resolution="high", mesh=TRUE)
#'                        
#'# Set the shape=FALSE 
#'open3d()
#'LiDARForestStand(crownshape = "paraboloid", CL = 10, CW =7, 
#'            HCB = 5, X =0, Y = 0, dbh = 0.4, crowncolor = "forestgreen", 
#'            stemcolor = "chocolate4", resolution="high", mesh=FALSE)
#' 
#'#=======================================================================#                                           
#'#EXAMPLE 02: Plotting a forest plantation stand in virtual 3-D space
#'#=======================================================================#
#' 
#'# Set the dimensions of the displayed forest stand
#'xlength<-30 # x length
#'ylength<-20 # y length
#'
#'# Set the space between trees
#'sx<-3 # x space length
#'sy<-2 # y space length
#'
#'# Tree location grid
#'XYgrid <- expand.grid(x = seq(1,xlength,sx), y = seq(1,ylength,sy))
#'
#'# Get the number of trees
#'Ntrees<-nrow(XYgrid)
#'
#'# Plot a virtual Eucalyptus forest plantation stand using the halfellipsoid tree crown shape
#'
#'# Set stand trees parameters
#'meanHCB<-5  # mean of the height at canopy base
#'sdHCB<-0.1  # standard deviation of the height at canopy base
#'HCB<-rnorm(Ntrees, mean=meanHCB, sd=sdHCB) # height at canopy base
#'CL<-HCB     # tree crown height
#'CW<-HCB*0.6 # tree crown diameter
#'
#'open3d()    # open a rgl window
#'
#'# Plotting the stand
#'for( i in 1:Ntrees){
#'  LiDARForestStand(crownshape = "halfellipsoid", CL = CL[i], CW = CW[i], 
#'              HCB = HCB[i], X = XYgrid[i,1], Y = XYgrid[i,2], dbh = 0.4, 
#'              crowncolor = "forestgreen", stemcolor = "chocolate4", 
#'              resolution="high", mesh=TRUE) 
#'}
#'                            
#'# Add other plot parameters
#'axes3d(c("x-", "x-", "y-", "z"), col="gray")       # axes
#'title3d(xlab = "X Coord", ylab = " Y Coord", zlab = "Height", col="red") # title
#'planes3d(0, 0, -1, 0.001, col="gray", alpha=0.7)   # set a terrain plane
#'
#'
#'# Plotting a virtual single-species forest plantation stand using "cone" tree crown shape
#'
#'# Set parameters f trees growing within the virtual stand
#'meanHCB<-3  # mean of the height at canopy base
#'sdHCB<-0.1  # standard deviation of the height at canopy base
#'HCB<-rnorm(Ntrees, mean=meanHCB, sd=sdHCB) # height at canopy base
#'CL<-HCB*2.0 # tree crown height
#'CW<-HCB*1.3 # tree crown diameter
#'
#'open3d() # open a rgl window
#'# Plot stand
#'for( i in 1:Ntrees){
#'  LiDARForestStand(crownshape = "cone", CL = CL[i], CW = CW[i], 
#'              HCB = HCB[i], X = XYgrid[i,1], Y = XYgrid[i,2], dbh = 0.4, 
#'              crowncolor = "forestgreen", stemcolor = "chocolate4", 
#'              resolution="high", mesh=TRUE) 
#'}
#'                            
#'# Add other plot parameters
#'axes3d(c("x-", "x-", "y-", "z"), col="gray")       # axes
#'title3d(xlab = "X Coord", ylab = " Y Coord", zlab = "Height", col="red") # title
#'planes3d(0, 0, -1, 0.001, col="gray", alpha=0.7)   # set a terrain plane
#'
#'#=======================================================================#
#'# EXAMPLE 03: Plotting a virtual mixed forest stand
#'#=======================================================================#
#'
#'# 01. Plot different trees species in the stand with different crown shapes 
#'
#'# Set the number of trees
#'Ntrees<-80 
#'
#'# Set the trees locations
#'xcoord<-sample(1:100, Ntrees)  # x coord
#'ycoord<-sample(1:100, Ntrees)  # y coord
#'
#'# Set a location grid of trees 
#'XYgrid<-cbind(xcoord,ycoord)
#'
#'# Plot the location of the trees
#'plot(XYgrid, main="Tree location")
#'
#'meanHCB<-7 # mean of the height at canopy base
#'sdHCB<-3   # standard deviation of height at canopy base
#'HCB<-rnorm(Ntrees, mean=meanHCB, sd=sdHCB) # height at canopy base
#'crownshape<-sample(c("cone", "ellipsoid","halfellipsoid", 
#'                   "paraboloid"), Ntrees, replace=TRUE) # tree crown shape 
#'CL<-HCB*1.3 # tree crown height
#'CW<-HCB     # tree crown diameter
#'
#'open3d() # open a rgl window
#'# Plot stand
#'
#'for( i in 1:Ntrees){
#'  LiDARForestStand(crownshape = crownshape[i], CL = CL[i], CW = CW[i], 
#'              HCB = HCB[i], X = as.numeric(XYgrid[i,1]), Y = as.numeric(XYgrid[i,2]), 
#'              dbh = 0.4, crowncolor = "forestgreen", stemcolor = "chocolate4", 
#'              resolution="high", mesh=TRUE)
#'}
#'                          
#'# Add other plot parameters
#'axes3d(c("x-", "x-", "y-", "z"), col="gray")       # axes
#'title3d(xlab = "X Coord", ylab = " Y Coord", zlab = "Height", col="red") # title
#'planes3d(0, 0, -1, 0.001, col="gray", alpha=0.7)   # set a terrain plane
#'
#'
#'# 02. Plot different tree height in the stand using different crown colors
#'
#'# Set the number of trees
#'Ntrees<-80 
#'
#'# Set the tree locations
#'xcoord<-sample(1:100, Ntrees) # x coord
#'ycoord<-sample(1:100, Ntrees) # y coord
#'
#'# Set a location grid of trees 
#'XYgrid<-cbind(xcoord,ycoord)
#'
#'# plot the location of the trees
#'plot(XYgrid, main="Tree location")
#'
#'meanHCB<-7 # mean of the height at canopy base
#'sdHCB<-3   # standard deviation of the height at canopy base
#'HCB<-rnorm(Ntrees, mean=meanHCB, sd=sdHCB) # height at canopy base
#'crownshape<-sample(c("cone", "ellipsoid","halfellipsoid", "paraboloid"), 
#'                   Ntrees, replace=TRUE) # tree crown shape 
#'CL<-HCB*1.3 # tree crown height
#'CW<-HCB     # tree crown diameter
#'
#'# Plot tree height based on the HCB quantiles
#'HCBq<-quantile(HCB) # HCB quantiles
#'crowncolor<-NA*(1:Ntrees) # set an empty crowncolor vector
#'
#'# classify trees by HCB quantile
#'for (i in 1:Ntrees){
#'  if (HCB[i] <= HCBq[2]) {crowncolor[i]<-"red"}                        # group 1
#'  if (HCB[i] > HCBq[2] & HCB[i] <= HCBq[3] ) {crowncolor[i]<-"blue"}   # group 2
#'  if (HCB[i] > HCBq[3] & HCB[i] <= HCBq[4] ) {crowncolor[i]<-"yellow"} # group 3
#'  if (HCB[i] >= HCBq[4]) {crowncolor[i]<-"dark green"}                 # group 4
#'}
#'    
#'open3d() # open a rgl window
#'
#'# Plot stand
#'for(i in 1:Ntrees){  
#'  LiDARForestStand(crownshape = crownshape[i], CL = CL[i], CW = CW[i], 
#'    HCB = HCB[i], X = as.numeric(XYgrid[i,1]), Y = as.numeric(XYgrid[i,2]), 
#'    dbh = 0.4, crowncolor = crowncolor[i],stemcolor = "chocolate4", 
#'    resolution="high", mesh=TRUE) 
#'}
#'    
#'# Add other plot parameters
#'axes3d(c("x-", "x-", "y-", "z"), col="gray")       # axes
#'title3d(xlab = "X Coord", ylab = " Y Coord", zlab = "Height", col="red") # title
#'planes3d(0, 0, -1, 0.001, col="gray", alpha=0.7)   # set a terrain plane
#'
#'}
#'@export
#'@importFrom geometry convhulln
#'@importFrom rgl plot3d open3d bg3d rgl.triangles 
LiDARForestStand<-function(crownshape = c("cone", "ellipsoid","halfellipsoid", "paraboloid", "cylinder"), 
                        CL = 4, CW = 8,HCB = 10, X = 0, Y = 0, dbh = 0.3, crowncolor = "forestgreen", 
                       stemcolor = "chocolate4", resolution="high",mesh=TRUE) 
{
  
  if (crownshape!="cone"& crownshape!="ellipsoid"&crownshape!="halfellipsoid"&crownshape!="paraboloid"&crownshape!="cylinder") {stop("The crownshape parameter is invalid. Please, use one of this crownshape types: 'cone','ellipsoid','halfellipsoid','paraboloid','cylinder'")}
  if (class(HCB)!="numeric") {stop("The HCB parameter is invalid. It is not a numeric parameter")}
  if (class(X)!="numeric") {stop("The X parameter is invalid. It is not a numeric parameter")}
  if (class(Y)!="numeric") {stop("The Y parameter is invalid. It is not a numeric parameter")}
  if (class(dbh)!="numeric") {stop("The HCB parameter is invalid. It is not a numeric parameter")}
  if (resolution!="high" & resolution!="medium" & resolution!="low") {stop("The resolution parameter is invalid. It must to be 'high', 'median' or 'low'")}
  if (class(mesh)!="logical") {stop("The shape parameter is invalid. It must to be a TRUE or FALSE logical statement")}

  
  if (resolution=="low"){nz<-15;nalpha<-15}
  if (resolution=="medium"){nz<-25;nalpha<-25}
  if (resolution=="high"){nz<-40;nalpha<-40}
  
  if (mesh==TRUE) {
    
  shape <- match.arg(crownshape)

  H <- HCB + CL
  dbase <- dbh * (H/(H - 1.3))
  if (!is.finite(dbase)) 
    dbase <- dbh
  
  
  m1 <- coord3dshape(shape, CW = CW, CL = CL, z0 = HCB, x0 = X, 
                                        y0 = Y, nz = nz, nalpha = nalpha)
  m2 <- coord3dshape("cone", CW = dbase, CL = H, z0 = 0, x0 = X, 
                      y0 = Y, nz = nz, nalpha = nalpha)
  
  interpol(m1, col = crowncolor)
  interpol(m2, col = stemcolor)
  
  } else {
    
    TreesModel(crownshape=crownshape, CW = CW, CL = CL, z0 = 0,HCB=HCB, x0 = X, 
                     y0 = Y, nz = nz, nalpha = nalpha, dbh = dbh,crowncolor = crowncolor, 
               stemcolor = stemcolor)
  }
  
  
}

