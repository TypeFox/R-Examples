
#############  this demo shows how to
###########  use RFOC to interactively plot and rotate focal planes

library(tcltk) 
library(tkrplot) 
library(rpanel)
library(RFOC)

#####################  create a set of data points
########  these are sets of azimuth and dip for seismic radiation

ud = matrix(c(  20.2437,   80.3351,
  -47.0462,   38.9554,
  -43.1122,   76.5689,
  226.4546,   74.7748,
  163.0923 ,  85.8983,
  201.9288,   88.5597
), ncol=2, byrow=TRUE)


dd =  matrix(c(
  267.7700,   44.9155,
  233.3922 ,  33.5276,
  222.6540 ,  55.4004,
  190.4177 ,  47.6243,
  224.4884 ,  11.2741,
   97.1261,   73.6107,
  -83.8927,   88.6226
), ncol=2, byrow=TRUE)

ad = matrix(c(
255.8644,   76.2435,
  -83.2491,   78.4293
), ncol=2, byrow=TRUE)



#############################
##########  initially: plot  stereonet and add points on the plot:

MN = makenet()
pnet(MN)

 qp = qpoint(180+ud[,1] , ud[,2] ,  col =2 , pch=1, UP=TRUE, PLOT=TRUE)             

 qp = qpoint(180+dd[,1] , dd[,2] ,  col =4 , pch=1, UP=TRUE, PLOT=TRUE)             

 qp = qpoint(180+ad[,1] , ad[,2] ,  col =1 , pch=1, UP=TRUE, PLOT=TRUE)             

############# define function for animation.
###  this could be changed to adjust for different input ud, dd, ad, etc....


focal.draw <- function(panel)
{
  print(c(panel$az1, panel$dip1, panel$rake1))
  pnet(MN)
  qp = qpoint(180+ud[,1] , ud[,2] ,  col =2 , pch=1, UP=TRUE, PLOT=TRUE)             

 qp = qpoint(180+dd[,1] , dd[,2] ,  col =4 , pch=1, UP=TRUE, PLOT=TRUE)             

 qp = qpoint(180+ad[,1] , ad[,2] ,  col =1 , pch=1, UP=TRUE, PLOT=TRUE)             

  mec = CONVERTSDR(panel$az1, panel$dip1, panel$rake1)
  MEC = MRake(mec$M)
  PlotPlanes(MEC)
  
  panel
}
#############  define control sequence

panel <- rp.control(az1=0, dip1=40, rake1=132)

############## define and stack the slider panels: 
rp.slider(panel, az1, 0, 360,initval=0,   action = focal.draw)
rp.slider(panel, dip1, 0, 90,initval=0,   action = focal.draw)
rp.slider(panel, rake1, 0, 180,initval=0,   action = focal.draw)
