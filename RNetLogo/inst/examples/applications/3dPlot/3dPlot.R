library(RNetLogo)

# adapt the path
nl.path <- "C:/Program Files/NetLogo 5.3/app"
NLStart(nl.path, gui=FALSE)
model.path <- "models/Curricular Models/Urban Suite/Urban Suite - Sprawl Effect.nlogo"
NLLoadModel(paste(nl.path,model.path,sep="/"))
NLCommand("random-seed 32453433")
NLCommand("resize-world -20 20 -20 20 ")
NLCommand("set smoothness 10",
          "set max-attraction 5",
          "set population 500",
          "set seeker-search-angle 200",
          "set seeker-patience 15",
          "set wait-between-seeking 5")

NLCommand("setup")
NLDoCommand(150,"go")
attraction <- NLGetPatches("attraction",as.matrix=T)
pxcor <- NLReport(c("min-pxcor","max-pxcor"))
pycor <- NLReport(c("min-pycor","max-pycor"))
d <- list(x=seq(pxcor[[1]],pxcor[[2]]),
     y=seq(pycor[[1]],pycor[[2]]),
     z=attraction)



# adapted from http://addictedtor.free.fr/graphiques/RGraphGallery.php?graph=1
# author: Romain Francois
kde2dplot <- function(d,                # a 2d density computed by kde2D
                      ncol=50,          # the number of colors to use
                      zlim=c(0,max(z)), # limits in z coordinates
                      nlevels=20,       # see option nlevels in contour
                      theta=30,         # see option theta in persp
		                  phi=30)           # see option phi in persp
{
  z   <- d$z
  nrz <- nrow(z)
  ncz <- ncol(z)
  couleurs  <- tail(topo.colors(trunc(1.4 * ncol)),ncol)
  fcol      <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
  dim(fcol) <- c(nrz,ncz)
  fcol      <- fcol[-nrz,-ncz]
  
  par(mfrow=c(1,2),mar=c(0.5,0.5,0.5,0.5))
  persp(d,col=fcol,zlim=zlim,theta=theta,phi=phi,
        zlab="attraction",xlab="x",ylab="y")
  
  par(mar=c(2,2,2,2))
  image(d,col=couleurs)
  contour(d,add=T,nlevels=nlevels)
  box()
}
kde2dplot(d)

NLQuit()