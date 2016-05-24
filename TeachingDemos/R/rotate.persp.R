"rotate.persp" <-
function(x,y,z){

  if(!requireNamespace('tcltk', quietly = TRUE)){stop('The tcltk package is needed')}

  persp.refresh <- function(...){
    persp(x,y,z, theta=slider(no=1), phi=slider(no=2),
          r=slider(no=3), d=slider(no=4),
          ltheta=slider(no=5), lphi=slider(no=6),
          shade=slider(no=7),col='lightblue'
          )
  }

  slider( persp.refresh,
	c('theta','phi','r','d','ltheta','lphi','shade'),
	c(-360,-180,0,0,0,0,0),
	c(360,180, 10, 5, 360, 180, 1),
	c(5,5,.25,.1,5,5,.05),
	c(0,15,sqrt(3),1,120,15,.7),
	'PerspControl')

}

