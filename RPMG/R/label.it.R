`label.it` <-
function(a="", corn=1, ...)
  {	   
##X## put a label in a figure
##X##    ###  corn = corner 1 2 3 4 = Low-Left UP-Left UP-Right Low-Right
    if(missing(corn)) { corn = 1 }
    
	aa=a
	u=par("usr")
        oce=par("cex")
        par(cex=1.2)
	p1em=par('cin')
        shiftx= (u[2]-u[1])*0.02
        shifty=  (u[4]-u[3])*0.02 

	if(corn==1) {
#	text(u[1]+shiftx*p1em[1],u[3]+shifty*p1em[2]   ,labels=aa,  adj=0)
	text(u[1]+shiftx,u[3]+shifty   ,labels=aa, adj=c(0 ,0.5), ...)
	}
 
	if(corn==2) {
#	text(u[2]-shiftx*p1em[1],u[3]+shifty*p1em[2]   ,labels=aa,  adj=1)
	text(u[2]-shiftx,u[3]+shifty   ,labels=aa, adj=c(1, 0), ... )
	}

	if(corn==3) {
#	text(u[2]-shiftx*p1em[1],u[4]- shifty*p1em[2]  ,labels=aa,  adj=1)
        text(u[2]-shiftx,u[4]-2*shifty  ,labels=aa,  adj=c(1, 0.5), ... )
	}
	if(corn==4) {
#	text(u[1]+shiftx*p1em[1],u[4]- shifty*p1em[2] ,labels=aa,  adj=0)
	text(u[1]+shiftx,u[4]-shifty ,labels=aa,  adj=c(0 , 0.5 ), ... )
	}
	par(cex=oce)
}

