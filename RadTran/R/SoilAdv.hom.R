SoilAdv.hom <-
function (lx, ly, nx, ny, e, k_soil, miu, dp_bot, dp_top,solution, ...)
{ 


Lx	=	lx			
Ly	=	-ly				
Nx	=	nx			
Ny	=	ny				
dx	=	Lx/Nx				
dy	=	Ly/Ny				
X	=	seq(dx/2,by=dx,len=Nx)		
Y	=	seq(dy/2,by=dy,len=Ny)		
A	=	Lx*(-Ly)				

xgrid=setup.grid.1D(N=Nx,L=Lx)
ygrid=setup.grid.1D(N=Ny,L=-Ly)
xygrid=setup.grid.2D(x.grid=xgrid,y.grid=ygrid)

x.axis 	 	= X
y.axis.press 	= Y
y.axis.flux 	= -ygrid$x.int

poro		=	e
k		=	k_soil
vis		=	miu
dp_bot2	=	dp_bot
dp_top2	=	dp_top
D		=	k/vis
			
Dgrid=setup.prop.2D(value=D,grid=xygrid,y.value=D)
VFgrid=setup.prop.2D(value=poro,grid=xygrid,y.value=poro)
Agrid=setup.prop.2D(value=1,grid=xygrid,y.value=1)

Diff2D=function(t,y,parms){

PRESS=matrix(nrow=Nx,ncol=Ny,y)

Tran=tran.2D(PRESS,C.x.up=0,C.x.down=0,C.y.up=-dp_top2,
	C.y.down=dp_bot2,flux.x.up=0,flux.x.down=0,
	a.bl.x.up=0,a.bl.x.down=0,D.grid=Dgrid,
	grid=xygrid,A.grid=Agrid,full.output=TRUE)
	   
	
dPRESS=Tran$dC
xFlux=Tran$x.flux		
yFlux=Tran$y.flux		
yTopFlux=Tran$flux.y.up		
yBottFlux=Tran$flux.y.down	

return(list(as.vector(dPRESS),yFlux=yFlux,xFlux=xFlux,
	yTopFlux=yTopFlux,yBottFlux=yBottFlux))

}

y=runif(Nx*Ny) #condição inicial
std2=steady.2D(func=Diff2D,y=as.vector(y),time=0,
	positive=TRUE,parms=NULL,lrw=9e7,dimens=c(Nx,Ny))

mat=matrix(nrow=Nx,ncol=Ny,-std2$y) 

mat2=matrix(nrow=Nx,ncol=Ny+1,std2$yFlux)


	SoilAdv.hom=list()
	SoilAdv.hom$conc=mat
	SoilAdv.hom$flux=mat2


if (solution == "steady") {

		return(list(x.axis=x.axis,y.axis.press=y.axis.press,
			y.axis.flux=y.axis.flux,
			press=SoilAdv.hom$conc,flux=SoilAdv.hom$flux))

		}

}
