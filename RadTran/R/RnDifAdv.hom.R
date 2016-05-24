RnDifAdv.hom <-
function (lx, ly, nx, ny, e, m,
bdc_top, rn_lam, rn_sol, k_soil, d_bulk, miu, dp, solution, ...)
{ 


Lx	=	lx				
Ly	=	ly			
Nx	=	nx			
Ny	=	ny			
dx	=	Lx/Nx				
dy	=	Ly/Ny				
X	=	seq(dx/2,by=dx,len=Nx)	
Y	=	seq(dy/2,by=dy,len=Ny)		
A	=	Lx*Ly				

xgrid=setup.grid.1D(N=Nx,L=Lx)
ygrid=setup.grid.1D(N=Ny,L=Ly)
xygrid=setup.grid.2D(x.grid=xgrid,y.grid=ygrid)

x.axis 	 	= X
y.axis.conc 	= Y
y.axis.flux 	= -ygrid$x.int

poro		=	e				
hum		=	m				
Bc		=	bdc_top				
lambda	=	rn_lam			
L		=	rn_sol			
k		=	k_soil
D		=	d_bulk					
vis		=	miu				
deltaP	=	dp				

beta	=	(1-hum+L*hum)*poro			
G	=	lambda*10000
v	=	(k/vis)*deltaP/Ly		


Dgrid=setup.prop.2D(value=0,grid=xygrid,y.value=D)
VFgrid=setup.prop.2D(value=poro,grid=xygrid,y.value=poro)
Agrid=setup.prop.2D(value=1,grid=xygrid,y.value=1)
AFDWgrid=setup.prop.2D(value=1,grid=xygrid,y.value=1)
vgrid=setup.prop.2D(value=0,grid=xygrid,y.value=v)


Diff2D=function(t,y,parms){

CONC=matrix(nrow=Nx,ncol=Ny,y)

Tran=tran.2D(CONC,C.x.up=0,C.x.down=0,C.y.up=Bc,
	C.y.down=0,flux.x.up=0,flux.x.down=0,a.bl.x.up=0,
	a.bl.x.down=0,D.grid=Dgrid,v.grid=vgrid,
	AFDW.grid=AFDWgrid,VF.grid=VFgrid,grid=xygrid,
	A.grid=Agrid,full.output=TRUE)
	   

	dCONC=Tran$dC-lambda*beta*CONC+poro*G
	xFlux=Tran$x.flux*A			
	yFlux=Tran$y.flux*A
	yTopFlux=Tran$flux.y.up		
	yBottFlux=Tran$flux.y.down

return(list(as.vector(dCONC),yFlux=yFlux,xFlux=xFlux,
	yTopFlux=yTopFlux,yBottFlux=yBottFlux))

}


y=runif(Nx*Ny) #condição inicial
std2=steady.2D(func=Diff2D,y=as.vector(y),time=0,
	positive=TRUE,parms=NULL,lrw=9e7,dimens=c(Nx,Ny))

mat=matrix(nrow=Nx,ncol=Ny,std2$y)

mat2=matrix(nrow=Nx,ncol=Ny+1,std2$yFlux)

RnDifAdv.hom=list()
RnDifAdv.hom$conc=mat
RnDifAdv.hom$flux=mat2

if (solution == "steady") {

	return(list(x.axis=x.axis,y.axis.conc=y.axis.conc,
			y.axis.flux=y.axis.flux,conc=RnDifAdv.hom$conc,
			flux=RnDifAdv.hom$flux))


	}

}
