RnDif.hom <-
function (lx, ly, nx, ny, e, m, rho, bdc_top, rn_lam, rn_ema, rn_dif, rn_sol, solution, ...) 
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

xgrid		=	setup.grid.1D(N=Nx,L=Lx)
ygrid		=	setup.grid.1D(N=Ny,L=-Ly)
xygrid	=	setup.grid.2D(x.grid=xgrid,y.grid=ygrid)

x.axis 	 	= X
y.axis.conc 	= Y
y.axis.flux 	= -ygrid$x.int

poro		=	e			
hum		=	m
rho_g		=	rho				
Bc		=	bdc_top				
lambda	=	rn_lam			
E		=	rn_ema				
D0		=	rn_dif				
L		=	rn_sol				
	
beta		=	(1-hum+L*hum)*poro			
De		=	D0*poro*exp(-6*hum*poro-6*hum^(14*poro)) 	
D		=	beta*De
G		=	lambda*E*((1-poro)/poro)*rho_g

	
Dgrid=setup.prop.2D(value=D,grid=xygrid,y.value=D)

VFgrid=setup.prop.2D(value=poro,grid=xygrid,y.value=poro)

Agrid=setup.prop.2D(value=1,grid=xygrid,y.value=1)


Diff2D = function(t,y,parms){
		
CONC=matrix(nrow=Nx,ncol=Ny,y)

	   
Tran=tran.2D(CONC,C.x.up=0,C.x.down=0,C.y.up=Bc,
flux.x.up=0,flux.x.down=0,a.bl.x.up=0,a.bl.x.down=0,
D.grid=Dgrid,VF.grid=VFgrid,grid=xygrid,A.grid=Agrid,
full.output=TRUE)
		

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


RnDif.hom=list()
RnDif.hom$conc=mat
RnDif.hom$flux=mat2

	
if (solution == "steady") {

	return(list(x.axis=x.axis,y.axis.conc=y.axis.conc,
			y.axis.flux=y.axis.flux,conc=RnDif.hom$conc,
			flux=RnDif.hom$flux))
	}



}

