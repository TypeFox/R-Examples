RnDif.het <-
function (lx, ly, nx, ny, rho, x.poro, x.hum, y.poro,y.hum, bdc_top, rn_lam, rn_ema, rn_dif, rn_sol, solution, ...)
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

rho_g		=	rho
Bc		=	bdc_top				
lambda	=	rn_lam		
E		=	rn_ema 			
D0		=	rn_dif				
L		=	rn_sol

x.poro2D	=	function(X,Y,a=x.poro)
			return(rep(a,length(X)))


x.hum2D	=	function(X,Y,a=x.hum)
			return(rep(a,length(X)))


x.beta2D	=	function(X,Y,L=rn_sol)
			return(rep((1-x.hum2D(X,Y,a=x.hum)+L*x.hum2D				(X,Y,a=x.hum))*x.poro2D(X,Y,a=x.poro)))

y.beta2D	=	function(X,Y,L=rn_sol)
			return((1-y.hum(X,Y)+L*y.hum(X,Y))*y.poro(X,Y))


x.De2D	=	function(X,Y,D0=rn_dif,a=x.hum,b=x.poro)
			return(rep(D0*b*exp(-6*a*b-6*a^(14*b)),
			length(X)))

y.De2D	=	function(X,Y,D0=rn_dif)
			return(D0*y.poro(X,Y)*exp(-6*y.hum(X,Y)*y.poro				(X,Y)-6*y.hum(X,Y)^(14*y.poro(X,Y))))

x.D2D		=	function(X,Y,a=x.hum,b=x.poro,L=rn_sol,
			D0=rn_dif)
			return(rep(((1-a+L*a)*b)*(D0*b*exp(-6*a*b-6*a^				(14*b))),length(X)))

y.D2D		=	function(X,Y,D0=rn_dif,L=rn_sol)
			return(((1-y.hum(X,Y)+L*y.hum(X,Y))*y.poro(X,Y))			*(D0*y.poro(X,Y)*exp(-6*y.hum(X,Y)*y.poro(X,Y)-6			*y.hum(X,Y)^(14*y.poro(X,Y)))))


x.G2D		=	function(X,Y,lambda=rn_lam,E=rn_ema,rho_g=rho,
			b=x.poro)
			return(rep(lambda*E*((1-b)/b)*rho_g,length(X)))


y.G2D		=	function(X,Y,lambda=rn_lam,E=rn_ema,rho_g=rho)
			return(lambda*E*((1-y.poro(X,Y))/y.poro(X,Y))				*rho_g)




VFgrid	=	setup.prop.2D(func=x.poro2D,grid=xygrid,
			y.func=y.poro)

Hgrid		=	setup.prop.2D(func=x.hum2D,grid=xygrid,
			y.func=y.hum)

Bgrid		=	setup.prop.2D(func=x.hum2D,grid=xygrid,
			y.func=y.beta2D)

Degrid	=	setup.prop.2D(func=x.De2D,grid=xygrid,
			y.func=y.De2D)

Dgrid		=	setup.prop.2D(func=x.D2D,grid=xygrid,
			y.func=y.D2D)

Ggrid		=	setup.prop.2D(func=x.G2D,grid=xygrid,
			y.func=y.G2D)

Agrid		=	setup.prop.2D(value=1,grid=xygrid,y.value=1)



Diff2D=function(t,y,parms){

CONC=matrix(nrow=Nx,ncol=Ny,y)
	   
Tran=tran.2D(CONC,C.x.up=0,C.x.down=0,C.y.up=bdc_top,
flux.x.up=0,flux.x.down=0,a.bl.x.up=0,a.bl.x.down=0,
D.grid=Dgrid,VF.grid=VFgrid,A.grid=Agrid,grid=xygrid,
full.output=TRUE)

dCONC=Tran$dC-lambda*Bgrid$y.mid*CONC+VFgrid$y.mid*Ggrid$y.mid

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


RnDif.het=list()
RnDif.het$conc=mat
RnDif.het$flux=mat2


if (solution == "steady") {

	return(list(x.axis=x.axis,y.axis.conc=y.axis.conc,
			y.axis.flux=y.axis.flux,conc=RnDif.het$conc,
			flux=RnDif.het$flux))
	}



}
