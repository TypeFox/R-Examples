
##==============================================================================
## Transport in a two-dimensional finite difference grid
## KS: Made the grid a matrix rather than a vector...
##==============================================================================

tran.2D <- function(C, C.x.up=C[1,], C.x.down=C[nrow(C),],
  C.y.up=C[,1],  C.y.down=C[,ncol(C)],
  flux.x.up=NULL, flux.x.down=NULL, flux.y.up=NULL, flux.y.down=NULL,
  a.bl.x.up=NULL, a.bl.x.down=NULL, a.bl.y.up=NULL, a.bl.y.down=NULL, 
  D.grid=NULL, D.x=NULL, D.y=D.x, v.grid=NULL, v.x=0, v.y=0, 
  AFDW.grid=NULL, AFDW.x=1, AFDW.y=AFDW.x,
  VF.grid=NULL,VF.x=1, VF.y=VF.x,
  A.grid=NULL, A.x=1, A.y=1,
  grid=NULL,dx=NULL, dy=NULL,
  full.check = FALSE, full.output = FALSE)
											
{
  if (is.null(grid))
   if (is.null(dx) | is.null(dy))
      stop("error: either grid or dx and dy should be specified  ")

  Nx <- nrow(C)
  Ny <- ncol(C)

# DEFAULT INFILLING OF GRID PARAMETERS

#==============================================================================
# infilling of 2D numerical grid
#==============================================================================
  if (is.null(grid)) {
     DX    <-  if (is.list(dx)) dx$dx else rep(dx,length.out=Nx)
     DXaux <-  if (is.list(dx)) dx$dx.aux else 0.5*(c(0,rep(dx,length.out=Nx))+
                                  c(rep(dx,length.out=Nx),0))
     DY    <-  if (is.list(dy)) dy$dx else rep(dy,length.out=Ny)
     DYaux <-  if (is.list(dy)) dy$dx.aux else 0.5*(c(0,rep(dy,length.out=Ny))+
                                  c(rep(dy,length.out=Ny),0))

    grid <- list(
      dx    = matrix(nrow=Nx, ncol=Ny, data = DX),  
      dx.aux= matrix(nrow=Nx+1, ncol=Ny, data = DXaux),
      dy    = matrix(nrow=Nx, ncol=Ny, data = DY, byrow=TRUE),  
      dy.aux= matrix(nrow=Nx, ncol=Ny+1, data = DYaux, byrow=TRUE)
                 )
   } else if (class(grid) !="grid.2D") {
   # KS: changed the next statement, -> everywhere dx, dy is a MATRIX
    if (length(grid$dx) != Nx*Ny)         
       grid$dx <- matrix(nrow=Nx, ncol=Ny, grid$dx)
    if (length(grid$dx.aux) != (Nx+1)*Ny) 
       grid$dx.aux <- matrix(nrow=Nx+1, ncol=Ny, grid$dx.aux)
    if (length(grid$dy) != Nx*Ny)         
       grid$dy <- matrix(nrow=Nx, ncol=Ny, data = grid$dy, byrow=TRUE)
    if (length(grid$dy.aux) != Nx*(Ny+1)) 
       grid$dy.aux <- matrix(nrow=Nx, ncol=Ny+1, grid$dy.aux, byrow=TRUE)
   }
   
#==============================================================================
# infilling of grids with x.int, y.int, x.mid, y.mid needed
#==============================================================================
  gridFill <- function(G.x,G.y,Name)    
  {
  
    # check if G.x and G.y is not NULL
    if (is.null(G.x) | is.null(G.y))
      stop( (paste("error: ",Name,"and (",Name,".x and", Name,".y) cannot be NULL at the same time", del="")))
    G.grid <- list()
    # infilling of x-matrix
    if (is.matrix(G.x)) {
      if (sum(abs(dim(G.x) - c(Nx+1,Ny)))!=0)
        stop (paste("error: ",Name,".x matrix not of correct (Nx+1) Ny dimensions", del=""))
      G.grid$x.int <- G.x
      G.grid$x.mid <- 0.5*(G.x[1:Nx,]+G.x[2:(Nx+1),])
    } else if (class(G.x)=="prop.1D") {
      G.grid$x.int <- matrix(data=G.x$int,nrow=(Nx+1),ncol=Ny)
      G.grid$x.mid <- matrix(data=G.x$mid, nrow=Nx, ncol=Ny)
    } else if (length(G.x) == 1) {
      G.grid$x.int <- matrix(data=G.x,nrow=(Nx+1),ncol=Ny)
      G.grid$x.mid <- matrix(data=G.x,nrow=Nx,ncol=Ny)
    } else if (length(G.x) != Nx+1) {
        stop (paste("error: ",Name,".x should be a vector of length 1 or Nx+1", del=""))
    } else {  # correct length
      G.grid$x.int <- matrix(data=G.x,nrow=(Nx+1),ncol=Ny)
      G.grid$x.mid <- matrix(data=0.5*(G.x[1:Nx]  +G.x[2:(Nx+1)]),
                     nrow=Nx, ncol=Ny)
    }
    # infilling of y-matrix
    if (is.matrix(G.y)) {
      if (sum(abs(dim(G.y) - c(Nx,Ny+1)))!=0)
        stop (paste("error: ",Name,".y matrix not of correct Nx(Ny+1)dimensions", del=""))
      G.grid$y.int <- G.y
      G.grid$y.mid <- 0.5*(G.y[,1:Ny]+G.y[,2:(Ny+1)])
    } else if (class(G.y)=="prop.1D") {
      G.grid$y.int <- matrix(data=G.y$int,nrow=Nx,ncol=(Ny+1))
      G.grid$y.mid <- matrix(data=G.y$mid, nrow=Nx, ncol=Ny)
    } else if (length(G.y) == 1) {
      G.grid$y.int <- matrix(data=G.y,nrow=Nx,ncol=(Ny+1))
      G.grid$y.mid <- matrix(data=G.y,nrow=Nx,ncol=Ny)
    } else if (length(G.y) != Ny+1) {
        stop (paste("error: ",Name,".y should be a vector of length 1 or Ny+1", del=""))
    } else {  # correct length
      G.grid$y.int <- matrix(data=G.y,nrow=Nx,ncol=(Ny+1))
      G.grid$y.mid <- matrix(data=0.5*(G.y[1:Nx]  +G.y[2:(Nx+1)]),
                     nrow=Nx, ncol=Ny)
    }
    G.grid
  }

# Need this for VF and A (volume fraction and surface

  if (is.null(VF.grid)) VF.grid <- gridFill(VF.x,VF.y,"VF")
  if (is.null(A.grid))  A.grid  <- gridFill(A.x,A.y,"A")

#==============================================================================
# infilling of other grids with only  x.int and y.int needed
#==============================================================================

  gridInt <- function(G.x,G.y,Name)     # define a function first
  {
    # check if G.x and G.y is not NULL
    if (is.null(G.x) | is.null(G.y))
      stop( (paste("error: ",Name,"and (",Name,".x and", Name,".y) cannot be NULL at the same time", del="")))
    G.grid <- list()
    # infilling of x-matrix
    if (is.matrix(G.x)) {
      if (sum(abs(dim(G.x) - c(Nx+1,Ny)))!=0)
        stop (paste("error: ",Name,".x matrix not of correct (Nx+1) Ny dimensions", del=""))
      G.grid$x.int <- G.x
    } else if (class(G.x)=="prop.1D") {
      G.grid$x.int <- matrix(data=G.x$int,nrow=(Nx+1),ncol=Ny)
    } else if (length(G.x) == 1) {
      G.grid$x.int <- matrix(data=G.x,nrow=(Nx+1),ncol=Ny)
    } else if (length(G.x) != Nx+1) {
        stop (paste("error: ",Name,".x should be a vector of length 1 or Nx+1", del=""))
    } else {  # correct length
      G.grid$x.int <- matrix(data=G.x,nrow=(Nx+1),ncol=Ny)
    }
    # infilling of y-matrix
    if (is.matrix(G.y)) {
      if (sum(abs(dim(G.y) - c(Nx,Ny+1)))!=0)
        stop (paste("error: ",Name,".y matrix not of correct Nx(Ny+1)dimensions", del=""))
      G.grid$y.int <- G.y
    } else if (class(G.y)=="prop.1D") {
      G.grid$y.int <- t(matrix(data=G.y$int,ncol=Nx,nrow=(Ny+1)))     #Karline: changed 007/2011
    } else if (length(G.y) == 1) {
      G.grid$y.int <- matrix(data=G.y,nrow=Nx,ncol=(Ny+1))
    } else if (length(G.y) != Ny+1) {
        stop (paste("error: ",Name,".y should be a vector of length 1 or Ny+1", del=""))
    } else {  # correct length
      G.grid$y.int <- matrix(data=G.y,nrow=Nx,ncol=(Ny+1))
    }
    G.grid
  }

# Need this for AFDW , D and v

  if (is.null(AFDW.grid)) AFDW.grid <- gridInt(AFDW.x,AFDW.y,"AFDW")
  if (is.null(D.grid))    D.grid <- gridInt(D.x,D.y,"D")
  if (is.null(v.grid))    v.grid <- gridInt(v.x,v.y,"v")

#==============================================================================
# INPUT CHECKS  
#==============================================================================


  if (full.check) {

## check dimensions of input concentrations

    if (!is.null(C.x.up)) {
      if (!((length(C.x.up)==1) || (length(C.x.up)==(Ny))))
        stop("error: C.x.up should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(C.x.down)) {
      if (!((length(C.x.down)==1) || (length(C.x.down)==(Ny))))
        stop("error: C.x.down should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(C.y.up)) {
      if (!((length(C.y.up)==1) || (length(C.y.up)==(Nx))))
        stop("error: C.y.up should be a vector of length 1 or nrow(C)")
    }
    if (!is.null(C.y.down)) {
      if (!((length(C.y.down)==1) || (length(C.y.down)==(Nx))))
        stop("error: C.y.down should be a vector of length 1 or nrow(C)")
    }

# check dimensions of input fluxes

    if (!is.null(flux.x.up)) {
      if (!((length(flux.x.up)==1) || (length(flux.x.up)==(Ny))))
        stop("error: flux.x.up should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(flux.x.down)) {
      if (!((length(flux.x.down)==1) || (length(flux.x.down)==(Ny))))
        stop("error: flux.x.down should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(flux.y.up)) {
      if (!((length(flux.y.up)==1) || (length(flux.y.up)==(Nx))))
        stop("error: flux.y.up should be a vector of length 1 or nrow(C)")
    }

    if (!is.null(flux.y.down)) {
      if (!((length(flux.y.down)==1) || (length(flux.y.down)==(Nx))))
        stop("error: flux.y.down should be a vector of length 1 or nrow(C)")
    }


## check input of grid

    if (is.null(dx) && is.null(dy) && is.null(grid))
      stop("error: dx, dy, and grid cannot be NULL at the same time")

    gn <- names(grid)
    if (! "dx" %in% gn)
      stop("error: grid should be a list that contains 'dx' ")
    if (! "dx.aux" %in% gn)
    	stop("error: grid should be a list that contains 'dx.aux' ")
    if (! "dy" %in% gn)
      stop("error: grid should be a list that contains 'dy' ")
    if (! "dy.aux" %in% gn)
    	stop("error: grid should be a list that contains 'dy.aux' ")
    if (is.null(grid$dx) || is.null(grid$dx.aux))
    	stop("error: the grid should be a list with (numeric) values for 'dx' and 'dx.aux' ")
    if (is.null(grid$dy) || is.null(grid$dy.aux))
    	stop("error: the grid should be a list with (numeric) values for 'dy' and 'dy.aux' ")
    if (any(grid$dx <= 0) || any(grid$dx.aux <= 0) )
    	stop("error: the grid distances dx and dx.aux should always be positive")
    if (any(grid$dy <= 0) || any(grid$dy.aux <= 0) )
    	stop("error: the grid distances dy and dy.aux should always be positive")


## check input of AFDW.grid

    gn <- names(AFDW.grid)
    if (! "x.int" %in% gn)
      stop("error: AFDW.grid should be a list that contains 'x.int', the AFDW values at the interfaces of the grid cells in x-direction")
    if (! "y.int" %in% gn)
      stop("error: AFDW.grid should be a list that contains 'y.int', the AFDW values at the interfaces of the grid cells in y-direction")
    if (is.null(AFDW.grid$x.int))
      stop("error: AFDW.grid$x.int should be a list with (numeric) values")
    if (is.null(AFDW.grid$y.int))
      stop("error: AFDW.grid$y.int should be a list with (numeric) values")
    if (any (AFDW.grid$x.int < 0)||any (AFDW.grid$x.int > 1))
    	stop("error: the AFDW should range between 0 and 1")
    if (any (AFDW.grid$y.int < 0)||any (AFDW.grid$y.int > 1))
	    stop("error: the AFDW should range between 0 and 1")

## check input of D.grid

    gn <- names(D.grid)
    if (! "x.int" %in% gn)
      stop("error: D.grid should be a list that contains 'x.int', the D values at the interfaces of the grid cells in x-direction")
    if (! "y.int" %in% gn)
      stop("error: D.grid should be a list that contains 'y.int', the D values at the interfaces of the grid cells in y-direction")
    if (is.null(D.grid$x.int))
      stop("error: D.grid$x.int should be a list with (numeric) values")
    if (is.null(D.grid$y.int))
      stop("error: D.grid$y.int should be a list with (numeric) values")
    if (any (D.grid$x.int < 0)||any (D.grid$y.int < 0))
    	stop("error: the diffusion coefficient should always be positive")

## check input of v.grid

    gn <- names(v.grid)
    if (! "x.int" %in% gn)
      stop("error: v.grid should be a list that contains 'x.int', the velocity values at the interfaces of the grid cells in x-direction")
    if (! "y.int" %in% gn)
      stop("error: v.grid should be a list that contains 'y.int', the velocity values at the interfaces of the grid cells in y-direction")
    if (is.null(v.grid$x.int))
      stop("error: the advective velocity v.grid$x.int should be a list with (numeric) values")
    if (is.null(v.grid$y.int))
      stop("error: the advective velocity v.grid$y.int should be a list with (numeric) values")

## check input of VF.grid

    gn <- names(VF.grid)
    if (! "x.int" %in% gn)
      stop("error: VF.grid should be a list that contains 'x.int'")
    if (! "y.int" %in% gn)
      stop("error: VF.grid should be a list that contains 'y.int'")
    if (! "x.mid" %in% gn)
      stop("error: VF.grid should be a list that contains 'x.mid'")
    if (! "y.mid" %in% gn)
      stop("error: VF.grid should be a list that contains 'y.mid'")
    if (is.null(VF.grid$x.int) || is.null(VF.grid$y.int) || is.null(VF.grid$x.mid) || is.null(VF.grid$y.mid))
     stop("error: VF should contain (numeric) values")
    if (any (VF.grid$x.int < 0) || any (VF.grid$y.int < 0) || any (VF.grid$x.mid < 0) || any (VF.grid$y.mid < 0))
      stop("error: the VF values should always be positive")

## check input of A.grid
    gn <- names(A.grid)
    if (! "x.int" %in% gn)
      stop("error: A.grid should be a list that contains 'x.int'")
    if (! "y.int" %in% gn)
      stop("error: A.grid should be a list that contains 'y.int'")
    if (! "x.mid" %in% gn)
      stop("error: A.grid should be a list that contains 'x.mid'")
    if (! "y.mid" %in% gn)
      stop("error: A.grid should be a list that contains 'y.mid'")
    if (is.null(A.grid$x.int) || is.null(A.grid$y.int) || is.null(A.grid$x.mid) || is.null(A.grid$y.mid))
     stop("error: the VF should contain (numeric) values")
    if (any (A.grid$x.int < 0) || any (A.grid$y.int < 0) || any (A.grid$x.mid < 0) || any (A.grid$y.mid < 0))
      stop("error: the A values should always be positive")

  }
      
## FUNCTION BODY: CALCULATIONS

## Impose boundary flux at upstream x-boundary when needed
## Default boundary condition is no gradient
  if (! is.null (flux.x.up[1])) {
    nom <- flux.x.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1,] +
           (1-AFDW.grid$x.int[1,])*v.grid$x.int[1,])*C[1,]
    denom <- VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1,]+
             AFDW.grid$x.int[1,]*v.grid$x.int[1,])
    C.x.up <- nom/denom
  }

## Impose boundary flux at downstream x-boundary when needed
## Default boundary condition is no gradient
  if (! is.null (flux.x.down[1])) {
  	nom <- flux.x.down - VF.grid$x.int[(Nx+1),]*(D.grid$x.int[(Nx+1),]/
            grid$dx.aux[Nx+1,] + AFDW.grid$x.int[(Nx+1),]*v.grid$x.int[(Nx+1),])*C[Nx,]
    denom <- -VF.grid$x.int[(Nx+1),]*(D.grid$x.int[(Nx+1),]/grid$dx.aux[Nx+1,]+
            (1-AFDW.grid$x.int[(Nx+1),])*v.grid$x.int[(Nx+1),])
    C.x.down <- nom/denom
  }

# Impose boundary flux at upstream y-boundary when needed
# Default boundary condition is no gradient
  if (! is.null (flux.y.up[1])) {
    nom <- flux.y.up + VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[,1] +
           (1-AFDW.grid$y.int[,1])*v.grid$y.int[,1])*C[,1]
    denom <- VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[,1]+
             AFDW.grid$y.int[,1]*v.grid$y.int[,1])
    C.y.up <- nom/denom
  }

# Impose boundary flux at downstream y-boundary when needed
# Default boundary condition is no gradient
  if (! is.null (flux.y.down[1]))  {
	  nom <- flux.y.down - VF.grid$y.int[,(Ny+1)]*(D.grid$y.int[,(Ny+1)]/
           grid$dy.aux[,(Ny+1)] + AFDW.grid$y.int[,(Ny+1)]*v.grid$y.int[,(Ny+1)])*C[,Ny]
    denom <- -VF.grid$y.int[,(Ny+1)]*(D.grid$y.int[,(Ny+1)]/grid$dy.aux[,Ny+1]+
             (1-AFDW.grid$y.int[,(Ny+1)])*v.grid$y.int[,(Ny+1)])
    C.y.down <- nom/denom
  }

## when upper boundary layer is present, calculate new C.x.up
  if (!is.null(a.bl.x.up) & !is.null(C.x.up[1])) {
	  nom <- a.bl.x.up*C.x.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/
           grid$dx.aux[1,] + (1-AFDW.grid$x.int[1,])*v.grid$x.int[1,])*C[1,]
    denom <- a.bl.x.up + VF.grid$x.int[1,]*(D.grid$x.int[1,]/grid$dx.aux[1,]+
             AFDW.grid$x.int[1,]*v.grid$x.int[1,])
	  C.x.up <- nom/denom
  }

## when lower boundary layer is present, calculate new C.x.down
  if (!is.null(a.bl.x.down) & !is.null(C.x.down[1])) {
	  nom <- a.bl.x.down*C.x.down + VF.grid$x.int[(Nx+1),]*(D.grid$x.int[(Nx+1),]/
           grid$dx.aux[(Nx+1),] + (1-AFDW.grid$x.int[(Nx+1),])*
           v.grid$x.int[(Nx+1),])*C[Nx,]
    denom <- a.bl.x.down + VF.grid$x.int[(Nx+1),]*(D.grid$x.int[(Nx+1),]/
             grid$dx.aux[(Nx+1),]+ AFDW.grid$x.int[(Nx+1),]*v.grid$x.int[(Nx+1),])
	  C.x.down <- nom/denom
  }

## when upper y boundary layer is present, calculate new C.y.up
  if (!is.null(a.bl.y.up) & !is.null(C.y.up[1])) {
	  nom <- a.bl.y.up*C.y.up + VF.grid$y.int[,1]*(D.grid$y.int[,1]/
           grid$dy.aux[,1] + (1-AFDW.grid$y.int[,1])*v.grid$y.int[,1])*C[,1]
    denom <- a.bl.y.up + VF.grid$y.int[,1]*(D.grid$y.int[,1]/grid$dy.aux[,1]+
             AFDW.grid$y.int[,1]*v.grid$y.int[,1])
	  C.y.up <- nom/denom
  }

## when lower y boundary layer is present, calculate new C.y.down
  if (!is.null(a.bl.y.down) & !is.null(C.y.down[1]))   {
	  nom <- a.bl.y.down*C.y.down + VF.grid$y.int[,(Ny+1)]*
           (D.grid$y.int[,(Ny+1)]/grid$dy.aux[,(Ny+1)] +
           (1-AFDW.grid$y.int[,(Ny+1)])*v.grid$y.int[,(Ny+1)])*C[,Ny]
    denom <- a.bl.y.down + VF.grid$y.int[,(Ny+1)]*(D.grid$y.int[,(Ny+1)]/
             grid$dy.aux[,(Ny+1)]+ AFDW.grid$y.int[,(Ny+1)]*v.grid$y.int[,(Ny+1)])
	  C.y.down <- nom/denom
  }

## Calculate diffusive part of the flux
  x.Dif.flux <- as.matrix(-VF.grid$x.int * D.grid$x.int *
                diff(rbind(C.x.up, C, C.x.down, deparse.level = 0))/
                grid$dx.aux)
  y.Dif.flux <- as.matrix(-VF.grid$y.int * D.grid$y.int *
                t(diff(t(cbind(C.y.up,C,C.y.down,deparse.level = 0))))/
                grid$dy.aux)

## Calculate advective part of the flux
  x.Adv.flux <- 0
  
  if (any(v.grid$x.int >0) ) {
    vv <- v.grid$x.int
    vv[vv<0]<-0
    x.Adv.flux <-  x.Adv.flux + as.matrix(VF.grid$x.int * vv * (
                 AFDW.grid$x.int * rbind(C.x.up,C,deparse.level = 0)
                 + (1-AFDW.grid$x.int) * rbind(C,C.x.down,deparse.level = 0)))
  }
  if (any (v.grid$x.int < 0))  {
    vv <- v.grid$x.int
    vv[vv>0]<-0
    x.Adv.flux <-  x.Adv.flux + as.matrix(VF.grid$x.int * vv * (
                    (1-AFDW.grid$x.int) * rbind(C.x.up,C,deparse.level = 0)
                 +   AFDW.grid$x.int * rbind(C,C.x.down,deparse.level = 0)))

  }
  y.Adv.flux <- 0
  if (any(v.grid$y.int >0) ) {
    vv <- v.grid$y.int
    vv[vv<0]<-0
    y.Adv.flux <-  y.Adv.flux + as.matrix(VF.grid$y.int * vv * (
                 AFDW.grid$y.int * cbind(C.y.up,C,deparse.level = 0)
                 + (1-AFDW.grid$y.int) * cbind(C,C.y.down,deparse.level = 0)))
  }
  if (any (v.grid$y.int < 0)) {
    vv <- v.grid$y.int
    vv[vv>0]<-0
    y.Adv.flux <-  y.Adv.flux + as.matrix(VF.grid$y.int * vv * (
                    (1-AFDW.grid$y.int) * cbind(C.y.up,C,deparse.level = 0)
                 +  AFDW.grid$y.int * cbind(C,C.y.down,deparse.level = 0)))
  }

  x.flux <- x.Dif.flux + x.Adv.flux
  y.flux <- y.Dif.flux + y.Adv.flux

## Impose boundary fluxes when needed
## Default boundary condition is no gradient
  if (! is.null (flux.x.up[1]))
    x.flux[1,]   <- flux.x.up
  if (! is.null (flux.x.down[1]))
    x.flux[nrow(x.flux),] <- flux.x.down
    
  if (! is.null (flux.y.up[1]))
    y.flux[,1]   <- flux.y.up
  if (! is.null (flux.y.down[1]))
    y.flux[,ncol(y.flux)] <- flux.y.down

## Calculate rate of change = flux gradient
  dFdx <- - (diff(A.grid$x.int*x.flux)   / A.grid$x.mid)  /grid$dx  / VF.grid$x.mid
  dFdy <- -t(diff(t(A.grid$y.int*y.flux))/t(A.grid$y.mid))/grid$dy / VF.grid$y.mid


  if (!full.output) {
    return (list (dC = dFdx + dFdy,                    # Rate of change due to advective-diffuisve transport in each grid cell
                  flux.x.up = x.flux[1,],                # flux across lower boundary interface; positive = IN
                  flux.x.down = x.flux[nrow(x.flux),],   # flux across lower boundary interface; positive = OUT
                  flux.y.up = y.flux[,1],                # flux across lower boundary interface; positive = IN
                  flux.y.down = y.flux[,ncol(y.flux)]))  # flux across lower boundary interface; positive = OUT

  } else {
    return (list (dC = dFdx + dFdy,                    # Rate of change in the centre of each grid cells
                  C.x.up = C.x.up,                     # concentration at upper interface
                  C.x.down = C.x.down,                 # concentration at upper interface
                  C.y.up = C.y.up,                     # concentration at upper interface
                  C.y.down = C.y.down,                 # concentration at upper interface
                  x.flux = x.flux,                     # flux across at the interface of each grid cell
                  y.flux = y.flux,                     # flux across at the interface of each grid cell
                  flux.x.up = x.flux[1,],               # flux across lower boundary interface; positive = IN
                  flux.x.down = x.flux[nrow(x.flux),],  # flux across lower boundary interface; positive = OUT
                  flux.y.up = y.flux[,1],               # flux across lower boundary interface; positive = IN
                  flux.y.down = y.flux[,ncol(y.flux)])) # flux across lower boundary interface; positive = OUT
  }
} # end tran.2D

