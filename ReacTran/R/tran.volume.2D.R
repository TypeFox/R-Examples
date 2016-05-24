
##==============================================================================
## Transport in a two-dimensional finite volume grid
##==============================================================================

tran.volume.2D <- function(C, 
  C.x.up = C[1,], C.x.down = C[nrow(C),],
  C.y.up = C[,1], C.y.down = C[,ncol(C)],
  C.z = C, masscons = TRUE,
  F.x.up = NULL, F.x.down = NULL, 
  F.y.up = NULL, F.y.down = NULL,
  Disp.grid = NULL, Disp.x = NULL, Disp.y = Disp.x, 
  flow.grid = NULL, flow.x = NULL, flow.y = NULL,
  AFDW.grid = NULL, AFDW.x = 1, AFDW.y = AFDW.x,
  V = NULL,  full.check = FALSE, full.output = FALSE)
											
{
  DD <- dim(C)
  Nx <- nrow(C)
  Ny <- ncol(C)
  if (is.null(V))
      stop("'V' should be specified  ")
      

#==============================================================================
# infilling of grids with only  x.int and y.int needed
#==============================================================================

  gridInt <- function(G.x, G.y, Name)    {   # define a function first

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

# Need this for AFDW , Disp and flow

  if (is.null(AFDW.grid)) AFDW.grid <- gridInt(AFDW.x, AFDW.y, "AFDW")
  if (is.null(Disp.grid)) Disp.grid <- gridInt(Disp.x, Disp.y, "Disp")
  if (is.null(flow.grid)) flow.grid <- gridInt(flow.x, flow.y, "flow")
  
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

    if (!is.null(F.x.up)) {
      if (!((length(F.x.up)==1) || (length(F.x.up)==(Ny))))
        stop("error: F.x.up should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(F.x.down)) {
      if (!((length(F.x.down)==1) || (length(F.x.down)==(Ny))))
        stop("error: F.x.down should be a vector of length 1 or ncol(C)")
    }
    if (!is.null(F.y.up)) {
      if (!((length(F.y.up)==1) || (length(F.y.up)==(Nx))))
        stop("error: F.y.up should be a vector of length 1 or nrow(C)")
    }

    if (!is.null(F.y.down)) {
      if (!((length(F.y.down)==1) || (length(F.y.down)==(Nx))))
        stop("error: F.y.down should be a vector of length 1 or nrow(C)")
    }


## check input of volumes
    if (any(V <= 0))
    	stop("error: the volumes should always be positive")
    if (nrow(V) != Nx || ncol(V) != Ny)
    	stop("error: the dimension of 'V' should be = dimension of 'C'")

## If mass should be conserved...
    if (masscons) {
    if (nrow(C.z) != Nx || ncol(C.z) != Ny)
    	stop("error: the dimension of 'C.z' should be = dimension of 'C' if 'masscons' = TRUE")
   }
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

## check input of Disp.grid

    gn <- names(Disp.grid)
    if (! "x.int" %in% gn)
      stop("error: Disp.grid should be a list that contains 'x.int', the Disp values at the interfaces of the grid cells in x-direction")
    if (! "y.int" %in% gn)
      stop("error: Disp.grid should be a list that contains 'y.int', the Disp values at the interfaces of the grid cells in y-direction")
    if (is.null(Disp.grid$x.int))
      stop("error: Disp.grid$x.int should be a list with (numeric) values")
    if (is.null(Disp.grid$y.int))
      stop("error: Disp.grid$y.int should be a list with (numeric) values")
    if (any (Disp.grid$x.int < 0)||any (Disp.grid$y.int < 0))
    	stop("error: the diffusion coefficient, Disp, should always be positive")

## check input of flow.grid

    gn <- names(flow.grid)
    if (! "x.int" %in% gn)
      stop("error: flow.grid should be a list that contains 'x.int', the flows at the interfaces of the grid cells in x-direction")
    if (! "y.int" %in% gn)
      stop("error: flow.grid should be a list that contains 'y.int', the flows at the interfaces of the grid cells in y-direction")
    if (is.null(flow.grid$x.int))
      stop("error: the flows flow.grid$x.int should be a list with (numeric) values")
    if (is.null(flow.grid$y.int))
      stop("error: the flows flow.grid$y.int should be a list with (numeric) values")
  }
      
## FUNCTION BODY: CALCULATIONS
## Calculate diffusive part of the flow
  x.Dif.flow <- as.matrix(-Disp.grid$x.int *
                diff(rbind(C.x.up, C, C.x.down, deparse.level = 0)))
  y.Dif.flow <- as.matrix(-Disp.grid$y.int *
                t(diff(t(cbind(C.y.up,C,C.y.down,deparse.level = 0)))))

## Calculate advective part of the flow
  x.Adv.flow <- 0
  
  if (any(flow.grid$x.int >0) ) {
    vv <- flow.grid$x.int
    vv[vv < 0]<-0
    x.Adv.flow <-  x.Adv.flow + as.matrix(vv * (
                 AFDW.grid$x.int * rbind(C.x.up,C,deparse.level = 0)
                 + (1-AFDW.grid$x.int) * rbind(C,C.x.down,deparse.level = 0)))
  }
  if (any (flow.grid$x.int < 0))  {
    vv <- flow.grid$x.int
    vv[vv>0]<-0
    x.Adv.flow <-  x.Adv.flow + as.matrix(vv * (
                    (1-AFDW.grid$x.int) * rbind(C.x.up,C,deparse.level = 0)
                 +   AFDW.grid$x.int * rbind(C,C.x.down,deparse.level = 0)))

  }
  y.Adv.flow <- 0
  if (any(flow.grid$y.int >0) ) {
    vv <- flow.grid$y.int
    vv[vv<0]<-0
    y.Adv.flow <-  y.Adv.flow + as.matrix(vv * (
                 AFDW.grid$y.int * cbind(C.y.up,C,deparse.level = 0)
                 + (1-AFDW.grid$y.int) * cbind(C,C.y.down,deparse.level = 0)))
  }
  if (any (flow.grid$y.int < 0)) {
    vv <- flow.grid$y.int
    vv[vv>0]<-0
    y.Adv.flow <-  y.Adv.flow + as.matrix(vv * (
                    (1-AFDW.grid$y.int) * cbind(C.y.up,C,deparse.level = 0)
                 +  AFDW.grid$y.int * cbind(C,C.y.down,deparse.level = 0)))
  }

  x.flow <- x.Dif.flow + x.Adv.flow
  y.flow <- y.Dif.flow + y.Adv.flow

  z.Adv.flow <- 0
  ww <- 0
  if (masscons) {
   ww <- -(flow.grid$x.int[-1, ] - flow.grid$x.int[-nrow(flow.grid$x.int), ] + 
           flow.grid$y.int[ ,-1] - flow.grid$y.int[ ,-ncol(flow.grid$y.int)] )
   AFDW.z <- AFDW.grid$y.int[1]
   if (any(ww >0) ) {
    vv <- ww
    vv[vv<0]<-0
    z.Adv.flow <-  z.Adv.flow + as.matrix(vv * (
                 AFDW.z * C + (1-AFDW.z) * C.z))
  }
  if (any (ww < 0)) {
    vv <- ww
    vv[vv>0]<-0
    z.Adv.flow <-  z.Adv.flow + as.matrix(vv * (
                    (1-AFDW.z) * C +  AFDW.z * C.z))
    }
  } 

## Impose boundary fluxes when needed
## Default boundary condition is no gradient
  if (! is.null (F.x.up[1]))
    x.flow[1,]   <- F.x.up
  if (! is.null (F.x.down[1]))
    x.flow[nrow(x.flow),] <- F.x.down
    
  if (! is.null (F.y.up[1]))
    y.flow[,1]   <- F.y.up
  if (! is.null (F.y.down[1]))
    y.flow[,ncol(y.flow)] <- F.y.down

## Calculate rate of change = flux gradient
  dFdx <- - diff(x.flow)/ V  
  dFdy <- -t(diff(t(y.flow))/t(V))
  dFdz <- - z.Adv.flow/V

  if (!full.output) {
    return (list (dC = dFdx + dFdy + dFdz,            # Rate of change due to advective-diffuisve transport in each grid cell
                  F.x.up = x.flow[1,],                # flux across upper boundary interface; positive = IN
                  F.x.down = x.flow[nrow(x.flow),],   # flux across lower boundary interface; positive = OUT
                  F.y.up = y.flow[,1],                
                  F.y.down = y.flow[,ncol(y.flow)],  
                  F.z  = z.Adv.flow                    ))  

  } else {
    return (list (dC = dFdx + dFdy,                    # Rate of change in the centre of each grid cells
                  x.flow = x.flow,                     # flow across at the interface of each grid cell
                  y.flow = y.flow,                     # flow across at the interface of each grid cell
                  F.x.up = x.flow[1,],              
                  F.x.down = x.flow[nrow(x.flow),], 
                  F.y.up = y.flow[,1],              
                  F.y.down = y.flow[,ncol(y.flow)], 
                  F.z = z.Adv.flow ))
  }
} # end tran.2D

