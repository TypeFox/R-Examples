
##==============================================================================
## Take differences of an array
##==============================================================================

diff3D <- function(Array, along) {
  dimens <- dim(Array)
  if (length(dimens) != 3)
    stop("'Array' should be an array with dimension = 3")

  if (along == 1)
    return(Array[2:dimens[1], , ] - Array[1:(dimens[1]-1), , ])
  else if (along == 2)
    return(Array[, 2:dimens[2], ] - Array[, 1:(dimens[2]-1), ])
  else if (along == 3)
    return(Array[, , 2:dimens[3]] - Array[, , 1:(dimens[3]-1)])
  else  
    stop("'along' should be 1, 2, or 3")
}

# function to bind two matrices to an array, to left and right
Mbind <- function (Matleft = NULL, Array, Matright = NULL, along = 1)  {
  dimens <- dim(Array)
  if (length(dimens) != 3)
    stop("'Array' should be an array with dimension = 3")
  
  if (along < 1 | along > 3) 
    stop("'along' should be 1, 2, or 3")

  if (! is.null(Matright)) {
    if (sum(abs(dimens[-along] - dim(Matright))) != 0) 
    stop("'Matright' not compatible with Array")
    dimens[along] <- dimens[along] + 1
    Array <- mbindr(Array, Matright, along)
  }
  if (! is.null(Matleft)) {
   if (sum(abs(dimens[-along] - dim(Matleft))) != 0) 
    stop("'Matleft' not compatible with Array")
    dimens[along] <- dimens[along] + 1
    Array <- mbindl(Matleft, Array, along)
  }
  Array  
}

# function to bind two matrices to an array, to left and right
# NO error checking
mbind <- function (Mat1, Array, Mat2, along = 1)  {
  dimens <- dim(Array)
  
  dimens[along] <- dimens[along] + 2
     if (along == 3)
       array(dim = dimens, data = c(Mat1, Array, Mat2))
     else if (along == 1)
       aperm(array(dim = dimens[c(3, 2, 1)],
         data = c(t(Mat1), aperm(Array, c(3, 2, 1)), t(Mat2))),
         c(3, 2, 1))
     else if (along == 2)
       aperm(array(dim = dimens[c(1, 3, 2)],
         data = c(Mat1, aperm(Array, c(1, 3, 2)), Mat2)),
         c(1, 3, 2))
}

# function to bind a matrix to an array on the left
mbindl <- function (Mat1, Array, along = 1)  {
  dimens <- dim(Array)

  dimens[along] <- dimens[along]+1
  if (along == 3)
       array(dim = dimens, data = c(Mat1, Array))
  else if (along == 1)
       aperm(array(dim = dimens[c(3, 2, 1)],
         data = c(t(Mat1), aperm(Array, c(3, 2, 1)))),
         c(3, 2, 1))
  else if (along == 2)
       aperm(array(dim = dimens[c(1, 3, 2)],
         data = c(Mat1, aperm(Array, c(1, 3, 2)))),
         c(1, 3, 2))
   }
# function to bind a matrix to an array on the right
mbindr <- function (Array, Mat2, along = 1)  {
  dimens <- dim(Array)

  dimens[along] <- dimens[along]+1
  if (along == 3)
       array(dim = dimens, data = c(Array, Mat2))
  else if (along == 1)
       aperm(array(dim = dimens[c(3, 2, 1)],
         data = c(aperm(Array, c(3, 2, 1)), t(Mat2))),
         c(3, 2, 1))
  else if (along == 2)
       aperm(array(dim = dimens[c(1, 3, 2)],
         data = c(aperm(Array, c(1, 3, 2)), Mat2)),
         c(1, 3, 2))
}

##==============================================================================
## Transport in a three-dimensional finite difference grid
##==============================================================================

tran.volume.3D <- function(C, 
  C.x.up = C[1, , ], C.x.down = C[dim(C)[1], ,],
  C.y.up = C[, 1, ], C.y.down = C[,dim(C)[2],],
  C.z.up = C[, , 1], C.z.down = C[,,dim(C)[3]],
  F.x.up = NULL, F.x.down = NULL,
  F.y.up = NULL, F.y.down = NULL,
  F.z.up = NULL, F.z.down = NULL,
  Disp.grid = NULL, Disp.x = NULL, Disp.y = Disp.x, Disp.z = Disp.x,
  flow.grid = NULL, flow.x = 0, flow.y = 0, flow.z = 0,
  AFDW.grid = NULL, AFDW.x = 1, AFDW.y = AFDW.x, AFDW.z = AFDW.x,
  V = NULL, full.check = FALSE, full.output = FALSE)
											
{
  Nx <- dim(C)[1]
  Ny <- dim(C)[2]
  Nz <- dim(C)[3]

  if(length (C.x.up)   == 1) C.x.up   <- matrix(nrow=Ny,ncol=Nz,C.x.up)
  if(length (C.x.down) == 1) C.x.down <- matrix(nrow=Ny,ncol=Nz,C.x.down)
  if(length (C.y.up)   == 1) C.y.up   <- matrix(nrow=Nx,ncol=Nz,C.y.up)
  if(length (C.y.down) == 1) C.y.down <- matrix(nrow=Nx,ncol=Nz,C.y.down)
  if(length (C.z.up)   == 1) C.z.up   <- matrix(nrow=Nx,ncol=Ny,C.z.up)
  if(length (C.z.down) == 1) C.z.down <- matrix(nrow=Nx,ncol=Ny,C.z.down)

  if (is.null(V))
      stop("volume of each grid cell, 'V' should be specified  ")

#==============================================================================
# infilling of grids with only  x.int and y.int needed
#==============================================================================

  gridInt <- function(G.x,G.y,G.z,Name)   {    # define a function first

    # check if G.x and G.y and G.z is not NULL
    if (is.null(G.x) | is.null(G.y) | is.null(G.z))
      stop( (paste("error: ",Name,"and (",Name,".x and", Name,".y and",
             Name,".z) cannot be NULL at the same time", del="")))
    G.grid <- list()
    # infilling of x-array
    if (is.array(G.x)) {
      if (sum(abs(dim(G.x) - c(Nx+1,Ny,Nz)))!=0)
        stop (paste("error: ",Name,".x array not of correct (Nx+1) Ny ,Nz dimensions", del=""))
      G.grid$x.int <- G.x
    } else if (length(G.x) == 1) {
      G.grid$x.int <- array(data=G.x,dim=c(Nx+1,Ny,Nz))
    } else if (length(G.x) == Nx+1) {
      G.grid$x.int <- array(data=G.x,dim=c(Nx+1,Ny,Nz))
    } else
      stop (paste("error: ",Name,".x should be one element, a vector, or an array", del=""))

    # infilling of y-array
    if (is.array(G.y)) {
      if (sum(abs(dim(G.y) - c(Nx,Ny+1,Nz)))!=0)
        stop (paste("error: ",Name,".y array not of correct Nx, Ny+1 ,Nz dimensions", del=""))
      G.grid$y.int <- G.y
    } else if (length(G.y) == 1) {
      G.grid$y.int <- array(data=G.y,dim=c(Nx,Ny+1,Nz))
    } else if (length(G.y) == Ny+1) {
      G.grid$y.int <- aperm(array(data=G.y,dim=c(Ny+1,Nx,Nz)),c(2,1,3))
    } else
      stop (paste("error: ",Name,".y should be one element, a vector or an array", del=""))

    # infilling of z-array
    if (is.array(G.z)) {
      if (sum(abs(dim(G.z) - c(Nx,Ny,Nz+1)))!=0)
        stop (paste("error: ",Name,".z array not of correct Nx, Ny ,Nz+1 dimensions", del=""))
      G.grid$z.int <- G.z
    } else if (length(G.z) == 1) {
      G.grid$z.int <- array(data=G.z,dim=c(Nx,Ny,Nz+1))
    } else if (length(G.z) == Nz+1) {
      G.grid$z.int <- aperm(array(data=G.z,dim=c(Nz+1,Ny,Nx)),c(3,2,1))
    } else
      stop (paste("error: ",Name,".z should be one element, a vector or an array", del=""))

    G.grid
  }

# Need this for AFDW , D and v

  if (is.null(AFDW.grid)) AFDW.grid <- gridInt(AFDW.x,AFDW.y,AFDW.z,"AFDW")
  if (is.null(Disp.grid)) Disp.grid <- gridInt(Disp.x,Disp.y,Disp.z,"Disp")
  if (is.null(flow.grid)) flow.grid <- gridInt(flow.x,flow.y,flow.z,"flow")


#==============================================================================
# INPUT CHECKS  
#==============================================================================


  if (full.check) {

## check dimensions of input concentrations

    if (!is.null(C.x.up)) {
      if (!((length(C.x.up)==1) || (sum(abs(dim(C.x.up) - c(Ny,Nz)))==0)))
        stop("error: C.x.up should be of length 1 or a matrix with dim Ny, Nz")
    }
    if (!is.null(C.x.down)) {
      if (!((length(C.x.down)==1) || (sum(abs(dim(C.x.down) - c(Ny,Nz)))==0)))
        stop("error: C.x.down should be of length 1 or a matrix with dim Ny, Nz")
    }
    if (!is.null(C.y.up)) {
      if (!((length(C.y.up)==1) || (sum(abs(dim(C.y.up) - c(Nx,Nz)))==0)))
        stop("error: C.y.up should be of length 1 or a matrix with dim Nx, Nz")
    }
    if (!is.null(C.y.down)) {
      if (!((length(C.y.down)==1) || (sum(abs(dim(C.y.down) - c(Nx,Nz)))==0)))
        stop("error: C.y.down should be of length 1 or a matrix with dim Nx, Nz")
    }
    if (!is.null(C.z.up)) {
      if (!((length(C.z.up)==1) || (sum(abs(dim(C.z.up) - c(Nx,Ny)))==0)))
        stop("error: C.z.up should be of length 1 or a matrix with dim Nx, Ny")
    }
    if (!is.null(C.z.down)) {
      if (!((length(C.z.down)==1) || (sum(abs(dim(C.z.down) - c(Nx,Ny)))==0)))
        stop("error: C.z.down should be of length 1 or a matrix with dim Nx, Ny")
    }

# check dimensions of input fluxes
    if (!is.null(F.x.up)) {
      if (!((length(F.x.up)==1) || (sum(abs(dim(F.x.up) - c(Ny,Nz)))!=0)))
        stop("error: F.x.up should be of length 1 or a matrix with dim Ny, Nz")
    }
    if (!is.null(F.x.down)) {
      if (!((length(F.x.down)==1) || (sum(abs(dim(F.x.down) - c(Ny,Nz)))!=0)))
        stop("error: F.x.down should be of length 1 or a matrix with dim Ny, Nz")
    }
    if (!is.null(F.y.up)) {
      if (!((length(F.y.up)==1) || (sum(abs(dim(F.y.up) - c(Nx,Nz)))!=0)))
        stop("error: F.y.up should be of length 1 or a matrix with dim Nx, Nz")
    }
    if (!is.null(F.y.down)) {
      if (!((length(F.y.down)==1) || (sum(abs(dim(F.y.down) - c(Nx,Nz)))!=0)))
        stop("error: F.y.down should be of length 1 or a matrix with dim Nx, Nz")
    }
    if (!is.null(F.z.up)) {
      if (!((length(F.z.up)==1) || (sum(abs(dim(F.z.up) - c(Nx,Ny)))!=0)))
        stop("error: F.z.up should be of length 1 or a matrix with dim Nx, Ny")
    }
    if (!is.null(F.z.down)) {
      if (!((length(F.z.down)==1) || (sum(abs(dim(F.z.down) - c(Nx,Ny)))!=0)))
        stop("error: F.z.down should be of length 1 or a matrix with dim Nx, Ny")
    }

## check input of volumes
    if (any(V <= 0))
    	stop("error: the volumes should always be positive")
    if (dim(V)[1] != Nx || dim(V)[2] != Ny || dim(V)[3] != Nz)
    	stop("error: the dimension of 'V' should be = dimension of 'C'")

## check input of AFDW.grid

    if (is.null(AFDW.x) && is.null(AFDW.y) && is.null(AFDW.z) && is.null(AFDW.grid))
      stop("error: AFDW.x, AFDW.y, AFDW.z and AFDW.grid cannot be NULL at the same time")

    gn <- names(AFDW.grid)
    if (! "x.int" %in% gn)
      stop("error: AFDW.grid should be a list that contains 'x.int', the values at the interfaces of the grid cells in x-direction")
    if (! "y.int" %in% gn)
      stop("error: AFDW.grid should be a list that contains 'y.int', the values at the interfaces of the grid cells in y-direction")
    if (! "z.int" %in% gn)
      stop("error: AFDW.grid should be a list that contains 'z.int', the values at the interfaces of the grid cells in z-direction")
    if (is.null(AFDW.grid$x.int))
      stop("error: AFDW.grid$x.int should be a list with (numeric) values")
    if (is.null(AFDW.grid$y.int))
      stop("error: AFDW.grid$y.int should be a list with (numeric) values")
    if (is.null(AFDW.grid$z.int))
      stop("error: AFDW.grid$z.int should be a list with (numeric) values")
    if (any (AFDW.grid$x.int < 0)||any (AFDW.grid$x.int > 1))
    	stop("error: the AFDW should range between 0 and 1")
    if (any (AFDW.grid$y.int < 0)||any (AFDW.grid$y.int > 1))
	    stop("error: the AFDW should range between 0 and 1")
    if (any (AFDW.grid$z.int < 0)||any (AFDW.grid$z.int > 1))
	    stop("error: the AFDW should range between 0 and 1")

## check input of Disp.grid

    if (is.null(Disp.x) && is.null(Disp.y) && is.null(Disp.grid))
      stop("error: Disp.x, Disp.y, and Disp.grid cannot be NULL at the same time")

    gn <- names(Disp.grid)
    if (! "x.int" %in% gn)
      stop("error: DDisp.grid should be a list that contains 'x.int', the Disp values at the interfaces of the grid cells in x-direction")
    if (! "y.int" %in% gn)
      stop("error: Disp.grid should be a list that contains 'y.int', the Disp values at the interfaces of the grid cells in y-direction")
    if (! "z.int" %in% gn)
      stop("error: Disp.grid should be a list that contains 'z.int', the Disp values at the interfaces of the grid cells in z-direction")
    if (is.null(Disp.grid$x.int))
      stop("error: Disp.grid$x.int should be a list with (numeric) values")
    if (is.null(Disp.grid$y.int))
      stop("error: Disp.grid$y.int should be a list with (numeric) values")
    if (is.null(Disp.grid$z.int))
      stop("error: Disp.grid$z.int should be a list with (numeric) values")
    if (any (Disp.grid$x.int < 0)||any (Disp.grid$y.int < 0)||any (Disp.grid$z.int < 0))
    	stop("error: the diffusion coefficient should always be positive")

## check input of flow.grid

    if (is.null(flow.x) && is.null(flow.y) && is.null(flow.grid))
      stop("error: flow.x, flow.y, and flow.grid cannot be NULL at the same time")

    gn <- names(flow.grid)
    if (! "x.int" %in% gn)
      stop("error: flow.grid should be a list that contains 'x.int', the flow values at the interfaces of the grid cells in x-direction")
    if (! "y.int" %in% gn)
      stop("error: flow.grid should be a list that contains 'y.int', the flow values at the interfaces of the grid cells in y-direction")
    if (! "z.int" %in% gn)
      stop("error: flow.grid should be a list that contains 'z.int', the flow values at the interfaces of the grid cells in z-direction")
    if (is.null(flow.grid$x.int))
      stop("error: flow.grid$x.int should be a list with (numeric) values")
    if (is.null(flow.grid$y.int))
      stop("error: flow.grid$y.int should be a list with (numeric) values")
    if (is.null(flow.grid$z.int))
      stop("error: flow.grid$z.int should be a list with (numeric) values")


  }
## FUNCTION BODY: CALCULATIONS

## Calculate diffusive part of the flow
  x.Dif.flow <- -Disp.grid$x.int *
                diff3D(mbind(C.x.up, C, C.x.down, along=1), along=1)
  y.Dif.flow <-  -Disp.grid$y.int *
                diff3D(mbind(C.y.up, C, C.y.down, along=2), along=2)
  z.Dif.flow <-  -Disp.grid$z.int *
                diff3D(mbind(C.z.up, C, C.z.down, along=3), along=3)

## Calculate advective part of the flow
  x.Adv.flow <- 0
  
  if (any(flow.grid$x.int > 0) ) {
    vv <- flow.grid$x.int
    vv[vv < 0] <- 0
    x.Adv.flow <-  x.Adv.flow + vv * (
                      AFDW.grid$x.int  * mbindl (C.x.up,   C, along = 1)
                 + (1-AFDW.grid$x.int) * mbindr (C, C.x.down, along = 1))
  }
  if (any (flow.grid$x.int < 0))  {
    vv <- flow.grid$x.int
    vv [vv > 0] <- 0
    x.Adv.flow <-  x.Adv.flow + vv * (
                (1-AFDW.grid$x.int) * mbindl(C.x.up,   C, along = 1)
                 + AFDW.grid$x.int  * mbindr(C, C.x.down, along = 1))
  }
  y.Adv.flow <- 0
  if (any(flow.grid$y.int >0) ) {
    vv <- flow.grid$y.int
    vv[vv<0]<-0
    y.Adv.flow <-  y.Adv.flow + vv * (
                      AFDW.grid$y.int  * mbindl(C.y.up  , C, along = 2)
                 + (1-AFDW.grid$y.int) * mbindr(C, C.y.down, along = 2))
  }
  if (any (flow.grid$y.int < 0))  {
    vv <- flow.grid$y.int
    vv[vv>0]<-0
    y.Adv.flow <-  y.Adv.flow + vv * (
                (1-AFDW.grid$y.int) * mbindl(C.y.up  , C, along = 2)
                  + AFDW.grid$y.int * mbindr(C, C.y.down, along = 2))
  }
  z.Adv.flow <- 0
  if (any(flow.grid$z.int > 0 ) ) {
    vv <- flow.grid$z.int
    vv[vv<0]<-0
    z.Adv.flow <-  z.Adv.flow + vv * (
                       AFDW.grid$z.int * mbindl(C.z.up,   C, along = 3)
                 + (1-AFDW.grid$z.int) * mbindr(C, C.z.down, along = 3))
  }
  if (any (flow.grid$z.int < 0))  {
    vv <- flow.grid$z.int
    vv[vv>0]<-0
    z.Adv.flow <-  z.Adv.flow + vv * (
                  (1-AFDW.grid$z.int) * mbindl(C.z.up,   C, along = 3)
                 +    AFDW.grid$z.int * mbindr(C, C.z.down, along = 3))

  }

  x.flow <- x.Dif.flow + x.Adv.flow
  y.flow <- y.Dif.flow + y.Adv.flow
  z.flow <- z.Dif.flow + z.Adv.flow

## Impose boundary fluxes when needed
## Default boundary condition is no gradient
  if (! is.null (F.x.up[1]))
    x.flow[1,,]   <- F.x.up
  if (! is.null (F.x.down[1]))
    x.flow[dim(x.flow)[1],,] <- F.x.down
    
  if (! is.null (F.y.up[1]))
    y.flow[,1,]   <- F.y.up
  if (! is.null (F.y.down[1]))
    y.flow[,dim(y.flow)[2],] <- F.y.down

  if (! is.null (F.z.up[1]))
    z.flow[,,1]   <- F.z.up
  if (! is.null (F.z.down[1]))
    z.flow[,,dim(z.flow)[3]] <- F.z.down

## Calculate rate of change = flow gradient     NOG DOEN
  dFdx <- - (diff3D(x.flow,along=1))/ V
  dFdy <- - (diff3D(y.flow,along=2))/ V
  dFdz <- - (diff3D(z.flow,along=3))/ V

  if (!full.output) {
    return (list (dC = dFdx + dFdy + dFdz,               # Rate of change due to advective-diffuisve transport in each grid cell
                  F.x.up = x.flow[1,,],                  # flow across lower boundary interface; positive = IN
                  F.x.down = x.flow[dim(x.flow)[1],,],   # flow across lower boundary interface; positive = OUT
                  F.y.up = y.flow[,1,],                  # flow across lower boundary interface; positive = IN
                  F.y.down = y.flow[,dim(y.flow)[2],],   # flow across lower boundary interface; positive = OUT
                  F.z.up = z.flow[,,1],                  # flow across lower boundary interface; positive = IN
                  F.z.down = z.flow[,,dim(z.flow)[3]]))  # flow across lower boundary interface; positive = OUT

  } else {
    return (list (dC = dFdx + dFdy + dFdz,               # Rate of change due to transport in each grid cell
                  x.flow = x.flow,                       # flow across at the interface of each grid cell
                  y.flow = y.flow,                       # flow across at the interface of each grid cell
                  z.flow = z.flow,                       # flow across at the interface of each grid cell
                  F.x.up = x.flow[1,,],                  # flow across lower boundary interface; positive = IN
                  F.x.down = x.flow[dim(x.flow)[1],,],   # flow across lower boundary interface; positive = OUT
                  F.y.up = y.flow[,1,],                  # flow across lower boundary interface; positive = IN
                  F.y.down = y.flow[,dim(y.flow)[2],],   # flow across lower boundary interface; positive = OUT
                  F.z.up = z.flow[,,1],                  # flow across lower boundary interface; positive = IN
                  F.z.down = z.flow[,,dim(z.flow)[3]]))  # flow across lower boundary interface; positive = OUT
  }
} # end tran.3D

