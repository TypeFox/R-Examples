## -----------------------------------------------------------------------------
## GRADIENT functions for internal use of the ReacTran functions
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## Binding one or two matrices to an array, to left and right
## -----------------------------------------------------------------------------


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


# Differences of a 3-D array, no checks, internal use only
Diff3D <- function(Array, along) {
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

##==============================================================================
## Transport in a three-dimensional finite difference grid
##==============================================================================

tran.3D <- function(C, 
  C.x.up = C[1, , ], C.x.down = C[dim(C)[1], ,],
  C.y.up = C[, 1, ], C.y.down = C[,dim(C)[2],],
  C.z.up = C[, , 1], C.z.down = C[,,dim(C)[3]],
  flux.x.up = NULL, flux.x.down = NULL,
  flux.y.up = NULL, flux.y.down = NULL,
  flux.z.up = NULL, flux.z.down = NULL,
  a.bl.x.up = NULL, a.bl.x.down = NULL, 
  a.bl.y.up = NULL, a.bl.y.down = NULL, 
  a.bl.z.up = NULL, a.bl.z.down = NULL, 
  D.grid = NULL, D.x = NULL, D.y = D.x, D.z = D.x,
  v.grid = NULL, v.x = 0, v.y = 0, v.z = 0,
  AFDW.grid = NULL, AFDW.x = 1, AFDW.y = AFDW.x, AFDW.z = AFDW.x,
  VF.grid = NULL, VF.x = 1, VF.y = VF.x, VF.z = VF.x,
  A.grid = NULL, A.x = 1, A.y = 1, A.z = 1,
  grid = NULL, dx = NULL, dy = NULL, dz = NULL,
  full.check = FALSE, full.output = FALSE)
											
{
  if (is.null(grid))
   if (is.null(dx) | is.null(dy) | is.null(dz))
      stop("error: either 'grid' or ('dx' and 'dy' and 'dz') should be specified  ")

  Nx <- dim(C)[1]
  Ny <- dim(C)[2]
  Nz <- dim(C)[3]

  if(length (C.x.up)   == 1) C.x.up   <- matrix(nrow=Ny,ncol=Nz,C.x.up)
  if(length (C.x.down) == 1) C.x.down <- matrix(nrow=Ny,ncol=Nz,C.x.down)
  if(length (C.y.up)   == 1) C.y.up   <- matrix(nrow=Nx,ncol=Nz,C.y.up)
  if(length (C.y.down) == 1) C.y.down <- matrix(nrow=Nx,ncol=Nz,C.y.down)
  if(length (C.z.up)   == 1) C.z.up   <- matrix(nrow=Nx,ncol=Ny,C.z.up)
  if(length (C.z.down) == 1) C.z.down <- matrix(nrow=Nx,ncol=Ny,C.z.down)

# DEFAULT INFILLING OF GRID PARAMETERS

# infilling of 3D grid
  if (is.null(grid)) {
     DX    <-  if (is.list(dx)) dx$dx else rep(dx,length.out=Nx)
     DXaux <-  if (is.list(dx)) dx$dx.aux else 0.5*(c(0,rep(dx,length.out=Nx))+
                                  c(rep(dx,length.out=Nx),0))
     DY    <-  if (is.list(dy)) dy$dx else rep(dy,length.out=Ny)
     DYaux <-  if (is.list(dy)) dy$dx.aux else 0.5*(c(0,rep(dy,length.out=Ny))+
                                  c(rep(dy,length.out=Ny),0))
     DZ    <-  if (is.list(dz)) dz$dx else rep(dz,length.out=Nz)
     DZaux <-  if (is.list(dz)) dz$dx.aux else 0.5*(c(0,rep(dz,length.out=Nz))+
                                  c(rep(dz,length.out=Nz),0))

    grid <- list(
      dx    = DX,  dx.aux= DXaux,
      dy    = DY,  dy.aux= DYaux,
      dz    = DZ,  dz.aux= DZaux
                 )
   }

#==============================================================================
# infilling of grids with x.int, y.int, z.int, x.mid, y.mid, z.mid needed
#==============================================================================
  gridFill <- function(G.x,G.y,G.z,Name)    # define a function first
  {
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
      G.grid$x.mid <- 0.5*(G.x[1:Nx,,]+G.x[2:(Nx+1),,])
    } else if (length(G.x) == 1) {
      G.grid$x.int <- array(data=G.x,dim=c(Nx+1,Ny,Nz))
      G.grid$x.mid <- array(data=G.x,dim=c(Nx,Ny,Nz))
    } else
      stop (paste("error: ",Name,".x should be one element or an array", del=""))

    # infilling of y-array
    if (is.array(G.y)) {
      if (sum(abs(dim(G.y) - c(Nx,Ny+1,Nz)))!=0)
        stop (paste("error: ",Name,".x array not of correct Nx, Ny+1 ,Nz dimensions", del=""))
      G.grid$y.int <- G.y
      G.grid$y.mid <- 0.5*(G.y[,1:Ny,]+G.y[,2:(Ny+1),])
    } else if (length(G.y) == 1) {
      G.grid$y.int <- array(data=G.y,dim=c(Nx,Ny+1,Nz))
      G.grid$y.mid <- array(data=G.y,dim=c(Nx,Ny,Nz))
    } else
      stop (paste("error: ",Name,".y should be one element or an array", del=""))

    # infilling of z-array
    if (is.array(G.z)) {
      if (sum(abs(dim(G.z) - c(Nx,Ny,Nz+1)))!=0)
        stop (paste("error: ",Name,".z array not of correct Nx, Ny ,Nz+1 dimensions", del=""))
      G.grid$z.int <- G.z
      G.grid$z.mid <- 0.5*(G.z[,1:Ny,]+G.z[,2:(Ny+1),])
    } else if (length(G.z) == 1) {
      G.grid$z.int <- array(data=G.z,dim=c(Nx,Ny,Nz+1))
      G.grid$z.mid <- array(data=G.z,dim=c(Nx,Ny,Nz))
    } else
      stop (paste("error: ",Name,".z should be one element or an array", del=""))

    G.grid
  }

# Need this for VF and A (volume fraction and surface

  if (is.null(VF.grid)) VF.grid <- gridFill(VF.x,VF.y,VF.z,"VF")
  if (is.null(A.grid)) A.grid <- gridFill(A.x,A.y,A.z,"A")

#==============================================================================
# infilling of other grids with only  x.int and y.int needed
#==============================================================================

  gridInt <- function(G.x,G.y,G.z,Name)     # define a function first
  {
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
  if (is.null(D.grid)) D.grid <- gridInt(D.x,D.y,D.z,"D")
  if (is.null(v.grid)) v.grid <- gridInt(v.x,v.y,v.z,"v")


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
    if (!is.null(flux.x.up)) {
      if (!((length(flux.x.up)==1) || (sum(abs(dim(flux.x.up) - c(Ny,Nz)))!=0)))
        stop("error: flux.x.up should be of length 1 or a matrix with dim Ny, Nz")
    }
    if (!is.null(flux.x.down)) {
      if (!((length(flux.x.down)==1) || (sum(abs(dim(flux.x.down) - c(Ny,Nz)))!=0)))
        stop("error: flux.x.down should be of length 1 or a matrix with dim Ny, Nz")
    }
    if (!is.null(flux.y.up)) {
      if (!((length(flux.y.up)==1) || (sum(abs(dim(flux.y.up) - c(Nx,Nz)))!=0)))
        stop("error: flux.y.up should be of length 1 or a matrix with dim Nx, Nz")
    }
    if (!is.null(flux.y.down)) {
      if (!((length(flux.y.down)==1) || (sum(abs(dim(flux.y.down) - c(Nx,Nz)))!=0)))
        stop("error: flux.y.down should be of length 1 or a matrix with dim Nx, Nz")
    }
    if (!is.null(flux.z.up)) {
      if (!((length(flux.z.up)==1) || (sum(abs(dim(flux.z.up) - c(Nx,Ny)))!=0)))
        stop("error: flux.z.up should be of length 1 or a matrix with dim Nx, Ny")
    }
    if (!is.null(flux.z.down)) {
      if (!((length(flux.z.down)==1) || (sum(abs(dim(flux.z.down) - c(Nx,Ny)))!=0)))
        stop("error: flux.z.down should be of length 1 or a matrix with dim Nx, Ny")
    }

## check input of grid

    if (is.null(dx) && is.null(dy) && is.null(dz) && is.null(grid))
      stop("error: dx, dy, dz and grid cannot be NULL at the same time")

    gn <- names(grid)
    if (! "dx" %in% gn)
      stop("error: grid should be a list that contains 'dx' ")
    if (! "dx.aux" %in% gn)
    	stop("error: grid should be a list that contains 'dx.aux' ")
    if (! "dy" %in% gn)
      stop("error: grid should be a list that contains 'dy' ")
    if (! "dy.aux" %in% gn)
    	stop("error: grid should be a list that contains 'dy.aux' ")
    if (! "dz" %in% gn)
      stop("error: grid should be a list that contains 'dz' ")
    if (! "dz.aux" %in% gn)
    	stop("error: grid should be a list that contains 'dz.aux' ")

    if (is.null(grid$dx) || is.null(grid$dx.aux))
    	stop("error: the grid should be a list with (numeric) values for 'dx' and 'dx.aux' ")
    if (is.null(grid$dy) || is.null(grid$dy.aux))
    	stop("error: the grid should be a list with (numeric) values for 'dy' and 'dy.aux' ")
    if (is.null(grid$dz) || is.null(grid$dz.aux))
    	stop("error: the grid should be a list with (numeric) values for 'dz' and 'dz.aux' ")
    if (any(grid$dx <= 0) || any(grid$dx.aux <= 0) )
    	stop("error: the grid distances dx and dx.aux should always be positive")
    if (any(grid$dy <= 0) || any(grid$dy.aux <= 0) )
    	stop("error: the grid distances dy and dy.aux should always be positive")
    if (any(grid$dz <= 0) || any(grid$dz.aux <= 0) )
    	stop("error: the grid distances dz and dz.aux should always be positive")

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

## check input of D.grid

    if (is.null(D.x) && is.null(D.y) && is.null(D.grid))
      stop("error: D.x, D.y, and D.grid cannot be NULL at the same time")

    gn <- names(D.grid)
    if (! "x.int" %in% gn)
      stop("error: D.grid should be a list that contains 'x.int', the D values at the interfaces of the grid cells in x-direction")
    if (! "y.int" %in% gn)
      stop("error: D.grid should be a list that contains 'y.int', the D values at the interfaces of the grid cells in y-direction")
    if (! "z.int" %in% gn)
      stop("error: D.grid should be a list that contains 'z.int', the D values at the interfaces of the grid cells in z-direction")
    if (is.null(D.grid$x.int))
      stop("error: D.grid$x.int should be a list with (numeric) values")
    if (is.null(D.grid$y.int))
      stop("error: D.grid$y.int should be a list with (numeric) values")
    if (is.null(D.grid$z.int))
      stop("error: D.grid$z.int should be a list with (numeric) values")
    if (any (D.grid$x.int < 0)||any (D.grid$y.int < 0)||any (D.grid$z.int < 0))
    	stop("error: the diffusion coefficient should always be positive")

## check input of v.grid

    if (is.null(v.x) && is.null(v.y) && is.null(v.grid))
      stop("error: v.x, v.y, and v.grid cannot be NULL at the same time")

    gn <- names(v.grid)
    if (! "x.int" %in% gn)
      stop("error: v.grid should be a list that contains 'x.int', the velocity values at the interfaces of the grid cells in x-direction")
    if (! "y.int" %in% gn)
      stop("error: v.grid should be a list that contains 'y.int', the velocity values at the interfaces of the grid cells in y-direction")
    if (! "z.int" %in% gn)
      stop("error: v.grid should be a list that contains 'z.int', the velocity values at the interfaces of the grid cells in z-direction")
    if (is.null(v.grid$x.int))
      stop("error: the advective velocity v.grid$x.int should be a list with (numeric) values")
    if (is.null(v.grid$y.int))
      stop("error: the advective velocity v.grid$y.int should be a list with (numeric) values")
    if (is.null(v.grid$z.int))
      stop("error: the advective velocity v.grid$z.int should be a list with (numeric) values")

## check input of VF.grid

    gn <- names(VF.grid)
    if (! "x.int" %in% gn)
      stop("error: VF.grid should be a list that contains 'x.int'")
    if (! "y.int" %in% gn)
      stop("error: VF.grid should be a list that contains 'y.int'")
    if (! "z.int" %in% gn)
      stop("error: VF.grid should be a list that contains 'z.int'")
    if (! "x.mid" %in% gn)
      stop("error: VF.grid should be a list that contains 'x.mid'")
    if (! "y.mid" %in% gn)
      stop("error: VF.grid should be a list that contains 'y.mid'")
    if (! "z.mid" %in% gn)
      stop("error: VF.grid should be a list that contains 'z.mid'")
    if (is.null(VF.grid$x.int) || is.null(VF.grid$y.int)
     || is.null(VF.grid$x.mid) || is.null(VF.grid$y.mid)
     || is.null(VF.grid$z.mid) || is.null(VF.grid$z.int))
     stop("error: VF should contain (numeric) values")
    if (any (VF.grid$x.int < 0) || any (VF.grid$y.int < 0)
     || any (VF.grid$x.mid < 0) || any (VF.grid$y.mid < 0)
     || any (VF.grid$z.mid < 0) || any (VF.grid$z.int < 0))
      stop("error: the VF values should always be positive")

## check input of A.grid
    gn <- names(A.grid)
    if (! "x.int" %in% gn)
      stop("error: A.grid should be a list that contains 'x.int'")
    if (! "y.int" %in% gn)
      stop("error: A.grid should be a list that contains 'y.int'")
    if (! "z.int" %in% gn)
      stop("error: A.grid should be a list that contains 'z.int'")
    if (! "x.mid" %in% gn)
      stop("error: A.grid should be a list that contains 'x.mid'")
    if (! "y.mid" %in% gn)
      stop("error: A.grid should be a list that contains 'y.mid'")
    if (! "z.mid" %in% gn)
      stop("error: A.grid should be a list that contains 'z.mid'")
    if (is.null(A.grid$x.int) || is.null(A.grid$y.int)
     || is.null(A.grid$x.mid) || is.null(A.grid$y.mid)
     || is.null(A.grid$z.mid) || is.null(A.grid$z.int))
     stop("error: the VF should contain (numeric) values")
    if (any (A.grid$x.int < 0) || any (A.grid$y.int < 0)
     || any (A.grid$x.mid < 0) || any (A.grid$y.mid < 0)
     || any (A.grid$z.mid < 0) || any (A.grid$z.int < 0))
      stop("error: the A values should always be positive")

  }
## FUNCTION BODY: CALCULATIONS

## Impose boundary flux at upstream x-boundary when needed
## Default boundary condition is no gradient
  if (! is.null (flux.x.up[1])) {
    nom <- flux.x.up + VF.grid$x.int[1,,]*(D.grid$x.int[1,,]/grid$dx.aux[1] +
           (1-AFDW.grid$x.int[1,,])*v.grid$x.int[1,,])*C[1,,]
    denom <- VF.grid$x.int[1,,]*(D.grid$x.int[1,,]/grid$dx.aux[1]+
             AFDW.grid$x.int[1,,]*v.grid$x.int[1,,])
    C.x.up <- nom/denom
  }

## Impose boundary flux at downstream x-boundary when needed
## Default boundary condition is no gradient
  if (! is.null (flux.x.down[1])) {
  	nom <- flux.x.down - VF.grid$x.int[(Nx+1),,]*(D.grid$x.int[(Nx+1),,]/
            grid$dx.aux[Nx+1] + AFDW.grid$x.int[(Nx+1),,]*v.grid$x.int[(Nx+1),,])*C[Nx,,]
    denom <- -VF.grid$x.int[(Nx+1),,]*(D.grid$x.int[(Nx+1),,]/grid$dx.aux[Nx+1]+
            (1-AFDW.grid$x.int[(Nx+1),,])*v.grid$x.int[(Nx+1),,])
    C.x.down <- nom/denom
  }

# Impose boundary flux at upstream y-boundary when needed
# Default boundary condition is no gradient
  if (! is.null (flux.y.up[1])) {
    nom <- flux.y.up + VF.grid$y.int[,1,]*(D.grid$y.int[,1,]/grid$dy.aux[1] +
           (1-AFDW.grid$y.int[,1,])*v.grid$y.int[,1,])*C[,1,]
    denom <- VF.grid$y.int[,1,]*(D.grid$y.int[,1,]/grid$dy.aux[1]+
             AFDW.grid$y.int[,1,]*v.grid$y.int[,1,])
    C.y.up <- nom/denom
  }

# Impose boundary flux at downstream y-boundary when needed
# Default boundary condition is no gradient
  if (! is.null (flux.y.down[1]))  {
	  nom <- flux.y.down - VF.grid$y.int[,(Ny+1),]*(D.grid$y.int[,(Ny+1),]/
           grid$dy.aux[Ny+1] + AFDW.grid$y.int[,(Ny+1),]*v.grid$y.int[,(Ny+1),])*C[,Ny,]
    denom <- -VF.grid$y.int[,(Ny+1),]*(D.grid$y.int[,(Ny+1),]/grid$dy.aux[Ny+1]+
             (1-AFDW.grid$y.int[,(Ny+1),])*v.grid$y.int[,(Ny+1),])
    C.y.down <- nom/denom
  }

# Impose boundary flux at upstream z-boundary when needed
# Default boundary condition is no gradient
  if (! is.null (flux.z.up[1])) {
    nom <- flux.z.up + VF.grid$z.int[,,1]*(D.grid$z.int[,,1]/grid$dz.aux[1] +
           (1-AFDW.grid$z.int[,,1])*v.grid$z.int[,,1])*C[,,1]
    denom <- VF.grid$z.int[,,1]*(D.grid$z.int[,,1]/grid$dz.aux[1]+
             AFDW.grid$z.int[,,1]*v.grid$z.int[,,1])
    C.z.up <- nom/denom
  }

# Impose boundary flux at downstream z-boundary when needed
# Default boundary condition is no gradient
  if (! is.null (flux.z.down[1]))  {
	  nom <- flux.z.down - VF.grid$z.int[,,(Nz+1)]*(D.grid$z.int[,,(Nz+1)]/
           grid$dz.aux[Nz+1] + AFDW.grid$z.int[,,(Nz+1)]*v.grid$z.int[,,(Nz+1)])*C[,,Nz]
    denom <- -VF.grid$z.int[,,(Nz+1)]*(D.grid$z.int[,,(Nz+1)]/grid$dz.aux[Nz+1]+
             (1-AFDW.grid$z.int[,,(Nz+1)])*v.grid$z.int[,,(Nz+1)])
    C.z.down <- nom/denom
  }

## when upper boundary layer is present, calculate new C.x.up
  if (!is.null(a.bl.x.up[1]) & !is.null(C.x.up[1])) {
	  nom <- a.bl.x.up*C.x.up + VF.grid$x.int[1,,]*(D.grid$x.int[1,,]/
           grid$dx.aux[1] + (1-AFDW.grid$x.int[1,,])*v.grid$x.int[1,,])*C[1,,]
    denom <- a.bl.x.up + VF.grid$x.int[1,,]*(D.grid$x.int[1,,]/grid$dx.aux[1]+
             AFDW.grid$x.int[1,,]*v.grid$x.int[1,,])
	  C.x.up <- nom/denom
  }

## when lower boundary layer is present, calculate new C.x.down
  if (!is.null(a.bl.x.down[1]) & !is.null(C.x.down[1])) {
	  nom <- a.bl.x.down*C.x.down + VF.grid$x.int[(Nx+1),,]*(D.grid$x.int[(Nx+1),,]/
           grid$dx.aux[(Nx+1)] + (1-AFDW.grid$x.int[(Nx+1),,])*
           v.grid$x.int[(Nx+1),,])*C[Nx,,]
    denom <- a.bl.x.down + VF.grid$x.int[(Nx+1),,]*(D.grid$x.int[(Nx+1),,]/
             grid$dx.aux[(Nx+1)]+ AFDW.grid$x.int[(Nx+1),,]*v.grid$x.int[(Nx+1),,])
	  C.x.down <- nom/denom
  }

## when upper y boundary layer is present, calculate new C.y.up
  if (!is.null(a.bl.y.up[1]) & !is.null(C.y.up[1])) {
	  nom <- a.bl.y.up*C.y.up + VF.grid$y.int[,1,]*(D.grid$y.int[,1,]/
           grid$dy.aux[1] + (1-AFDW.grid$y.int[,1,])*v.grid$y.int[,1,])*C[,1,]
    denom <- a.bl.y.up + VF.grid$y.int[,1,]*(D.grid$y.int[,1,]/grid$dy.aux[1]+
             AFDW.grid$y.int[,1,]*v.grid$y.int[,1,])
	  C.y.up <- nom/denom
  }

## when lower y boundary layer is present, calculate new C.y.down
  if (!is.null(a.bl.y.down[1]) & !is.null(C.y.down[1]))   {
	  nom <- a.bl.y.down*C.y.down + VF.grid$y.int[,(Ny+1),]*
           (D.grid$y.int[,(Ny+1),]/grid$dy.aux[(Ny+1)] +
           (1-AFDW.grid$y.int[,(Ny+1),])*v.grid$y.int[,(Ny+1),])*C[,Ny,]
    denom <- a.bl.y.down + VF.grid$y.int[,(Ny+1),]*(D.grid$y.int[,(Ny+1),]/
             grid$dy.aux[(Ny+1)]+ AFDW.grid$y.int[,(Ny+1),]*v.grid$y.int[,(Ny+1),])
	  C.y.down <- nom/denom
  }

## when upper z boundary layer is present, calculate new C.z.up
  if (!is.null(a.bl.z.up[1]) & !is.null(C.z.up[1])) {
	  nom <- a.bl.z.up*C.z.up + VF.grid$z.int[,1,]*(D.grid$z.int[,1,]/
           grid$dz.aux[1] + (1-AFDW.grid$z.int[,1,])*v.grid$z.int[,1,])*C[,1,]
    denom <- a.bl.z.up + VF.grid$z.int[,1,]*(D.grid$z.int[,1,]/grid$dz.aux[1]+
             AFDW.grid$z.int[,1,]*v.grid$z.int[,1,])
	  C.z.up <- nom/denom
  }

## when lower z boundary layer is present, calculate new C.z.down
  if (!is.null(a.bl.z.down[1]) & !is.null(C.z.down[1]))   {
	  nom <- a.bl.z.down*C.z.down + VF.grid$z.int[,(Nz+1),]*
           (D.grid$z.int[,,(Nz+1)]/grid$dz.aux[(Nz+1)] +
           (1-AFDW.grid$z.int[,,(Nz+1)])*v.grid$z.int[,,(Nz+1)])*C[,,Nz]
    denom <- a.bl.z.down + VF.grid$z.int[,,(Nz+1)]*(D.grid$z.int[,,(Nz+1)]/
             grid$dz.aux[(Nz+1)]+ AFDW.grid$z.int[,,(Nz+1)]*v.grid$z.int[,,(Nz+1)])
	  C.z.down <- nom/denom
  }

## Calculate diffusive part of the flux
# Nx = 10, Ny = 2, Nz = 3
#  DX <- array(dim = c(Nx, Ny, Nz), dx.aux)
#  DY <- aperm(array(dim = c(Ny,Nx,Nz), dy.aux),c(2,1,3))
#  DZ <- aperm(array(dim=c(Nz,Ny,Nx),dz.aux),c(3,2,1))
  x.Dif.flux <- -VF.grid$x.int * D.grid$x.int *
                Diff3D(mbind(C.x.up, C, C.x.down, along=1), along=1)/
                array(dim=c(Nx+1,Ny,Nz),data=grid$dx.aux)
  y.Dif.flux <-  -VF.grid$y.int * D.grid$y.int *
                Diff3D(mbind(C.y.up, C, C.y.down, along=2), along=2)/
                aperm(array(data=grid$dy.aux,dim=c(Ny+1,Nx,Nz)),c(2,1,3))
  z.Dif.flux <-  -VF.grid$z.int * D.grid$z.int *
                Diff3D(mbind(C.z.up, C, C.z.down, along=3), along=3)/
                aperm(array(data=grid$dz.aux,dim=c(Nz+1,Ny,Nx)),c(3,2,1))

## Calculate advective part of the flux
  x.Adv.flux <- 0
  
  if (any(v.grid$x.int >0) ) {
    vv <- v.grid$x.int
    vv[vv<0]<-0
    x.Adv.flux <-  x.Adv.flux + VF.grid$x.int * vv * (
                 AFDW.grid$x.int * mbindl (C.x.up,  C,along=1)
                 +  (1-AFDW.grid$x.int)  * mbindr (C,C.x.down,along=1))
  }
  if (any (v.grid$x.int < 0))  {
    vv <- v.grid$x.int
    vv[vv>0]<-0
    x.Adv.flux <-  x.Adv.flux + VF.grid$x.int * vv * (
                   (1- AFDW.grid$x.int) * mbindl(C.x.up,C,along=1)
                 + AFDW.grid$x.int * mbindr(C,C.x.down,along=1))

  }
  y.Adv.flux <- 0
  if (any(v.grid$y.int >0) ) {
    vv <- v.grid$y.int
    vv[vv<0]<-0
    y.Adv.flux <-  y.Adv.flux + VF.grid$y.int * vv * (
                 AFDW.grid$y.int * mbindl(C.y.up,C,along=2)
                 + (1-AFDW.grid$y.int) * mbindr(C,C.y.down,along=2))
  }
  if (any (v.grid$y.int < 0))  {
    vv <- v.grid$y.int
    vv[vv>0]<-0
    y.Adv.flux <-  y.Adv.flux + VF.grid$y.int * vv * (
                   (1- AFDW.grid$y.int) * mbindl(C.y.up  ,C,along=2)
                 + AFDW.grid$y.int * mbindr(C,C.y.down,along=2))
  }
  z.Adv.flux <- 0
  if (any(v.grid$z.int >0) ) {
    vv <- v.grid$z.int
    vv[vv<0]<-0
    z.Adv.flux <-  z.Adv.flux + VF.grid$z.int * vv * (
                   AFDW.grid$z.int * mbindl(C.z.up,C,along=3)
                 + (1-AFDW.grid$z.int) * mbindr(C,C.z.down,along=3))
  }
  if (any (v.grid$z.int < 0))  {
    vv <- v.grid$z.int
    vv[vv>0]<-0
    z.Adv.flux <-  z.Adv.flux + VF.grid$z.int * vv * (
                    (1-AFDW.grid$z.int) * mbindl(C.z.up,C,along=3)
                 +  AFDW.grid$z.int * mbindr(C,C.z.down,along=3))

  }

  x.flux <- x.Dif.flux + x.Adv.flux
  y.flux <- y.Dif.flux + y.Adv.flux
  z.flux <- z.Dif.flux + z.Adv.flux

## Impose boundary fluxes when needed
## Default boundary condition is no gradient
  if (! is.null (flux.x.up[1]))
    x.flux[1,,]   <- flux.x.up
  if (! is.null (flux.x.down[1]))
    x.flux[dim(x.flux)[1],,] <- flux.x.down
    
  if (! is.null (flux.y.up[1]))
    y.flux[,1,]   <- flux.y.up
  if (! is.null (flux.y.down[1]))
    y.flux[,dim(y.flux)[2],] <- flux.y.down

  if (! is.null (flux.z.up[1]))
    z.flux[,,1]   <- flux.z.up
  if (! is.null (flux.z.down[1]))
    z.flux[,,dim(z.flux)[3]] <- flux.z.down

## Calculate rate of change = flux gradient     NOG DOEN
  dFdx <- - (Diff3D(A.grid$x.int*x.flux,along=1)/ A.grid$x.mid/grid$dx) / VF.grid$x.mid
  dFdy <- - (Diff3D(A.grid$y.int*y.flux,along=2)/ A.grid$y.mid/grid$dy) / VF.grid$y.mid
  dFdz <- - (Diff3D(A.grid$z.int*z.flux,along=3)/ A.grid$z.mid/grid$dz) / VF.grid$z.mid

  if (!full.output) {
    return (list (dC = dFdx + dFdy + dFdz,                  # Rate of change due to advective-diffuisve transport in each grid cell
                  flux.x.up = x.flux[1,,],                  # flux across lower boundary interface; positive = IN
                  flux.x.down = x.flux[dim(x.flux)[1],,],   # flux across lower boundary interface; positive = OUT
                  flux.y.up = y.flux[,1,],                  # flux across lower boundary interface; positive = IN
                  flux.y.down = y.flux[,dim(y.flux)[2],],   # flux across lower boundary interface; positive = OUT
                  flux.z.up = z.flux[,,1],                  # flux across lower boundary interface; positive = IN
                  flux.z.down = z.flux[,,dim(z.flux)[3]]))  # flux across lower boundary interface; positive = OUT

  } else {
    return (list (dC = dFdx + dFdy + dFdz,                  # Rate of change due to advective-diffuisve transport in each grid cell
                  C.x.up = C.x.up,                     # concentration at upper interface
                  C.x.down = C.x.down,                 # concentration at upper interface
                  C.y.up = C.y.up,                     # concentration at upper interface
                  C.y.down = C.y.down,                 # concentration at upper interface
                  C.z.up = C.z.up,                     # concentration at upper interface
                  C.z.down = C.z.down,                 # concentration at upper interface
                  x.flux = x.flux,                     # flux across at the interface of each grid cell
                  y.flux = y.flux,                     # flux across at the interface of each grid cell
                  z.flux = z.flux,                     # flux across at the interface of each grid cell
                  flux.x.up = x.flux[1,,],                  # flux across lower boundary interface; positive = IN
                  flux.x.down = x.flux[dim(x.flux)[1],,],   # flux across lower boundary interface; positive = OUT
                  flux.y.up = y.flux[,1,],                  # flux across lower boundary interface; positive = IN
                  flux.y.down = y.flux[,dim(y.flux)[2],],   # flux across lower boundary interface; positive = OUT
                  flux.z.up = z.flux[,,1],                  # flux across lower boundary interface; positive = IN
                  flux.z.down = z.flux[,,dim(z.flux)[3]]))  # flux across lower boundary interface; positive = OUT
  }
} # end tran.3D

