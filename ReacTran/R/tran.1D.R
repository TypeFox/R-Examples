
##==============================================================================
## Transport in a one-dimensional finite difference grid
##==============================================================================

tran.1D <- function(C, C.up=C[1], C.down=C[length(C)],
     flux.up=NULL, flux.down=NULL, a.bl.up=NULL, 
		 a.bl.down=NULL, D=0, v=0, AFDW=1, VF=1, A=1,
     dx, full.check = FALSE, full.output = FALSE) {


### INPUT CHECKS
  if (is.null(dx))
    stop("error: dx should be inputted ")

  N <- length(C)
  if (N == 0)
    stop("C should be a vector with numeric values")

  if (is.null(C.up))   C.up <- C[1]
  if (is.null(C.down)) C.down <- C[N]

  ## default infilling of grid parameters

  if (!is.list(AFDW)) AFDW <- list(int=AFDW)
  if (!is.list(D)) D <- list(int=D)
  if (!is.list(v)) v <- list(int=v)

  if (!is.list(VF))
    VF <- list(int=rep(VF,length.out=(N+1)),
                mid=0.5*(rep(VF,length.out=(N+1))[1:N]+
                rep(VF,length.out=(N+1))[2:(N+1)]))
  if (!is.list(A))
    A <- list(int=rep(A,length.out=(N+1)),
              mid=0.5*(rep(A,length.out=(length(C)+1))[1:N]+
              rep(A,length.out=(N+1))[2:(N+1)]))

  if (is.list(dx)) grid <- dx
  if (!is.list(dx))
    grid <- list(dx=rep(dx,length.out=N),
                 dx.aux=0.5*(c(0,rep(dx,length.out=N))+
                 c(rep(dx,length.out=N),0)))

  ## check dimensions of input arguments

  if (full.check) {

    ## check input of AFDW
    gn <- names(AFDW)
    if (! "int" %in% gn)
      stop("error: AFDW should be a list that contains 'int', the AFDW values at the interface of the grid cell ")
    if (is.null(AFDW$int))
      stop("error: AFDW is NULL, should contain (numeric) values")
    if (!is.null(AFDW$int)) {
      if (!((length(AFDW$int)==1) || (length(AFDW$int)==(N+1))))
        stop("error: AFDW should be a vector of length 1 or N+1")
    }
    if (any(AFDW$int < 0)||any(AFDW$int > 1))
	    stop("error: the AFDW should always range between 0 and 1")

    ## check input of D
    gn <- names(D)
    if (! "int" %in% gn)
      stop("error: D should be a list that contains 'int', the D values at the interface of the grid cell ")
    if (is.null(D$int))
      stop("error: D is NULL, should contain (numeric) values")
    if (!is.null(D$int)) {
      if (!((length(D$int)==1) || (length(D$int)==(N+1))))
        stop("error: D should be a vector of length 1 or N+1")
    }
    if (any(D$int < 0))
	    stop("error: the diffusion coefficient should always be positive")

    ## check input of v
    gn <- names(v)
    if (! "int" %in% gn)
      stop("error: v should be a list that contains 'int', the v values at the interface of the grid cell ")
    if (is.null(v$int))
      stop("error: the advective velocity v is NULL, should contain (numeric) values")
    if (!is.null(v$int)) {
      if (!((length(v$int)==1) || (length(v$int)==(N+1))))
        stop("error: v should be a vector of length 1 or N+1")
    }

    ## check input of VF
    gn <- names(VF)
    if (! "int" %in% gn)
      stop("error: VF should be a list that contains 'int', the area values at the interface of the grid cell ")
    if (! "mid" %in% gn)
      stop("error: VF should be a list that contains 'mid', the area at the middle of the grid cells")
    if (is.null(VF$int) || is.null(VF$mid))
      stop("error: the volume fraction VF should contain (numeric) values")
    if (!is.null(VF$int)) {
      if (!((length(VF$int)==1) || (length(VF$int)==(N+1))))
        stop("error: VF$int should be a vector of length 1 or N+1")
    }
    if (!is.null(VF$mid)) {
      if (!((length(VF$mid)==1) || (length(VF$mid)==(N))))
        stop("error: VF$mid should be a vector of length 1 or N")
    }
    if (any(VF$int < 0) || any(VF$mid < 0) || any(VF$int > 1) ||any(VF$mid > 1))
      stop("error: the volume fraction should range between 0 and 1")

    ## check input of A
    gn <- names(A)
    if (! "int" %in% gn)
      stop("error: A should be a list that contains 'int', the area values at the interface of the grid cell ")
    if (! "mid" %in% gn)
      stop("error: A should be a list that contains 'mid', the area at the middle of the grid cells")
    if (is.null(A$int) || is.null(A$mid))
      stop("error: the surface area A is NULL, should contain (numeric) values")
    if (!is.null(A$int)) {
      if (!((length(A$int)==1) || (length(A$int)==(N+1))))
        stop("error: A$int should be a vector of length 1 or N+1")
    }
    if (!is.null(A$mid)) {
      if (!((length(A$mid)==1) || (length(A$mid)==(N))))
        stop("error: A$mid should be a vector of length 1 or N")
    }
    if (any (A$int < 0) || any (A$mid < 0))
      stop("error: the area A should always be positive")

    ## check input of grid
    gn <- names(grid)
    if (! "dx" %in% gn)
      stop("error: grid should be a list that contains 'dx' ")
    if (! "dx.aux" %in% gn)
	    stop("error: grid should be a list that contains 'dx.aux' ")
    if (is.null(grid$dx) || is.null(grid$dx.aux))
    	stop("error: the grid should be a list with (numeric) values for 'dx' and 'dx.aux' ")
    if (!is.null(grid$dx)) {
      if (!((length(grid$dx)==1) || (length(grid$dx)==N)))
        stop("error: dx should be a vector of length 1 or N")
    }
    if (!is.null(grid$dx.aux)) {
      if (!((length(grid$dx.aux)==1) || (length(grid$dx.aux)==(N+1))))
        stop("error: dx.aux should be a vector of length 1 or N+1")
    }
    if (any(grid$dx <= 0) || any(grid$dx.aux <= 0) )
    	stop("error: the grid distances dx and dx.aux should always be positive")

    ## check input of boundary layer parameters
    if (!is.null(a.bl.up) & !is.null(C.up)) {
    	if (a.bl.up < 0)
        stop("error: the boundary layer transfer coefficient should be positive")
    }

    if (!is.null(a.bl.down) & !is.null(C.down)) {
	    if (a.bl.down < 0)
        stop("error: the boundary layer transfer coefficient should be positive")
    }
  } # end full.check


### FUNCTION BODY: CALCULATIONS
  if (full.output) {
    ## Impose boundary flux at upper boundary when needed
    ## Default boundary condition is zero gradient
    if (! is.null (flux.up)) {
      if (v$int[1] >= 0) {
      ## advection directed downwards
        nom <- flux.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               (1-AFDW$int[1])*v$int[1])*C[1]
        denom <- VF$int[1]*(D$int[1]/grid$dx.aux[1]+
                 AFDW$int[1]*v$int[1])
	    } else	{
      ## advection directed upwards
        nom <- flux.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               AFDW$int[1]*v$int[1])*C[1]
        denom <- VF$int[1]*(D$int[1]/grid$dx.aux[1]+
                (1-AFDW$int[1])*v$int[1])
	    }
      C.up <- nom/denom
    }
    ## Impose boundary flux at lower boundary when needed
    ## Default boundary condition is no gradient
    if (! is.null (flux.down)) {
      if (v$int[N+1] >= 0) {
        ## advection directed downwards
	      nom <- flux.down - VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
               AFDW$int[N+1]*v$int[N+1])*C[N]
        denom <- -VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1]+
                (1-AFDW$int[N+1])*v$int[N+1])
      } else {
        ## advection directed downwards
	      nom <- flux.down - VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
               (1-AFDW$int[N+1])*v$int[N+1])*C[N]
        denom <- -VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1]+
                 AFDW$int[N+1]*v$int[N+1])
      }
      C.down <- nom/denom
    }

    ## when upstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.up) & !is.null(C.up)) {
      if (v$int[1] >= 0)	{
        ## advection directed downwards
        nom <- a.bl.up*C.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
              (1-AFDW$int[1])*v$int[1])*C[1]
        denom <- a.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
                 AFDW$int[1]*v$int[1])
	    } else	{
        ## advection directed upwards
        nom <- a.bl.up*C.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               AFDW$int[1]*v$int[1])*C[1]
        denom <- a.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
                (1-AFDW$int[1])*v$int[1])
	    }
      C.up <- nom/denom
    }

    ## when downstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.down) & !is.null(C.down)) {
      if (v$int[N+1] >= 0)	{
        ## advection directed downwards
       nom <- a.bl.down*C.down + VF$int[N+1]*(D$int[N+1]/
              grid$dx.aux[N+1] + (1-AFDW$int[N+1])*v$int[N+1])*C[N]
       denom <- a.bl.down + VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
                AFDW$int[N+1]*v$int[N+1])
     	} else	{
        ## advection directed upwards
        nom <- a.bl.down*C.down + VF$int[N+1]*(D$int[N+1]/
               grid$dx.aux[N+1] + AFDW$int[N+1]*v$int[N+1])*C[N]
        denom <- a.bl.down + VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
                 (1-AFDW$int[N+1])*v$int[N+1])
	    }
      C.down <- nom/denom
    }

    ## Calculate diffusive part of the flux

    dif.flux <- as.vector(-VF$int*D$int*diff(c(C.up,C,C.down))/
                          grid$dx.aux)
    adv.flux <- rep(0,length.out=length(dif.flux))

    ## Add advective part of the flux if needed
      if (any (v$int > 0)) {   # advection directed downwards
       vv <- v$int
       vv[v$int<0] <- 0
	     conc <- AFDW$int*c(C.up,C)
	     if (any (AFDW$int < 1))
         conc <- conc +(1-AFDW$int)*c(C,C.down)
	     adv.flux <- adv.flux + as.vector(VF$int*vv*conc)
    }
    if (any (v$int < 0)) {   # advection directed upwards
       vv <- v$int
       vv[v$int>0] <- 0
       conc <- AFDW$int*c(C,C.down)
	     if (any (AFDW$int < 1))
         conc <- conc +(1-AFDW$int)*c(C.up,C)
	     adv.flux <- adv.flux + as.vector(VF$int*vv*conc)
    }

    flux <- dif.flux + adv.flux

  } else { # not full.output

    ## when upstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.up) & !is.null(C.up))  {
      if (v$int[1] >= 0) { # advection directed downwards
        nom <- a.bl.up*C.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               (1-AFDW$int[1])*v$int[1])*C[1]
        denom <- a.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               AFDW$int[1]*v$int[1])
	    } else	{  # advection directed upwards
        nom <- a.bl.up*C.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
               AFDW$int[1]*v$int[1])*C[1]
        denom <- a.bl.up + VF$int[1]*(D$int[1]/grid$dx.aux[1] +
                 (1-AFDW$int[1])*v$int[1])
	    }
      C.up <- nom/denom
    }

    ## when upstream boundary layer is present, calculate new C.up
    if (!is.null(a.bl.down) & !is.null(C.down)) {
      if (v$int[N+1] >= 0)	{   # advection directed downwards
        nom <- a.bl.down*C.down + VF$int[N+1]*(D$int[N+1]/
               grid$dx.aux[N+1] + (1-AFDW$int[N+1])*v$int[N+1])*C[N]
        denom <- a.bl.down + VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
                 AFDW$int[N+1]*v$int[N+1])
	    } else	{  # advection directed upwards
        nom <- a.bl.down*C.down + VF$int[N+1]*(D$int[N+1]/
                 grid$dx.aux[N+1] + AFDW$int[N+1]*v$int[N+1])*C[N]
        denom <- a.bl.down + VF$int[N+1]*(D$int[N+1]/grid$dx.aux[N+1] +
                (1-AFDW$int[N+1])*v$int[N+1])
	    }
      C.down <- nom/denom
    }

    ## Calculate diffusive part of the flux
	  flux <- as.vector(-(VF$int)*D$int*
                      diff(c(C.up,C,C.down))/grid$dx.aux)
    ## Add advective part of the flux if needed
    if (any (v$int > 0)) {   # advection directed downwards
       vv <- v$int
       vv[v$int<0] <- 0
	     conc <- AFDW$int*c(C.up,C)
	     if (any (AFDW$int < 1))
         conc <- conc +(1-AFDW$int)*c(C,C.down)
	     flux <- flux + as.vector(VF$int*vv*conc)
    }
    if (any (v$int < 0)) {   # advection directed upwards
       vv <- v$int
       vv[v$int>0] <- 0
       conc <- AFDW$int*c(C,C.down)
	     if (any (AFDW$int < 1))
         conc <- conc +(1-AFDW$int)*c(C.up,C)
	     flux <- flux + as.vector(VF$int*vv*conc)
    }

  }

  if (! is.null (flux.up)) flux[1] <- flux.up
  if (! is.null (flux.down)) flux[N+1] <- flux.down

    
## Calculate rate of change = Flux gradient       
  dC <- -diff(A$int*flux)/A$mid/VF$mid/grid$dx

  if (!full.output){
    return (list (dC = dC,                   # Rate of change due to advective-diffusive transport in each grid cell
  					flux.up = flux[1],               # Flux across lower boundary interface; positive = IN
	  				flux.down = flux[length(flux)])) # Flux across lower boundary interface; positive = OUT

  } else {
  	return (list (dC = dC,                   # Rate of change due to advective-diffusive transport in each grid cell
	  				C.up = C.up,                     # Concentration at upper interface
		  			C.down = C.down,                 # Concentration at lower interface
			  		dif.flux = dif.flux,             # Diffusive flux across at the interface of each grid cell
				  	adv.flux = adv.flux,             # Advective flux across at the interface of each grid cell
					  flux = flux,                     # Flux across at the interface of each grid cell
  					flux.up = flux[1],               # Flux across lower boundary interface; positive = IN
	  				flux.down = flux[length(flux)])) # Flux across lower boundary interface; positive = OUT
  }

} # end tran.1D

