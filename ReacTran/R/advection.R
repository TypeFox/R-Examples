
# januari 2012: solved bug in advection.1D, when negative velocities
# changed default C.up and C.down to NULL; when NULL the previous default
# is used - but the type of boundary is set to zero-gradient

advection.1D <- function (C, C.up = NULL, C.down = NULL,
                       flux.up = NULL, flux.down=NULL,
                       v, VF=1, A =1,  dx, dt.default = 1,
                       adv.method = c("muscl","super","quick","p3","up"),
                       full.check = FALSE) {

  # number of compartments
  n <- as.integer(length(C))

  adv.method <- match.arg(adv.method)
  advmet <- pmatch(adv.method,c("muscl","super","quick","p3","up"))
  if (is.na(advmet))
    stop ("'adv.method' not known: ",adv.method)
    
  # types of boundaries and inputs
  if (!is.null(flux.up)) {
    upbnd  <- 1
    upval  <- flux.up
  } else if (is.null (C.up)) {
    upbnd  <- 3
    upval  <- C[1]
  } else {
    upbnd  <- 2
    upval  <- C.up
  }
  if (!is.null(flux.down)) {
    dwnbnd  <- 1
    dwnval  <- flux.down
  } else if (is.null (C.down)){
    dwnbnd  <- 3
    dwnval  <- C[length(C)]
  } else {
    dwnbnd  <- 2
    dwnval  <- C.down
  }

  # timestep this works only for latest version of deSolve  !
  dt <- timestep(prev=FALSE)
#  print(dt)
  if (is.nan(dt)) dt <- dt.default
  if (dt < 1e-30) dt <- dt.default   # if dt is = 0 or ~0
  if (dt > 1e30)  dt <- dt.default

  # velocity, grid sizes, volume fractions, surface areas
  v     <- rep(v, length.out = n+1)
  
  if (is.list(dx)) {
    dx.aux <- dx$dx.aux
    dx <- dx$dx
  } else { 
    dx    <- rep(dx,length.out = n)
#    dx.aux<- 0.5*(c(0,rep(dx,length.out=n))+
#                  c(rep(dx,length.out=n),0))
    dx.aux <- rep(dx,length.out = n+1)

  }
  
  if (is.list(VF)) {
    VFint <- VF$int
    VFmid <- VF$mid
  } else { 
    VFint <- rep(VF,length.out=(n+1))
    VFmid <- 0.5*(rep(VF,length.out=(n+1))[1:n]+
                  rep(VF,length.out=(n+1))[2:(n+1)])
  }
  if (is.list(A)) {
    Aint <- A$int
    Amid <- A$mid
  } else { 
    Aint  <- rep(A,length.out=(n+1))
    Amid  <- 0.5*(rep(A,length.out=(length(C)+1))[1:n]+
                  rep(A,length.out=(n+1))[2:(n+1)])
  }
  storage.mode(dx) <- storage.mode(dx.aux) <- "double"
  storage.mode(VFint) <- storage.mode(VFmid) <- "double"
  storage.mode(Aint) <- storage.mode(Amid) <- "double"
  storage.mode(v) <- "double"

  # The advection routine ...
  O<-.Fortran("advection",n,as.double(C), as.double(dt),
   dx, dx.aux, v,
   as.integer(upbnd), as.integer(dwnbnd), as.double(upval), as.double (dwnval),
   VFint, VFmid, Aint, Amid,
   as.integer(advmet), as.integer(1), as.integer(1), dy=as.double(rep(0.,n)),
   cu=as.double(rep(0.,n+1)), it=as.integer(0), package = "ReacTran")
   
  # return the "rate of change" and the fluxes
   list(dC=O$dy, adv.flux=O$cu, flux.up =O$cu[1], flux.down=O$cu[n+1], 
     it=O$it)
}

## =============================================================================
## volumetric advective transport
## =============================================================================

advection.volume.1D <- function (C, C.up = C[1], C.down = C[length(C)],
                       F.up = NULL, F.down=NULL,
                       flow, V, dt.default = 1,
                       adv.method = c("muscl","super","quick","p3","up"),
                       full.check = FALSE) {

  # number of compartments
  n <- as.integer(length(C))
  adv.method <- match.arg(adv.method)
  advmet <- pmatch(adv.method,c("muscl","super","quick","p3","up"))
  if (is.na(advmet))
    stop ("'adv.method' not known: ",adv.method)
    
  # types of boundaries and inputs
  if (is.null(F.up)) {
    upbnd  <- 2
    upval  <- C.up
  } else {
    upbnd  <- 1
    upval  <- F.up
  }
  if (is.null(F.down)) {
    dwnbnd  <- 2
    dwnval  <- C.down
  } else {
    dwnbnd  <- 1
    dwnval  <- F.down
  }

  # timestep - this works only for latest version of deSolve  !
  dt <- timestep(prev = FALSE)
  if (dt == 0)   dt <- dt.default
  if (dt > 1e30) dt <- dt.default
  
  # velocity, grid sizes, volume fractions, surface areas
  flow  <- rep(flow, length.out = n+1)
  V     <- rep(V,length.out = n)
  V.aux <- c(V[1],V)

  storage.mode(V) <- storage.mode(V.aux) <- "double"
  storage.mode(flow) <- "double"

  # The advection routine ...
  O<-.Fortran("advectvol",n,as.double(C), as.double(dt),
   V, V.aux, flow,
   as.integer(upbnd), as.integer(dwnbnd), 
   as.double(upval), as.double (dwnval),
   as.integer(advmet), as.integer(1), dy=as.double(rep(0.,n)),
   cu=as.double(rep(0.,n+1)), it=as.integer(0), package = "ReacTran")

  # return the "rate of cange" and the fluxes
   list(dC=O$dy, F=O$cu, F.up =O$cu[1], F.down=O$cu[n+1], it=O$it)
}
