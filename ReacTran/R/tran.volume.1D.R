##==============================================================================
## Advective-diffusive transport in a river
##==============================================================================

tran.volume.1D <- function(C, C.up=C[1], C.down=C[length(C)],
   C.lat=C, F.up=NULL, F.down=NULL, F.lat=NULL,
	 Disp, flow=0, flow.lat=NULL, AFDW = 1, V=NULL,
   full.check = FALSE, full.output = FALSE)

{
  Rep <- function(x, length.out) 
      if (! is.null(x)) 
        return(rep(x, length.out = length.out))
      else
        return(NULL)
        
## INPUT CHECKS

  N <- length(C)
  if (N == 0)
    stop("Error: C should be a vector with numeric values")

  if (!is.list(AFDW))
    AFDW <- list(int=AFDW)

  if (!is.list(F.lat))
    F.lat <- list(mid=Rep(F.lat,length.out=N))

  if (!is.list(Disp))
    Disp <- list(int=rep(Disp,length.out=N+1))

  if (!is.list(flow))
    flow<-list(int=flow)  # check later if one value/vector

  if (!is.list(flow.lat))
    flow.lat <-list(mid=flow.lat) # check later if one value/vector

  if (!is.list(V))
    V <- list(mid=rep(V,length.out=N))

  if (!is.list(C.lat))
    C.lat <- list(mid=rep(C.lat,length.out=length(C)))

## Start full check when requested

  if (full.check) {
## check input of grids
    gn <- names(AFDW)
    if (! "int" %in% gn)
      stop("error: AFDW should be a list that contains 'int', the AFDW values at the interface of the grid cells")
    if (is.null(AFDW$int) || length(AFDW$int) == 0)
      stop("error: AFDW should contain (numeric) values")
    if (any(AFDW$int < 0)||any(AFDW$int > 1))
	    stop("error: AFDW values should always range between 0 and 1")

## check input of V
    gn <- names(V)
    if (! "mid" %in% gn)
      stop("error: V should be a list that contains 'mid'")
    if (is.null(V$mid))
	    stop("error: the argument V should contain (numeric) values")
    if (any(V$mid < 0))
	    stop("error: teh volume V should always be larger than zero")

## check input of Disp

    gn <- names(Disp)
    if (! "int" %in% gn)
      stop("error: Disp should be a list that contains 'int', the dispersion coefficient at the grid cell interfaces")
    if (is.null(Disp$int))
      stop("error: the argument Disp should be contain (numeric) values")
    if (any (Disp$int < 0))
	    stop("error: the bulk dispersion coefficient Disp should always be positive")

## check input of flow
    gn <- names(flow)
    if (! "int" %in% gn)
      stop("error: flow should be a list that contains 'int'")
    if (is.null(flow$int))
      stop("error: the argument flow should contain (numeric) values")
#    if (any (flow$int < 0) & any (flow$int > 0))
#  	  stop("error: the discharge flow cannot be both positive and negative within the same domain")

## check input of flow.lat
    gn <- names(flow.lat)
    if (! "mid" %in% gn)
      stop("error: flow.lat should be a list that contains 'mid'")
#    if (is.null(flow.lat$mid))
#      stop("error: the argument flow.lat should contain (numeric) values")

## check input of flow
    gn <- names(C.lat)
    if (! "mid" %in% gn)
      stop("error: the argument C.lat should be a list that contains 'mid'")
    if (any (C.lat$int < 0))
	    stop("error: the concentration C.lat should always be positive")

## check input of flow
    gn <- names(F.lat)
    if (! "mid" %in% gn)
      stop("error: the argument F.lat should be a list that contains 'mid'")

  } # end full.check

## FUNCTION BODY: CALCULATIONS

## Calculate the discharge in each box (using the water balance) or the upstream
## flow velocity
  ## 1. flow.lat is NULL -> estimate flow.lat based on flow gradient (diff(flow))
  if (is.null(flow.lat$mid)) {
    flow$int <- rep(flow$int,length.out=N+1)
    flow.lat$mid <- diff(flow$int)
  } else

  ## 2. flow.lat contains numeric values - check length; should be N;
  ##    flow should contain upstream flow (one value); create flows at interface
  { if (!is.vector(flow.lat$mid) || length(flow.lat$mid) == 1) flow.lat$mid <- rep(flow.lat$mid,length.out=N)

    if (length(flow.lat$mid) != N)
      stop ("flow.lat should be one number or a vector of length = N")
    if (!(length(flow$int) = 1))
      stop ("flow should be of length = 1 if flow.lat is specified")
    f1 <- flow$int[1]
    flow$int <- c(f1,f1+cumsum(flow.lat$mid))
  }

## Calculate the lateral mass input  - IF OUTPUT: C is transported instead.

  if (is.null (F.lat$mid)) {
    if (is.null (C.lat$mid))
      stop ("C.lat and F.lat cannot be both NULL")
    v1 <- flow.lat$mid
    F.lat$mid <- 0

    if (any(flow.lat$mid>0)) {
      v1[v1<0] <- 0
      F.lat$mid <- C.lat$mid * v1
    }
    if (any(flow.lat$mid < 0)) {
      ii <- which (flow.lat$mid<0)
      F.lat$mid[ii] <- C[ii]*flow.lat$mid[ii]
    }
  } else {
    C.lat$mid <- F.lat$mid / flow.lat$mid
  }
 
## Calculate diffusive part of the mass flow F

  Dif.F <- as.vector(-Disp$int*diff(c(C.up,C,C.down)))

## Calculate advective part of the mass flow F

## positive flows first

  v1 <- flow$int
	Adv.F <- 0

  if ( any (v1 > 0 )) {
    v1[v1<0] <- 0
    conc <- AFDW$int *c(C.up,C)
    if (any (AFDW$int < 1))
         conc <- conc +(1-AFDW$int)*c(C,C.down)
	  Adv.F <- Adv.F + as.vector(v1 * conc)
  } else Adv.F <- 0
  
## If there are negative flows:
  if ( any (flow$int < 0 )) {
    v1 <- flow$int
    v1[v1>0] <- 0
    conc <- AFDW$int*c(C,C.down)
	  if (any (AFDW$int < 1))
      conc <- conc +(1-AFDW$int)*c(C.up,C)
  	Adv.F <- Adv.F + as.vector(v1 * conc)
  }

## Assemble the total mass flow

  F <- as.vector(Dif.F + Adv.F)

## Impose boundary fluxes when needed
  if (! is.null (F.up))
    F[1]   <- F.up
  if (! is.null (F.down))
    F[length(F)] <- F.down
    
## Calculate rate of change = Flux gradient + lateral input

  dC <- -diff(F)/V$mid + F.lat$mid/V$mid

    if (!full.output){
    return (list (dC = dC,               # Rate of change due to advective-diffuisve transport in each grid cell
                F.up = F[1],             # Flux across upstream boundary ; positive = IN
                F.down = F[length(F)]    # Flux across downstream boundary ; positive = OUT
                ))
  } else {
    return (list (dC = dC,                 # Rate of change in the centre of each grid cells
                flow = flow$int,
                flow.up = flow$int[1],
                flow.down = flow$int[N+1],
                flow.lat = flow.lat$mid,
                F = F,                   # Flux across at the interface of each grid cell
                F.up = F[1],             # Flux across upstream boundary ; positive = IN
                F.down = F[length(F)],   # Flux across downstream boundary ; positive = OUT
				  	    F.lat = F.lat$mid)) # Flux lateral ; positive = IN
  }
}
