
##==============================================================================
## Calculates the advective velocities assuming steady state compaction
##==============================================================================

setup.compaction.1D <- function(v.0 = NULL, v.inf = NULL,
	 por.0, por.inf, por.grid) {

## check input
  gn <- names(por.grid)
  if (! "mid"  %in% gn)
    stop("error in setup.prop: porosity should be a list that contains mid")
  if (! "int"  %in% gn)
    stop("error in setup.prop: porosity should be a list that contains int")

  if (is.null(v.0) && is.null(v.inf))
    stop("error in setup.advection: either the sedimentation velocity <v.0> or the burial velocity <v.inf> should be specified")

## calculate velocities at infinity

  if (is.null(v.inf)) { # sedimentation velocity is specified
  	v.factor <-  v.0*(1-por.0)
  	v.inf <- v.factor/(1-por.inf)
  	u.inf <- v.inf
  	u.factor <-  u.inf*por.inf
  } else { # burial velocity is specified
  	v.factor <-  v.inf*(1-por.inf)
	  v.0 <- v.factor/(1-por.0)
  	u.inf <- v.inf
	  u.factor <-  u.inf*por.inf
	}

  v.mid <- v.factor/(1-por.grid$mid)
  v.int <- v.factor/(1-por.grid$int)

  u.mid <- u.factor/por.grid$mid
  u.int <- u.factor/por.grid$int

  Res <- list(u = list(mid=u.mid,int=u.int), v=list(mid=v.mid,int=v.int))
  class(Res) <- "compaction.1D"
  return(Res)
}
