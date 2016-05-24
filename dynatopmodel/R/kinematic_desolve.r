# downslope distribution function resulting from exponential transmissivity
# assumption
dqdt <- function(q, t, w, m, r, a)
{
#	dddt <- q- (t(w)%*%(a*q))/a-r
#	dqdt <- dddt*-q/m
	A <- diag(1/a) %*% t(w) %*% diag(a) -identity.matrix(nrow(w))
  # throttle input
  #q <- pmax(q, qmax)
	res <- q/m * (A %*%q + r)
    # impose
	return(res)
#	return(res)
}

dqdt.ode <- function(t, q, parms, qmax,
                      ...)
{
 	w <- parms$w
 #	m <- parms$m
 	r <- parms$r
 	a <- parms$groups$area
    fun <- parms$dqds
    S <- parms$S
 #   qmax <- parms$qmax

	#	dddt <- q- (t(w)%*%(a*q))/a-r
	#	dqdt <- dddt*-q/m
	A <- diag(1/a) %*% t(w) %*% diag(a) -identity.matrix(nrow(w))

	res <- -fun(q, S, params=parms$groups) * (A %*%q + r)

  # impose maximum q and ensure non-negative
  #res <- pmax(pmin(res, parms$groups$qbmax, na.rm=T),0)
	return(list(res))
	#	return(res)
}

require(deSolve)

# ******************************************************************************
# kinematic.r kinematic wave routes flow downslope through groupings from
# 4-point solution using input flows and output flows from previous and current
# time steps
# ==============================================================================
# Inputs
# ------------------------------------------------------------------------------
# time step dt in hrs
# flowst1: group flows and storages at previous time step
# flows: groups flows at current time step- Qbf to be determined
# w: time-stepping weighting coefficient: zero for a totally implicit solution
# (depends only on flows for previous steps). w=0.5
# niter: max number of iterations in iterative scheme
# ==============================================================================
# Returns
# ------------------------------------------------------------------------------
# updated storage deficts, estimates for base flows at this time step
# ==============================================================================
# References
# ------------------------------------------------------------------------------
# Beven and Freer (2001). A Dynamic TOPMODEL
# Beven (1981). Kinematic subsurface stormflow
# Li el al (1975).  Li, Simons and Stevens 1975 - Nonlinear Kinematic Wave
# Approximation for Water Routing
# Beven (2012). Rainfall runoff modelling, chapter 5 pp141-150, pp.180-183
# ******************************************************************************
# note: we now exclude lateral input from land hsus from the flux inpu
route.kinematic.euler <- function(groups,
                            flows,        # fluxes at previous time step (prediction time)
                            stores,       # current storage
                            dtt,
                            ichan,
                            time,
							              w,
                            nstep=1,
				                    dqds,
              method="lsoda"              # Livermore solver
)
{

  # UZ: recharge is drainage from unsaturated into saturated zone - assumed constant over time steps
  r <- flows$uz #

 	qb0 <- flows$qbf
 	qb0[ichan] <- 0
 	res <- ode(y=qb0, times=seq(0, dtt, length.out=nstep+1),
                func=dqdt.ode,
                method=method,
                parms=list(w=w, r=r,
                           groups=groups,
#                            a=groups$area,
#                            m=groups$m,
#                            z.drain=groups$z.drain,
#                            ln_t0_plus=groups$
                           dqds=dqds,
                           S=stores$sd))
  # final row gives baseflows at each stage
 	qbf <- res[nstep+1,-1]

  res <- matrix(res[-1,-1], nrow=nstep)
  # excess base flow
  ex <- t(apply(res, MARGIN=1, function(x)ifelse(x>groups$qbf_max,
                                                    x-groups$qbf_max, 0)))
  # remove and add to saturation excess
  res <- res-ex

  # specific storage (defict) change over this step, assumming the drainage constant
#  sd.add <- -t(apply(res[,-1], MARGIN=1,
#                     function(x)((x*groups$area)%*%w-x*groups$area)/groups$area-r)*dtt/nstep)

#  sd.add <- t(apply(res, MARGIN=1,
#                   function(x)x-r-flows$qin/groups$area)*dtt/nstep)

  flows$ex <- colMeans(ex)

 # limit output to storage available (although more available from upslope)
#  qbf <- pmin(qbf, stores$sd-groups)
  # max saturated flux is calculated from topography of groups, setting SD=0
  # if we exceed this then set the flux to that value
# throttle outlet discharge to max given gradient and conductivity and ensure >0
 flows$qbf <- qbf
# stores$sd <- stores$sd + colSums(sd.add)
#stores$ex <- stores$ex +colSums(ex*dtt/nstep)
#    flows$ex <- pmax(qbf-groups$qbmax, 0, na.rm=T)
#    flows$qbf <- pmax(qbf-flows$ex,0)
  #  flows$qbf[ichan] <- 0
# updated base flow. storage calculated in update.storages
  return(list("flows"=flows, "stores"=stores))
}





