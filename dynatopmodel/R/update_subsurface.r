#source("kinematic_desolve.r")
#c.route.kinematic.euler <- cmpfun(route.kinematic.euler)
################################################################################
# run inner loop Solve kinematic equation for base flow at given base flows at
# previous steps and input at this one
# Input (all)
# w               : weighting matrix
# pe              : potential evap
# tm,             : current simulation time
# ntt,            : no. inner time steps
# dt,             : main time step
# dqds,           : gradient functions
# ------------------------------------------------------------------------------
# Input (each HSU)
# ------------------------------------------------------------------------------
# flows$qbf       : base (specific) flow  at previous time step
# flows$qin       : input (total) flow at previous time step
# flows$uz        : drainage recharge from unsaturated zone (m/hr)
# groups$area : plan area
#
# stores$sd       : specific storage deficits
# ichan           : channel identifiers
# ------------------------------------------------------------------------------
# Input
# ------------------------------------------------------------------------------
# Returns (each HSU)
# ------------------------------------------------------------------------------
# flows$qbf       : specific base flow at time step, per grouping
# flows$qin       : total input (upslope) input, per grouping
# stores$sd       : updated specific storage deficits: unsat drainage and qin inputs, qbf out
# flows$qof       : updated specific overland flux per areal grouping
#
# weighting matrix
################################################################################
update.subsurface <- function (groups, flows, stores,
						  w,
                          pe=0, # potential evap
                          tm,   # current simulation time
                          ntt,  # no. inner time steps
                          dt,   # main time step
                          dqds,   # gradient functions
                          ichan=1)
{
  # save storages
  stores1 <- stores

  dtt <- dt/ntt
  # intial river flux is given by the saturation excess storage redistributed from
  Qriv <- 0

  # subsurface flows for subsurface and channels at inner time steps
  qriv.in <- matrix(0,ncol=length(ichan), nrow=ntt)
  colnames(qriv.in) <- groups[ichan,"id"]
  qriv.out <- qriv.in

  # base flow excess (per inner time step)
 # qb.ex <- matrix(0,ncol=nrow(groups), nrow=ntt)
  # record of actual evapotranpiration over each inner step
  ae.dtt <- matrix(0,ncol=nrow(groups), nrow=ntt)
  timei<-tm
  for(inner in 1:ntt)
  {
    iter <- 1
  	# apply rain input and evapotranspiration (*rates*) across the inner time step
    # note that rain and actual evap are maintained in flows
    updated <- root.zone(groups, flows, stores, pe,
                    dtt, timei, ichan)        #
    # ae is removed from root zone only - note that original DynaTM has code to remove evap
    # at max potential rate from unsat zone
    flows <- updated$flows  # includes ae and rain
  	stores <- updated$stores  #  storage excess

    # route excess flow from root zone into unsaturated zone and  drainage into water table
    updated <- unsat.zone(groups, flows, stores, dtt, timei, ichan)
  	flows <- updated$flows
    stores <- updated$stores

    # Distribute baseflows downslope through areas using precalculated inter-group
    # flow weighting matrix to give input flows for at next time step - required for ann implicit soln
    # note total input qin returned andconverted within kinematic routine
#     if(any(flows$qbf > groups$qbmax*0.75, na.rm=T))
#     {
#     #  iter <- round(20/ntt)
#    #   browser()
#     }
#     dtt.ode <- dtt/iter
#     for(i in 1:iter)
#     {
   #   message("Increasing no. iterations")
    # solution of ODE system. Now uses the Livermore solver by default
      updated <- route.kinematic.euler(groups, flows, stores, dtt,
    								ichan=ichan, w=w, time=timei,
                    dqds=dqds)

			flows <- updated$flows

     # update stores and route any excess flow
     updated <- update.storages(groups, flows, stores, dtt, ichan, tm=time)
     stores <- updated$stores

#    if(any(flows$ex>0, na.rm=T))
#    {
 #     browser()
#    }
     # stores updated by net baseflow and drainage from unsat zone across time step
#    }

  	# distribute fluxes to give new estimate for qin(t) given qbf(t) determined above
  	# base flux transferred from other areas across inner time step
  	flows$qin <- as.vector(dist.flux(groups, flows$qbf,
  												 ichan = ichan,
  												 W=w))

 		# channel flow into input
    qriv.in[inner,] <- flows[ichan,]$qin
    # actual ae at this time step
    ae.dtt[inner,] <- flows$ae
    # total excess over inner time step: sat excess surface storage and base flow excess
    # record base flow into and out of all river reaches
#    qriv.out[inner,] <- flows[ichan,"qbf"]

    flows$ex <- 0

    timei <- timei + dtt*3600
  }

  # average out total ae
  flows$ae <- colMeans(ae.dtt)  #Sums(ae.dtt)*dtt/dt

	# channel flows are rain in over time step minus evapotranspiration, which doesn't vary according
	# to root zone storage, only whether rain is falling at the time. take mean of total input across inner loop
	flows[ichan,]$qin <- colMeans(qriv.in) + (flows[ichan,]$rain-flows[ichan,]$ae)*groups[ichan,]$area

  Qriv <- colMeans(qriv.out)#+stores[ichan,]$ex

  # specific overland flow (rate)
#  flows$qof <- stores$ex/dt

	# ############################## water balance check ###########################
# 	stores.diff <- stores - stores0
# 	store.gain <- stores.diff$ex + stores.diff$srz + stores.diff$suz-stores.diff$sd
# 	net.in <- (flows$rain- flows$ae)*dtt
# 	bal <- store.gain-net.in
	# ############################## water balance check ###########################
  # return updated fluxes and storages, total discharge into river
  return(list("flows"=flows, "stores"=stores))
}

# adjust storages given updated base flow and inputs over this time step
# if max storage deficit or saturation recahed due to net inflow (outflow) then
# route the excess overland
update.storages <- function(groups, flows, stores, dtt, ichan, tm)
{
  # initially add any excess baseflow to the excess(surface) storage
  stores$ex <- stores$ex + flows$ex*dtt

  noflow <- setdiff(which(stores$sd>=groups$sd_max), ichan)
  if(length(noflow)>0)
  {
    LogEvent("SD > SDmax")  #, tm=tm)  #, paste0(groups[noflow,]$id, collapse=",")), tm=tm)
    stores[noflow,]$sd<- groups[noflow,]$sd_max
    #cat("Maximum storage deficit reached. This usually indicates pathological behaviour leading to extreme performance degradation. Execution halted")
   # stop()
    flows[noflow,]$qbf <- 0  # flows[noflow,]$qbf - (stores[noflow,]$sd - groups[noflow,]$sd_max) / dtt
  }
	# check for max saturated flow for grouping, in which case set SD to zero and
	# route excess as overland flow update storage deficits using the inflows and
	# outflows (fluxes) - note that sd is a deficit so is increased by base flow
	# out of groups and increased by flow from unsaturated zone and upslope areas
	# inc drainage from unsat zone
	# net outflow from land zone - adds to sd.
	bal <- dtt*(flows$qbf - flows$uz -flows$qin/groups$area)

	# inflow / outflow for channel is handled by the time delay histogram so storage isn't relevant here
	bal[ichan]<-0
	stores$sd <- stores$sd + bal


	# do not allow base flow to reduce storage  below zero. This should be routed overland
	saturated <- which(stores$sd<=0)  #    #setdiff(which(stores$sd<=0), ichan)
	if(length(saturated)>0)
	{
		#LogEvent(paste("Base flow saturation in zones ", paste(groups[saturated,]$tag, collapse=",")), tm=time)
		#  browser()
		# transfer negative deficit to excess store
		stores[saturated,]$ex <- stores[saturated,]$ex  -  stores[saturated,]$sd
		#
		#     # add to overland storage
		#     flows[saturated,]$ex <-
		stores[saturated,]$sd<- 0
	}

	# check channels
	return(list("flows"=flows, "stores"=stores))
}


