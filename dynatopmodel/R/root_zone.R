

# root zone storage, excess flow into unsaturated zone and actual evapotranspiration
# time step dt in hrs
# ==================================================
# Inputs
# --------------------------------------------------
# time step dt in hrs
# group flows and storages
# p.e and rainfall over thjis time step
# ==================================================
# Returns
# --------------------------------------------------
# updated root zone storage, excess flow into unsaturated zone and actual evapotranspiration
# precipitaion excess, if any, from root zone into unsaturated zone
# ==================================================

# areal input rainfall = flows$rain
# pe = potential evapotranspiration per unit time (mm/hr)
# flow: current root zone storage
root.zone <- function(groups, flows, stores, pe, dt, time,
                     ichan)
{
  # precipitation excess initially zero
  flows$pex <- 0
  rain <- flows$rain
  # extract land hsus
 # iland <- setdiff(1:nrow(groups), ichan)
  # determine storage change across this time step
  # This is limited by the max rz storage.
  # rainfall and pe in storage units per unit plan area per hour
  # if rain above a certain threshold then p.e. is zero

  # add rainfall input over this time step  rain is an hourly *rate*
  # any excess flow (from e.g surface storage) will have already been added
  stores$srz <- stores$srz + dt*rain
  # any rainfall on channel immediately becomes storage excess
#  stores[ichan, ]$ex <- rain[ichan]
  # only consider evap if rainfall minimal (Beven), there is some actual evap, and also moisture to lose!
 # evaporating <- which(rain < 0.001 & stores$srz > 0 & pe>0)  # & groups$srz_max>0)
  #fact<- rep(0, nrow(groups))
 # pe[which(rain > 1e-6)]<- 0


#  fact[ichan] <- 1
  # actual ET for land groups is calculated using ratio of storage to max storage
  fact <- stores$srz/groups$srz_max

  # pe (and ae) is a rate not absolute amount. (Evapotranspiration only during dry periods)
  ae <- pe * fact# * as.numeric(rain < 1e-8) #* dt

	# for evap from land units, limit to amount of storage remaining in the root zone - ensures that store doesn't fall below zero.
#	ae[-ichan] <- pmin(ae[-ichan], stores[-ichan,]$srz/dt)  # convert evap rate to amount over (inner) time step
	ae <- pmin(ae, stores$srz/dt)  # convert evap rate to amount over (inner) time step

  # remove ae from the root zone, making sure that the storage never drops below zero
  stores$srz <- stores$srz - ae*dt    #pmax(, 0)

	# evaporation removed at maximum allowed from the channel
	ae[ichan] <- pe[ichan]

	# for channel, rain -> "root zone" -> excess
#	stores[ichan,c("srz", "suz")] <- 0

	#
# Route any excess over max storage into unsaturated zone - ignore channel and reservoirs
  full <- setdiff(StoresFull(groups, stores), ichan)
  if(length(full)>0)
  {
  #  cat(paste(time, ": root zone(s) full "), groups[full,]$id, "\n")
    # storage excess is routed into the unsat zone where it drains into water table
    flows[full,]$pex <- (stores[full,]$srz -groups[full,]$srz_max)/dt   # note pex is a rate
    # root zone full
    stores[full,]$srz <- groups[full,]$srz_max
  }
  #

  flows$ae <- ae
  flows[ichan] <- NA
  return(list("flows"=flows, "stores"=stores))
}

# return indexes of full & non-empty root zones
StoresFull <- function(groups,stores)
{
  return(which(stores$srz>=groups$srz_max)) #  & stores$srz>0))
}

StoresEmpty <- function(groups, stores)
{
  return(which(stores$srz<=0 & groups$SRmax>0))
}
