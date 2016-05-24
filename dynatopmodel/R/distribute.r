# ==============================================================================
# Distribute flux from areal group to  downslope
# areas according to weighting matrix
# ==================================================
# groups$area: plan areas for each group
# flows$qbf: existing specific base flow for grouping
# stores : output saturated zone (base) flow for each areal subdivision
# W: weightings for flows out of the subdivisions into other 
# dtt: time step (likely to be subdivison of main time step)
# ==================================================
# returns
# flows$qin: updated (total) input flux for groupings
# ==================================================
dist.flux <- function(groups, qbf, 
                      W,  # flow distribution matrix 
                      #dtt, # time step, inner or outer (needed?)
                      ichan=1)  # channel identifiers
{  
    # total baseflow
    QBF <- qbf*groups$area  # total base flow
    # ensure river fluxes doesn't get recycled to itself (why?)
    diag.chan <- matrix(rep(ichan, 2), ncol=2)
    W[diag.chan]<- rep(0, length(ichan)) #
    
    # all land - land, land-riv and riv - riv transfers
    # distribute base flux to inputs of other HSUs; note, total (not specific) discharge
    #flows$qin <-    
    return(QBF %*% W)  # as vector?
    
}


dist.ex.dummy <- function(groups, flows, stores, w, dt, ichan=1, nstep=10, debug=F)
{
	ex <- stores$ex
	# transfer everything into the river
	if(any(ex>0) )
	{
		# shift all flux into channel
		ex.chan <- sum(ex*groups$area)
		ex[] <- 0
		# route everything into the channel all at once
		# convert to specific flow
		ex[ichan] <- ex.chan/sum(groups$area[ichan])

		
	}
	return(ex)
}
	

dist.ex.2 <- function(groups, flows, stores, w, dt, ichan=1, nstep=10, debug=F)
{
	ex <- stores$ex
	if(any(ex>0) )
	{
		dtt <- dt/nstep
		#		ex[ichan]<- 0 # 
		#	if(sum(ex[ichan])>0) {stop("excess channel storage!")}
		w.in <- ex*groups$area
		# channel flux remains there until return		
		groups[ichan,]$vof <- 0
		ex.dtt <- ex
		flows.dtt <- flows
		for(i in 1:nstep)
		{ 
			# total flow out of all groups, limited to storage available
			flows$qbf <- pmin(ex*groups$vof, ex/dtt)
			
			# total flow into each group from elsewhere (including recycled flux)
			# channel recyles flux into itsel
			flows$qin <- dist.flux(groups, flows$qbf, w, ichan)
			
			flows.dtt <- rbind(flows.dtt, flows)
			
			# difference in specific flux over time step
			ex <- ex + dtt*( flows$qin/groups$area- flows$qbf)
			# ex reduced below zero: 
			ex.dtt <- rbind(ex.dtt, ex) 
		}
		
		# convert to a specific storage
		#	ex <- ex/groups$area
		# water balance check
		w.out <- ex*groups$area
		bal <- sum(w.out-w.in)
		if(bal > 0.01){warning("Water balance check error")}
	}
	
	return(ex)  
}





# Distribute rainfall (pe) input for this time step and across HRU groupings 
# ==============================================================================
# Inputs
# ==============================================================================
# rain: rainfall or pe in time step
# ==============================================================================
# returns
# updated input flows at time step
# adding recharge to river hsus 
# ==============================================================================

allocate.PE <- function(pe, groups, it)
{
    # currently as for rainfall 
    # partition rain equally across groups 
    return(allocate.rain(pe, groups, it))  
}

is.datetime <- function(obj)
{
    return(inherits(obj, "POSIXt"))
    
}

# return a vector with the correct rainfall allocated from multiple gauges to appropriate area,
# identified by gauge.id 
# this needs optimisation: the rainfall vectors could be allocated at initialisation, one for each
# HRU (which might lead to them being duplicated...) and the index automatically 
allocate.rain <- function(rain, groups, it)
{
    # can index by time instead of integer (performance implications?)
    if(is.datetime(it)){it <- which(index(rain)==it)}
    # 	if(all(colnames(rain)==groups$id))
    # 	{
    # 		# already split the rainfall input up by response unit
    # 		rain.dist <- rain[it,]   
    # 	}
    # 	else
    # 	{	
    # 	
    # 	  if(it>nrow(rain))
    # 	  {
    # 	    warning("Rainfall record truncated at ", tail(index(rain), 1))
    # 	    it <- nrow(rain)
    # 	  }
    # 	
    #ncol <- ncol(rain)
    
    rain.dist <- rep(0, nrow(groups))
#     if(length(it)==0)  # || it>nrow(rain))
#     {
#         warning("Exceeded data range in allocate rainfall, returning zero for all groups")
#         return()
#     }
    # gauge id, limited to nunber of data columns	
    try(rain.dist <- rain[it,pmin(groups$gauge.id, ncol(rain))])   
    
    return(as.vector(rain.dist))
}

