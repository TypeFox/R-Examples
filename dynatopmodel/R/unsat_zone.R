
# Route any excess storage (effective rainfall) pex from the root zone into the
# unsat zone and calculate vertical drainage out into the water table
unsat.zone <- function(groups, flows, stores, dt, time, ichan)
{
  # initialise flow into saturated zone as zero
  flows$uz <- 0
  flows$exus <- 0
  flows$ex <- 0
  # pex is excess drainage (rate) from root zone (over time step): add to unsat storage suz.
  stores$suz <- stores$suz + flows$pex*dt # check that this is total over time step
   
  # route excess drainage (can't be stored by storage deficit) into excess store
  unsatExcess <- which(stores$suz >= stores$sd & stores$sd>0)     #which( & stores$suz>stores$sd)
  if(length(setdiff(unsatExcess, ichan))>0)
  {
#    unsat.land <- setdiff(unsatExcess, ichan)
#    if(length(unsat.land)>0)
#    { LogEvent(paste0("sat. excess flow in ", paste(groups[unsat.land,]$tag, collapse=",")),tm=time)}
  	#reached max storage in usz
  	stores[unsatExcess,]$suz <- stores[unsatExcess,]$suz - stores[unsatExcess,]$sd   # exus is routed overland stores[unsatExcess,]$sd	
  	# allocate sufficient drainage to fill remaining deficit and route the excess overland
  	stores[unsatExcess,]$sd <- 0
  	flows[unsatExcess,]$uz <- stores[unsatExcess,]$sd/dt
    # allocate storage to excess that exceeds that required to fill the storage deficit
  	stores[unsatExcess,]$ex <-  stores[unsatExcess,]$ex + stores[unsatExcess,]$suz
    # remove unsat storage
  	stores[unsatExcess,]$suz <- 0
  }  

  #   otherwise unsaturated regions with +ve unsat storage drain into sat zone as
  #   enough storage left to accommodate
  draining <- which(stores$sd>0 &stores$suz >0)  # & groups$td>0)
  if(length(draining)>0)
  {  	
    #browser()
    # assume a proportion of the storage flows into the unsat zone controlled by
    # unsaturated time delay for this group td (time in hrs per vertical depth of
    # infiltration; 1/td = infiltration velocity). larger storage = higher drainage rate
  	# lower sd = higher soil moisture content = higher "
    # see beven and wood (1983), Beven (2012)
    uz <- stores[draining,]$suz /(stores[draining,]$sd * groups[draining,]$td)  # m/hr
    
    if(any((dt*uz)>stores[draining,]$suz))
    {
 #     browser()
    }
    # cap drainage so that storage never falls below zero "channel" hsus have
    # very low time delay so all the storage drain into the SZ in one time step
    drainage <- pmin(dt*uz,stores[draining,]$suz)
        
	#     reduce storage by drainage out of zone over time step, not allowing this 
	#     to fall below zero
    stores[draining,]$suz <- stores[draining,]$suz-drainage
    # 
    # vertical drainage into saturated zone ("recharge") - will be routed in kinematic wave routine
    # river channels will also be recharged along their lengths by baseflow from the land - also
    # overland flow

    flows[draining,]$uz <- drainage/dt # recharge rate into saturated zone
  }
  
  # input to river is rainfall - evap (can be -ve). note that pex is a storage
  flows[ichan,]$uz <- flows[ichan,]$pex
  
  # remove rainfall input as has been allocated either to excess or recharge
  flows$pex <- 0
        
  return(list("flows"=flows, "stores"=stores))    
}

# return a list of saturated zones (zero storage deficit), excluding the channel
GetSaturated <- function(stores, ichan)
{	
	stores[ichan,]<-0.1
  sat <- which(stores$sd<=0)
  return(sat)
}