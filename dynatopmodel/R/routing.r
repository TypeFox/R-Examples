# ################################################################################
# routing of subsurface,  overland flow to channel and channel flows to the
# catchment outlet
###############################################################################


#-------------------------------------------------------------------------------
# input
# groups      : average effective overland
# w                 : flux distribution matrix
# ichan             : channel identifiers
# vof               ; overall overland flow velocity
#-------------------------------------------------------------------------------

get.routing <- function(groups,
                        w,
												routing,   # channel routing table, first column is average distance to outlet, second the proportion of flow entering at this upstream distance
                        ichan,      # channel identifiers
                        dt       # time step (hr)
                       )
{
	if(length(routing) == 0)
	{
		# everthing moved to outlet in one time step
		return(matrix(1))
	}
	# apply average channel wave velocity??? need a reach" table
	vchan <- mean(groups[ichan,]$vchan, na.rm=T)
  chan.dists <- routing[,1]
  groups <- groups[-ichan,]

  # time delay for flow once in channel(s) until reaching outlet
  # assumming a constant channel velocity
  dtichan <- round(chan.dists/vchan+0.5) # vector of length nchan

  # number of outlet considered
  nout <- ncol(routing)-1 #length(ichan)

  # time shift (in simulation steps) for channel flow to reach outlet
	# create a matrix that can hold delays up to the maximum time
  nshift.chan <- max(dtichan)
  chan.routing <- matrix(ncol=nout,
                        nrow= nshift.chan) # enough to distribute all channel flows into future time steps

  chan.routing[] <- 0

  for(iout in 1:nout)
  {
    props <- routing[,iout+1]
  	for(i in 1:length(props))
  	{
  		# sum proportions of channel flow entering at future time step
  		chan.routing[dtichan[i], iout] <- chan.routing[dtichan[i], iout]+props[i]
  	}
    # route all flow to a particular future time step, but could distribute between steps for longer intervals and reach lengths
   # chan.routing[dtichan[chan], chan] <- routing[,2] # second column is proprtion of flow delayed by this time step
  }


  # mulitply the channel flows at a particular time step by the cahhnel routing matrix to get a time-shifted version
  # rewlative to current stimulation step
  # for two channel system this could be something like
  # 0   0
  # 1   0
  # 0   1
  # 0   0
  # ..
  # of flows moved through the putlet that can be added to existing river flows

  return(chan.routing)
}


# route both overland flow and channel flows to outlet
route.channel.flows <- function (groups,
								flows,  # updated flows
                                delays = NULL, # precalculated time delay matrix, 1 column per channel
								w,        # flux distribution matrix
								Qr,      # time series of simulated total discharge at catchment outlet at simulation time
								it,     # current time step
                dt,      # time step in hours
								ichan=1,
								i.out=1,
								chan.dist.matrix=NULL  # option distribution of land
                         		)      # channel id(s)
{
	if(is.null(delays))
	{
		#  total flux transferred out of channel across time step
		Qr[it,] <- flows[ichan[i.out],]$qbf *groups[ichan[i.out],]$area
	}
	else
	{
	#	if(sum(flows[ichan,]$qof)>0){browser()}
		# input flows into channel: subsurface & surface (overland), converted to a total flow
		qchan <- flows[ichan,]$qin + flows[ichan,]$qof*groups[ichan,]$area

		if(!is.finite(qchan))
		{
			warning("Bad in channel flow ")
		}
		else if(sum(qchan)>0)
  	{
      if(!is.null(chan.dist.matrix))
      {
        Qchan <- w[,ichan] * flows$qbf*groups$area
        # split input into channel according to distribution matrix
        Qreach <- t(Qchan[-ichan]) %*% chan.dist.matrix
        Qr[it,2:ncol(Qr)] <- Qreach
      }

  		# determine fluxes at future times passing through each outlet reach specified. colsums may not all = 1
      # as much of the flow will not pass through points upstream. The path length analysis will have taken this account
      # when constructing the distance table supplied vi athe routing parameter
  		nout <- ncol(delays)   # number of outlet reaches specified
  		q.shift <- delays * qchan
  		# how many time steps will it take all the flux entering the channel to reach the outlet (x dt to get actual times)
  		n.shift <- nrow(delays)
       # add time shifted fluxes for every channel to the time series. (flow delayed by at least one time step?)
	    Qr[it:(it+n.shift-1),1:nout] <- Qr[it:(it+n.shift-1),1:nout] + q.shift
  	}

  }

  return(Qr)
}


# route.channel.flux.to.outlet <- function(Qreach, drn, pt, dt)
# {
#   # Qreach is total flux over time step into each reach
#   # specific flux (per m) per hour
#   reach.lengths <- gLength(drn, byid=T)
#   qreach <- t(apply(Qreach, MARGIN=1, function(x)x/reach.lengths))
#   qreach <- xts(qreach, order.by=index(run$Qr)[1:nrow(Qreach)])
#
#
# }





