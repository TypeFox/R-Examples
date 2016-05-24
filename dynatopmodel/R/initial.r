#source("update.fluxes.r")
#source("DistFlux.r")

################################################################################
# initialisation and input checks for Dynamic TOPMODEL
#-------------------------------------------------------------------------------
ResetFluxes<- function(flows, ichan)
{
  # tidy up
  # remove base flows from channel  now they have been allocated
#  flows[ichan, ]$qin <- 0
  # remove overland flows - also allocated in RouteFlow routine
  flows[, c("qof","pex", "ex", "exus")]<- 0
  return(flows)
}

ResetStorages <- function(stores, ichan)
{
  # tidy up

  # remove overland flows - also allocated in RouteFlow routine
  stores[, "ex"]<- 0
  return(stores)
}

# check that tabular object has given column names
check.cols <- function(obj, names, stop=T)
{
  msg <- NULL
  for(nm in names)
  {
    if(!(nm %in% colnames(obj)))
    {
      msg <- paste(msg, "Missing column name ", nm, "\n")
    }
  }
  if(length(msg)>0 & stop){stop(msg)}
  return(msg)
}

check.and.fix.input <- function(envir)
{

  errs <- NULL
  warns <- NULL
  if(is.null(envir$qt0) | !is.finite(envir$qt0))
  {
    warns <- c(warns, "Initial specific discharge qt0 required, setting to 1e-5")
    envir$qt0 <- 1e-5
  }
  if(is.null(envir$qmax) | !is.finite(envir$qt0))
  {
    warns <- c(warns, "Max displayed specific discharge required: setting to 5/1000")
    envir$qmax <- 5/100
  }
  if(!(is.data.frame(envir$groups)|is.matrix(envir$groups))){
    errs <- c(errs, "Areal grouping information required as data frame or matrix")
  }
  ngroups <- nrow(envir$groups)

  if(!is.matrix(envir$weights)){
    errs <- c(errs, "Flow weightings - matrix required")
  }
  if(nrow(envir$weights)!= ncol(envir$weights))
  {
    errs <- c(errs, "Square input weighting matrix required")
  }
  if(nrow(envir$weights)!=ngroups | ncol(envir$weights) != nrow(envir$weights))
  {
    errs <- c(errs, "Flow weighting matrix incorrect size - require square matrix of side equal to # areal groupings")
  }
  # ensure we have the expected col names and minimal data requirements
  check.cols(envir$groups, c("id", "area", "atb.bar"))

  if(length(envir$rain)==0)
  {
    errs <- c(errs, "Rainfall input required. Check start / end dates")
  }
  else
  {
    # check that input rain and pe (output)
    if(nrow(envir$rain) != nrow(envir$pe))
    {
      min.len <- min(nrow(envir$rain), nrow(envir$pe))
      warns <- c(warns, "rainfall and pe series differ in length, trimming to shortest")
      envir$rain <- envir$rain[1:min.len,]
      envir$pe <- envir$pe[1:min.len,]
    }
  }
  # each row of weighting matrix is total flow out of the group so should add to 1
  dist.sum <- round(rowSums(envir$weights), 1)
  not.unity <- which(dist.sum !=1)
  if(length(not.unity>0))
  {
    envir$weights <-normalise.rows(envir$weights)
  #  warns <- c(warns, paste("Group ", groups[not.unity]$tag, " weights should add to 1", sep=""))
    # normalise
   # envir$weights[not.unity,] <- normalise.vector(weights[not.unity,])
  }

  if(length(warns)>0){lapply(warns, warning)}
  if(length(errs)>0)(stop(errs))

#      cur.env <-  environment()
#      # copy values to the specified environment, typicaly the  current environment in
#      # the calling context
#      if(!is.null(calling.env))
#      {
#        for(n in ls(cur.env, all.names=TRUE))
#        {
#          # assign name-value in this envir to that specified
#          assign(n, get(n, cur.env), calling.env)
#        }
#      }
  return(envir)
}



# determine base flows assuming input output mass balance in steady state
init.base.flows <- function(groups,
                          W,     # flux distribution matrix
                          r,     # specific recharge e.g rainfall and discharge (assummed equal)
                          ichan)
{
  # across each time step
  # input:
  # external recharge e.g. rain, applied equally over each HSU, entering water
  # table as gravity drainage from unsaturated zone
  # base flow distributed out of other HSUs
  # Qint1 <- dt*(r + t(flows$qbf) %*% W)*Ag
  # output at next time:
  # base flow
  # effective evap (ignore)
  # channel discharge - ignore for time being by assumming everything distributed to
  # channel reaches outlet within one step
  # Qout <- dt*flows$qbf*groups$area
  # in steady state these must be equal. Equating
  # Qout = Qin
  # (r + qbf %*% W)*Ag = qbf*Ag

  # t(W)*qbf - qbf = -r
  # (I-t(W))*qbf=r
  # this gives a system of n equations for n unknowns qbf1...qbf which can be
  # easily solved used base::solve. It also gives a quick way to check that the
  # supplied downslope weighting
  # ignore base flow in channel HSUs as these are distributed by the channel
  # delay histogram
  W[ichan,] <- 0 # recharge to river is by redistributed baseflow  only?

  n <- nrow(W)
  ItW <- identity.matrix(n)-t(W)

  ra<-groups$area * r    # rep(r, n)*g
  ra[ichan]  <- 0
  # default
  Qbf0 <- r * groups$area
#  ra[ichan] <- 0  # # recharge to river is by redistributed baseflow  only?
  Qbf0 <- tryCatch(Qbf0 <- solve(ItW, ra), silent =T)
  qbf0 <- Qbf0/groups$area
  # if matrix singular somthing is wrong with the way the matrix has been set-up
  # need to check in CheckInput

  return(qbf0)
}


# inialise base flows assumming an initial steady state
init.fluxes <- function(groups,      # hsu definitions
                        W,          # flux distribution matrix (not required anymore)
                        ichan,      # channel identifiers
                        qt0,        # initial (specific) catchment discharge / recharge (m/hr)
                        dtt=1)      # inner time step

{

  # ex = saturation excess overland flow
  # uz = gravity drainage from unsaturated into saturated zone
  # qof = saturation excess overland flow
  flows <- data.frame("id"=groups$id, "tag"=groups$tag,
                      "qin"=0,    # qin = total input flow into areal groupings (m^3/hr)
                      "qbf"=0,    # qbf = base flow out of saturated zone (m^3/hr/m^2)
                      "exus"=0,   # exus = root zone excess flow into saturated zone

                      "ex"=0,     # base flow excess
                      "pex"=0,    # precipitation excess
                      "qof"=0,    #
                      "ae"=0,     # actual evap
                      "sat.evap"=0,
                      "rain"=0,
                      "uz"=0)     # uz drainage into water table

  # theoretical max for each group
  qbfmax <- exp(groups$ln_t0-groups$atb.bar)

  # solve for base flows assumming initial steady state  recharge / discharge
  flows$qbf <- init.base.flows(groups,  W, r=as.numeric(qt0), ichan)

  # total upslope area for all elements in each group
  # assumme uniform recharge over all upslope area ai for element i
  # then base flow qi for elemnet is rai and total for area is r.sigma(ai)
  # at steady state specific discharge qbf0 is equivalent to rainfall
  # a.bar is specific discharge for entire group per unit recharge
 # flows$qbf <- qt0*groups$sigma.a/15  # aaargh initialisation depends on times steps as these are aggregated over inner loop when distributed to channel
  # initilise to max value and let the loop sort it out

#  flows[ichan,]$qbf <- 0#qt0/15

  # mass balance implies that input is output minus recharge
  # recharge = r + flux distributed from other areas using weighting matrix
 # r <- DistributeDownslopeFluxes(flows$qbf, w)$qin+qt0
#   qin <- dist.flux(groups, flows, W)$qin - qt0*groups$area
#
#   # something awry with init base flows and / or calculated upslope areas
#   if(any(qin<0))
#   {
#     warning("Given initial specific recharge / discharge qt0 and accumulated uplsope areas sigma.a give -ve input fluxes - check")
#   }
  flows$qin <- 0 #pmax(qin, 0)

  # in steady state drainage to unsat zone is equal to recharge (this *is* correct)
  flows$uz <- as.numeric(qt0)

  return(flows)
}

# transmissivity as a function of depth to water table and limiting (saturated) transmissivity
td <- function(z, m, T0)
{
  return(T0*exp(-z/m))
}

# quick estimate for initial deficits given base flows and groups' hydrological properties
sd0 <- function(groups, flows)
{
  stores$sd <- groups$sd_max

  # q0  saturated specific areal base flows  estimated from areal average
  # of topographic index and saturated transmissivity
  q0 <- exp(groups$ln_t0-groups$atb.bar)

  # initial base flows estimated from initial recharge assumming steady state
  # see Beven 2012 equation 6.1.21
  sdbar<- -groups$m * log(flows$qbf/q0)

  return(max(sdbar, 0))
}

# #############################################################################
# see Beven (2012) eqns B6.1.18 and B6.1.21, pp213-215, for use of areal mean of
# soil-topographic index (Beven, 1986b) to initialise (mean) storage deficts
# -----------------------------------------------------------------------------
init.stores <- function(groups,
                         flows,
                         ichan, qt0)
{
  # storages and deficits for each areal group ([L] per plan area)
  # suz = unsaturated zone storage
  # srz = root zone storage
  # sd = storage deficit (depth to water table)
  # storages all initially zero. Root zone fills from srz0 until reaching
  # "field capacity" (srz_max) - the max RZ storage for the areal group
  stores <- data.frame("id"=groups$id, "tag"=groups$tag, "suz"=0, "srz"=0,
                       "ex"=0,    # excess storage in e.g. overland flow
                       "bf.ex"=0,
                       "sd"=0, "wetting"=TRUE)
  # initialise storage deficit to maximum allowed - allows init base flow to remove some before recharge reduces deficit again
  stores$sd <- groups$sd_max

  # q0  saturated specific areal base flows  estimated from areal average
  # of topographic index and saturated transmissivity
  q0 <- exp(groups$ln_t0-groups$atb.bar)

  # initial base flows estimated from initial recharge assumming steady state
  # see Beven 2012 equation 6.1.21
  sdbar<- pmax(-groups$m * log(flows$qbf/q0), 0, na.rm=T)

  # Assign mean value to areal effective storage deficits (no -ve storages)
  if(any(sdbar<=0))
  {
    LogEvent("Calculated initial storage deficit -ve, setting to 0 (saturation)")
  }
  stores$sd <- sdbar
  if(any(stores$sd>groups$sd_max))
  {
    LogEvent("Calculated initial SD below maximum allowed for groups, set to max")
    stores$sd <- pmin(sdbar, groups$sd_max)
  }

  if(any(groups$srz0 > 1)){warning("Initial root zone storage proportion > 1: ignored")}

  # srz0 is proportion 0 - 1 of root zoone initially filled
  stores$srz <- pmin(abs(groups$srz0*groups$srz_max), groups$srz_max)


  # initialise unsat storage now estimates for initial deficit obtained
  # see formulations from Beven and Wood (1983) relating drainage flux (which we take as equal to
  # steady-state recharge) to unsat storage, defict and time delay parameter
  # qv = suz/td.sd
  # what is the maximum unsat storage given a particular sd?
  stores$suz <- qt0*groups$td*stores$sd

  return(stores)
}



InitialiseChannelFlows <- function(flowlengths, dt, deltaR, q0, chanV)
{
  # channel flow along river reaches - initially zero - use length factors supplied with
  # channel delays table
  qriver <- flowlengths[1,]

  # qchan is flux per downstream length for the river reaches
  # Q0 is the total flow out of the catchment in 1hr.
  # the river volume that flows out of the catchment in 1 hr is qchan[1]*chanV*dt
  # As steady-state solve and apply to entire river network - means that narrower reaches
  # have higher fluxes?
  qr <- q0/chanV
  qriver[] <- qr

  # expand into a vector of flows per 1m lengths of the channel
  qriver <- as.vector(matrix(rep(qriver, deltaR),nrow=deltaR,byrow=TRUE))
  #return(qriver)
}

# use the specified soil /hydraulic properties if supplied or assign from defaults
init.hru<- function(groups,
                            params=def.hsu.par(),
                            chan.par=def.chan.par(),
							ichan=1, rain=NULL)
{
  # groups is table of area groupings and params is list of default parameter values
  # Search for columns for saturation transmissivity ln_t0, max root zone storage srMax
  # and exponential transmissivity profile factor m

  # add missing params using defaults. groups' parameters override defaults
	nms <- setdiff(colnames(groups), names(params))
  #groups[igroup,] <- as.vector(merge.lists(params, def.hsu.par()))
#    }
	if(length(ichan)>0)
	{
	  # the "river" has zero root zone storage and there is no time delay in drainage reaching baseflow
	  # Storage is limited only by bankfull level - set SD and SDMax to large values to simulate
	  # Could apply physically realistic max storage (i.e average river depth at bankfull level) to simulate overbank events
	  groups[ichan,]$srz0 <- chan.par$srz0
	  groups[ichan,]$srz_max <- chan.par$srz_max   # rainfall goes straight to the unsat zone infiltration
	  groups[ichan,]$sd_max <- chan.par$sd_max   # maximum allowable SD
	  groups[ichan,]$vof <-  groups[ichan,]$vchan   # routing velocity

	#  groups[1:nchan,]$SD <- 3  # this is a huge value reflecting the channel's storage capacity
	  groups[ichan,]$td <- chan.par$td  # all drainage will be routed directly to base flow
	}
	# checj taht groups wetness is in decreasing order
	#if(!all(order(groups[-ichan,$atb.bar, decreasing=T)==groups[-ichan,]$order)){warning("Wetness index not in decreasing order of HRU number")}
  # labels for groups default to ids (shoudl these be sequential?)
  if(is.null(groups$tag)){groups$tag<-groups$id}
	# max base flow when deficit =0
	try(groups$qbmax <- exp(groups$ln_t0-groups$atb.bar), silent=T)
	n.gauge <- ifelse(is.null(rain),
										1,
										ncol(rain))
  # ensure gauge reading supplied is withing range
  groups$gauge.id <- pmin(groups$gauge.id, n.gauge)

  # any drainage params
  groups <- init.drainage.params(groups)
  return(groups)
}

# add column for parameter unless it already exists
AddParam <- function(groups, name, val)
{
  if(!name %in% names(groups))
  {
    # add column with default values if not found
    groups <- cbind(groups, list(val))
    names(groups)[length(groups)]<-name
    return(groups)
  }

  # just return the input
  return(groups)

}

# returns the input ts subsetted by a time range
# either bound or neither can be supplied, if a bound is not supplied then the corresponding bound of the
# input series is used
# dt = time step in hours. If the data intervals are greater than these then
# disaggregate or aggregate if they are larger
GetTimeSeriesInputRange <- function(series, start, endt,
                                    cap="rainfall", cols=1:ncol(series),
                                    dt=1, verbose=TRUE)
{
  if(is.null(series)){return(series)}
  #if(length(series)==0){browser()}
	# subset of columns if specified
	series <- series[,cols]
	res <- series[index(series)>=start & index(series)<=endt]
	nas <- which(is.na(res))
  if(length(nas)>0)
  {
    ina <- as.numeric(which(apply(res, MARGIN=1, function(x){any(is.na(x))})))
  	if(length(ina)>0)
  	{
  		if(verbose==TRUE)
  		{
  			cat(paste(length(nas), " ", cap, " NA data encountered: replaced with 0\n"))
  		}
  		res[ina,] <- 0
  	}
  }

	if(verbose==TRUE)
	{
	  # Determine if this value is in mm or m by looking at magnitude
	  # 1m of rainfall or pe is impossible so divide through by 1000!
	  if(max(res,na.rm=TRUE)>10)
	  {
	    message("Large values encountered in series: are data in mm/hr?")
	    #res <- res /1000
	  }
    msg <- paste("Using ", length(res), " ", cap, " observations from ", start, " to ", endt, sep="")
		cat(msg, "\n")
	}
	return(res)
}


require(xts)

# if the max value of the given series or value is greater than 1m it implies that the values are in mm/hr.
# ensure consistency by dividing by the given factor.
convert.values <- function(x, msg="Value(s) converted to m/hr", max=1, fact=1000)
{
  if(max(x, na.rm = T)> max)
  {
    message(msg)
    x <- x/fact
  }
  return(x)
}

init.input <- function(groups, dt, ntt,
                      weights, rain, pe, routing=NULL,
					            ichan=1, i.out=1,
                      qobs=NULL,   # specific or absolute discharge
                      qt0=1e-4,    # initial discharge in m/hr
                      dqds=NULL,
                      disp.par, run.par,
                      sim.start, sim.end,
                      calling.env = NULL)
{
  if(is.null(i.out)){i.out<-1}
  # combine lists, supplied attributes overriding the defaults and missing values being supplied
  # by defaults
  disp.par <- merge.lists(disp.par(), disp.par)
  run.par <- merge.lists(def.run.par(), run.par)
    # need this for $ notation
  groups<-data.frame(groups)
  runid <- run.par$id
  qmax <- disp.par$max.q

  if(length(rain)==0)
  {
    stop("Rainfall input required")
  }

  tz <- ""
  # use the rainfall data to infer the time step and run start and end times, if not otherwise specified
  if(is.null(dt))
  {
    dt <- as.numeric(difftime(index(rain)[2], index(rain)[1], units="hour"))
  }

  # default limits of run period inferred from time series input
  sim.start <- max(start(rain), as.POSIXct(sim.start), na.rm=T)
  sim.end <- min(end(rain), as.POSIXct(sim.end), na.rm=T)

  dispStart <- sim.start + disp.par$"graphics.delay"*3600 # hours -> seconds

  # check that the river reach distribution is of the correct dimensions -
  # should be ngroup x nreach, where nth row determines the distribution of the
  # nth group into the river network
  ngroups <- nrow(groups)
  # number of river reaches determined from flow length distribution - each row
  # represents a reach length and proportion of channel within that reach
  # nreach <- nrow(routing)

  # inner time step
  dtt <- dt/ntt

  nchan <- length(ichan)
  # channel identifiers
#  ichan <- 1:nchan
  weights <- as.matrix(weights)

  # logging output (defaults to console if not specified)
  if(nchar(run.par$log.out)>0)
  {
 #   glogout <- file(run.par$log.out)
    # route messages to a file, otherwise sent straight to console
   # sink(run.par$log.out, type="message", append=T)
  }

  # check consistent matrix dims, array sizes and values
  # check times series input in correct format and duration

  LogEvent("Dynamic TOPMODEL run started")

  # times for output series, dt in hours
  times <- seq(sim.start, sim.end, by = 3600*dt)
  nmax <- length(times)

  # create dummy pe and qobs if none supplied, Q) can be used to set absolute discharge max limit
  if(is.null(pe)){
    pe <- xts(rep(0,nmax), order.by=times)
  }
  # locate and fix any NA records
  pe[which(is.na(pe))] <-0

  # ensure everything in same units
  rain <- convert.values(rain)
  qt0 <- convert.values(qt0)

  pe <- convert.values(pe)

  qobs <- GetTimeSeriesInputRange(qobs, sim.start, sim.end, verbose=FALSE)

  # dummy series if nothing supplied or outside range of simulation
  if(length(qobs)>0)
  {
		qobs <- xts(qobs)

		# want a specific input, observed values often (always?) in cubic m^3/s
		qobs <- convert.values(qobs, msg="Observed flows to m/hr: check that specfic discharge in m/hr has been supplied")
  }

  # use the overland velocity as proxy for channel routing velocity
 # def.chan.par$vchan <- vchan
  # initialise storages and deficits for each group add static parameters to
  # areal groups if a soils map provided then fix up the soil parameters from
  # here, otherwise continue to use default
  groups <- init.hru(groups, ichan=ichan, rain=rain)
  # Input rainfall, observed flows and p.e.
  rain <- GetTimeSeriesInputRange(rain, sim.start, sim.end, verbose=FALSE)
  pe <- GetTimeSeriesInputRange(pe, sim.start, sim.end, verbose=FALSE)
  # add the actual evap
  evap <- cbind("pe"=pe, "ae"=NA)

  try(check.and.fix.input(environment()))

	# calculated time delay matrix for channel routing. removes initial column with flow lengths
	routing <- get.routing(groups,
											 weights,
											 routing,
											 ichan,	dt)
  # max number of time steps that river flux can be shifted
  # Time series of total output flows from output reach(s). add enough extra elements to hold time delayed flows
  n.shift <- nrow(routing)
  times.ex <-seq(sim.start, sim.end+nrow(routing)*dt*3600, by = 3600*dt)

  nms <- groups[ichan,]$id

  # numbet of outlet reaches to produce flow output
  nout <- ncol(routing)
  Qr <- xts(matrix(0,nrow=nmax+nrow(routing), ncol=nout), order.by=times.ex)

  #colnames(Qr)<- nms
# total overland flow contribution
  Qof <- xts(matrix(0,nrow=nmax,ncol=1), order.by=times)


  #*****************************************************************************
  # initialisation of flows, storages and deficits for each areal group flows
  # (at time step t, per plan area except where stated)
#  message("... initialising fluxes and storages...\n")

  # can determine a theorectical max for discharge: set all base flow to max ,
  # controlled by saturation transmissivity (z=0)  add in any rainfall that
  # could be transferred from overland flow, distribute to channel by flux routing matrix,
  # then route to outlet by time delay table
  # see discussion of TOPMODEL theory in Beven (2012), Chapter 6, Box 6.1, p210
  # gamma is the average soil-topographic index defined in  Beven (1986a)
  # typical values for ln_t0 are in range (-7,8).
  # Could examine catchment and its observed discharges to calibrate this
  groups$qbf_max <- exp(groups$ln_t0-groups$atb.bar)
  # huge value ensures no excess in river (or set to zero to ensure that ay base flow in river produces "surface" run-off)
  groups$qbf_max[ichan]<-1e6

 # weights <- spam::as.spam(weights)

  # in the summer, effective rainfall recharge will be affected significantly be ea.
  # as we don't know the storages yet, apply a

  # estimate for initial base fluxes given a steady state
  flows <- init.fluxes(groups=groups,
                      weights,           # flux distribution matrix
                      ichan,       # channel identifiers
                      qt0,         # initial recharge / discharge, assummed equal in steady state, taking into account evapotranspiration
                      dtt)

  # storage deficits estimated from relationship with base flows estimated above
  stores <- init.stores(groups, flows, ichan,
                         qt0)


  # initial determined from storages vs. deficit

  # index in time series after which to start showing results
  i.start <- which(times>=sim.start)[1]
  if(is.na(i.start))
  {stop(paste("Error determining start of simulation period ", dispStart))}
  # 	LogEvent(paste("initialising time series: start ", as.character(start)))


  # structure to save the average storage deficit. ae, river input, water balance etc etc
  dat <- data.frame("ae"=rep(0,nmax))
                   # "wb"=rep(0,nmax),
                    #"qr"=rep(0,nmax),
                    #"qriv.in"=rep(0,nmax),   # riveer input fluxes
                    #"ae"=rep(0,nmax))
  data.out <- xts(dat, order.by=times)

  qriver <- xts(matrix(rep(0,nmax,each=nchan), ncol=nchan), order.by=times)
  hsus <- groups[-ichan,]
  nhsu <- nrow(hsus)

  fluxes <- array(dim=c(nmax, ngroups, 7), dimnames=c( "it", "hsu", "flux"))
  storages <- array(dim=c(nmax, ngroups, 4), dimnames=c( "it", "hsu", "store"))
  # graphics interval expressed in hours: convert to time steps

  # start of outer time stepping logic using a time step at each rainfall observation
  time<- sim.start
  it <- 1

  # initialise excess (overland) flow to 0
  q.ex <- rep(0, nrow(groups))

# initialisation
    # run this number of time steps to get an equilibrium with
  # rainfall and specified initial discharge
  # number of time steps to wait before start
#   n.step <- round(run.par$sim.delay/dt)
#
#   flows <- init$flows
#   stores <- init$stores
  # set  first discharge estimates to value specified (total, not specific)
  # extend into lag period whilst flux is stil in transit. If this isn't done then it appears that
  # there is a period when no flux enters and we see a initial dip in the hydrograph
 # Qr[it:(it+n.shift),]<- qt0*sum(groups$area)
  nroute <- length(routing)
  flows1 <- flows
  # at steady state the flows into the channel equal those at the outlet
  flows1$qin[ichan] <- qt0*sum(groups$area)

  Qr1 <- Qr
  # distribute the initial flow to future flows
  for(it1 in 1:nroute)
  {
    Qr1 <- route.channel.flows(groups, flows1, delays=routing,
                          weights, Qr1, it=it1, dt, ichan)
  }
  Qr[1:nroute,]<- qt0*sum(groups$area)

  # remove distributed flows from initial steps of river outlet series
  Qr[1:nroute,] <- Qr[1:nroute,] - Qr1[1:nroute,]

# set channel deficts to a typical value
	stores[ichan,]$sd <- groups[ichan, ]$sd_max/2
  run.par$sim.delay <- 0

	# initialise surface excess storage to zero
	surf.ex <- stores$ex
	pe.dist <- xts(order.by=index(pe), matrix(rep(pe, ngroups), ncol=ngroups))

  # open graphics windows, calculate interval etc
  disp.par <-setup.output(disp.par, run.par, rain=rain, qobs=qobs)


    # function specifying gradient of q - s function defaults to exponential
    if(is.null(dqds))
    {
        dqds <-
        function(q, S, params)
        {                   #
            return(-q/params$m)
        }
    }

  rain <- as.xts(rain)
  pe <- as.xts(pe)
  # ensure that if it's raining there is nil evapotranspiration
  tm.rain <- index(rain)[which(rain>0)]
  pe[tm.rain]<- 0

  sim.end <- as.POSIXct(min(max(index(pe)), max(index(rain)), sim.end))
  cur.env <-  environment()
  # copy values to the specified environment, typicaly the  current environment in
  # the calling context
  if(!is.null(calling.env))
  {
    for(n in ls(cur.env, all.names=TRUE))
    {
      # assign name-value in this envir to that specified
      assign(n, get(n, cur.env), calling.env)
    }
  }


  return(cur.env)
}


init.drainage.params <- function(groups)
{
  idrained <- which(!is.null(groups$drained) && groups$drained)
  if(length(idrained)>0)
  {
    # set the drained groups' parameters from the field capacity phi, drained thickness dz.drain and top of drained layer z.drain
    # drainage depth and field capacity
    groups[idrained,]$sd_drain_max <- groups[idrained,]$z_drain * groups[idrained,]$phi
    groups[idrained,]$sd_drain_min <- (groups[idrained,]$z_drain-groups[idrained,]$dz_drain)* groups[idrained,]$phi
  }
  return(groups)
}



