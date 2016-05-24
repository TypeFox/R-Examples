#' Discrete a catchment into hydrological response units (HRUs)
#'
#' @description Discrete a catchment into a set hydrological response units (HRUs) according to any number of landscape layers and cuts
#' @details This applies the given cuts to the supplied landscape layers to produce areal groupings of the catchment.
#' @export discretise
#' @import raster
#' @param layers A multi-band raster (stack) comprising the catchment data. This should be in a projected coordinate system (or none) and have reqular cells. The first layer should be the elevation raster, and subsequent (named) layers should supply the landscape data drawn in to create the discretisation
#' @param cuts A list of cuts of the form layer_name=number. Each name should correspond to a layer name in the layers parameter.
#' @param order.by Name of layer whose values will be use to sort the response units, in decreasing order. Defaults to the name of the first cut
#' @param area.thresh	Minimum area for response units, expressed as a percentage of the catchment plan area, excluding channel cells. Areas smaller than this are aggregate dwith adjacent areas until exceeding the threshold aea
#' @param chans	Raster containing channel reach locations, of the same dimensions and resolution of the DEM and other catchment layers. The reaches should be numbered sequentially and any areas not containing part of the channel should be NA. If a second band is supplied with values 0-1 then this is taken to be the proportion of the corresponding non-zero cell occuppied by the channel. If this layer is not present then the proportion is infered from the channel width as p = min(1, chan.width/xres(dem))
#' @param chan.width Channel width, in same units as DEM. Only used if chans doesn't contain a layer to specify the proportion of each river cell comprised of the channel.
#' @return A list comprising the following
#' @return weights	Flux distribution (weighting) matrix. A ngroup x ngroup matrix defining the downslope flux distributions between groups, between land and the channel, and between channel reaches, where nh is the number of land discretisations identified by applying cuts to the catchment layers and nc the numbbr of channel reaches defined. The nth row gives the proportions of flow out of HRU #n to other response units and the channel. Row sums should thus always add to 1. The m th column gives the proportion of flow from the other response units into the m th group.
#' @return groups A data frame whose rows comprising the names, plan area and model parameters of each response unit. See Beven and Freer (2001) and Metcalfe et al (2015) for a description of these parameters
#' @return hru	Multi-band raster comprising the original rasters that the specifed cuts were applied to produce the discretisation; the channel network;the resultant response unit locations
#' @examples
#' # Landcover and soils are fairly homogenous throughout the Brompton catchment;
#' # storm response of the appears to be mostly controlled by proximity to the
#' # channel network. A simple discretisation according to flow distance from the
#' # nearest channel thus appears to capture the dynamics during the  2012 event
#' # without introducing unnecessary complexity.
#'\dontrun{
#' require(dynatopmodel)
#'
#' data(brompton)
#'
#' chans <- build.chans(brompton$dem, drn=brompton$drn, chan.width=2)
#' # sort by distance but want areas closest the channel to come first
#' layers <- addLayer(brompton$dem, 2000-brompton$flowdists)

#' disc <- discretise(layers, cuts=c(flowdists=5), chans=chans, area.thresh=2/100)
#'
#' write.table(disc$groups, sep="\t", row.names=FALSE)
#'}

discretise <- function(layers,
                       chans,
                       cuts=list(a=10),
                       area.thresh=2/100,
                       order.by=names(cuts)[[1]],
                       chan.width=5)
{
  dem <- layers[[1]]
  catch <- layers

  compareRaster(dem, chans)

  if(is.null(area.thresh)){
    area.thresh<- get.defs()$area.thresh
  }

  if(is.null(chan.width)){
    chan.width<- get.defs()$chan.width
  }

  if(area.thresh>=1){area.thresh<-area.thresh/100}

  if(nlayers(chans)>1)
  {
    cellprops <- chans[[2]]
  }
  else
  {
    cellprops <- min(1, chan.width/xres(chans))
  }

  nchan <- length(unique(chans[[1]], na.rm=T))
  #  nchan <- nrow(drn)
  # channel identifiers
  ichan <- 1:nchan

  message("Combining layers...")
  cm <- combine.groupings(dem, catch=catch, chans=chans,
                          cuts=cuts, thresh=area.thresh)
  #nms <- names(cuts)
  # default sort order is using upslope area

  # reorder by upslope area
  zonal.vals <- raster::zonal(catch, cm[[1]])
  # create sequences for ordering the HRUs
  ord.args <- lapply(as.list(order.by),
                     function(x)
                     {
                       zonal.vals[, x]
                     })
  ords <- do.call(order, c(ord.args, decreasing=T))  # by default zone with higher averages appear first
  sub.df <- data.frame(cbind(zonal.vals[,"zone"], zonal.vals[ords,"zone"]))

  cm[[1]] <- subs(cm[[1]], sub.df)
  nms <- names(cm)
  nms[[1]]<- "HRU"
  names(cm)<-nms

  message("Building group info table....")

  # layers that went into discretisation are held in further layers of classification matrix:
  # pass into proc to produce group summary table
  groups <- data.frame(build.hru.table(cm, dem=dem,
                                       reaches=chans[[1]],
                                       cellareas=1-cellprops))

  groups[ichan,"chan.no"] <- ichan
  groups[ichan,"vof"] <- NA
  groups[-ichan,"vchan"] <- NA

  # add in zonal info
  #	zchan <- zonal.vals[1, names(cuts)]
  #  	zchan[] <- NA
  #  	groups <- cbind(groups, rbind(zchan, zonal.vals[ords, names(cuts)]))
  nms <- names(cm)
  # invalid cuts removed
  layer.nms <- nms[2:length(nms)]

  w<-get.flow.distribution.matrix(dem,
                                  cm=cm[[1]],
                                  reaches=addLayer(chans[[1]], cellprops))

  return(list(
    "groups"=groups,
    "hru"=cm[[1]],
    "weights"=w))
}





#load.source("dtm.main.r", chdir=T)

require(raster)
require(xts)

get.sim.range <- function(proj)
{
  #  require(intervals)
  obs <- proj$obs$obs
  pe <- proj$obs$pe
  qobs <- proj$obs$qobs
  try(proj$sim.start <- as.POSIXct(proj$sim.start), silent=T)
  try(proj$sim.end <- as.POSIXct(proj$sim.end), silent=T)

  if(length(proj$sim.start)==0)
  {
    s1 <- NA
    s2 <- NA
    s3 <- NA
    try(s1 <- start(obs), silent=T)
    try(s2 <- start(pe), silent=T)
    try(s3 <- start(qobs), silent=T)

    proj$sim.start <- max(c(s1, s2, s3), na.rm=T)
    if(!is.null(proj$sim.start))
    {
      message(paste("Start of simulation inferred from input as ", proj$sim.start))
    }
  }
  if(length(proj$sim.end)==0)
  {
    e1 <- NA
    e2 <- NA
    e3 <- NA
    try(e1 <- end(obs), silent=T)
    try(e2 <- end(pe), silent=T)
    try(e3 <- end(qobs), silent=T)
    proj$sim.end <- min(c(e1, e2, e3), na.rm=T)
    if(!is.null(proj$sim.end))
    {
      message(paste("End of simulation inferred from input as ", proj$sim.end))
    }
    if(proj$sim.end < proj$sim.start)
    {
      stop("Error: sim.start after end. Check supplied data and values")
    }
  }
  return(proj)
}


exists.not.null <- function(obj.name, check.file=T, warn=NULL)
{
  res <- FALSE
  # calling frame / environment, up one level in calling stack, must be at least one
  p.env <- sys.frame(-1)
  # treat the argument as a charcater object name and look in the calling frame
  if(exists(obj.name, where=p.env))
  {
    obj <- get(obj.name, p.env)
    if(!is.null(obj))
    {
      res <- ifelse(check.file,
                    file.exists(obj),
                    TRUE)

    }
  }
  if(!res & !is.null(warn))
  {
    warning(warn)
  }
  return(res)
}


# use catchment area calculated from dem to determine specific discharges from the
#input given in cu.m/sec
convert.to.specific.discharges <- function(proj, q)
{
	if(max(q))
	# catchment area
  a <- with(proj, sum(length(which(!is.na(dem[])))*xres(dem)*yres(dem)))
  # assumme in cu.m/s
  res <- 3600*q/a
  if(max(res, na.rm=T)>1){warning("Very large specific discharges calculated: check input not in mm")}

  return(res)
}

add.layers <- function(proj, reload=F)
{
	nms <-  dir(proj$dir, "\\.shp")
	if(!reload)
	{
		nms <- setdiff(nms, paste0(names(proj), ".shp"))
	}
	#  try(proj$dem <- rast.mem.copy(raster(file.path(data.dir, "dem.tif"))))
	#try(proj$drn <- readOGR(data.dir, "drn"))
	# iterate through any shape files located and add verbatim
	for(shp in nms)
	{
		# name without extention
		shp.nm <- sub("*.shp", "", shp)
		cat("Loading shape file ", shp.nm, "...")
		try(proj[[shp.nm]] <- readOGR(proj$dir, shp.nm))
		if(!is.null(proj[[shp.nm]])){cat("...done\n")}
		else{cat("...failed\n")}
	}
	return(proj)
}

check.time.intervals <- function(proj)
{
  int <- time.interval.intersection(proj$obs$rain, proj$sim.start, proj$sim.end)
  return(length(int)>0)
}

time.interval.intersection <- function(obs, sim.start, sim.end)
{
  int <- which(as.numeric(index(obs))<= as.numeric(sim.end) &
                 as.numeric(index(obs))>= as.numeric(sim.start))
  if(length(int)>0)
  {
    return(index(obs)[int])
  }
  return(NULL)
}

# given a set of observations and a specified run interval
# expand / contrcat the run times to accommodate the observations
fix.run.dates <- function(proj)
{
  obs <- proj$obs$rain

  if(!is.null(obs) & !check.time.intervals(proj))
  {

    warning("No rainfall data within specified run start/ end times: adjusting...")
    len.sim <- proj$sim.start-proj$sim.end

    start <- start(index(obs$rain))
    cat("Setting sim start to ", "\n")
    proj$sim.start <- start
    end <- min(start + len.sim, end(index(obs$rain)))
    cat("Setting sim end to", end, "\n")
    proj$sim.end <- as.POSIXct(end)
  }
  return(proj)
}

aggregate_observations <- function(proj)
{
  try(proj<-fix.run.dates(proj))
  obs <- proj$obs
  # check that the specified start end and end dates contain at least some rainfall data. Other data are
  # less important and take null defaults if not specified
  try(obs$pe <- disaggregate_xts(proj$obs$pe,
                                 ser.start=proj$sim.start,
                                 ser.end=proj$sim.end,
                                 dt=proj$dt, is.rate=T))
  try(obs$rain <- disaggregate_xts(proj$obs$rain,
                                   ser.start=proj$sim.start,
                                   ser.end=proj$sim.end,
                                   dt=proj$dt, is.rate=T))
  # note observed flows required in specific discharge m/hr
  try(obs$qobs <- disaggregate_xts(proj$obs$qobs,
                                   ser.start=proj$sim.start,
                                   ser.end=proj$sim.end,
                                   dt=proj$dt, is.rate=T))

  return(obs)
}

# simple output of results, allowing the simulated values to be supplied separately
plot.run <- function(run,
                     qsim=NULL,
                     cols=c("blue", "green"),
                     fn=NULL,
                     main="",
                     start = run$sim.start,
                     evap=run$evap[,"ae"]*1000,
                     end= run$sim.end,
                     par=disp.par(),
                     ...)
{
  sim.start <- as.POSIXct(start, tz="GMT")
  sim.end <- as.POSIXct(end, tz="GMT")

  sel <- paste0(sim.start, "::", sim.end)
  rain <- as.xts(run$rain)[sel]*1000
  qobs <- run$qobs[sel]*1000
  if(is.null(qsim))
  {
    qsim <- run$qsim*1000  #r[sel]/run$catch.area*1000
  }

  if(!is.null(evap))
  {
    ae <- evap[sel]
  }
  else
  {
    ae <- NULL
  }
  qsim <- qsim[sel]

  #  try(print(NSE(qobs[sel]), silent=T))

  # try(print(format("time at peak =%H:%M", format(time_at_peak(qsim))))
  #	nresp<- ncol(qresp)

  par$max.q <- max(c(disp.par$qmax, qobs[], qsim[]), na.rm=T)

  # par(mar=c(4,4,3,4))
  if(length(fn)>0)
  {
    jpeg(filename=fn, width=1024, height=768)
    # larger axes labels
    par("cex.lab"=1.5)
    on.exit(dev.off())
    par(mar=c(3,4,3,4.5))
    title <- ""
  }

  # everything
  disp.output(main=main,
              qsim=qsim,
                  #		start=sim.start,
                  #		end=sim.end,
              evap=ae,
              rain=rain,
              tm=NULL,
              qobs=qobs,
              par=par, ...)

}


# # plot of simulated discharges etc after
# disp.run.results <- function(run,
#                              qmax=NULL, legend=F,
#                              title = "",
#                              disp.par = disp.par(),
#                              ...)
# {
#   qobs <- run$qobs
#   layout(matrix(1))
#   qmax <- max(run$qsim, run$qobs, na.rm=T)*1000*1.25
#   qsim <- run
#   pe <- run$evap
#   par(family="serif")
#   disp.par$legend.show <- legend
#   disp.par$title.main<- title
#   dDischargeSelection(qsim=run$qsim, evap=run$evap[,"ae"], rain=run$rain, qobs= run$qobs,
#                             qmax=qmax,disp.par=disp.par,...)
#   #, run.par=run.par)
#
#
# }





# apply the given parameters to groups in all discretisations
apply.params <- function(proj, params, which=1:length(proj$disc))
{
  if(length(proj$disc)==0){return(proj)}
  params <- params[which(names(params) %in% colnames(proj$disc[[1]]$groups))]
  if(length(params)==0){return(proj)}

  proj$disc[which] <- lapply(proj$disc[which],
                             function(disc)
                             {

                               vals <- matrix(rep(unlist(params),nrow(disc$groups)), nrow=nrow(disc$groups), byrow=T)
                               disc$groups[,names(params)]<- vals

                               return(disc)
                             }
  )
  return(proj)
}


# graphics and text output
gr.on <- function(proj, spatial=F)
{
  return(graphics.on.off(proj,T, spatial))
}

gr.off <- function(proj, spatial=F)
{
  return(graphics.on.off(proj, F, spatial))
}

graphics.on.off <- function(proj, val=T, spatial=F)
{
  proj$disp.par$graphics.show <- val
  if(spatial)
  {
    proj$disp.par$graphics.spatial.show <- val

  }
  return(proj)
}

output.off <- function(proj)
{
  proj$disp.par$text.out <- NULL
  return(proj)
}

output.on <- function(proj)
{
  proj$disp.par$text <- stdout()

}




# plot a set of results associated with a topmodel run
plot.run.response <- function(run,   # dtm run output
                               qresp, # response series
                               evap=NULL,
                               ymax=NULL,
                                show.qobs=F,
                               fn=NULL,
                               lwd=1,
                               lty=1,
                               cols=c("blue", rainbow(ncol(qresp)-1)),
                               title="",
                              start=index(qresp)[1],
                              end=index(qresp)[length(index(qresp))],

                              ...)
{
  sel <- paste0(start, "::", end)
  run$qsim <- run$qsim[sel]
  qresp <- qresp[sel]
  if(length(evap)>0)
  {
    run$evap <- evap
    names(run$evap)<- "ae"
  }
  else
  {
    run$evap <- NULL
  }
  nresp<- ncol(qresp)
  if(!show.qobs)
  {
    run$qobs<-NULL
  }

  plot.run(run, qsim=qresp,
           cols=cols, fn=fn,
           title=title,
           ymax=ymax,
           lwd=lwd,
           lty=lty,
           start=start,
           end=end,
           #legend=F,
         #  legend.col=cols,
           ...)


}


# plot the results of a run using another matrix of q and a source project
plot.q <- function(proj,
                   qresp,
                   evap=NULL,
                   ymax=NULL,
                   fn=NULL,
                   lwd=2,
                   lty=1,
                   cols=rainbow(ncol(qresp)),
                   title="", ...)
{
  run <- proj
  run$disc <- NULL
  run$proj <- proj
  run$qobs <- NULL
  run$rain <- proj$obs$rain

  if(evap==F)
  {
    evap <- NULL
  }
  else if(is.null(evap) & !is.null(proj$obs$pe))
  {
    run$evap <- proj$obs$pe
    names(run$evap)<- "ae"
  }
  nresp<- ncol(qresp)

  plot.run(run, qsim=qresp*1000,
           cols=cols, fn=fn,
           title=title,
           ymax=ymax,
           lwd=lwd,
           lty=lty,
           ...)

}

get.calib.dir <- function(proj.ex)
{

  s <- format(proj.ex$sim.start, "s=%Y-%m-%d")
  e <- format(proj.ex$sim.end, "e=%Y-%m-%d")
  return(paste0(s, e, collapse=","))
}

# return the goodness of fit for a simulation run
gof.run <- function(run, pos=1) # specify column for multiple results
                    #s=),
                    #e=)
{
  sel <- paste0(first(index(run$qsim)), "::", last(index(run$qsim)))
  pos <- min(pos, ncol(run$obs$qobs))
  qobs <- run$qobs[sel, pos]
  if(!is.null(run$qobs) & nrow(qobs)==nrow(run$qsim))
  {

    res <- run.gof(run$qsim, qobs)
    #     res <- run.gof(as.vector(subset_zoo(run$qsim, s, e)),
    #                    as.vector(subset_zoo(run$qobs[,pos], s, e)))
    #
    return(res)
  }
}

