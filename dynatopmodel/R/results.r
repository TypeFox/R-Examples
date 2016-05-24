require(lubridate)
update.output <- function(output, groups, stores, it, ichan,
                          items=NULL)
{
  # weighted average of storage deficits
 # output[time, "sd"]<-
  # record river input flow
 # output[time, "qriv.in"] <- sum(Qchan)

  # water balance: rain in, evap and specific discharge out
#  output[it, "wb"]  <- output[it, "wb"] + as.numeric(current.input(groups, rain[it,], output[it, "ae"], Qr[it,]/catchArea))
#  output[time, "ae"] <- weighted.average(ae, groups$area)
#  storage.deficits[time,]<- stores[-ichan, "sd"]
#  base.flows[time,]<- flows[-ichan, "qbf"]
  for(nm in names(items))
  {
    if(!nm %in% names(output))
    {
      output <- cbind(output, rep(0, nrow(output)))
      names(output)<-c(names(output)[1:(length(names(output))-1)], nm)
    }
    output[it, nm] <- items[[nm]]

  }

    return(output)
}


#
RunSummary <- function(groups, stores, storage.in, ichan, qsim, qobs, start.time,
                       rain, ae=0, text.out)
                     #  sat.ex =0,
                    # qbf.ex = 0, rain=0,
                    #   ae=0)
{
  if(!is.null(text.out)){return()}
  run.time <- difftime(Sys.time(), start.time, units="auto")
  # depending on length of run could be mins, secs (or hours!)
  units <- units(run.time)
  # difference between net input and net storage gain
  storage.gain <- current.storage(groups, stores, ichan)-storage.in
  wb <- water.balance(groups, stores, storage.in, ichan, qsim,
                       rain, ae)

  #sum(rain-ae)-sum(qsim)-storage.gain

  LogEvent(paste("Run completed in ", round(run.time), units))
  # water balance etc
  cat("Total river discharge      = ", round(sum(qsim),2), "m\n")
#  cat("   saturated excess        = ", round(sum(sat.ex),3), "m)\n")
#  cat("   base flow excess        = ", round(sum(qbf.ex),3), "m)\n")
  cat("Actual  evapotranspiration = ", round(sum(ae),2), "m\n")
  cat("Total rain input           = ", round(sum(rain),2), "m\n")
  cat("Net storage gain           = ", round(sum(storage.gain),2), "m\n")

  cat("Water balance              = ", round(wb,2), "m\n")

}



# # initialise the directory for run info - call before run to initialise
#
#
#
# SaveGraphicsOutput <- function(save.dir, tm,
#                                title, runid, qr,
#                                evap, rain,
#                                groups, qobs,
#                                wb,
#                                eff,
#                                ymax,
#                                disp.par,
#                                run.par,
#                                buff,ichan)
# {
#   # save a snapshot of observed / simulated discharges
# #   fn <- paste0(format(tm, disp.par$graphics.fn.fmt), ".jpg")
# #
# #   jpeg(1024, 768, file=file.path(save.dir, fn))
#  #  browser()
# #   disp.output(title, runid, qr,
# #                     evap, rain, tm,
# #                     groups, qobs,
# #                     wb,
# #                     eff,
# #                     ymax,
# #                     disp.par=disp.par,
# #                     run.par=run.par,
# #                     buff,
# #                   ichan=ichan)
#
#   # fmt by default include extension
#   fn <- format(tm, disp.par$graphics.fn.fmt)
#   width <- disp.par$graphics.save.width
#   height <- disp.par$graphics.save.height
#   devn <- dev.copy(device=jpeg, width=width, height=height,
#                       filename=file.path(save.dir, fn))
#   dev.off()
#
#   #
#   #       disp.output(title, runid, qsim, pe, rain, tm,
#   #                        groups, qobs, bal, eff=eff,
#   #                        timeBuffer)
#   #       dev.off()
# }




# identify the first window of default type on list with hydrograph output
# call once to estalish titled output windows
setup.output <- function(disp.par,
                        run.par,
                        rain,
                        qobs=NULL,
												sim.start=NULL, dt=1, pe=NULL)

{
	if(disp.par$"graphics.show"==T)
	{
    disp.par <- merge.lists(disp.par(), disp.par)
    # simStart <-  start + run.par$"sim.delay"*3600  # add simulation delay (expressed in hours)
    disp.par$disp.start <- sim.start + disp.par$"graphics.delay"*3600 # hours -> seconds

    if(!(disp.par$text.out == "N"))
    {
        disp.par$text.out <- stdout()
    }
    else
    {
        disp.par$text.out <- NULL
    }
    if(nchar(disp.par$text.out)>0)
    {
        #  file specified
        #   text.out <- file(disp.par$text.out)
    }

    if(disp.par$"graphics.save"==T & !file.exists(disp.par$"graphics.out"))
    {
      # mmmm
        dir.create(disp.par$"graphics.out", recursive=T)
    }
    # graphics interval is in hours but needs to be converted if time step is not
    disp.par$graphics.interval <- max(round(disp.par$graphics.interval/dt), 1)

    disp.par$time.tick.int <- "day"

    # set the tick interval according to the size of the window
    if(disp.par$graphics.window.length > 3*7*24)
    {
      disp.par$time.tick.int <- "week"
    }

    # ensure suffcient windows open to display results
    disp.par$winid <- open.devices(1, title=run.par$id,
                width=12, #disp.par$graphics.save.width,
                height=8) #disp.par$graphics.save.height)

    # for display should be in mm/hr
 #   disp.par$max.q <- convert.vals(disp.par$max.q)
#    disp.par$max.rain <- convert.vals(disp.par$max.rain)
	}
  # return the updated display parameters
  return(disp.par)
}



SetFigWindow <- function(xmn=0, xmx=1, ymn=0, ymx=1)
{
  # set the fg parameter to place a window at the given position expressed in
  # proportion of the screen. Call with no parameters to reset to default window
  par(fig=c(xmn,xmx,ymn,ymx))
}

# Plot discharge predictions, rainfall and potential (actual?) evapotranspiration,
# storage deficits and root and unsat zone storages.
# time buffer is in hours
disp.results <- function (it,  # current time step
                          tm,  # current simulation time
                          qr,  # calculated discharge at outlet
                          rain,
                          evap,  # actula evapotranspiration
                          groups,
                          flows,
                          stores,
                          qobs=NULL,
                          wb=NULL,
                          ichan=1,   # channel indicators
                          text.out=stdout(),
                          log.msg = "",
                          start, end,
                          disp.par,
                          run.par, sf=2)
{
    fmt <-  disp.par$"time.text.fmt"
    # sim.delay expressed as hours, convert to seconds
    if(tm >= start & !is.null(text.out))
    {
        # send to console by default,  disp.par$output.text.flows
        #  txt.out <- paste(flows[,disp.par$output.text.flows], sep="\t", )
        cat(format(tm, fmt), "\t", signif(qr[it,ichan],sf), file=text.out)   #
#         if(length(qobs)>0)
#         {
#           cat("\t\t", signif(qobs[it,], sf), file=text.out)
#           cat("\t\t", signif((qr[it,ichan]-qobs[it,])/qobs[it,]*100, 2), "%")
#         }
		if(any(flows$qof>0))
		{
			# overland flow distributed to channel over this time step
            cat("\t\t", signif(flows[ichan, "qof"]*1000*groups[ichan,]$area/sum(groups$area),2), file=text.out)   #
		}
		#cat("\t\tmm/hr\n")   #
		cat("\n")   #
    }
    else{
        # waiting...
    #    cat(".")
    }

    disp.par$disp.start <- start + run.par$"sim.delay"*3600

    if(!disp.par$"graphics.show" | tm < disp.par$disp.start){return()}

  title <- disp.par$"main"

  # render either graphical output at this time step?
  show.graphics <-  it%%disp.par$"graphics.interval"==0

  # time buffer, the length of time around the current observation that should be displayed
  buff <- disp.par$"graphics.window.length"

  # graphical output, if requested, after specified time interval
  if(show.graphics)
  {
   # dev.set(disp.par$winid)
 		activate.device(1)
    # buffer output
    dev.hold()
    on.exit(dev.flush())

    qr <- GetTimeSeriesInputRange(qr, start, tm, verbose=FALSE)

    # determine y axis range from
    # (a) actual evapotranspiration
    # (b) any observed discharges
    # (c) max base flow routed through catchment outlet
    # (d) explicitly specified limits
#    ymax <- max(max(evap, disp.par$max.q, na.rm=T), qobs[], na.rm=T)
#    ymax <- max(c(disp.par$max.q, qobs[], 0.01), na.rm=T)
    par("mar"=c(5, 4.5, 5, 4.5))
    evap <- evap[,"ae"]

    # flows observed and simulated
    disp_output(main=title,
                qsim=qr,
                evap=evap,
                rain=rain,
                tm=tm,
                qobs=qobs,
                par=disp.par)
  }
}



disp.discharge.selection <- function(groups,
                                      qsim,  # in m/hr
                                      qobs,  #  "
                                      rain, evap,
                                      sel=NULL,
                                      max.q=0,   # upper bound for discharge display, in mm/hr
                                      ichan=1,
                                      disp.par=disp.par,
                                      title=disp.par$title.main,...
                                      )

{
  # plot full range
#  if(e<=s){stop("invalid time selection")}
  if(is.null(sel))
  {
    sel<-range(index(qsim))
  }

  buf <- as.numeric(difftime(sel[2], sel[1], units="hours")) #[2]-sel[1]
  #  difftime(end, start, end,units="hours")
  # calculate R for this period
#  sel <- sel<-paste(s, "::", e, sep="")
  #qobs <- qobs[sel]
  ae <- xts(1000*subset_zoo(evap, sel[1], sel[2]))
  if(!is.null(qobs))
  {
    qobs <- xts(1000*subset_zoo(qobs, sel[1], sel[2]))}
  qsim <- xts(1000*subset_zoo(qsim, sel[1], sel[2]))

  mar <- c(4, 4, 5, 5)
  if(!disp.par$legend.show)
  {
    # less space at base
  #  mar[1]<-2
  }
  if(title=="")
  {
    # don't need so much space at top
    mar[3]<-3
  }
    par("mar"=mar)

  if(!is.null(qobs))
  {
    #names(qobs)<-"Observed discharges"

  }
  #ae <- 1000*evap   #[,"ae"]

  disp_output(title=title, disp.par=disp.par, runid=NULL,
                    qsim, ae, 1000*rain, tm=sel[1],
                    groups, qobs, #eff=eff,
                    ymax=max.q,

                                 # show to end of simulation period
                    timeBuffer= buf,
                    ichan=ichan,...)

}

# the beginning of the epoch is 1970, so coercing a POSIXct to numeric gives
# no. secs since 01-01-1970:00:00
asPOSIXct <- function(x, origin = "1970-01-01", tz = "GMT")
{
  return(as.POSIXct(x, origin = origin, tz = tz))

}


# returns the size of an time based object in seconds or ghours
length.period <- function(qr, units="hr")
{

	res <- diff(as.numeric(range(index(qr))))

	if(units=="hr")
	{
		res <- res/3600
	}
	return(res)
}

# Output results of a Dynamic TOPMODEL run. Wrapper to disp.output
plot.run <- function(run,
                     qsim=NULL,
                     fn=NULL,
                     start = run$sim.start,
                     end = run$sim.end,
                     par=disp.par(),
                     ...)
{
  sel <- paste0(start, "::", end)
  rain <- run$rain[sel]*1000
  evap <- run$ae[sel]*1000

  # observations, if supplied
  qobs <- run$qobs

  if(is.null(qsim))
  {
    # use the run output - can over
    qsim <- run$qsim[sel]*1000
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
  par$max.q <- max(c(par$max.q, qsim[]), na.rm=T)

  if(!is.null(qobs))
  {
    qobs <- qobs[sel]*1000
    try(print(NSE(qsim, qobs), silent=T))
    # update the limits
    par$max.q <- max(c(par$max.q, qobs), na.rm=T)
    cat("Observed time at peak = ", format(time_at_peak(qobs)), "\n")

    cat("Observed peak discharge = ", round(max(qobs),2), " mm/hr\n")
  }
  cat("Time at peak = ", format(time_at_peak(qsim)), "\n")
  cat("Peak discharge = ", round(max(qsim),2), " mm/hr\n")
  #	nresp<- ncol(qresp)

  # par(mar=c(4,4,3,4))
  if(length(fn)>0)
  {
    jpeg(filename=fn, ...)
    # larger axes labels
    par("cex.lab"=1.5)
    on.exit(dev.off())
    par(mar=c(3,4,3,4.5))
    title <- ""
  }

  #
  disp_output(qsim=qsim,
              #		start=sim.start,
              #		end=sim.end,
              evap=evap,
              rain=rain,
              tm=NULL,
              start=start,
              end=end,
              qobs=qobs,
              par=par, ...)

}

#' Display output of a Dynamic TOPMODEL  run
#'
#' @description Simple output of the results of a simulation.
#'
#' @details This will render the hydrograph, any observations, actual evapotranpiration, if supplied, and the rainfall hyetograph.
#' @export disp_output
#' @param qsim Time series of simulated discharges.
#' @param rain Time series of rainfall, at same interval as simulated values.
#' @param evap Time series of evapotranspiration (optional), at same interval as simulated values.
#' @param qobs Time series of evapotranspiration (optional), at same interval as simulated values.
#' @param start Start time for plot in a format interpretable as POSIXct date time. Defaults start of simulated discharges.
#' @param end End time for plot in a format interpretable as POSIXct date time. Defaults to end of simulated discharges.
#' @param tm Display a vertical line at this time in the simulation. If NULL no line will be drawn.
#' @param par Parameters controlling display output. A default set may be obtained through a call to disp.par.
#' @param ... Any further named parameters will be treated as graphics parameters and applied throughout the plot.
#' @author Peter Metcalfe
#' @seealso disp.par
#' @examples
#' \dontrun{
#' # Show the output of the storm simulation, overriding label colours and vertical axis limits.
#' require(dynatopmodel)
#'
#' data(brompton)
#'
#' x11()
#' with(brompton$storm.run, disp_output(evap=ae*1000,qobs=qobs*1000,
#'                                      qsim=qsim*1000, rain=rain*1000,
#'                                      max.q=2, cex.main=1, col.axis="slategrey", las.time=1))
#'}
disp_output <- function(qsim,
                        rain, evap=NULL,
                        qobs = NULL, tm=NULL,
                        par = NULL,
                        start=min(index(qsim)),
                        end=max(index(qsim)),...)
{
  if(is.null(par))
  {
    par <- disp.par()
  }
  par <- merge.lists(par, list(...))
  # save and restore display parameters
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))

  # any relevant graphical parameters supplied in disp.par will applied
  par(delete.NULLs(par[names(old.par)]))

  # add more space to right for axis titles, trim left margin a bit and shrink bottom margin
 # par("mgp"=c(2,1,0))
  par("xpd"=F)

  # if the time is specified then extend the time bounds around it
  if(!is.null(tm))
  {
    timeBuffer <- par$graphics.window.length # in hours
    # Window onto the desired period in the simulation, or all of it if limits not specified. Assume buffer in hours
    # rainfall and evapotranspiration on same plot (conversion to seconds for POSix time)
    start  <- max(c(tm-timeBuffer*3600/2, start, na.rm=T))
    end <- max(c(start+timeBuffer*3600, end, na.rm=T))

    # if past the end of the rainfall stop there
    if(end>max(index(rain)))
    {
      end <- max(index(rain))
    }
  }

  # inverted plot of rainfall
  # reversed limits, plus a buffer
 # ylim2 <- c(max(c(2*rain),na.rm=TRUE),min(Qr, 0, na.rm=T))
  sel <- paste0(start, "::", end)

  rain <- rain[sel]
  qriv <- qsim[sel]

  nobs <- 0
  cols <- par$col.qsim

  lty <- par$lty.qsim
  lwd <- par$lwd.qsim

  leg.txt <- expression(q[sim]) #"Simulated"
 	if(length(qobs)>0)
 	{
 	  # if observations add to output and colours plus vector of line widths and types used in disp.qsim
 		qobs <- qobs[sel]
 		qriv <- cbind(qriv, qobs)
 		nobs <- ncol(qobs)
 		cols <- c(cols, par$col.qobs[1:nobs])
 		lty <- c(lty, rep(par$lty.qobs, nobs))
 		lwd <- c(lwd, rep(par$lwd.qobs, nobs))
 		for(i in 1:nobs)
 		{
 		  leg.txt <- c(leg.txt, expression(q[obs]))
 	  }# "Observed"
 	}

  leg.txt <- c(leg.txt, "Precipitation")

  xlim <- range(index(qriv))

  # rain at top, inverted
  layout(matrix(1:2), heights=c(par$prop, (1-par$prop)))
  on.exit(layout(matrix(1)))
  # no margin at base
  with(par, par("mar"=c(0,xmar,ymar,xmar)))

  # inverted plot of rainfall (note limits inverted)
  disp.rain(rain,
            xlim = xlim,
            ylim = c(par$max.rain, 0),
            lwd=2,
            main=title)

  # axis is labelled if the position specified is at the top
  with(par, add_time_axis(side=3,
                las=las.time,
                labels=time.axis.side=="top",
                cex=cex,
                fmt=time.fmt,
                time.int=par$int.time))

  # border up the sides and top
  axis(side=3, col.ticks = NA, labels = F)
  axis(side=2, col.ticks = NA, labels = F)

  disp.pos <- !is.null(tm) && (tm>start) && (tm<end)

  if(disp.pos)
  {
    disp.sim.time(tm, label=F)
  }

  title(main=par$main)
  # calculate y limits from maximum of specified values evap, observed flows plotted on LH axis
  ylim <- c(min(qsim, 0, na.rm=T), max(c(0.01, par$max.q), na.rm=T))

#   par("cex.axis"=cex)
#   par("cex.lab"=cex)
  with(par, par("mar"=c(ymar,xmar,0,xmar)))

  evap <- evap[sel]
	if(length(evap)>0 && length(which(is.finite(evap)))>0)
	{
	  leg.txt <- c(leg.txt, expression(E[a]))
	  disp.evap(evap, ylim=ylim, col=par$col.evap)
	  par(new=T)
	}

	disp.qsim(qriv, xlim=xlim, ylim=ylim,
	          cols=cols,
					  qlab=par$lab.q,
						lty=lty,
						lwd=lwd,
						cex=par$cex,
						qint=as.numeric(par$int.q))

  # line showing current time
	if(disp.pos)
	{
	    disp.sim.time(tm, fmt=par$time.pos.fmt)
	}
	# U-shaped borders
  box(bty="u")

  # time axis lines
  with(par, add_time_axis(side=1,
                               las=las.time,
                               cex=cex,
                               # axis is labelled if the position specified is at the bottom
                               labels=time.axis.side=="bottom",
                               fmt=time.fmt,
                               time.int=par$int.time))

  if(par$legend.show)
  {
    legend.col <- c(cols, "black", par$col.evap)

    add.legend(legend.col=legend.col,
               legend=leg.txt,
               cex=par$cex,
               yoff=par$legend.yoff)

  }


}

disp.sim.time <- function(tm, label=T, fmt="%d-%b-%y")
{
    if(label)
    {
        time.str <- format(tm, fmt )
        mtext(side=1, at=tm, text=time.str, las=3,cex=0.8, line=0.25)
    }
    abline(v=tm, col="red",lwd=2)
}

add_time_axis <- function(side=1,
                          las=par("las"),
                          labels=T,
													time.int="week",
                          col="slategray", lty=2,
                          cex=par("cex.lab"),
													fmt="%d-%b-%y",...)
{
	#grid(col="slategray", nx=NA)
	#  horz axis above with perpendicular labels - more ticks than default- compute tick locations
	par("xaxp"=c(par("usr")[1:2], 1))
	# time axis at top, add formatted labels and lines

	# round_date is a lubridate method
	tm.range <- lubridate::round_date(as.POSIXct(par("xaxp")[1:2], origin="1970-01-01"), "day")
	tms <- seq(tm.range[1], tm.range[2], by=3600)		# hourly intervals

	if(is.numeric(time.int))
	{
	  i.at <- which(lubridate::hour(tms) %% time.int==0)
	}
	else
	{
  	# more detail. get pretty breaks based on day numbers (date-times in seconds from t.origin)
  	# use more lubridate methods
  	i.at <- switch(time.int,
  	"month"= which(lubridate::mday(tms)==1 & lubridate::hour(tms)==0),
  	"week"= which(lubridate::wday(tms)==1 & lubridate::hour(tms)==0),
  	"day"=which(lubridate::hour(tms)==0),
      "hour"=which(lubridate::second(tms)==0)
  	)
	}

	at <- tms[i.at]

	if(labels)
	{

		labs <- at  # as.POSIXct(at, origin="1970-01-01")
		labs <- format(labs, format=fmt)
        for(iside in side)
        {
		      axis(side=iside, at = at, las=las, labels=labs, cex.lab=cex, ...)
        }
	}
	abline(v=at, col=col, lty=lty)

}


# plot simulated discharges
disp.qsim <- function(qsim, xlim=range(index(qsim)),
											prop=1,
											ylim=c(0, max(qsim)),
											cols=qsim.cols(qsim),
											qint=0.5,
											qlab="Specific discharge (mm/hr)", ...)
{
  if(prop<1)
  {
      # fit to base of screen
      xflims <- par("fig")[1:2]
      par(fig=c(xflims, 0, prop))
      on.exit(par(fig=c(xflims,0,1)))
  }

  plot.zoo(qsim, main="", plot.type="single",
           xlim=xlim,
           ylim=ylim,
           ylab="", #Specific discharge / evap. (mm/hr)",
           col=cols,
           xaxt="n",
           yaxt="n",
           bty="n", xlab="", ...)

  add.q.axis(side=2, qsim, ylim, int=qint, lab=qlab)
  # horizontal grid lines up to max discharge
 #	ylim <- range(qsim)

  at <- seq(floor(min(qsim, 0, na.rm=T)), ceiling(ylim[2]), by=qint)#  pretty(par("usr")[3:4])
#  at <- at[-length(at)]
  abline(h=at, col="gray", lty=2)
  abline(h=0, col="gray")
  axis(2, at=at)
  # neater
  mtext(text=qlab, las=3, side=2, line=2.5,
        font=par("font.lab"),
  			cex=par("cex.lab"))
 #add_time_axis(side=1, las=2)

}

add.q.axis <- function(side=2, qsim, ylim=range(qsim), int=0.1, lab="", line=2.5, cex=par("cex.axis"),...)
{
  at <- seq(floor(min(qsim, 0, na.rm=T)), ceiling(ylim[2]), by=as.numeric(int))#  pretty(par("usr")[3:4])
  #  at <- at[-length(at)]
  abline(h=at, col="gray", lty=2)
  abline(h=0, col="gray")
  axis(side=side, at=at, cex.axis=cex)
  # neater
  mtext(text=lab, las=3, side=side, line=line,
        font=par("font.lab"),
        cex=cex) #par("cex.lab"))

}

# display actual evapotranspiration, in mm/hr. Axis is at least ymin high
disp.evap <- function(evap,
                      ylim=par("yaxp")[1:2],
                      ymin=1,    # min axis value in mm/hr
                      prop=1,text=expression(E[a]*"(mm/hr)"),
                      cex=par("cex.axis"),
											col="brown",
											...)
{
#	evap[evap==0]<- NA
  if(prop<1)
  {
	  par(fig=c(0,1, 0, prop))
	  on.exit(par(fig=c(0,1,0,1)))
  }
#	if(length(which()))
    plot.zoo(evap, type="h", ylim=ylim, col=col, #make.transparent("brown"),
             xaxt="n", bty="n",
             yaxt="n", main="", ann=F, lty=1, lwd=1, ylab="",xlab="", ...)
#	axis.range <- range(evap, na.rm=T)
  # ensure axis is always at least ymin high, greater if the evap exceeds this
	labs <- seq(0, ceiling(max(c(ymin, evap), na.rm=T)), by=0.2)
	axis(side=4, labels=labs, at=labs, cex.lab=cex)
	mtext(text=text, las=3,
				side=4, line=2.5, adj=0,
				cex=cex)
}

# inverted plot of rainfall
disp.rain <- function(rain, prop=1, col=rain.cols(rain),
                      xlim=range(index(rain)),
                      ylim=c(max(rain,na.rm=TRUE),0),
                      text="Precipitation (mm/hr)", axes=T, side=4,
                      cex=par("cex.axis"), line=2.5*cex,
                      ...)
{
    if(prop<1)
    {
        # place in upper part of screen
        xflims <- par("fig")[1:2]
        par(fig=c(xflims, 1-prop, 1))   #
        on.exit(par(fig=c(xflims,0,1)))
    }

  rain[rain==0] <- NA
  plot.zoo(rain,
           xlim=xlim,
           ylim=ylim, type="h", plot.type="single",
           ylab="",
           xaxt="n",
           yaxt="n",
           #cex.main=1,
           bty="n",
           lwd=2,
           col=col, #rain.cols(rain),
           xlab="")


	if(axes)
	{
		# right hand axis - precipitation and pe
		at <- pretty(range(rain, na.rm=T))
		axis(side=side, at=at, labels=at, cex.lab=cex)

		# grid lines
		abline(h=at, lty=3, col="gray")
		#     add_time_axis(side=3, labels=F)

		mtext(text=text, las=3,
					side=side, line=line,
			#		font=par("font.lab"),
					cex=cex)
	}

}

rain.cols <- function(rain)
{
  cols <- c("#000000","slategray","gray")
  return(cols[1:ncol(rain)])
}


qsim.cols <- function(qsim)
{
  cols <-   rev(c("#00FFFF","#00BFDF","#007FBF","#003F9F","#000080")) #colorRampPalette(c("navy", "cyan"), ncol(qsim))

  return(cols[1:ncol(qsim)])
}

add.legend <- function(nrow=2,
                       legend = expression("Simulated",
                                            "Observed",  "Precipitation",
                                            E[a]),
                       cex=par("cex.lab"),
                      legend.col = c("blue", "green", "black", "brown"),
                       legend.title=NULL,yoff=-0.05,...
                       )
{
  if(length(legend)>0)
  {
    #cols <- get.qobs.cols(qobs)
#     titles <- c(paste("Simulated: ", colnames(qsim)),
#                 paste("Rain:", colnames(rain)),
#                 paste("Observed:", colnames(qobs)),
#                 "Evapotranspiration")




    #  determine plot limits
    xlim <- par("usr")[1:2]
    # xjust and yjust controls how legend justified wrt x and ycoord: 2=right / top justified (doc is wrong)
    legend(x=xlim[2], y=yoff, legend=legend, #ncol=length(titles),
           title=legend.title,cex=cex,
           xpd=T,  # needed in order to plot outside figure margin
           yjust=2, xjust=1, horiz=T,
          # ncol=max(2, round(length(titles)/nrow+0.5)),
           col=legend.col, bg="white", lty=1, lwd=2)
    #bg="#FFFFFFCC")

  }
}



# barplot showing channel and reservoir storages
disp.chan.stores <- function(groups, stores, ichan)
{
	chan.stores <- stores[ichan,]
	chan.groups <- groups[ichan,]
	nchan <- length(ichan)
	# scale to percentage?
	chan.groups$sd_max <-0.01
	dat <- rbind(chan.groups$sd_max-chan.stores$sd, chan.stores$sd)
	cols <- rbind(rep("blue", nchan), rep("white", nchan), names.arg = chan.groups$id)
	barplot(horiz=T, dat, col=cols, xlab="Average depth (m)")

}

disp.stores <- function(groups,stores,ichan=1, nlev=5)
{
  # colours for unsat zone
  uzcols <- colorRampPalette(c("#A0D600FF", "blue"))(nlev)

  # move axis title
  par("mgp"=c(1,1,0))
  # add a bit more space at the bottom for the legend and HSU IDs
  par("mar"=c(4, 2, 1, 3))

  groups <- groups[-ichan,]  #rbind(groups[-ichan,], dumgr)
  stores <- stores[-ichan,]  #rbind(, dumst)

  ngroup <- nrow(groups)

  # stacked plot of:
  # max and actual root zone storage
  # unsat zone storage - shown by shading of unsat zone
  # max and actual storage deficit
  # colours for storage diagram
  # tan, brown, light blue, dark blue
  # the uz zone can contain at max the remaining storage deficit?
 # browser()
#   sds <- stores$sd
#   sds[sds==0]<-1
#   suzProp <- stores$srz/sds
  suzProp <- rep(1, ngroup)
  unsat <- which(stores$sd>0)

  uzcol <- uzcols[suzProp]
#   col <- c("red",     # storage excess
#         #   "#A6611A", # dark brown - wetted root zone
#            "#DFC27D", # tan - dried out root zone
#            "#80CDC1", # green - unsat
           col<-c("gray", # bed rock
                  colours()[125],  # gunmetal blue - saturated zone
                  colours()[85],   # khaki - unsat zone
                  "tan",    # dried root zone
                  "#A6611A",  # wetted rooT zone
                  "blue")   # surface storage
          # "#018571")  # blue - sat zone (water table)
 # col <- rev(col) # now show the bar the correct way up
  # brewer.pal(5, "BrBG")
  cols <-  rbind(col[1],
                col[2],
                uzcol,  # colour the uz according to amount of storage
                col[4]) #matrix(rep(col, ngroup),nrow=nstores)

#  densities <- rbind(0, 0, round(suzProp), 0)

  rz <- rbind(groups$srz_max-stores$srz,stores$srz)
  # if drying out then zero storage deficit
  drying <- which(stores$wetting==FALSE)
#   if(length(drying)>0)
#   {
#    #    browser("Stores drying out")
#     # if drying then the root zone dries from the top so wet areas are at bottom,
#     # if wetting-up then the base of  the zone wets last. Wett
#     rz[1:2,drying]<- rz[2:1,drying]
#
#     # if the region is "wetting" up then darker shade increases from the top downwards,
#     # if drying out, lighter shade increases from top
#     cols[1:2,drying] <- cols[2:1,drying]
#   }
  # max height of sat zone, including space for "bedrock"
  max.sz <- max(0.1+groups$sd_max)      #   max(rz + groups$sd_max, na.rm=T)
  # excess (overland) storage
  ex <- stores$ex

  # make columns the same height, draw a dotted line at maximum sd
  #  dat<-t(cbind(rz, stores$suz,stores$sd-stores$suz,groups$sd_max-stores$sd))
  dat<- rbind(as.vector(max.sz-groups$sd_max), as.vector(groups$sd_max-stores$sd),
              as.vector(stores$sd), rz, ex)  #  maxh--)


  ylim <- max(c(0,colSums(dat)),na.rm=T)

  widths <- sqrt(groups$area/sum(groups$area, na.rm=T))

  widths[is.na(widths)]<-0.2
  # shade unsaturated zone according to ammount of storage. value gives no. shading
  # lines per inch
  #densities <- rbind(NULL,NULL,20*stores$suz/stores$sd,NULL)
  barplot(dat, col=col,
         # main="Subsurface storages",
          #sub ="Specific storages",
          beside=FALSE,
          las=2,
          space=0,
         # density=densities,
          names.arg=groups$tag,
          xlab ="",
         ylab="",
          width=widths,
         # ylab="Storage (m)",
         cex.main=1, axes=F)
      #    ylim=c(ylim,0))

  mtext(side=2, text="Storage (m)", line=0)
 # axis(side=2, at=pretty(c(0, max(groups$srz_max))))

#   axis(side=2,
#        at=pretty(c(max(groups$srz_max), ylim)))  #,
#        labels=pretty(c(max(groups$srz_max), ylim))-max(groups$srz_max))
  # legend in bottom margin,centered vertically
  leg.y <- grconvertY(par("fin")[2] + par("mai")[3]/2, from="inches")  # par("usr")[4]
#leg.y <- par("usr")[4]
  legend(x=0, y=leg.y,
         ncol=4,
         bty="n",
      #   bg="white",
         xpd=T,
         legend=c("Sat. RZ", "Unsat RZ", "Unsat Zone", "Max SD"),cex=0.8, fill=col)
  #browser()
}

require(fields)


# returns the input ts subsetted by a time range
# either bound or neither can be supplied, if a bound is not supplied then the corresponding bound of the
# input series is used
# dt = time step in hours. If the data intervals are greater than these then
# disaggregate or aggregate if they are larger
GetTimeSeriesInputRange <- function(series, start, endt,
                                    cap="rainfall", cols=1,
                                    dt=1, verbose=TRUE)
{
    if(is.null(series)){return(series)}
    #if(length(series)==0){browser()}
    # subset of columns if specified
    if(length(which(is.finite(series[])))==0){return(series)}
    series <- series[,cols]
    res <- series[index(series)>=start & index(series)<=endt]
    nas <- which(is.na(res))
    if(length(nas)>0)
    {
        if(verbose==TRUE)
        {
            cat(paste(length(nas), " ", cap, " NA data encountered: replaced with 0\n"))
        }
       # res[nas] <- 0
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
    return(xts(res))
}

# make lists of time series for each hru and flux
collect_flux_ts <- function(groups, fluxes, tms)
{
    catch.area <- sum(groups$area)
    # convert fluxes to a named list
    fluxes <- apply(fluxes, MARGIN=3, function(x){list(x)})
    fluxes <- lapply(fluxes, function(x)x[[1]])
    names(fluxes) <- c("qbf", "qin", "uz", "rain", "ae", "ex", "qof")
    fluxes <- lapply(fluxes,
                     function(flux)
                     {
                         flux <- flux[1:length(tms),]
                         xts(cbind(flux, "mean"=rowSums(flux*groups$area)/catch.area),
                             order.by=tms)
                     })

    return(fluxes)
}


collect_storage_ts <- function(groups, storages, tms)
{
    catch.area <- sum(groups$area)

    storages <- apply(storages, MARGIN=3, function(x){list(x)})
    storages <- lapply(storages, function(x)x[[1]])
    names(storages) <- c("srz", "suz", "sd")


    storages <- lapply(storages,
                       function(storage)
                       {
                           storage <- storage[1:length(tms),]
                           xts(cbind(storage, "mean"=rowSums(storage*groups$area)/catch.area),
                               order.by=tms)

                       })
    return(storages)

}
