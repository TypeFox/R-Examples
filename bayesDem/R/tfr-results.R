TFRresults.group <- function(g, main.win, parent) {
	e <- new.env()
								
	e$sim.dir <- parent$sim.dir
	graph.defaults <- formals(png)
	mcmc.defaults <- formals(run.tfr.mcmc)
	nb <- bDem.gnotebook(container=g, expand=TRUE)
	show.traj.g <- ggroup(label="<span color='#0B6138'>TFR trajectories</span>", 
							markup=TRUE, horizontal=FALSE, container=nb)
	create.trajectories.group(show.traj.g, e, main.win)	
	
	map.g <- ggroup(label="<span color='#0B6138'>TFR world maps</span>", markup=TRUE, 
					horizontal=FALSE, container=nb)
	create.maps.group(map.g, e, main.win)

	dl.curve.g <- ggroup(label="<span color='#0B6138'>DL curve</span>", markup=TRUE, 
							horizontal=FALSE, container=nb)
	create.dlcurves.group(dl.curve.g, e, main.win)
	
	traces.g <- ggroup(label="<span color='#0B6138'>Parameter traces</span>", markup=TRUE, horizontal=FALSE, container=nb)
	create.partraces.group.all(traces.g, e, main.win)
	
	############################################
	# Convergence Diagnostics
	############################################
	convergence.g <- ggroup(label="<span color='#0B6138'>Convergence</span>", markup=TRUE, horizontal=FALSE, container=nb)
	create.convergence.group.all(convergence.g, e$sim.dir, main.win=main.win)	
	svalue(nb) <- 1
	return(e)
}

.create.trajectories.settings.group <- function(g, e, defaults, l=1) {
	leftcenter <- c(-1,0)
	show.traj.tfr.f <- gframe("<span color='blue'>Trajectories settings</span>", markup=TRUE, 
								horizontal=FALSE, container=g)
	lo <- glayout(container=show.traj.tfr.f) 
	lo[l,1, anchor=leftcenter] <- glabel('CI (%):', container=lo)
	lo[l,2] <- e$pi <- gedit('80, 95', width=7, container=lo)
	lo[l+1,1, anchor=leftcenter] <- glabel('# trajectories:', container=lo)
	lo[l+1,2] <- e$nr.traj <- gedit(20, width=6, container=lo)
	lo[l,3, anchor=leftcenter] <- 	glabel('From year:', container=lo)
	lo[l,4] <- e$start.year <- gedit(width=4, container=lo)
	lo[l+1,3, anchor=leftcenter] <- glabel('To year:', container=lo)
	lo[l+1,4] <- e$end.year <- gedit(width=4, container=lo)
	lo[l,5] <- glabel('     ', container=lo)
	lo[l,6:7] <- e$half.child.variant <- gcheckbox('+/- 0.5 child', checked=defaults$half.child.variant, 
								container=lo)
	lo[l+1,6:7] <- e$typical.trajectory <- gcheckbox('Typical trajectory', checked=defaults$typical.trajectory, 
								container=lo)
	return(lo)
}
#################################################################################

create.trajectories.group <- function(g, parent.env, main.win) {
	############################################
	# TFR Trajectories
	############################################
	e <- new.env()
	e$sim.dir <- parent.env$sim.dir
	e$pred.type <- 'tfr'

	defaults.pred <- formals(tfr.predict)
	defaults.traj <- formals(tfr.trajectories.plot)
	defaults.traj.all <- formals(tfr.trajectories.plot.all)
	mcmc.defaults <- formals(run.tfr.mcmc)
	
	addSpace(g, 10)
	show.traj.country.f <- gframe("<span color='blue'>Country settings</span>", markup=TRUE, 
									horizontal=FALSE, container=g)
	e$show.traj.country <- create.country.widget(show.traj.country.f, defaults.traj.all, 
									main.win, prediction=TRUE, parent.env=e)
		
	addSpace(g, 10)
	.create.trajectories.settings.group(g, e, defaults=list(start.year=mcmc.defaults$start.year, end.year=defaults.pred$end.year,
															half.child.variant=defaults.traj$half.child.variant,
															typical.trajectory=defaults.traj$typical.trajectory))
	addSpace(g, 10)
	
	show.traj.graph.f <- gframe("<span color='blue'>Advanced graph parameters</span>", markup=TRUE, 
									horizontal=FALSE, container=g)
	e$graph.pars <- create.graph.pars.widgets(show.traj.graph.f, main.win=main.win)
	addSpring(g)
	show.traj.bg <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='tfr.trajectories.plot', package='bayesTFR', parent.group=show.traj.bg,
						parent.window=main.win)
	addSpring(show.traj.bg)
	create.generate.script.button(handler=show.e0.traj, action=list(mw=main.win, env=e, type='plot', 
									script=TRUE, pred.type='tfr', package='bayesTFR'),
								container=show.traj.bg)
	addSpace(show.traj.bg, 5)
	TableB.show.traj.act <- gaction(label='Table', icon='dataframe', handler=show.e0.traj, 
						action=list(mw=main.win, env=e, type='table', script=FALSE, pred.type='tfr', package='bayesTFR'))
	#GraphB.show.traj.act <- gaction(label='Graph', icon='lines', handler=showTFRtraj, 
	#					action=list(mw=main.win, env=e, type='plot', script=FALSE))
	GraphB.show.traj.act <- gaction(label='Graph', icon='lines', handler=show.e0.traj, 
						action=list(mw=main.win, env=e, type='plot', script=FALSE, pred.type='tfr', package='bayesTFR'))
	e$TableB.show.traj <- bDem.gbutton(action=TableB.show.traj.act, container=show.traj.bg)
	bDem.gbutton(action=GraphB.show.traj.act, container=show.traj.bg)
}

.create.map.settings.group <- function(g, e, measures=c('TFR', 'lambda', tfr.parameter.names.cs.extended())) {
	leftcenter <- c(-1,0)
	map.set.f <- gframe("<span color='blue'>Map settings</span>", markup=TRUE, 
									horizontal=FALSE, container=g)
	mlo <- glayout(container=map.set.f)
	mlo[1,1, anchor = leftcenter] <- glabel('Percentile:', container=mlo)
	e$percentiles <- list('median'=0.5, 'lower 80'=0.1, 'upper 80'=0.9, 'lower 90'=0.05, 'upper 90'=0.95,
						'lower 95'=0.025, 'upper 95'=0.975, 'lower 60'=0.2, 'upper 60'=0.8,
						'lower 50'=0.25, 'upper 50'=0.75, 'lower 40'=0.3, 'upper 40'=0.7, 
						'lower 20'=0.4, 'upper 20'=0.6
						)
	mlo[1,2] <- e$map.percentile <- bDem.gdroplist(names(e$percentiles), container=mlo)
	mlo[1,3] <- '    ' # add some space between the two groups
	mlo[2,1, anchor = leftcenter] <- glabel('Bounds:    ', container=mlo)
	mlo[2,2] <- bounds.g <- ggroup(horizontal=TRUE, container=mlo)
	e$map.bounds <- bDem.gdroplist(c(80, 90, 95, 60, 50, 40, 20), container=bounds.g)
	glabel('%', container=bounds.g)	
	mlo[3,1, anchor = leftcenter] <- glabel('Measure:', container=mlo)
	mlo[3,2] <- e$map.measure <- bDem.gdroplist(measures, container=mlo)	
	mlo[1,4, anchor = leftcenter] <- glabel('Use R package:', container=mlo)
	mlo[1:2,5] <- e$map.package <- gradio(c('rworldmap', 'googleVis'), horizontal = FALSE, 
						handler=function(h, ...) {
							enabled(e$map.bounds) <- svalue(h$obj) == 'googleVis';
							enabled(e$same.scale) <- svalue(h$obj) == 'rworldmap'}, 
						container=mlo)
	mlo[3,4:5] <- e$same.scale <- gcheckbox('Same scale for all maps', checked=TRUE, container=mlo)
	enabled(e$map.bounds) <- svalue(e$map.package) == 'googleVis'
	enabled(e$same.scale) <- svalue(e$map.package) == 'rworldmap'
	return(mlo)
}

create.maps.group <- function(g, e, main.win) {
	############################################
	# TFR World Maps
	############################################
	addSpace(g, 10)
	.create.map.settings.group(g, e)
	addSpring(g)
	map.bg <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='tfr.map', package='bayesTFR', parent.group=map.bg,
						parent.window=main.win)
	addSpring(map.bg)
	create.generate.script.button(handler=showMap, action=list(mw=main.win, env=e, script=TRUE),
								container=map.bg)
	addSpace(map.bg, 5)
	GraphB.map <- gaction(label=' Show Map ', handler=showMap, 
						action=list(mw=main.win, env=e, script=FALSE))
	bDem.gbutton(action=GraphB.map, container=map.bg)
	addHandlerChanged(e$map.measure, 
					handler=function(h,...) enabled(e$same.scale) <- svalue(h$obj) == 'TFR')
}

create.dlcurves.group <- function(g, parent.env, main.win) {
	e <- new.env()
	e$sim.dir <- parent.env$sim.dir
	e$pred.type <- 'tfr'
	leftcenter <- c(-1,0)
	############################################
	# DL Curves
	############################################
	defaults.dl <- formals(DLcurve.plot)
	defaults.dl.all <- formals(DLcurve.plot.all)
	addSpace(g, 10)
	dlc.country.f <- gframe("<span color='blue'>Country settings</span>", markup=TRUE, 
							horizontal=FALSE, container=g)
	e$dlc.country <- create.country.widget(dlc.country.f, defaults.dl.all, main.win, prediction=FALSE, 
											parent.env=e, disable.table.button=FALSE)
	addSpace(g, 10)
	dlc.dl.f <- gframe("<span color='blue'>DL curve settings</span>", markup=TRUE, 
							horizontal=FALSE, container=g)
	dllo <- glayout(horizontal=TRUE, container=dlc.dl.f)
	dllo[1,1, anchor=leftcenter] <- glabel('CI (%):', container=dllo)
	dllo[1,2] <- e$pi <- gedit('80, 95', width=7, container=dllo)
	dllo[1,3, anchor=leftcenter] <- glabel('Burnin:', container=dllo)
	dllo[1,4] <- e$burnin <- gedit(defaults.dl$burnin, width=5, container=dllo)
	dllo[2,3, anchor=leftcenter] <- glabel('Maximum TFR:', container=dllo)
	dllo[2,4] <- e$tfr.max <- gedit(defaults.dl$tfr.max, width=2, container=dllo)
	dllo[2,1, anchor=leftcenter] <- glabel('# curves:', container=dllo)
	dllo[2,2] <- e$nr.curves <- gedit(20, width=6, container=dllo)
	
	dllo[1,5] <- e$predictive.distr <- gcheckbox('Predictive distribution', 
							checked=defaults.dl$predictive.distr, container=dllo)
	addSpace(g, 10)			
	dlc.graph.f <- gframe("<span color='blue'>Advanced graph parameters</span>", markup=TRUE, 
						horizontal=FALSE, container=g)
	e$graph.pars <- create.graph.pars.widgets(dlc.graph.f, main.win=main.win)
	addSpring(g)
	dlc.bg <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='DLcurve.plot', package='bayesTFR', parent.group=dlc.bg,
						parent.window=main.win)
	addSpring(dlc.bg)
	create.generate.script.button(handler=showDLcurve, action=list(mw=main.win, env=e, script=TRUE),
								container=dlc.bg)
	addSpace(dlc.bg, 5)
	GraphB.dlc <- gaction(label='Graph', icon='lines', handler=showDLcurve, 
						action=list(mw=main.win, env=e, script=FALSE))
	bDem.gbutton(action=GraphB.dlc, container=dlc.bg)
}

.create.partraces.settings.group <- function(g, e, main.win, par.names, par.names.cs) {
	leftcenter <- c(-1,0)
	addSpace(g, 10)	
	f <- gframe("<span color='blue'>Parameter traces settings</span>", markup=TRUE, 
							horizontal=FALSE, container=g)
	tlo <- glayout(container=f)
	tlo[1,1, anchor=leftcenter] <- glabel('Parameters:', container=tlo)
	tlo[1,2] <- e$pars.chb <- gcheckbox("all", container=tlo, checked=TRUE, 
							handler=function(h,...) {
								if(svalue(e$cs.chb, index=TRUE)==2) {
									enabled(e$par.cs.dl)<-!svalue(h$obj)
									enabled(e$par.dl)<-FALSE
								} else {
									enabled(e$par.dl)<-!svalue(h$obj)
									enabled(e$par.cs.dl)<-FALSE
									}})
	tlo[2:3,1:2] <- e$cs.chb <- gradio(c('World parameters', 'Country specific'), horizontal = FALSE,
						container=tlo, 
						handler=function(h,...) {
								if (svalue(h$obj, index=TRUE)==2) {
									enabled(e$par.cs.dl)<-!svalue(e$pars.chb)
									enabled(e$par.dl)<-FALSE
								} else {
									enabled(e$par.dl)<-!svalue(e$pars.chb)
									enabled(e$par.cs.dl)<-FALSE
									}
								enabled(e$country$country.w) <- svalue(h$obj, index=TRUE)==2
								enabled(e$country$country.select.b) <- svalue(h$obj, index=TRUE)==2
								}
							)
	tlo[2,3:4] <- e$par.dl <- bDem.gdroplist(par.names, container=tlo)
	enabled(e$par.dl) <- FALSE
	tlo[3,3:4] <- e$par.cs.dl <- bDem.gdroplist(par.names.cs, container=tlo)
	enabled(e$par.cs.dl) <- FALSE
	tlo[4,1:4] <- cw <- ggroup(horizontal=TRUE, container=tlo)
	e$country <- create.country.widget(cw,  main.win=main.win, show.all=FALSE, prediction=FALSE, 
											parent.env=e)
	enabled(e$country$country.w) <- svalue(e$cs.chb, index=TRUE)==2
	enabled(e$country$country.select.b) <- svalue(e$cs.chb, index=TRUE)==2
	tlo[1,5] <- '    '
	tlo[2,6, anchor=leftcenter] <- glabel('# points:', container=tlo)
	tlo[2,7] <- e$nr.points <- gedit(100, width=5, container=tlo)
	tlo[3,6, anchor=leftcenter] <- glabel("Burnin:", container=tlo)
	tlo[3,7] <- e$burnin <- gedit(0, width=5, container=tlo)
	tlo[4,6, anchor=leftcenter] <- glabel("Thin:", container=tlo)
	tlo[4,7] <- e$thin <- gedit(1, width=5, container=tlo)
}

create.partraces.group.all <- function(g, parent.env, main.win) {
	type.nb <- gnotebook(container=g, expand=TRUE)
	set.widget.bgcolor(type.nb, color.main)
	set.widget.basecolor(type.nb, color.nb.inactive)
	phaseII.g <- ggroup(label="<span color='darkred'>Phase II</span>", markup=TRUE, horizontal=FALSE, container=type.nb)
	addSpace(phaseII.g, 10)
	create.partraces.group(phaseII.g, parent.env, main.win)
	phaseIII.g <- ggroup(label="<span color='darkred'>Phase III</span>", markup=TRUE, horizontal=FALSE, container=type.nb)
	addSpace(phaseIII.g, 10)
	create.partraces.group(phaseIII.g, parent.env, main.win, type='tfr3')
	svalue(type.nb) <- 1
}

create.partraces.group <- function(g, parent.env, main.win, type='tfr') {
	############################################
	# Parameter Traces
	############################################
	e <- new.env()
	e$sim.dir <- parent.env$sim.dir
	e$pred.type <- type
	.create.partraces.settings.group(g, e, main.win, 
			par.names=do.call(paste(type, '.parameter.names', sep=''), list()), 
			par.names.cs=do.call(paste(type, '.parameter.names.cs', sep=''), list()))
	addSpring(g)
	traces.bg <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic=paste(type, '.partraces.plot', sep=''), package='bayesTFR', parent.group=traces.bg,
						parent.window=main.win)	
	addSpring(traces.bg)
	SummaryB.traces <- gaction(label='Show summary', icon='dataframe', handler=showParTraces, 
						action=list(mw=main.win, env=e, print.summary=TRUE))
	bDem.gbutton(action=SummaryB.traces, container=traces.bg)
	GraphB.traces <- gaction(label='Graph', icon='lines', handler=showParTraces, 
						action=list(mw=main.win, env=e, print.summary=FALSE))
	bDem.gbutton(action=GraphB.traces, container=traces.bg)
}

create.convergence.group.all <- function(g, sim.dir, main.win) {
	type.nb <- gnotebook(container=g, expand=TRUE)
	set.widget.bgcolor(type.nb, color.main)
	set.widget.basecolor(type.nb, color.nb.inactive)
	phaseII.g <- ggroup(label="<span color='darkred'>Phase II</span>", markup=TRUE, horizontal=FALSE, container=type.nb)
	addSpace(phaseII.g, 10)
	create.convergence.tab(phaseII.g, sim.dir, main.win=main.win)
	phaseIII.g <- ggroup(label="<span color='darkred'>Phase III</span>", markup=TRUE, horizontal=FALSE, container=type.nb)
	addSpace(phaseIII.g, 10)
	e3 <- create.convergence.tab(phaseIII.g, sim.dir, type='tfr3', main.win=main.win)
	svalue(e3$keep.thin.mcmc) <- FALSE
	enabled(e3$keep.thin.mcmc) <- FALSE
	svalue(type.nb) <- 1
}

create.convergence.tab <- function(parent, sim.dir, type='tfr', package='bayesTFR', main.win=NULL) {
	defaults <- formals(paste(type,'.diagnose', sep=''))
	e <- new.env()
	e$sim.dir <- sim.dir
	leftcenter <- c(-1,0)
	addSpace(parent, 10)
	#g <- ggroup(horizontal=FALSE, container=parent, expand=FALSE)
	g <- glayout(container=parent)
	#g1 <- ggroup(horizontal=TRUE, container=g, expand=TRUE)
	g[1,1] <- g2 <- gframe("<span color='blue'>Diagnostics settings</span>", markup=TRUE, horizontal=FALSE, container=g)
	mclo <- glayout(container=g2)
	mclo[1,1, anchor=leftcenter] <- glabel('Thin:', container=mclo)
	mclo[1,2] <- e$thin <- gedit(defaults$thin, width=5, container=mclo)
	mclo[2,1, anchor=leftcenter] <- glabel('Burnin:', container=mclo)
	mclo[2,2] <- e$burnin <- gedit(defaults$burnin, width=5, container=mclo)
	mclo[3,1:2] <- e$keep.thin.mcmc <- gcheckbox('Keep thinned MCMCs', checked = defaults$keep.thin.mcmc, container=mclo)
	#g3  <- gframe("<span color='blue'>Optional settings</span>", markup=TRUE, horizontal=FALSE, container=g1)
	#oplo <- glayout(container=g3)
	mclo[1,3] <- e$express <- gcheckbox('Express', checked=defaults$express, container=mclo,
						handler=function(h,...) enabled(e$country.sampling.prop) <- !svalue(h$obj))
	mclo[1,4] <- e$verbose <- gcheckbox('Verbose', checked = defaults$verbose, container=mclo)
	
	mclo[2,3, anchor=leftcenter] <- glabel('Proportion of countries included (0-1):', container=mclo)
	mclo[2,4] <- e$country.sampling.prop <- gedit(1, width=5, container=mclo)
	
	#addSpace(parent, 10)
	g[2, 1] <- '  '
	g[3,1] <- bDem.gbutton('    Show available convergence diagnostics    ', container=g, handler=showConvergenceDiag,
					action=list(mw=main.win, env=e, type=type), fill=TRUE)
	# This line is commented out because (for some reason) on Windows OS it causes the main window to go out of whack.
	#glo[l,2] <- (e$country.sampling.prop <- gslider(from=0, to=1, by=1/200, value=1, container=g4))
	addSpace(parent, 10)
	.create.status.label(parent, e)
	enabled(e$country.sampling.prop) <- !defaults$express
	
	addSpring(parent)
	butg <- ggroup(horizontal=TRUE, container=parent)
	create.help.button(topic=paste(type,'.diagnose', sep=''), package=package, parent.group=butg,
						parent.window=main.win)	
	addSpring(butg)
	create.generate.script.button(handler=computeConvergenceDiag, 
					action=list(mw=main.win, env=e, type=type, package=package, script=TRUE),
								container=butg)
	addSpace(butg, 5)
	ComputeDiag <- gaction(label=' Compute New Diagnostics ', icon='execute', handler=computeConvergenceDiag, 
						action=list(mw=main.win, env=e, type=type, package=package, script=FALSE))
	bDem.gbutton(action=ComputeDiag, container=butg)
	return(e)
}

create.country.widget <- function(parent, defaults=NULL, main.win=NULL, glo = NULL, start.row=1, show.all=TRUE, 
									prediction=FALSE, parent.env=NULL, disable.table.button=TRUE) {
	e <- new.env()
	e$parent.env <- parent.env
	e$prediction <- prediction
	leftcenter <- c(-1, 0)
	rightcenter <- c(1, 0)
	g1 <- if(is.null(glo)) glayout(container=parent) else glo
	l <- start.row
	g1[l,1, anchor=leftcenter] <- glabel('Country:', container=g1)
	g1[l,2] <- e$country.w <- gedit(width=20, container=g1)
	g1[l,3] <- e$country.select.b <- bDem.gbutton('Select', container=g1, handler=selectCountryMenu,
								action=list(mw=main.win, text.widget=e$country.w, env=e))
	if(show.all) {
		l <- l+1
		g1[l,1:3] <- gseparator(container=g1)
		l <- l+1
		g1[l,1] <- e$all.countries.chb <- gcheckbox('All countries', checked=FALSE, container=g1,
									handler=function(h,...){
										enabled(e$country.w) <- !svalue(h$obj)
										enabled(e$country.select.b) <- !svalue(h$obj)
										enabled(e$all.type) <- svalue(h$obj)
										enabled(e$all.output) <- svalue(h$obj)
										if(disable.table.button) enabled(parent.env$TableB.show.traj) <- !svalue(h$obj)
										})
		g1[l,2, anchor=rightcenter] <- glabel("Output type:", container=g1)
		g1[l,3] <- e$all.type <- bDem.gdroplist(c('png', 'jpeg', 'pdf', 'tiff', 'bmp', 'postscript'), container=g1)
		enabled(e$all.type) <- FALSE
		l <- l+1
		g1[l,1, anchor=leftcenter] <- glabel("Output directory:", container=g1)
		g1[l,2:3] <- e$all.output <- bDem.gfilebrowse(eval(defaults$output.dir), type='selectdir', 
					  width=39, quote=FALSE, container=g1)
		enabled(e$all.output) <- FALSE
	}
	e$country.lo <- g1
	return(e)	
}



create.graph.pars.widgets <- function (parent, main.win=NULL) {
	g <- ggroup(horizontal=FALSE, container=parent)
	glabel("Comma-separated parameters in R format (see '?par'), e.g. ylim=c(0,6), xlab='Time'", container=g)
	pars.w <- gedit('', width=40, container=g)
	return(pars.w)
}

get.table.of.countries.from.meta <- function(sim.dir, prediction=FALSE, sorted=TRUE, 
										pred.type='tfr', env=NULL) {
	if(!is.null(env$prior.select.countries.function)) # function to run prior the selection. Can be used to set something to the env.
		do.call(env$prior.select.countries.function, list(env))
	if(prediction) {
		pred.call <- paste('get.', pred.type, '.prediction', sep='')
		args <- formals(pred.call)
		args$sim.dir <- NULL
		lenv <- as.list(env)
		add.args <- lenv[names(args)[!sapply(lenv[names(args)], is.null)]]
		pred <- do.call(pred.call, c(list(sim.dir=sim.dir), add.args))
		if(is.null(pred)) {
			gmessage('Simulation directory contains no valid predictions.', 
					title='Input Error', icon='error')
			return(NULL)
		}
		loc.data <- get.countries.table(pred)
	} else { #simulation
		mcmc.set <- do.call(paste('get.', pred.type, '.mcmc', sep=''), list(sim.dir=sim.dir))
		if(is.null(mcmc.set)) {
			gmessage('Simulation directory contains no valid MCMC results.', title='Input Error',
					icon='error')
			return(NULL)
		}
		loc.data <- get.countries.table(mcmc.set)	}
	if(sorted) {
		ord.idx <- order(loc.data[,'name'])
		loc.data <- loc.data[ord.idx,]
	}
	return(loc.data)
}

draw.new.country.select <- function(used, env) {
	for(item in c('sim.dir', env$new.country.select.if.changed))
		if (svalue(env[[item]]) != used[[item]]) return (TRUE)
	return(FALSE)
}
	
set.used.items <- function(env.used, env) {
	for(item in c('sim.dir', env$new.country.select.if.changed))
		env.used[[item]] <- svalue(env[[item]])
}

selectCountryMenu <- function(h, ...) {
	country.selected <- function(h1, ...) {
		selected.country <- as.numeric(svalue(h$action$env$selcountry.gt))
		selected.country <- get.country.object(selected.country, 
								country.table=h$action$env$country.table)
		if (length(selected.country) > 0) {
			svalue(h$action$text.widget) <- selected.country$name
		} 
		visible(h$action$env$country.sel.win) <- FALSE
	}
	new.window <- TRUE
	if (!is.null(h$action$env$country.sel.win)) {
		# if anything has changed (sim.dir or the data), the window needs to be re-built
		if(draw.new.country.select(h$action$env$used, h$action$env$parent.env)) {
			dispose(h$action$env$country.sel.win)
			new.window <- TRUE
		} else {
			country.table <- get.table.of.countries.from.meta(svalue(h$action$env$parent.env$sim.dir), 
								prediction=h$action$env$prediction, 
								pred.type=if(is.null(h$action$env$parent.env$pred.type)) 'tfr' 
											else h$action$env$parent.env$pred.type, 
								env=h$action$env$parent.env)
			if(is.null(country.table)) {
				dispose(h$action$env$country.sel.win)
				return(NULL)
			}
			if(dim(country.table)[1] != dim(h$action$env$country.table)[1]) {
				dispose(h$action$env$country.sel.win)
				new.window <- TRUE
			} else {
				new.window <- FALSE
				visible(h$action$env$country.sel.win) <- TRUE
			}
		}
	}
	if(new.window) {
		sim.dir.used <- svalue(h$action$env$parent.env$sim.dir)
		country.table <- get.table.of.countries.from.meta(sim.dir.used, prediction=h$action$env$prediction,
										pred.type=if(is.null(h$action$env$parent.env$pred.type)) 'tfr' 
											else h$action$env$parent.env$pred.type, 
										env=h$action$env$parent.env)
		if (is.null(country.table)) return(NULL)
		h$action$env$used <- new.env()
		set.used.items(h$action$env$used, h$action$env$parent.env)
		h$action$env$country.table <- country.table
		h$action$env$country.sel.win <- win <- gwindow('Select country', parent=h$action$mw, height=450,
					handler=function(h, ...) {
						h$action$env$country.sel.win<-NULL;
						h$action$env$selcountry.ok.handler <- NULL;
						h$action$env$selcountry.gt.handler <- NULL
					},
					action=list(env=h$action$env))
		t.group <- ggroup(horizontal=FALSE, container=win)
		h$action$env$selcountry.gt <- gtable(h$action$env$country.table, container=t.group, expand=TRUE,				handler=country.selected)
		b.group <- ggroup(horizontal=TRUE, container=t.group)
		gbutton('Cancel', container=b.group, handler=function(h, ...) 
					visible(win) <- FALSE)
		addSpring(b.group)
		h$action$env$selcountry.okbutton <- gbutton('OK', container=b.group)
	}
	if(!is.null(h$action$env$selcountry.ok.handler)) 
		removehandler(h$action$env$selcountry.okbutton, h$action$env$selcountry.ok.handler)
	h$action$env$selcountry.ok.handler <- addhandlerclicked(h$action$env$selcountry.okbutton, 
												handler=country.selected)
	if(!is.null(h$action$env$selcountry.gt.handler)) 
		removehandler(h$action$env$selcountry.gt, h$action$env$selcountry.gt.handler)
	h$action$env$selcountry.gt.handler <- addhandlerdoubleclick(h$action$env$selcountry.gt, 
												handler=country.selected)
}
	


get.country.code.from.widget <- function(country.widget, env, force.country.spec=FALSE, allow.null.country=FALSE) {
	country <- svalue(country.widget)
	country.selected <- TRUE
	if (!is.null(env$all.countries.chb)) { 
		if(svalue(env$all.countries.chb)) country.selected <- FALSE
	}
	if (force.country.spec) country.selected <- TRUE
	if (country.selected) {
		if (nchar(country)==0) {
			if(!allow.null.country)
				gmessage('Country must be specified.', title='Input Error', icon='error')
			return(NULL)
		}
		warn <- getOption('warn')
		options(warn=-1)
		country.code <- as.numeric(country)
		if (!is.na(country.code)) country <- country.code
		options(warn=warn)
		country <- get.country.object(country, country.table=env$country.table)
		if(is.null(country$name)) {
			gmessage('Country does not exist.', title='Input Error',
					icon='error')
			return(NULL)
		}
		return(country)
	}
	return(list(code=NULL, output.dir=svalue(env$all.output), output.type=svalue(env$all.type)))
}
	
get.additional.tfr.param <- function(e, ...) {
	hchv <- svalue(e$half.child.variant)
	return(list(add=list(half.child.variant=hchv), 
				plot=c('pi', 'xlim', 'nr.traj', 'half.child.variant', 'typical.trajectory'), 
				table=c('pi', 'country', 'half.child.variant'), table.decimal=2))
}
	
assemble.tfr.plot.cmd <- function(param, e, all=FALSE) {
	all.suffix <- if(all) '.all' else ''
	return(paste('tfr.trajectories.plot',all.suffix, '(pred,', assemble.arguments(param, svalue(e$graph.pars)), ')', sep=''))
}

get.tfr.table.title <- function(country, pred, ...) 
	return (country)

tfr.get.trajectories.table.values <- function(pred, param, ...) {
	t <- do.call('tfr.trajectories.table', c(list(pred), param))
	# change the column names, otherwise gtable will prefix an 'X'
	tend <- ncol(t)
	if(param$half.child.variant) tend <- tend-2
	colnames(t)[2:tend] <- paste('q', colnames(t)[2:tend], sep='')
	if(param$half.child.variant)
		colnames(t)[(tend+1):(tend+2)] <- c('minus0.5child', 'plus0.5child')
	return(t)
}

showMap <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	percentile <- svalue(e$map.percentile)
	quantile <- e$percentiles[[percentile]]
	param.env <-list(sim.dir=svalue(e$sim.dir), quantile=quantile)
	param.names1 <- list(text='sim.dir')
	param.pred <- get.parameters(param.names1, env=param.env, quote=h$action$script, retrieve.from.widgets=FALSE)
	same.scale <- svalue(e$same.scale)
	par.name <- svalue(e$map.measure)
	bounds <- svalue(e$map.bounds)
	package <- svalue(e$map.package)
	map.function <- if(package == 'rworldmap') 'tfr.map' else 'tfr.map.gvis'
	if(h$action$script) {
		cmd <- paste('pred <- get.tfr.prediction(', assemble.arguments(param.pred), 
						')\n', sep='')
		if (par.name == 'TFR') {
			 if(package == 'rworldmap') {
				cmd <- paste(cmd, "param.map <- get.tfr.map.parameters(pred, same.scale=", same.scale,
					", quantile=", quantile, ")\n", sep="")
				cmd <- paste(cmd, 'do.call("', map.function, '", param.map)', sep='')
			} else {
				cmd <- paste(cmd, map.function, '(pred, quantile=', quantile, ', pi=', bounds, ')', sep='')
			}
		} else {
			cmd <- paste(cmd, map.function, '(pred, quantile=', quantile, ', par.name="', par.name, '"', sep='')
			cmd <- paste(cmd, if (package == 'googleVis') paste(', pi=', bounds, sep='') else '', sep='')
			cmd <- paste(cmd, if (par.name == 'lambda' && package == 'rworldmap') 
						', catMethod="pretty",  numCats=20' else '', ')', sep='')
		}
		create.script.widget(cmd, h$action$mw, package="bayesTFR")
	} else {
		pred <- do.call('get.tfr.prediction', param.pred)
		if (par.name == 'TFR' && package == 'rworldmap') {
			param.map <-  get.tfr.map.parameters(pred, same.scale=same.scale, quantile=quantile)
		} else {
			param.map <- list(pred=pred, quantile=quantile)
			if (par.name != 'TFR')
				param.map[['par.name']]<- par.name
				if(par.name=='lambda' && package == 'rworldmap') 
					param.map <- c(param.map, list(catMethod='pretty',  numCats=20))
		}
		if(package == 'rworldmap') param.map[['device']] <- 'dev.new'
		if (package == 'googleVis') param.map[['pi']] <- bounds
		g <- create.graphics.map.window(parent=h$action$mw, pred=pred, params=param.map, percentile=percentile, 
										is.gvis= package == 'googleVis', title="World Map", 
										cw.main=paste(c(par.name, percentile), collapse=', '))
	}
}

tfr.get.time.info <- function(pred) {
	meta <- pred$mcmc.set$meta
	return(list(est.periods=bayesTFR:::get.tfr.periods(meta), 
				proj.periods=bayesTFR:::get.prediction.periods(meta, pred$nr.projections+1),
				proj.ind.func=bayesTFR:::get.predORest.year.index,
				present.year=meta$present.year,
				est.years=bayesTFR:::get.estimation.years(meta),
				proj.years=bayesTFR:::get.all.prediction.years(pred)
				))
	
}

e0.get.time.info <- function(pred) return(tfr.get.time.info(pred))
	
create.graphics.map.window <- function(parent, pred, params, percentile,  is.gvis=FALSE, title='', type='tfr', 
											cw.main='', dpi=80) {
	time.info <- do.call(paste(type, '.get.time.info', sep=''), list(pred))
	est.periods <- time.info$est.periods
	proj.periods <- time.info$proj.periods
	
	newMap <- function(h, ...) {
		if (!is.null(h$action$dev)) dev.set(h$action$dev)
		if(!is.null(h$action$map.pars$device)) h$action$map.pars$device <- "dev.cur"
		do.show.map(as.numeric(svalue(proj.year)), h$action$map.pars)
	}
	do.show.map <- function(projection.year, map.pars, update.control.win=TRUE) {
		#is.median <- percentile == 'median'
		ind.proj <- do.call(time.info$proj.ind.func, list(pred, projection.year))
		#projection.index <- ind.proj['index']
		is.projection <- ind.proj['is.projection']
		if(update.control.win)
			svalue(year.label) <- if(is.projection) 'Projection year:' else 'Estimation year:'
		do.call(paste(type, '.map', if(is.gvis) '.gvis' else '', sep=''), 
				c(map.pars, list(year=projection.year #, main=main
					)))
	}
	close.map <- function(h, ...) dev.off(h$action$dev)
	
	if(is.gvis && !is.null(params[['par.name']])) {
		do.show.map(time.info$present.year, params, update.control.win=FALSE)
		return(NULL)
	}
	lest.periods <- length(est.periods)
	periods <- c(est.periods[-lest.periods], # remove the present period, otherwise doubled 
				 proj.periods)
	est.years <- time.info$est.years
	years <- c(est.years[-lest.periods], time.info$proj.years)
	e <- new.env()
	win <- bDem.gwindow(paste(title, 'Control Panel'), height=70, parent=parent, horizontal=FALSE)
	g <- ggroup(container=win, horizontal=FALSE, expand=TRUE)
	glabel(cw.main, container=g)
	g1 <- ggroup(container=g, horizontal=TRUE)
	year.label <- glabel("Projection year:", container=g1)
	proj.year <- gspinbutton(from= min(years), to=max(years), by=5, value=years[lest.periods], container=g1)
	if(!is.null(params$par.name)) enabled(proj.year) <- FALSE
	do.show.map(time.info$present.year, params)
	if (!is.gvis) {
		addSpring(g1)
		glabel("Output type:", container=g1)
		e$type <- bDem.gdroplist(c("pdf", "postscript", "png", "jpeg", "tiff", "bmp"), container=g1)
		height <- list()
		height[['png']] <- height[['jpeg']] <- height[['tiff']] <- height[['bmp']] <- 500
		height[['pdf']] <- height[['postscript']] <- 7
		e$height <- height
		e$width <- 'default'
		gb <- bDem.gbutton('Save', container=g1)
		addHandlerClicked(gb, handler=saveGraph, action=list(mw=win, env=e, dpi=dpi, dev=dev.cur()))
		addHandlerChanged(proj.year, handler=newMap, action=list(dev=dev.cur(), map.pars=params))
		addHandlerDestroy(win, handler=close.map, action=list(dev=dev.cur()))
	} else {
		addHandlerChanged(proj.year, handler=newMap, action=list(map.pars=params))
	}
	return(g)
}
	

showDLcurve <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	country.pars <- get.country.code.from.widget(e$dlc.country$country.w, e$dlc.country)
	if(is.null(country.pars)) return(NULL)
	param.names.all <- list(text='sim.dir', numvector=c('pi'),
							numeric=c('nr.curves', 'burnin', 'tfr.max'),
							logical='predictive.distr')
	param.env <- get.parameters(param.names.all, env=e, quote=h$action$script)
	param.env.rest <- list(country=country.pars$code, output.dir=country.pars$output.dir,
							output.type=country.pars$output.type, verbose=TRUE)
	param.env <- c(param.env, 
					get.parameters(list(text=c('output.dir', 'output.type'), 
										logical='verbose', numeric='country'), 
									param.env.rest, quote=TRUE,
									retrieve.from.widgets=FALSE))

	param.mcmc <- param.env['sim.dir']
	if(h$action$script) {
		cmd <- paste('m <- get.tfr.mcmc(', assemble.arguments(param.mcmc), ')\n', sep='')
	} else {
		m <- do.call('get.tfr.mcmc', param.mcmc)
		cmd <- ''
	}
	pars.value <- svalue(e$graph.pars)
	param.plot1c <- list()
	for (par in c('pi', 'nr.curves', 'tfr.max', 'country', 'burnin', 'predictive.distr')) 
		if(is.element(par, names(param.env))) param.plot1c <- c(param.plot1c, param.env[par])

	if(!is.null(country.pars$code)) { # one country
		cmd <- paste(cmd, 'DLcurve.plot(mcmc.list=m, ', assemble.arguments(param.plot1c, pars.value), ')', sep='')
		if (h$action$script) {
			create.script.widget(cmd, h$action$mw, package="bayesTFR")
		} else {
			create.graphics.window(parent=h$action$mw, title=paste("Double Logistic Curves for", country.pars$name))
			eval(parse(text=cmd))
		}
	} else { # all countries
		param.plot.allc <- param.env[c(names(param.plot1c), 'output.dir', 'output.type',  'verbose')]
		cmd <- paste(cmd, 'DLcurve.plot.all(mcmc.list=m, ', assemble.arguments(param.plot.allc, pars.value), ')', sep='')
		if (h$action$script) {
			create.script.widget(cmd, h$action$mw, package="bayesTFR")
		} else {
			eval(parse(text=cmd))
		}
	}
}

showParTraces <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	param.names <- list(text='sim.dir', numeric=c('nr.points', 'burnin', 'thin'))
	params <- get.parameters(param.names, env=e, quote=FALSE)
	cs <- svalue(e$cs.chb, index=TRUE)
	all.pars <- svalue(e$pars.chb)
	print.summary <- h$action$print.summary
	if (cs==2) {
		country.pars <- get.country.code.from.widget(e$country$country.w, e$country)
		if(is.null(country.pars)) return(NULL)
	}
	type <- e$pred.type
	if(print.summary) {
		mc.summary <- c()
		warn <- getOption('warn')
		options(warn=-1) # disable warning messages
		mcmc.set <- do.call(paste('get.', type, '.mcmc', sep=''), params['sim.dir'])
		options(warn=warn)
		con <- textConnection("mc.summary", "w", local=TRUE)
		mc.exist <- TRUE
		sink(con)
		if (is.null(mcmc.set)) {
			cat('No simulation available in this directory.')
			mc.exist <- FALSE
		}
	} else create.graphics.window(parent=h$action$mw, title="Parameter traces", dpi=100)
	if (cs==2) { # country-specific parameters
		if (!all.pars) {
			pars <- svalue(e$par.cs.dl)
			if(print.summary) {if (mc.exist) print(summary(mcmc.set, country=country.pars$code, par.names.cs=pars, par.names=NULL, 
											burnin=params[['burnin']], thin=params[['thin']]))
			} else 
			do.call(paste(type, '.partraces.cs.plot', sep=''), c(list(country=country.pars$code, par.names=pars), params))
		} else {
			if(print.summary){if (mc.exist) print(summary(mcmc.set, country=country.pars$code, par.names=NULL, 
											burnin=params[['burnin']], thin=params[['thin']]))
			} else 
			do.call(paste(type, '.partraces.cs.plot', sep=''), c(list(country=country.pars$code), params))
		}
	} else { # World-parameters
		if (!all.pars) { # selected pars
			pars <- svalue(e$par.dl)
			if(print.summary) {if (mc.exist) print(summary(mcmc.set, par.names.cs=NULL, par.names=pars, 
											burnin=params[['burnin']], thin=params[['thin']]))
			} else 
			do.call(paste(type, '.partraces.plot', sep=''), c(list(par.names=pars), params))
		} else { # all pars
			if(print.summary) {if (mc.exist) print(summary(mcmc.set, par.names.cs=NULL, 
											burnin=params[['burnin']], thin=params[['thin']]))
			} else 
			do.call(paste(type, '.partraces.plot', sep=''), params)
		}
	}
	if(print.summary) {
		sink()
		close(con)
		sum.win <- gwindow('MCMC summary', parent=h$action$mw, width=500, height=400)
		set.widget.bgcolor(sum.win, "white")
		gtext(mc.summary, container=sum.win)
	}
}

computeConvergenceDiag <- function(h, ...) {
	e <- h$action$env
	type <- h$action$type
	if(!has.required.arguments(list(sim.dir='Simulation directory', burnin='Burnin'), env=e)) return()
	param.names <- list(numeric=c('burnin', 'thin', 'country.sampling.prop'),
						text=c('sim.dir'),
						logical=c('express', 'verbose', 'keep.thin.mcmc'))
	params <- get.parameters(param.names, e, quote=h$action$script)
	if(params$express || params$country.sampling.prop >= 1) params$country.sampling.prop <- NULL
	if (h$action$script) {
		cmd <- paste(type, '.diagnose(', assemble.arguments(c(params, e$params)), ')',sep='')
		create.script.widget(cmd, h$action$mw, package=h$action$package)
	} else {
		run <- FALSE
		mcmc.set <- do.call(paste('get.', type, '.mcmc', sep=''), list(sim.dir=params$sim.dir))
		iter <- get.total.iterations(mcmc.set$mcmc.list, burnin=params$burnin)
		if (iter < 0) gmessage('Number of iterations is smaller than burnin. Change the value of burnin.',
							container=h$action$mw)
		else {
			if (iter > 10000 && !params$express) {
				gconfirm('Computing convergence diagnostics with these settings can take a very long time. Do you want to continue?',
					icon='question', parent=h$action$mw,
					handler=function(h1, ...) run <<- TRUE)
			} else run <- TRUE
			if(run) 
				.run.diagnostics(e, type=h$action$package, handler=get.diagnostics.status, 
								option=paste('bDem', h$action$package, 'diagnose', sep='.'), 
								call=paste(type, 'diagnose', sep='.'), params=params, 
								sim.name=paste(h$action$package, 'diagnose'), main.win=h$action$mw,
								action=list(sb=e$statuslabel, package=h$action$package),
								interval=1000)
		}
	}
}

showConvergenceDiag <- function(h, ...) {
	e <- h$action$env
	type <- h$action$type
	dir <- svalue(e$sim.dir)
	diag.all <- do.call(paste('get.', type, '.convergence.all', sep=''), list(dir))
	ldiag <- length(diag.all)
	if(ldiag <=0) {
		gmessage(paste('There is no available convergence diagnostics in', dir), container=h$action$mw)
		return()
	}
	path = system.file("images",package="bayesDem")
	win <- bDem.gwindow('Available Convergence Diagnostics', parent=h$action$mw)
	g <- ggroup(horizontal=FALSE, container=win)
	verb <- 'are'
	noun.postfix <- 's'
	if(ldiag == 1) {
		verb <- 'is'
		noun.postfix <- ''
	}
	glabel(paste('There ', verb, ' ', ldiag, ' set', noun.postfix, ' of diagnostics.', sep=''),
			container=g)
	glabel('Click on the traffic lights to get a report.', container=g)
	addSpace(g, 10)
	glo <- glayout(container=g)
	glo[4,1] <- glabel('Needs at least', container=glo)
	l <- 1
	for(i in 1:ldiag) {
		diag <- diag.all[[i]]
		light <- names(diag$status)[diag$status]
		if(length(light) > 1) light <- 'green-yellow'
		image.name <- file.path(path, 'traffic_light', paste(light, 'png', sep='.'))
		glo[l,i+1] <- glabel(paste('Burnin =', diag$burnin), container=glo)
		glo[l+1,i+1] <- glabel(paste('Thin =', diag$thin), container=glo)
		glo[l+2,i+1] <- gimage(image.name, container=glo, handler=showDiagInfo, 
								action=list(mw=win, env=e, diag=diag))
		glo[l+3,i+1] <- glabel(diag$iter.needed, container=glo)
	}
	glo[4,ldiag+2] <- glabel('additional iterations', container=glo)
	addSpace(g,20)
}

showDiagInfo <- function(h, ...) {
	diag <- h$action$diag
	conv.diag <- c()
	con <- textConnection("conv.diag", "w", local=TRUE)
	sink(con)
	summary(diag, expand=TRUE)
	sink()
	close(con)
	win <- gwindow(paste('Convergence Diagnostics for burnin=', diag$burnin, sep=''), 
						parent=h$action$mw, width=500, height=400)
	gtext(conv.diag, container=win)
}

.run.diagnostics <- function(type='bayesTFR', ...) {
	statusopt <- paste('bDem.', type, '.diagnose.status', sep='')
	opt <- list()
	opt[statusopt] <- list(NULL)
	options(opt)
	res <- .run.simulation(...)
	options(opt)
	return(res)
}

get.diagnostics.status <- function(h, ...) 
	.update.status(h$action$sb, paste('bDem', h$action$package, 'diagnose.status', sep='.'), 'Running diagnostics ...')
