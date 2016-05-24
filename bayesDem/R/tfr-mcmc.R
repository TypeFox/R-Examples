TFRrunMCMCgroup <- function(g, main.win, parent) {
	# parent needed for wpp set in the parent function
	nb <- gnotebook(container=g, expand=TRUE)
	set.widget.bgcolor(nb, color.main)
	set.widget.basecolor(nb, color.nb.inactive)
	all.c.g <- ggroup(label="  <span color='#0B6138'>All Countries</span>  ", markup=TRUE, 
						horizontal=FALSE, container=nb)
	mcmc.all.countries.group(all.c.g, main.win, parent)
	extra.c.g <- ggroup(label="  <span color='#0B6138'>Extra Areas &amp; Regions</span>  ", 
						markup=TRUE, horizontal=FALSE, container=nb)
	mcmc.extra.countries.group(extra.c.g, main.win, parent)
	svalue(nb) <- 1
}

.enable.auto.run <- function(enable, e) {
	enabled(e$auto.conf.b) <- enable
	enabled(e$nr.chains) <- !enable
	enabled(e$iter) <- !enable 
}

.create.process.group <- function(g, e, defaults, show.buffer.size=TRUE, label='Process control') {
	leftcenter <- c(-1,0)
	paral.g <- gframe(paste("<span color='blue'>", label, "</span>"), markup=TRUE, horizontal=TRUE, spacing=10, container=g)
	pclo <- glayout(container=paral.g)
	pclo[1,1] <- e$verbose <- gcheckbox("Verbose", checked=defaults$verbose)
	pclo[1,2, anchor=leftcenter] <- "Verbose iter:"
	pclo[1,3] <- e$verbose.iter <- gedit(defaults$verbose.iter, width=4, container=pclo)
	tooltip(e$verbose.iter) <- 'Frequency of outputting the process status.'
	pclo[2,1] <- e$parallel <- gcheckbox("Parallel", checked=defaults$parallel, container=pclo)
	tooltip(e$parallel) <- 'Chains will be processed in parallel.'
	pclo[2,2, anchor=leftcenter] <- "Number of nodes:"
	pclo[2,3] <- e$nr.nodes <- gedit(defaults$nr.chains, width=2, container=pclo)
	tooltip(e$nr.nodes) <- 'Used if Parallel is checked. Set it to number of chains.'
	if(show.buffer.size) {
		pclo[3,2, anchor=leftcenter] <- "Buffer size:"
		pclo[3,3] <- e$buffer.size <- gedit(defaults$buffer.size, width=4, container=pclo)
		tooltip(e$buffer.size) <- "Number of iterations from which data are kept in memory before writing to disk."
	}
}

.create.mcmc.process.group <- function(g, e, main.win, defaults, type='tfr', advance.settigs.function=mcmc.advance.settings,
											mcmc.label='MCMC', process.label='Process control') {
	leftcenter <- c(-1,0)
	g2 <- ggroup(horizontal=TRUE, container=g)
	mcmc.g <- gframe(paste("<span color='blue'>", mcmc.label, "</span>"), markup=TRUE, horizontal=FALSE, spacing=10, container=g2)
	mclo <- glayout(container=mcmc.g)
	e$run.prediction <- FALSE
	mclo[1,1] <- e$run.auto <- gcheckbox("Auto simulation", checked=defaults$iter=='auto', container=mclo,
							handler=function(h,...){.enable.auto.run(svalue(h$obj), e)})
	mclo[1,2:4] <- e$auto.conf.b <- bDem.gbutton(' Configure auto run ', container=mclo, handler=configure.auto.run, 
				action=list(mw=main.win, env=e, cont.run=FALSE, type=type))
	mclo[2,1, anchor=leftcenter] <- glabel("Number of chains:")
	mclo[2,2] <- e$nr.chains <- gedit(defaults$nr.chains, width=2)
	mclo[3,1, anchor=leftcenter] <- glabel("Number of iterations:")
	mclo[3,2] <- e$iter <- gedit(defaults$iter, width=7)
	tooltip(e$iter) <- "How many iterations per chain to be simulated."
	mclo[2,3, anchor=leftcenter] <- glabel("Thin:")
	mclo[2,4] <- e$thin <- gedit(defaults$thin, width=2)
	tooltip(e$thin) <- "Thinning for storing MCMCs on disk."
	.enable.auto.run(defaults$iter=='auto', e)
	mclo[3,3, anchor=leftcenter] <- glabel("RNG seed:")
	mclo[3,4] <- e$seed <- gedit(defaults$seed, width=4)
	mclo[4,1, anchor=leftcenter] <- glabel("File compression:")
	mclo[4,2] <- e$compression.type <- bDem.gdroplist(c('None', 'xz', 'bz', 'gz'), container=mclo)
	mclo[5,1:4] <- bDem.gbutton('Priors & Advanced MCMC Settings',  handler=advance.settigs.function, 
				action=list(mw=main.win, env=e))	
	addSpring(g2)
	.create.process.group(g2, e, defaults, label=process.label)
}
mcmc.all.countries.group <- function(g, main.win, parent) {
	e <- new.env()
	defaults <- formals(run.tfr.mcmc) # default argument values
	defaults3 <- formals(run.tfr3.mcmc)
	e$output.dir <- parent$sim.dir
	leftcenter <- c(-1,0)	
	addSpace(g, 10)
	time.g <- gframe("<span color='blue'>TFR time series</span>", markup=TRUE, horizontal=FALSE, spacing=10, container=g)
	timelo <- glayout(container=time.g)
	timelo[1,1, anchor=leftcenter] <- glabel("Start year:", container=timelo)
	timelo[1,2] <- e$start.year <- gedit(defaults$start.year, width=4, container=timelo)
	tooltip(e$start.year) <- 'Historical data prior to this year will be ignored.'
	timelo[1,3, anchor=leftcenter] <- glabel("Present year:", container=timelo)
	timelo[1,4] <- e$present.year <- gedit(defaults$present.year, width=4, container=timelo)
	tooltip(e$present.year) <- 'Historical data after this year will be ignored.'
	timelo[1,5, anchor=leftcenter] <- glabel("WPP year:", container=timelo)
	timelo[1,6, anchor=c(1,0)] <- wpp <- glabel(parent$wpp.year, container=timelo)
	tooltip(wpp) <- 'To change this start bayesDem with wpp.year as argument.'
	timelo[2,1:2, anchor=leftcenter] <- glabel("User-defined TFR file:", container=timelo)
	timelo[2,3:6] <- e$my.tfr.file <- bDem.gfilebrowse(eval(defaults$my.tfr.file), type='open', 
					  width=30, quote=FALSE, container=timelo)
	tooltip(e$my.tfr.file) <- 'Overwrites default wpp data.'
					  
	addSpace(g, 10)
	phase.g <- gframe("", markup=TRUE, horizontal=FALSE, spacing=10, container=g)
	phase.g1 <- ggroup(horizontal=TRUE, container=phase.g)
	glabel("Run TFR MCMC for ", container=phase.g1)
	e$mcmc.type <- gcheckboxgroup(c('Phase II', 'Phase III'), horizontal=TRUE, checked=TRUE, container=phase.g1,
								handler=function(h,...){
									enabled(phaseII.g) <- is.element(1, svalue(h$obj, index=TRUE))
									enabled(phaseIII.g) <- is.element(2, svalue(h$obj, index=TRUE))
								})
	type.nb <- gnotebook(container=phase.g, expand=TRUE)
	set.widget.bgcolor(type.nb, color.main)
	set.widget.basecolor(type.nb, color.nb.inactive)
	phaseII.g <- ggroup(label="<span color='darkred'>Phase II</span>", markup=TRUE, horizontal=FALSE, container=type.nb)
	addSpace(phaseII.g, 10)
	.create.mcmc.process.group(phaseII.g, e, main.win, defaults, mcmc.label="MCMC Phase II", process.label="Process control for Phase II")
	#addSpace(g, 10)
	phaseIII.g <- ggroup(label="<span color='darkred'>Phase III</span>", markup=TRUE, horizontal=FALSE, container=type.nb)
	addSpace(phaseIII.g, 10)
	e$phase3 <- new.env()
	.create.mcmc.process.group(phaseIII.g, e$phase3, main.win, defaults3, type='tfr3', mcmc.label="MCMC Phase III", 
								advance.settigs.function=mcmc3.advance.settings,
								process.label="Process control for Phase III")
	svalue(type.nb) <- 1
	addSpace(g, 10)
	.create.status.label(g, e)
	addSpring(g)
	adv.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic=c('run.tfr.mcmc', 'run.tfr3.mcmc', 'bayesTFR-package'), package='bayesTFR', parent.group=adv.g,
						parent.window=main.win)
				
	addSpring(adv.g)
	create.generate.script.button(handler=mcmc.run, action=list(mw=main.win, env=e, script=TRUE, wpp.year=parent$wpp.year),
								container=adv.g)
	bDem.gbutton(action=gaction(label=' Run MCMC ', icon='execute', handler=mcmc.run, 
				action=list(mw=main.win, env=e, script=FALSE, wpp.year=parent$wpp.year, parent.group=g)), 
				container=adv.g)
	return(e)
	}


mcmc.run <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(output.dir='Simulation directory'), env=e)) return()
	param.names <- list(numeric=c('buffer.size', 'nr.nodes', 'iter', 'thin', 'nr.chains', 'start.year', 
									'present.year', 'seed', 'verbose.iter'),
						text=c('output.dir', 'my.tfr.file', 'compression.type'),
						logical=c('verbose', 'parallel'))
	param.names.p3 <- list(numeric=c('buffer.size', 'nr.nodes', 'iter', 'thin', 'nr.chains', 'seed', 'verbose.iter'),
						text=c('compression.type'),
						logical=c('verbose', 'parallel'))
	params <- get.parameters(param.names, e, quote=h$action$script)
	params.p3 <- get.parameters(param.names.p3, e$phase3, quote=h$action$script)
	params[['wpp.year']] <- h$action$wpp.year
	params.p3[['my.tfr.file']] <- params[['my.tfr.file']]
	
	mcmc.type <- svalue(e$mcmc.type, index=TRUE)
	if(length(mcmc.type)<=0) {
		gmessage('Eiter Phase II or Phase III must be checked.', title='Input Error', icon='error')
		return()
	}
	envs <- list(e, e$phase3)
	parlist <- list(params, params.p3)
	auto <- list()
	for(i in 1:2) {
		auto[[i]] <- list()
		auto[[i]]$run.auto <- FALSE
		if(!is.element(i, mcmc.type)) next
		auto[[i]]$run.auto <- svalue(envs[[i]]$run.auto)
		if (auto[[i]]$run.auto) {
			parlist[[i]][['auto.conf']] <- envs[[i]]$auto.conf
			op <- options("useFancyQuotes")
			options(useFancyQuotes = FALSE)
			parlist[[i]][['iter']] <- if (h$action$script) sQuote('auto') else 'auto'
			options(op)
		}
	}
	outdir <- get.parameters(list(text='output.dir'), e, quote=FALSE)$output.dir # to avoid double quoting if script is TRUE
	if(!is.element(1, mcmc.type)) {
		m <- get.tfr.mcmc(outdir)
		if(is.null(m)) {
			gmessage('Simulation directory contains no valid MCMC results. Phase II must be run before phase III.', 
							title='Input Error', icon='error')
			return()
		}
	}
	if(file.exists(outdir) && (length(list.files(outdir)) > 0) && is.element(1, mcmc.type)) {
		params[['replace.output']] <- FALSE
		if (gconfirm(paste('Non-empty directory', outdir, 
								'already exists.\nDo you want to overwrite existing results?'),
				icon='question', parent=h$action$mw))
			params[['replace.output']] <- TRUE
		else return(NULL)
	}
	run.prediction <- (auto[[1]]$run.auto || auto[[2]]$run.auto) && (e$run.prediction || e$phase3$run.prediction)
	use.diag <- if(auto[[1]]$run.auto && e$run.prediction) 'TRUE' else 'FALSE'
	if (h$action$script) {
		commands <- ''
		if(is.element(1, mcmc.type)) 
			commands <- paste(commands, 'm <- run.tfr.mcmc(', assemble.arguments(c(parlist[[1]], e$params)), ')', sep=' ')
		if(is.element(2, mcmc.type)) 
			commands <- paste(commands, '\n\nm3 <- run.tfr3.mcmc(', parlist[[1]]$output.dir, ', ', 
									assemble.arguments(c(parlist[[2]], e$phase3$params)), ')', sep=' ')
		if(run.prediction) 
			commands <- paste(commands, '\n\ntfr.predict(sim.dir=,', parlist[[1]]$output.dir, 
							', use.diagnostics=',use.diag, ', replace.output=TRUE)', sep='')
		create.script.widget(commands, h$action$mw, package="bayesTFR")
	} else { # run a simulation
		run <- FALSE
		long.run <- FALSE
		for (i in 1:2) {
			if (is.element(i, mcmc.type) && ((parlist[[i]][['iter']] == 'auto' && ((!is.null(parlist[[i]][['auto.conf']]) 
				&& parlist[[i]][['auto.conf']]$iter > 100) || is.null(parlist[[i]][['auto.conf']]))) 
					|| (parlist[[i]][['iter']] != 'auto' && parlist[[i]][['iter']] > 100))) {
						long.run <- TRUE
						break
			}
		}
		if(long.run) {
			gconfirm('Running MCMC with these settings can take a very long time. Do you want to continue?',
					icon='question', parent=h$action$mw,
					handler=function(h, ...) run <<- TRUE)
		} else run <- TRUE
		if(run) {
			if(is.element(1, mcmc.type)) { # run phase II MCMCs
				if(file.exists(params[['output.dir']])) unlink(params[['output.dir']], recursive=TRUE)
				m <- .run.simulation(e, handler=get.tfr.simulation.status, option='bDem.TFRmcmc', 
								call='run.tfr.mcmc', params=parlist[[1]], 
								sim.name='TFR MCMC Phase II', main.win=h$action$mw,
								action=list(sb=e$statuslabel, sim.dir=params[['output.dir']]),
								interval=5000)
			} 
			if(is.element(2, mcmc.type)) {
				m3 <- .run.simulation(e, handler=get.tfr.simulation.status, option='bDem.TFRmcmc', 
								call='run.tfr3.mcmc', params=c(list(sim.dir=params[['output.dir']]), parlist[[2]]), 
								sim.name='TFR MCMC Phase III', main.win=h$action$mw,
								action=list(sb=e$statuslabel, sim.dir=params[['output.dir']]),
								interval=5000)
			}
			if(run.prediction) {
				.run.prediction(e, handler=get.tfr.prediction.status, option='bDem.TFRpred', 
								call='tfr.predict', params=list(sim.dir=params[['output.dir']], use.diagnostics=eval(use.diag), replace.output=TRUE), 
								sim.name='TFR prediction', main.win=h$action$mw,
								action=list(sb=e$statuslabel),
								interval=1000)
			}
		}
	}
}

.run.simulation <- function(e, handler, option, call, params, sim.name, main.win, action=list(), interval=1000) {
	svalue(e$statuslabel) <- paste('Starting', sim.name, '...')
	handler.id <- addHandlerIdle(e$statuslabel, interval=interval, handler = handler, action=action)
	opt <- list()
	opt[[option]] <- TRUE
	options(opt)
	m <- try(do.call(call, params))
	opt[[option]] <- FALSE
	options(opt)
	if(inherits(m, "try-error")) svalue(e$statuslabel) <- paste(sim.name, 'failed.')
	else svalue(e$statuslabel) <- paste(sim.name, 'finished.')
	removeHandler(main.win, handler.id)
	gSourceRemove(handler.id)
	return(m)
}

get.tfr.simulation.status <- function(h, ...) {
	sb <- h$action$sb
	sim.dir <- h$action$sim.dir
	warn <- getOption('warn')
	options(warn=-1) # disable warning messages
	mcmc.set <- get.tfr.mcmc(sim.dir)
	options(warn=warn)
	if(is.null(mcmc.set)) return()
	get.finished <- function(x) return(x$finished.iter)
	get.iter <- function(x) return(x$iter)
	finished <- sapply(mcmc.set$mcmc.list, get.finished)
	total <- sapply(mcmc.set$mcmc.list, get.iter)
	svalue(sb) <- paste('Running TFR simulation ... Iterations - ', 
				 	paste('chain ', 1:length(finished), ': ', finished, ' (', total, ')', sep='', collapse=', '))
	
}

get.tfr.simulation.extra.status <- function(h, ...) {
	sb <- h$action$sb
	# We don't have any status info for the extra run, so just keep a static label
	svalue(sb) <- 'Running TFR extra simulation ... '
	
}


configure.auto.run <- function(h, ...) {
	cont.run <- h$action$cont.run
	type <- if(is.null(h$action$type)) 'tfr' else h$action$type
	if (!cont.run)
		defaults <- eval(formals(paste("run.", type, ".mcmc", sep=''))$auto.conf)
	else {
		defaults <- do.call(paste('.get.defaults.for.auto.cont.', type, sep=''), list(h$action$env))
		if(is.null(defaults)) return()
	}
	set.defaults <- function(h2, ...) {
		for (par in names(defaults)) {
			svalue(h$action$env$auto.conf.env[[par]]) <- defaults[[par]]
		}
		svalue(h$action$env$auto.conf.env[['run.prediction']]) <- FALSE
	}
	set.auto.conf <- function(h2, ...) {
		params <- get.parameters(list(numvector=names(defaults)), env=h$action$env$auto.conf.env)
		h$action$env$auto.conf <- params
		h$action$env$run.prediction <- svalue(h$action$env$auto.conf.env$run.prediction)
		visible(h$action$env$auto.conf.win) <- FALSE
	}
	
	if (!is.null(h$action$env$auto.conf.win) && !h$action$env$auto.conf.env$window.destroyed) { # window exists
		if(!is.null(h$action$env$auto.conf)) { # OK button previously clicked 
			for (par in names(defaults)) 
				svalue(h$action$env$auto.conf.env[[par]]) <- h$action$env$auto.conf[[par]]
			svalue(h$action$env$auto.conf.env$run.prediction) <- h$action$env$run.prediction
		} else { # OK button not clicked yet, values are set to defaults
			for (par in names(defaults)) 
				svalue(h$action$env$auto.conf.env[[par]]) <- defaults[[par]]
				#svalue(h$action$env$auto.conf.env[[par]]) <- h$action$env$auto.conf.env$defaults[[par]]
			svalue(h$action$env$auto.conf.env[['run.prediction']]) <- FALSE
		}
		visible(h$action$env$auto.conf.win) <- TRUE
	} else { # create the window
		
	h$action$env$auto.conf.win <- auto.conf.win <- 
					bDem.gwindow('Configuration of Auto Run',
						parent=h$action$mw, visible=FALSE,
						handler=function(h, ...) {
							h$action$env$auto.run.okhandler <- NULL
						})
	e <- new.env()
	g <- ggroup(container=auto.conf.win, horizontal=FALSE)
	mcmc.g <- gframe("<span color='blue'>MCMC</span>", markup=TRUE, horizontal=FALSE, container=g)
	#mcmc.g1 <- ggroup(container=mcmc.g, horizontal=TRUE)
	mclo <- glayout(container=mcmc.g)
	mclo[1,1] <- glabel("Number of chains:", container=mclo)
	mclo[1,2] <- e$nr.chains <- gedit(defaults$nr.chains, width=2, container=mclo)
	mclo[2,1] <- glabel("Number of iterations:", container=mclo)
	mclo[2,2] <- e$iter <- gedit(defaults$iter, width=7, container=mclo)
	tooltip(e$iter) <- "Number of iterations to run inititally after which the covergence is checked."
	enabled(e$nr.chains) <- !cont.run

	#mcmc.g2 <- ggroup(container=mcmc.g, horizontal=TRUE)
	mclo[3,1] <- glabel("Iteration increments: ", container=mclo)
	mclo[3,2] <- e$iter.incr <- gedit(defaults$iter.incr, width=7, container=mclo)
	tooltip(e$iter.incr) <- "Number of iterations to run in the following loops after each of which the convergence is checked."
	
	#conv.g <- gframe("<span color='blue'>Convergence diagnostics</span>", markup=TRUE, horizontal=TRUE, container=g)
	mclo[1,3:4] <- 'Convergence diagnostics:'
	mclo[2,3] <- glabel("Burnin:", container=mclo)
	mclo[2,4] <- e$burnin <- gedit(defaults$burnin, width=7, container=mclo)
	tooltip(e$burnin) <- "Burnin used when checking convergence."
	mclo[3,3] <- glabel("Thin:", container=mclo)
	mclo[3,4] <- e$thin <- gedit(defaults$thin, width=7, container=mclo)
	tooltip(e$thin) <- "Thin used when checking convergence."
	
	setup.g <- gframe("<span color='blue'>Run setup</span>", markup=TRUE, horizontal=TRUE, container=g)
	glabel("Maximum loops:", container=setup.g)
	e$max.loops <- gedit(defaults$max.loops, width=2, container=setup.g)
	tooltip(e$max.loops) <- "How many times should the convergence be checked?"
	addSpace(setup.g, 10)
	e$run.prediction <- gcheckbox("Make predictions", checked=FALSE, container=setup.g)
	tooltip(e$run.prediction) <- "If checked the prediction is generated using settings from the converged MCMCs."
	addSpring(g)
	# Buttons
	button.g <- ggroup(container=g, horizontal=TRUE)
	bDem.gbutton('Cancel', container=button.g, handler=function(h, ...) 
					visible(auto.conf.win) <- FALSE)
	addSpring(button.g)
	bDem.gbutton(action=gaction(label='  Set to Default Values  ', icon='refresh', handler=set.defaults), container=button.g)
	e$auto.conf.okbutton <- bDem.gbutton('OK', container=button.g)
	
		if(!is.null(h$action$env$auto.conf)) { # OK button previously clicked 
			for (par in names(defaults)) 
				svalue(e[[par]]) <- h$action$env$auto.conf[[par]]
			svalue(e$run.prediction) <- h$action$env$run.prediction
		}
			
	visible(auto.conf.win) <- TRUE
	e$window.destroyed <- FALSE
	h$action$env$auto.conf.env <- e
	addHandlerDestroy(auto.conf.win, 
							handler=function(h1, ...) h$action$env$auto.conf.env$window.destroyed <- TRUE)

	}
	if(!is.null(h$action$env$auto.conf.okhandler)) 
		removehandler(h$action$env$auto.conf.env$auto.conf.okbutton, h$action$env$auto.conf.okhandler)
	h$action$env$auto.conf.okhandler <- addhandlerclicked(h$action$env$auto.conf.env$auto.conf.okbutton, 
												handler=set.auto.conf)

}

mcmc.advance.settings <- function(h, ...) {
	param.names <- list('Triangle_c4.low', 'Triangle_c4.up', 'Triangle_c4.trans.width', 'Triangle4.0',
						'mean.eps.tau0', 'sd.eps.tau0', 'delta4.0',
						'delta0', 'nu.delta0', 'd.low', 'd.up', 'd.trans.width',
						'dl.p1', 'dl.p2', 'U.c.low', 'U.up', 'U.width', 'S.low', 'S.up', 'S.width', 
						'a.low', 'a.up', 'a.width', 'b.low', 'b.up', 'b.width', 
						'const.low', 'const.up', 'const.width', 
						'sigma0.low', 'sigma0.up', 'sigma0.width', 'sigma0.min',
						'chi0', 'psi0', 'nu.psi0', 'alpha0.p', 'nu4', 'nu.tau0',
						'S.ini', 'a.ini', 'b.ini', 'const.ini', 'gamma.ini', 'sigma0.ini', 'Triangle_c4.ini'
						)
	get.defaults <- function() {
		all.defaults <- formals(run.tfr.mcmc) # default argument values
		defaults <- list()
		for (par in param.names) { 
			defaults[[par]] <- all.defaults[[par]]
			if (length(defaults[[par]]) > 1) {
				if (defaults[[par]][[1]] == '-')  # handle negative values
					defaults[[par]] <- paste(defaults[[par]], sep='', collapse='')
			}
		}
		# special cases
		defaults$b.low <- defaults$a.low
		defaults$b.up <- defaults$a.up
		defaults$alpha0.p <- '-1, 0.5, 1.5'
		return(defaults)
	}

	set.defaults <- function(h2, ...) {
		defaults <- get.defaults()
		for (par in param.names) {
			svalue(h$action$env$adv.set.env[[par]]) <- defaults[[par]]
		}
		for (par in names(linked.pars.list)) svalue(h$action$env$adv.set.env[[par]]) <- defaults[[linked.pars.list[[par]]]]
	}
	
	set.advance.pars <- function(h2, ...) {
		params <- get.parameters(list(numvector=param.names), env=h$action$env$adv.set.env)
		#params <- list()
		#for (par in param.names) {
		#	params[[par]] <- as.numeric(strsplit(svalue(h$action$env$adv.set.env[[par]]), ',')[[1]])
		#	}
		h$action$env$params <- params
		visible(h$action$env$adv.set.win) <- FALSE
	}
		
	if (!is.null(h$action$env$adv.set.win) && !h$action$env$adv.set.env$window.destroyed) { #Advanced Parameters window exists
		if(!is.null(h$action$env$params)) { # OK button previously clicked 
			for (par in param.names) {
				#print (par)
				#print(h$action$env$params[[par]])
				#print(svalue(h$action$env$adv.set.env[[par]]))
				svalue(h$action$env$adv.set.env[[par]]) <- paste(h$action$env$params[[par]], collapse=', ')
			}
		} else { # OK button not clicked yet, values are set to defaults
			for (par in param.names) {
				svalue(h$action$env$adv.set.env[[par]]) <- h$action$env$adv.set.env$defaults[[par]]
			}
		}
		visible(h$action$env$adv.set.win) <- TRUE
	} else { # create the Advanced Parameters window
		h$action$env$adv.set.win <- adv.set.win <- 
					bDem.gwindow('Settings for Bayesian Hierarchical TFR Model ',
						parent=h$action$mw, visible=FALSE,
						handler=function(h, ...) {
							h$action$env$adv.set.okhandler <- NULL
						})
		e <- new.env()
		e$defaults <- defaults <- get.defaults()
		e$adv.g <- ggroup(container=adv.set.win, horizontal=FALSE)
	
	linked.pars.list <- list()
	
	modelpars.f <- gframe("<span color='blue'>Model parameters</span>", markup=TRUE, container=e$adv.g, horizontal=TRUE)
	modelpars.flo <- glayout(container=modelpars.f)
	l <- 1 # row 1
	modelpars.flo[l,1] <- ''
	modelpars.flo[l,2] <- '   lower bound   '
	modelpars.flo[l,3] <- '   upper bound   '
	modelpars.flo[l,4] <- 'width (trans)'
		
	l <- l+1 # new row
	#modelpars.flo[l,1] <- glabel('<span>&#916;<sub>c4</sub>:</span>', markup=TRUE, container=modelpars.flo)
	modelpars.flo[l,1] <- glabel('<span>Triangle<sub>c4</sub>:</span>', markup=TRUE, container=modelpars.flo)
	modelpars.flo[l,2] <- e$Triangle_c4.low <- gedit(defaults$Triangle_c4.low, width=5, container=modelpars.flo)
	modelpars.flo[l,3] <- e$Triangle_c4.up <- gedit(defaults$Triangle_c4.up, width=5, container=modelpars.flo)
	modelpars.flo[l,4] <- e$Triangle_c4.trans.width <- gedit(defaults$Triangle_c4.trans.width, width=5, container=modelpars.flo)
	#modelpars.flo[l,5] <- '      '
	
	
	l <- l+1 # new row
	modelpars.flo[l,1] <- 'd:'
	modelpars.flo[l,2] <- e$d.low <- gedit(defaults$d.low, width=5, container=modelpars.flo)
	modelpars.flo[l,3] <- e$d.up <- gedit(defaults$d.up, width=5, container=modelpars.flo)
	modelpars.flo[l,4] <- e$d.trans.width <- gedit(defaults$d.trans.width, width=5, container=modelpars.flo)
	#modelpars.flo[l,5] <- '(trans)'
	
	addSpace(modelpars.f, 30)
	modelpars.g2 <- glayout(container=modelpars.f)
	l <- 1
	modelpars.g2[l,1] <- ''
	l <- l+1
	modelpars.g2[l,1] <- glabel('p1:', container=modelpars.g2)
	modelpars.g2[l,2] <- e$dl.p1 <- gedit(defaults$dl.p1, width=5, container=modelpars.g2)
	#addSpace(modelpars.g2, 10)
	l <- l+1
	modelpars.g2[l,1] <- glabel('p2:', container=modelpars.g2)
	modelpars.g2[l,2] <- e$dl.p2 <- gedit(defaults$dl.p2, width=5, container=modelpars.g2)
	
	# Priors
	priors.f <- gframe("<span color='blue'>Prior distributions</span>", markup=TRUE, container=e$adv.g, horizontal=TRUE)
	
	priors.uniform.f <- gframe("<span  color='#0B6138'>Uniform Priors</span>", markup=TRUE, container=priors.f)
	uniform.flo <- glayout(container=priors.uniform.f)
	
	l <- 1 # row 1
	uniform.flo[l,1] <- ''
	uniform.flo[l,2] <- '   lower   '
	uniform.flo[l,3] <- '   upper   '
	uniform.flo[l,4] <- ' width '
	uniform.flo[l,5] <- ' min '
	
	l <- l+1 # new row
	uniform.flo[l,1] <- 'U:'
	uniform.flo[l,2] <- e$U.c.low <- gedit(defaults$U.c.low, width=5, container=uniform.flo)
	uniform.flo[l,3] <- e$U.up <- gedit(defaults$U.up, width=5, container=uniform.flo)
	uniform.flo[l,4] <- e$U.width <- gedit(defaults$U.width, width=5, container=uniform.flo)

	l <- l+1 # new row
	uniform.flo[l,1] <- 'S:'
	uniform.flo[l,2] <- e$S.low <- gedit(defaults$S.low, width=5, container=uniform.flo)
	uniform.flo[l,3] <- e$S.up <- gedit(defaults$S.up, width=5, container=uniform.flo)
	uniform.flo[l,4] <- e$S.width <- gedit(defaults$S.width, width=5, container=uniform.flo)
	
	l <- l+1 # new row
	uniform.flo[l,1] <- 'a:'
	uniform.flo[l,2] <- e$a.low <- gedit(defaults$a.low, width=5, container=uniform.flo)
	uniform.flo[l,3] <- e$a.up <- gedit(defaults$a.up, width=5, container=uniform.flo)
	uniform.flo[l,4] <- e$a.width <- gedit(defaults$a.width, width=5, container=uniform.flo)

	l <- l+1 # new row
	uniform.flo[l,1] <- 'b:'
	uniform.flo[l,2] <- e$b.low <- gedit(defaults$b.low, width=5, container=uniform.flo)
	uniform.flo[l,3] <- e$b.up <- gedit(defaults$b.up, width=5, container=uniform.flo)
	uniform.flo[l,4] <- e$b.width <- gedit(defaults$b.width, width=5, container=uniform.flo)

	l <- l+1 # new row
	uniform.flo[l,1] <- glabel('<span>sigma<sub>0</sub>:</span>', markup=TRUE, container=uniform.flo)
	uniform.flo[l,2] <- e$sigma0.low <- gedit(defaults$sigma0.low, width=5, container=uniform.flo)
	uniform.flo[l,3] <- e$sigma0.up <- gedit(defaults$sigma0.up, width=5, container=uniform.flo)
	uniform.flo[l,4] <- e$sigma0.width <- gedit(defaults$sigma0.width, width=5, container=uniform.flo)
	uniform.flo[l,5] <- e$sigma0.min <- gedit(defaults$sigma0.min, width=5, container=uniform.flo)

	l <- l+1 # new row
	uniform.flo[l,1] <- 'const (c):'
	uniform.flo[l,2] <- e$const.low <- gedit(defaults$const.low, width=5, container=uniform.flo)
	uniform.flo[l,3] <- e$const.up <- gedit(defaults$const.up, width=5, container=uniform.flo)
	uniform.flo[l,4] <- e$const.width <- gedit(defaults$const.width, width=5, container=uniform.flo)
	
	normal.gamma.g <- ggroup(horizontal=TRUE, container=priors.f)
	priors.normal.f <- gframe("<span  color='#0B6138'>Normal Priors</span>", markup=TRUE, container=normal.gamma.g,
								horizontal=FALSE)
	normal.flo <- glayout(container=priors.normal.f)
	
	l <- 1 # row 1
	normal.flo[l,1] <- ''
	normal.flo[l,2] <- ' mean '
	normal.flo[l,3] <- ' sd '
		
	l <- l+1 # new row
	normal.flo[l,1] <- glabel('<span>chi:</span>', markup=TRUE, container=normal.flo)
	normal.flo[l,2] <- e$chi0 <- gedit(defaults$chi0, width=5, container=normal.flo)
	normal.flo[l,3] <- e$l.psi0 <- glabel(width=5, container=normal.flo)

	l <- l+1 # new row
	normal.flo[l,1] <- glabel('<span>alpha<sub>i</sub>:</span>', markup=TRUE, container=normal.flo)
	normal.flo[l,2] <- e$alpha0.p <- gedit(defaults$alpha0.p, width=8, container=normal.flo)
	normal.flo[l,3] <- e$l.delta0 <- glabel(width=5, container=normal.flo)
	
	l <- l+1 # new row
	normal.flo[l,1] <- glabel('<span>Triangle<sub>4</sub>:</span>', markup=TRUE, container=normal.flo)
	normal.flo[l,2] <- e$Triangle4.0 <- gedit(defaults$Triangle4.0, width=5, container=normal.flo)
	normal.flo[l,3] <- e$l.delta4.0 <- glabel(width=5, container=normal.flo)
	
	l <- l+1 # new row
	normal.flo[l,1] <- glabel('<span>m<sub>tau</sub>:</span>', markup=TRUE, container=normal.flo)
	normal.flo[l,2] <- e$mean.eps.tau0 <- gedit(defaults$mean.eps.tau0, width=5, container=normal.flo)
	normal.flo[l,3] <- e$l.sd.eps.tau0 <- glabel(width=5, container=normal.flo)

	addSpace(priors.normal.f, 10)
	normal.g2 <- ggroup(container=priors.normal.f)
	glabel('sd given by s_0 in Gamma prior', container=normal.g2)

	priors.gamma.f <- gframe("<span  color='#0B6138'>Gamma Priors</span>", markup=TRUE, container=priors.f,
								horizontal=FALSE)
	gamma.flo <- glayout(container=priors.gamma.f)
	
	l <- 1 # row 1
	gamma.flo[l,1] <- ''
	gamma.flo[l,2] <- ' nu_0 '
	gamma.flo[l,3] <- ' s_0  '

	l <- l+1 # new row
	gamma.flo[l,1] <- glabel('<span>1/psi<sup>2</sup>:</span>', markup=TRUE, container=gamma.flo)
	gamma.flo[l,2] <- e$nu.psi0 <- gedit(defaults$nu.psi0, width=5, container=gamma.flo)
	gamma.flo[l,3] <- e$psi0 <- gedit(defaults$psi0, width=5, container=normal.flo)
	svalue(e$l.psi0) <- svalue(e$psi0)
	addHandlerChanged(e$psi0, handler=function(h,...) svalue(e$l.psi0) <- svalue(e$psi0))
	linked.pars.list[['l.psi0']] <- 'psi0'
	
	l <- l+1 # new row
	gamma.flo[l,1] <- glabel('<span>1/delta<sub>i</sub><sup>2</sup>:</span>', markup=TRUE, container=gamma.flo)
	gamma.flo[l,2] <- e$nu.delta0 <- gedit(defaults$nu.delta0, width=5, container=gamma.flo)
	gamma.flo[l,3] <- e$delta0 <- gedit(defaults$delta0, width=5, container=gamma.flo)
	svalue(e$l.delta0) <- svalue(e$delta0)
	addHandlerChanged(e$delta0, handler=function(h,...) svalue(e$l.delta0) <- svalue(e$delta0))
	linked.pars.list[['l.delta0']] <- 'delta0'
	
	l <- l+1 # new row
	gamma.flo[l,1] <- glabel('<span>1/delta<sub>4</sub><sup>2</sup>:</span>', markup=TRUE, container=gamma.flo)
	gamma.flo[l,2] <- e$nu4 <- gedit(defaults$nu4, width=5, container=gamma.flo)
	gamma.flo[l,3] <- e$delta4.0 <- gedit(defaults$delta4.0, width=5, container=gamma.flo)
	svalue(e$l.delta4.0) <- svalue(e$delta4.0)
	addHandlerChanged(e$delta4.0, handler=function(h,...) svalue(e$l.delta4.0) <- svalue(e$delta4.0))
	linked.pars.list[['l.delta4.0']] <- 'delta4.0'
	
	l <- l+1 # new row
	gamma.flo[l,1] <- glabel('<span>1/s<sub>tau</sub><sup>2</sup>:</span>', markup=TRUE, container=gamma.flo)
	gamma.flo[l,2] <- e$nu.tau0 <- gedit(defaults$nu.tau0, width=5, container=gamma.flo)
	gamma.flo[l,3] <- e$sd.eps.tau0 <- gedit(defaults$sd.eps.tau0, width=5, container=gamma.flo)
	svalue(e$l.sd.eps.tau0) <- svalue(e$sd.eps.tau0)
	addHandlerChanged(e$sd.eps.tau0, handler=function(h,...) svalue(e$l.sd.eps.tau0) <- svalue(e$sd.eps.tau0))
	linked.pars.list[['l.sd.eps.tau0']] <- 'sd.eps.tau0'
	
	addSpace(priors.gamma.f, 10)
	gamma.g2 <- glayout(container=priors.gamma.f, horizontal=FALSE)
	gamma.g2[1,1, expand=FALSE, anchor=c(-1,0)] <- glabel('Gamma shape: nu_0/2', container=gamma.g2)
	gamma.g2[2,1, expand=FALSE] <- glabel('Gamma rate: nu_0/2*s_0^2', container=gamma.g2)
	
	# Starting values
	start.f <- gframe("<span color='blue'>Starting values</span>", markup=TRUE, container=e$adv.g, horizontal=FALSE)
	start.g1 <- ggroup(container=start.f, horizontal=TRUE)
	glabel('Leave empty for starting values being equally-spaced between lower and upper bound.', container=start.g1)
	
	start.flo <- ggroup(container=start.f, horizontal=TRUE)
	glabel('S:', container=start.flo)
	e$S.ini <- gedit(defaults$S.ini, width=5, container=start.flo)
	glabel('a:', container=start.flo)
	e$a.ini <- gedit(defaults$a.ini, width=5, container=start.flo)
	glabel('b:', container=start.flo)
	e$b.ini <- gedit(defaults$b.ini, width=5, container=start.flo)
	glabel('<span>sigma<sub>0</sub>:</span>', markup=TRUE, container=start.flo)
	e$sigma0.ini <- gedit(defaults$sigma0.ini, width=5, container=start.flo)
	glabel('const (c):', container=start.flo)
	e$const.ini <- gedit(defaults$const.ini, width=5, container=start.flo)
	glabel('<span>Triangle<sub>c4</sub>:</span>', markup=TRUE, container=start.flo)
	e$Triangle_c4.ini <- gedit(defaults$Triangle_c4.ini, width=5, container=start.flo)
	addSpace(start.flo, 20)
	glabel('gamma:', container=start.flo)
	e$gamma.ini <- gedit(defaults$gamma.ini, width=5, container=start.flo)

	# Buttons
	button.g <- ggroup(container=e$adv.g, horizontal=TRUE)
	bDem.gbutton('Cancel', container=button.g, handler=function(h, ...) 
					visible(adv.set.win) <- FALSE)
	addSpring(button.g)
	e$adv.set.defaultbutton <- bDem.gbutton(action=gaction(label='  Set to Default Values  ', icon='refresh', handler=set.defaults), 
								container=button.g)
	e$adv.set.okbutton <- bDem.gbutton('OK', container=button.g)
	
		if(!is.null(h$action$env$params)) { # OK button previously clicked 
			for (par in param.names) 
				svalue(e[[par]]) <- paste(h$action$env$params[[par]], collapse=', ')
			
		}
		
	visible(adv.set.win) <- TRUE
	e$window.destroyed <- FALSE
	h$action$env$adv.set.env <- e
	
	addHandlerDestroy(adv.set.win, 
					handler=function(h1, ...) h$action$env$adv.set.env$window.destroyed <- TRUE)

	}
	if(!is.null(h$action$env$adv.set.okhandler)) 
		removehandler(h$action$env$adv.set.env$adv.set.okbutton, h$action$env$adv.set.okhandler)
	h$action$env$adv.set.okhandler <- addhandlerclicked(h$action$env$adv.set.env$adv.set.okbutton, 
												handler=set.advance.pars)
}

.create.extra.TS.group <- function(g, main.win, e, defaults, my.file.item='my.tfr.file', type='tfr', label.type='TFR') {
	g2 <- ggroup(horizontal=TRUE, container=g)
	tfr.g <- gframe(paste("<span color='blue'>", label.type, " time series</span>", sep=''), markup=TRUE, 
					horizontal=FALSE, container=g2)
	tlo <- glayout(container=tfr.g)
	tlo[1,1:2] <- bDem.gbutton(paste("  Select countries/regions from the UN ", label.type, " dataset  ", sep=''), 
				container=tlo, handler=multiSelectCountryMenu,
				action=list(mw=main.win, env=e, type=type, label.widget.name='extra.country.label'))
	# For showing selected countries
	tlo[2,1:2, anchor=c(-1,0)] <- e$extra.country.label <- glabel('', container=tlo)
	tlo[3,1, anchor=c(-1,0)] <- glabel(paste("User-defined ", label.type, " file:", sep=''), container=tlo)
	tlo[3,2] <- e[[my.file.item]] <- bDem.gfilebrowse(eval(defaults[[my.file.item]]), type='open', 
					  width=30, quote=FALSE, container=tlo)
	tooltip(e[[my.file.item]]) <- "Overwrites default WPP data."
}

.create.extra.mcmc.process <- function(g, e, defaults) {
	leftcenter <- c(-1,0)
	g2 <- ggroup(horizontal=TRUE, container=g)		  
	iter.g <- gframe("<span color='blue'>MCMC</span>", markup=TRUE, horizontal=TRUE, container=g2)
	itlo <- glayout(container=iter.g)
	itlo[1,1] <- glabel("Number of iterations:", container=itlo)
	itlo[1,2] <- e$iter <- gedit(defaults$iter, width=7, container=itlo)
	tooltip(e$iter) <- "Leave empty if the same as the main simulation."
	itlo[2,1] <- glabel("Thin:", container=itlo)
	itlo[2,2] <- e$thin <- gedit(defaults$thin, width=2, container=itlo)
	tooltip(e$thin) <- "Thinning interval for sampling from the hyperparameters."
	itlo[3,1] <- glabel("Burnin:", container=itlo)
	itlo[3,2] <- e$burnin <- gedit(defaults$burnin, width=4, container=itlo)
	tooltip(e$burnin) <- "Burnin for sampling from the hyperparameters."
		
	addSpace(g2,10)			  					  
	paral.g <- gframe("<span color='blue'>Process control</span>", markup=TRUE, horizontal=TRUE, container=g2)
	pclo <- glayout(container=paral.g)
	pclo[1,1] <- e$verbose <- gcheckbox("Verbose", checked=defaults$verbose, container=pclo)
	#addSpace(e$paral.g, 10)
	pclo[1,2, anchor=leftcenter] <- "Verbose iter:"
	pclo[1,3] <- e$verbose.iter <- gedit(defaults$verbose.iter, width=4, container=pclo)
	tooltip(e$verbose.iter) <- 'Frequency of outputting the process status.'
	pclo[2,1] <- e$parallel <- gcheckbox("Parallel", checked=defaults$parallel, container=pclo)
	tooltip(e$parallel) <- 'Chains will be processed in parallel.'
	pclo[2,2, anchor=leftcenter] <- glabel("Number of nodes:", container=pclo)
	pclo[2,3] <- e$nr.nodes <- gedit(defaults$nr.nodes, width=2, container=pclo)
	tooltip(e$nr.nodes) <- 'Used if Parallel is checked. Set it to number of chains in main simulation.'
}

.create.status.label <- function(g, e) {
	statusg <- ggroup(container=g, horizontal=TRUE)
	e$statuslabel <- glabel('', container=statusg)
	addSpring(statusg)
}

mcmc.extra.countries.group <- function(g, main.win, parent) {
	e <- new.env()
	defaults <- formals(run.tfr.mcmc.extra) # default argument values
	e$sim.dir <- parent$sim.dir
	addSpace(g, 10)
	.create.extra.TS.group(g, main.win, e, defaults)			
	addSpace(g, 10)
	.create.extra.mcmc.process(g, e, defaults)
	addSpace(g, 10)
	.create.status.label(g, e)
	addSpring(g)
	button.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='run.tfr.mcmc.extra', package='bayesTFR', parent.group=button.g,
						parent.window=main.win)	
	addSpring(button.g)
	create.generate.script.button(handler=mcmc.run.extra, action=list(mw=main.win, env=e, script=TRUE), container=button.g)
	bDem.gbutton(action=gaction(label=' Run MCMC ', icon='execute', handler=mcmc.run.extra, 
				action=list(mw=main.win, env=e, script=FALSE)), container=button.g)
	return(e)
}

mcmc.run.extra <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	param.names <- list(numeric=c('nr.nodes', 'iter', 'thin', 'burnin', 'verbose.iter'),
						text=c('sim.dir', 'my.tfr.file'),
						logical=c('verbose', 'parallel'))
	params <- get.parameters(param.names, e, quote=h$action$script)
	params[['countries']] <- e$selected.extra.countries
	if (h$action$script) {
		cmd <- paste('run.tfr.mcmc.extra(', assemble.arguments(params), ')',sep=' ')
		create.script.widget(cmd, h$action$mw, package="bayesTFR")
	} else {
		.run.simulation(e, handler=get.tfr.simulation.extra.status, option='bDem.TFRmcmcExtra', 
								call='run.tfr.mcmc.extra', params=params, 
								action=list(sb=e$statuslabel),
								sim.name='TFR MCMC extra simulation', main.win=h$action$mw,
								interval=1000)
	}
}


get.table.of.countries.from.locfile <- function(sim.dir, sorted=TRUE, type='tfr') {
	mcmc.set <- do.call(paste('get.', type, '.mcmc', sep=''), list(sim.dir=sim.dir))
	if(is.null(mcmc.set)) {
		gmessage('Simulation directory contains no valid MCMC results.', title='Input Error',
					icon='error')
		return(NULL)
	}
	wpp.year <- mcmc.set$meta$wpp.year
	# get the UN dataset
	tfr.data <- do.call(paste('get.', type, '.UN.data', sep=''), list(meta=mcmc.set$meta))
	codes <- tfr.data[,'country_code']
	# filter out countries used already for an estimation
	used.codes <- mcmc.set$meta$regions$country_code[1:(get.nr.countries.est(mcmc.set$meta))]
	codes <- codes[!is.element(codes, used.codes)] # codes not used in the estimation
			
	# get the UN location dataset
	loc.data <- bayesTFR:::load.bdem.dataset('UNlocations', wpp.year)
	loc.data <- loc.data[,c("country_code", "name")]
	colnames(loc.data)[1] <- 'code'
	
	#include only those that are contained in codes
	loc.data <- loc.data[is.element(loc.data[,'code'], codes),]
	loc.data[,'name'] <- gsub(' +$', '', loc.data[,'name']) # remove trailing space
	loc.data[,'name'] <- gsub('^ +', '', loc.data[,'name']) # remove spaces at the beginning
	if(sorted) {
		ord.idx <- order(loc.data[,'name'])
		loc.data <- loc.data[ord.idx,]
	}
	return(loc.data)
}

multiSelectCountryMenu <- function(h, ...) {
	country.selected <- function(h1, ...) {
		h$action$env$selected.extra.countries <- svalue(h$action$env$sel.extra.country.gt)
		visible(h$action$env$extra.country.sel.win) <- FALSE
		if(!is.null(h$action$label.widget.name) && !is.null(h$action$env[[h$action$label.widget.name]])) {
			svalue(h$action$env[[h$action$label.widget.name]]) <- paste(h$action$env$selected.extra.countries, collapse=',')
		}
		if(length(h$action$env$selected.extra.countries) == 0) h$action$env$selected.extra.countries <- NULL
	}
	new.window <- TRUE
	if (!is.null(h$action$env$extra.country.sel.win)) {
		# if anything has changed (sim.dir or the data), the window needs to be re-built
		if (svalue(h$action$env$sim.dir) != h$action$env$sim.dir.used) {
			dispose(h$action$env$extra.country.sel.win)
			new.window <- TRUE
		} else {
			extra.country.table <- get.table.of.countries.from.locfile(
												sim.dir=svalue(h$action$env$sim.dir),
												sorted=FALSE, 
												type=if(is.null(h$action$type)) 'tfr' else h$action$type)
			if(dim(extra.country.table)[1] != dim(h$action$env$extra.country.table)[1]) {
				dispose(h$action$env$extra.country.sel.win)
				new.window <- TRUE
			} else {
				new.window <- FALSE
				visible(h$action$env$extra.country.sel.win) <- TRUE
			}
		}
	}
	if(new.window) {
		sim.dir.used <- svalue(h$action$env$sim.dir)
		country.table <- get.table.of.countries.from.locfile(sim.dir=sim.dir.used, sorted=FALSE,
								type=if(is.null(h$action$type)) 'tfr' else h$action$type)
		if (is.null(country.table)) return(NULL)
		h$action$env$sim.dir.used <- sim.dir.used
		h$action$env$extra.country.table <- country.table
		h$action$env$extra.country.sel.win <- win <- gwindow('Select countries & regions', 
							parent=h$action$mw, height=450,
							handler=function(h, ...) {
								h$action$env$extra.country.sel.win<-NULL;
								h$action$env$sel.extra.country.ok.handler <- NULL
							},
							action=list(env=h$action$env))
		t.group <- ggroup(horizontal=FALSE, container=win)
		h$action$env$sel.extra.country.gt <- gtable(h$action$env$extra.country.table, container=t.group, 
					expand=TRUE, multiple=TRUE, handler=country.selected)
		b.group <- ggroup(horizontal=TRUE, container=t.group)
		gbutton('Cancel', container=b.group, handler=function(h, ...) 
					visible(win) <- FALSE)
		addSpring(b.group)
		h$action$env$sel.extra.country.okbutton <- gbutton('OK', container=b.group)
	}
	if(!is.null(h$action$env$sel.extra.country.ok.handler)) 
		removehandler(h$action$env$sel.extra.country.okbutton, h$action$env$sel.extra.country.ok.handler)
	h$action$env$sel.extra.country.ok.handler <- addhandlerclicked(
						h$action$env$sel.extra.country.okbutton, handler=country.selected)
}

mcmc3.advance.settings <- function(h, ...) {
	base.par.names <- c('mu', 'rho', 'sigma.mu', 'sigma.rho', 'sigma.eps')
	param.names <- c(paste(base.par.names, 'ini', sep='.'), 
					paste(base.par.names, 'ini.low', sep='.'),
					paste(base.par.names, 'ini.up', sep='.'),
					paste(base.par.names, 'prior.low', sep='.'),
					paste(base.par.names, 'prior.up', sep='.'))

	get.defaults <- function() {
		all.defaults <- formals(run.tfr3.mcmc) # default argument values
		defaults <- list()
		for (par in param.names) {
			if(is.element(par, names(all.defaults)))
				defaults[[par]] <- eval(all.defaults[[par]])
		}
		for(par in base.par.names) {
			defaults[[paste(par, 'ini.low', sep='.')]] <- eval(all.defaults[[paste(par, 'prior.range', sep='.')]])[1]
			defaults[[paste(par, 'ini.up', sep='.')]] <- eval(all.defaults[[paste(par, 'prior.range', sep='.')]])[2]
			defaults[[paste(par, 'prior.low', sep='.')]] <- eval(all.defaults[[paste(par, 'prior.range', sep='.')]])[1]
			defaults[[paste(par, 'prior.up', sep='.')]] <- eval(all.defaults[[paste(par, 'prior.range', sep='.')]])[2]
		}
		return(defaults)
	}

	set.defaults <- function(h2, ...) {
		defaults <- get.defaults()
		for (par in param.names) {
			if(!is.null(h$action$env$adv.set.env[[par]])) {
				if(length(defaults[[par]])>1) {
					for (i in 1:length(defaults[[par]])) svalue(h$action$env$adv.set.env[[par]][[i]]) <- defaults[[par]][i]
				} else 
					svalue(h$action$env$adv.set.env[[par]]) <- defaults[[par]]
			}
		}
		widget.defaults <- h$action$env$adv.set.env$widget.defaults
		for (par in names(widget.defaults)) {
			if(length(widget.defaults[[par]])>1) {
				for (i in 1:length(widget.defaults[[par]])) svalue(h$action$env$adv.set.env[[par]][[i]]) <- widget.defaults[[par]][i]
			} else svalue(h$action$env$adv.set.env[[par]]) <- widget.defaults[[par]]
		}
	}
	
	set.advance.pars <- function(h2, ...) {
		# The order in which the widgets are handled must be the same as in function init.widget.value.pairs.
		params <- widget.value <- list()
		counter <- 1		
		linked.pars.list <- h$action$env$adv.set.env$linked.pars.list
		for(par in names(linked.pars.list)) {
			if(is.list(linked.pars.list[[par]])) {
				for(i in 1:length(linked.pars.list[[par]])) {
					value <- svalue(linked.pars.list[[par]][[i]])
					params[[par]] <- c(params[[par]], if(nchar(value)==0) NULL else as.numeric(value))
					h$action$env$adv.set.env$widget.value.pairs[[counter]][[2]] <- value
					#widget.value <- c(widget.value, list(c(linked.pars.list[[par]][[i]], value)))
					counter <- counter + 1
				} 
			} else {
				value <- svalue(linked.pars.list[[par]])
				params[[par]] <- if(nchar(value)==0) NULL else as.numeric(value)
				h$action$env$adv.set.env$widget.value.pairs[[counter]][[2]] <- value
				counter <- counter + 1
				#widget.value <- c(widget.value, list(c(linked.pars.list[[par]], value)))
			}
		}
		h$action$env$params <- params
		#h$action$env$adv.set.env$widget.value.pairs <- widget.value
		visible(h$action$env$adv.set.win) <- FALSE
	}
	
	init.widget.value.pairs <- function(e) {
		# The order in which the widgets are handled must be the same as in function set.advance.pars.
		widget.value <- list()		

		for(par in names(e$linked.pars.list)) {
			if(is.list(e$linked.pars.list[[par]])) {
				for(i in 1:length(e$linked.pars.list[[par]])) 
					widget.value <- c(widget.value, list(c(e$linked.pars.list[[par]][[i]], 0)))
			} else {
				widget.value <- c(widget.value, list(c(e$linked.pars.list[[par]], 0)))
			}
		}
		e$widget.value.pairs <- widget.value
	}

		
	if (!is.null(h$action$env$adv.set.win) && !h$action$env$adv.set.env$window.destroyed) { #Advanced Parameters window exists
		if(!is.null(h$action$env$params)) { # OK button previously clicked 
			for (i in 1:length(h$action$env$adv.set.env$widget.value.pairs)) 
				svalue(h$action$env$adv.set.env$widget.value.pairs[[i]][[1]]) <- h$action$env$adv.set.env$widget.value.pairs[[i]][[2]]
		} else { # OK button not clicked yet, values are set to defaults
			#for (par in param.names) {
			#	svalue(h$action$env$adv.set.env[[par]]) <- h$action$env$adv.set.env$defaults[[par]]
			#}
			set.defaults(h)
		}
		visible(h$action$env$adv.set.win) <- TRUE
	} else { # create the Advanced Parameters window
		h$action$env$adv.set.win <- adv.set.win <- bDem.gwindow('Settings for Bayesian Hierarchical Model Phase III TFR',
						parent=h$action$mw, visible=FALSE,
						handler=function(h, ...) {
							h$action$env$adv.set.okhandler <- NULL
						})
		e <- new.env()
		e$defaults <- defaults <- get.defaults()
		e$adv.g <- ggroup(container=adv.set.win, horizontal=FALSE)
	
	linked.pars.list <-  widget.defaults <- list()

	priors.f <- gframe("<span color='blue'>Uniform prior parameters and initial values</span>", markup=TRUE, 
						container=e$adv.g, horizontal=FALSE)
	uni.g <- ggroup(horizontal=TRUE, container=priors.f)
	uni.flo <- glayout(container=uni.g)
	
	l <- 1 # row 1
	uni.flo[l,1] <- ''
	uni.flo[l,2] <- '   prior\nlower b.'
	uni.flo[l,3] <- '   prior\nupper b.'
	uni.flo[l,4] <- ' init\nlower'
	uni.flo[l,5] <- ' init\nupper'
	uni.flo[l,6] <- 'init values'
	lower <- c(defaults$mu.ini.low, defaults$rho.ini.low, defaults$sigma.mu.ini.low, 
					defaults$sigma.rho.ini.low,defaults$sigma.eps.ini.low)
	upper <- c(defaults$mu.ini.up, defaults$rho.ini.up, defaults$sigma.mu.ini.up, 
					defaults$sigma.rho.ini.up,defaults$sigma.eps.ini.up)
	ini <- c(defaults$mu.ini, defaults$rho.ini, defaults$sigma.mu.ini, 
					defaults$sigma.rho.ini,defaults$sigma.eps.ini)
	prior.low <- c(defaults$mu.prior.low, defaults$rho.prior.low, defaults$sigma.mu.prior.low, 
					defaults$sigma.rho.prior.low,defaults$sigma.eps.prior.low)
	prior.up <- c(defaults$mu.prior.up, defaults$rho.prior.up, defaults$sigma.mu.prior.up, 
					defaults$sigma.rho.prior.up,defaults$sigma.eps.prior.up)

	e$lower <- e$upper <- e$prior.low <- e$prior.up <-e$init <- NULL
	labels <- paste('<span>', c('mu', 'rho', 'sigma<sub>mu</sub>', 'sigma<sub>rho</sub>', 'sigma<sub>eps</sub>'), ':</span>', sep='')
	for (i in 1:5) {
		row <- l + i
		uni.flo[row,1] <- glabel(labels[i], markup=TRUE, container=uni.flo)
		e$prior.low <- c(e$prior.low, uni.flo[row,2] <- gedit(prior.low[i], width=5, container=uni.flo))
		e$prior.up <- c(e$prior.up, uni.flo[row,3] <- gedit(prior.up[i], width=5, container=uni.flo))
		e$lower <- c(e$lower, uni.flo[row,4] <- gedit(lower[i], width=5, container=uni.flo))
		e$upper <- c(e$upper, uni.flo[row,5] <-  gedit(upper[i], width=5, container=uni.flo))
		e$init  <- c(e$init, uni.flo[row,6] <- gedit(ini[i], width=10, container=uni.flo))
		addHandlerChanged(e$init[[i]], action=list(idx=i), handler=function(h1,...) {
						isempty <- nchar(svalue(e$init[[h1$action$idx]]))==0
						enabled(e$lower[[h1$action$idx]]) <- isempty
						enabled(e$upper[[h1$action$idx]]) <- isempty
						})
	}
	widget.defaults[['prior.low']] <- prior.low
	widget.defaults[['prior.up']] <- prior.up
	widget.defaults[['lower']] <- lower
	widget.defaults[['upper']] <- upper
	widget.defaults[['init']] <- if(is.null(ini)) rep('', 5) else ini
	
	ipar <- 1
	for(par in base.par.names) {
		linked.pars.list[[paste(par, 'ini', sep='.')]] <- e$init[[ipar]]
		linked.pars.list[[paste(par, 'ini.range', sep='.')]] <- c(e$lower[[ipar]], e$upper[[ipar]])
		linked.pars.list[[paste(par, 'prior.range', sep='.')]] <- c(e$prior.low[[ipar]], e$prior.up[[ipar]])
		ipar <- ipar+1
	}		

	addSpring(priors.f)
	info.ini.group <- ggroup(horizontal=FALSE, container=priors.f)
	glabel('NOTE:', container=info.ini.group, anchor=c(-1,0))
	glabel('Leave "init values" blank for the starting values being equally', container=info.ini.group, anchor=c(-1,0))
	glabel('distributed between "init lower" and "init upper". For specific', container=info.ini.group, anchor=c(-1,0))
	glabel('initial values enter one value per chain separated by commas.', container=info.ini.group, anchor=c(-1,0))
	
	# Buttons
	button.g <- ggroup(container=e$adv.g, horizontal=TRUE)
	bDem.gbutton('Cancel', container=button.g, handler=function(h, ...) 
					visible(adv.set.win) <- FALSE)
	addSpring(button.g)
	e$adv.set.defaultbutton <- bDem.gbutton(action=gaction('  Set to Default Values  ', icon='refresh', handler=set.defaults),
										container=button.g)
	e$adv.set.okbutton <- bDem.gbutton('OK', container=button.g)
	
	e$linked.pars.list <- linked.pars.list
	e$widget.defaults <- widget.defaults
	init.widget.value.pairs(e)
	
		if(!is.null(h$action$env$params)) { # OK button previously clicked 
			for (i in 1:length(h$action$env$adv.set.env$widget.value.pairs)) {
				svalue(e$widget.value.pairs[[i]][[1]]) <- h$action$env$adv.set.env$widget.value.pairs[[i]][[2]]
				e$widget.value.pairs[[i]][[2]] <- h$action$env$adv.set.env$widget.value.pairs[[i]][[2]]
			}
		} 
		
	visible(adv.set.win) <- TRUE
	e$window.destroyed <- FALSE
	h$action$env$adv.set.env <- e
	
	addHandlerDestroy(adv.set.win, 
					handler=function(h1, ...) h$action$env$adv.set.env$window.destroyed <- TRUE)
	
	}
	if(!is.null(h$action$env$adv.set.okhandler)) 
		removehandler(h$action$env$adv.set.env$adv.set.okbutton, h$action$env$adv.set.okhandler)
	h$action$env$adv.set.okhandler <- addhandlerclicked(h$action$env$adv.set.env$adv.set.okbutton, 
											handler=set.advance.pars)

}