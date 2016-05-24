TFRcontinueMCMCgroup <- function(g, main.win, parent=NULL) {
	e <- new.env()
	defaults <- formals(continue.tfr.mcmc) # default argument values
	defaults3 <- formals(continue.tfr3.mcmc)
	e$output.dir <- parent$sim.dir
	e$run.prediction <- FALSE
	addSpace(g, 10)
	phase.g <- gframe("", markup=TRUE, horizontal=FALSE, spacing=10, container=g)
	phase.g1 <- ggroup(horizontal=TRUE, container=phase.g)
	glabel("Run TFR MCMC for ", container=phase.g1)
	e$mcmc.type <- gradio(c('Phase II', 'Phase III'), horizontal=TRUE, container=phase.g1,
								handler=function(h,...){
									val <- svalue(h$obj, index=TRUE)
									enabled(phaseII.g) <- val == 1
									enabled(phaseIII.g) <- val == 2
									svalue(type.nb) <- val
								})
	type.nb <- gnotebook(container=phase.g, expand=TRUE)
	set.widget.bgcolor(type.nb, color.main)
	set.widget.basecolor(type.nb, color.nb.inactive)
	phaseII.g <- ggroup(label="<span color='darkred'>Phase II</span>", markup=TRUE, horizontal=FALSE, container=type.nb)
	addSpace(phaseII.g, 10)
	.create.cont.mcmc.process.group(phaseII.g, e, main.win, defaults, mcmc.label="MCMC Phase II", process.label="Process control for Phase II")
	phaseIII.g <- ggroup(label="<span color='darkred'>Phase III</span>", markup=TRUE, horizontal=FALSE, container=type.nb)
	addSpace(phaseIII.g, 10)
	e$phase3 <- new.env()
	e$phase3$output.dir <- e$output.dir
	.create.cont.mcmc.process.group(phaseIII.g, e$phase3, main.win, defaults3, type='tfr3', mcmc.label="MCMC Phase III",
								process.label="Process control for Phase III")
	svalue(type.nb) <- 1
	enabled(phaseIII.g) <- svalue(e$mcmc.type, index=TRUE) == 2
	addSpace(g, 10)
	.create.status.label(g, e)

	addSpring(g)
	cont.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='continue.tfr.mcmc', package='bayesTFR', parent.group=cont.g,
						parent.window=main.win)
	addSpring(cont.g)
	create.generate.script.button(handler=mcmc.continue, action=list(mw=main.win, env=e, script=TRUE),
								container=cont.g)
	bDem.gbutton(action=gaction(label=' Continue MCMC ', icon='execute', handler=mcmc.continue, 
				action=list(mw=main.win, env=e, script=FALSE)), container=cont.g)
	return(e)
}

.create.cont.mcmc.process.group <- function(g, e, main.win, defaults, type='tfr', mcmc.label="MCMC", process.label="Process control") {
	g2 <- ggroup(horizontal=TRUE, container=g)
	.create.autoconf.cont.group(g2, e, main.win, defaults, type=type, label=mcmc.label)
	.create.process.group(g2, e, defaults, show.buffer.size=FALSE, label=process.label)
}

.create.autoconf.cont.group <- function(g, e, main.win, defaults, type='tfr', label='MCMC') {
	leftcenter <- c(-1,0)
	mcmc.g <- gframe(paste("<span color='blue'>", label, "</span>"), markup=TRUE, horizontal=FALSE, spacing=10, container=g)
	mclo <- glayout(container=mcmc.g)
	mclo[1,1] <- e$run.auto <- gcheckbox("Auto simulation", checked=FALSE, 
							container=mclo, handler=function(h,...){.enable.auto.cont(svalue(h$obj), e)})
	mclo[1,2] <- e$auto.conf.b <- bDem.gbutton(' Configure auto run ', container=mclo, handler=configure.auto.run, 
				action=list(mw=main.win, env=e, cont.run=TRUE, type=type))
	mclo[2,1, anchor=leftcenter] <- glabel("Number of iterations:", container=mclo)
	mclo[2,2] <- e$iter <- gedit(defaults$iter, width=7, container=mclo)
	tooltip(e$iter) <- "How many iterations to add to each chain."
	mclo[3,1, anchor=leftcenter] <- glabel("Chain ids:", container=mclo)
	mclo[3,2] <- e$chain.ids <- gedit(defaults$chain.ids, width=10, container=mclo)
	tooltip(e$chain.ids) <- "Comma-separated ids of the chains. Leave empty if all chains should be continued."
	.enable.auto.cont(FALSE, e)
}

.enable.auto.cont <- function(enable, e) {
	enabled(e$auto.conf.b) <- enable
	enabled(e$chain.ids) <- !enable
	enabled(e$iter) <- !enable 
}

.get.defaults.for.auto.cont.tfr <- function(e) {
	mcmc.set <- get.tfr.mcmc(sim.dir=svalue(e$output.dir))
	if(is.null(mcmc.set)) {
		gmessage('Simulation directory contains no valid phase II TFR MCMCs.', title='Input Error',
					icon='error')
		return(NULL)
	}
	return(mcmc.set$meta$auto.conf)
}

.get.defaults.for.auto.cont.tfr3 <- function(e) {
	mcmc.set <- get.tfr3.mcmc(sim.dir=svalue(e$output.dir))
	if(is.null(mcmc.set)) {
		gmessage('Simulation directory contains no valid phase III TFR MCMCs.', title='Input Error',
					icon='error')
		return(NULL)
	}
	return(mcmc.set$meta$auto.conf)
}

mcmc.continue <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(output.dir='Simulation directory'), env=e)) return()
	params.sim.dir <- get.parameters(list(text=c('output.dir')), e, quote=h$action$script)
	outdir <- get.parameters(list(text='output.dir'), e, quote=FALSE)$output.dir # to avoid double quoting if script is TRUE
	
	param.names <- list(numeric=c('iter', 'nr.nodes'),
						logical=c('verbose', 'parallel'),
						numvector=c('chain.ids'))
	mcmc.type <- svalue(e$mcmc.type, index=TRUE)
	env <- if(mcmc.type==1) e else e$phase3
	params <- get.parameters(param.names, env, quote=h$action$script)
	run.auto <- svalue(env$run.auto)
	if (!run.auto && !has.required.arguments(list(iter='Number of iterations'), env=env)) return()
	if (run.auto) {
		params[['auto.conf']] <- env$auto.conf
		params[['iter']] <- if (h$action$script) sQuote('auto') else 'auto'
	}
	if(mcmc.type==1) { # phase II
		m <- get.tfr.mcmc(outdir)
		if(is.null(m)) {
			gmessage('Simulation directory contains no valid MCMC results.', 
							title='Input Error', icon='error')
			return()
		}
	} else { # phase III
		m <- get.tfr3.mcmc(outdir)
		if(is.null(m)) {
			gmessage('Simulation directory contains no valid phase III MCMC results.', 
							title='Input Error', icon='error')
			return()
		}
	}
	if (h$action$script) {
		if(mcmc.type==1)
			commands <- paste('m <- continue.tfr.mcmc(', assemble.arguments(c(params, params.sim.dir)), ')',sep=' ')
		else
			commands <- paste('m <- continue.tfr3.mcmc(', params.sim.dir$output.dir, ', ', assemble.arguments(params), ')', sep=' ')
		if(run.auto && env$run.prediction) 
			commands <- paste(commands, '\n\ntfr.predict(sim.dir=,', params.sim.dir$output.dir, ', use.diagnostics=TRUE, replace.output=TRUE)', sep='')
		create.script.widget(commands, h$action$mw, package="bayesTFR")
	} else {  # run a simulation
		if(mcmc.type==1) # phase II
			m <- .run.simulation(e, handler=get.tfr.simulation.status, option='bDem.TFRmcmc', 
								call='continue.tfr.mcmc', params=c(params, params.sim.dir), 
								sim.name='TFR MCMC Phase II', main.win=h$action$mw,
								action=list(sb=e$statuslabel, sim.dir=params.sim.dir$output.dir),
								interval=5000)
		else # phase III
			m <- .run.simulation(e, handler=get.tfr.simulation.status, option='bDem.TFRmcmc', 
								call='continue.tfr3.mcmc', params=c(list(sim.dir=params.sim.dir$output.dir), params), 
								sim.name='TFR MCMC Phase III', main.win=h$action$mw,
								action=list(sb=e$statuslabel, sim.dir=params.sim.dir$output.dir),
								interval=5000)
		if(run.auto && env$run.prediction)
			.run.prediction(e, handler=get.tfr.prediction.status, option='bDem.TFRpred', 
								call='tfr.predict', params=list(sim.dir=params.sim.dir$output.dir, use.diagnostics=TRUE, replace.output=TRUE), 
								sim.name='TFR prediction', main.win=h$action$mw,
								action=list(sb=e$statuslabel),
								interval=1000)
	}
}
