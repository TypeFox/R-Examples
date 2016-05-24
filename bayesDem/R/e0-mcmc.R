e0RunMCMCgroup <- function(g, main.win, parent) {
	# parent needed for wpp set in the parent function
	e <- new.env()
	nb <- bDem.gnotebook(container=g, expand=TRUE)
	all.c.g <- ggroup(label="<span color='#0B6138'>All Countries</span>", markup=TRUE, 
						horizontal=FALSE, container=nb)
	e0mcmc.all.countries.group(all.c.g, main.win, parent)
	extra.c.g <- ggroup(label="<span color='#0B6138'>Extra Areas &amp; Regions</span>", 
						markup=TRUE, horizontal=FALSE, container=nb)
	e0mcmc.extra.countries.group(extra.c.g, main.win, parent)
	svalue(nb) <- 1
}

e0mcmc.all.countries.group <- function(g, main.win, parent) {
	e <- new.env()
	defaults <- formals(run.e0.mcmc) # default argument values
	e$output.dir <- parent$sim.dir
	addSpace(g, 10)
	time.g <- gframe("<span color='blue'>e0 time series</span>", markup=TRUE, horizontal=FALSE, container=g, expand=FALSE)
	timelo <- glayout(container=time.g)
	timelo[1,1, anchor=c(-1,0)] <- glabel("Sex:", container=timelo)
	timelo[1,2] <- e$sex <- bDem.gdroplist(c('Female', 'Male'), container=timelo, selected=1)
	timelo[1,3, anchor=c(-1,0)] <- glabel("Start year:", container=timelo)
	timelo[1,4] <- e$start.year <- gedit(defaults$start.year, width=4, container=timelo)
	tooltip(e$start.year) <- 'Historical data prior to this year will be ignored.'
	timelo[1,5, anchor=c(-1,0)] <- glabel("Present year:", container=timelo)
	timelo[1,6] <- e$present.year <- gedit(defaults$present.year, width=4, container=timelo)
	tooltip(e$present.year) <- 'Historical data after this year will be ignored.'
	timelo[1,7, anchor=c(-1,0)] <- glabel("WPP year:", container=timelo)
	timelo[1,8, anchor=c(-1,0)] <- wpp <- glabel(parent$wpp.year, container=timelo)
	tooltip(wpp) <- 'To change this start bayesDem with wpp.year as second argument.'
	timelo[2,1:2, anchor=c(-1,0)] <- glabel("User-defined e0 file:", container=timelo)
	timelo[2,3:8] <- e$my.e0.file <- bDem.gfilebrowse(eval(defaults$my.e0.file), type='open', 
					  width=30, quote=FALSE, container=timelo)
	tooltip(e$my.e0.file) <- 'Overwrites default wpp data.'
	addSpace(g, 10)
	.create.mcmc.process.group(g, e, main.win, defaults, type='e0', advance.settigs.function=e0mcmc.advance.settings)
	addSpace(g, 10)
	.create.status.label(g, e)
	addSpring(g)
	adv.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic=c('run.e0.mcmc', 'bayesLife-package'), package='bayesLife', parent.group=adv.g,
						parent.window=main.win)
	addSpring(adv.g)
	create.generate.script.button(handler=e0mcmc.run, action=list(mw=main.win, env=e, script=TRUE, wpp.year=parent$wpp.year),
								container=adv.g)
	bDem.gbutton(action=gaction(label=' Run MCMC ', icon='execute', handler=e0mcmc.run, 
				action=list(mw=main.win, env=e, script=FALSE, wpp.year=parent$wpp.year, parent.group=g)), 
				container=adv.g)
	#e$statusbar <- gstatusbar() # don't display now
	return(e)
	}

e0mcmc.extra.countries.group <- function(g, main.win, parent) {
	e <- new.env()
	defaults <- formals(run.e0.mcmc.extra) # default argument values
	e$sim.dir <- parent$sim.dir
	addSpace(g, 10)
	.create.extra.TS.group(g, main.win, e, defaults, my.file.item='my.e0.file', type='e0', label.type='e0')
	addSpace(g, 10)
	.create.extra.mcmc.process(g, e, defaults)
	addSpace(g, 10)
	.create.status.label(g, e)
	addSpring(g)
	button.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='run.e0.mcmc.extra', package='bayesLife', parent.group=button.g,
						parent.window=main.win)	
	addSpring(button.g)
	create.generate.script.button(handler=e0mcmc.run.extra, action=list(mw=main.win, env=e, script=TRUE), container=button.g)
	bDem.gbutton(action=gaction(label=' Run MCMC ', icon='execute', handler=e0mcmc.run.extra, 
				action=list(mw=main.win, env=e, script=FALSE)), container=button.g)

	return(e)
}

e0mcmc.run <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(output.dir='Simulation directory'), env=e)) return()
	param.names <- list(numeric=c('buffer.size', 'nr.nodes', 'iter', 'thin', 'nr.chains', 'start.year', 
									'present.year', 'seed'),
						text=c('output.dir', 'my.e0.file', 'sex', 'compression.type'),
						logical=c('verbose', 'parallel'))
	params <- get.parameters(param.names, e, quote=h$action$script)
	params[['wpp.year']] <- h$action$wpp.year
	run.auto <- svalue(e$run.auto)
	if (run.auto) {
		params[['auto.conf']] <- e$auto.conf
		op <- options("useFancyQuotes")
		options(useFancyQuotes = FALSE)
		params[['iter']] <- if (h$action$script) sQuote('auto') else 'auto'
		options(op)
	}
	outdir <- get.parameters(list(text='output.dir'), e, quote=FALSE)$output.dir # to avoid double quoting if script is TRUE
	if(file.exists(outdir) && (length(list.files(outdir)) > 0)) {
		params[['replace.output']] <- FALSE
		if (gconfirm(paste('Non-empty directory', outdir, 
								'already exists.\nDo you want to overwrite existing results?'),
				icon='question', parent=h$action$mw))
			params[['replace.output']] <- TRUE
		else return(NULL)
	}

	if (h$action$script) {
		commands <- paste('m <- run.e0.mcmc(', assemble.arguments(c(params, e$params)), ')',sep=' ')
		if(run.auto && e$run.prediction) {
			commands <- paste(commands, '\n\ne0.predict(m, use.diagnostics=TRUE)', sep='')
			#commands <- paste(commands, '\n\ne0.predict(m)', sep='')
			}
		create.script.widget(commands, h$action$mw, package="bayesLife")
	} else {
		run <- FALSE
		if ((params[['iter']] == 'auto' && ((!is.null(params[['auto.conf']]) 
				&& params[['auto.conf']]$iter > 5000) || is.null(params[['auto.conf']]))) 
					|| (params[['iter']] != 'auto' && params[['iter']] > 5000)) {
			gconfirm('Running MCMC with these settings can take a very long time. Do you want to continue?',
					icon='question', parent=h$action$mw,
					handler=function(h, ...) run <<- TRUE)
		} else run <- TRUE
		if(run) {
			if(file.exists(params[['output.dir']])) unlink(params[['output.dir']], recursive=TRUE)
			m <- .run.simulation(e, handler=get.e0.simulation.status, option='bDem.e0mcmc', 
								call='run.e0.mcmc', params=c(params, e$params), 
								sim.name='e0 MCMC simulation', main.win=h$action$mw,
								action=list(sb=e$statuslabel, sim.dir=params[['output.dir']]),
								interval=5000)
			if(run.auto && e$run.prediction) {
				.run.prediction(e, type='e0', handler=get.e0.prediction.status, option='bDem.e0pred', 
								call='e0.predict', params=list(m, use.diagnostics=TRUE), 
								sim.name='e0 prediction', main.win=h$action$mw,
								action=list(sb=e$statuslabel),
								interval=1000)
			}
		}
	}
}

get.e0.simulation.status <- function(h, ...) {
	sb <- h$action$sb
	sim.dir <- h$action$sim.dir
	warn <- getOption('warn')
	options(warn=-1) # disable warning messages
	mcmc.set <- get.e0.mcmc(sim.dir)
	options(warn=warn)
	if(is.null(mcmc.set)) return()
	get.finished <- function(x) return(x$finished.iter)
	get.iter <- function(x) return(x$iter)
	finished <- sapply(mcmc.set$mcmc.list, get.finished)
	total <- sapply(mcmc.set$mcmc.list, get.iter)
	svalue(sb) <- paste('Running e0 simulation ... Iterations - ', 
				 	paste('chain ', 1:length(finished), ': ', finished, ' (', total, ')', sep='', collapse=', '))
	
}

get.e0.simulation.extra.status <- function(h, ...) {
	sb <- h$action$sb
	# We don't have any status info for the extra run, so just keep a static label
	svalue(sb) <- 'Running e0 extra simulation ... '
	
}


e0mcmc.run.extra <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	param.names <- list(numeric=c('nr.nodes', 'iter', 'thin', 'burnin'),
						text=c('sim.dir', 'my.e0.file'),
						logical=c('verbose', 'parallel'))
	params <- get.parameters(param.names, e, quote=h$action$script)
	params[['countries']] <- e$selected.extra.countries
	if (h$action$script) {
		cmd <- paste('run.e0.mcmc.extra(', assemble.arguments(params), ')', sep=' ')
		create.script.widget(cmd, h$action$mw, package="bayesLife")
	} else {
		.run.simulation(e, handler=get.e0.simulation.extra.status, option='bDem.e0mcmcExtra', 
								call='run.e0.mcmc.extra', params=params, 
								action=list(sb=e$statuslabel),
								sim.name='e0 MCMC extra simulation', main.win=h$action$mw,
								interval=1000)
	}
}

e0mcmc.advance.settings <- function(h, ...) {
	param.names <- list('a', 'delta', 'tau', 'Triangle.ini', 'k.ini', 'z.ini',
					'lambda.ini', 'lambda.k.ini', 'lambda.z.ini', 'omega.ini',
					'Triangle.ini.low', 'Triangle.ini.up', 
					'k.ini.low', 'z.ini.low', 'k.ini.up', 'z.ini.up',
					'Triangle.prior.low', 'Triangle.prior.up',
					'k.prior.low', 'z.prior.low', 'k.prior.up', 'z.prior.up',
					'lambda.ini.low', 'lambda.k.ini.low', 'lambda.z.ini.low', 
					'lambda.ini.up', 'lambda.k.ini.up', 'lambda.z.ini.up',
					'omega.ini.low', 'omega.ini.up',
					'Triangle.c.prior.low', 'Triangle.c.prior.up',
					'k.c.prior.low', 'z.c.prior.low', 'k.c.prior.up', 'z.c.prior.up',
					 'nu', 'dl.p1', 'dl.p2'
						)
	get.defaults <- function() {
		all.defaults <- formals(run.e0.mcmc) # default argument values
		defaults <- list()
		for (par in param.names) { 
			defaults[[par]] <- eval(all.defaults[[par]])
		}
		# special cases
		defaults$Triangle.c.ini.norm <- list(round(defaults$Triangle.ini.low + 
			(defaults$Triangle.ini.up - defaults$Triangle.ini.low)/2),c(2,2,2,2))
		defaults$k.c.ini.norm <- c(round(defaults$k.ini.low + (defaults$k.ini.up - defaults$k.ini.low)/2),2)
		defaults$z.c.ini.norm <- c(round(defaults$z.ini.low + (defaults$z.ini.up - defaults$z.ini.low)/2, 2), 0.2)
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
		array.pars <- c('a', 'delta', 'tau')
		counter <- 1
		for (i in 1:6) {
			for(par in array.pars) {
				value <- svalue(h$action$env$adv.set.env[[par]][[i]])
				params[[par]] <- c(params[[par]], as.numeric(value))
				h$action$env$adv.set.env$widget.value.pairs[[counter]][[2]] <- value
				#widget.value <- c(widget.value, list(c(h$action$env$adv.set.env[[par]][[i]], value)))
				counter <- counter + 1
			}
		}
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
		linked.pars.tuple <- h$action$env$adv.set.env$linked.pars.tuple
		for(par in names(linked.pars.tuple)) {
			for(item in 1:2) {
				if(is.list(linked.pars.tuple[[par]][[item]])) {
					if(item==1) params[[par]] <- list()
					params[[par]][[item]] <- rep(NA, length(linked.pars.tuple[[par]][[item]]))
					for(i in 1:length(linked.pars.tuple[[par]][[item]])) {
						value <- svalue(linked.pars.tuple[[par]][[item]][[i]])
						params[[par]][[item]][i] <-if(nchar(value)==0) NULL else as.numeric(value)
						h$action$env$adv.set.env$widget.value.pairs[[counter]][[2]] <- value
						counter <- counter + 1
						#widget.value <- c(widget.value, list(c(linked.pars.tuple[[par]][[item]][[i]], value)))
					} 
				} else {
					if(item==1) params[[par]] <- c()
					value <- svalue(linked.pars.tuple[[par]][[item]])
					params[[par]][item] <- if(nchar(value)==0) NULL else as.numeric(value)
					h$action$env$adv.set.env$widget.value.pairs[[counter]][[2]] <- value
					counter <- counter + 1
					#widget.value <- c(widget.value, list(c(linked.pars.tuple[[par]][[item]], value)))
				}
			}
		}
		rest.par.names <- setdiff(param.names, 
								union(union(array.pars, names(linked.pars.list)), names(linked.pars.tuple)))
		params <- c(params, get.parameters(list(numvector=rest.par.names), env=h$action$env$adv.set.env))
		for (par in rest.par.names) {
			h$action$env$adv.set.env$widget.value.pairs[[counter]][[2]] <- if(is.null(params[[par]])) '' else params[[par]]
			counter <- counter + 1
			#widget.value <- c(widget.value, list(c(h$action$env$adv.set.env[[par]], 
			#		if(is.null(params[[par]])) '' else params[[par]])))
		}
		h$action$env$params <- params
		#h$action$env$adv.set.env$widget.value.pairs <- widget.value
		visible(h$action$env$adv.set.win) <- FALSE
	}
	
	init.widget.value.pairs <- function(e) {
		# The order in which the widgets are handled must be the same as in function set.advance.pars.
		widget.value <- list()		
		array.pars <- c('a', 'delta', 'tau')
		for (i in 1:6) {
			for(par in array.pars) {
				widget.value <- c(widget.value, list(c(e[[par]][[i]], 0)))
			}
		}
		for(par in names(e$linked.pars.list)) {
			if(is.list(e$linked.pars.list[[par]])) {
				for(i in 1:length(e$linked.pars.list[[par]])) 
					widget.value <- c(widget.value, list(c(e$linked.pars.list[[par]][[i]], 0)))
			} else {
				widget.value <- c(widget.value, list(c(e$linked.pars.list[[par]], 0)))
			}
		}
		for(par in names(e$linked.pars.tuple)) {
			for(item in 1:2) {
				if(is.list(e$linked.pars.tuple[[par]][[item]])) {
					for(i in 1:length(e$linked.pars.tuple[[par]][[item]])) 
						widget.value <- c(widget.value, list(c(e$linked.pars.tuple[[par]][[item]][[i]], 0)))
				} else 
					widget.value <- c(widget.value, list(c(e$linked.pars.tuple[[par]][[item]], 0)))
			}
		}
		rest.par.names <- setdiff(param.names, 
								union(union(array.pars, names(e$linked.pars.list)), names(e$linked.pars.tuple)))
		for (par in rest.par.names) widget.value <- c(widget.value, list(c(e[[par]], 0)))
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
		h$action$env$adv.set.win <- adv.set.win <- bDem.gwindow('Settings for Bayesian Hierarchical Model of Life Expectancy',
						parent=h$action$mw, visible=FALSE,
						handler=function(h, ...) {
							h$action$env$adv.set.okhandler <- NULL
						})
		e <- new.env()
		e$defaults <- defaults <- get.defaults()
		e$adv.g <- ggroup(container=adv.set.win, horizontal=FALSE)
	
	linked.pars.list <- linked.pars.tuple <- widget.defaults <- list()

	priors.f <- gframe("<span color='blue'>Prior parameters and initial values</span>", markup=TRUE, container=e$adv.g, horizontal=FALSE)
	normal.g <- ggroup(horizontal=TRUE, container=priors.f)
	priors.normal.f <- gframe("<span  color='#0B6138'>Normal priors for world parameters</span>", markup=TRUE, container=normal.g,
								horizontal=FALSE)
	normal.flo <- glayout(container=priors.normal.f)
	
	l <- 1 # row 1
	normal.flo[l,1] <- ''
	normal.flo[l,2] <- 'prior mean\n     (a)'
	normal.flo[l,3] <- 'prior sd\n (delta)'
	normal.flo[l,4] <- '   prior\nlower b.'
	normal.flo[l,5] <- '   prior\nupper b.'
	normal.flo[l,6] <- ' init\nlower'
	normal.flo[l,7] <- ' init\nupper'
	normal.flo[l,8] <- 'init values'
	lower <- c(defaults$Triangle.ini.low, defaults$k.ini.low, defaults$z.ini.low)
	upper <- c(defaults$Triangle.ini.up, defaults$k.ini.up, defaults$z.ini.up)
	ini <- c(unlist(defaults$Triangle.ini), defaults$k.ini, defaults$z.ini)
	prior.low <- c(defaults$Triangle.prior.low, defaults$k.prior.low, defaults$z.prior.low)
	prior.up <- c(defaults$Triangle.prior.up, defaults$k.prior.up, defaults$z.prior.up)

	e$a <- e$delta <- e$lower <- e$upper <- e$prior.low <- e$prior.up <-e$init <- NULL
	labels <- paste('<span>', c(paste('Triangle<sub>', 1:4, '</sub>', sep=''), 'k', 'z'), ':</span>', sep='')
	for (i in 1:6) {
		row <- l + i
		normal.flo[row,1] <- glabel(labels[i], markup=TRUE, container=normal.flo)
		e$a <- c(e$a, normal.flo[row,2] <- gedit(round(defaults$a[i],4), width=7, container=normal.flo))
		e$delta <- c(e$delta, normal.flo[row,3] <- gedit(round(defaults$delta[i],3), width=5, container=normal.flo))
		e$prior.low <- c(e$prior.low, normal.flo[row,4] <- gedit(prior.low[i], width=5, container=normal.flo))
		e$prior.up <- c(e$prior.up, normal.flo[row,5] <- gedit(prior.up[i], width=5, container=normal.flo))
		e$lower <- c(e$lower, normal.flo[row,6] <- gedit(lower[i], width=5, container=normal.flo))
		e$upper <- c(e$upper, normal.flo[row,7] <-  gedit(upper[i], width=5, container=normal.flo))
		e$init  <- c(e$init, normal.flo[row,8] <- gedit(ini[i], width=10, container=normal.flo))
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
	widget.defaults[['init']] <- if(is.null(ini)) rep('', 6) else ini
	
	linked.pars.list[['Triangle.ini']] <- e$init[1:4]
	linked.pars.list[['k.ini']] <- e$init[[5]]
	linked.pars.list[['z.ini']] <- e$init[[6]]
	linked.pars.list[['Triangle.ini.low']] <- e$lower[1:4]
	linked.pars.list[['Triangle.ini.up']] <- e$upper[1:4]
	linked.pars.list[['k.ini.low']] <- e$lower[[5]]
	linked.pars.list[['k.ini.up']] <- e$upper[[5]]
	linked.pars.list[['z.ini.low']] <- e$lower[[6]]
	linked.pars.list[['z.ini.up']] <- e$upper[[6]]
	linked.pars.list[['Triangle.prior.low']] <- e$prior.low[1:4]
	linked.pars.list[['Triangle.prior.up']] <- e$prior.up[1:4]
	linked.pars.list[['k.prior.low']] <- e$prior.low[[5]]
	linked.pars.list[['k.prior.up']] <- e$prior.up[[5]]
	linked.pars.list[['z.prior.low']] <- e$prior.low[[6]]
	linked.pars.list[['z.prior.up']] <- e$prior.up[[6]]
		
	ini.normal.f <- gframe("<span  color='#0B6138'>Normal priors for country specific parameters</span>", markup=TRUE, container=normal.g,
								horizontal=FALSE)
	nini.flo <- glayout(container=ini.normal.f)
	
	l <- 1 # row 1
	nini.flo[l,1] <- ''
	nini.flo[l,2] <- '   prior\nlower b.'
	nini.flo[l,3] <- '   prior\nupper b.'
	nini.flo[l,4] <- ' init\nmean'
	nini.flo[l,5] <- 'init\n sd'

	prior.c.low <- c(defaults$Triangle.c.prior.low, defaults$k.c.prior.low, defaults$z.c.prior.low)
	prior.c.up <- c(defaults$Triangle.c.prior.up, defaults$k.c.prior.up, defaults$z.c.prior.up)
	ini.c.means <- c(defaults$Triangle.c.ini.norm[[1]], 
					defaults$k.c.ini.norm[1], defaults$z.c.ini.norm[1])
	ini.c.sd <- c(defaults$Triangle.c.ini.norm[[2]], 
					defaults$k.c.ini.norm[2], defaults$z.c.ini.norm[2])
	e$c.low <- e$c.up <- e$c.means <- e$c.sd <- NULL
	c.labels <- paste('<span>', c(paste('Triangle<sup>c</sup><sub>', 1:4, '</sub>', sep=''), 
						'k<sup>c</sup>', 'z<sup>c</sup>'), ':</span>', sep='')
	for (i in 1:6) {
		row <- l + i
		nini.flo[row,1] <- glabel(c.labels[i], markup=TRUE, container=nini.flo)
		e$c.low <- c(e$c.low, nini.flo[row,2] <- gedit(prior.c.low[i], width=5, container=nini.flo))
		e$c.up <- c(e$c.up, nini.flo[row,3] <- gedit(prior.c.up[i], width=5, container=nini.flo))
		e$c.means <- c(e$c.means, nini.flo[row,4] <- gedit(ini.c.means[i], width=5, container=nini.flo))
		e$c.sd <- c(e$c.sd, nini.flo[row,5] <- gedit(ini.c.sd[i], width=5, container=nini.flo))
	}
	widget.defaults[['c.low']] <- prior.c.low
	widget.defaults[['c.up']] <- prior.c.up
	widget.defaults[['c.means']] <- ini.c.means
	widget.defaults[['c.sd']] <- ini.c.sd

	linked.pars.list[['Triangle.c.prior.low']] <- e$c.low[1:4]
	linked.pars.list[['Triangle.c.prior.up']] <- e$c.up[1:4]
	linked.pars.list[['k.c.prior.low']] <- e$c.low[[5]]
	linked.pars.list[['k.c.prior.up']] <- e$c.up[[5]]
	linked.pars.list[['z.c.prior.low']] <- e$c.low[[6]]
	linked.pars.list[['z.c.prior.up']] <- e$c.up[[6]]
	linked.pars.tuple[['Triangle.c.ini.norm']] <- list(e$c.means[1:4], e$c.sd[1:4])
	linked.pars.tuple[['k.c.ini.norm']] <- c(e$c.means[[5]], e$c.sd[[5]])
	linked.pars.tuple[['z.c.ini.norm']] <- c(e$c.means[[6]], e$c.sd[[6]])

	gamma.uniform.g <- ggroup(horizontal=TRUE, container=priors.f)
	priors.gamma.f <- gframe("<span  color='#0B6138'>Gamma priors</span>", markup=TRUE, container=gamma.uniform.g,
								horizontal=FALSE)
	shape.g <- ggroup(horizontal=TRUE, container=priors.gamma.f)
	glabel('shape: nu', container=shape.g)
	e$nu <- gedit(defaults$nu, width=3, container=shape.g)
	glabel('/2', container=shape.g)
	
	gamma.flo <- glayout(container=priors.gamma.f)
	
	l <- 1 # row 1
	gamma.flo[l,1] <- ''
	gamma.flo[l,2] <- 'sqrt{rate} (tau)'
	gamma.flo[l,3] <- 'center'
	gamma.flo[l,4] <- 'init lower'
	gamma.flo[l,5] <- 'init upper'
	gamma.flo[l,6] <- 'init values'
	
	
	centers <- round(sqrt(2*defaults$tau/defaults$nu),2)
	lambda.lower <- c(defaults$lambda.ini.low, defaults$lambda.k.ini.low, defaults$lambda.z.ini.low)
	lambda.upper <- c(defaults$lambda.ini.up, defaults$lambda.k.ini.up, defaults$lambda.z.ini.up)
	lambda.ini <- c(unlist(defaults$lambda.ini), defaults$lambda.k.ini, defaults$lambda.z.ini)

	g.labels <- paste('<span>', paste('lambda<sub>', c(1:4, 'k', 'z'), '</sub>', sep=''), 
						':</span>', sep='')
	e$tau <- e$center <- NULL
	for (i in 1:6) {
		row <- l + i
		gamma.flo[row,1] <- glabel(g.labels[i], markup=TRUE, container=gamma.flo)
		e$tau <- c(e$tau, gamma.flo[row,2] <- gedit(defaults$tau[i], width=5, container=gamma.flo))
		e$center <- c(e$center, gamma.flo[row,3] <- gedit(centers[i], width=5, container=gamma.flo))
		addHandlerChanged(e$tau[[i]], action=list(idx=i), handler=function(h1,...) {
			blockHandler(e$center[[h1$action$idx]])
			svalue(e$center[[h1$action$idx]]) <- round(sqrt(2*as.numeric(svalue(e$tau[[h1$action$idx]]))/as.numeric(svalue(e$nu))),2)
			unblockHandler(e$center[[h1$action$idx]])
			})
		addHandlerChanged(e$center[[i]], action=list(idx=i), handler=function(h1,...) {
			blockHandler(e$tau[[h1$action$idx]])
			svalue(e$tau[[h1$action$idx]]) <- as.numeric(svalue(e$center[[h1$action$idx]]))^2 * as.numeric(svalue(e$nu))/2
			unblockHandler(e$tau[[h1$action$idx]])
			})
		e$lambda.lower <- c(e$lambda.lower, gamma.flo[row,4] <- gedit(lambda.lower[i], width=7, container=gamma.flo))
		e$lambda.upper <- c(e$lambda.upper, gamma.flo[row,5] <-  gedit(lambda.upper[i], width=7, container=gamma.flo))
		e$lambda.init  <- c(e$lambda.init, gamma.flo[row,6] <- gedit(lambda.ini[i], width=12, container=gamma.flo))
		addHandlerChanged(e$lambda.init[[i]], action=list(idx=i), handler=function(h1,...) {
						isempty <- nchar(svalue(e$lambda.init[[h1$action$idx]]))==0
						enabled(e$lambda.lower[[h1$action$idx]]) <- isempty
						enabled(e$lambda.upper[[h1$action$idx]]) <- isempty
						})

	}
	addHandlerChanged(e$nu, handler=function(h1,...) for(i in 1:6) {
			blockHandler(e$center[[i]])
			svalue(e$center[[i]]) <- round(sqrt(2*as.numeric(svalue(e$tau[[i]]))/as.numeric(svalue(e$nu))),2)
			unblockHandler(e$center[[i]])
			})
	widget.defaults[['lambda.lower']] <- lambda.lower
	widget.defaults[['lambda.upper']] <- lambda.upper
	widget.defaults[['lambda.init']] <- if(is.null(lambda.ini)) rep('', 6) else lambda.ini
	
	linked.pars.list[['lambda.ini']] <- e$lambda.init[1:4]
	linked.pars.list[['lambda.k.ini']] <- e$lambda.init[[5]]
	linked.pars.list[['lambda.z.ini']] <- e$lambda.init[[6]]
	linked.pars.list[['lambda.ini.low']] <- e$lambda.lower[1:4]
	linked.pars.list[['lambda.ini.up']] <- e$lambda.upper[1:4]
	linked.pars.list[['lambda.k.ini.low']] <- e$lambda.lower[[5]]
	linked.pars.list[['lambda.k.ini.up']] <- e$lambda.upper[[5]]
	linked.pars.list[['lambda.z.ini.low']] <- e$lambda.lower[[6]]
	linked.pars.list[['lambda.z.ini.up']] <- e$lambda.upper[[6]]

	info.gamma.glo <- glayout(container=priors.gamma.f)
	info.gamma.glo[1,1, expand=FALSE, anchor=c(-1,0)] <- glabel('"center" is the prior center for sigma: sqrt(2*tau/nu)', 
					container=  info.gamma.glo)
	info.gamma.glo[2,1, expand=FALSE] <- glabel('After editing a cell, hit Enter and wait for the other cells to be updated.', 
				container=  info.gamma.glo)
	#info.gamma.glo[4,1, expand=FALSE] <- glabel('', container=  info.gamma.glo)
	
	uniform.dlf.g <- ggroup(horizontal=FALSE, container=gamma.uniform.g, expand=FALSE)

	
	priors.unif.f <- gframe("<span  color='#0B6138'>Uniform priors</span>", markup=TRUE, container=uniform.dlf.g,
								horizontal=TRUE)
	unif.flo <- glayout(container=priors.unif.f)
	l <- 1
	unif.flo[l,1] <- ''
	unif.flo[l,2] <- 'init lower'
	unif.flo[l,3] <- 'init upper'
	unif.flo[l,4] <- 'init values'
	
	l<-l+1
	unif.flo[l,1] <- 'omega'
	e$omega.ini.low <- unif.flo[l,2] <- gedit(defaults$omega.ini.low, width=5, container=unif.flo)
	e$omega.ini.up <- unif.flo[l,3] <- gedit(defaults$omega.ini.up, width=5, container=unif.flo)
	e$omega.ini <- unif.flo[l,4] <- gedit(defaults$omega.ini, width=10, container=unif.flo)
	addHandlerChanged(e$omega.ini, handler=function(h1,...) {
						isempty <- nchar(svalue(e$omega.ini))==0
						enabled(e$omega.ini.low) <- isempty
						enabled(e$omega.ini.up) <- isempty
						})
	widget.defaults[['omega.ini']] <- if(is.null(defaults$omega.ini)) '' else defaults$omega.ini
	
	addSpace(uniform.dlf.g, 10)
	dlf.f <- gframe("<span  color='blue'>Double-logistic parameters</span>", markup=TRUE, container=uniform.dlf.g,
								horizontal=TRUE)
	glabel('p1:', container=dlf.f)
	e$dl.p1 <- gedit(defaults$dl.p1, width=5, container=dlf.f)
	addSpace(dlf.f, 10)
	glabel('p2:', container=dlf.f)
	e$dl.p2 <- gedit(defaults$dl.p2, width=5, container=dlf.f)

	addSpring(uniform.dlf.g)
	info.ini.group <- ggroup(horizontal=FALSE, container=uniform.dlf.g)
	glabel('NOTE:', container=info.ini.group, anchor=c(-1,0))
	glabel('Leave "init values" blank for the starting values', container=info.ini.group, anchor=c(-1,0))
	glabel('being equally distributed between "init lower" and', container=info.ini.group, anchor=c(-1,0))
	glabel('"init upper". For specific initial values enter', container=info.ini.group, anchor=c(-1,0))
	glabel('one value per chain separated by commas.', container=info.ini.group, anchor=c(-1,0))

	
	# Buttons
	button.g <- ggroup(container=e$adv.g, horizontal=TRUE)
	bDem.gbutton('Cancel', container=button.g, handler=function(h, ...) 
					visible(adv.set.win) <- FALSE)
	addSpring(button.g)
	e$adv.set.defaultbutton <- bDem.gbutton(action=gaction('  Set to Default Values  ', icon='refresh', handler=set.defaults),
										container=button.g)
	e$adv.set.okbutton <- bDem.gbutton('OK', container=button.g)
	
	e$linked.pars.list <- linked.pars.list
	e$linked.pars.tuple <- linked.pars.tuple
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