e0Results.group <- function(g, main.win, parent) {
	e <- new.env()
								
	e$sim.dir <- parent$sim.dir
	graph.defaults <- formals(png)

	nb <- bDem.gnotebook(container=g, expand=TRUE)
	
	traj.g <- ggroup(label="<span color='#0B6138'>e0 trajectories</span>", 
							markup=TRUE, horizontal=FALSE, container=nb)
	traj.env <- e0.show.trajectories.group(traj.g, main.win, e)
	map.g <- ggroup(label="<span color='#0B6138'>e0 world maps</span>", 
							markup=TRUE, horizontal=FALSE, container=nb)
	map.env <- e0.show.map.group(map.g, main.win, e)
	dl.g <- ggroup(label="<span color='#0B6138'>DL curve</span>", 
							markup=TRUE, horizontal=FALSE, container=nb)
	dl.env <- e0.show.dl.group(dl.g, main.win, e)
	traces.g <- ggroup(label="<span color='#0B6138'>Parameter traces</span>", 
							markup=TRUE, horizontal=FALSE, container=nb)
	traces.env <- e0.show.traces.group(traces.g, main.win, e)
	convergence.g <- ggroup(label="<span color='#0B6138'>Convergence</span>", markup=TRUE, horizontal=FALSE, container=nb)
	create.convergence.tab(convergence.g, e$sim.dir, type='e0', package='bayesLife', main.win)

	svalue(nb) <- 1
}

e0.show.trajectories.group <- function(g, main.win, parent.env) {
	leftcenter <- c(-1,0)
	e <- new.env()
	e$sim.dir <- parent.env$sim.dir
	e$pred.type <- 'e0'
	defaults.pred <- formals(e0.predict)
	defaults.traj <- formals(e0.trajectories.plot)
	defaults.traj.all <- formals(e0.trajectories.plot.all)
	addSpace(g, 10)	
	country.f <- gframe("<span color='blue'>Country settings</span>", markup=TRUE, 
									horizontal=FALSE, container=g)
	e$show.traj.country <- create.country.widget(country.f, defaults.traj.all, 
									main.win, prediction=TRUE, parent.env=e)
		
	addSpace(g, 10)
	traj.settings.f <- gframe("<span color='blue'>Trajectories settings</span>", markup=TRUE, 
								horizontal=TRUE, container=g)
	lo <- glayout(container=traj.settings.f) 
	lo[1,1, anchor=leftcenter] <- glabel('CI (%):', container=lo)
	lo[1,2] <- e$pi <- gedit('80, 95', width=7, container=lo)
	lo[2,1, anchor=leftcenter] <- glabel('# trajectories:', container=lo)
	lo[2,2] <- e$nr.traj <- gedit(20, width=6, container=lo)
	lo[1,3, anchor=leftcenter] <- 	glabel('From year:', container=lo)
	lo[1,4] <- e$start.year <- gedit('', width=4, container=lo)
	lo[2,3, anchor=leftcenter] <- glabel('To year:', container=lo)
	lo[2,4] <- e$end.year <- gedit('', width=4, container=lo)
	lo[1,5] <- glabel('     ', container=lo)
	lo[1,6] <- e$typical.trajectory <- gcheckbox('Typical trajectory', checked=defaults.traj$typical.trajectory, 
								container=lo)
	lo[3,1, anchor=leftcenter] <- "Type:"
	lo[3,2] <- e$sex <- bDem.gdroplist(c('Female', 'Male', 'Average', 'Both Marginal', 'Both Joint', 'Gap'), container=lo, selected=1,
				handler=function(h,...) {
					if(svalue(h$obj) == 'Both Marginal') {svalue(e$nr.traj) <- 0; svalue(e$pi) <- 95}
					if(svalue(h$obj) == 'Both Joint') {svalue(e$nr.traj) <- 500; svalue(e$pi) <- 95}
					if(svalue(h$obj) == 'Gap') {svalue(e$nr.traj) <- 0; svalue(e$pi) <- '80, 95'}
					if(is.element(svalue(h$obj), c('Female', 'Male', 'Average'))) {svalue(e$nr.traj) <- 20; svalue(e$pi) <- '80, 95'}
					enabled(e$TableB.show.traj) <- is.element(svalue(h$obj), c('Female', 'Male', 'Average'))
					enabled(e$years) <- svalue(h$obj) == 'Both Joint'
					enabled(e$typical.trajectory) <- is.element(svalue(h$obj), c('Female', 'Male', 'Average', 'Both Marginal'))
				})
	lo[3,3, anchor=leftcenter] <- glabel('Years:', container=lo)
	lo[3,4] <- e$years <- gedit('2013, 2048, 2098', width=10, container=lo)
	enabled(e$years) <- FALSE

	addSpace(g, 10)
	graph.f <- gframe("<span color='blue'>Advanced graph parameters</span>", markup=TRUE, 
									horizontal=FALSE, container=g)
	e$graph.pars <- create.graph.pars.widgets(graph.f, main.win=main.win)
	addSpring(g)
	button.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic=c('e0.trajectories.plot', 'e0.joint.plot', 'e0.gap.plot'), package='bayesLife', parent.group=button.g,
						parent.window=main.win)
	addSpring(button.g)
	create.generate.script.button(handler=show.e0.traj, action=list(mw=main.win, env=e, type='plot', script=TRUE),
								container=button.g)
	addSpace(button.g, 5)
	TableB.show.traj.act <- gaction(label='Table', icon='dataframe', handler=show.e0.traj, 
						action=list(mw=main.win, env=e, type='table', script=FALSE))
	GraphB.show.traj.act <- gaction(label='Graph', icon='lines', handler=show.e0.traj, 
						action=list(mw=main.win, env=e, type='plot', script=FALSE))
	e$TableB.show.traj <- bDem.gbutton(action=TableB.show.traj.act, container=button.g)
	bDem.gbutton(action=GraphB.show.traj.act, container=button.g)
	return(e)
}

get.additional.e0.param <- function(e, ...) {
	sex <- svalue(e$sex)
	param <- list(both.sexes= if(sex=='Average') "'A'" else sex=='Both Marginal', joint.male= sex == 'Male')
	return(list(add=param, plot=c('pi', 'xlim', 'nr.traj', 'both.sexes', 'typical.trajectory'), 
					pred='joint.male', pred.unquote=param['joint.male'],
					table=c('pi', 'country', 'both.sexes'), table.decimal=2))
}
	

assemble.e0.plot.cmd <- function(param, e, all=FALSE) {
	plot.type <- svalue(e$sex)
	all.suffix <- if(all) '.all' else ''
	if(is.element(plot.type, c('Female', 'Male', 'Average', 'Both Marginal'))) 
		return(paste('e0.trajectories.plot',all.suffix, '(pred,',
				assemble.arguments(param, svalue(e$graph.pars)), ')', sep=''))
	
	if(plot.type == 'Both Joint') {
		parlist <- list(country=param$country, years=as.numeric(strsplit(svalue(e$years), ',')[[1]]), 
						nr.points=param$nr.traj, pi=param$pi)
		return(paste('e0.joint.plot',all.suffix, '(pred,', assemble.arguments(parlist, svalue(e$graph.pars)), ')', sep=''))
	}
	# gap plot
	parlist <- list(country=param$country, nr.traj=param$nr.traj, pi=param$pi, xlim=param$xlim)
	return(paste('e0.gap.plot',all.suffix, '(pred,', assemble.arguments(parlist, svalue(e$graph.pars)), ')', sep=''))
}

e0.get.trajectories.table.values <- function(pred, param, ...) {
	t <- do.call('e0.trajectories.table', c(list(pred), param))
	# change the column names, otherwise gtable will prefix an 'X'
	colnames(t)[2:ncol(t)] <- paste('q', colnames(t)[2:ncol(t)], sep='')
	return(t)
}


show.e0.traj <- function(h, ...) {
	e <- h$action$env
	pred.type <- if(is.null(h$action$pred.type)) 'e0' else h$action$pred.type
	package <- if(is.null(h$action$package)) 'bayesLife' else h$action$package
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	show.type <- h$action$type
	allow.null.country <- if(is.null(h$action$allow.null.country)) FALSE else h$action$allow.null.country
	country.pars <- NULL
	if(is.null(h$action$env$set.country.to.null) || !h$action$env$set.country.to.null)
		country.pars <- get.country.code.from.widget(e$show.traj.country$country.w, e$show.traj.country, 
							force.country.spec = show.type!='plot', allow.null.country=allow.null.country)
	if(is.null(country.pars)) {
		if(!allow.null.country) return(NULL)
		else if(!do.call(paste('.', pred.type, '.traj.country.check', sep=''), list(e)))
				 return(NULL)
	}
	draw.one.country <- is.null(country.pars$output.dir)
	param.names.all <- list(text='sim.dir', numvector='pi', 
							numeric=c('nr.traj', 'start.year', 'end.year'),
							logical='typical.trajectory')
	param.env <- get.parameters(param.names.all, env=e, quote=h$action$script)
	param.env.rest <- list(country=country.pars$code, output.dir=country.pars$output.dir,
							output.type=country.pars$output.type, verbose=TRUE)
	param.env <- c(param.env, 
					get.parameters(list(text=c('output.dir', 'output.type'), 
										logical='verbose', numeric='country'), 
									param.env.rest, quote=TRUE,
									retrieve.from.widgets=FALSE))
	add.param.names <- do.call(paste('get.additional.',pred.type, '.param', sep=''), 
									list(e, script=h$action$script, type=show.type))
	param.env <- c(param.env, add.param.names[['add']])
	param.pred <- param.env[c('sim.dir', 
					add.param.names[['pred']][is.element(add.param.names[['pred']], names(param.env))])]
	#param.pred.ev <- param.pred
	 # get it now unquoted (to avoid double quotes if script is TRUE)
	param.pred.ev <- c(get.parameters(list(text='sim.dir'), env=e, quote=FALSE),
						add.param.names[['pred.unquote']])
	pred <- do.call(paste('get.', pred.type, '.prediction', sep=''), param.pred.ev)
	if(is.null(pred)) {
		gmessage('Simulation directory contains no prediction of the specified type.', 
					title='Input Error', icon='error')
		return(NULL)
	}
	if(h$action$script) {
		cmd <- paste('pred <- get.', pred.type, '.prediction(', 
					assemble.arguments(param.pred), ')\n', sep='')
	} else {	
		cmd <- ''
	}
	if(!is.null(add.param.names$delete))
		for(par in add.param.names$delete) param.env[[par]] <- NULL
	xmin <- param.env$start.year
	xmax <- param.env$end.year
	if(!is.null(xmin) || !is.null(xmax)) {
		param.env.xlim <- list(xlim=paste(if(!is.null(xmin)) xmin else pred$mcmc.set$meta$start.year,
									 if(!is.null(xmax)) xmax else pred$end.year, sep=', '))
		param.env <- c(param.env, get.parameters(list(numvector='xlim'), param.env.xlim, quote=h$action$script,
									retrieve.from.widgets=FALSE))
	}

	pars.value <- svalue(e$graph.pars)
	if (show.type == 'plot') {
		param.plot1c <- param.env[add.param.names[['plot']][is.element(add.param.names[['plot']], names(param.env))]]
		if(is.element('country', names(param.env))) param.plot1c <- c(param.plot1c, param.env['country'])
		if(draw.one.country) { # one country
			cmd <- paste(cmd, do.call(paste('assemble.', pred.type, '.plot.cmd', sep=''), list(param.plot1c, e)), sep='')
			if (h$action$script) {
				create.script.widget(cmd, h$action$mw, package=package)
			} else {
				create.graphics.window(parent=h$action$mw, title=paste("Trajectories", 
							if(!is.null(country.pars$name)) paste("for", country.pars$name) else ""))
				eval(parse(text=cmd))
			}
		} else { # all countries
			param.plot.allc <- param.env[c(names(param.plot1c), 'output.dir', 'output.type',  'verbose')]
			cmd <- paste(cmd, do.call(paste('assemble.', pred.type, '.plot.cmd', sep=''), list(param.plot.allc, e, all=TRUE)), sep='')
			# cmd <- paste(cmd, paste(pred.type, '.trajectories.plot', 
						# if(pred.type=='pop') 'All' else '.all',
						# '(pred, ', sep=''), 
					# paste(paste(names(param.plot.allc), param.plot.allc, sep='='), collapse=', '), sep='')
			# if(!is.null(pars.value)) {
				# if(nchar(pars.value)>0)
					# cmd <- paste(cmd, ',', pars.value)
			# }
			# cmd <- paste(cmd, ')', sep='')
			if (h$action$script) {
				create.script.widget(cmd, h$action$mw, package=package)
			} else {
				eval(parse(text=cmd))
			}
		}
	} else {
		# Table
		param.table <- param.env[add.param.names[['table']][is.element(add.param.names[['table']], names(param.env))]]
		table.values <- do.call(paste(pred.type, '.get.trajectories.table.values', sep=''), 
							list(pred, param.table, e))
		table.values <- round(table.values[!apply(is.na(table.values), 1, all),],
								add.param.names[['table.decimal']])
		table.values <- cbind(rownames(table.values), table.values)
		colnames(table.values)[1] <- 'year'
		win <- bDem.gwindow(do.call(paste('get.', pred.type, '.table.title', sep=''), 
						list(country.pars$name, pred, e)),
					parent=h$action$mw, height=max(min(22.2*(dim(table.values)[1]+1),600), 100))
		g <- ggroup(container=win, horizontal=FALSE, expand=TRUE)
		gt <- gtable(table.values, container=g, expand=TRUE)
		bDem.gbutton('Print to R Console', container=g, handler=function(h,...){
										print(do.call(paste(pred.type, '.trajectories.table', sep=''), 
												c(list(pred), param.table)))})
	}
}

get.e0.table.title <- function(country, pred, ...) 
	return (paste(country, '-', bayesLife:::get.sex.label(pred$mcmc.set$meta)))
	
e0.show.map.group <- function(g, main.win, parent.env) {
	e <- new.env()
	e$sim.dir <- parent.env$sim.dir
	addSpace(g, 10)
	lo <- .create.map.settings.group(g, e, measures=c('e0', e0.parameter.names.cs.extended()))
	lo[4, 1, anchor=c(-1,0)] <- "Sex:"
	lo[4, 2] <- e$sex <- bDem.gdroplist(c('Female', 'Male'), container=lo, selected=1)
	addSpring(g)
	bg <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='e0.map', package='bayesLife', parent.group=bg,
						parent.window=main.win)
	addSpring(bg)
	create.generate.script.button(handler=e0.showMap, action=list(mw=main.win, env=e, script=TRUE),
								container=bg)
	addSpace(bg, 5)
	GraphB.map <- gaction(label=' Show Map ', handler=e0.showMap, 
						action=list(mw=main.win, env=e, script=FALSE))
	bDem.gbutton(action=GraphB.map, container=bg)
}

e0.showMap <- function(h, ...) {
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
	sex <- svalue(e$sex)
	param.pred$joint.male <- sex == 'Male'
	map.function <- if(package == 'rworldmap') 'e0.map' else 'e0.map.gvis'
	if(h$action$script) {
		cmd <- paste('pred <- get.e0.prediction(', assemble.arguments(param.pred), ')\n', sep='')
		if (par.name == 'e0') {
			if(package == 'rworldmap') {
				cmd <- paste(cmd, "param.map <- get.e0.map.parameters(pred, same.scale=", same.scale,
					", quantile=", quantile, ")\n", sep="")
				cmd <- paste(cmd, 'do.call("', map.function, '", param.map)', sep='')
			} else {
				cmd <- paste(cmd, map.function, '(pred, quantile=', quantile, ', pi=', bounds, ')', sep='')
			}
		} else {
			cmd <- paste(cmd, map.function, '(pred, quantile=', quantile, ', par.name="', par.name, '"', sep='')
			cmd <- paste(cmd, if (package == 'googleVis') paste(', pi=', bounds, sep='') else '', sep='')
			cmd <- paste(cmd, ')', sep='')
		}
		create.script.widget(cmd, h$action$mw, package="bayesLife")
	} else {
		pred <- do.call('get.e0.prediction', param.pred)
		if (par.name == 'e0' && package == 'rworldmap') {
			param.map <- get.e0.map.parameters(pred, same.scale=same.scale, quantile=quantile)
		} else {
			param.map <- list(pred=pred, quantile=quantile)
			if (par.name != 'e0')
				param.map[['par.name']]<- par.name
		}
		if(package == 'rworldmap') param.map[['device']] <- 'dev.new'
		if (package == 'googleVis') param.map[['pi']] <- bounds
		g <- create.graphics.map.window(parent=h$action$mw, pred=pred, params=param.map, percentile=percentile, 
										is.gvis= package == 'googleVis', title="World Map", type='e0', 
										cw.main=paste(c(par.name, sex, percentile), collapse=', '))
	}
}

e0.show.dl.group <- function(g, main.win, parent.env) {
	e <- new.env()
	leftcenter <- c(-1,0)
	e$sim.dir <- parent.env$sim.dir
	e$pred.type <- 'e0'
	defaults.dl <- formals(e0.DLcurve.plot)
	defaults.dl.all <- formals(e0.DLcurve.plot.all)
	addSpace(g, 10)
	country.f <- gframe("<span color='blue'>Country settings</span>", markup=TRUE, 
							horizontal=FALSE, container=g)
	e$dlc.country <- create.country.widget(country.f, defaults.dl.all, main.win, prediction=FALSE, 
											parent.env=e, disable.table.button=FALSE)
	addSpace(g, 10)
	dl.f <- gframe("<span color='blue'>DL curve settings</span>", markup=TRUE, 
							horizontal=FALSE, container=g)
	dllo <- glayout(horizontal=TRUE, container=dl.f)
	dllo[1,1, anchor=leftcenter] <- glabel('CI (%):', container=dllo)
	dllo[1,2] <- e$pi <- gedit('80, 95', width=7, container=dllo)
	dllo[1,3, anchor=leftcenter] <- glabel('Burnin:', container=dllo)
	dllo[1,4] <- e$burnin <- gedit(defaults.dl$burnin, width=5, container=dllo)
	dllo[2,3, anchor=leftcenter] <- glabel('e0 min, max:', container=dllo)
	dllo[2,4] <- e$e0.lim <- gedit(defaults.dl$tfr.max, width=2, container=dllo)
	dllo[2,1, anchor=leftcenter] <- glabel('# curves:', container=dllo)
	dllo[2,2] <- e$nr.curves <- gedit(20, width=6, container=dllo)
	dllo[1,5] <- e$predictive.distr <- gcheckbox('Predictive distribution', 
							checked=defaults.dl$predictive.distr, container=dllo)
	addSpace(g, 10)
	graph.f <- gframe("<span color='blue'>Advanced graph parameters</span>", markup=TRUE, 
						horizontal=FALSE, container=g)
	e$graph.pars <- create.graph.pars.widgets(graph.f, main.win=main.win)
	addSpring(g)
	button.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='e0.DLcurve.plot', package='bayesLife', parent.group=button.g,
						parent.window=main.win)
	addSpring(button.g)
	create.generate.script.button(handler=e0.showDLcurve, action=list(mw=main.win, env=e, script=TRUE),
								container=button.g)
	addSpace(button.g, 5)
	GraphB.dlc <- gaction(label='Graph', icon='lines', handler=e0.showDLcurve, 
						action=list(mw=main.win, env=e, script=FALSE))
	bDem.gbutton(action=GraphB.dlc, container=button.g)
}

e0.showDLcurve <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	country.pars <- get.country.code.from.widget(e$dlc.country$country.w, e$dlc.country)
	if(is.null(country.pars)) return(NULL)
	param.names.all <- list(text='sim.dir', numvector=c('pi', 'e0.lim'),
							numeric=c('nr.curves', 'burnin'),
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
		cmd <- paste('m <- get.e0.mcmc(', assemble.arguments(param.mcmc), ')\n', sep='')
	} else {
		m <- do.call('get.e0.mcmc', param.mcmc)
		cmd <- ''
	}
	pars.value <- svalue(e$graph.pars)
	param.plot1c <- list()
	for (par in c('pi', 'nr.curves', 'e0.lim', 'country', 'burnin', 'predictive.distr')) 
		if(is.element(par, names(param.env))) param.plot1c <- c(param.plot1c, param.env[par])

	if(!is.null(country.pars$code)) { # one country
		cmd <- paste(cmd, 'e0.DLcurve.plot(mcmc.list=m, ', assemble.arguments(param.plot1c, pars.value), ')', sep='')
		if (h$action$script) {
			create.script.widget(cmd, h$action$mw, package="bayesLife")
		} else {
			create.graphics.window(parent=h$action$mw, title=paste("Double Logistic Curves for", country.pars$name))
			eval(parse(text=cmd))
		}
	} else { # all countries
		param.plot.allc <- param.env[c(names(param.plot1c), 'output.dir', 'output.type',  'verbose')]
		cmd <- paste(cmd, 'e0.DLcurve.plot.all(mcmc.list=m, ', 
						assemble.arguments(param.plot.allc, pars.value), ')', sep='')
		if (h$action$script) {
			create.script.widget(cmd, h$action$mw, package="bayesLife")
		} else {
			eval(parse(text=cmd))
		}
	}
}

e0.show.traces.group <- function(g, main.win, parent.env) {
	e <- new.env()
	e$sim.dir <- parent.env$sim.dir
	e$pred.type <- 'e0'
	.create.partraces.settings.group(g, e, main.win, par.names=e0.parameter.names(), par.names.cs=e0.parameter.names.cs())
	addSpring(g)
	button.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='e0.partraces.plot', package='bayesLife', parent.group=button.g,
						parent.window=main.win)	
	addSpring(button.g)
	SummaryB.traces <- gaction(label='Show summary', icon='dataframe', handler=e0.showParTraces, 
						action=list(mw=main.win, env=e, print.summary=TRUE))
	bDem.gbutton(action=SummaryB.traces, container=button.g)
	GraphB.traces <- gaction(label='Graph', icon='lines', handler=e0.showParTraces, 
						action=list(mw=main.win, env=e, print.summary=FALSE))
	bDem.gbutton(action=GraphB.traces, container=button.g)
}

e0.showParTraces <- function(h, ...) {
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
	if(print.summary) {
		mc.summary <- c()
		warn <- getOption('warn')
		options(warn=-1) # disable warning messages
		mcmc.set <- get.e0.mcmc(params[['sim.dir']])
		options(warn=warn)
		con <- textConnection("mc.summary", "w", local=TRUE)
		mc.exist <- TRUE
		sink(con)
		if (is.null(mcmc.set)) {
			cat('No simulation available in this directory.')
			mc.exist <- FALSE
		}
	} else create.graphics.window(parent=h$action$mw, title="Parameter traces", dpi=150, height=4*150)
	
	if (cs==2) { # country-specific parameters
		if (!all.pars) {
			pars <- svalue(e$par.cs.dl)
			if(print.summary) {if (mc.exist) print(summary(mcmc.set, country=country.pars$code, par.names.cs=pars, par.names=NULL, 
											burnin=params[['burnin']], thin=params[['thin']]))
			} else e0.partraces.cs.plot(sim.dir=params[['sim.dir']], country=country.pars$code, par.names=pars, 
											nr.points=params[['nr.points']], 
											burnin=params[['burnin']], thin=params[['thin']])
		} else {
			if(print.summary){if (mc.exist) print(summary(mcmc.set, country=country.pars$code, par.names=NULL, 
											burnin=params[['burnin']], thin=params[['thin']]))
			} else e0.partraces.cs.plot(sim.dir=params[['sim.dir']], country=country.pars$code, nr.points=params[['nr.points']], 
											burnin=params[['burnin']], thin=params[['thin']])
		}
	} else { # World-parameters
		if (!all.pars) { # selected pars
			pars <- svalue(e$par.dl)
			if(print.summary) {if (mc.exist) print(summary(mcmc.set, par.names.cs=NULL, par.names=pars, 
											burnin=params[['burnin']], thin=params[['thin']]))
			} else e0.partraces.plot(sim.dir=params[['sim.dir']], par.names=pars, nr.points=params[['nr.points']], 
											burnin=params[['burnin']], thin=params[['thin']])
		} else { # all pars
			if(print.summary) {if (mc.exist) print(summary(mcmc.set, par.names.cs=NULL, 
											burnin=params[['burnin']], thin=params[['thin']]))
			} else e0.partraces.plot(sim.dir=params[['sim.dir']], nr.points=params[['nr.points']], 
											burnin=params[['burnin']], thin=params[['thin']])
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


