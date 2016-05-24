popResults.group <- function(g, main.win, parent) {
	e <- new.env()
								
	e$sim.dir <- parent$sim.dir
	
	graph.defaults <- formals(png)
	
	gg <-  ggroup(horizontal=FALSE, container=g, expand=TRUE)
	aggr.g <- ggroup(horizontal=TRUE, container=gg)
	addSpring(aggr.g)
	nb <- bDem.gnotebook(container=gg, expand=TRUE)
	traj.g <- ggroup(label="<span color='#0B6138'>Population trajectories</span>", 
							markup=TRUE, horizontal=FALSE, container=nb)
	traj.env <- pop.show.trajectories.group(traj.g, main.win, e)
	cohort.g <- ggroup(label="<span color='#0B6138'>Cohort trajectories</span>", 
							markup=TRUE, horizontal=FALSE, container=nb)
	cohort.env <- pop.show.cohort.trajectories.group(cohort.g, main.win, e)
	pyr.g <- ggroup(label="<span color='#0B6138'>Population pyramid</span>", 
							markup=TRUE, horizontal=FALSE, container=nb)
	pyr.env <- pop.show.pyramid.group(pyr.g, main.win, e)
	map.g <- ggroup(label="<span color='#0B6138'>Pop world maps</span>", 
							markup=TRUE, horizontal=FALSE, container=nb)
	map.env <- pop.show.map.group(map.g, main.win, e)

	svalue(nb) <- 1
}

pop.show.trajectories.group <- function(g, main.win, parent.env) {
	leftcenter <- c(-1,0)
	e <- new.env()
	e$sim.dir <- parent.env$sim.dir
	e$pred.type <- 'pop'
	e$new.country.select.if.changed <- c('aggregation.dl')
	defaults.pred <- formals(pop.predict)
	defaults.traj <- formals(pop.trajectories.plot)
	defaults.traj.all <- formals(pop.trajectories.plotAll)
	addSpace(g, 10)	
	country.f <- gframe("<span color='blue'>Country settings</span>", markup=TRUE, 
									horizontal=FALSE, container=g)
	glo <- glayout(container=country.f)
	glo[1,1, anchor=leftcenter] <- 'Aggregation:'
	glo[1,2] <- e$aggregation.dl <- bDem.gbutton(" None ", container=glo,
				handler=selectAggregationMenuPop,
				action=list(mw=main.win, env=e, label.widget.name='aggregation.dl'))
	e$aggregation <- NULL
	e$show.traj.country <- create.country.widget(country.f, defaults.traj.all, 
									main.win, glo=glo, start.row=2, prediction=TRUE, parent.env=e)
	l <- 2+4
	exp.g <- e$show.traj.country$country.lo
	exp.g[l,1:3] <- gseparator(container=exp.g)
	exp.g[l+1,1] <- e$use.expression <- gcheckbox('Expression: ', container=exp.g, checked=FALSE, 
						handler=function(h,...){
							use.expr <- svalue(h$obj)
							enabled(e$expression) <- use.expr
							enabled(e$show.traj.country$country.w) <- !use.expr
							enabled(e$show.traj.country$country.select.b) <- !use.expr
							enabled(e$sex) <- !use.expr
							enabled(e$sum.over.ages) <- !use.expr
							enabled(e$age.gb) <- !use.expr
							enabled(e$half.child.variant) <- !use.expr
							e$set.country.to.null <- use.expr
							})
	exp.g[l+1,2:3] <- e$expression <- gedit('', container=exp.g, fill=TRUE, expand=TRUE)
	tooltip(e$expression) <- "See ?pop.expressions. Use XXX as country code if 'All countries' is checked."
	enabled(e$expression) <- FALSE
	addSpace(g, 10)
	
	lo <- .create.trajectories.settings.group(g, e, 
				defaults=list(start.year=defaults.pred$start.year, end.year=defaults.pred$end.year,
							half.child.variant=defaults.traj$half.child.variant,
							typical.trajectory=defaults.traj$typical.trajectory), l=2)	
	lo[1,1, anchor=leftcenter] <- glabel('Sex:', container=lo)
	lo[1,2] <- e$sex <- bDem.gdroplist(c('both', 'female', 'male'), container=lo, selected=1)
	lo[1,3:4] <- e$sum.over.ages <- gcheckbox("Sum over ages", container=lo, checked=TRUE,
							handler=function(h,...){
								enabled(e$TableB.show.traj) <- svalue(h$obj)})
	lo[1,6] <- e$age.gb <- bDem.gbutton(" Age ", container=lo,
				handler=selectAgeMenuPop,
				action=list(mw=main.win, env=e, multiple=TRUE, label.widget.name='age.label'))
	lo[4,1:2] <- e$plot.by.age <- gcheckbox("Plot by age", container=lo, checked=FALSE,
							handler=function(h,...){
								by.age <- svalue(h$obj)
								enabled(e$sum.over.ages) <- !by.age && !svalue(e$use.expression)
								enabled(e$age.gb) <- !by.age && !svalue(e$use.expression)
								enabled(e$start.year) <- !by.age
								enabled(e$end.year) <- !by.age
								enabled(e$year) <- by.age
								if(by.age) enabled(e$TableB.show.traj) <- TRUE
								else enabled(e$TableB.show.traj) <- svalue(e$sum.over.ages)
								})
	lo[4,3, anchor=leftcenter] <- glabel('Year:', container=lo)
	lo[4,4] <- e$year <- gedit('', width=4, container=lo)
	enabled(e$year) <- FALSE
	e$age <- 'all'
	lo[1,7, anchor=leftcenter] <- e$age.label <- glabel('', container=lo)
	addSpace(g, 10)
	graph.f <- gframe("<span color='blue'>Advanced graph parameters</span>", markup=TRUE, 
									horizontal=FALSE, container=g)
	e$graph.pars <- create.graph.pars.widgets(graph.f, main.win=main.win)
	addSpring(g)
	button.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic=c('pop.trajectories.plot', 'pop.expressions'), package='bayesPop', parent.group=button.g,
						parent.window=main.win)
	addSpring(button.g)
	create.generate.script.button(handler=show.e0.traj, 
				action=list(mw=main.win, env=e, type='plot', script=TRUE,
							pred.type='pop', package='bayesPop', allow.null.country=TRUE),
								container=button.g)
	addSpace(button.g, 5)
	TableB.show.traj.act <- gaction(label='Table', icon='dataframe', handler=show.e0.traj, 
						action=list(mw=main.win, env=e, type='table', script=FALSE,
											pred.type='pop', package='bayesPop', allow.null.country=TRUE))
	GraphB.show.traj.act <- gaction(label='Graph', icon='lines', handler=show.e0.traj, 
						action=list(mw=main.win, env=e, type='plot', script=FALSE,
											pred.type='pop', package='bayesPop', allow.null.country=TRUE))
	e$TableB.show.traj <- bDem.gbutton(action=TableB.show.traj.act, container=button.g)
	bDem.gbutton(action=GraphB.show.traj.act, container=button.g)
	enabled(e$TableB.show.traj) <- svalue(e$sum.over.ages)
	return(e)

}


get.additional.pop.param <- function(e, script, type, ...) {
	param.names <- list(text=c('sex', 'expression'), logical=c('sum.over.ages', 'half.child.variant'), numeric='year')
	non.widget.pars <- list(text='aggregation')
	if (is.element(0, e$age) || e$age=='all') {
		non.widget.pars$text <- c(non.widget.pars$text, 'age')
		e$age <- 'all'
	}
	quote <- if(type=='table') script else TRUE
	
	param <- c(get.parameters(param.names, env=e, quote=quote), 
			get.parameters(non.widget.pars,
					env=e, quote=quote, retrieve.from.widgets=FALSE))
	if(!is.element('age', names(param))) param[['age']] <- e$age # if age is not all it is a vector
	delete.pars <- c()
	if(!svalue(e$use.expression)) delete.pars <- c(delete.pars, 'expression')
	else  # disable if using expression
		delete.pars <- c(delete.pars, c('sex', 'sum.over.ages', 'half.child.variant', 'age'))
	if(svalue(e$plot.by.age)) delete.pars <- c(delete.pars, c('sum.over.ages', 'age', 'start.year', 'end.year'))
	else delete.pars <- c(delete.pars, 'year')
	return(list(add=param, 
			plot=c('pi', 'xlim', 'nr.traj', 'sex', 'age', 'sum.over.ages', 'half.child.variant', 
						'typical.trajectory', 'expression', 'year'), 
			 table=c('pi', 'country', 'sex', 'age', 'half.child.variant', 'expression', 'year'),
			 pred=c('aggregation'),
			 delete=delete.pars,
			 pred.unquote=get.parameters(list(text='aggregation'), 
					 				env=e, quote=FALSE, retrieve.from.widgets=FALSE),
			table.decimal=2))
}

.pop.traj.country.check <- function(e) {
	if(nchar(svalue(e$expression))<=0) {
		gmessage('Country or expression must be specified.', title='Input Error',
					icon='error')
		return(FALSE)
	}
	return(TRUE)
}

assemble.pop.plot.cmd <- function(param, e, all=FALSE) {
	all.suffix <- if(all) 'All' else ''
	cmd <- if(svalue(e$plot.by.age)) 'pop.byage.plot' else 'pop.trajectories.plot'
	return(paste(paste(cmd, all.suffix, sep=''), '(pred,',
				assemble.arguments(param, svalue(e$graph.pars)), ')'))
}

get.pop.table.title <- function(country, pred, e) {
	title <- country
	if(svalue(e$sex) != 'both') title <- paste(title, svalue(e$sex))
	if(nchar(svalue(e$age.label)) > 0) title <- paste(title, 'age:', svalue(e$age.label))
	return (title)
}

pop.get.trajectories.table.values <- function(pred, param, e, ...) {
	cmd <- if(svalue(e$plot.by.age)) 'pop.byage.table' else 'pop.trajectories.table'
	t <- do.call(cmd, c(list(pred), param))
	# change the column names, otherwise gtable will prefix an 'X'
	colnames(t)[2:ncol(t)] <- paste('q', colnames(t)[2:ncol(t)], sep='')
	return(t)
}

pop.show.cohort.trajectories.group <- function(g, main.win, parent.env) {
	leftcenter <- c(-1, 0)
	e <- new.env()
	e$sim.dir <- parent.env$sim.dir
	e$pred.type <- 'pop'
	e$new.country.select.if.changed <- c('aggregation.dl')
	defaults.pred <- formals(pop.predict)
	defaults.cohort <- formals(bayesPop::pop.cohorts.plot)
		
	addSpace(g, 10)
	country.f <- gframe("<span color='blue'>Country settings</span>", markup=TRUE, 
									horizontal=FALSE, container=g)
	glo <- glayout(container=country.f)
	glo[1,1, anchor=leftcenter] <- 'Aggregation:'
	glo[1,2] <- e$aggregation.dl <- bDem.gbutton(" None ", container=glo,
				handler=selectAggregationMenuPop,
				action=list(mw=main.win, env=e, label.widget.name='aggregation.dl'))
	e$show.cohort.country <- create.country.widget(country.f, defaults.cohort, show.all=FALSE,
									main.win, glo=glo, start.row=2, prediction=TRUE, parent.env=e, disable.table.button=FALSE)
	addSpace(g, 10)
	l <- 3
	exp.g <- e$show.cohort.country$country.lo
	exp.g[l,1:3] <- gseparator(container=exp.g)
	exp.g[l+1,1] <- e$use.expression <- gcheckbox('Expression: ', container=exp.g, checked=FALSE, 
						handler=function(h,...){
							use.expr <- svalue(h$obj)
							enabled(e$expression) <- use.expr
							enabled(e$show.cohort.country$country.w) <- !use.expr
							enabled(e$show.cohort.country$country.select.b) <- !use.expr
							e$set.country.to.null <- use.expr
							})
	exp.g[l+1,2:3] <- e$expression <- gedit('', container=exp.g, fill=TRUE, expand=TRUE)
	tooltip(e$expression) <- "See ?pop.expressions. Use XXX as country code if 'All countries' is checked."
	enabled(e$expression) <- FALSE
	addSpace(g, 10)

	type.settings.f <- gframe("<span color='blue'>Cohort trajectories settings</span>", markup=TRUE, 
								horizontal=FALSE, container=g)
	type.g1 <- glayout(horizontal=TRUE, container=type.settings.f)
	type.g1[1,1, anchor=leftcenter] <- glabel('CI (%):', container=type.g1)
	type.g1[1,2] <- e$pi <- gedit('80, 95', width=7, container=type.g1)
	type.g1[2,1] <- cohorts.gb <- bDem.gbutton(" Cohorts ", container=type.g1,
				handler=selectYearsMenuPop,
				action=list(mw=main.win, env=e, widget='cohorts', multiple=TRUE))
	type.g1[2,2] <- e$cohorts <- gedit('', width=10, container=type.g1)
	tooltip(e$cohorts) <- "There will be one plot per cohort."
	type.g1[1,3, anchor=leftcenter] <- glabel('Legend position:', container=type.g1)
	type.g1[1,4] <- e$legend.pos <- gedit(defaults.cohort$legend.pos, width=12, container=type.g1)
	tooltip(e$legend.pos) <- 'Position of the legend passed to the "legend" function.'
	type.g1[2,3, anchor=leftcenter] <- glabel('# dev columns:', container=type.g1)
	type.g1[2,4] <- e$dev.ncol <- gedit(defaults.cohort$dev.ncol, width=3, container=type.g1)
	tooltip(e$dev.ncol) <- "Number of columns in the plot matrix."
	type.g1[1,5] <- e$show.legend <- gcheckbox('Show legend', container=type.g1, checked=defaults.cohort$show.legend)


	addSpace(g, 10)
	graph.f <- gframe("<span color='blue'>Advanced graph parameters</span>", markup=TRUE, 
									horizontal=FALSE, container=g)
	e$graph.pars <- create.graph.pars.widgets(graph.f, main.win=main.win)
	addSpring(g)
	button.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic=c('pop.cohorts.plot', 'pop.expressions'), package='bayesPop', parent.group=button.g,
						parent.window=main.win)
	addSpring(button.g)
	create.generate.script.button(handler=show.pop.cohorts, 
							action=list(mw=main.win, env=e, script=TRUE, type='plot', allow.null.country=TRUE), container=button.g)
	addSpace(button.g, 5)
	TableB.show.cohort.act <- gaction(label='Table', icon='dataframe', handler=show.pop.cohorts, 
						action=list(mw=main.win, env=e, type='table', script=FALSE, allow.null.country=TRUE))
	GraphB.show.cohort.act <- gaction(label='Graph', icon='lines', handler=show.pop.cohorts, 
						action=list(mw=main.win, env=e, type='plot', script=FALSE, allow.null.country=TRUE))
	bDem.gbutton(action=TableB.show.cohort.act, container=button.g)
	bDem.gbutton(action=GraphB.show.cohort.act, container=button.g)
	return(e)
}



show.pop.cohorts <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	show.type <- h$action$type
	allow.null.country <- if(is.null(h$action$allow.null.country)) FALSE else h$action$allow.null.country
	country.pars <- get.country.code.from.widget(e$show.cohort.country$country.w, e$show.cohort.country, 
							 allow.null.country=allow.null.country)
	if(is.null(country.pars)) {
		if(!allow.null.country) return(NULL)
		else if(!.pop.traj.country.check(e))
				 return(NULL)
	}
	param.names.all <- list(text=c('sim.dir'), numvector=c('pi',  'cohorts'), 
							numeric=c('dev.ncol'), logical=c('show.legend'))
	param.env <- get.parameters(param.names.all, env=e, quote=h$action$script)
	param.env <- c(param.env, get.parameters(list(text=c('legend.pos', 'expression')), env=e, quote=show.type=="plot" || h$action$script==TRUE))
	param.env.rest <- list(country=country.pars$code, aggregation=e$aggregation)
	param.env <- c(param.env, get.parameters(list(text=c('aggregation'), numeric='country'), 
									param.env.rest, quote=TRUE,
									retrieve.from.widgets=FALSE))
	if(!svalue(e$use.expression)) param.env[['expression']] <- NULL
	else param.env[['country']] <- NULL
	param.pred <- param.env['sim.dir']
	if(is.element('aggregation', names(param.env))) param.pred <- c(param.pred, param.env['aggregation'])
	# get it now unquoted (to avoid double quotes if script is TRUE)
	param.pred.ev <- c(get.parameters(list(text=c('sim.dir')), env=e, quote=FALSE),
			get.parameters(list(text=c('aggregation')), env=e, quote=FALSE, retrieve.from.widgets=FALSE))		 
	pred <- do.call('get.pop.prediction', param.pred.ev)
	if(is.null(pred)) {
		gmessage('Simulation directory contains no prediction of the specified type.', 
					title='Input Error', icon='error')
		return(NULL)
	}
	if(h$action$script) {
		cmd <- paste('pred <- get.pop.prediction(', assemble.arguments(param.pred), ')\n', sep='')
	} else {	
		cmd <- ''
	}
	param.plot1c <- list()
	for(item in names(param.env)) {
		if(!item %in% c('sim.dir', 'aggregation'))
			param.plot1c[[item]] <- param.env[[item]]
	}
	pars.value <- svalue(e$graph.pars)
	if (show.type == 'plot') {
		#param.plot1c <- param.env #[c('expression', 'pi', 'cohorts', 'nr.traj', 'show.legend', 'dev.ncol', 'legend.pos')]
		#if(is.element('country', names(param.env))) param.plot1c <- c(param.plot1c, param.env['country'])
		cmd <- paste(cmd, 'pop.cohorts.plot', '(pred,', assemble.arguments(param.plot1c, pars.value), ')\n')
		if (h$action$script) {
				create.script.widget(cmd, h$action$mw, package="bayesPop")
			} else {
				create.graphics.window(parent=h$action$mw, title=paste("Trajectories", 
							if(!is.null(country.pars$name) && !is.null(param.env[['country']])) paste("for", country.pars$name) else ""))
				eval(parse(text=cmd))
			}
	} else {
		# Table
		param.table <- param.env['pi']
		param.table <- c(param.table, if('country' %in% names(param.env)) param.env['country'] else param.env['expression'])
		table.values <- do.call('cohorts', c(list(pred), param.table))
		table.values[['last.observed']] <- NULL
		# select cohorts
		if(!is.null(param.env$cohorts)) {
			.round.to.lower5 <- function(x) 5*floor(x/5) 
			from.cohorts <- .round.to.lower5(param.env$cohorts)
			cohorts <- paste(from.cohorts, '-', from.cohorts+5, sep="")
			table.values <- table.values[cohorts]
		}
		# print
		cohort.table <- c()
		con <- textConnection("cohort.table", "w", local=TRUE)
		sink(con)
		print(table.values)
		close(con)
		create.script.widget(cohort.table, h$action$mw, 'bayesPop', title="Cohort data")
	}
}


pop.show.pyramid.group <- function(g, main.win, parent.env) {
	leftcenter <- c(-1, 0)
	e <- new.env()
	e$sim.dir <- parent.env$sim.dir
	e$pred.type <- 'pop'
	e$new.country.select.if.changed <- c('aggregation.dl')
	defaults.pred <- formals(pop.predict)
	defaults.pyr <- formals(bayesPop:::pop.pyramid.bayesPop.prediction)
	defaults.pyr.all <- formals(pop.pyramidAll)
	defaults.traj.pyr <- formals(bayesPop:::pop.trajectories.pyramid.bayesPop.prediction)
		
	addSpace(g, 10)
	country.f <- gframe("<span color='blue'>Country settings</span>", markup=TRUE, 
									horizontal=FALSE, container=g)
	glo <- glayout(container=country.f)
	glo[1,1, anchor=leftcenter] <- 'Aggregation:'
	glo[1,2] <- e$aggregation.dl <- bDem.gbutton(" None ", container=glo,
				handler=selectAggregationMenuPop,
				action=list(mw=main.win, env=e, label.widget.name='aggregation.dl'))
	e$show.pyr.country <- create.country.widget(country.f, defaults.pyr.all, 
									main.win, glo=glo, start.row=2, prediction=TRUE, parent.env=e, disable.table.button=FALSE)
	addSpace(g, 10)
	type.settings.f <- gframe("<span color='blue'>Pyramid settings</span>", markup=TRUE, 
								horizontal=FALSE, container=g)
	type.g1 <- glayout(horizontal=TRUE, container=type.settings.f)
	type.g1[1,1, anchor=leftcenter] <- glabel('CI (%):', container=type.g1)
	type.g1[1,2] <- e$pi <- gedit('80, 95', width=7, container=type.g1)
	type.g1[2,1:2] <- e$is.traj.pyr <- gcheckbox("Trajectory pyramid", container=type.g1, checked=FALSE,
						handler=function(h,...) enabled(e$nr.traj) <- svalue(h$obj))
	type.g1[3,1:2] <- e$proportion <- gcheckbox("x-axis as proportion", container=type.g1, checked=defaults.pyr$proportion)
	type.g1[2,3, anchor=leftcenter] <- glabel('# trajectories:', container=type.g1)
	type.g1[2,4] <- e$nr.traj <- gedit(20, width=5, container=type.g1)
	type.g1[3,3, anchor=leftcenter] <- glabel('Highest age category:', container=type.g1)
	type.g1[3,4] <- e$age <- gedit(defaults.pyr$age[length(defaults.pyr$age)], width=3, container=type.g1)
	tooltip(e$age) <- "21: 100-104, ..., 27: 130+"
	type.g1[2,5] <- year.gb <- bDem.gbutton(" Years ", container=type.g1,
				handler=selectYearsMenuPop,
				action=list(mw=main.win, env=e, widget='year', multiple=TRUE))
	type.g1[2,6] <- e$year <- gedit(defaults.pred$present.year, width=10, container=type.g1)
	tooltip(e$year) <- "There will be one plot per year."
	type.g1[3,5] <- year.comp.gb <- bDem.gbutton(" Additional years ", container=type.g1,
				handler=selectYearsMenuPop,
				action=list(mw=main.win, env=e, widget='year.comp', multiple=TRUE))
	type.g1[3,6] <- e$year.comp <- gedit('', width=10, container=type.g1)
	tooltip(e$year.comp) <- "All will be shown in the same plot."
	enabled(e$nr.traj) <- svalue(e$is.traj.pyr)
	addSpring(g)
	button.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='pop.pyramid', package='bayesPop', parent.group=button.g,
						parent.window=main.win)
	addSpring(button.g)
	create.generate.script.button(handler=show.pop.pyramid, 
							action=list(mw=main.win, env=e, script=TRUE), container=button.g)
	addSpace(button.g, 5)
	GraphB.show.pyr.act <- gaction(label='Graph', icon='lines', handler=show.pop.pyramid, 
						action=list(mw=main.win, env=e, script=FALSE))
	bDem.gbutton(action=GraphB.show.pyr.act, container=button.g)
	return(e)
}

show.pop.pyramid <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	country.pars <- get.country.code.from.widget(e$show.pyr.country$country.w, e$show.pyr.country, 
							force.country.spec=FALSE)
	if(is.null(country.pars)) return(NULL)
	traj.pyramid <- svalue(e$is.traj.pyr)
	param.names.all <- list(text='sim.dir', numvector=c('pi', 'year', 'year.comp'), 
							numeric='age', logical=c('proportion'))
	if(traj.pyramid) param.names.all$numeric <- c(param.names.all$numeric, 'nr.traj')
	param.env <- get.parameters(param.names.all, env=e, quote=h$action$script)
	param.env.rest <- list(country=country.pars$code, output.dir=country.pars$output.dir,
							output.type=country.pars$output.type, verbose=TRUE, aggregation=e$aggregation)
	param.env <- c(param.env, get.parameters(list(text=c('output.dir', 'output.type', 'aggregation'), 
										logical='verbose', numeric='country'), 
									param.env.rest, quote=TRUE,
									retrieve.from.widgets=FALSE))
	param.env$age <- if(!is.null(param.env$age)) 1:param.env$age else NULL
	param.pred <- param.env['sim.dir']
	if(is.element('aggregation', names(param.env))) param.pred <- c(param.pred, param.env['aggregation'])
	#param.pred.ev <- param.pred
	# get it now unquoted (to avoid double quotes if script is TRUE)
	param.pred.ev <- c(get.parameters(list(text=c('sim.dir')), env=e, quote=FALSE),
			get.parameters(list(text=c('aggregation')), env=e, quote=FALSE, retrieve.from.widgets=FALSE))		 
	pred <- do.call('get.pop.prediction', param.pred.ev)
	if(is.null(pred)) {
		gmessage('Simulation directory contains no prediction of the specified type.', 
					title='Input Error', icon='error')
		return(NULL)
	}
	if(h$action$script) {
		cmd <- paste('pred <- get.pop.prediction(', assemble.arguments(param.pred), ')\n', sep='')
	} else {	
		cmd <- ''
	}
	param.plot1c <- param.env[c('pi', 'age', 'year', 'proportion')]
	if(traj.pyramid) {
		param.plot1c <- c(param.plot1c, param.env['nr.traj'])
		pyr.command <- 'pop.trajectories.pyramid'
	} else {
		pyr.command <- 'pop.pyramid'
	}
	comp.years <- if(is.null(param.env$year.comp)) c() else param.env$year.comp
	if(is.element('country', names(param.env))) param.plot1c <- c(param.plot1c, param.env['country'])

	if(!is.null(param.env$country)) { # one country
		years <- param.env$year
		if (h$action$script) {
			for(year in years) {
				param.plot1c$year <- c(year, comp.years)
				cmd <- paste(cmd, pyr.command, '(pred,', assemble.arguments(param.plot1c), ')\n')
			}
			create.script.widget(cmd, h$action$mw, package='bayesPop')
		} else {
			base.cmd <- cmd
			for(year in years) {
				param.plot1c$year <- c(year, comp.years)
				cmd <- paste(base.cmd, pyr.command, '(pred,', assemble.arguments(param.plot1c), ')\n')
				create.graphics.window(parent=h$action$mw, title=paste("Pyramid for", country.pars$name), 
										ps=9, width=700, height=500)
				eval(parse(text=cmd))
			}
		}
	} else { # all countries
		years <- lapply(param.env$year, append, comp.years)
		year.str <- paste('year=list(',paste(years, collapse=', '), '), ', sep='')
		param.plot1c$year <- NULL
		param.plot.allc <- param.env[c(names(param.plot1c), 'output.dir', 'output.type',  'verbose')]
		cmd <- paste(cmd, pyr.command, 'All(pred, ', year.str, assemble.arguments(param.plot.allc), ')', sep='')
		if (h$action$script) {
			create.script.widget(cmd, h$action$mw, package='bayesPop')
		} else {
			eval(parse(text=cmd))
		}
	}
}

selectAgeMenuPop <- function(h, ...) {
	age.selected <- function(h1, ...) {
		age <- svalue(h$action$env$age.gt)
		h$action$env$age <- if (is.element(0, age)) 'all' else age 
		visible(h$action$env$age.sel.win) <- FALSE
		if(!is.null(h$action$label.widget.name) && !is.null(h$action$env[[h$action$label.widget.name]])) 
			svalue(h$action$env[[h$action$label.widget.name]]) <- if (h$action$env$age[1] == 'all') '' 
										else paste(bayesPop:::get.age.labels(h$action$env$age, collapse=TRUE, age.is.index=TRUE), 
																		collapse=',')
	}
	if (!is.null(h$action$env$age.sel.win)) 
		visible(h$action$env$age.sel.win) <- TRUE
	else {
		ages <- c('All', bayesPop:::get.age.labels(seq(0, by=5, length=27)))
		h$action$env$age.table <- data.frame(Index=seq(0,length=length(ages)), Age=ages)
		h$action$env$age.sel.win <- win <- gwindow('Select ages', 
							parent=h$action$mw, height=450, width=200,
							handler=function(h, ...) {
								h$action$env$age.sel.win<-NULL;
								h$action$env$age.ok.handler <- NULL
							},
							action=list(env=h$action$env))
		t.group <- ggroup(horizontal=FALSE, container=win)
		h$action$env$age.gt <- gtable(h$action$env$age.table, container=t.group, 
					expand=TRUE, multiple=h$action$multiple, handler=age.selected)
		b.group <- ggroup(horizontal=TRUE, container=t.group)
		gbutton('Cancel', container=b.group, handler=function(h, ...) 
					visible(win) <- FALSE)
		addSpring(b.group)
		h$action$env$age.okbutton <- gbutton('OK', container=b.group)
	}
	if(!is.null(h$action$env$age.ok.handler)) 
		removehandler(h$action$env$age.okbutton, h$action$env$age.ok.handler)
	h$action$env$age.ok.handler <- addhandlerclicked(
						h$action$env$age.okbutton, handler=age.selected)

}
selectYearsMenuPop <- function(h, ...) {
	w <- h$action$widget
	year.selected <- function(h1, ...) {
		years <- svalue(h$action$env$year.gt[[w]])
		svalue(h$action$env[[w]]) <- paste(years, collapse=',')
		visible(h$action$env$year.sel.win[[w]]) <- FALSE
	}
	if (is.null(h$action$env$year.sel.win)) h$action$env$year.sel.win <- list()
	if (is.null(h$action$env$year.gt)) h$action$env$year.gt <- list()
	if (is.null(h$action$env$year.ok.handler)) h$action$env$year.ok.handler <- list()
	if (!is.null(h$action$env$year.sel.win[[w]])) 
		visible(h$action$env$year.sel.win[[w]]) <- TRUE
	else {
		if(!has.required.arguments(list(sim.dir='Simulation directory'), env=h$action$env)) return()
		param <-list(sim.dir=svalue(h$action$env$sim.dir))
		pred <- do.call('get.pop.prediction', param)
		if(is.null(pred)) {
			gmessage('Simulation directory contains no valid population projection.', 
						title='Input Error', icon='error')
        	return(NULL)
		}
		h$action$env$year.table <- data.frame(Mid.year=c(as.integer(colnames(pred$inputs$pop.matrix[['male']])), pred$proj.years[-1]))
		h$action$env$year.sel.win[[w]] <- win <- gwindow('Select years', 
							parent=h$action$mw, height=450, width=100,
							handler=function(h1, ...) {
								h$action$env$year.sel.win[[w]]<-NULL;
								h$action$env$year.ok.handler[[w]] <- NULL
							},
							action=list(env=h$action$env))
		t.group <- ggroup(horizontal=FALSE, container=win)
		h$action$env$year.gt[[w]] <- gtable(h$action$env$year.table, container=t.group, 
					expand=TRUE, multiple=h$action$multiple, handler=year.selected)
		b.group <- ggroup(horizontal=TRUE, container=t.group)
		gbutton('Cancel', container=b.group, handler=function(h, ...) 
					visible(win) <- FALSE)
		addSpring(b.group)
		if (is.null(h$action$env$year.okbutton)) h$action$env$year.okbutton <- list()
		h$action$env$year.okbutton[[w]] <- gbutton('OK', container=b.group)
	}
	if(!is.null(h$action$env$year.ok.handler[[w]])) 
		removehandler(h$action$env$year.okbutton[[w]], h$action$env$year.ok.handler[[w]])
	h$action$env$year.ok.handler[[w]] <- addhandlerclicked(
						h$action$env$year.okbutton[[w]], handler=year.selected)

}

selectAggregationMenuPop <- function(h, ...) {
	aggr.selected <- function(h1, ...) {
		aggr <- h$action$env$aggregation <- as.character(svalue(h$action$env$aggr.gt))	
		visible(h$action$env$aggr.sel.win) <- FALSE
		if(length(aggr) == 0) aggr <- 'None' # OK button was pressed without making selection
		if(!is.null(h$action$label.widget.name) && !is.null(h$action$env[[h$action$label.widget.name]])) 
			svalue(h$action$env[[h$action$label.widget.name]]) <- aggr
		if(aggr == 'None') h$action$env$aggregation <- NULL
	}
	get.avail.aggregations <- function() {
		if(!has.required.arguments(list(sim.dir='Simulation directory'), env=h$action$env)) return()
		param <-list(sim.dir=svalue(h$action$env$sim.dir))
		pred <- do.call('get.pop.prediction', param)
		if(is.null(pred)) {
			gmessage('Simulation directory contains no valid population projection.', 
						title='Input Error', icon='error')
        	return(NULL)
		}
		return(data.frame(Name=c('None', bayesPop:::available.pop.aggregations(pred))))
	}
	aggrs <- get.avail.aggregations()
	if (!is.null(h$action$env$aggr.sel.win) && setequal(aggrs$Name, h$action$env$aggr.table$Name))
		visible(h$action$env$aggr.sel.win) <- TRUE
	else {		
		h$action$env$aggr.table <- aggrs
		h$action$env$aggr.sel.win <- win <- gwindow('Select aggregation', 
							parent=h$action$mw, height=450, width=200,
							handler=function(h, ...) {
								h$action$env$aggr.sel.win<-NULL;
								h$action$env$aggr.ok.handler <- NULL
							},
							action=list(env=h$action$env))
		t.group <- ggroup(horizontal=FALSE, container=win)
		h$action$env$aggr.gt <- gtable(h$action$env$aggr.table, container=t.group, 
					expand=TRUE, multiple=FALSE, handler=aggr.selected)
		b.group <- ggroup(horizontal=TRUE, container=t.group)
		gbutton('Cancel', container=b.group, handler=function(h, ...) 
					visible(win) <- FALSE)
		addSpring(b.group)
		h$action$env$aggr.okbutton <- gbutton('OK', container=b.group)
	}
	if(!is.null(h$action$env$aggr.ok.handler)) 
		removehandler(h$action$env$aggr.okbutton, h$action$env$aggr.ok.handler)
	h$action$env$aggr.ok.handler <- addhandlerclicked(
						h$action$env$aggr.okbutton, handler=aggr.selected)

}


pop.show.map.group <- function(g, main.win, parent.env) {
	leftcenter <- c(-1, 0)
	e <- new.env()
	e$sim.dir <- parent.env$sim.dir
	addSpace(g, 10)
	lo <- .create.map.settings.group(g, e, measures=c('Population', 'Expression'))
	lo[4, 1, anchor=leftcenter] <- "Sex:"
	lo[4, 2] <- e$sex <- bDem.gdroplist(c('Both', 'Female', 'Male'), container=lo, selected=1)
	lo[4, 4] <- age.gb <- bDem.gbutton(" Age ", container=lo,
				handler=selectAgeMenuPop,
				action=list(mw=main.win, env=e, multiple=TRUE, label.widget.name='age.label'))
	e$age <- 'all'
	lo[4,5, anchor=leftcenter] <- e$age.label <- glabel('', container=lo)
	lo[5, 1, anchor=leftcenter] <- "Expression:"
	lo[5, 2:5] <- e$expression <- gedit('', container=lo)
	tooltip(e$expression) <- "See ?pop.expressions. Use XXX as country code."
	addHandlerChanged(e$map.measure, function(h, ...) {
					is.expression <- svalue(h$obj) == 'Expression'
					enabled(e$expression) <- is.expression
					enabled(e$sex) <- !is.expression
					enabled(age.gb) <- !is.expression
					})
	enabled(e$expression) <- FALSE
	addSpring(g)
	bg <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='pop.map', package='bayesPop', parent.group=bg,
						parent.window=main.win)
	addSpring(bg)
	create.generate.script.button(handler=pop.showMap, action=list(mw=main.win, env=e, script=TRUE),
								container=bg)
	addSpace(bg, 5)
	GraphB.map <- gaction(label=' Show Map ', handler=pop.showMap, 
						action=list(mw=main.win, env=e, script=FALSE))
	bDem.gbutton(action=GraphB.map, container=bg)
}

pop.showMap <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	param.names1 <- list(text='sim.dir')
	param.pred <- get.parameters(param.names1, env=e, quote=h$action$script)
	param.names.rest <- list(text=c('sex', 'expression'))
	param.rest <- get.parameters(param.names.rest, env=e, quote=h$action$script)
	param.rest$age <- if (is.element(0, e$age) || e$age=='all') 
		get.parameters(list(text='age'), env=e, retrieve.from.widgets=FALSE, quote=h$action$script)$age
					 else e$age 
	param.rest2 <- get.parameters(list(logical=c('same.scale')), env=e, quote=h$action$script)
	param.rest$sex <- tolower(param.rest$sex)
	if(svalue(e$map.measure) != 'Expression') param.rest$expression <- NULL
	else for(par in c('sex', 'age')) param.rest[[par]] <- NULL
	
	percentile <- svalue(e$map.percentile)
	param.rest$quantile <- e$percentiles[[percentile]]
	bounds <- svalue(e$map.bounds)
	package <- svalue(e$map.package)
	map.function <- if(package == 'rworldmap') 'pop.map' else 'pop.map.gvis'
	if(h$action$script) {
		cmd <- paste('pred <- get.pop.prediction(', assemble.arguments(param.pred), ')\n', sep='')
		if(package == 'rworldmap') {
			cmd <- paste(cmd, "param.map <- get.pop.map.parameters(pred, ", 
						assemble.arguments(c(param.rest, param.rest2)), ")\n", sep="")
			cmd <- paste(cmd, 'do.call("', map.function, '", param.map)', sep='')
		} else {
			cmd <- paste(cmd, map.function, '(pred, ', assemble.arguments(c(param.rest, list(pi=bounds))), ')', sep='')
		}
		create.script.widget(cmd, h$action$mw, package="bayesPop")
	} else {
		pred <- do.call('get.pop.prediction', param.pred)
		param.map <-  if (package == 'rworldmap') do.call('get.pop.map.parameters', c(list(pred), param.rest, param.rest2))
						else c(list(pred=pred), param.rest)
		if(package == 'rworldmap') param.map[['device']] <- 'dev.new'
		if (package == 'googleVis') param.map[['pi']] <- bounds
		g <- create.graphics.map.window(parent=h$action$mw, pred=pred, params=param.map, percentile=percentile, 
										is.gvis= package == 'googleVis', title="World Map", type='pop', 
										cw.main=paste(.map.main.default(pred, param.map), param.rest$quantile))
	}
}

pop.get.time.info <- function(pred) {
	return(list(est.periods=bayesPop:::get.pop.observed.periods(pred), 
				proj.periods=bayesPop:::get.pop.prediction.periods(pred),
				proj.ind.func=bayesPop:::get.predORobs.year.index,
				present.year=pred$proj.years[1],
				est.years=pred$estim.years,
				proj.years=pred$proj.years
				))
	
}