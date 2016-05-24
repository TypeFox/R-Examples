popNewPred.group <- function(g, main.win, parent) {
	nb <- bDem.gnotebook(container=g, expand=TRUE)
	c.g <- ggroup(label="  <span color='#0B6138'>Countries</span>  ", 
							markup=TRUE, horizontal=FALSE, expand=TRUE, container=nb)
	pop.pred.countries.group(c.g, main.win, parent)
	a.g <- ggroup(label="  <span color='#0B6138'>Aggregations</span>  ", 
							markup=TRUE, horizontal=FALSE, expand=TRUE, container=nb)
	pop.pred.aggregation.group(a.g, main.win, parent)
	svalue(nb) <- 1
}

pop.pred.countries.group <- function(g, main.win, parent) {
	e <- new.env()
	defaults <- formals(pop.predict) # default argument values
	e$output.dir <- e$sim.dir <- parent$sim.dir
	leftcenter <- c(-1,0)
	addSpace(g, 10)
	pred.g <- gframe("<span color='blue'>Prediction settings</span>", markup=TRUE, horizontal=FALSE, container=g)
	pred.g1 <- glayout(container=pred.g)
	pred.g1[1,1, anchor=leftcenter] <- glabel("End year:", container=pred.g1)
	pred.g1[1,2] <- e$end.year <- gedit(defaults$end.year, width=4, container=pred.g1)
	tooltip(e$end.year) <- 'End year of the prediction.'
	pred.g1[1,3, anchor=leftcenter] <- glabel("     WPP year:", container=pred.g1)
	pred.g1[1,4, anchor=leftcenter] <- wpp <- glabel(parent$wpp.year, container=pred.g1)
	tooltip(wpp) <- 'To change this start bayesDem with wpp.year as (third) argument.'
	pred.g1[2,1, anchor=leftcenter] <- glabel("Start year:", container=pred.g1)
	pred.g1[2,2] <- e$start.year <- gedit(defaults$start.year, width=4, container=pred.g1)
	tooltip(e$start.year) <- 'Historical data on mortality prior to this year will be ignored.'
	pred.g1[2,3, anchor=leftcenter] <- glabel("Present year:", container=pred.g1)
	pred.g1[2,4] <- e$present.year <- gedit(defaults$present.year, width=4, container=pred.g1)
	tooltip(e$present.year) <- 'Initial year of the population data.'
	pred.g1[1,5] <- '    '
	pred.g1[3,1, anchor=leftcenter] <- glabel("Nr. trajectories:", container=pred.g1)
	pred.g1[3,2] <- e$nr.traj <- gedit(defaults$nr.traj, width=5, container=pred.g1)
	tooltip(e$nr.traj) <- 'If left empty, #trajectories of TFR or e0 is used.'
	pred.g1[1,6:7] <- e$keep.vital.events <- gcheckbox("Keep vital events", checked=defaults$keep.vital.events, container=pred.g1)
	pred.g1[2,6:7] <- e$fixed.mx <- gcheckbox("Fixed mortality", checked=defaults$fixed.mx, container=pred.g1)
	tooltip(e$fixed.mx) <- 'If selected, future mx is used instead of e0.'
	pred.g1[3,6:7] <- e$fixed.pasfr <- gcheckbox("Fixed percent fert.", checked=defaults$fixed.pasfr, container=pred.g1)
	tooltip(e$fixed.pasfr) <- 'If not selected, future percent fertility is computed on the fly.'
	pred.g1[1,8:9] <- e$verbose <- gcheckbox("Verbose", checked=defaults$verbose, container=pred.g1)
	
	addSpace(g, 10)
	countries.g <- gframe("<span color='blue'>Countries selection</span>", markup=TRUE, 
							horizontal=FALSE, container=g)
	countries.g1 <- glayout(horizontal=TRUE, container=countries.g)
	countries.g1[1,1] <- e$all.countries <- gcheckbox("All countries", checked=TRUE, container=countries.g1,
									handler=function(h,...){
										enabled(countries.gb) <- !svalue(h$obj)
										enabled(e$all.remaining.countries) <- !svalue(h$obj)
										})
	countries.g1[1,2] <- countries.gb <- bDem.gbutton("  Select countries  ", container=countries.g1,
				handler=selectCountryMenuPop,
				action=list(mw=main.win, env=e, multiple=TRUE, wpp.year=parent$wpp.year, label.widget.name='country.label'))
	countries.g1[2,1] <- e$all.remaining.countries <- gcheckbox("All countries without prediction", 
									checked=FALSE, container=countries.g1,
									handler=function(h,...){
										enabled(countries.gb) <- !svalue(h$obj)
										enabled(e$all.countries) <- !svalue(h$obj)
										})
	countries.g1[2,2, anchor=leftcenter] <- e$country.label <- glabel('', container=countries.g1)
	enabled(countries.gb) <- !svalue(e$all.countries)
	enabled(e$all.remaining.countries) <- !svalue(e$all.countries)
	
	e$inputs <- list()				
	addSpace(g, 20)
	bDem.gbutton("Projection inputs", markup=TRUE, container=g,
				handler=OptInputFilesPopPred,
				action=list(mw=main.win, env=e, defaults=defaults))
	addSpace(g, 10)
	.create.status.label(g, e)

	addSpring(g)
	button.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic=c('pop.predict', 'bayesPop-package'), package='bayesPop', 
				parent.group=button.g,
						parent.window=main.win)
	addSpring(button.g)
	create.generate.script.button(handler=run.pop.prediction, 
							action=list(mw=main.win, env=e, script=TRUE, wpp.year=parent$wpp.year), container=button.g)
	bDem.gbutton(action=gaction(label=' Make Prediction ', icon='evaluate', handler=run.pop.prediction, 
				action=list(mw=main.win, env=e, script=FALSE, wpp.year=parent$wpp.year)), 
				container=button.g)
	return(e)					
}

pop.pred.aggregation.group <- function(g, main.win, parent) {
	enable.prob.inputs <- function(enable) {
		enabled(e$inputs$tfr.sim.dir) <- enable
		enabled(e$inputs$e0F.sim.dir) <- enable
		enabled(e$inputs$e0M.sim.dir) <- enable
		enabled(e$inputs$e0M.joint) <- enable
		enabled(e$inputs$e0M.sim.dir) <- !svalue(e$inputs$e0M.joint)
	}
	e <- new.env()
	defaults <- formals(pop.aggregate) # default argument values
	e$sim.dir <- parent$sim.dir
	leftcenter <- c(-1,0)
	addSpace(g, 10)
	aggr.f <- gframe("<span color='blue'>Aggregation settings</span>", markup=TRUE, horizontal=FALSE, container=g)
	lo <- glayout(container=aggr.f)
	lo[1,1, anchor=leftcenter] <- 'Input Type:'
	lo[1,2] <- e$input.type <- gradio(c('country', 'region'), selected=1, container=lo, horizontal=TRUE,
									handler=function(h, ...) {
										enable <- svalue(h$obj) == 'region'
										enable.prob.inputs(enable)
									})
	lo[2,1, anchor=leftcenter] <- 'Name:'
	lo[2,2] <- e$name <- gedit('', container=lo)
	lo[3,1:2] <- bDem.gbutton('   Select regions   ', container=lo, handler=selectRegionMenuPop,
				action=list(mw=main.win, env=e, wpp.year=parent$wpp.year, multiple=TRUE, label.widget.name='region.label'))
	lo[3,3:4, anchor=leftcenter] <- e$region.label <- glabel('', container=lo)
	lo[1,3] <- '     ' 
	lo[1,4] <- e$verbose <- gcheckbox("Verbose", checked=defaults$verbose, container=lo)
	addSpace(g, 10)
	e$inputs <- list()
	inp.f <- gframe("<span color='blue'>Probabilistic inputs</span>", markup=TRUE, horizontal=FALSE, container=g)
	ilo <- glayout(container=inp.f)
	ilo[1,1, anchor=leftcenter] <- glabel("TFR sim directory:", container=ilo)
	ilo[1,2:4] <- e$inputs$tfr.sim.dir <- bDem.gfilebrowse(eval(defaults$inputs$tfr.sim.dir), type='selectdir', 
					  width=20, quote=FALSE, container=ilo)
	ilo[2,1, anchor=leftcenter] <- glabel("Female e0 sim directory:", container=ilo)
	ilo[2,2:4] <- e$inputs$e0F.sim.dir <- bDem.gfilebrowse(eval(defaults$inputs$e0F.sim.dir), type='selectdir', 
					  width=20, quote=FALSE, container=ilo)
	ilo[3,1, anchor=leftcenter] <- glabel("Male e0 sim directory:", container=ilo)
	ilo[3,2] <- e$inputs$e0M.joint <- gcheckbox("Joint with Female", checked=defaults$inputs$e0M.sim.dir == 'joint_', 
						container=ilo,
							handler=function(h,...) {
								enabled(e$inputs$e0M.sim.dir) <- !svalue(h$obj)
							})
	ilo[3,3, anchor=leftcenter] <- glabel("Sim directory:", container=ilo)
	ilo[3,4] <- e$inputs$e0M.sim.dir <- bDem.gfilebrowse('', type='selectdir', 
					  width=15, quote=FALSE, container=ilo)
	enable.prob.inputs(FALSE)
	addSpace(g, 10)
	.create.status.label(g, e)
	
	addSpring(g)
	button.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic=c('pop.aggregate', 'bayesPop-package'), package='bayesPop', 
				parent.group=button.g,
						parent.window=main.win)
	addSpring(button.g)
	create.generate.script.button(handler=run.pop.aggregation, 
							action=list(mw=main.win, env=e, script=TRUE), container=button.g)
	bDem.gbutton(action=gaction(label=' Make Prediction ', icon='evaluate', handler=run.pop.aggregation, 
				action=list(mw=main.win, env=e, script=FALSE)), 
				container=button.g)
}

run.pop.prediction <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(output.dir='Simulation directory', end.year='End year', 
									start.year='Start year', present.year='Present year'), env=e)) return()
	param.names <- list(numeric=c('end.year', 'start.year', 'present.year', 'nr.traj'), 
						text=c('output.dir'),
						logical=c('verbose', 'keep.vital.events', 'fixed.mx', 'fixed.pasfr')
						)
	param.input.names.opt <- list(
		text=c('tfr.sim.dir', 'tfr.file', 'e0M.sim.dir', 'e0M.file', 'e0F.sim.dir', 'e0F.file', 'migMtraj', 'migFtraj',
				'popM', 'popF', 'mxM', 'mxF', 'srb', 'migM', 'migF', 'mig.type'))
	params <- get.parameters(param.names, e, h$action$script)
	params.inp <- list()
	if(!is.null(e$opt) && e$inputs.modified) {
		opt.names <- list(text=c())
		for (par in param.input.names.opt$text)	
			if (is.element(par, names(e$inputs)) && !is.null(e$inputs[[par]])) 
				opt.names$text <- c(opt.names$text, par)
		if(e$inputs$e0M.joint && is.element('e0F.sim.dir', opt.names$text)) {
			e$inputs$e0M.sim.dir <- 'joint_'
			if(!is.element('e0M.sim.dir', opt.names$text))
				opt.names$text <- c(opt.names$text, 'e0M.sim.dir')	
		}
		params.inp <- c(params.inp, get.parameters(opt.names, e$inputs, h$action$script, 
							retrieve.from.widgets=FALSE))
	}
	params$wpp.year <- h$action$wpp.year
	params$countries <- NULL
	if (!svalue(e$all.countries)) {
		if (svalue(e$all.remaining.countries)) params$countries <- NA
		else params$countries <- e$selected.countries
	}
	simdir <- get.parameters(list(text='output.dir'), e, quote=FALSE)$output.dir # to avoid double quoting if script is TRUE
	if(has.pop.prediction(sim.dir=simdir)) {
		params[['replace.output']] <- FALSE
		if (gconfirm(paste('Prediction for', simdir, 
								'already exists.\nDo you want to overwrite existing results?'),
				icon='question', parent=h$action$mw))
			params[['replace.output']] <- TRUE
		else return(NULL)
	}

	if (h$action$script) {
		cmd <- paste('pop.predict(', assemble.arguments(params), sep=' ')
		if (length(params.inp) > 0)
			cmd <- paste(cmd, ', inputs = list(', assemble.arguments(params.inp),')', sep=' ')
		cmd <- paste(cmd, ')', sep='')
		create.script.widget(cmd, h$action$mw, package="bayesPop")
	} else {
		.run.prediction(e, type='Pop', handler=get.pop.prediction.status, option='bDem.Poppred', 
								call='pop.predict', 
								params=c(params, if (length(params.inp) > 0) list(inputs=params.inp) else NULL), 
								sim.name='Population prediction', main.win=h$action$mw,
								action=list(sb=e$statuslabel),
								interval=1000)
	}

}

run.pop.aggregation <- function(h, ...) {
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	param.names <- list(text=c('sim.dir', 'name', 'input.type'),
						logical=c('verbose')
						)
	params <- get.parameters(param.names, e, h$action$script)
	params.inputs <- list()
	if(svalue(e$input.type) == 'region') {
		param.input.names.opt <- list(text=c('tfr.sim.dir', 'e0M.sim.dir', 'e0F.sim.dir'))
		params.inputs <- c(params.inputs, get.parameters(param.input.names.opt, e$inputs, h$action$script))
		if(svalue(e$inputs$e0M.joint)) 
			params.inputs$e0M.sim.dir <- if (h$action$script) "'joint_'" else 'joint_'
	}
	params$regions <- e$selected.regions
	if(is.null(params$regions)) {
		gmessage("No region selected.")
		return(NULL)
	}
	simdir <- get.parameters(list(text='sim.dir'), e, quote=FALSE)$sim.dir # to avoid double quoting if script is TRUE
	if(!has.pop.prediction(sim.dir=simdir)) {
		gmessage("No prediction exists in the simulation directory.")
		return(NULL)
	}
	pop.pred <- get.pop.prediction(simdir)
	if (h$action$script) {
		cmd <- paste('pred <- get.pop.prediction(sim.dir=', params[['sim.dir']], ')', sep='')
		params[['sim.dir']] <- NULL
		cmd <- paste(cmd, '\npop.aggregate(pred,', assemble.arguments(params), sep=' ')
		if (length(params.inputs) > 0)
			cmd <- paste(cmd, ', inputs = list(', assemble.arguments(params.inputs),')', sep=' ')
		cmd <- paste(cmd, ')', sep='')
		create.script.widget(cmd, h$action$mw, package="bayesPop")
	} else {
		params[['sim.dir']] <- NULL
		.run.prediction(e, type='PopAg', handler=get.pop.aggregation.status, option='bDem.PopAgpred', 
								call='pop.aggregate', 
								params=c(list(pop.pred), params, if (length(params.inputs) > 0) list(inputs=params.inputs) else NULL), 
								sim.name='Population aggregations', main.win=h$action$mw,
								action=list(sb=e$statuslabel),
								interval=1000)
	}

}

get.pop.prediction.status <- function(h, ...) 
	.update.status(h$action$sb, 'bDem.Poppred.status', 'Running Pop prediction ...')
	
get.pop.aggregation.status <- function(h, ...) 
	.update.status(h$action$sb, 'bDem.PopAgpred.status', 'Running Pop aggregation ...')


OptInputFilesPopPred <- function(h, ...) {
	input.names <- c('tfr.sim.dir', 'tfr.file', 'e0M.sim.dir', 'e0M.file', 'e0M.joint', 'e0F.sim.dir', 'e0F.file', 'migMtraj', 'migFtraj',
				'popM', 'popF', 'mxM', 'mxF', 'srb', 'migM', 'migF', 'mig.type')
	defaults <- h$action$defaults
	defaults$inputs$e0M.joint <- TRUE
	
	setOptionalInputsPop <- function(h1, ...) {
		for(par in input.names) {
			value <- svalue(h$action$env$opt[[par]])
			if (nchar(value) > 0 || ((nchar(value)==0) && is.element(par, names(h$action$env$inputs)) && nchar(h$action$env$inputs[[par]]) > 0))
				h$action$env$inputs[[par]] <- value
		}
		h$action$env$inputs.modified <- TRUE
		visible(h$action$env$inputs.win) <- FALSE
	}
	if (!is.null(h$action$env$inputs.win) && !h$action$env$opt$window.destroyed) { #Window exists
		if(h$action$env$inputs.modified) { # OK button previously clicked
			for(par in input.names) 
				if(is.element(par, names(h$action$env$inputs))) 
					svalue(h$action$env$opt[[par]]) <- h$action$env$inputs[[par]]
		} else # OK button not clicked yet, values are set to defaults
			for(par in input.names) svalue(h$action$env$opt[[par]]) <- defaults$inputs[[par]]
		visible(h$action$env$inputs.win) <- TRUE
	} else { # create the inputs window
		h$action$env$inputs.win <- win <- bDem.gwindow('Select optional input files and directories', 
							parent=h$action$mw)
		h$action$env$inputs.modified <- FALSE
		e <- new.env() # h$action$env
		g <- ggroup(horizontal=FALSE, container=win)
		g.prob <- gframe("<span color='blue'>Probabilistic Inputs (trajectories)</span>", markup=TRUE, 
							horizontal=FALSE, container=g)
		g.tfr <- gframe("<span color='#0B6138'>Total Fertility Rate (select one)</span>", markup=TRUE, 
							horizontal=TRUE, container=g.prob)
		glabel("bayesTFR sim folder:", container=g.tfr)
		e$tfr.sim.dir <- bDem.gfilebrowse(eval(defaults$inputs$tfr.sim.dir), type='selectdir', 
					  width=20, quote=FALSE, container=g.tfr)
		glabel("CSV file:", container=g.tfr)
		e$tfr.file <- bDem.gfilebrowse(eval(defaults$inputs$tfr.sim.dir), type='open', 
					  width=20, quote=FALSE, container=g.tfr)
		
		g.e0f <- gframe("<span color='#0B6138'>Female Life Expectancy (select one)</span>", markup=TRUE, 
							horizontal=TRUE, container=g.prob)
		glabel("bayesLife sim folder: ", container=g.e0f)
		e$e0F.sim.dir <- bDem.gfilebrowse(eval(defaults$inputs$e0F.sim.dir), type='selectdir', 
					  width=20, quote=FALSE, container=g.e0f)
		glabel("CSV file:", container=g.e0f)
		e$e0F.file <- bDem.gfilebrowse(eval(defaults$inputs$e0F.sim.dir), type='open', 
					  width=20, quote=FALSE, container=g.e0f)
					  			  
		g.e0m <- gframe("<span color='#0B6138'>Male Life Expectancy (select one)</span>", markup=TRUE, 
							horizontal=TRUE, container=g.prob)
		e$e0M.joint <- gcheckbox("Joint with Female", checked=defaults$inputs$e0M.joint, container=g.e0m,
							handler=function(h1,...) {
								enabled(e$e0M.sim.dir) <- !svalue(h1$obj)
								enabled(e$e0M.file) <- !svalue(h1$obj)
							})
		addSpace(g.e0m, 5)
		glabel("Sim folder:", container=g.e0m)
		e$e0M.sim.dir <- bDem.gfilebrowse(eval(defaults$inputs$e0M.sim.dir), type='selectdir', 
					  width=15, quote=FALSE, container=g.e0m)
		glabel("CSV file:", container=g.e0m)
		e$e0M.file <- bDem.gfilebrowse(eval(defaults$inputs$e0M.sim.dir), type='open', 
					  width=15, quote=FALSE, container=g.e0m)
		enabled(e$e0M.sim.dir) <- !svalue(e$e0M.joint)
		enabled(e$e0M.file) <- !svalue(e$e0M.joint)
		
		g.migtraj <- gframe("<span color='#0B6138'>Net Age-Specific Migration Counts</span>", markup=TRUE, 
							horizontal=TRUE, container=g.prob)
		glabel("CSV file for Male:      ", container=g.migtraj)
		e$migMtraj <- bDem.gfilebrowse(eval(defaults$inputs$migMtraj), type='open', 
					  width=20, quote=FALSE, container=g.migtraj)
		glabel(" Female:", container=g.migtraj)
		e$migFtraj <- bDem.gfilebrowse(eval(defaults$inputs$migFtraj), type='open', 
					  width=20, quote=FALSE, container=g.migtraj)
							
		g.other <- gframe("<span color='blue'>Deterministic Inputs (optional files)</span>", markup=TRUE, 
							horizontal=FALSE, container=g)
		#g.other <- gframe("<span color='#0B6138'>Other Optional Files</span>", markup=TRUE, 
		#					horizontal=FALSE, container=g)
		glo <- glayout(container=g.other)
		glo[1,1] <- glabel("Initial Male Population:", container=glo)
		glo[1,2] <- e$popM <- bDem.gfilebrowse(eval(defaults$inputs$popM), type='open', 
					  width=50, quote=FALSE, container=glo)
		glo[2,1] <- glabel("Initial Female Population:", container=glo)
		glo[2,2] <- e$popF <- bDem.gfilebrowse(eval(defaults$inputs$popF), type='open', 
					  width=50, quote=FALSE, container=glo)
					  
		glo[3,1] <- glabel("Mortality Rate Male:", container=glo)
		glo[3,2] <- e$mxM <- bDem.gfilebrowse(eval(defaults$inputs$mxM), type='open', 
					  width=50, quote=FALSE, container=glo)
		glo[4,1] <- glabel("Mortality Rate Female:", container=glo)
		glo[4,2] <- e$mxF <- bDem.gfilebrowse(eval(defaults$inputs$mxF), type='open', 
					  width=50, quote=FALSE, container=glo)
		glo[5,1] <- glabel("Sex Ratio at Birth:", container=glo)
		glo[5,2] <- e$srb <- bDem.gfilebrowse(eval(defaults$inputs$srb), type='open', 
					  width=50, quote=FALSE, container=glo)
		glo[6,1] <- glabel("% Age-specific Fertility Ratio:", container=glo)
		glo[6,2] <- e$pasfr <- bDem.gfilebrowse(eval(defaults$inputs$pasfr), type='open', 
					  width=50, quote=FALSE, container=glo)
		glo[7,1] <- glabel("Migration Male:", container=glo)
		glo[7,2] <- e$migM <- bDem.gfilebrowse(eval(defaults$inputs$migM), type='open', 
					  width=50, quote=FALSE, container=glo)
		glo[8,1] <- glabel("Migration Female:", container=glo)
		glo[8,2] <- e$migF <- bDem.gfilebrowse(eval(defaults$inputs$migF), type='open', 
					  width=50, quote=FALSE, container=glo)
		glo[9,1] <- glabel("Migration Type:    ", container=glo)
		glo[9,2] <- e$mig.type <- bDem.gfilebrowse(eval(defaults$inputs$mig.type), type='open', 
					  width=50, quote=FALSE, container=glo)

		b.group <- ggroup(horizontal=TRUE, container=g)
		bDem.gbutton('Cancel', container=b.group, handler=function(h1, ...) 
					visible(win) <- FALSE)
		addSpring(b.group)
		e$okbutton <- bDem.gbutton('OK', container=b.group)
		e$window.destroyed <- FALSE
		e$inputs.modified <- FALSE
		h$action$env$opt <- e
		addHandlerDestroy(win, handler=function(h1, ...) h$action$env$opt$window.destroyed <- TRUE)
	}
	if(!is.null(h$action$env$opt.okhandler)) 
		removehandler(h$action$env$opt$okbutton, h$action$env$opt.okhandler)
	h$action$env$opt.okhandler <- addhandlerclicked(h$action$env$opt$okbutton, 
											handler=setOptionalInputsPop)
}



selectCountryMenuPop <- function(h, ...) {
	country.selected <- function(h1, ...) {
		h$action$env$selected.countries <- svalue(h$action$env$country.gt)
		visible(h$action$env$country.sel.win) <- FALSE
		if(!is.null(h$action$label.widget.name) && !is.null(h$action$env[[h$action$label.widget.name]])) 
			svalue(h$action$env[[h$action$label.widget.name]]) <- paste('Selected countries:', 
								paste(h$action$env$selected.countries, collapse=', '))
	}
	if (!is.null(h$action$env$country.sel.win)) 
			visible(h$action$env$country.sel.win) <- TRUE
	else {
		SRB <- bayesPop:::load.wpp.dataset('sexRatio', h$action$wpp.year)
		country.table <- SRB[,c("country_code", "country")]
		h$action$env$country.table <- country.table
		h$action$env$country.sel.win <- win <- gwindow('Select countries', 
							parent=h$action$mw, height=450,
							handler=function(h, ...) {
								h$action$env$country.sel.win<-NULL;
								h$action$env$country.ok.handler <- NULL
							},
							action=list(env=h$action$env))
		t.group <- ggroup(horizontal=FALSE, container=win)
		h$action$env$country.gt <- gtable(h$action$env$country.table, container=t.group, 
					expand=TRUE, multiple=h$action$multiple, handler=country.selected)
		b.group <- ggroup(horizontal=TRUE, container=t.group)
		gbutton('Cancel', container=b.group, handler=function(h, ...) 
					visible(win) <- FALSE)
		addSpring(b.group)
		h$action$env$country.okbutton <- gbutton('OK', container=b.group)
	}
	if(!is.null(h$action$env$country.ok.handler)) 
		removehandler(h$action$env$country.okbutton, h$action$env$country.ok.handler)
	h$action$env$country.ok.handler <- addhandlerclicked(
						h$action$env$country.okbutton, handler=country.selected)

}

selectRegionMenuPop <- function(h, ...) {
	region.selected <- function(h1, ...) {
		h$action$env$selected.regions <- svalue(h$action$env$region.gt)
		visible(h$action$env$region.sel.win) <- FALSE
		if(!is.null(h$action$label.widget.name) && !is.null(h$action$env[[h$action$label.widget.name]])) 
			svalue(h$action$env[[h$action$label.widget.name]]) <- paste('Selected regions:', 
								paste(h$action$env$selected.regions, collapse=', '))
	}
	if (!is.null(h$action$env$region.sel.win)) 
			visible(h$action$env$region.sel.win) <- TRUE
	else {
		LOCATIONS <- bayesTFR:::load.bdem.dataset('UNlocations', h$action$wpp.year)		
		region.table <- LOCATIONS[is.element(LOCATIONS[,'location_type'], c(0,2,3)), c("country_code", "name")]
		names(region.table) <- c('code', 'name')
		h$action$env$region.table <- region.table
		h$action$env$region.sel.win <- win <- gwindow('Select regions', 
							parent=h$action$mw, height=450,
							handler=function(h, ...) {
								h$action$env$region.sel.win<-NULL;
								h$action$env$region.ok.handler <- NULL
							},
							action=list(env=h$action$env))
		t.group <- ggroup(horizontal=FALSE, container=win)
		h$action$env$region.gt <- gtable(h$action$env$region.table, container=t.group, 
					expand=TRUE, multiple=h$action$multiple, handler=region.selected)
		b.group <- ggroup(horizontal=TRUE, container=t.group)
		gbutton('Cancel', container=b.group, handler=function(h, ...) 
					visible(win) <- FALSE)
		addSpring(b.group)
		h$action$env$region.okbutton <- gbutton('OK', container=b.group)
	}
	if(!is.null(h$action$env$region.ok.handler)) 
		removehandler(h$action$env$region.okbutton, h$action$env$region.ok.handler)
	h$action$env$region.ok.handler <- addhandlerclicked(
						h$action$env$region.okbutton, handler=region.selected)
}

