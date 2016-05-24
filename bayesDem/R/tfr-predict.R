TFRnewPred.group <- function(g, main.win, parent) {
	nb <- bDem.gnotebook(container=g, expand=TRUE)
	all.c.g <- ggroup(label="  <span color='#0B6138'>All Countries</span>  ", 
							markup=TRUE, horizontal=FALSE, expand=TRUE, container=nb)
	pred.all <- tfr.pred.all.countries.group(all.c.g, main.win, parent)
	extra.c.g <- ggroup(label="  <span color='#0B6138'>Extra Areas &amp; Regions</span>  ", 
							markup=TRUE, horizontal=FALSE, expand=TRUE, container=nb)
	pred.extra <- tfr.pred.extra.countries.group(extra.c.g, main.win, parent)
	edit.g <- ggroup(label="  <span color='#0B6138'>Edit Predictions</span>  ", 
							markup=TRUE, horizontal=FALSE, expand=TRUE, container=nb)
	pred.edit <- edit.predictions.group(edit.g, main.win, parent)
	svalue(nb) <- 1
	return(pred.all)
}

tfr.pred.all.countries.group <- function(g, main.win, parent) {
	e <- new.env()
	defaults <- formals(tfr.predict) # default argument values
	e$sim.dir <- parent$sim.dir
	addSpace(g, 10)
	lo <- .create.prediction.setting.group(g, e, defaults)
	lo[2,6] <- e$use.correlation <- gcheckbox("Correlation between countries", checked=defaults$use.correlation, container=lo)
	addSpace(g, 10)
	ar.g <- gframe("<span color='blue'>Phase III Model</span>", markup=TRUE, horizontal=FALSE, container=g)
	p3.g1 <- ggroup(horizontal=TRUE, container=ar.g)
	glabel("Phase III modeled using ", container=p3.g1)
	e$p3.type <- gradio(c('BHM', 'Classic AR(1)'), horizontal=TRUE, container=p3.g1,
								handler=function(h,...){
									val <- svalue(h$obj, index=TRUE)
									enabled(p3.g2a) <- val == 1
									enabled(p3.g2b) <- val == 2
								})
	addSpace(ar.g, 10)
	p3.g2 <- ggroup(horizontal=TRUE, container=ar.g)
	p3.g2a <- ggroup(horizontal=TRUE, container=p3.g2)
	glabel("Burnin:", container=p3.g2a)
	e$burnin3 <- gedit(defaults$burnin3, width=10, container=p3.g2a)
	addSpace(p3.g2, 10)
	p3.g2b <- ggroup(horizontal=FALSE, container=p3.g2)
	p3.g2b1 <- ggroup(horizontal=TRUE, container=p3.g2b)
	glabel("mu:", container=p3.g2b1)
	e$mu <- gedit(defaults$mu, width=10, container=p3.g2b1)
	addSpace(p3.g2b1, 10)
	glabel("rho:<sup>*</sup>", markup=TRUE, container=p3.g2b1)
	e$rho <- gedit(defaults$rho, width=10, container=p3.g2b1)
	addSpace(p3.g2b1, 10)
	glabel("sigma:<sup>*</sup>", markup=TRUE, container=p3.g2b1)
	e$sigmaAR1 <- gedit(defaults$sigmaAR1, width=10, container=p3.g2b1)
	p3.g2b2 <- ggroup(horizontal=TRUE, container=p3.g2b)
	glabel('<sup>*</sup>If left blank values are re-estimated from the data.', markup=TRUE, container=p3.g2b2)
	enabled(p3.g2b) <- FALSE
	addSpace(g, 10)
	.create.status.label(g, e)
	addSpring(g)
	predict.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='tfr.predict', package='bayesTFR', 
				parent.group=predict.g,
						parent.window=main.win)
	addSpring(predict.g)
	create.generate.script.button(handler=run.tfr.prediction, action=list(mw=main.win, env=e, script=TRUE),
								container=predict.g)
	bDem.gbutton(action=gaction(label=' Make Prediction ', icon='evaluate', handler=run.tfr.prediction, 
				action=list(mw=main.win, env=e, script=FALSE)), container=predict.g)
	return(e)

}

.create.prediction.setting.group <- function(g, e, defaults) {
	enable.pred.settings <- function(not.use.diag) {
		enabled(e$burnin) <- not.use.diag
		enabled(e$nr.traj) <- not.use.diag
		enabled(e$thin) <- not.use.diag	
	}
	leftcenter <- c(-1,0)
	pred.g <- gframe("<span color='blue'>Prediction settings</span>", markup=TRUE, horizontal=FALSE, container=g)
	#pred.g1 <- ggroup(horizontal=TRUE, container=pred.g)
	plo <- glayout(container=pred.g)
	plo[1,1, anchor=leftcenter] <- glabel("End year:", container=plo)
	plo[1,2] <- e$end.year <- gedit(defaults$end.year, width=4, container=plo)
	tooltip(e$end.year) <- 'End year of the prediction.'
	l <- 2
	plo[l,1, anchor=leftcenter] <- glabel("Start year:", container=plo)
	plo[l,2] <- e$start.year <- gedit(defaults$start.year, width=4, container=plo)
	tooltip(e$start.year) <- 'Start year of the prediction. Leave empty if the same as present.year in estimation.'
	plo[l,4, anchor=leftcenter] <- 	glabel("Burnin:", container=plo)
	plo[l,5] <- e$burnin <- gedit(defaults$burnin, width=7, container=plo)
	tooltip(e$burnin) <- 'Must be smaller than #iter in MCMCs.'

	l <- l+1
	plo[l,1, anchor=leftcenter] <- glabel("Nr. trajectories:", container=plo)
	plo[l,2] <- e$nr.traj <- gedit(defaults$nr.traj, width=5, container=plo)
	tooltip(e$nr.traj) <- "How many trajectories to generate. Can be alternatively determined by Thin or by settings in Auto simulation."
	plo[l,3, anchor=leftcenter] <- 'OR'
	plo[l,4, anchor=leftcenter] <- glabel("Thin:", container=plo)
	plo[l,5] <- e$thin <- gedit(defaults$thin, width=5, container=plo)
	tooltip(e$thin) <- "Determines how many trajectories to generate. Can be alternatively determined by Nr. trajectories or by settings in Auto simulation."
	l <- l+1
	plo[l,1, anchor=leftcenter] <- glabel("Nr. ascii trajectories:", container=plo)
	plo[l,2] <- e$save.as.ascii <- gedit(defaults$save.as.ascii, width=5, container=plo)
	tooltip(e$save.as.ascii) <- 'Number of trajectories to be exported into an ascii file. Set 0 if no export is desired.'
	plo[l,4, anchor=leftcenter] <- glabel("RNG seed:", container=plo)
	plo[l,5] <- e$seed <- gedit(defaults$seed, width=4, container=plo)
	
	l <- 3
	plo[l,6:7] <- e$use.diagnostics <- gcheckbox("Use diagnostics", checked = defaults$use.diagnostics, 
									handler=function(h, ...) enable.pred.settings(!svalue(h$obj)), 
									container=plo)
	tooltip(e$use.diagnostics) <- "Check this if an 'auto' simulation or convergence diagnostics completed successfully. Nr. trajectories and Thin will be ignored - Settings is taken from converged MCMCs."
	plo[l+1,6:7] <- e$verbose <- gcheckbox("Verbose", checked=defaults$verbose, container=plo)
	enable.pred.settings(!defaults$use.diagnostics)
	return(plo)
}

run.tfr.prediction <- function(h, ...)
{
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory',
									end.year='End year', burnin='Burnin'), env=e)) return()
	param.names <- list(numeric=c('mu', 'rho', 'end.year', 'start.year', 'burnin', 'seed', 'nr.traj', 'thin', 'burnin3'), 
						numvector=c('sigmaAR1'),
						text=c('sim.dir'),
						logical=c('verbose', 'use.diagnostics', 'use.correlation'),
						numtext=c('save.as.ascii') #can be both - numeric or text
						)
	params <- get.parameters(param.names, e, h$action$script)
	params[['use.tfr3']] <- svalue(e$p3.type, index=TRUE)==1
	if(params$use.tfr3) 
		params$rho <- params$sigmaAR1 <- params$mu <- NULL
	else {
		if(is.null(params$rho)) params$rho <- NA
		if(is.null(params$sigmaAR1)) params$sigmaAR1 <- NA
		params$burnin3 <- NULL
	}
	if(params$use.diagnostics) {
		params[['burnin']] <- NULL
		params[['thin']] <- NULL
		params[['nr.traj']] <- NULL
	}
	
	simdir <- get.parameters(list(text='sim.dir'), e, quote=FALSE)$sim.dir # to avoid double quoting if script is TRUE
	if(has.tfr.prediction(sim.dir=simdir)) {
		params[['replace.output']] <- FALSE
		if (gconfirm(paste('Prediction for', simdir, 
								'already exists.\nDo you want to overwrite existing results?'),
				icon='question', parent=h$action$mw))
			params[['replace.output']] <- TRUE
		else return(NULL)
	}

	if (h$action$script) {
		cmd <- paste('tfr.predict(', assemble.arguments(params), ')', sep=' ')
		create.script.widget(cmd, h$action$mw, package="bayesTFR")
	} else {
		if(params$use.diagnostics) {
			mcmc.set <- get.tfr.mcmc(params[['sim.dir']])
			diag.list <- get.tfr.convergence.all(mcmc.set$meta$output.dir)
			if(length(diag.list) == 0) {
				gmessage(paste('There is no diagnostics available for', params[['sim.dir']],
							'. Use manual settings for Nr. of trajectories or Thin.'))
				return()
			}
		}
		.run.prediction(e, type='TFR', handler=get.tfr.prediction.status, option='bDem.TFRpred', 
								call='tfr.predict', params=params, 
								sim.name='TFR prediction', main.win=h$action$mw,
								action=list(sb=e$statuslabel),
								interval=1000)
	}
}

.run.prediction <- function(type='TFR', ...) {
	statusopt <- paste('bDem.', type, 'pred.status', sep='')
	opt <- list()
	opt[statusopt] <- list(NULL)
	options(opt)
	res <- .run.simulation(...)
	options(opt)
	return(res)
}
.update.status <- function(sb, status.option, prefix) {
	predstatus <- getOption(status.option, default=NULL)
	status <- prefix
	if(!is.null(predstatus))
		status <- paste(status, predstatus)
	svalue(sb) <- status
}

get.tfr.prediction.status <- function(h, ...) 
	.update.status(h$action$sb, 'bDem.TFRpred.status', 'Running TFR prediction ...')


.create.extra.countries.group <- function(g, e, main.win, defaults) {
	addSpace(g, 10)
	countries.g <- gframe("<span color='blue'>Countries/Regions selection</span>", markup=TRUE, 
							horizontal=FALSE, container=g)
	countries.g1 <- glayout(container=countries.g)
	countries.g1[1,1] <- e$all.countries <- gcheckbox("All without prediction", checked=TRUE, container=countries.g1,
									handler=function(h,...){
										enabled(e$e.countries.gb) <- !svalue(h$obj)
										})
	countries.g1[1,2] <- '      '
	countries.g1[1,3] <- e$e.countries.gb <- bDem.gbutton("  Select specific countries/regions  ", container=countries.g1,
				handler=selectCountryMenuPred,
				action=list(mw=main.win, env=e, not.predicted=TRUE, multiple=TRUE, sorted=FALSE, 
							type=e$pred.type, label.widget.name='extra.pred.country.label'))
	enabled(e$e.countries.gb) <- !svalue(e$all.countries)
	countries.g1[2,3, anchor=c(-1,0)] <- e$extra.pred.country.label <- glabel('', container=countries.g1)
	addSpace(g, 10)
	sim.g <- gframe("<span color='blue'>Output</span>", markup=TRUE, horizontal=FALSE, container=g)
	sim.g2 <- ggroup(horizontal=TRUE, container=sim.g)
	glabel("Nr. of ascii trajectories:", container=sim.g2)
	e$save.as.ascii <- gedit(defaults$save.as.ascii, width=5, container=sim.g2)
	tooltip(e$save.as.ascii) <- 'Number of trajectories to be exported into an ascii file. Set 0 if no export is desired.'
	addSpace(sim.g2, 30)
	e$verbose <- gcheckbox("Verbose", checked=defaults$verbose, container=sim.g2)
	addSpace(g, 10)
	.create.status.label(g, e)
}
	
tfr.pred.extra.countries.group <- function(g, main.win, parent) {
	e <- new.env()
	defaults <- formals(tfr.predict.extra) # default argument values
	e$sim.dir <- parent$sim.dir
	e$pred.type <- 'tfr'
	.create.extra.countries.group(g, e, main.win, defaults)
	addSpring(g)
	predict.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='tfr.predict.extra', package='bayesTFR', 
				parent.group=predict.g, parent.window=main.win)
	addSpring(predict.g)
	create.generate.script.button(handler=run.tfr.prediction.extra, action=list(mw=main.win, env=e, script=TRUE),
								container=predict.g)
	bDem.gbutton(action=gaction(label=' Make Prediction ', icon='evaluate', handler=run.tfr.prediction.extra, 
				action=list(mw=main.win, env=e, script=FALSE)), container=predict.g)
	return(e)
		  
}

run.tfr.prediction.extra <- function(h, ...)
{
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	param.names <- list(text=c('sim.dir'),
						logical=c('verbose'),
						numtext=c('save.as.ascii') #can be both - numeric or text
						)
	params <- get.parameters(param.names, e, h$action$script)
	params[['countries']] <- if(svalue(e$all.countries)) NULL else e$selected.countries
	
	if (h$action$script) {
		cmd <- paste('tfr.predict.extra(', assemble.arguments(params), ')',sep=' ')
		create.script.widget(cmd, h$action$mw, package="bayesTFR")
	} else {
		.run.prediction(e, type='TFR', handler=get.tfr.prediction.status, option='bDem.TFRpred', 
								call='tfr.predict.extra', params=params, 
								sim.name='TFR extra prediction', main.win=h$action$mw,
								action=list(sb=e$statuslabel), interval=1000)
	}
}

edit.predictions.group <- function(g, main.win, parent, type='tfr') {
	e <- new.env()
	e$sim.dir <- parent$sim.dir
	addSpace(g, 10)
	e$edit.frame <- gframe("<span color='blue'>Edit Median</span>", markup=TRUE, 
						horizontal=FALSE, container=g, expand=TRUE)
	bDem.gbutton("  Select specific country/region  ", container=e$edit.frame,
				handler=selectCountryMenuPred,
				action=list(mw=main.win, env=e, not.predicted=FALSE, multiple=FALSE,
							edit.median=TRUE, sorted=TRUE, type=type))
	#addSpring(g)
	button.g <- ggroup(horizontal=TRUE, container=g)
	addSpring(button.g)
	e$restoreb <- bDem.gbutton(action=gaction(label='Restore BHM medians', 
				handler=restore.bhm.medians, action=list(mw=main.win, env=e, type=type)), container=button.g)
	enabled(e$restoreb) <- FALSE
	e$applyb <- bDem.gbutton(action=gaction(label='Save', handler=edit.prediction, 
				action=list(mw=main.win, env=e, type=type)), container=button.g)
	enabled(e$applyb) <- FALSE
	
	addSpring(g)
	convert.frame <- gframe("<span color='blue'>Convert trajectories to ASCII</span>", markup=TRUE, 
						horizontal=TRUE, container=g)
	glabel("Nr. of ascii trajectories:", container=convert.frame)
	e$save.as.ascii <- gedit(formals(paste(type,'.predict', sep=''))$save.as.ascii, width=5, container=convert.frame)
	e$write.summary <- gcheckbox("Write summary files", checked=TRUE, container=convert.frame)
	e$verbose <- gcheckbox("Verbose", checked=FALSE, container=convert.frame)
	addSpring(convert.frame)
	bDem.gbutton(action=gaction(label='Convert', handler=.convert.trajectories, 
				action=list(mw=main.win, env=e, type=type)), container=convert.frame)
}

edit.prediction <- function(h, ...) {
	e <- h$action$env
	values <- unlist(e$median.df[1,])
	where.modified <- values != e$medians
	do.call(paste(h$action$type,'.median.set', sep=''), list(e$sim.dir.value, e$edit.country.obj$code, 
					values=values[where.modified], years=as.numeric(names(values)[where.modified])))
	h$action$env$medians <- values
	enabled(e$applyb) <- FALSE
}

restore.bhm.medians <- function(h, ...) {
	e <- h$action$env
	new.pred <- do.call(paste(h$action$type, '.median.shift', sep=''), list(e$sim.dir.value, 
						e$edit.country.obj$code, reset=TRUE))
	data <- .get.data.for.median.editor(new.pred, e$edit.country.obj)
	e$median.df[,] <- data
}

.convert.trajectories <- function(h, ...) {
	e <- h$action$env
	sim.dir <- svalue(e$sim.dir)
	nr.traj <- as.numeric(svalue(e$save.as.ascii))
	do.call(paste('convert.', h$action$type, '.trajectories', sep=''), 
				list(dir=sim.dir, n=nr.traj, verbose=svalue(e$verbose)))
	if(svalue(e$write.summary)) {
		if(h$action$type == 'tfr')
			write.projection.summary(dir=sim.dir)
		else write.e0.projection.summary(dir=sim.dir)
	}
}

get.table.of.countries.from.prediction <- function(sim.dir, not.predicted=TRUE, sorted=TRUE, 
					type='tfr', env=NULL) {
	loc.data.pred <- get.table.of.countries.from.meta(sim.dir, prediction=TRUE, sorted=sorted, pred.type=type, env=env)
	if(is.null(loc.data.pred)) return(NULL)
	if(!not.predicted) return(loc.data.pred)
	loc.data.sim <- get.table.of.countries.from.meta(sim.dir, prediction=FALSE, sorted=sorted, pred.type=type, env=env)
	if(is.null(loc.data.sim)) return(NULL)
	mcmc.set <- do.call(paste('get.', type, '.mcmc', sep=''), list(sim.dir=sim.dir))
	if (get.nr.countries(mcmc.set$meta) <= get.nr.countries.est(mcmc.set$meta)) {
		gmessage("No countries/regions without a prediction available.")
		return(NULL)
	}
	loc.data.sim.extra <- loc.data.sim[(get.nr.countries.est(mcmc.set$meta)+1):(get.nr.countries(mcmc.set$meta)),]
	# find countries without a prediction
	is.predicted <- is.element(loc.data.sim.extra[,'code'], loc.data.pred[,'code'])
	not.predicted.idx <- (1:dim(loc.data.sim.extra)[1])[!is.predicted]
	pr <- rep('yes', dim(loc.data.sim.extra)[1])
	pr[not.predicted.idx] <- 'no'
	loc.data<-cbind(loc.data.sim.extra, predicted=pr)
	if(sorted) {
		ord.idx <- order(loc.data[,'name'])
		loc.data <- loc.data[ord.idx,]
	}
	return(loc.data)
}

.get.data.for.median.editor <- function(pred, country.obj) {
	pred.years <- dimnames(pred$quantiles)[[3]]
	medians <- get.median.from.prediction(pred, country.obj$index, country.obj$code)
	data <- matrix(medians, ncol=length(medians))
	colnames(data) <- pred.years
	rownames(data) <- 'medians'
	return(data)
}

show.median.editor <- function(e, type='tfr'){
	sim.dir <- svalue(e$sim.dir)
	pred <- do.call(paste('get.', type, '.prediction', sep=''), list(sim.dir))
	country.obj <- get.country.object(e$selected.countries, meta=pred$mcmc.set$meta)
	if(!is.null(e$country.label)) svalue(e$country.label) <- country.obj$name
	else {
		addSpring(e$edit.frame)	
		e$country.label <- glabel(country.obj$name, container=e$edit.frame)
	}
	e$sim.dir.value <- sim.dir
	e$edit.country.obj <- country.obj
	data <- .get.data.for.median.editor(pred, country.obj)
	e$medians <- data[1,]
	if(!is.null(e$median.df)) e$median.df[,] <- data
	else {
		f <- function(env) {
				enabled(env$applyb) <- TRUE; enabled(env$restoreb) <- TRUE
				}
		df.view <- makeDFView(data, e$edit.frame, f, e)
		e$median.df <- df.view$model
		#e$median.df <- gdf(data, container=e$edit.frame)
		#addhandlerchanged(e$median.df, handler=function(h, ...) {
		#		enabled(e$applyb) <- TRUE; enabled(e$restoreb) <- TRUE
		#		})
	}
}

selectCountryMenuPred <- function(h, ...) {
	country.selected <- function(h1, ...) {
		h$action$env$selected.countries <- svalue(h$action$env$sel.extra.country.gt)
		if(!is.null(h$action$edit.median) && h$action$edit.median) show.median.editor(h$action$env, h$action$type) 
		visible(h$action$env$extra.country.sel.win) <- FALSE
		if(!is.null(h$action$label.widget.name) && !is.null(h$action$env[[h$action$label.widget.name]])) 
			svalue(h$action$env[[h$action$label.widget.name]]) <- paste(h$action$env$selected.countries, collapse=',')
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
												sorted=h$action$sorted, 
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
		country.table <- get.table.of.countries.from.prediction(sim.dir=sim.dir.used, 
							not.predicted=h$action$not.predicted, sorted=h$action$sorted,
							type=if(is.null(h$action$type)) 'tfr' else h$action$type, env=h$action$env)
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
					expand=TRUE, multiple=h$action$multiple, handler=country.selected)
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
