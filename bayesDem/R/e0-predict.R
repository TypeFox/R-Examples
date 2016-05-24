e0NewPred.group <- function(g, main.win, parent) {
	nb <- bDem.gnotebook(container=g, expand=TRUE)
	all.c.g <- ggroup(label="<span color='#0B6138'>All Countries</span>", markup=TRUE, 
					horizontal=FALSE, expand=TRUE, container=nb)
	pred.all <- e0.pred.all.countries.group(all.c.g, main.win, parent)
	extra.c.g <- ggroup(label="<span color='#0B6138'>Extra Areas &amp; Regions</span>",
					 markup=TRUE, horizontal=FALSE, expand=TRUE, container=nb)
	pred.extra <- e0.pred.extra.countries.group(extra.c.g, main.win, parent)
	joint.male.g <- ggroup(label="<span color='#0B6138'>Joint Male</span>",
					 markup=TRUE, horizontal=FALSE, expand=TRUE, container=nb)
	joint.male <- e0.joint.male.group(joint.male.g, main.win, parent)
	edit.g <- ggroup(label="<span color='#0B6138'>Edit Predictions</span>", markup=TRUE, horizontal=FALSE, expand=TRUE, container=nb)
	pred.edit <- edit.predictions.group(edit.g, main.win, parent, type='e0')
	svalue(nb) <- 1
	return(pred.all)
}

e0.pred.all.countries.group <- function(g, main.win, parent) {
	e <- new.env()
	defaults <- formals(e0.predict) # default argument values
	e$sim.dir <- parent$sim.dir
	addSpace(g, 10)
	lo <- .create.prediction.setting.group(g, e, defaults)
	lo[2,6] <- e$predict.jmale <- gcheckbox("Predict joint male", checked=defaults$predict.jmale, container=lo)
	addSpace(g, 10)
	.create.status.label(g, e)
	addSpring(g)
	predict.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='e0.predict', package='bayesLife', 
				parent.group=predict.g,
						parent.window=main.win)
	addSpring(predict.g)
	create.generate.script.button(handler=run.e0.prediction, action=list(mw=main.win, env=e, script=TRUE),
								container=predict.g)
	bDem.gbutton(action=gaction(label=' Make Prediction ', icon='evaluate', handler=run.e0.prediction, 
				action=list(mw=main.win, env=e, script=FALSE)), container=predict.g)
	return(e)

}

run.e0.prediction <- function(h, ...)
{
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory',
									end.year='End year', burnin='Burnin'), env=e)) return()
	param.names <- list(numeric=c('end.year', 'burnin', 'seed', 'nr.traj', 'thin'), 
						text=c('sim.dir'),
						logical=c('verbose', 'use.diagnostics', 'predict.jmale'),
						numtext=c('save.as.ascii') #can be both - numeric or text
						)
	params <- get.parameters(param.names, e, h$action$script)
	if(params$use.diagnostics) {
		params[['burnin']] <- NULL
		params[['thin']] <- NULL
		params[['nr.traj']] <- NULL
	}

	simdir <- get.parameters(list(text='sim.dir'), e, quote=FALSE)$sim.dir # to avoid double quoting if script is TRUE
	if(has.e0.prediction(sim.dir=simdir)) {
		params[['replace.output']] <- FALSE
		if (gconfirm(paste('Prediction for', simdir, 
								'already exists.\nDo you want to overwrite existing results?'),
				icon='question', parent=h$action$mw))
			params[['replace.output']] <- TRUE
		else return(NULL)
	}

	if (h$action$script) {
		cmd <- paste('e0.predict(', paste(paste(names(params), params, sep='='), collapse=', '),
											 ')', sep=' ')
		create.script.widget(cmd, h$action$mw, package="bayesLife")
	} else {
		if(params$use.diagnostics) {
			mcmc.set <- get.e0.mcmc(params[['sim.dir']])
			diag.list <- get.e0.convergence.all(mcmc.set$meta$output.dir)
			if(length(diag.list) == 0) {
				gmessage(paste('There is no diagnostics available for', params[['sim.dir']],
							'. Use manual settings for Nr. of trajectories or Thin.'))
				return()
			}
		}
		.run.prediction(e, type='e0', handler=get.e0.prediction.status, option='bDem.e0pred', 
								call='e0.predict', params=params, 
								sim.name='e0 prediction', main.win=h$action$mw,
								action=list(sb=e$statuslabel),
								interval=1000)
	}
}

get.e0.prediction.status <- function(h, ...) 
	.update.status(h$action$sb, 'bDem.e0pred.status', 'Running e0 prediction ...')


e0.pred.extra.countries.group <- function(g, main.win, parent) {
	e <- new.env()
	defaults <- formals(e0.predict.extra) # default argument values
	e$sim.dir <- parent$sim.dir
	e$pred.type <- 'e0'
	.create.extra.countries.group(g, e, main.win, defaults)
	addSpring(g)
	predict.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic='e0.predict.extra', package='bayesLife', 
				parent.group=predict.g,
						parent.window=main.win)
	addSpring(predict.g)
	create.generate.script.button(handler=run.e0.prediction.extra, action=list(mw=main.win, env=e, script=TRUE),
								container=predict.g)
	bDem.gbutton(action=gaction(label=' Make Prediction ', icon='evaluate', handler=run.e0.prediction.extra, 
				action=list(mw=main.win, env=e, script=FALSE)), container=predict.g)

	return(e)		  
}

run.e0.prediction.extra <- function(h, ...)
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
		cmd <- paste('e0.predict.extra(', paste(paste(names(params), params, sep='='), collapse=', '),
											 ')',sep=' ')
		create.script.widget(cmd, h$action$mw, package="bayesLife")
	} else {
		.run.prediction(e, type='e0', handler=get.e0.prediction.status, option='bDem.e0pred', 
								call='e0.predict.extra', params=params, 
								sim.name='e0 extra prediction', main.win=h$action$mw,
								interval=1000)
	}
}

e0.joint.male.group <- function(g, main.win, parent) {
	e <- new.env()
	e$sim.dir <- parent$sim.dir
	defaults <- formals(e0.jmale.predict) # default argument values
	defaults.est <- formals(e0.jmale.estimate)
	addSpace(g, 10)
	estim.f <- gframe("<span color='blue'>Estimation settings</span>", markup=TRUE, 
							horizontal=FALSE, container=g)
	g1 <- ggroup(horizontal=TRUE, container=estim.f)
	est.lo <- glayout(container=g1)
	est.lo[1,1] <- ''
	est.lo[1,2] <- 'Estimate DoF'
	est.lo[1,3] <- 'DoF'
	est.lo[1,4] <- 'Min. e0'
	est.lo[2,1] <- 'Equation 1'
	est.lo[2,2, anchor=c(0,0)] <- e$estDof.eq1 <- gcheckbox('', checked=defaults.est$estDof.eq1, container=est.lo)
	est.lo[2,3] <- e$eq1.dof <- gedit(defaults.est$start.eq1$dof, container=est.lo, width=5)
	est.lo[3,1] <- 'Equation 2'
	est.lo[3,2, anchor=c(0,0)] <- e$estDof.eq2 <- gcheckbox('', checked=defaults.est$estDof.eq2, container=est.lo)
	est.lo[3,3] <- e$eq2.dof <- gedit(defaults.est$start.eq2$dof, container=est.lo, width=5)
	est.lo[3,4] <- e$min.e0.eq2 <- gedit(defaults.est$min.e0.eq2, container=est.lo, width=5)
	
	g2 <- ggroup(horizontal=TRUE, container=estim.f)
	glabel("User-defined male e0 file:", container=g2)
	e$my.e0.file <- bDem.gfilebrowse(eval(defaults$my.e0.file), type='open', 
					  width=40, quote=FALSE, container=g2)
	
	addSpace(g, 10)				  
	pred.f <- gframe("<span color='blue'>Prediction settings</span>", markup=TRUE, 
							horizontal=FALSE, container=g)
	g3 <- ggroup(horizontal=TRUE, container=pred.f)
	glabel('Female-Male Gap limits:', container=g3)
	e$gap.lim <- gedit("0, 18", container=g3, width=10)
	addSpace(g3, 20)
	e$verbose <- gcheckbox("Verbose", checked=defaults$verbose, container=g3)
	addSpace(g, 10)
	.create.status.label(g, e)
	addSpring(g)
	predict.g <- ggroup(horizontal=TRUE, container=g)
	create.help.button(topic=c('e0.jmale.predict', 'e0.jmale.estimate'), package='bayesLife', 
				parent.group=predict.g, parent.window=main.win)
	addSpring(predict.g)
	create.generate.script.button(handler=joint.male.prediction, 
				action=list(mw=main.win, env=e, script=TRUE), container=predict.g)
	bDem.gbutton(action=gaction(label=' Make Prediction ', icon='evaluate', handler=joint.male.prediction, 
				action=list(mw=main.win, env=e, script=FALSE)), container=predict.g)
	return(e)
}

joint.male.prediction <- function(h, ...)
{
	e <- h$action$env
	if(!has.required.arguments(list(sim.dir='Simulation directory'), env=e)) return()
	param.names <- list(numeric=c('eq1.dof', 'eq2.dof', 'min.e0.eq2'), 
						numvector=c('gap.lim'),
						text=c('sim.dir', 'my.e0.file'),
						logical=c('verbose', 'estDof.eq1', 'estDof.eq2')
						)
	params <- get.parameters(param.names, e, h$action$script)
	if(!is.null(params$eq1.dof)) {
		params$start.eq1 <- list(dof=params$eq1.dof)
		params$eq1.dof <- NULL
	}
	if(!is.null(params$eq2.dof)) {
		params$start.eq2 <- list(dof=params$eq2.dof)
		params$eq2.dof <- NULL
	}
	
	if (h$action$script) {
		cmd <- paste('pred <- get.e0.prediction(sim.dir=', params$sim.dir, 
						')\n', sep='')
		params$sim.dir <- NULL
		cmd <- paste(cmd, 'e0.jmale.predict(pred,', paste(paste(names(params), params, sep='='), collapse=', '),
											 ')',sep=' ')
		create.script.widget(cmd, h$action$mw, package="bayesLife")
	} else {
		pred <- get.e0.prediction(params$sim.dir)
		params$sim.dir <- NULL
		.run.prediction(e, type='e0', handler=get.e0.prediction.status, option='bDem.e0pred', 
								call='e0.jmale.predict', params=c(list(e0.pred=pred), params), 
								sim.name='e0 joint male prediction', main.win=h$action$mw,
								action=list(sb=e$statuslabel),
								interval=1000)
	}
}