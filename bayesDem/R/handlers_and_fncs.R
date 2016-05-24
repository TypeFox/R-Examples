selectDir <- function(h, ...) {
	text.widget <- h$action$text.widget
	gf <- gfile(inititalfilename=svalue(text.widget), 
			type='selectdir',
			handler=function(h2,...){svalue(text.widget)<-h2$file;
				dispose(gf)})
	
}

set.widget.basecolor <- function(widget, color, state="normal") {
	gtk.widget <- getToolkitWidget(widget)
	gtk.widget$modifyBase(state, color)
} 

set.widget.bgcolor <- function(widget, color, state="normal") {
	gtk.widget <- getToolkitWidget(widget)
	gtk.widget$modifyBg(state, color)
} 

create.help.button <- function(topic, package, parent.group, parent.window) {
	helpaction <- gaction(label='Help', icon='help', 
		handler=function(h, ...) {
			oldhtmloption <- options(htmlhelp=FALSE);
			oldchmoption <- options(chmhelp=FALSE);
			ltop <- length(topic)
			gh <- ghelp(topic=topic[ltop], package=package, 
					container=gwindow(parent=parent.window, width=700, height=700));
			if(ltop > 1) {
				for(i in (ltop-1):1) # create the pages in reverse order in order to have focus on the first one
					add(gh, list(topic=topic[i], package=package))
			}
			options(htmlhelp=oldhtmloption);
			options(chmhelp=oldchmoption)
		})
	bDem.gbutton(action=helpaction, container=parent.group)
}

create.generate.script.button <- function(handler, action, container) {
	bDem.gbutton(action=gaction(label=' Generate Script ', icon='justify-fill', handler=handler,
							action=action), container=container)
}

bDem.gwindow <- function(...) {
	win <- gwindow(...)
	set.widget.bgcolor(win, color.main)
	return(win)
}

bDem.gdroplist <- function(items, ...) {
	dl <- gdroplist(items, ...)
	cr <- gtkCellRendererCombo()
	model <- rGtkDataFrame(items)
	cr['model'] <- model
	cr['text-column'] <- 0
	cr['cell-background'] <- color.button
	combo <- getToolkitWidget(dl)
	combo$clear()
	combo$packStart(cr)
	combo$addAttribute(cr, 'text', 0)
	#combo$setActive(selected)
	return(dl)
}

bDem.gfilebrowse <- function(...) {
	fb <- gfilebrowse(...)
	get.filebrowse.button(fb)$modifyBg("normal", color.button)
	return(fb)
}

bDem.gbutton <- function(...) {
	b <- gbutton(...)
	set.widget.bgcolor(b, color.button)
	set.widget.bgcolor(b, color.button, state='insensitive')
	return(b)
}

bDem.gnotebook <- function(...) {
	nb <- gnotebook(...)
	set.widget.bgcolor(nb, color.main)
	return(nb)
}


create.sim.dir.widget <- function(env, parent, main.win, default, dir.widget.name='sim.dir', ...) {
	sim.g <- ggroup(horizontal=TRUE, container=parent)
	glabel("Simulation directory:", container=sim.g)
	#glabel("<span color='red'>*</span>", markup=TRUE, container=sim.g)
	env[[dir.widget.name]] <- bDem.gfilebrowse(default, type='selectdir', width=40, quote=FALSE, container=sim.g)
	create.info.button(dir.widget.name, sim.g, main.win, env, ...)
}

get.filebrowse.button <- function(widget) {
	# got this hack from John Verzani
	return(widget@widget@block@widget@widget[[1]][[2]][[1]])
}

create.info.button <- function(dir.widget.name, parent.group, parent.window, env, type, no.mcmc=FALSE) {
	infoaction <- gaction(label='Info', icon='info', handler=show.summary,
					action=list(mw=parent.window, env=env, dir.widget.name=dir.widget.name,
								type=type, no.mcmc=no.mcmc))
	bDem.gbutton(action=infoaction, container=parent.group)
}

create.graphics.window <- function(parent, title='', dpi=80, ps=10, ...) {
	e <- new.env()
	win <- gwindow(title, parent=parent, horizontal=FALSE)
	g <- ggroup(container=win, horizontal=FALSE, expand=TRUE)
	g1 <- ggroup(container=g, horizontal=TRUE)
	glabel("Output type:", container=g1)
	#types <- formals(savePlot)$type
	e$type <- gdroplist(c("pdf", "postscript", "png", "jpeg", "tiff", "bmp"), container=g1)
	gb <- gbutton('Save', container=g1)
	g2 <- ggroup(container=g, horizontal=TRUE, expand=TRUE)
	ggraphics(container=g2, ps=ps, dpi=dpi, ...)
	addHandlerClicked(gb, handler=saveGraph, action=list(mw=win, env=e, dpi=dpi, dev=dev.cur()))
	Sys.sleep(1)
	return(g)
}

create.script.widget <- function(script, parent, package, title=paste(package, 'commands')) {
	script.widget <- gwindow(title, parent=parent, visible=FALSE, width=600, height=150)
	set.widget.bgcolor(script.widget, "white")
	gt <- gtext("", container=script.widget)
	insert(gt, script)
	#gtext(script, container=script.widget)
	addHandlerFocus(script.widget, handler=function(h, ...) focus(gt) <- TRUE)
	visible(script.widget) <- TRUE
	focus(script.widget) <- TRUE
}

saveGraph <- function(h, ...){
	e <- h$action$env
	type <- svalue(e$type)
	postfix <- list(png='png', jpeg='jpg', pdf='pdf', tiff='tiff', bmp='bmp', postscript='ps')
	filter <- list()
	filter[[paste(type, 'files')]] <- list(patterns=paste('*', postfix[[type]], sep='.'))
	gfile(type='save', quote=FALSE,
					filter=c(filter, list("All files" = list(patterns = c("*")))),
					handler=function(h1,...){
						filename <- h1$file
						pattern <- paste('[.]', postfix[[type]], '$|[.]', type, '$', sep='')
						if(length(grep(pattern, filename))==0)
							filename <- paste(filename,  postfix[[type]], sep='.')
						size <- list()
						for(measure in c('width', 'height')) {
							size[[measure]]  <- NULL
							if(is.null(e[[measure]])) {
								if(type != 'pdf' & type != 'postscript') {
									size[[measure]] <- h$action$dpi*6
								}
							} else {
								if(e[[measure]] != 'default') {
									if(!is.null(e[[measure]][[type]]))
										size[[measure]] <- e[[measure]][[type]]
									else
										size[[measure]] <- e[[measure]]
								}
							}
						}
						dev.set(h$action$dev)
						do.call('dev.print', c(list(file=filename, device=eval(parse(text=type))),
												size))
					}
			)
	
}

show.summary <- function(h, ...) {
	e <- h$action$env
	dir <- svalue(e[[h$action$dir.widget.name]])
	type <- h$action$type
	warn <- getOption('warn')
	options(warn=-1) # disable warning messages
	mcmc.sets <- list()
	if(!h$action$no.mcmc) {
		for(ityp in 1:length(type)) {
			mcmc.sets[[ityp]] <- do.call(paste('get.', type[ityp], '.mcmc', sep=''), list(dir))
		}
	} 
	# get prediction
	pred <- do.call(paste('get.', type[1], '.prediction', sep=''), list(sim.dir=dir))
	options(warn=warn)
	info <- c()
	con <- textConnection("info", "w", local=TRUE)
	sink(con)
	if(!h$action$no.mcmc) {
		if (is.null(mcmc.sets[[1]])) {
			cat('No simulation available in this directory.')
		} else { 
			cat('Simulation results\n')
			cat('********************\n')
			for(i in 1:length(mcmc.sets)) {
				print(summary(mcmc.sets[[i]], meta.only=TRUE))
				cat('===============================\n')
			}
		}
	}
	if (is.null(pred)) {
		cat('\nNo projections available in this directory.\n')	} else {	
		cat('\nProjections')
		if(!is.null(pred$end.year)) cat(' for end year:', pred$end.year) else cat(':')
		cat('\n------------------')
		print(summary(pred))
	}
	sink()
	close(con)
	info.win <- gwindow('Directory Info', parent=h$action$mw, width=500, height=400)
	gt <- gtext("", container=info.win)
	insert(gt, info)
	#gtext(info, container=info.win)
}

get.parameters <- function(par.names, env, quote=FALSE, retrieve.from.widgets=TRUE) {
	# Get values from a bunch of widgets, separated into groups by the value types
	# par.names is a list of vectors with elements: 'numeric', 'numvector', 'text', 'numtext' etc.
	# Each vector contains the names of the corresponding widgets in the environment 'env'. 
	params <- list()
	op <- options("useFancyQuotes")
	options(useFancyQuotes = FALSE)
	for (par in unlist(par.names)) {
		params[[par]] <- if(retrieve.from.widgets) svalue(env[[par]]) else env[[par]]
		if(is.null(params[[par]])) next
		if (nchar(params[[par]])==0) {params[[par]] <- NULL; next}
		if(any(par == par.names$numeric)) params[[par]] <- as.numeric(params[[par]])
		if(any(par == par.names$numvector)) params[[par]] <- as.numeric(strsplit(params[[par]], ',')[[1]])
		if (quote & any(par == par.names$text)) params[[par]] <- sQuote(params[[par]])
		# convert back- to forward-slashes in path names
		if(any(par == par.names$text)) params[[par]] <- gsub("\\\\", "/", params[[par]]) 
		if(any(par == par.names$numtext)) {
			warn <- getOption('warn')
			options(warn=-1)
			numpar <- as.numeric(params[[par]])
			options(warn=warn)
			if (!is.na(numpar)) {params[[par]] <- numpar}
			else {if (quote) params[[par]] <- sQuote(params[[par]])}
		}
	}
	options(op)
	return(params)
}

assemble.arguments <- function(params, add.params=NULL) {
	if(is.null(params)) return ('')
	argline <- paste(paste(names(params), params, sep='='), collapse=', ')
	if(!is.null(add.params) && nchar(add.params) > 0)
		argline <- paste(argline, ',', add.params)
	return(argline)
}

has.required.arguments <- function(par.names, env) {
	for (par in names(par.names)) {
		value <- svalue(env[[par]])
		if (nchar(value)==0) {
			gmessage(paste('Argument', par.names[[par]], 'is required.'))
			return(FALSE)
		}
	}
	return(TRUE)
}

makeDFView <- function(df, container, f=NULL, ...) {
	# from John Verzani
	model <- rGtkDataFrame(df)
	view <- gtkTreeView(model)
	## Michael Lawrence's trick
	mapply(view$insertColumnWithAttributes,-1, colnames(model), 
		lapply(1:ncol(model), function(i) gtkCellRendererText()), text = seq_len(ncol(model)) - 1)

	sw <- gtkScrolledWindow()
	sw$add(view)
	
	coerceVar <- function(x, value) UseMethod("coerceVar")
	coerceVar.default <- function(x, value) value
	coerceVar.numeric <- function(x, value) as.numeric(value)
	coerceVar.integer <- function(x, value) as.integer(value)
	coerceVar.logical <- function(x, value) as.logical(value)
	coerceVar.factor <- function(x, value) {
		ind <- pmatch(value, levels(x))
		if(is.na(ind))
			ind
		else
			levels(x)[ind]
	}

	sapply(1:ncol(model), function(j) {
		col <- view$getColumn(j-1)
		#print(col)
		cr <- col$getCellRenderers()[[1]]
		cr['editable'] <- TRUE
		gSignalConnect(cr, "edited", 
			f=function(cr, path, newtext, user.data) {
				if(!is.null(f)) f(...) else NULL
				curRow <- as.numeric(path) + 1
				curCol <- user.data$column
				model <- user.data$model
				## coerce newtext from character to desired type
				## otherwise this coerces to character
				model[curRow, curCol] <- coerceVar(model[,curCol], newtext)
			}, data=list(model=model, column=j))
	})

	add(container, sw, expand=TRUE)
	list(model=model, view=view)
}

get.data.path <- function(type) {
	package <- NULL
	if(type=='tfr') package <- 'bayesTFR'
	else{if(type=='e0') package <- 'bayesLife'}
	if(is.null(package)) stop('Wrong type ', type, '. Only "tfr" and "e0" allowed.')
	return(file.path(find.package(package), "data"))
}

get.tfr.UN.data <- function(meta) {
	return(bayesTFR:::load.bdem.dataset('tfr', meta$wpp.year))
}

get.e0.UN.data <- function(meta) {
	return(bayesTFR:::load.bdem.dataset(paste('e0', meta$sex, sep=''), meta$wpp.year))
}
