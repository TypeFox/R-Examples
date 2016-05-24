bayesDem.go <- function(wpp.year.tfr=wpp.year.default, wpp.year.e0=wpp.year.tfr,
						wpp.year.pop=wpp.year.tfr) {
			
	quit.bayesdem <- function(h, ...) {
		dispose(main.win)
		#detach(gtkdb)
	}
	options(guiToolkit=guiToolkit.default)
	# 
	path = system.file("images",package="bayesDem")
	wait.window <- gwindow('Bayesian Demographer Initialization', width=400, height=100,
						parent=c(500, 300), visible=FALSE)
	set.widget.bgcolor(wait.window, "white")
	#glabel('Starting Bayesian Demographer ...', container=wait.window)
	gimage(file.path(path, 'startpyramid.png'), , container=wait.window)

	visible(wait.window) <- TRUE
	
	# main window
	main.win <<- gwindow(paste('Bayesian Demographer  v.', 
			packageVersion("bayesDem")), visible=FALSE, parent=c(400,50))
	main.g <- ggroup(horizontal=FALSE, container=main.win)
	
	# notebook with tabs
	main.notebook <- bDem.gnotebook(container=main.g, expand=TRUE)

	# TFR prediction tab
	tfr.pred <- ggroup(label="<span weight='bold' color='blue'>Projection of Total Fertility Rate</span>", 
		markup=TRUE, horizontal=FALSE, container=main.notebook)

	tfrPredTab(tfr.pred, main.win, wpp.year=wpp.year.tfr)

	# Life expectancy
	e0w <- ggroup(label="<span weight='bold' color='blue'>Projection of Life Expectancy</span>", 
		markup=TRUE, horizontal=FALSE, container=main.notebook)
	e0PredTab(e0w, main.win, wpp.year=wpp.year.e0)
	
	# Population Prediction
	pop.w <- ggroup(label="<span weight='bold' color='blue'>Population Projection</span>", 
		markup=TRUE, horizontal=FALSE, container=main.notebook)
	popPredTab(pop.w, main.win, wpp.year=wpp.year.pop)

	svalue(main.notebook)<- 1
	
	# Quit Button
	button.group <- ggroup(container=main.g, expand=TRUE)
	gbutton('Quit', handler=quit.bayesdem, container=button.group)
	addSpring(button.group)
	label <- glabel('BayesPop group\nUniversity of Washington', container=button.group)
	font(label) <- c(style='normal', family='serif', size='xx-small')
	#addHandlerRightclick(main.win, handler=function(h, ...) focus(h$obj) <- TRUE)
	addHandlerFocus(main.win, handler=function(h, ...) focus(h$obj) <- TRUE)
	#addHandlerDragMotion(main.win, handler=function(h, ...) focus(h$obj) <- TRUE)
	dispose(wait.window)
	visible(main.win) <- TRUE
	# This is necessary in order for bayesTFR etc to find RGtk2 functions
	# gtkdb <- list(gtkEventsPending=RGtk2::gtkEventsPending, gtkMainIteration=RGtk2::gtkMainIteration)
	# attach(gtkdb)
	do.call('require', list('RGtk2'))
}

