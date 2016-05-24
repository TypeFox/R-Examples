navGraph <- function(data, graph = NULL, viz = NULL, settings = NULL) {
	## omit no visible binding note
	t.vizcounter <- NULL
	g <- NULL
	Glist <- NULL


	##
	## removed GGobi
	##
	## Start Rgobi if there is an NG_Scatterplot object in the vizList
	if(!is.null(settings) && is.list(settings) && !is.null(settings$defaultDisplay)){
		if(tolower(settings$defaultDisplay)=="ggobi"){
			stop("navGraph: as from RnavGraph version 0.1.6 ggobi is not supported anymore!")					
		}
	}



	##
	## RnavGraph can be started by
	## - either passing by a navGraph handler
	## - or by passing by: graph, viz and settings objects
	## - or by passing only by a data object
	## - or by passing data and variable graph
	###################################################################
	if(identical(data, 'tclreset')) {
		cat('Currently disabled. Restart your R session please to get the same effect.\n\n')
		##cat('Note that any running navGraph session will not work anymore.\n')
		##.Tcl('set ng_data ""')
		##.Tcl('set ng_windowManager ""')
		
	}else if(is(data, "NavGraph_handler")){
		
		## if navgraph handler gets passed by
		graphList <- data@graphs
		dataList <- data@data
		
		settings <- data@settings
		paths <- data@paths
		activePath <- tclVar(data@activePath)
		activePathGraph <- tclVar(data@activePathGraph)
		tk2dcolors <- data@tk2dcolors
		
		## take image ids away
		vizList <- sapply(data@viz,FUN=function(viz){
					if(is(viz,"NG_Viztk2d")){
						if(!is.null(viz@imgIds)){
							warning(paste("[navGraph] restarting session with navGraph handler. Images from ",viz@viz_name," visualization instructions won't show."))								
						}
						viz@imgIds <- NULL
					}
					return(viz)
				})
	}else {
		tk2dcolors <- list()
		
		in_args <- c(!is.null(data),!is.null(graph),!is.null(viz))
		if( all(in_args)) {
			## data, graph and viz objects were passed by
			
			##  First check whether graph, data and viz are from the correct format.
			##  The rest of navGraph works with objects called
			##  graphList, dataList, vizList
			
			for(arg in c("graph", "data", "viz")){
				if(arg == "viz"){
					arg_class <- "NG_Visualization"
				}else{
					arg_class <- paste("NG_",arg,sep = '')
				}	
				if(is(get(arg),arg_class)){
					assign(paste(arg,"List", sep = ''),list(get(arg)))	
				} else if(is.list(get(arg))){
					if(all(sapply(get(arg),function(x){is(x,arg_class)}))){
						assign(paste(arg,"List", sep = ''),get(arg))	
					}else{
						stop(paste("[navGraph] elements in list ",arg," are not from class ",arg_class, sep = ''))
					}
				} else {
					stop(paste("[navGraph] argument ",arg," is neither a list of NG_", arg," objects nor a NG_", arg," object",sep = ''))
				}	
			}
		} else if((in_args[1] == TRUE) && (in_args[2] == TRUE) && (in_args[3] == FALSE)) {
			## Data and variable graph
			n <- length(data)
			
			if(is.list(data)) {
				dataList <- data
			} else {
				dataList <- list(data)
			}
			if(is.list(graph)) {
				GList <- graph
			} else {
				GList <- list(graph)
			}
			if(!all(sapply(dataList, FUN=function(x){is(x,"NG_data")}))) {
				stop("[navGraph] elements in list data are not from class NG_data")
			}
			
			if(!all(sapply(GList, FUN=function(x){is(x,"graph")}))) {
				stop("[navGraph] elements in list graph are not from class graph")
			}
			
			
			if((length(graph) > 1) && (n == 1)) {
				## multiple graphs and one data set
				## check that all graphs relate to data set
				dnames <- names(data) 
				snames <- shortnames(data)
				
				graphMatchdata <- sapply(GList, FUN=function(G){
							all(nodes(G) %in% dnames) || all(nodes(G) %in% snames)
						})
				if(!all(graphMatchdata)) {
					stop("[navGraph] you passed on graphs that don't match your data.")
				}
				graphList <- vector("list", length = 2*n)
				vizList <- vector("list", length = 2*n)
				
				if(is.null(names(GList))){
					graph_names <- paste('g',1:length(Glist), sep = '')
				} else {
					graph_names <- names(GList)
				}
				sep = ":"
				for(i in 1:length(GList)) {
					G = GList[[i]]
					if(any(grepl(sep, nodes(G), fixed = TRUE))) {
						notanswered <- TRUE
						while(notanswered){
							sep <- readline(paste("Seperator '",sep,"' exists in node names. Choose different separator: ", sep=''))
							if(sep == '' || any(grepl(sep, nodes(G), fixed = TRUE))){
								sep <- readline(paste("Seperator '",sep,"' exists in node names. Choose different separator: ", sep=''))
							}else {
								notanswered <- FALSE
							}
						}
					}
					ii <- (i-1)*2+1
					
					LG <- linegraph(G, sep = sep)
					LGnot <- complement(LG)
					graphList[[ii]] <- ng_graph(paste(graph_names[i],"3d",data@name),LG,sep = sep)
					graphList[[ii+1]] <- ng_graph(paste(graph_names[i],"4d",data@name),LGnot,sep = sep)
					
					vizList[[ii]] <- ng_2d(data, graphList[[ii]])
					vizList[[ii+1]] <- ng_2d(data, graphList[[ii+1]])
				}
				
				
			} else if(n == length(graph)) {
				## one graph per one data
				graphList <- vector("list", length = 2*n)
				vizList <- vector("list", length = 2*n)
				
				sep = ":"
				for(i in 0:(n-1)) {
					G = GList[[i+1]]
					if(any(grepl(sep, nodes(G), fixed = TRUE))) {
						notanswered <- TRUE
						while(notanswered){
							sep <- readline(paste("Seperator '",sep,"' exists in node names. Choose different separator: ", sep=''))
							if(sep == '' || any(grepl(sep, nodes(G), fixed = TRUE))){
								sep <- readline(paste("Seperator '",sep,"' exists in node names. Choose different separator: ", sep=''))
							}else {
								notanswered <- FALSE
							}
						}
					}
					ii <- i*2+1
					
					LG <- linegraph(G, sep = sep)
					LGnot <- complement(LG)
					graphList[[ii]] <- ng_graph(paste("3d",dataList[[i+1]]@name),LG,sep = sep)
					graphList[[ii+1]] <- ng_graph(paste("4d",dataList[[i+1]]@name),LGnot,sep = sep)
					
					vizList[[ii]] <- ng_2d(dataList[[i+1]], graphList[[ii]])
					vizList[[ii+1]] <- ng_2d(dataList[[i+1]], graphList[[ii+1]])
				}
			} else {
				stop("[navGraph] if only data and (variable!)graph arguments are submitted, then they must be of the same length.")				
			}
			
			
			
			
		}	else if( (in_args[1] == TRUE) && (in_args[2] == FALSE) && (in_args[3] == FALSE) ) {
			## only look at data object
			if(!is.list(data)) {
				dataList <- list(data)
			} else {
				dataList <- data
			}
			if(!all(sapply(dataList, FUN=function(x){is(x,"NG_data")}))) {
				stop("[navGraph] some elements in list data are nGlistot from class NG_data")
			}
			
	
			
			
			out <- unlist(lapply(dataList, FUN=function(data){
								V <- shortnames(data)
								if(length(V) == 0) {
									V <- names(data)
								}
								G <- completegraph(V)
								LG <- linegraph(G)
								LGnot <- complement(LG)
								ng.lg <- ng_graph(name = paste(data@name,': 3D'), graph = LG, layout = 'circle')
								ng.lgnot <- ng_graph(name = paste(data@name,': 4D'), graph = LGnot, layout = 'circle')
								
								viz1 <- ng_2d(data,ng.lg)
								viz2 <- ng_2d(data,ng.lgnot)
			
								return(list(ng.lg,ng.lgnot,viz1,viz2))
							}))
			
			is.graph <- sapply(out, FUN=function(element){
						is(element,"NG_graph")
					})
			
			graphList <- out[is.graph]
			vizList <- out[!is.graph]
		} else {
			## weird combination of arguments passed by	
			stop("[navGraph] if the graph argument is used, then viz argument must be used too and vice versa.")
		}
		
		
		
		## check if data and graph names are unique
		data_names <- sapply(dataList,FUN=function(x)x@name)
		graph_names <- sapply(graphList,FUN=function(x)x@name)
		
		if(length(unique(data_names)) != length(data_names)) {
			stop("[navGraph] the names of the NG_data objects passed by are not unique")
		}
		if(length(unique(graph_names)) != length(graph_names)) {
			stop("[navGraph] the names of the NG_graph objects passed by are not unique")
		}
		
		
		## check if the objects in the visualization instruction (data and graph) have been passed by
		if(!all(sapply(vizList,FUN=function(x)x@graph) %in% graph_names)) {
			stop("[navGraph] some of the visualization instructions contain connect to graphs that werent passed to navGraph")
		}
		
		if(!all(sapply(vizList,FUN=function(x)x@data) %in% data_names)) {
			stop("[navGraph] some of the visualization instructions contain connect to data that werent passed to navGraph")
		}
		
		
		##  create settings object
		##
		##  The program creats an object containing all the program settings. 
		## 	The user can pass on a list with custom settings.
		
		if(is.null(settings) == FALSE){
			if(is.list(settings)==FALSE){
				stop("[navGraph] argument settings in navGraph has to be a list")
			}else{
				if(is.null(settings$color) == FALSE){
					s_col <- do.call('new',append('ColorSettings', settings$color))
				}else{
					s_col <- new('ColorSettings')
				}
				
				if(is.null(settings$interaction) == FALSE){
					s_inter <- do.call('new',append('InteractionSettings', settings$interaction))
				}else{
					s_inter <- new('InteractionSettings')
				}
				
				if(is.null(settings$display) == FALSE){
					s_disp <- do.call('new',append('DisplaySettings', settings$display))
				}else{
					s_disp <- new('DisplaySettings')
				}
				if(is.null(settings$tk2d) == FALSE){
					s_tk2d <- do.call('new',append('Tk2dDisplay', settings$tk2d))
				}else{
					s_tk2d <- new('Tk2dDisplay')
				}
			}
		}else{
			s_col <- new('ColorSettings')
			s_inter <- new('InteractionSettings')
			s_disp <- new('DisplaySettings')
			s_tk2d <- new('Tk2dDisplay')
		}
		
		## Settings object
		settings <- new('NG_Settings', color = s_col,
				interaction = s_inter, display = s_disp, tk2d = s_tk2d)
		
		## Path related stuff
		paths <- new("NG_path")
		activePath <- tclVar("")
		activePathGraph <- tclVar("")
	}		
	
	
	
	
	## loacl variables
	## ---------------
	## I'll often just pass by the ngEnv variable and access the
	## graphs, settings, etc... vie the ngEnv variable
	###################################################################
	
	ngEnv <- environment()		# navGraph enivronment
	
	
	startNode <- nodes(graphList[[1]]@graph)[1]	
	bulletState <- list(from = startNode, to = "", percentage = 0)
	windowHistory <- list(width = 400, height = 400)
	mouseHistory <- list(x=0, y=0, xClick = 0, yClick = 0, Click = FALSE)
	windowManager <- list() ## all windows that are open
	
	
	
	
	## name Data list object
	dataNames <- sapply(dataList,function(x)x@name)
	names(dataList) <- dataNames
	
	
	
	## Graphs & Visualizations
	## -----------------------
	## Graphs are unlisted into the objects: graph1, graph2, ...
	## This makes it easier to access them later
	## ###########################################################
	
	## adj() does not work with graphAM object, thats why I use the coerce function
	nGraphs <- length(graphList)
	for(i in 1:nGraphs){
		slot(graphList[[i]],'graph') <- as(graphList[[i]]@graph,'graphNEL')
		assign(paste('graph',i,sep=''), graphList[[i]], envir = ngEnv)
	}
	cGraph <- 'graph1' ## active graph
	assign('graph', get(cGraph)) ## physical active graph	

	
	
	
	
	## TCLTK
	## ------
	## Here we create the canvas window and define it's bindings
	## 
	###################################################################
	
	
	
	
	
	
	
	## TODO: document linked state in vignette
	if(as.numeric(.Tcl('info exists ng_windowManager("ngInstance")'))) {
		ng_instance <- as.numeric(.Tcl('incr ng_windowManager("ngInstance") 1'))
	} else {
		## create global array windowManager
		.Tcl('set ng_windowManager("ngInstance") 1')
		ng_instance <- 1
	}
	
	if(ngEnv$settings@tk2d@linked) {
		ng_LinkedInstance <- 0
	}else {
		ng_LinkedInstance <- ng_instance
	}
	
	
	
	## Progress bar
	ttp <- tktoplevel()
	tktitle(ttp) <- paste("Session", ng_instance)	
	tkgrab(ttp)
	tkpack(f.progress <- tkframe(ttp, padx = 5, pady = 5))
	
	progress <- tkwidget(f.progress,"ttk::progressbar", orient = 'horizontal', length=200, mode='determinate')
	tkpack(progress, side = "top", anchor = "w")
	progress.label <- tklabel(f.progress, text = "Graph Display")
	tkpack(progress.label, side = "top", anchor = "w")
	
	tkconfigure(progress, value = "10")
	
	
	tt <- tktoplevel()
	tktitle(tt) <- paste("Session ", ng_instance,", RnavGraph Version ", utils::packageDescription("RnavGraph", field="Version"), sep = '')
	
	canvas <- tkcanvas(tt, width = windowHistory$width,
			height = windowHistory$height,
			background = settings@color@background)
	tkpack(canvas, side = "top", fill = "both", expand = TRUE)
	
	tkbind(tt, '<Destroy>', function(){
				.closenavGraph(ngEnv)
			})
	
	
	
	##
	## tk menu
	##
	#####################################
	
	topMenu <- tkmenu(tt)
	tkconfigure(tt,menu=topMenu)
	
	fileMenu <- tkmenu(topMenu, tearoff = FALSE)
	tkadd(topMenu,"cascade",label="File",menu=fileMenu)
	
	tkadd(fileMenu,"command", label="Settings", command=function(){.settingsMenue(ngEnv)})
	tkadd(fileMenu,"command",label="Reset Edges Seen",command=function(){
				tkitemconfigure(ngEnv$canvas,'edge && graph', fill=ngEnv$settings@color@notVisitedEdge)
				tkdtag(ngEnv$canvas, 'all', 'visited')
			})
#	tkadd(fileMenu,"command", label="Save navGraph Handler", command=function(){				
#				ttwin <- tktoplevel()
#				tktitle(ttwin) <- "Save a navGraph handler"
#				loc <- .tcl2str(tkwm.geometry(ngEnv$tt))
#				loc1 <- unlist(strsplit(loc,split = '+',  fixed = TRUE))
#				tkwm.geometry(ttwin,paste("+",as.numeric(loc1[2])+30,"+",as.numeric(loc1[3])+150,sep=""))			
#				tkgrab(ttwin)
#				tkpack(tkframe(ttwin, height = 10),side = "top")
#				l <- tklabel(ttwin, text = "Save the navGraph handler in the global environment as:")
#				tkpack(l,side = "top", anchor = "w", padx = 5)
#				entry <- tkentry(ttwin)
#				tkpack(entry, side = "top", fill = "x",padx = 5, pady = 5)
#				f.b <- tkframe(ttwin)
#				tkpack(f.b, side = "top", padx = 5)
#				tkpack(tkframe(ttwin, height = 5),side = "top")
#				bok <- tkbutton(f.b,text="OK", command=function(){.save()})
#				bcan <- tkbutton(f.b,text="Cancel", command=function(){tkdestroy(ttwin)})
#				tkpack(bok,bcan, side = "left", anchor = "center", padx = 10)
#				tkfocus(entry)
#				tkbind(entry,"<Return>",function(){.save()})
#				.save <- function(){
#					varName <- .tcl2str(tkget(entry))
#					if (length(varName) != 0) {
#						if (varName %in% ls(.GlobalEnv)) {
#							if (tk_messageBox("yesno",paste('Variable "',varName,'" already exists. Overwrite it?',sep=""), parent = ttwin) == "yes") {
#								cat(paste("Session", ngEnv$ng_instance,"navGraph handler saved (overwritten) as:",varName,'\n'))
#								assign(varName, new("NavGraph_handler", env = ngEnv,
#												graphs = graphList, data = dataList, viz = vizList,
#												settings = settings, paths = paths, activePath = tclvalue(activePath),
#												activePathGraph = tclvalue(activePathGraph), dateCreated = date(),
#												dateUpdated = "not"), envir = .GlobalEnv)
#								tkdestroy(ttwin)
#							}
#						} else {
#							cat(paste("Session", ngEnv$ng_instance, "navGraph handler saved as:",varName,'\n'))
#							assign(varName, new("NavGraph_handler", env = ngEnv,
#											graphs = graphList, data = dataList, viz = vizList,
#											settings = settings, paths = paths, activePath = tclvalue(activePath),
#											activePathGraph = tclvalue(activePathGraph), dateCreated = date(),
#											dateUpdated = "not"), envir = .GlobalEnv)
#							tkdestroy(ttwin)
#						}
#					}
#				}
#			})
	
	tkadd(fileMenu,"separator")
	tkadd(fileMenu,"command",label="Quit",command=function(){
				tkdestroy(ngEnv$tt)
			})
	
	
	graphMenu <- tkmenu(topMenu, tearoff = FALSE)
	tkadd(topMenu, "cascade", label = "Graph", menu = graphMenu)
	sapply(1:nGraphs,function(graphNr){.switchGraph(ngEnv, graphNr)})
	
	toolsMenu <- tkmenu(topMenu, tearoff = FALSE)
	tkadd(topMenu, "cascade", label = "Tools", menu = toolsMenu)
	tkadd(toolsMenu,"command",label="Paths",command=function(){.pathGUI(ngEnv)})
	
	## TODO: reset graph layout
#	tkadd(toolsMenu,"separator")
	
#	tkadd(toolsMenu,"command",label="Reset Graph Layout",command=function(){})
	
	
	
	
	
	##
	## tk binding
	##
	#####################################
	tkconfigure(progress.label, text = "Define Bindings")
	tkconfigure(progress, value = "20")
	
	## Buttonclick on every canvas element
	## register it
	tkitembind(canvas, 'all', '<Button-1>', function(x,y) {
				.registerButtonClick(ngEnv,x,y) 
			})
	
	
	## move bullet
	tkitembind(canvas, 'bullet', '<B1-Motion>', function(x,y) {
				.moveBullet(ngEnv,x,y)
			})
	
	## necessay to break after this 
	## (in order to not call Button-1 if control is pressed)
#tkitembind(canvas, 'node', '<Control-Button-1>', function(x,y) {
#			##	.registerButtonClick(ngEnv,x,y)
#			##tcl('break')
#		})
	
	
	
	tkitembind(canvas, 'node', '<Control-B1-Motion>', function(x,y) {
				.moveNode(ngEnv,x,y)
				tcl('break')
			})
	
	tkitembind(canvas, 'label', '<Control-B1-Motion>', function(x,y) {
				.moveLabel(ngEnv,x,y)
				tcl('break')
			})
	
	
	## Resizing the window and canvas objects
	tkbind(tt,'<Configure>', function(W,w,h){
				##w <- as.numeric(w)
				##h <- as.numeric(h)
				## need exact canvs size because of menue
				w <- .tcl2num(tcl('winfo','width',canvas))
				h <- .tcl2num(tcl('winfo','height',canvas))
				xf <- w/ngEnv$windowHistory$width
				yf <- h/ngEnv$windowHistory$height
				tkitemscale(canvas,'edge || label',0,0,xf,yf)
				tcl('scaleNoArea', ngEnv$canvas, 'node || bullet', 0,0,xf,yf)
				ngEnv$windowHistory$width <- w
				ngEnv$windowHistory$height <- h
				
			})
	
	
	## select adjacent node
	tkitembind(canvas,'node || label','<Button-1>', function(){
				
				node <- .tcl2str(tkgettags(ngEnv$canvas,'current'))[3]
				if(!(node %in% c(ngEnv$bulletState$from,ngEnv$bulletState$to))){
					if((node %in% adj(ngEnv$graph@graph,ngEnv$bulletState$from)[[ngEnv$bulletState$from]]) && 
							(ngEnv$bulletState$percentage == 0)){
						## select path						
						ngEnv$bulletState$to <- node
						.leaveNode(ngEnv)
						
					}else {
						## jump to node
						ngEnv$bulletState$from <- node
						ngEnv$bulletState$percentage <- 0 
						.updatePlots(ngEnv)
						.arriveAtNode(ngEnv)
						xynode <- .tcl2xy(tkcoords(ngEnv$canvas,paste('node && ',node)))
						xybullet <- .tcl2xy(tkcoords(ngEnv$canvas,'bullet'))
						
						dxy <- xynode-xybullet
						tkmove(ngEnv$canvas, 'bullet',dxy[1],dxy[2])
					}
				}
			})
	## look that the above event does not get executed with ctr-B-1
	tkitembind(canvas,'node || label','<Control-Button-1>', function(){})
	
	
	## animate edge
	tkitembind(canvas, 'node || label', '<Double-Button-1>', function(){
				node <- .tcl2str(tkgettags(ngEnv$canvas,'current'))[3]
				if(node == ngEnv$bulletState$to) {
					## walk forward
					.walkEdge(ngEnv)					
				}else if(node == ngEnv$bulletState$from){
					if(ngEnv$bulletState$percentage != 0) {
						.walkEdge(ngEnv,TRUE)	
					}
				}
				
			})
	## look that the above event does not get executed with ctr-B-1
	tkitembind(canvas,'node || label','<Control-Double-1>', function(){})
	
	
	
	sysname <- Sys.info()[1]
	if(sysname == "Windows") {
		tkbind(tt, '<MouseWheel>', function(D){
					if(ngEnv$bulletState$to != ''){
						xyfrom <- .tcl2xy(tkcoords(ngEnv$canvas,paste('node &&', ngEnv$bulletState$from)))
						xyto <- .tcl2xy(tkcoords(ngEnv$canvas,paste('node &&', ngEnv$bulletState$to)))
						dxy <- (xyto - xyfrom)*1/ngEnv$settings@interaction@NSteps	
						if(D > 0) {
							if(ngEnv$bulletState$percentage == 0) {
								ngEnv$mouseHistory$Click <- TRUE ## if you want to leave the node							
							}else {
								ngEnv$mouseHistory$Click <- FALSE
							}
							.moveBullet(ngEnv, ngEnv$mouseHistory$x + dxy[1], ngEnv$mouseHistory$y + dxy[2])						
						}else {
							.moveBullet(ngEnv, ngEnv$mouseHistory$x - dxy[1], ngEnv$mouseHistory$y - dxy[2])						
						}
					}
				})
	} else {
		## scroll up
		tkbind(tt, '<Button-4>', function(){
					if(ngEnv$bulletState$to != ''){
						if(ngEnv$bulletState$percentage == 0) {
							ngEnv$mouseHistory$Click <- TRUE ## if you want to leave the node							
						}else {
							ngEnv$mouseHistory$Click <- FALSE
						}
						xyfrom <- .tcl2xy(tkcoords(ngEnv$canvas,paste('node &&', ngEnv$bulletState$from)))
						xyto <- .tcl2xy(tkcoords(ngEnv$canvas,paste('node &&', ngEnv$bulletState$to)))
						dxy <- (xyto - xyfrom)*1/ngEnv$settings@interaction@NSteps	
						.moveBullet(ngEnv, ngEnv$mouseHistory$x + dxy[1], ngEnv$mouseHistory$y + dxy[2])						
					}
				})
		## scroll back to 'from' node
		tkbind(tt, '<Button-5>', function(){
					if(ngEnv$bulletState$to != ''){
						xyfrom <- .tcl2xy(tkcoords(ngEnv$canvas,paste('node &&', ngEnv$bulletState$from)))
						xyto <- .tcl2xy(tkcoords(ngEnv$canvas,paste('node &&', ngEnv$bulletState$to)))
						dxy <- (xyto - xyfrom)*1/ngEnv$settings@interaction@NSteps	
						.moveBullet(ngEnv, ngEnv$mouseHistory$x - dxy[1], ngEnv$mouseHistory$y - dxy[2])
					}
				})
	}
	
	
	
	## paths: add edge - intialize path
	isKeyShift <- FALSE
	tkbind(tt,'<KeyPress-Shift_L>', function(){
				## Careful, When shifts is beeing kept pressed in windows it fire repeatedly events
				if(!ngEnv$isKeyShift) {
					if(ngEnv$bulletState$to == "") {
						tclObj(ngEnv$activePath) <- ngEnv$bulletState$from
						tclvalue(ngEnv$activePathGraph) <- ngEnv$graph@name
					}
					ngEnv$isKeyShift <- TRUE
				}
			})
	
	
	tkbind(tt,'<KeyRelease-Shift_L>', function(){
				ngEnv$isKeyShift <- FALSE
				.normalState(ngEnv)
			})
	tkbind(tt,'<FocusOut>', function(){
				ngEnv$isKeyShift <- FALSE
				.normalState(ngEnv)
			})
	tkbind(tt,'<FocusIn>', function(){
				.normalState(ngEnv)
			})
	
	
	
#	
	tkbind(tt,'<Enter>',function()tkfocus(tt))
	tkitembind(canvas,'node || label','<Shift-Button-1>', function(){})	
	tkitembind(canvas,'node || label || bullet','<Shift-Button-1>', function(){
				
			#	tkfocus(tt)	
				if(!ngEnv$isKeyShift) {
					if(ngEnv$bulletState$to == "") {
						tclObj(ngEnv$activePath) <- ngEnv$bulletState$from
						tclvalue(ngEnv$activePathGraph) <- ngEnv$graph@name
					}
					ngEnv$isKeyShift <- TRUE
				}
				
				
				if(ngEnv$bulletState$to == "") {
					
					path <- .parsePath2Vec(tclvalue(activePath))
					## get last node of activePath						
					from <- tail(path,1)
					## get node
					tags <- .tcl2str(tkgettags(ngEnv$canvas,'current'))
					if('bullet' %in% tags) {
						to <- ngEnv$bulletState$from
					} else {
						to <- tags[3]
					}
					
					#cat(paste('from',from,'to',to,'\n'))
					## is the node adjacent?
					adjEdge <- .tcl2str(.Tcl(paste(canvas$ID," find withtag {edge&&",from,'&&',to,'}', sep = '')))
					
					if(length(adjEdge) != 0) {
						if(from != to) {
							##add to active Path
							#cat('add it\n')
							tclvalue(ngEnv$activePath) <- paste(tclvalue(activePath),to)
							## Highlight Path
							.highlightPath(ngEnv,from,to)
							## Highlight adjoining edges
						}
					}
					
				}

				#else {
				#tk_messageBox(message = "The window containing the graph must be the needs Please first focus on the window containing the graph.\nAn ordinary mouse-click should be sufficient.", parent=tt)
				#}
				
			})
	
	## walk path
	tkitembind(canvas,'node || label || bullet','<Shift-Double-1>', function(){
				path <- .parsePath2Vec(tclvalue(activePath))
				## get last node of activePath						
				from <- tail(path,1)
				## get node
				tags <- .tcl2str(tkgettags(ngEnv$canvas,'current'))
				if('bullet' %in% tags) {
					to <- ngEnv$bulletState$from
				} else {
					to <- tags[3]
				}
				
				## is the node adjacent?
				adjEdge <- .tcl2str(.Tcl(paste(canvas$ID," find withtag {edge&&",from,'&&',to,'}', sep = '')))
				
				if(length(adjEdge) != 0) {
					.walkPath(ngEnv, tclvalue(ngEnv$activePath))				
				}
			})
	
	
	
	
	
	
	##
	## highlight labels and nodes
	## 
	
	## node
	tkitembind(ngEnv$canvas,'node', '<Any-Enter>', function() {
				node <- .tcl2str(tkgettags(ngEnv$canvas,'current'))[3]
				col <- .tcl2str(tkitemcget(ngEnv$canvas, paste('node && ',node), activefill = NULL))
				ngEnv$tempColor <- tkitemcget(ngEnv$canvas,paste('label &&', node), fill = NULL)
				tkitemconfigure(ngEnv$canvas,paste('label &&', node), fill = col)
			})
	tkitembind(ngEnv$canvas,'node', '<Any-Leave>', function() {
				node <- .tcl2str(tkgettags(ngEnv$canvas,'current'))[3]
				tkitemconfigure(ngEnv$canvas,paste('label &&', node), fill = ngEnv$tempColor)
			})
	## label
	tkitembind(ngEnv$canvas,'label', '<Any-Enter>', function() {
				label <- .tcl2str(tkgettags(ngEnv$canvas,'current'))[3]
				col <- .tcl2str(tkitemcget(ngEnv$canvas, paste('label && ',label), activefill = NULL))
				ngEnv$tempColor <- tkitemcget(ngEnv$canvas,paste('node &&', label),fill = NULL)
				tkitemconfigure(ngEnv$canvas,paste('node &&', label), fill = col)
			})
	tkitembind(ngEnv$canvas,'label', '<Any-Leave>', function() {
				label <- .tcl2str(tkgettags(ngEnv$canvas,'current'))[3]
				tkitemconfigure(ngEnv$canvas,paste('node &&', label), fill = ngEnv$tempColor)
			})
	
	## bullet
	tkitembind(ngEnv$canvas,'bullet', '<Any-Enter>', function() {
				if(ngEnv$bulletState$percentage == 0) {
					ngEnv$tempColor <- tkitemcget(ngEnv$canvas,paste('label &&', ngEnv$bulletState$from),fill = NULL)
					tkitemconfigure(ngEnv$canvas,paste('label &&', ngEnv$bulletState$from), fill = ngEnv$settings@color@bulletActive)					
				}
			})
	tkitembind(ngEnv$canvas,'bullet', '<Any-Leave>', function() {
				if(ngEnv$bulletState$percentage == 0) {
					tkitemconfigure(ngEnv$canvas,paste('label &&', ngEnv$bulletState$from), fill = ngEnv$tempColor)					
				}
			})
	## Edge
	tkitembind(ngEnv$canvas,'edge', '<Any-Enter>', function() {
				tags <- .tcl2str(tkgettags(ngEnv$canvas,'current'))
				node1 <- tags[3]
				node2 <- tags[4]
				ngEnv$tempColor  <- tkitemcget(ngEnv$canvas,paste('node &&', node1),fill = NULL)
				ngEnv$tempColor2 <- tkitemcget(ngEnv$canvas,paste('node &&', node2),fill = NULL)
				ngEnv$tempColor3  <- tkitemcget(ngEnv$canvas,paste('label &&', node1),fill = NULL)
				ngEnv$tempColor4 <- tkitemcget(ngEnv$canvas,paste('label &&', node2),fill = NULL)
				tkitemconfigure(ngEnv$canvas,paste('(node || label) && (', node1, ' || ', node2,')'),
						fill = ngEnv$settings@color@edgeActive)
			})
	
	tkitembind(ngEnv$canvas,'edge', '<Any-Leave>', function() {
				tags <- .tcl2str(tkgettags(ngEnv$canvas,'current'))
				node1 <- tags[3]
				node2 <- tags[4]
				tkitemconfigure(ngEnv$canvas,paste('node &&', node1), fill = ngEnv$tempColor)
				tkitemconfigure(ngEnv$canvas,paste('node &&', node2), fill = ngEnv$tempColor2)
				tkitemconfigure(ngEnv$canvas,paste('label &&', node1),fill = ngEnv$tempColor3)
				tkitemconfigure(ngEnv$canvas,paste('label &&', node2),fill = ngEnv$tempColor4)
			})
	
	
	tkconfigure(progress.label, text = "Load Graph")
	tkconfigure(progress, value = "40")
	
	
	##
	## initialize first graph in tk canvas
	##
	#####################################
	tkdelete(canvas,'all')
	.graph2canvas(ngEnv)
	.arriveAtNode(ngEnv)	
	
	
	
	
	
	
	
	
	## 
	##  Visualization of data
	##
	## #################################
	tkconfigure(progress.label, text = "Visualize Data")
	tkconfigure(progress, value = "50")
	

	
	tkconfigure(progress.label, text = "Scale Data")
	tkconfigure(progress, value = "70")
	
	
	## Tk Canvas Scatterplot
	whichScaledData <- sapply(vizList,function(x){is(x,"NG_Viztk2d")|is(x,"NG_Viz2DAxis")})
	
	if(any(whichScaledData)) {
		ngEnv$scaledData <- list()
		
		## vector with names of the data sets
		t.data <- unique(unlist(sapply(vizList,function(viz){
									if(is(viz,"NG_Viztk2d")|| is(viz,"NG_Viz2DAxis") ){
										return(viz@data)
									}else{
										return()
									}
								})))
		
		for(data in t.data) {
			## scale data
			X <- as.data.frame(lapply(dataList[[data]]@data,function(var){
								(var - mean(var))/sd(var)
							}))
			b <- sqrt(2)*max(abs(X), na.rm = TRUE) ##max(abs(min(X,na.rm = TRUE)),max(X,na.rm = TRUE))	
			Xs <- as.data.frame(lapply(X,function(var){
								var/b
							}))
			ngEnv$scaledData[[data]] <- Xs		
		}
		
		
		
		#initialize the tk scatterplot stuff
		
		
		ngEnv$t.vizcounter <- 1
		vizList <- lapply(vizList, FUN = function(viz){
					if(is(viz,"NG_Viztk2d")) {
						
						viz@viz_name <- paste("viz",t.vizcounter,sep = '') 
						vizN <- viz@viz_name
						
						ngEnv$t.vizcounter <- ngEnv$t.vizcounter + 1	
						
						
						
						
						if(!is.null(viz@glyphVarOrder)) {
							## glyph_seq
							#ng_windowManager("$ngInstance\.$viz\.glyph_seq") 
							tcl('set',paste('ng_windowManager("',ng_instance,'.',viz@viz_name,'.','glyph_seq','")',sep = ''), viz@glyphVarOrder)
							
							## glyphs coordinate
							#ng_windowManager("$ngInstance\.$viz\.glyphs")
							
							t.X <- sapply(dataList[[viz@data]]@data, FUN = function(var){(var-min(var))/diff(range(var-min(var)))})
							
							t.glyphs <- t.X[,viz@glyphVarOrder]
							
							t.array <- paste('ng_windowManager("',ng_instance,'.',viz@viz_name,'.','glyphs','")',sep = '')
							tcl('set',t.array, '')
							apply(t.glyphs,1,FUN = function(row){
										tcl('lappend',t.array, row)
									})
							
							## glyph alpha
							#ng_windowManager("$ngInstance\.$viz\.glyph_alpha")
							t.n <- length(viz@glyphVarOrder)
							tcl('set',paste('ng_windowManager("',ng_instance,'.',viz@viz_name,'.','glyph_alpha','")',sep = ''), 2*pi/t.n *(0:(t.n-1)))
						}
						
						
						
						if(!is.null(viz@imgIds)) {
							## images (ids for a particular instance)
							##set ng_windowManager("$ngInstance\.$viz\.images") {}
							##set ng_windowManager("$ngInstance\.$viz\.images_orig") {}
							
#							ng_windowManager("$ngInstance\.$viz\.images_orig")
#							ng_windowManager("$ngInstance\.$viz\.images")						
							tcl('set',paste('ng_windowManager("',ng_instance,'.',viz@viz_name,'.','images_orig','")',sep = ''), viz@imgIds)	
							
							
							t.img <- sapply(viz@imgIds, FUN = function(img){
										t.im <- as.character(tclvalue(tkimage.create('photo')))
										tcl('image_scale',img,50,t.im)
										return(t.im)
									})
							tcl('set',paste('ng_windowManager("',ng_instance,'.',viz@viz_name,'.','images','")',sep = ''), t.img)
							
							## image_w
#							ng_windowManager("$ngInstance\.$viz\.image_w2")
							t.w2 <- sapply(t.img, FUN = function(img.w){
										.tcl2num(tcl('image','width',img.w))/2
									})
							tcl('set',paste('ng_windowManager("',ng_instance,'.',viz@viz_name,'.','image_w2','")',sep = ''), t.w2)
							
							
							## image_h
#							ng_windowManager("$ngInstance\.$viz\.image_h2")
							t.h2 <- sapply(t.img, FUN = function(img){
										as.numeric(tclvalue(tcl('image','height',img)))/2	
									})
							tcl('set',paste('ng_windowManager("',ng_instance,'.',viz@viz_name,'.','image_h2','")',sep = ''), t.h2)
							
							## image_halo width (initialize 2)
							tcl('set',paste('ng_windowManager("',ng_instance,'.',viz@viz_name,'.','image_halo','")',sep = ''), rep(2,length(viz@imgIds)))
							tcl('set',paste('ng_windowManager("',ng_instance,'.',viz@viz_name,'.','image_diag_old','")',sep = ''), rep(50,length(viz@imgIds)))		
							
						}
						
						
						
						return(viz)
						
					} else {
						return(viz)
					}
					
					
				})
		
		
		defaulTkColors <- ngEnv$settings@tk2d@brush_colors
		additionalColors <- c(defaulTkColors,colors()[-match(defaulTkColors,colors())])
		ncolors <- length(additionalColors)
		
		for(dataName in t.data) {
			
			if(!.tcl2num(tcl('info','exists',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.bg")',sep = '')))) {
				cat(paste("Session ", ngEnv$ng_instance,", Data ",dataName," is new!\n",sep=''))
				
				t.ndata <- dim(ngEnv$scaledData[[dataName]])[1]
				
				## color of points (groups)
				## ng_data("$ngLinkedInstance\.$dataName\.color")
				##
				x <- which(dataNames == dataName)
				group <- dataList[[x]]@group
				if(length(group) == 0) {
					## No group variable submitted everything is red
					pointColors <- rep(defaulTkColors[1],t.ndata)
					pointSizes <- rep(8,t.ndata)
				} else {
					## coerce variable to color and size
					
					if(all(grepl("^c.+;s.+$",group))){ ## restore old one
						# restore old values
						tmp <- sapply(group, FUN = function(str) {
									unlist(strsplit(unlist(strsplit(str,"^c")),";s",fixed=TRUE))									
								})
						pointColors <- tmp[1,]
						pointSizes <- tmp[2,]
						
					} else if(all(group %in% colors())){
						pointColors <- group
						tk2dcolors[[dataName]] <- unique( c(group,defaulTkColors))[1:9]
						pointSizes <- rep(5,t.ndata)
					} else {
						## create new coloring
						pointColors <- rep(defaulTkColors[1],t.ndata)
						group.u <- unique(group)
						
						if(length(group.u)>ncolors) {
							warning(paste("Data", dataName, "has more than",ncolors,"classes. Only", ncolors, "classes will be distingushed by colors."))
							pointColors <- rep(defaulTkColors[1],t.ndata)
							for(i in 2:ncolors) {
								pointColors[group == group.u[i-1]] <- additionalColors[i]
							}	
						} else {
							i = 1
							for(group.i in group.u) {
								pointColors[group == group.i] <- additionalColors[i]	
								i = i + 1
							}
						}
						pointSizes <- rep(5,t.ndata)
						
					}
					
				}
				
				tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','color','")',sep = ''), pointColors)
				tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','size','")',sep = ''), pointSizes)
				if(!is.null(tk2dcolors[[dataName]])){
					tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','brush_colors','")',sep = ''), tk2dcolors[[dataName]])
				} else {
					tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','brush_colors','")',sep = ''), defaulTkColors)
				}
				
				
				if(is.null(tk2dcolors[[paste('sel:',dataName,sep='')]])){
					tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','brush_color','")',sep = ''),ngEnv$settings@tk2d@brush_color)
				} else {
					tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','brush_color','")',sep = ''), tk2dcolors[[paste('sel:',dataName,sep='')]] )
				}
				
				if(is.null(tk2dcolors[[paste('bg:',dataName,sep='')]])){
					tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','bg','")',sep = ''), ngEnv$settings@tk2d@bg)
				} else {
					tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','bg','")',sep = ''), tk2dcolors[[paste('bg:',dataName,sep='')]] )
				}
				
				
				## text (labels)
				## ng_data("$ngLinkedInstance\.$dataName\.text")
				if(length(dataList[[dataName]]@labels)){
					tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','text','")',sep = ''), dataList[[dataName]]@labels)
				}else{
					tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','text','")',sep = ''), '')
				}	
				
				## selected points (initialize none)
				##ng_data("$ngLinkedInstance\.$dataName\.selected")
				tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','selected','")',sep = ''), rep(0,t.ndata))
				tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','deactivated','")',sep = ''), rep(0,t.ndata))
				tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','anyDeactivated','")',sep = ''), 0)
				
				## total_brushed (initialize none)
				## ng_data("$ngLinkedInstance\.$dataName\.total_brushed") ""
				tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','total_brushed','")',sep = ''), "")
				
				## radius of dots (for the moment initialize with 8)
				## ng_data("$ngLinkedInstance\.$dataName\.radius") 
				tcl('set',paste('ng_data("',ng_LinkedInstance,'.',dataName,'.','size','")',sep = ''), rep(5,t.ndata))
				
				
			} else {
				cat(paste("Session ", ngEnv$ng_instance,", Data ",dataName," is linked!\n",sep=''))
			}
		}
	}
	
	tkconfigure(progress.label, text = "Vizlist")
	tkconfigure(progress, value = "90")
	
	
	## Make new lists: viz1, viz2, ... vizn (n = number of graphs)
	## each list contains the visualizations for a graph
	for(i in 1:nGraphs){
		t.graph <- get(paste('graph',i, sep = ''))@name
		t.vizList <- sapply(vizList,
				function(viz){
					if(viz@graph == t.graph){
						return(viz)
					}else{
					}
				})
		
		## Delete NULL entries
		t.viz <- t.vizList[unlist(lapply(t.vizList, length) != 0)]
		
		assign(paste('viz',i,sep = ''),t.viz, envir = ngEnv)
	}
	selViz <- 'viz1'
	
	
	
	
	.initializePlots(ngEnv)
	
	
# Create navGraph handler
	if(!is(data, "NavGraph_handler")){
		## Create and return the navGraph Handler
		ng <- new("NavGraph_handler",
				env = ngEnv,
				graphs = graphList,
				data = dataList,
				viz = vizList,
				settings = settings,
				paths = paths,
				activePath = "",
				activePathGraph = "",
				dateCreated = date(),
				dateUpdated = "not")
		
		tkconfigure(progress.label, text = "Vizlist")
		tkconfigure(progress, value = "100")
		tkdestroy(ttp)
		
	}else{
		ng <- data
		ng@env <- ngEnv
		ng@dateUpdated <- date()
	}
	
	isDestoying <- FALSE
	return(ng)
}

