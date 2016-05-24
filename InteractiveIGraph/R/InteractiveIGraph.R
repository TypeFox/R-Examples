
InteractiveIGraph <- function(g, ...){
	
	### susikuriam ir nusipiesiam
	GName = deparse(substitute(g))

	SavedGraph <- InteractiveIGraph.Constructor(g, GName=GName, ...)
	g <- plot.InteractiveIGraph(SavedGraph)
	LastGraph <- g


	### apsibreziam kintamuosius
	ButtonsDown <- NULL	
	StartX <- 0
	StartY <- 0	
	MD <- NULL	
	MenuSelected <- FALSE
	LogMsg <- ""
	
	
    devset <- function()
        if (dev.cur() != eventEnv$which) dev.set(eventEnv$which)

		
	mousedown <- function(buttons, x, y) {
		### issisaugom pries pat pakeitimus
		LastGraph <<- g		
	 
		### issivalom
		g <<- GraphClear(mode="mousedown", g=g)
	
	
		# devset()
		x = grconvertX(x, "ndc", "user")
		y = grconvertY(y, "ndc", "user")
		
		LogMsg <<- ""
		Msg = paste("ButtonsDown ", paste(buttons, collapse=" "), "(at ", round(x,3), round(y,3), ") ")
		LogMsg <<- paste(LogMsg, Msg, sep = "")
		
		
		### nusirodom kintamuosius
		ButtonsDown <<- buttons		
		StartX <<- x
		StartY <<- y			
		MD <<- NULL	
		MenuSelected <<- FALSE
		
		### Pasiziurim ar buvo nuspaustas koks nors Menu, jei taip tai ji ivykdome (aisku jei aktyvus)
		MenuItem = GetMenuItem(x, y, g$Menu$ActiveMenuList)
		if(!is.null(MenuItem)){
			if(!is.null(MenuItem$FUN)){
				res = do.call(MenuItem$FUN, args=as.list(MenuItem$FUNParams))
			}	
			#gr123 = g
			if(!is.null(MenuItem$command)){
				# res = try(eval(parse(text=MenuItem$command), envir = sys.parent(n = 3)),silent=TRUE)
				res = try(eval(parse(text=MenuItem$command), envir = as.environment(-1)),silent=TRUE)
			}	
			MenuSelected <<- TRUE
			Msg = " -> Menu("%.% MenuItem$label %.% ")"
			LogMsg <<- paste(LogMsg, Msg, sep = "")
			return(NULL)
		}
	
	
		### pasizuirim ar pazymejo koki nors verteksta
		MD <<- MinDistGraph(StartX, StartY, g)	
			
		if(!(MD$selected & V(g)[MD$vid]$selected) & buttons[1]!=1){
			eventEnv$onMouseMove <- dragmousemove		
		}
		
			
		NULL
	}
	
	dragmousemove <- function(buttons, x, y) {
		x = grconvertX(x, "ndc", "user")
		y = grconvertY(y, "ndc", "user")

		polygon(c(StartX,x,x,StartX), c(StartY,StartY,y,y), col=rgb(0,100,0,10,maxColorValue=255), border = NA)
		NULL
	}	
	
	mouseup <- function(buttons, x, y) {
		eventEnv$onMouseMove <- NULL		
		
		### issivalom
		g <<- GraphClear(mode="mouseup", g=g)
		
		x = grconvertX(x, "ndc", "user")
		y = grconvertY(y, "ndc", "user")		
		
		if(length(ButtonsDown)==1 & !MenuSelected){
			g <<- SelectEdge(eids=NULL, vids=NULL, g=g)	
			if(ButtonsDown==0){	
				DragDist = sqrt((x-StartX)^2 + (y-StartY)^2)
				if(DragDist<3/200){ #jei nebuvo peles judesio, tai zymim arba atzymim konkrecia virsune  
					if(MD$selected){
					    g <<- SelectVertex(vids=MD$vid, g=g)
						Msg = " -> " %.% GetObjectLabel(ids=MD$vid, type="V", g=g)
						LogMsg <<- paste(LogMsg, Msg, sep = "")					   
					}else{
					    g <<- SelectVertex(vids=NULL, g=g)
					}				
				}else{ # jei buvo judejimas 
					if(MD$selected & V(g)[MD$vid]$selected ){ # jei pagavom jau pazymeta, tai judiname visus pazymetus
						V(g)[selected]$x <<- V(g)[selected]$x - StartX + x
						V(g)[selected]$y <<- V(g)[selected]$y - StartY + y
						g <<- XY.norm(g)
						LogMsg <<- paste(LogMsg, "-> Move selected" , sep = "")	
					} else { # zymime regijona
						XX = range(c(StartX, x))
						YY = range(c(StartY, y))
						XSelect = XX[1] <= V(g)$x & V(g)$x <= XX[2]
						YSelect = YY[1] <= V(g)$y & V(g)$y <= YY[2]
						vids = V(g)[active & !hidden & XSelect & YSelect]
						g <<- SelectVertex(vids=vids, g=g)
						LogMsg <<- paste(LogMsg, "-> select region" , sep = "")	
					}
				}					
			}
			if(ButtonsDown==2){
				g$Menu$MenuList$MainMenu$active <<- TRUE		
				g$Menu$MenuList$MainMenu$Params$x <<- x
				g$Menu$MenuList$MainMenu$Params$y <<- y
				LogMsg <<- paste(LogMsg, "-> Menu" , sep = "")	
			} 
			if(ButtonsDown==1){
				g$PlotParam$xlim <<- g$PlotParam$xlim + StartX-x
				g$PlotParam$ylim <<- g$PlotParam$ylim + StartY-y
				LogMsg <<- paste(LogMsg, "-> Move plot" , sep = "")	
			}		
		}

		
		g <<- MsgToLogObj(LogMsg, g=g)
		g <<- plot.InteractiveIGraph(g)	
		
		
		NULL
	}

	keybd <- function(key) {
		
	    eventEnv$onMouseMove <- NULL
	   	
		# issivalom
		g <<- GraphClear(mode="keybd", g=g)
		LogMsg <<- ""

		
		### escape - cia reikia realizuoti visapusiska issivalima
		if(key=="ctrl-["){
			g <<- GraphClear(mode="esc", g=g)
			g <<- MsgToLogObj(Msg="<Esc>", g=g, add=FALSE)
			g <<- plot.InteractiveIGraph(g)
			return(NULL)
		}	
		
		### enter
	    if(key=="ctrl-J"){
		
			if(g$mode$input){
				if(g$mode$InputType=="Command"){				
					res = try(eval(parse(text = g$ExtraParam$input)), silent=TRUE)					
				}
				if(g$mode$InputType=="Attributes"){				
					g <<- AttributeChange(AttrCommand=g$ExtraParam$input, g=g)
				}
				cat("\n")
				g$ExtraParam$input <<- ""
				g$mode$input <<- FALSE
				g$Menu$MenuList$BottomMenu$MenuLines$Input$active <<- FALSE				
			}else{
				AO = GetActiveObject(g=g)
				g <<- PutViewObject(oProgId=AO$ProgId, type=AO$type, g=g)
			}			
			g <<- plot.InteractiveIGraph(g)
			return(NULL)
		}	  

		### jei input reikalaujantys klavisai, tai perkraunam input parametrus
		if(key %in% c("ctrl-C", "ctrl-F")){
			ClearBottomMenuPlot(g)
			g$mode$input <<- FALSE
			g$mode$InputType <<- "none"
			g$ExtraParam$input <<- ""
		}
	
		
		### toliau pildom inputa
		if(g$mode$input){
		    ClearBottomMenuPlot(g)
			if(key=="ctrl-V")		
				key = readClipboard()				
			if(key=="ctrl-H"){	# backspace	 
				g$ExtraParam$input <<- substr(g$ExtraParam$input, start=1, stop=nchar(g$ExtraParam$input)-1)
				cat("\nCommand: "%.%g$ExtraParam$input)
			}else{
				g$ExtraParam$input <<- g$ExtraParam$input %.% key
				cat(key)
			}
			g$Menu$MenuList$BottomMenu$MenuLines$Input$active <<- TRUE
			g$Menu$MenuList$BottomMenu$MenuLines$Input$MenuItems[[1]]$label <<- g$mode$InputType %.%": " %.% g$ExtraParam$input
			g$Menu$ActiveMenuList <<- list(BottomMenu=plot.Menu(g$Menu$MenuList$BottomMenu, x = par("usr")[1], y = par("usr")[3]))
			return(NULL)
		}

		
		### gyzimas atgal
		if(key == "ctrl-Z"){
			g <<- LastGraph
			g <<- plot.InteractiveIGraph(g)
			return(NULL)
		}
		### jei jau atejom iki cia, tai issisaukom	
		LastGraph <<- g	
	
		
		
	
		### vykdom normalias komandas
		LogMsg <<- paste(LogMsg, "Key <", key, ">", sep = "")
		if(!is.null(CommandList[[key]])){
			LogMsg <<- paste(LogMsg, "(", CommandList[[key]]$label, ")", sep = "")
			g <<- MsgToLogObj(LogMsg, g=g, add=FALSE)
			g <<- CommandList[[key]]$FUN(g)	
		} else {
			LogMsg <<- paste(LogMsg, "(nothing)", sep = "")
			g <<- MsgToLogObj(LogMsg, g=g, add=FALSE)
		}
	
		
		
		### isejimas
		if(g$ExtraParam$quite){
			ClearBottomMenuPlot(g)
			g$ExtraParam$quite <<- FALSE
			gr <- MsgToLogObj("bye.", g=g, add=FALSE)
			gr <- plot.Menu(gr$Menu$MenuList$BottomMenu, x = par("usr")[1], y = par("usr")[3])
			return(invisible(g))
		}
		### issaugojimas
		if(g$ExtraParam$save){
			g$ExtraParam$save = FALSE
			SavedGraph <<- g	
		}
		
		if(g$ExtraParam$replot){
			g <<- plot.InteractiveIGraph(g)		
			
		} else {
		    ClearBottomMenuPlot(g)
			g$Menu$ActiveMenuList <<- list(BottomMenu=plot.Menu(g$Menu$MenuList$BottomMenu, x = par("usr")[1], y = par("usr")[3]))
		}

		NULL
	}

 	
	setGraphicsEventHandlers(prompt="Pleasy, enjoy.",
				 onMouseDown = mousedown,
				 onMouseUp = mouseup,
				 onMouseMove = NULL,
				 onKeybd = keybd)
    eventEnv <- getGraphicsEventEnv()
	getGraphicsEvent()	
}

plot.InteractiveIGraph <- function(x, ...){
	g=x
		
	### pasirupiname, kad objektas tikrai butu reikiamos kalses, o jei jau toks yra, kad grafiniai parametrai butu overridinti	
	g <- InteractiveIGraph.Constructor(g, ...) 
	
	### apskritai pasaliname pasleptus kintamuosius	
	gplot = induced.subgraph(g, vids=V(g)[!hidden])
	gplot = delete.edges(gplot, edges=E(gplot)[hidden])

	### apsibreziame dirbtinus kintamuosius, kad jie islikytu teisingas kordinates. Jas pridesime tik pries plotinima
	DummyVertices = vertices(name=c("dummy1","dummy2"), label=c(NA, NA), x=range(V(g)$x), y=range(V(g)$y)
		, active=FALSE, hidden=TRUE, ProgType = "dummy", shape = "circle", size=1)
		
		
	# gplotList <- g$functions$PlotAdjustment(gplot)
	gplotList <- do.call(g$functions$PlotAdjustment, args=list(gplot))
	
	if(length(gplotList)==0){
		warning("PlotAdjustment return object of length 0.")
		plot(0,0, type="n")
	}	

	
	### breziam (kolkas juodrastinis)
	OptParSave = par(mar = gplot$ExtraParam$mar)	
	
	### pirmasis lygis
	gi = gplotList[[1]] 
	PlotParam = gi$PlotParam
	PlotParam$x = gi + DummyVertices
	PlotParam=PlotParam[-(which(sapply(PlotParam,is.null),arr.ind=TRUE))]
	do.call("plot.igraph", args=PlotParam)
	

	### visi like lygiai
	if(length(gplotList)>1){
		for(gi in gplotList[-1]){
			PlotParam = gi$PlotParam 
			PlotParam$x = gi + DummyVertices
			PlotParam$add = TRUE
			PlotParam=PlotParam[-(which(sapply(PlotParam,is.null),arr.ind=TRUE))]
			do.call("plot.igraph", args=PlotParam)		
		}	
	}
	
	### Menu
	gplot <- do.call(g$functions$BottomMenuAdjustment, args=list(gplot))
	
	gplot$Menu$MenuList$BottomMenu$Params$x <- par("usr")[1]
	gplot$Menu$MenuList$BottomMenu$Params$y <- par("usr")[3]		
	g$Menu$ActiveMenuList <- lapply(gplot$Menu$MenuList, plot.Menu)		
	
	
	par(OptParSave)
	
	g$ExtraParam$replot <- FALSE	
	
	invisible(g)
}

GraphClear <- function(mode=c("esc", "mousedown", "mouseup", "keybd"), g){
	
	mode = match.arg(mode)

	if(mode=="esc"){
		g$mode$input = FALSE
		g$mode$InputType = "none"		
		g$ExtraParam$input <- ""
		g$mode$select = FALSE
		g$mode$AllEdges = FALSE
		g$Menu$MenuList$BottomMenu$MenuLines$Input$active <- FALSE
		g$Menu$MenuList$MainMenu$active <- FALSE
		g$Menu$MenuList$BottomMenu$active <- TRUE
		V(g)$selected = FALSE
		E(g)$selected = FALSE
		g <- PutActiveObject(id=NA, type="none", g=g)
		g <- PutViewObject(id=NA, type="none", g=g)	
	}	
	if(mode=="mousedown"){
		g$mode$input <- FALSE
		g$ExtraParam$input <- ""
		g$Menu$MenuList$BottomMenu$MenuLines$Input$active <- FALSE
	}
	if(mode=="mouseup"){
		g$Menu$MenuList$MainMenu$active <- FALSE	 ### pagrindini Menu panaikiname	
	}
	
	if(mode=="keybd"){
		g$Menu$MenuList$MainMenu$active <- FALSE	 ### pagrindini Menu panaikiname
	}
	
	g
}	

InteractiveIGraph.Constructor <- function(g, ...){

	dots = list(...)

	if(class(g)[1]=="InteractiveIGraph"){
		
		if(length(dots)==0) return(g)
		
		### cia tik overridinam kaikurias reikmes
		
		### PlotParam
		NameList = c("xlim", "ylim", "main", "sub", "axes", "xlab", "ylab" , "add") 
		TakeNames = intersect(NameList, names(dots))
		g$PlotParam[TakeNames] = dots[TakeNames]
		dots[TakeNames] = NULL

		### ExtraParam
		if(!is.null(dots[["mar"]])){
			g$ExtraParam$mar = dots[["mar"]]
			dots[["mar"]] = NULL
		}

		g$ExtraParam$GName = "Interactive Graph"
		if(!is.null(dots[["GName"]])) g$ExtraParam$GName = dots[["GName"]]
		if(is.null(g$PlotParam$main)) g$PlotParam$main = g$ExtraParam$GName
		dots[["GName"]] = NULL
		
		### jei jau cia pateko uztikriname svarbus parametrai butu geri
		g$ExtraParam$replot = FALSE	
		# g$ExtraParam$save = FALSE	
		
	} else {
	
		### layout: 
		if(!is.null(dots[["layout"]])) g$layout = dots[["layout"]]
		if(!is.null(g$layout)){ 
			mat = layout.norm(g$layout, -1, 1, -1, 1)
		} else {
			if(!is.null(V(g)$x) & !is.null(V(g)$y)){
				if(any(c(is.na(V(g)$x),c(is.na(V(g)$y))))){
				    g = remove.vertex.attribute(g, "x")
				    g = remove.vertex.attribute(g, "y")
					g$layout = NULL
					g$layout = layout.auto(g)
					mat = layout.norm(g$layout, -1, 1, -1, 1)
				} else {
					mat = layout.norm(cbind(V(g)$x, V(g)$y), -1, 1, -1, 1)
				}
			} else {
				g$layout = layout.auto(g)
				mat = layout.norm(g$layout, -1, 1, -1, 1)
			}		
		}
		V(g)$x = mat[,1]
		V(g)$y = mat[,2]
		dots[["layout"]] =  NULL
		g = remove.graph.attribute(g, "layout")

	
		### Grafiniai parametrai
		NameList = c("xlim", "ylim", "main", "sub", "axes", "xlab", "ylab" , "add") 
		PlotParam = vector("list", length(NameList))
		names(PlotParam) = NameList
		g = set.graph.attribute(g, "PlotParam", PlotParam)
		# butinuju parametru defoltai
		g$PlotParam$xlim = c(-1, 1)
		g$PlotParam$xlim = c(-1, 1)
		g$PlotParam$ylim = c(-1.1, 1)
		g$PlotParam$add = FALSE		
		# overridinimas
		TakeNames = intersect(NameList, names(dots))
		g$PlotParam[TakeNames] = dots[TakeNames]
		dots[TakeNames] = NULL
	
		### reikalingi atributai		
		# butinieji ir programiniai
		g = TestAndSet.vertex.params(graph=g, name="active", default=TRUE)
		g = TestAndSet.vertex.params(graph=g, name="hidden", default=FALSE)
		g = TestAndSet.vertex.params(graph=g, name="selected", default=FALSE)
		g = TestAndSet.vertex.params(graph=g, name="block", default="")
		g = TestAndSet.vertex.params(graph=g, name="ProgType", default="normal")
		g = TestAndSet.vertex.params(graph=g, name="ProgId", default=seq_along(V(g)))
		g = TestAndSet.vertex.params(graph=g, name="name", default=as.character(V(g)$ProgId))
		g = TestAndSet.vertex.params(graph=g, name="label", default=V(g)$name)
		
		g = TestAndSet.edge.params(graph=g, name="ProgId", default=seq_along(E(g)))
		g = TestAndSet.edge.params(graph=g, name="name", default=as.character(E(g)$ProgId))
		g = TestAndSet.edge.params(graph=g, name="active", default=TRUE)
		g = TestAndSet.edge.params(graph=g, name="hidden", default=FALSE)
		g = TestAndSet.edge.params(graph=g, name="selected", default=FALSE)
		g = TestAndSet.edge.params(graph=g, name="ProgType", default="normal")
	

		# graifniai butinieji atributai / defoltai, siap tai reiktu padaryti, kad is dots irgi sugebetu pasikrauti. Kolkas praleidziam, tegu koreguoja per atributus.
		g = TestAndSet.vertex.params(graph=g, name="frame.color", default=NA)
		g = TestAndSet.vertex.params(graph=g, name="label.color", default="black")
		g = TestAndSet.vertex.params(graph=g, name="color", default="SkyBlue2")
		g = TestAndSet.vertex.params(graph=g, name="shape", default="circle")
		g = TestAndSet.vertex.params(graph=g, name="size", default=15)
		
		g = TestAndSet.edge.params(graph=g, name="lty", default=1)
		g = TestAndSet.edge.params(graph=g, name="width", default=1)
		g = TestAndSet.edge.params(graph=g, name="color", default="darkgrey")
		g = TestAndSet.edge.params(graph=g, name="arrow.size", default=0.7) ### nebutinas, bet norimas

		
		### parametai, ir saugomi kintamieji
		NameList = c("input", "replot", "quite", "save", "mar", "generated", "type" , "GName") 
		ExtraParam = vector("list", length(NameList))
		names(ExtraParam) = NameList
		g = set.graph.attribute(g, "ExtraParam", ExtraParam)
		
		g$ExtraParam$input = ""
		g$ExtraParam$replot = FALSE			
		g$ExtraParam$quite = FALSE
		g$ExtraParam$save = FALSE
		g$ExtraParam$mar = c(0,0,1,0)
		g$ExtraParam$generated = "igraph"	# is ko sudeneruota
		g$ExtraParam$type = "none"	# tipas - kad ir kam kokia prasme turetu
		
		if(!is.null(dots[["mar"]])) g$ExtraParam$mar = dots[["mar"]]
		dots[["mar"]] = NULL

		g$ExtraParam$GName = "Interactive Graph"
		if(!is.null(dots[["GName"]])) g$ExtraParam$GName = dots[["GName"]]
		if(is.null(g$PlotParam$main)) g$PlotParam$main = g$ExtraParam$GName
		dots[["GName"]] = NULL
		
		
		### grupes
		g$groups = list()
		if(!is.null(dots$groups)){
			for(GroupsInfoList in dots$groups){
				g <- GroupCreate(g=g, GroupsInfoList=GroupsInfoList)
			}		
		} else {
			if(!is.null(dots$mark.groups)){
				for(vProgIds in dots$mark.groups){
					g <- GroupCreate(g=g, vProgIds=vProgIds)
					
				}
			}
		}
		dots[c("groups","mark.groups")] = NULL

		### busenu nustatymai 
		NameList = c("select", "AllEdges", "ActiveObjectType", "ActiveObjectProgId", 
			"ViewObjectType", "ViewObjectProgId", "input" , "InputType") 
		mode = vector("list", length(NameList))
		names(mode) = NameList
		g = set.graph.attribute(g, "mode", mode)
		
		
		g$mode$select = FALSE
		g$mode$select = FALSE
		g$mode$AllEdges = FALSE
				
		g$mode$ActiveObjectType = "none" # "none", "V", "E", "G", "B"
		g$mode$ActiveObjectProgId = 0
		
		g$mode$ViewObjectType = "none" # "none", "V", "E", "G", "B"
		g$mode$ViewObjectProgId = 0

		g$mode$input = FALSE
		g$mode$InputType = "none" # "none", "Command", "Attributes"
		
		

		### apdorojimo funckijos:
		NameList = c("PlotAdjustment", "BottomMenuAdjustment", "ExtraInfo") 
		functions = vector("list", length(NameList))
		names(functions) = NameList
		g = set.graph.attribute(g, "functions", functions)

		g$functions$PlotAdjustment = "PlotAdjustment.default" # nezinau kodel, bet reikia dvieju eiluciu - tuomet gerai
		g$functions$PlotAdjustment = "PlotAdjustment.default"
		g$functions$BottomMenuAdjustment = "BottomMenuAdjustment.default"
		g$functions$ExtraInfo = "ExtraInfo.default"	
		# Overridinimas	
		TakeNames = intersect(NameList, names(dots))
		g$functions[TakeNames] = dots[TakeNames]
		dots[TakeNames] = NULL

		### Menu 
		
		### Apatinis Menu, cia nekoreguojasm, bet per funkcijas galima nurodyti ka rodyti. Cia galima tik ji isjungti
		InputItem = list(label="Input: ", RegionParams=list(YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8))
		InputLine = list(active = FALSE, MenuItems=list(InputItem),RecParams=list(lwd = NA, border=NA))
		InfoItem = list(label="Info:", RegionParams=list(YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8))
		InfoLine = list(active = FALSE, MenuItems=list(InfoItem),RecParams=list(lwd = NA, border=NA))
		AItem = list(label="A: ", RegionParams=list(YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8))
		ALine = list(active = FALSE, MenuItems=list(AItem),RecParams=list(lwd = NA, border=NA))
		GItem = list(label="G: ", RegionParams=list(YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8))
		GLine = list(active = FALSE, MenuItems=list(GItem),RecParams=list(lwd = NA, border=NA))

		LogItem = list(label="Log: ", RegionParams=list(YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8))
		LogLine = list(active = TRUE, MenuItems=list(LogItem),RecParams=list(lwd = NA, border=NA))
		BottomMenu = list(active = TRUE, MenuLines = list(Input=InputLine, Info=InfoLine, A=ALine, G=GLine, Log=LogLine)
			,Params=list(pos = 3), RecParams=list(col=rgb(230,230,230,255,maxColorValue=255), lwd = NA, border=NA))
		
		if(!is.null(dots[["BottomMenu.active"]])) BottomMenu$active = dots[["BottomMenu.active"]]
		dots[["BottomMenu.active"]] = NULL
		
		if(is.null(dots[["MainMenu"]])){	
			### pagrindinis Menu, cia sutvarkomas defoltas, jei useris nori kito metiu tegul juo ir pasirupina
			MenuItem = list(label="Print all commads", command = "keybd(key='m')", RecParams=list(col=rgb(220,220,220,255,maxColorValue=255), lwd = 1), TextParams=list(col=1))
			MenuLine1 = list(MenuItems=list(MenuItem), RecParams=list(col=rgb(220,220,220,255,maxColorValue=255), border=1))
			MenuItem = list(label="Groupe", command = "keybd(key='g')" , RecParams=list(col=rgb(220,220,220,255,maxColorValue=255), lwd = 1), TextParams=list(col=1))
			MenuLine2 = list(MenuItems=list(MenuItem), RecParams=list(col=rgb(220,220,220,255,maxColorValue=255), border=1))
			MainMenu <- list(active = FALSE, MenuLines = list(MenuLine1,MenuLine2),RecParams=list(col=rgb(220,220,220,255,maxColorValue=255), lwd = 1, border=1))
		} else {
			MainMenu = dots[["MainMenu"]]
		}
		# g$Menu$MenuList = list(BottomMenu=BottomMenu, MainMenu=MainMenu)	
		g = set.graph.attribute(g, "Menu", list(MenuList=list(), ActiveMenuList=list()))
		g$Menu$MenuList = list(BottomMenu=BottomMenu, MainMenu=MainMenu)
		g$Menu$MenuList = list(BottomMenu=BottomMenu, MainMenu=MainMenu)

		if(class(g)[1]!="InteractiveIGraph") class(g) <- c("InteractiveIGraph", class(g))	
	}
	
	### galiausiai
	if(length(dots)>0){
		warning("Some of the paramyters was not used. Namely: " %.% names(dots))
	}
	
	return(g)
}	

TestAndSet.vertex.params <- function(graph, name, default=NA){
	value = get.vertex.attribute(graph=graph, name)
	if(is.null(value)){
		g = set.vertex.attribute(graph=graph, name, value=default)
	} else {
		value[is.na(value)] = default
		g = set.vertex.attribute(graph=graph, name, value=value)
	}
	g
}

TestAndSet.edge.params <- function(graph, name, default=NA){
	value = get.edge.attribute(graph=graph, name)
	if(is.null(value)){
		g = set.edge.attribute(graph=graph, name, value=default)
	} else {
		value[is.na(value)] = default
		g = set.edge.attribute(graph=graph, name, value=value)
	}
	g
}

as.igraph.InteractiveIGraph <- function(x, KeepAttr = NULL, ...){

	x = induced.subgraph(x, vids=V(x)[ProgType=="normal"])
	x = subgraph.edges(x, eids=E(x)[ProgType=="normal"], delete.vertices = FALSE)
	
	GDeleteAttr = c("PlotParam", "ExtraParam", "groups", "mode", "functions", "Menu")
	GDeleteAttr = setdiff(GDeleteAttr, KeepAttr)
	GDeleteAttr = intersect(GDeleteAttr, list.graph.attributes(x))
	for(atr in GDeleteAttr){
		x[[atr]] = NULL
	}	
	
	
	VDeleteAttr = c("ProgId", "active", "hidden", "selected", "ProgType", "block")
	VDeleteAttr = setdiff(VDeleteAttr, KeepAttr)
	VDeleteAttr = intersect(VDeleteAttr, list.vertex.attributes(x))
	for(atr in VDeleteAttr){
		x = remove.vertex.attribute(graph=x, name=atr)
	}	
	
	EDeleteAttr = c("ProgId", "active", "hidden", "selected", "ProgType", "block")
	EDeleteAttr = setdiff(EDeleteAttr, KeepAttr)
	EDeleteAttr = intersect(EDeleteAttr, list.edge.attributes(x))
	for(atr in EDeleteAttr){
		x = remove.edge.attribute(graph=x, name=atr)
	}	

	x
}



#####################################################################################################
#####################################################################################################
#####################################################################################################

##### spec functions #####



PrintCommandList <- function(){	
	for(i in seq_along(CommandList)){
		cat(paste(format(names(CommandList[i]), width=6) , " : ",  CommandList[[i]]$label, " -- ", CommandList[[i]]$description,  ".\n", sep=""))
	}	
}

GetShortestPath <- function(mode="all", g){

	if(mode=="all") directed=FALSE else directed=TRUE

	eProgIds = NULL

	gi = induced.subgraph(g, vids=V(g)[active & !hidden])
	gi = delete.edges(gi, edges=E(gi)[!active])
	AO = GetActiveObject(gi)
	if(AO$type=="V"){
		AOvid = as.numeric(V(gi)[ProgId==AO$ProgId])
		vids = setdiff(as.numeric(V(gi)[selected & ProgType=="normal"]),AOvid)
		if(length(vids)>0){
		
			VPathsList = list()
			msg = tryCatch({
				VPathsList = get.shortest.paths(gi, from=AOvid, to=vids, output="vpath", mode=mode)
				msg = ""	
			}, error = function(msg) msg, warning=function(msg) msg)			
						

			if(length(VPathsList)>0){		
				EPathsList <- vector("list", length(VPathsList))
				for(i in seq_along(VPathsList)){
					if(mode=="in") VPathsList[[i]]=rev(VPathsList[[i]])
					EPathsList[[i]] <- E(gi, path=VPathsList[[i]], directed=directed)
				}
				EPaths = unique(unlist(EPathsList))		
				eProgIds = E(gi)[EPaths]$ProgId	
			} else {
				eProgIds = vector()
			}
		}
	}
	list(eProgIds=eProgIds, msg=msg)
}

MinDistGraph <- function(x, y, g=g){ 	
	vec = (V(g)[!hidden]$x-x)^2 + (V(g)[!hidden]$y-y)^2
	minID = which.min(vec)
	vid = as.numeric(V(g)[!hidden])[minID]
	vProgId = V(g)[vid]$ProgId
	Dist = sqrt(vec[minID])
	return(list(vProgId=vProgId, vid = vid, value=Dist, selected = Dist<=V(g)[vid]$size/200))		
}

PlotRegionGeneral <- function(x, y, VRange=1, col=rgb(0,100,0,50,maxColorValue=255), border = NA){
	grad = seq(from = 0, to = 2*pi, length.out = 100)
	x = cos(grad)*VRange + x
	y = sin(grad)*VRange + y	
	polygon(x=x, y=y, col=col, border = border)
}

PlotActiveRegionsGraphs <- function(g=g, col=rgb(0,100,0,50,maxColorValue=255), border = NA){
	for(i in as.numeric(V(g)[!hidden])){
		PlotRegionGeneral(V(g)[i]$x, V(g)[i]$y, VRange=V(g)[i]$size/200, col=col, border=border)
	}
} 

XY.norm <- function(g=g){
	V(g)$x[is.na(V(g)$x)] <- 0
	V(g)$y[is.na(V(g)$y)] <- 0
	mat = layout.norm(cbind(V(g)$x,V(g)$y), -1, 1, -1, 1)
	V(g)$x = mat[,1]
	V(g)$y = mat[,2]
	g
}

MsgToLogObj <- function(Msg="", g=g, add=FALSE){
	if(add){
		g$Menu$MenuList$BottomMenu$MenuLines[['Log']]$MenuItems[[1]]$label <- g$Menu$MenuList$BottomMenu$MenuLines[['Log']]$MenuItems[[1]]$label %.% Msg
	} else {
		g$Menu$MenuList$BottomMenu$MenuLines[['Log']]$MenuItems[[1]]$label <- Msg
	}
	cat(g$Menu$MenuList$BottomMenu$MenuLines[['Log']]$MenuItems[[1]]$label)
	cat("\n")
	
	g
}


### Objects manipulations ###

SelectVertex <- function(vids=NULL, vProgIds=NULL, g=g){	

	### pasirupiname teisingais indeksais
	if(length(vProgIds)>0){
		vids = as.numeric(V(g)[ProgId %in% vProgIds])
	}


	if(g$mode$select){
		if(length(vids)>0){
			V(g)[vids]$selected = !V(g)[vids]$selected
		}
	}else{
		if(length(vids)>0){
			V(g)$selected = FALSE
			V(g)[vids]$selected = TRUE
		} else {
			V(g)$selected = FALSE
		}
	}
	
	### aktivuoti objekta
	if(length(vids)>0){
		if(all(V(g)[vids]$ProgType == "groupe")){ # taisytina(Test)
			if(all(V(g)[vids]$selected) & length(vids)==1){
				# bid = which(sapply(g$groups, function(G) G$ProgId==V(g)[vids]$ProgId))[1] # taisytina		
				# browser()
				g = PutActiveObject(oProgId=V(g)[vids]$ProgId, type="G", g=g)
			} else {
				g = PutActiveObject(oProgId=NA, type="G", g=g)
			}
		} else {
			if(all(V(g)[vids]$selected) & length(vids)==1) g = PutActiveObject(oProgId=V(g)[vids]$ProgId, type="V", g=g) # taisytina	(Test)
			if(!all(V(g)[vids]$selected)) g = PutActiveObject(oProgId=NA, type="V", g=g)
			if(length(vids)>1) 	g = PutActiveObject(oProgId=NA, type="V", g=g)
		}
	} else {
		if(!g$mode$select) 	g = PutActiveObject(oProgId=NA, type="V", g=g)
	}

	g
}

SelectEdge  <- function(eids=NULL, eProgIds=NULL, vids=NULL, vProgIds=NULL, g=g){	

	### pasirupiname teisingais indeksais
	if(length(eProgIds)>0){
		eids = as.numeric(E(g)[ProgId %in% eProgIds])
	}
	if(length(vProgIds)>0){
		vids = as.numeric(V(g)[ProgId %in% vProgIds])
	}


	if(g$mode$select){
		if(length(eids)>0){
			E(g)[eids]$selected = !E(g)[eids]$selected
		}
	} else {
		if(length(eids)>0){
			E(g)$selected = FALSE
			E(g)[eids]$selected = TRUE
		} else {
			E(g)$selected = FALSE
		}
	}
		
	
	if(length(vids)>0){
		E(g)[adj(V(g)[vids])]$selected = TRUE	
		if(!g$mode$AllEdges){
			E(g)[adj(V(g)[!active])]$selected = FALSE	
		}
	}
	
	### aktivuoti objekta
	if(length(vids)==0){
		if(length(eids)>0){
			if(all(E(g)[eids]$selected) & length(eids)==1) 	g = PutActiveObject(oProgId=E(g)[eids]$ProgId, type="E", g=g)
			if(!all(E(g)[eids]$selected)) g = PutActiveObject(oProgId=NA, type="E", g=g)
			if(length(eids)>1) 	g = PutActiveObject(oProgId=NA, type="V", g=g)
		} else {
			# g = PutActiveObject(oProgId=NA, type="E", g=g)
		}
	}	
	g

}

PutActiveObject <- function(id=NULL, oProgId=NULL,  type=c("none", "V", "E", "G"), g=g){ 

	type = match.arg(type)

	### pasirupiname indeksais
	if(length(id)>0 && length(oProgId)==0){
		if(all(!is.na(id))){
			oProgId = switch(EXPR=type
				, none = NA
				, V = V(g)[id]$ProgId[1]
				, E = E(g)[id]$ProgId[1]
				, G = as.numeric(names(g$groups)[id][1])	
			)
		} else {
			oProgId = NA
		}
	}
	
	
	if(length(oProgId)>0){
		if(any(is.na(oProgId))){
			g$mode$ActiveObjectType = "none"
			g$mode$ActiveObjectProgId = 0
		} else {
			g$mode$ActiveObjectType = type
			g$mode$ActiveObjectProgId = oProgId
		}	
	} else {
		stop("The length of paramyters is 'PutActiveObject' is zero. This means tha in some place proper 'ProgId' was not found.")
	}

	g
}

GetActiveObject <- function(g=g){  #taisytina(Test)
	list(ProgId=g$mode$ActiveObjectProgId, type=g$mode$ActiveObjectType)	
}

PutViewObject <- function(id=NULL, oProgId=NULL, type=c("none", "V", "E", "G"), g=g){  
	type = match.arg(type)
	
	
	### pasirupiname indeksais
	if(length(id)>0 && length(oProgId)==0){
		if(all(!is.na(id))){
			oProgId = switch(EXPR=type
				, none = NA
				, V = V(g)[id]$ProgId[1]
				, E = E(g)[id]$ProgId[1]
				, G = as.numeric(names(g$groups)[id][1])	
			)
		} else {
			oProgId = NA
		}
	}
	
	if(length(oProgId)>0){
		if(any(is.na(oProgId))){
			g$mode$ViewObjectType = "none"
			g$mode$ViewObjectProgId = 0
		} else {
			g$mode$ViewObjectType = type
			g$mode$ViewObjectProgId = oProgId
		}	
	} else {
		stop("The length of paramyters is 'PutViewObject' is zero. This means tha in some place proper 'ProgId' was not found.")
	}
	g
}		
	
GetViewObject <- function(g=g){  #taisytina(Test)
	list(ProgId=g$mode$ViewObjectProgId, type=g$mode$ViewObjectType)	
}
	
GetObjectAttribute <- function(ids=NULL, oProgIds=NULL, Attr="name", type=c("none", "V", "E", "G"), g=g){  
	type = match.arg(type)	
	if(length(oProgIds)>0){
		ids = switch(EXPR=type
			, none = vector()
			, V = as.numeric(V(g)[ProgId %in% oProgIds])
			, E = as.numeric(E(g)[ProgId %in% oProgIds])
			, G = which(names(g$groups) %in% oProgIds)	
		)
	}
	
	
	value = rep(NA, length(ids))
	value = rep(NA, length(ids))
	idsNA = is.na(ids)	
	ids = ids[!idsNA]
	val = rep(NA, length(ids))	
	
	if(type=="none" | length(ids)==0){
		return(value)
	} else {
		if(type=="V"){
			val = get.vertex.attribute(graph=g, name=Attr, index=V(g)[ids])			
		}
		if(type=="E"){
			val = get.edge.attribute(graph=g, name=Attr, index=E(g)[ids])
		}
		if(type=="G"){# cia idetas tiesioginis id - gal reiktu naudoti ProgId - butu tvarkingiau
			for(i in seq_along(ids)){
				if(length(g$groups)>=ids[i]){
					val[i] = g$groups[[ids[i]]][[Attr]]
				}
			}
		}
		if(length(val)>0) value[!idsNA] = val
	}
	
	value
}

 # GetObjectAttribute(ids=1, Attr="name", type="V", g=g)
 # GetObjectAttribute(ids=c(1, NA, 500), Attr="name", type="V", g=g)
 # GetObjectAttribute(ids=c(1, NA, 500), Attr="ProgId", type="V", g=g)
 # GetObjectAttribute(ids=c(1, NA, 500), Attr="ProgId", type="E", g=g)
 # GetObjectAttribute(ids=c(1, NA, 500), Attr="Petras", type="E", g=g)
 # GetObjectAttribute(ids=c(1, NA, 500), Attr="name", type="G", g=g)

 # GetObjectAttribute(oProgIds=c(1, NA, 500), Attr="name", type="V", g=g)
 # GetObjectAttribute(oProgIds=c(1, NA, 500), Attr="name", type="E", g=g)

GetObjectLabel <- function(ids=NULL, oProgIds=NULL, type=c("none", "V", "E", "G"), g=g){ #taisytina(Test)
	type = match.arg(type)	
	if(length(ids)>0 && length(oProgIds)==0){
		oProgIds = switch(EXPR=type
			, none = rep(NA, length(ids))
			, V = V(g)[ids]$ProgId
			, E = E(g)[ids]$ProgId
			, G = as.numeric(names(g$groups)[ids])	
		)		
	}		
	if(length(oProgIds)==0) oProgIds = rep(NA, length(ids))
	
	value = rep("", length(oProgIds))	
	idsNA = is.na(oProgIds)
	oProgIds = oProgIds[!idsNA]	
	val = rep("", length(oProgIds))
	
	if(length(oProgIds)>0){
	
		name = GetObjectAttribute(oProgIds=oProgIds, Attr="name", type=type, g=g)
		name[is.na(name)] = ""
		label = GetObjectAttribute(oProgIds=oProgIds, Attr="label", type=type, g=g)
		label[is.na(label)] = ""
		
		
		if(type!="none") val[!is.na(name)] = paste(type,"(", oProgIds[!is.na(name)], ")", sep="")

		PrintName = ""
		if(any(nchar(name)>0)){
			if(any(paste(name)!=paste(oProgIds))) {
				PrintName = paste("<", name, ">", sep="")
			} else {
				PrintName = ""
			}
		}
		
		PrintLabel = ""
		if(any(nchar(label)>0)){
			if(any(paste(label)!=paste(name))) {
				PrintLabel = paste("[", label, "]", sep="")
			} else {
				PrintLabel = ""
			}
		}
		
		val = val%.%PrintName%.%PrintLabel
		
		if(length(val)>0) value[!idsNA] = val
		
	} 
	
	return(value)
}	

### debug(GetObjectLabel)
# GetObjectLabel(ids=2, type="V", g=g)
# GetObjectLabel(ids=c(1,2), type="V", g=g)
# GetObjectLabel(ids=9, type="E", g=g)
# GetObjectLabel(ids=1, type="G", g=g)

# GetObjectLabel(oProgIds=1, type="V", g=g)
# GetObjectLabel(oProgIds=1, type="E", g=g)
# GetObjectLabel(oProgIds=1, type="G", g=g)
###########################################
