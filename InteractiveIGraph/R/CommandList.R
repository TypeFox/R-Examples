
CommandList <- list()


### Special functional commands: save, quite, reset, delete,... ### 
CommandList[["Del"]] <- list(label="Delete", description = "Delete any active object.")
CommandList[["Del"]]$FUN <- function(g){
	
	AO = GetActiveObject(g)
	
	if(AO$type!="none"){
		if(AO$type=="G" & AO$ProgId>0){
			g <- GroupDelete(gProgId=AO$ProgId, g=g)			
		}	
		if(AO$type=="V" & AO$ProgId>0){
			vid = as.numeric(V(g)[ProgId == AO$ProgId])
			g <- delete.vertices(g, vid)
			g <- XY.norm(g)
		}	
		if(AO$type=="E" & AO$ProgId>0){
			eid = as.numeric(E(g)[ProgId == AO$ProgId])
			g <- delete.edges(g, eid)
		}	
		
		
		
		if(class(g)[1]!="InteractiveIGraph") class(g) <- c("InteractiveIGraph", class(g))
		g <- PutActiveObject(id=NA, type="none", g=g)
		g <- PutViewObject(id=NA, type="none", g=g)			
	}
		
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 
CommandList[["q"]] <- list(label="Quite", description = "Quits the function, returns last value of graph.")
CommandList[["q"]]$FUN <- function(g){
	g$ExtraParam$quite = TRUE
	invisible(g)
} 
CommandList[["ctrl-S"]] <- list(label="Save", description = "Saves graph and dump graph to text file in active directory.")
CommandList[["ctrl-S"]]$FUN <- function(g){	
	# browser()
	dump(c("g"), file="_Graph_"%.%g$ExtraParam$GName%.%".R")
	g$ExtraParam$save = TRUE
	# assign("SavedGraph", value=g, envir=parent.frame(1))
	invisible(g)
} 
CommandList[["ctrl-R"]] <- list(label="Reset", description = "Reset graph from las saving.")
CommandList[["ctrl-R"]]$FUN <- function(g){
	#browser()
	g <- get("SavedGraph", envir=parent.frame(1))
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 
CommandList[["ctrl-C"]] <- list(label="Command", description = "Waits for input of the command.")
CommandList[["ctrl-C"]]$FUN <- function(g){
	g$mode$input <- TRUE
	g$mode$InputType <- "Command"
	g$ExtraParam$input <- ""
	g$Menu$MenuList$BottomMenu$MenuLines$Input$active <- TRUE
	g$Menu$MenuList$BottomMenu$MenuLines$Input$MenuItems[[1]]$label <- g$mode$InputType %.%": "

	ClearBottomMenuPlot(g)
	g$Menu$ActiveMenuList <- list(BottomMenu=plot.Menu(g$Menu$MenuList$BottomMenu, x = par("usr")[1], y = par("usr")[3]))

	cat(g$mode$InputType %.%": ")	
	invisible(g)
}
CommandList[["ctrl-F"]] <- list(label="Change attributes", description = "Waits for input of the attributes specification (color='green')")
CommandList[["ctrl-F"]]$FUN <- function(g){
		
	g$mode$input <- TRUE
	g$mode$InputType <- "Attributes"
	g$ExtraParam$input <- ""	
	g$Menu$MenuList$BottomMenu$MenuLines$Input$active <- TRUE
	g$Menu$MenuList$BottomMenu$MenuLines$Input$MenuItems[[1]]$label <- g$mode$InputType %.%": "		
		
	ClearBottomMenuPlot(g)
	g$Menu$ActiveMenuList <- list(BottomMenu=plot.Menu(g$Menu$MenuList$BottomMenu, x = par("usr")[1], y = par("usr")[3]))

	cat(g$mode$InputType %.%": ")	
	invisible(g)		
}
CommandList[["ctrl-B"]] <- list(label="SavePdf", description = "Saves picture as pdf.")
CommandList[["ctrl-B"]]$FUN <- function(g){
	
	# lets find unused name:	
	FileTemp = "_"%.%g$ExtraParam$GName
	Ind = TRUE
	i=0
	while(Ind &(i<1000)){
		i = i + 1
		Ind = file.exists(FileTemp %.% "_" %.% i %.% ".pdf")
	}
	
	if(Ind){
		print("Warning: program did not manage to find free name for graph file:
			it tested the range " %.% FileTemp %.% "_{1 to 1000}.")
		return(invisible(g))
	} else {	
		d2file(file=FileTemp %.% "_" %.% i %.% ".pdf")
	}
	
	invisible(g)		
}

### navigation ###
CommandList[["+"]] <- list(label="zoom(+)", description = "zoom(+)")
CommandList[["+"]]$FUN <- function(g){
	g$PlotParam$xlim <- g$PlotParam$xlim - c(-0.1,0.1)
	g$PlotParam$ylim <- g$PlotParam$ylim - c(-0.1,0.1)
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 
CommandList[["-"]] <- list(label="zoom(-)", description = "zoom(-)")
CommandList[["-"]]$FUN <- function(g){
	g$PlotParam$xlim <- g$PlotParam$xlim + c(-0.1,0.1)
	g$PlotParam$ylim <- g$PlotParam$ylim + c(-0.1,0.1)
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 
CommandList[["Left"]] <- list(label="Left", description = "Move to left")
CommandList[["Left"]]$FUN <- function(g){
	g$PlotParam$xlim <- g$PlotParam$xlim-0.1
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 
CommandList[["Right"]] <- list(label="Right", description = "Move to right")
CommandList[["Right"]]$FUN <- function(g){
	g$PlotParam$xlim <- g$PlotParam$xlim+0.1
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 
CommandList[["Up"]] <- list(label="Up", description = "Move up")
CommandList[["Up"]]$FUN <- function(g){
	g$PlotParam$ylim <- g$PlotParam$ylim+0.1
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 
CommandList[["Down"]] <- list(label="Down", description = "Move down")
CommandList[["Down"]]$FUN <- function(g){
	g$PlotParam$ylim <- g$PlotParam$ylim-0.1
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 


### active / not active ###
CommandList[["a"]] <- list(label="(de)activate", description = "Activates or activates vertices.")
CommandList[["a"]]$FUN <- function(g){

	VO = GetViewObject(g)
	if(VO$type=="none"){
		V(g)[selected]$active <- !V(g)[selected]$active
	} else {
		if(VO$type=="V"){
			V(g)[ProgId==VO$ProgId]$active <- !V(g)[ProgId==VO$ProgId]$active
		}
		if(VO$type=="G"){
			g$groups[[as.character(VO$ProgId)]]$active = !all(g$groups[[as.character(VO$ProgId)]]$active)
		}
	}
	
	
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 
CommandList[["A"]] <- list(label="Activate only", description = "Deactivate all vertices but selected.")
CommandList[["A"]]$FUN <- function(g){
	V(g)$active <- FALSE
	V(g)[selected]$active <- TRUE
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 
CommandList[["ctrl-A"]] <- list(label="Select all", description = "Select all vertices.")
CommandList[["ctrl-A"]]$FUN <- function(g){
	V(g)[!hidden]$selected <- TRUE
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 


### groups ### 
CommandList[["g"]] <- list(label="Groupe", description = "Create new group, or join group if in group view mode.")
CommandList[["g"]]$FUN <- function(g){
	vProgIds = V(g)[selected]$ProgId
	if(length(vProgIds)>0){
		VO = GetViewObject(g)
		if(VO$type=="none"){
			### sukuriam grupe 
			g <- GroupCreate(g=g, vProgIds=vProgIds, selected=TRUE)
		 	g <- PutActiveObject(id=length(g$groups), type="G", g=g)
		} else {
			if(VO$type=="G"){	
				# browser()
				### pridedam arba atmetam elementus is grupes
				add = setdiff(vProgIds, g$groups[[as.character(VO$ProgId)]]$MemberProgIds)
				drop = intersect(vProgIds, g$groups[[as.character(VO$ProgId)]]$MemberProgIds)
				g$groups[[as.character(VO$ProgId)]]$MemberProgIds = setdiff(union(g$groups[[as.character(VO$ProgId)]]$MemberProgIds, add),drop)				
			}
		}		
		g$ExtraParam$replot <- TRUE
	}
	invisible(g)
} 
CommandList[["G"]] <- list(label="Select active group members", description = "Select active group members")
CommandList[["G"]]$FUN <- function(g){
	AO = GetActiveObject(g)
	
	if(AO$type!="none"){
		if(AO$type=="G" & AO$ProgId>0){
			G = g$groups[[as.character(AO$ProgId)]]
			V(g)$selected = FALSE	
			V(g)[ProgId %in% G$MemberProgIds]$selected = TRUE			
		}
	}
	g$ExtraParam$replot <- TRUE
	invisible(g)		
}
CommandList[["ctrl-G"]] <- list(label="Menu all groups", description = "Put all groups in bottom menu.")
CommandList[["ctrl-G"]]$FUN <- function(g){

		g <- MenuClearExtraItems(MenuLine='G', g=g)	
	
		Items = list()
		for(i in seq_along(g$groups)){		
			Items[[as.character(i)]] = list(label="<"%.%GetObjectAttribute(ids=i, Attr="name", type="G", g=g)%.%">", command="g <<- MenuSelectGroup(id="%.%i%.%", g=g)", RegionParams=list(YBufCof=0.2), RecParams=list(lwd = NA, border=NA), TextParams=list(cex=0.8, col=2))
		}					
		g <- MenuAddItems(Items,  MenuLine='G', g=g)	
	

	invisible(g)		
}
CommandList[["b"]] <- list(label="Block", description = "Block(unblock) active group.")
CommandList[["b"]]$FUN <- function(g){

	AO = GetActiveObject(g)
	# browser()
	if(AO$type=="G"){
		if(g$groups[[as.character(AO$ProgId)]]$closed){ # atidarom
			g <- GroupOpen(gProgId=AO$ProgId, g=g)
		} else { # uzdarom				
			g <- GroupClose(gProgId=AO$ProgId, g=g)			
		}	
	}	
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 


### modes / Menu ###
CommandList[["s"]] <- list(label="Mode(select)", description = "Turn on(off) select mode.")
CommandList[["s"]]$FUN <- function(g){
	g$mode$select  <- !g$mode$select
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 
CommandList[["i"]] <- list(label="Mode(AllEdges)", description = "Turn on(off) AllEdges mode.")
CommandList[["i"]]$FUN <- function(g){
	g$mode$AllEdges = !g$mode$AllEdges
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 
CommandList[["M"]] <- list(label="Print mode list", description = "Prints the list of modes.")
CommandList[["M"]]$FUN <- function(g){
	print(g$mode)
	invisible(g)
} 
CommandList[["m"]] <- list(label="Print all commands", description = "Prints the list of commands.")
CommandList[["m"]]$FUN <- function(g){
	PrintCommandList()	
	invisible(g)
} 
CommandList[["h"]] <- list(label="Hide bottom meniu", description = "Hides(unhides) bottom meniu.")
CommandList[["h"]]$FUN <- function(g){
	g$Menu$MenuList$BottomMenu$active = !g$Menu$MenuList$BottomMenu$active	
	g$ExtraParam$replot <- TRUE
	invisible(g)
} 






### vertex manipulations ###
CommandList[["v"]] <- list(label="Select vertices by edge", description = "Select vertices that associates with selected edges.")
CommandList[["v"]]$FUN <- function(g){
	eids = E(g)[selected]
	V(g)[adj(eids)]$selected = TRUE	

	g$ExtraParam$replot <- TRUE
	invisible(g)		
}
CommandList[["V"]] <- list(label="Deselect all vertices", description = "Remove any selections of the vertices.")
CommandList[["V"]]$FUN <- function(g){
	V(g)$selected = FALSE
	g$ExtraParam$replot <- TRUE
	invisible(g)		
}
CommandList[["ctrl-V"]] <- list(label="Add vertex", description = "Create new vertex in the place of past click.")
CommandList[["ctrl-V"]]$FUN <- function(g){
	ProgId = max(V(g)$ProgId)+1
	StartX = get("StartX", envir=parent.frame(1))
	StartY = get("StartY", envir=parent.frame(1))
	attr = list(name="n("%.%ProgId%.%")", label="n("%.%ProgId%.%")", ProgId=ProgId, x=StartX, y=StartY
			, active=TRUE, hidden=FALSE, selected=TRUE
			, ProgType = "normal", shape = "circle", size=15, color='SkyBlue2'
			, frame.color = NA)
	g <- add.vertices(g, 1, attr=attr)		
	g <- XY.norm(g)
	
	
	if(class(g)[1]!="InteractiveIGraph") class(g) <- c("InteractiveIGraph", class(g))
	g$ExtraParam$replot <- TRUE
	invisible(g)		
}


### edge manipulations ###
CommandList[["e"]] <- list(label="Select edges by vertices", description = "Select edges that associates with selected vertices.")
CommandList[["e"]]$FUN <- function(g){
	vids = V(g)[selected]
	E(g)[adj(vids)]$selected = TRUE	

	g$ExtraParam$replot <- TRUE
	invisible(g)			
}
CommandList[["E"]] <- list(label="Deselect all edges", description = "Remove any selections of the edges.")
CommandList[["E"]]$FUN <- function(g){
	E(g)$selected = FALSE
	g$ExtraParam$replot <- TRUE
	invisible(g)		
}
CommandList[["ctrl-E"]] <- list(label="Add edge", description = "Create new edge, from selected vertices to active vertex.")
CommandList[["ctrl-E"]]$FUN <- function(g){
	### adds edges with active vertex (if any) fiwh all other selected vertices
	# browser()
	AO = GetActiveObject(g)
	if(AO$type=="V"){
		AOvid = as.numeric(V(g)[ProgId==AO$ProgId])
		vids = setdiff(as.numeric(V(g)[selected & ProgType=="normal"]), AOvid)
		
		if(length(vids)>0){
			### pridadame krastine
			EfromDF = as.data.frame(cbind(V1=vids, V2=AOvid))
			eProgId = max(E(g)$ProgId)+seq_along(vids)
			EfromDF$ProgId = eProgId
			EfromDF$name = eProgId
				
			EfromDF$hidden = FALSE	
			EfromDF$selected = FALSE	# jungciu pazymejimas
			EfromDF$ProgType = "normal"
			EfromDF$active = TRUE
		
			EfromDF$lty = 1
			EfromDF$width = 1
			EfromDF$color = "darkgrey" 
			EfromDF$arrow.size = 0.7 ### nebutinas, bet norimas
			for(i in 1:nrow(EfromDF)){
				g = add.edges(g, as.vector(as.matrix(EfromDF[i,1:2])), attr=as.list(EfromDF[i,-(1:2)]))
			}	
			
			## reikia susitvarkyti su grupem
		}
	}
	if(class(g)[1]!="InteractiveIGraph") class(g) <- c("InteractiveIGraph", class(g))
	g$ExtraParam$replot <- TRUE
	invisible(g)		
}
CommandList[["d"]] <- list(label="Activate selected edges", description = "Activate selected edges")
CommandList[["d"]]$FUN <- function(g){
	E(g)[selected]$active = TRUE	
	g$ExtraParam$replot <- TRUE
	invisible(g)			
}
CommandList[["D"]] <- list(label="Deactivate selected edges", description = "Deactivate selected edges")
CommandList[["D"]]$FUN <- function(g){
	E(g)[selected]$active = FALSE	
	g$ExtraParam$replot <- TRUE
	invisible(g)		
}


### trees / paths ###
CommandList[["t"]] <- list(label="minimum.spanning.tree", description = "Selects all edges allond minimum spanning tree.")
CommandList[["t"]]$FUN <- function(g){

	gi = induced.subgraph(g, vids=V(g)[active & !hidden])
	gitree = minimum.spanning.tree(gi)
	eProgIds = E(gitree)$ProgId
	E(g)$selected = FALSE
	E(g)[ProgId %in% eProgIds]$selected = TRUE


	g$ExtraParam$replot <- TRUE
	invisible(g)		
}
CommandList[["p"]] <- list(label="shortest.paths(form)", description = "Finds shortest path from from selected vertices to active vertex.")
CommandList[["p"]]$FUN <- function(g){
	
	if(GetActiveObject(g)$type=="V" & sum(V(g)$selected)>1){
		ShortestPath = GetShortestPath(mode="in", g=g)

		eProgIds = ShortestPath$eProgIds
		if(length(eProgIds)>0){
			E(g)$selected = FALSE
			E(g)[ProgId %in% eProgIds]$selected = TRUE	
			
		} else {
			g = MsgToLogObj(" - No path found - ", g=g, add=TRUE)
		}
	} else {
		g = MsgToLogObj("Which vertices should be joined? Program searches from selected vertices to active one. ", g=g, add=TRUE)
	}
	g$ExtraParam$replot <- TRUE	
	invisible(g)		
}
CommandList[["P"]] <- list(label="shortest.paths(to)", description = "Finds shortest path from from active vertex to selected vertices.")
CommandList[["P"]]$FUN <- function(g){
	
	ShortestPath = GetShortestPath(mode="out", g=g)
	eProgIds = ShortestPath$eProgIds
	if(length(eProgIds)>0){
		E(g)$selected = FALSE
		E(g)[ProgId %in% eProgIds]$selected = TRUE			
	} else {
		g = MsgToLogObj(" - No path found - ", g=g, add=TRUE)
	}
	g$ExtraParam$replot <- TRUE
	invisible(g)		
}
CommandList[["ctrl-P"]] <- list(label="shortest.paths(none)", description = "Finds shortest path regardless directions.")
CommandList[["ctrl-P"]]$FUN <- function(g){
	
	ShortestPath = GetShortestPath(mode="all", g=g)
	eProgIds = ShortestPath$eProgIds
	if(length(eProgIds)>0){
		E(g)$selected = FALSE
		E(g)[ProgId %in% eProgIds]$selected = TRUE	
	} else {
		g = MsgToLogObj(" - No path found - ", g=g, add=TRUE)
	}
	g$ExtraParam$replot <- TRUE
	invisible(g)		
}


### appierences ### 
CommandList[["w"]] <- list(label="Adjust whide", description = "Adjust vertex whide to fit label.")
CommandList[["w"]]$FUN <- function(g){

	V(g)[active & !hidden]$size = ceiling((strwidth("  ") + strwidth(V(g)[active & !hidden]$label))*100) + 1
	V(g)[active & !hidden]$size2 = ceiling((strheight(" ") + strheight(V(g)[active & !hidden]$label))*100) + 1

	g$ExtraParam$replot <- TRUE
	invisible(g)		
}
CommandList[["l"]] <- list(label="Layout", description = "Redo layout (layout.outo)")
CommandList[["l"]]$FUN <- function(g){

	gi = induced.subgraph(g, vids=V(g)[!hidden])
	gi = remove.vertex.attribute(gi, "x")
	gi = remove.vertex.attribute(gi, "y")

	coords <- layout.auto(gi, dim=2)
	vProgIds = V(gi)$ProgId
	
	V(g)[ProgId %in% vProgIds]$x = coords[,1]
	V(g)[ProgId %in% vProgIds]$y = coords[,2]
	
	g <- XY.norm(g)
	
	g$ExtraParam$replot <- TRUE
	invisible(g)		
}
CommandList[["o"]] <- list(label="Plot Active Regions", description = "Marks the regions around vertices that will initiates selection if clicked.")
CommandList[["o"]]$FUN <- function(g){
	PlotActiveRegionsGraphs(g=g)
	g$ExtraParam$replot <- FALSE
	invisible(g)
} 



