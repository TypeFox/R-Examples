## navGraph internal class for objects containing path information
setClass(
		Class = "NG_path",
		representation = representation(
				path = "character",
				info = "character",
				graph = "character"
		),
		validity = function(object){		
			#if(length(object@path) == length(object@info)){# && length(object@info) == length(object@graph)){
			#	return(TRUE)
			#}else{
			#	return(FALSE)
			#}
		})
#
#setMethod(f = "[",
#		signature = "NG_path",
#		definition = function(x,i,j,drop){
#			if(all(length(i)==1, j %in% c("path","info","comment","graph","all"))){
#				if( j == "all"){ ## delete element
#					return(new("NG_path", path = x@path[i], info = x@info[i], graphName = x@graphName[i]))					
#				}else if(all(i > 0,any(match(j, c("path","info","graph","comment"))))){
#					if(j == "path"){
#						return(unlist(strsplit(x@path[i],' ')))
#					}else if(j %in% c("info","comment")){
#						return(x@info[i])
#					}else if(j == "graph"){
#						return(x@graphName[i])
#					}
#					
#				}else{
#					stop("[NG_path:'[' ]: wrong indices i and j")					
#				}		
#			}
#			
#		}
#)

setMethod(f = "show",
		signature = "NG_path",
		definition = function(object){
			
			if(length(object@path) != 0) {
				
				cat("NG_path object\n")
				for(i in 1:length(object@path)) {
					cat(paste("Path",i,"-----------------\n"))
					cat(paste("Graph:",object@graphName[i],"\n"))
					cat(paste("Path:",object@path[i],"\n"))
					cat(paste("Comment:",object@info[i],"\n"))
				}
			}else {
				cat("NG_path object is empty.\n")
			}
		}
)




.parsePath <- function(path){
	return(gsub('\\s+',' ',sub("\\s*\n*$","",sub("^\\s+","",path))))
}

.parsePath2Vec <- function(path){
	return(unlist(strsplit(.parsePath(path), split = " ")))
}

.pathGUI <- function(ngEnv) {
	if(!is.null(ngEnv$windowManager$paths)) {
		## bring window to front
		tkraise(ngEnv$windowManager$paths$tt)
	}else {
		
		
		ttwin <- tktoplevel(borderwidth=5)
		tktitle(ttwin) <- paste("Session ", ngEnv$ng_instance,", RnavGraph Paths", sep = '')
		
		ngEnv$windowManager$paths$tt <- ttwin
		tkbind(ttwin, '<Destroy>', function(){
					ngEnv$ngEnv$windowManager$paths <- NULL
				})
		
		
		f.active <- tkwidget(ttwin,"labelframe",text = "Active Path:")
		f.paths <- tkwidget(ttwin,"labelframe",text = "Saved Paths:")
		f.comments <- tkwidget(ttwin,"labelframe",text = "Comments on the selected path:")
		
		tkpack(f.active, side = "top", fill = "x",pady=2)
		tkpack(f.paths, side = "top", fill = "both", expand = TRUE)
		tkpack(f.comments, side = "top", fill = "x",pady=2)
		
		## active Path
		entry.activePath <- tkentry(f.active, bg = "white", relief = "sunken", textvariable = ngEnv$activePath)
		tkpack(entry.activePath, side="left", expand = TRUE, fill = "x", anchor="w")
		tkfocus(entry.activePath)
		
		## walk
		
		l.walk <- tklabel(f.active, text = "walk")# , activebackground="darkgrey")
		normalcolor <- tkcget(l.walk,background=NULL)
		tkbind(l.walk, '<Enter>', function(){tkconfigure(l.walk,background="darkgrey")})
		tkbind(l.walk, '<Leave>', function(){tkconfigure(l.walk,background=normalcolor)})
		tkbind(l.walk, '<Button-1>', function().walkPath(ngEnv, tclvalue(ngEnv$activePath)))
		
		
		
		## view
		l.view <- tklabel(f.active, text = "view")#, activebackground="darkgrey")
		tkbind(l.view, '<Enter>', function(){tkconfigure(l.view,background="darkgrey")})
		tkbind(l.view, '<Leave>', function(){tkconfigure(l.view,background=normalcolor)})
		tkbind(l.view, '<Button-1>', function().showPath(ngEnv, tclvalue(ngEnv$activePath)))
		
		
		
		## save
		l.save <- tklabel(f.active, text = "save")# , activebackground="darkgrey")
		tkbind(l.save, '<Enter>', function(){tkconfigure(l.save,background="darkgrey")})
		tkbind(l.save, '<Leave>', function(){tkconfigure(l.save,background=normalcolor)})
		
		
		tkpack(l.walk,l.view,l.save, side = "right", padx = "2")
		
		
		
		if(length(ngEnv$paths@path) == 1) {
			savedPaths <- tclVar(paste('{',ngEnv$paths@path,'}'))
			savedGraphs <-  tclVar(paste('{',ngEnv$paths@graph,'}'))
		}else {
			savedPaths <- tclVar(ngEnv$paths@path)
			savedGraphs <- tclVar(ngEnv$paths@graph)
		}
		
		## Saved Paths
#tclvalue(savedPaths) <- paths
		scr0 <- tkscrollbar(f.paths, repeatinterval=5, command=function(...){tkyview(tl,...);tkyview(tl2,...)})
		tl<-tklistbox(f.paths,height=12, selectmode="browse", listvariable=savedPaths, yscrollcommand=function(...){tkset(scr0,...)}, background="white", exportselection=0)
		tkpack(tl, side="left", fill = "both", expand = TRUE)
		tkpack(scr0, side="left", fill= "y")
		

		tl2<-tklistbox(f.paths,height=12, listvariable=savedGraphs, width = 20, selectmode="browse", yscrollcommand=function(...){tkset(scr0,...)}, background="white", exportselection=0)
		tkpack(tl2, side="right", fill = "both")
		
		
		## saved comments
		
		scr1 <- tkscrollbar(f.comments, repeatinterval=5, command=function(...)tkyview(txt,...))
		txt <- tktext(f.comments,bg="white",font="courier",yscrollcommand=function(...)tkset(scr1,...), width=30, height = 12, wrap="word")
		tkpack(txt, side="left", fill = "x", expand = TRUE)
		tkpack(scr1, side="left",fill="y")
		
		
		## initialize
		## TODO match active path
		pathEnv <- environment()
		tcl(tl,'activate',0)
		tcl(tl,'selection','set',0)
		tcl(tl,'see',0)
		tcl(tl2,'activate',0)
		tcl(tl2,'selection','set',0)
		tcl(tl2,'see',0)
		if(length(ngEnv$paths@path)>0){
			tkdelete(txt, '0.0', 'end')
			tkinsert(txt, "end", ngEnv$paths@info[1])
		}
		isBrowsing <- FALSE
		
		## interaction
		
		tkbind(tl, "<<ListboxSelect>>", function(){
					if(length(ngEnv$paths@path)>0){
						new <- as.numeric(tkcurselection(tl))
						tcl(tl2,'selection','clear',0,'end')
						tcl(tl2,'activate',new)
						tcl(tl2,'selection','set',new)
						tcl(tl2,'see',new)
						if(!pathEnv$isBrowsing){
							pathEnv$isBrowsing <- TRUE
							old <- as.numeric(tcl(tl,'index','active'))
							pathEnv$ngEnv$paths@info[old+1] <- gsub("\n$","",tclvalue(tkget(txt,"0.0","end")))    
						}
						tkdelete(txt, '0.0', 'end')
						tkinsert(txt, "end", ngEnv$paths@info[new+1])
					}
				})
		
		tkbind(tl2, "<<ListboxSelect>>", function(){
					if(length(ngEnv$paths@path)>0){
						new <- as.numeric(tkcurselection(tl2))
						tcl(tl,'selection','clear',0,'end')
						tcl(tl,'activate',new)
						tcl(tl,'selection','set',new)
						tcl(tl,'see',new)
						if(!pathEnv$isBrowsing){
							pathEnv$isBrowsing <- TRUE
							old <- as.numeric(tcl(tl2,'index','active'))
							pathEnv$ngEnv$paths@info[old+1] <- gsub("\n$","",tclvalue(tkget(txt,"0.0","end")))        
						}
						tkdelete(txt, '0.0', 'end')
						tkinsert(txt, "end", ngEnv$paths@info[new+1])
					}
				})
		
		tkbind(tl, "<ButtonRelease-1>", function(){pathEnv$isBrowsing <- FALSE})
		tkbind(tl2, "<ButtonRelease-1>", function(){pathEnv$isBrowsing <- FALSE})
		
		
		tkbind(l.save, '<Button-1>', function(){
					.saveComment()
					aPath <- as.character(tclvalue(ngEnv$activePath))
					aGraph <- as.character(tclvalue(ngEnv$activePathGraph))
					
					isNew <- FALSE
					
					
					lP <- ngEnv$paths@path %in% aPath
					lG <- ngEnv$paths@graph %in% aGraph
					
					if(any(lP & lG)){
						i <- which((lP & lG) == TRUE)[1]
					}else {
						isNew <- TRUE
					}
					
					if(!isNew) {
						## old Path
						tcl(tl,'selection','clear',0,'end')
						tcl(tl,'activate',i-1)
						tcl(tl,'selection','set',i-1)
						tcl(tl,'see',i-1)
						tcl(tl2,'selection','clear',0,'end')
						tcl(tl2,'activate',i-1)
						tcl(tl2,'selection','set',i-1)
						tcl(tl2,'see',i-1)
						tkdelete(txt, '0.0', 'end')
						tkinsert(txt, "end", ngEnv$paths@info[i])
					}else {
						## add new path
						ngEnv$paths@path <- c(ngEnv$paths@path,aPath)
						ngEnv$paths@graph <- c(ngEnv$paths@graph,aGraph)
						ngEnv$paths@info <- c(ngEnv$paths@info,"")
						if(length(ngEnv$paths@path) == 1) {
							tclvalue(savedPaths) <- paste('{',ngEnv$paths@path,'}')
							tclvalue(savedGraphs) <-  paste('{',ngEnv$paths@graph,'}')
						}else {
							tclvalue(savedPaths) <- ngEnv$paths@path
							tclvalue(savedGraphs) <- ngEnv$paths@graph
						}
						tcl(tl,'selection','clear',0,'end')
						tcl(tl,'selection','set','end')
						tcl(tl,'activate','end')
						tcl(tl,'see','end')
						tcl(tl2,'selection','clear',0,'end')
						tcl(tl2,'selection','set','end')
						tcl(tl2,'activate','end')
						tcl(tl2,'see','end')
						tkdelete(txt, '0.0', 'end')
					}
				})
		
		
		
		.saveComment <- function(){
			i <- as.numeric(tkcurselection(tl))+1
			pathEnv$ngEnv$paths@info[i] <- gsub("\n$","",tclvalue(tkget(txt,"0.0","end")))
			return(i)
		}
		
		
		tkbind(tl,"<Key-Delete>",function().delCurrent())
		tkbind(tl2,"<Key-Delete>",function().delCurrent())
		
		tkbind(tl,"<Key-BackSpace>",function().delCurrent())
		tkbind(tl2,"<Key-BackSpace>",function().delCurrent())
		
		.delCurrent <- function(){
			i <- as.numeric(tkcurselection(tl))+1
			if(length(i)>0){
				ngEnv$paths@info <- ngEnv$paths@info[-i]
				ngEnv$paths@path <- ngEnv$paths@path[-i]		
				ngEnv$paths@graph <- ngEnv$paths@graph[-i]
				if(length(ngEnv$paths@path) == 1) {
					tclvalue(savedPaths) <- paste('{',ngEnv$paths@path,'}')
					tclvalue(savedGraphs) <-  paste('{',ngEnv$paths@graph,'}')
				}else {
					tclvalue(savedPaths) <- ngEnv$paths@path
					tclvalue(savedGraphs) <- ngEnv$paths@graph
				}
				
				if(i > length(ngEnv$paths@path)) {
					i <- length(ngEnv$paths@path )
				}
				tcl(tl,'selection','clear',0,'end')
				tcl(tl,'selection','set',i-1)
				tcl(tl,'activate',i-1)
				tcl(tl2,'selection','clear',0,'end')
				tcl(tl2,'selection','set',i-1)
				tcl(tl2,'activate',i-1)
				tkdelete(txt, '0.0', 'end')
				tkinsert(txt, "end", ngEnv$paths@info[i])
			} else {
				tkdelete(txt, '0.0', 'end')
			}
		}
		
		
		tkbind(tl,"<Double-Button-1>",function().save2active())
		tkbind(tl2,"<Double-Button-1>",function().save2active())
		
		tkbind(tl,"<KeyPress-KP_Enter>",function().save2active())
		tkbind(tl2,"<KeyPress-KP_Enter>",function().save2active())
		
		tkbind(tl,"<KeyPress-Return>>",function().save2active())
		tkbind(tl2,"<KeyPress-Return>>",function().save2active())
		
		
		
		.save2active <- function(){
			i <- as.numeric(tkcurselection(tl))
			if(length(i)>0){  
				print(i)
				tclvalue(ngEnv$activePath) <- ngEnv$paths@path[i+1]
				tclvalue(ngEnv$activePathGraph) <- ngEnv$paths@graph[i+1]
#		tkdelete(txt, '0.0', 'end')
#		tkinsert(txt, "end", comments[i+1])
			}
		}
		
		
		
		
		closeWin <- function(){
			.saveComment()
			tkdestroy(ttwin)
		}
		tcl("wm", "protocol", ttwin, "WM_DELETE_WINDOW", function()closeWin())
	}
	
}





## carfule what argument the function takes (vector or string?)
.isPath <- function(graph,path){
	
	t.nodes <- nodeNr(graph,path)
	
	edgeM <- matrix(tail(head(rep(t.nodes,each=2), n=-1L), n=-1L),ncol=2, byrow = 2)	
	
	l.p <- apply(edgeM,1, FUN = function(x){any(x[2] == adjacent(graph,x[1],'node',retNr=TRUE))})
	
	
	if(any(is.na(l.p))){
		tkmessageBox(message = "at least one node does not exist")
		stop("your path is not correct.\n")
	}else{
		if(all(l.p)){
			return(TRUE)
		}else{
			tkmessageBox(message = "path does not exist")
			stop("path does not exist.\n")
		}
	}
}


## path is a string or vector
.isPathOnCanvas <- function(ngEnv, path) {
	
	if(length(path) == 1) {
		nodes <- .parsePath2Vec(path)
	}else {
		nodes <- path
	}
	n <- length(nodes)
	edgeMatrix <- cbind(nodes[1:(n-1)],nodes[2:n])
	edgeExists <- apply(edgeMatrix,1, FUN = function(row){
				length(.tcl2str(tcl(ngEnv$canvas,'find','withtag',paste('edge && ',row[1],' && ', row[2]))))
			})
	if(0 %in% edgeExists){
		tkmessageBox(message = "path does not exist")
		stop("path does not exist.\n")
	}else{
		return(TRUE)
	}
}

## path is a string 
.walkPath <- function(ngEnv, path) {
	if(path != "") {
		.isPathOnCanvas(ngEnv, path)
		
		if(length(path) == 1) {
			nodes <- .parsePath2Vec(path)
		}else {
			nodes <- path
		}
		
		n <- length(nodes)
		#browser()
		
		
		## first color the total path
		for(i in 1:(n-1)) {
			tkitemconfigure(ngEnv$canvas, paste('edge && ', nodes[i], ' && ', nodes[i+1]),
					fill = ngEnv$settings@color@path)
		}
		
		## jump to first node
		ngEnv$bulletState$from <- nodes[1]
		ngEnv$bulletState$percentage <- 0 
		.updatePlots(ngEnv)
		.arriveAtNode(ngEnv)
		xynode <- .tcl2xy(tkcoords(ngEnv$canvas,paste('node && ',nodes[1])))
		xybullet <- .tcl2xy(tkcoords(ngEnv$canvas,'bullet'))
		dxy <- xynode-xybullet
		tkmove(ngEnv$canvas, 'bullet',dxy[1],dxy[2])
		
		print(nodes)
		for(i in 2:n) {
#			ngEnv$bulletState$from <- nodes[i]
			ngEnv$bulletState$to <- nodes[i]
			ngEnv$bulletState$percentage <- 0
			
			tkitemconfigure(ngEnv$canvas, 'edge', width = ngEnv$settings@display@lineWidth)
			
			tkitemconfigure(ngEnv$canvas, paste('edge && ', nodes[i-1], ' && ', nodes[i]),
					fill = ngEnv$settings@color@path,
					width = ngEnv$settings@display@highlightedLineWidth)
			
			
			.walkEdge(ngEnv, path = !(i == n))
		}		
	}
}


.showPath <- function(ngEnv,path) {
	.isPathOnCanvas(ngEnv, path)
	if(length(path) == 1) {
		nodes <- .parsePath2Vec(path)
	}else {
		nodes <- path
	}
	n <- length(nodes)
	edgeMatrix <- cbind(nodes[1:(n-1)],nodes[2:n])
	
	
	for(i in 1:nrow(edgeMatrix)) {
		tkitemconfigure(ngEnv$canvas, 'edge', width = ngEnv$settings@display@lineWidth)
		tkitemconfigure(ngEnv$canvas, paste('edge && ', edgeMatrix[i,1], ' && ', edgeMatrix[i,2]),
				fill = ngEnv$settings@color@path,
				width = ngEnv$settings@display@highlightedLineWidth)
		tcl('update','idletasks')
		Sys.sleep(0.6)
	}
}












setMethod(f = "ng_get",
		signature = "NG_path",
		definition = function(obj, what=NULL, ...){
			possibleOptions <- c("path","graph","comment")
			
			if(is.null(what)){
				cat("Get what? Possible options are: ")
				cat(paste(possibleOptions, collapse = ", "))
				cat("\n")
			}else{
				if(any(is.na(match(what,possibleOptions)))){
					stop(paste("[ng_get] object",what,"is not defined."))
				}
				
				if(what == "path"){
					return(obj@path)
				}else if(what == "graph"){
					return(obj@graph)
				}else if(what == "comment"){
					return(obj@info)
				}
			}
		})






