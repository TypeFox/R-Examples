
## Convert tcl values (strings) into R variables (vectors)
.tcl2num <- function(x) {
	as.numeric(unlist(strsplit(tclvalue(x), split = ' ')))
}
.tcl2xy <- function(x) {
	coord <-as.numeric(unlist(strsplit(tclvalue(x), split = ' ')))
	c(mean(coord[c(1,3)]), mean(coord[c(2,4)]))
}
.tcl2str <- function(x) {
	unlist(strsplit(tclvalue(x), split = ' '))
}


## event
.registerButtonClick <- function(ngEnv,x,y) {
	ngEnv$mouseHistory$x <- as.numeric(x)
	ngEnv$mouseHistory$y <- as.numeric(y)
	ngEnv$mouseHistory$xClick <- ngEnv$mouseHistory$x
	ngEnv$mouseHistory$yClick <- ngEnv$mouseHistory$y
	
	ngEnv$mouseHistory$Click <- TRUE 
}


## graph & data visualization connection
## switching a graph calls .closePlots and .initializePlots
## changing the bullet state invokes .updatePlots
.initializePlots <- function(ngEnv) {
	assign(ngEnv$selViz,sapply(get(ngEnv$selViz,envir = ngEnv),function(viz)initializeViz(viz,ngEnv)),envir = ngEnv)	
}

.updatePlots <- function(ngEnv) {
	ngEnv$changedPlots <- TRUE
	assign(ngEnv$selViz,sapply(get(ngEnv$selViz,envir = ngEnv),function(viz)updateViz(viz,ngEnv)),envir = ngEnv)
}

.closePlots <- function(ngEnv) {
	assign(ngEnv$selViz,sapply(get(ngEnv$selViz,envir = ngEnv),function(viz)closeViz(viz,ngEnv)),envir = ngEnv)
}


## Switch the graph via the menu entry
.switchGraph <- function(ngEnv,graphNr){
	
	graphName <- paste('graph',graphNr,sep = '')
	tkadd(ngEnv$graphMenu,"command", label = get(paste('graph',graphNr,sep = ''),envir = ngEnv)@name,
			command=function(){
				if(graphName != ngEnv$cGraph){

					.closePlots(ngEnv)
					## Graph
					## update graph modifications

					.canvas2graph(ngEnv)
					## save current graph


					assign(ngEnv$cGraph, ngEnv$graph, envir = ngEnv)
					

					## add new graph
					ngEnv$cGraph <- graphName
					ngEnv$graph <- get(graphName, envir = ngEnv)
					ngEnv$selViz <-  sub('graph','viz',graphName)


					## reset state
					if(numNodes(ngEnv$graph@graph)>0){
						
						if(length(get(ngEnv$selViz,envir = ngEnv)) == 0) {
							ngEnv$bulletState$from <- nodes(ngEnv$graph@graph)[1]
						} else { 
							ngEnv$bulletState$from <- get(ngEnv$selViz,envir = ngEnv)[[1]]@from
						}
						ngEnv$bulletState$to <- ''
						ngEnv$bulletState$percentage <- 0
						.graph2canvas(ngEnv)
						.highlightAdj(ngEnv)
						
					}else {
						w <- .tcl2num(tcl('winfo','width',ngEnv$canvas))
						h <- .tcl2num(tcl('winfo','height',ngEnv$canvas))
						
						tkdelete(ngEnv$canvas,'all')
						tkcreate(ngEnv$canvas, 'text',w/2,h/2, text = 'No Nodes!', fill = ngEnv$settings$color$labels)	
						
					}

					.initializePlots(ngEnv)						
					
				}
				
			})
}



## Settings Menue
.settingsMenue <- function(ngEnv){
	
	.updateValues <- function(ngEnv, obj){
		vname <- paste("v.",obj, sep="")  ## tcl variable name
		value <- as.numeric(tclvalue(ngEnv$windowManager$settings[[vname]]))
		
		if(obj == "NSteps") {
			ngEnv$settings@interaction@NSteps <- value
		}else if(obj == "animationTime") {
			ngEnv$settings@interaction@animationTime <- value
		}else if(obj == "dragSelectRadius") {
			ngEnv$settings@interaction@dragSelectRadius <- value 
		}else if(obj == "labelDistRadius") {
			ngEnv$settings@interaction@labelDistRadius <- value 
		}else if(obj == "bulletRadius") {
			ratio <- value/ngEnv$settings@display@bulletRadius
			xy <- .tcl2xy(tkcoords(ngEnv$canvas,'bullet'))		
			tkitemscale(ngEnv$canvas, 'bullet', xy[1], xy[2], ratio, ratio)
			ngEnv$settings@display@bulletRadius <- value
			
		}else if(obj == "nodeRadius") {
			ratio <- value/ngEnv$settings@display@nodeRadius		
			tcl('scaleNodes',ngEnv$canvas,ratio)
			ngEnv$settings@display@nodeRadius <- value
		}else if(obj == "lineWidth") {
			ngEnv$settings@display@lineWidth <- value
			tkitemconfigure(ngEnv$canvas, "edge", width = value)
			.normalState(ngEnv)
		}else if(obj == "highlightedLineWidth") {
			ngEnv$settings@display@highlightedLineWidth <- value
			.normalState(ngEnv)
		}else {
			
		}
	}
	
	if(!is.null(ngEnv$windowManager$settings)) {
		## bring window to front
		tkraise(ngEnv$windowManager$settings$tt)
	}else {
		boldfont <- tkfont.create(weight="bold")
		
		tt <- tktoplevel(, borderwidth = 5)  ## main tk window
		tktitle(tt) <- tktitle(tt) <- paste("Session ", ngEnv$ng_instance,", RnavGraph Settings", sep = '')
		
		ngEnv$windowManager$settings$tt <- tt
		
		tkbind(tt, '<Destroy>', function(){
					ngEnv$ngEnv$windowManager$settings <- NULL
				})
		
		colframe <- tkframe(tt, relief = "groove", borderwidth = 3)
		tkpack(colframe, side = "left", padx = 5, pady = 5, fill = "y")
		
		intframe <- tkframe(tt, relief = "groove", borderwidth = 3)
		tkpack(intframe, side = "left", padx = 5, pady = 5, fill = "y")
		
		dispframe <- tkframe(tt , relief = "groove", borderwidth = 3) 
		tkpack(dispframe, side = "left", padx = 5, pady = 5, fill = "y")
		
		sapply(slotNames(ngEnv$settings@color), FUN = function(obj){
					fname <- paste("f.",obj,sep = "")
					assign(fname, tkframe(colframe))
					tkpack(get(fname),side = "top", fill = "x")
					
					tkpack(tklabel(get(fname), text = obj), side = "left")
					bname <- paste("b.", obj, sep = "")
					
					assign(bname,
							tkbutton(get(fname), text = "",
									background = slot(ngEnv$settings@color,obj),
									activebackground = slot(ngEnv$settings@color,obj),
									width = 2, height = 1,
									command = function(){
										col <- as.character(.Tcl(paste("tk_chooseColor -initialcolor", slot(ngEnv$settings@color,obj))))
										if(length(col) != 0){
											slot(ngEnv$settings@color,obj) <- col
											tkconfigure(get(paste("b.",obj,sep = "")), background = col,
													activebackground = col)
											if(obj == "background"){
												tkconfigure(ngEnv$canvas, background = col)	
											} else if (obj == "bulletActive") {
												tkitemconfigure(ngEnv$canvas, "bullet", activefill = col)	
											}else {
												.normalState(ngEnv)
											}
										}
									}))
					tkpack(get(bname), side = "right")
				})
		## buffer frame
		tkpack(tkframe(colframe, width = 200), anchor = 'w')
		
		
		## interaction Settings
		
		tkpack(tklabel(intframe, text = "Interaction", font = boldfont, padx = 38), side = "top", fill = "x")
		
		sapply(slotNames(ngEnv$settings@interaction), FUN = function(obj){		
					## frame
					fname <- paste("fi.",obj,sep = "")
					assign(fname, tkframe(intframe))
					tkpack(get(fname),side = "top", anchor = "w", fill = "x")
					
					fname_t <- paste("fit.",obj,sep = "")
					fname_b <- paste("fib.",obj,sep = "")
					assign(fname_t, tkframe(get(fname)))
					assign(fname_b, tkframe(get(fname)))
					tkpack(get(fname_t), get(fname_b),side = "top", fill = "x")
					
					## label
					tkpack(tklabel(get(fname_t), text = obj), side = "left")
					
					## create tcl variable
					vname <- paste("v.",obj, sep="")
					ngEnv$windowManager$settings[[vname]] <- tclVar(as.character(slot(ngEnv$settings@interaction,obj)))
					
					
					## and label
					lname <- paste("vl.",obj,sep = "")
					assign(lname, tklabel(get(fname_t),text=as.character(tclvalue(ngEnv$windowManager$settings[[vname]]))))
					tkconfigure(get(lname),textvariable=ngEnv$windowManager$settings[[vname]])
					tkpack(get(lname), side = "right")
					
					## slider
					sname <- paste("s.",obj,sep = "")
					#	
					
					if(obj == "animationTime") {
						from <- 0
						to <- 15
						resolution <- 0.1
					}else {
						from <- 1
						to <- 100
						resolution <- 1
					}
					assign(sname,tkscale(get(fname_b), from=from, to=to, length = 200,# width = 200,
									showvalue=F, variable=ngEnv$windowManager$settings[[vname]],
									resolution=resolution, orient="horizontal"))
					
					
					tkpack(get(sname),side = "top", fill = "x", expand = TRUE)
					tkconfigure(get(sname), command = function(...){
								.updateValues(ngEnv, obj)
							})
				})
		
		
		
		## Display Settings
		tkpack(tklabel(dispframe, text = "Display", font = boldfont, padx = 56), side = "top", fill = "x")
		
		
		sapply(slotNames(ngEnv$settings@display), FUN = function(obj){
					## frame
					fname <- paste("fi.",obj,sep = "")
					assign(fname, tkframe(dispframe))
					tkpack(get(fname),side = "top", anchor = "w", fill = "x")
					
					fname_t <- paste("fit.",obj,sep = "")
					fname_b <- paste("fib.",obj,sep = "")
					assign(fname_t, tkframe(get(fname)))
					assign(fname_b, tkframe(get(fname)))
					tkpack(get(fname_t), get(fname_b),side = "top", fill = "x")
					
					## label
					tkpack(tklabel(get(fname_t), text = obj), side = "left")
					
					## create tcl variable
					vname <- paste("v.",obj, sep="")
					ngEnv$windowManager$settings[[vname]] <- tclVar(as.character(slot(ngEnv$settings@display,obj)))
					
					## and label
					lname <- paste("vl.",obj,sep = "")
					assign(lname, tklabel(get(fname_t),text=as.character(tclvalue(ngEnv$windowManager$settings[[vname]]))))
					tkconfigure(get(lname),textvariable=ngEnv$windowManager$settings[[vname]])
					tkpack(get(lname), side = "right")
					
					## slider
					sname <- paste("s.",obj,sep = "")
					#	
					assign(sname,tkscale(get(fname_b), from=1, to=100, length = 200,# width = 20,
									showvalue=F, variable=ngEnv$windowManager$settings[[vname]],
									resolution=1, orient="horizontal"))
					
					
					tkpack(get(sname),side = "top", fill = "x", expand = TRUE)
					tkconfigure(get(sname), command = function(...){
								.updateValues(ngEnv, obj)  
							})
				})
	}

}




## Close all navGraph related windows
.closenavGraph <- function(ngEnv){
	if(!ngEnv$isDestoying){
		ngEnv$isDestoying <- TRUE
		cat(paste("RnavGraph Session", ngEnv$ng_instance,"closed.\n"))
		.closePlots(ngEnv)
		
		if(!is.null(ngEnv$windowManager$settings$tt)) {
			tkdestroy(ngEnv$windowManager$settings$tt)
		}
		
		if(!is.null(ngEnv$windowManager$paths$tt)) {
			tkdestroy(ngEnv$windowManager$paths$tt)
		}
	}
}
