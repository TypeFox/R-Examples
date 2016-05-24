setClass(
		Class = "NG_Viz2DAxis",
		representation = representation(
				devName = "integer",
				FUN = "character",
				devType = "OptionalCharNULL",
				scaled = "logical"
		),
		contains = "NG_Visualization2d"
)

## object creator function
ng_2d_myplot <- function(data,graph,fnName, device = "base", scaled=TRUE){
	
	if(is(data,"NG_data") == FALSE){
		stop("data is no NG_data object.\n")
	}
	if(is(graph,"NG_graph") == FALSE){
		stop("graph is no NG_graph object.\n")
	}
	
	if(!is.null(device)) {
	 if(!(device %in% c("base", "grid", "ggplot2", "lattice", "rgl"))) {
	 	stop("argument device has to be either 'base', 'grid', 'ggplot2', 'lattice' or 'rgl'.")
	 }
	}
	
	## match Variable names with Graph nodes
	varNames <-	vizVarNames(graph,data)
	
	return( new("NG_Viz2DAxis",
					graph = graph@name,
					data = data@name,
					from = nodes(graph@graph)[1],
					to = "",
					varList = varNames,
					devType = device,
					scaled = scaled,
					mat = matrix(rep(0,length(varNames)*2),ncol = 2),
					FUN = fnName,
					transitionKind = 0)
	)
}



setMethod(
		f = "initializeViz",
		signature = "NG_Viz2DAxis",
		definition = function(viz,ngEnv){
			
			## open new device
		  	if(!is.null(viz@devType)) {
			if(any(viz@devType %in% c("base", "grid", "ggplot2", "lattice"))){
				dev.new(title= ngEnv$graph@name)
				viz@devName <- dev.cur()			
			}else if(viz@devType == "rgl"){
				rgl.open()
				viz@devName <- rgl.cur()							
			}
			}
			
			## initialize rotation matrix	
			viz <- initRotation(viz,ngEnv)
			
			init <- paste(viz@FUN,'.init',sep='')
			fun <- viz@FUN
			
			if(exists(init,envir = .GlobalEnv)) {
				if(is.function(get(init,envir=.GlobalEnv))) {
					fun <- init
				}
			} 
			
			ii <-  c("x","y","group","labels","order","from","to","percentage","data") %in% names(formals(get(fun)))
	
			## plot regular update function
			do.call(fun, list(x = ng_2d_xcoord(viz,ngEnv),
							y =  ng_2d_ycoord(viz,ngEnv),
							group = ngEnv$dataList[[viz@data]]@group,
							labels = ngEnv$dataList[[viz@data]]@labels,
							order =  ng_2d_dist(viz,ngEnv),
							from = ngEnv$bulletState$from,
							to = ngEnv$bulletState$to,
							percentage = ngEnv$bulletState$percentage,
							data = viz@data)[ii],
					envir = .GlobalEnv)
			
			return(viz)
		}
)




setMethod(
		f = "updateViz",
		signature = "NG_Viz2DAxis",
		definition = function(viz,ngEnv){
			
			
			
			
			if(!is.null(viz@devType)) {
			if(viz@devType == "rgl"){
				rgl.set(viz@devName)
			} else {
				dev.set(viz@devName)			
			}
			}
			viz <- ng_2dRotationMatrix(viz,ngEnv)		

			ii <-  c("x","y","group","labels","order","from","to","percentage","data") %in% names(formals(get(viz@FUN)))
			
			## plot regular update function
			do.call(viz@FUN, list(x = ng_2d_xcoord(viz,ngEnv),
							y =  ng_2d_ycoord(viz,ngEnv),
							group = ngEnv$dataList[[viz@data]]@group,
							labels = ngEnv$dataList[[viz@data]]@labels,
							order =  ng_2d_dist(viz,ngEnv),
							from = ngEnv$bulletState$from,
							to = ngEnv$bulletState$to,
							percentage = ngEnv$bulletState$percentage,
							data = viz@data)[ii],
					envir = .GlobalEnv)
			
			
			return(viz)
		})




setMethod(
		f = "closeViz",
		signature = "NG_Viz2DAxis",
		definition = function(viz,ngEnv){
if(!is.null(viz@devType)) {
			if(viz@devType == "rgl"){
				rgl.set(viz@devName)
				rgl.close()		
			} else {
				dev.off(viz@devName)			
			}
}
			return(viz)
		})
