setClass(
		Class = "NG_Viztk2d",
		representation = representation(
				viz_name = "character",
				toplevel = "OptionalTkwin",
				glyphVarOrder = "OptionalCharNumNULL",
				imgIds = "OptionalCharNULL",
				scaled = "logical",
				geometry = "character"
		),
		contains = "NG_Visualization2d"
)



ng_2d <- function(data, graph, images = NULL, glyphs = NULL) {
	
	if(!is(data,"NG_data")) {
		stop("[ng_2d] data argument is not an NG_data object.")
	}
	if(!is(graph,"NG_graph")) {
		stop("[ng_2d] graph argument is not an NG_graph object.")
	}
	if(!is.null(images)) {
		if(!is(images,"NG_image")) {
			stop("[ng_2d] image argument is not an NG_image object.")
		}
	}
	if(!(all(glyphs %in% shortnames(data)) || all(glyphs %in% names(data)))) {
		stop("[ng_2d] glyph argument contains names which are not present in the data.")
	}
	
	## check if variable- and node names of graph and data match
	varNames <-	vizVarNames(graph,data)
	
	if(is.null(images)) {
		ids <- NULL	
	}else{
		if(is(images,"NG_image")) {
			ids <- images@ids			
		}else {
			stop('[ng_2d] image argument not of class NG_image')
		}
	}
	
	## Glyph Variable order
	if(!is.null(glyphs)) {
		if(all(glyphs %in% names(data@data))) {	
			glyphVarOrder <- match(glyphs, names(data@data))
		} else if(all(glyphs %in% data@shortnames)) {
			glyphVarOrder <- match(glyphs, data@shortnames)
		} else {
			stop('[ng_2d] argument glyph does not match with the data.')
		}
	}else {
		glyphVarOrder <- NULL
	}
	
	
	return( new("NG_Viztk2d",
					graph = graph@name,
					data = data@name,
					from = nodes(graph@graph)[1],
					to = "",
					varList = varNames,
					mat = matrix(rep(0,length(varNames)*2),ncol = 2),
					transitionKind = 0,
					viz_name = "",
					glyphVarOrder = glyphVarOrder,
					imgIds = ids,
					scaled = TRUE
			))
}



## Initialize Plots
setMethod(
		f = "initializeViz",
		signature = "NG_Viztk2d",
		definition = function(viz,ngEnv){
			
			viz@toplevel <- tktoplevel()
			
			if(length(viz@geometry)>0) {
				tkwm.geometry(viz@toplevel,viz@geometry)
			}else {
				if(exists('tk2dDefaultSize',envir = ngEnv)){
					tkwm.geometry(viz@toplevel,ngEnv$tk2dDefaultSize)
				}
			}
			
			viz <- initRotation(viz,ngEnv)
			tcl('set', paste("ng_data(\"",ngEnv$ng_instance,".",viz@data,".xcoord\")", sep = ''), ng_2d_xcoord(viz,ngEnv))
			tcl('set', paste("ng_data(\"",ngEnv$ng_instance,".",viz@data,".ycoord\")", sep = ''), ng_2d_ycoord(viz,ngEnv))
			tcl('tk_2d_display', viz@toplevel, ngEnv$ng_instance, ngEnv$ng_LinkedInstance, viz@data, viz@viz_name, !is.null(viz@imgIds), !is.null(viz@glyphVarOrder), paste('Session ', ngEnv$ng_instance,', data: ',viz@data,", graph: ",viz@graph,sep=""))
			
			
			return(viz)
		})


setMethod(
		f = "updateViz",
		signature = "NG_Viztk2d",
		definition = function(viz,ngEnv){
			
			## TODO: check changed plots variable
			#if(ngEnv$changedPlots){				
			viz <- ng_2dRotationMatrix(viz,ngEnv)			
			tcl('set', paste("ng_data(\"",ngEnv$ng_instance,".",viz@data,".xcoord\")", sep = ''), ng_2d_xcoord(viz,ngEnv))
			tcl('set', paste("ng_data(\"",ngEnv$ng_instance,".",viz@data,".ycoord\")", sep = ''), ng_2d_ycoord(viz,ngEnv))
			#	ngEnv$changedPlots <- FALSE
			#}
			
			tcl('update_displays',viz@toplevel, ngEnv$ng_instance, viz@data, viz@viz_name)
			
			return(viz)
		})


setMethod(
		f = "closeViz",
		signature = "NG_Viztk2d",
		definition = function(viz,ngEnv){
			## first save window size
			viz@geometry <- .tcl2str(tkwm.geometry(viz@toplevel))
			ngEnv$tk2dDefaultSize <- viz@geometry 
			
			tkdestroy(viz@toplevel)
			
			t.array <- paste("ng_windowManager(\"",ngEnv$ng_LinkedInstance,".",viz@data,".ttID\")", sep = '')
			i <- .tcl2num(.Tcl(paste('lsearch -exact $',t.array, ' ', viz@toplevel, sep = '')))
			
			if(i != -1) {
				.Tcl(paste('set ', t.array, ' [lreplace $',t.array,' ',i,' ',i,']', sep=''))
			}else {
				warning("[closeViz] could not find toplevel window!")
			}
			viz@toplevel <- NULL
			
			return(viz)
		})
