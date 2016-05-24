setClass(
		Class = "NavGraph_handler",
		representation = representation(
				env = "environment",
				graphs = "list",
				data = "list",
				viz = "list",
				settings = "NG_Settings",
				paths = "NG_path",
				activePath = "character",
				activePathGraph = "character",
				tk2dcolors = "list",
				dateCreated = "character",
				dateUpdated = "character"
		),
		validity = function(object){
			goThrough <- TRUE
			
			if(!all(sapply(object@graphs,FUN = function(x){is(x,"NG_graph")}))) {
				cat("[validity RnavGraph handler] elements in graphs are not from the NG_graph class\n")
				goThrough <- FALSE	
			}
			if(!all(sapply(object@data,FUN = function(x){is(x,"NG_data")}))) {
				cat("[validity RnavGraph handler] elements in graphs are not from the NG_data class\n")
				goThrough <- FALSE	
			}
			if(!all(sapply(object@viz,FUN = function(x){is(x,"NG_Visualization")}))) {
				cat("[validity RnavGraph handler] elements in graphs are not from the NG_Visualization class\n")
				goThrough <- FALSE	
			}
			return(goThrough)
		}
)


## show Method
setMethod("show","NavGraph_handler",
		function(object){
			cat("RnavGraph handler:\n---\n")
			cat(paste("created       :",object@dateCreated,"\nlast updated  :",object@dateUpdated,"\n"))
			cat("---\n")
			cat(paste("graphs        :",paste(sapply(object@graphs,function(x)x@name),collapse=", "),"\n"))
			cat(paste("data          :",paste(sapply(object@data,function(x)x@name),collapse=", "),"\n"))
		}
)


## Walk a path
ng_walk <- function(nghandler,path){
	if(is(nghandler,"NavGraph_handler") == FALSE){
		stop("[ng_walk] nghandler is not of class NavGraph_handler")
	}
	if(length(path)>1){
		path <- paste(path, collapse = " ")
	}
	tclvalue(nghandler@env$activePath) <- path
	tclvalue(nghandler@env$activePathGraph) <- get(nghandler@env$cGraph, envir=nghandler@env)@name
	.walkPath(nghandler@env, path)
}

ng_update <- function(nghandler){
	if(is(nghandler,"NavGraph_handler")){
		
		## does the environment still exist ?
		if(TRUE){		
			nghandler@dateUpdated <- date()		
			
			## create a Graph list
			n <- length(nghandler@graphs)
			nghandler@graphs <- eval(parse(text = paste("list(",paste(paste('graph',1:n,sep = ""),collapse = ","),")")), envir = nghandler@env)
			
			nghandler@settings <- nghandler@env$settings
			nghandler@paths <- nghandler@env$paths
			nghandler@activePath <- tclvalue(nghandler@env$activePath)
			nghandler@activePathGraph <- tclvalue(nghandler@env$activePathGraph)
			
			
			## update the group slot in the ng_data objects
			for(i in 1:length(nghandler@data)) {
				dataName <- nghandler@data[[i]]@name
				classes <- sapply(nghandler@viz, FUN=function(viz){
							if(viz@data == dataName) {
								return(class(viz))
							} else {
								return("NONE")
							}
						})
				
				
				
				if ("NG_Viztk2d" %in% classes ){
					## tk2d
					
					## find data name and instance
					size <- .tcl2num(paste('ng_data("',nghandler@env$ng_LinkedInstance,".",dataName,'.size")',sep = ''))
					col <- .tcl2str(paste('ng_data("',nghandler@env$ng_LinkedInstance,".",dataName,'.color")',sep = ''))
					
					group <- paste('c',col,';s',size,sep='')
				}	

				nghandler@data[[i]]@group <- group
				
				cols <- .tcl2str(paste('ng_data("',nghandler@env$ng_LinkedInstance,".",dataName,'.brush_colors")',sep = ''))
				## seems to be a tcl error lset a 0 #1 makes a list within the list
				
				cols <- sapply(cols, FUN=function(str){
							t1 <- sub("^\\{","",str)
							return(sub("\\}$","",t1))
						})
				nghandler@tk2dcolors[[dataName]] <- cols
				nghandler@tk2dcolors[[paste('bg:',dataName,sep='')]] <-  .tcl2str(paste('ng_data("',nghandler@env$ng_LinkedInstance,".",dataName,'.bg")',sep = ''))
				nghandler@tk2dcolors[[paste('sel:',dataName,sep='')]] <- .tcl2str(paste('ng_data("',nghandler@env$ng_LinkedInstance,".",dataName,'.brush_color")',sep = ''))
			}
			
			return(nghandler)
		}else{		
			stop("[ng_update]: environment does not exist anymore")
		}
	}else{
		stop("[ng_update]: argument is not a navGraph handler")
	}
}


setMethod(f = "ng_get",
		signature = "NavGraph_handler",
		definition = function(obj, what=NULL, ...){
			possibleOptions <- c("graphs","paths","data","viz")
			
			if(is.null(what)){
				cat("possible options are: ")
				cat(paste(possibleOptions, collapse = ", "))
				cat("\n")
			}else{
				if(any(is.na(match(what,possibleOptions)))){
					stop(paste("[ng_get] object",what,"is not defined."))
				}
				
				if(what == "graphs"){
					if(length(obj@graphs) == 1) {
						return(obj@graphs[[1]])
					}else {
						return(obj@graphs)
					}					
				}else if(what == "paths"){
					return(obj@paths)
				}else if(what == "data"){
					if(length(obj@data) == 1) {
						return(obj@data[[1]])
					}else {
						return(obj@data)
					}		
				} else if(what == "viz"){
					if(length(obj@viz) == 1) {
						return(obj@viz[[1]])
					}else {
						return(obj@viz)
					}
				}
			}
		})


setReplaceMethod(
		f = "ng_set","NavGraph_handler",
		function(object,what,value){
			possibleOptions <- c("graphs","paths","data","viz")
			
			tmp <- object
			
			if(!(what %in% possibleOptions)) {
				stop(paste("Replace what? Possible options are: ", paste(possibleOptions, collapse = ", "),'\n'))				
			}	
			
			slot(object, what, check = TRUE) <- value 
			
			if(validObject(object)){
				return(object)
			}else {
				cat("[ng_set] assignment is wrong\n")
				return(tmp)
			} 
		})


## ng_get_color
ng_get_color <- function(obj, dataName) {

  if(!is(obj, "NavGraph_handler")) {
    stop('[ng_get_color]: first arguments needs to be a NavGraph_handler object')
  }

  dnames <- sapply(obj@data, function(x)x@name)
  
  isMissing <- FALSE
  if(missing(dataName)) {
    if(length(obj@data) > 1) {
      cat("Specify data name. Choose from:\n")
      cat(paste('   ', paste(sapply(obj@data, function(x)x@name), collapse = ', '), '\n'))
      isMissing <- TRUE
    } else {
      dataName <- dnames[1]
    }
  } else {
    ## does dataName exist
    if(!(dataName %in% dnames)) {
      stop(paste('[ng_get_color]: data name "',dataName,'" does not exist in your NavGraph_handler', sep=''))
    }
  }

  if(!isMissing) {
    return(as.character(tcl('set',paste('ng_data("',obj@env$ng_LinkedInstance,'.',dataName,'.','color','")',sep = ''))))
  }
  
}


setReplaceMethod(
                 f = "ng_set_color",
                 signature = signature(obj = "NavGraph_handler"),
                 function(obj,dataName,value){

                   dnames <- sapply(obj@data, function(x)x@name)

                   isMissing <- FALSE
                   if(missing(dataName)) {
                     if(length(obj@data) > 1) {
                       cat("Specify data name. Choose from:\n")
                       cat(paste('   ', paste(sapply(obj@data, function(x)x@name), collapse = ', '), '\n'))
                       isMissing <- TRUE
                     } else {
                       dataName <- dnames[1]
                     }
                   } else {
                     ## does dataName exist
                     if(!(dataName %in% dnames)) {
                       stop(paste('[ng_get_color]: data name "',dataName,'" does not exist in your NavGraph_handler', sep=''))
                     }
                   }

                   if(!isMissing){
                     ## where is the data stored
                     ind <- match(dataName,dnames)[1]
                     n <- dim(obj@data[[ind]]@data)[1]
                     
                     if(length(value) == 1) {
                       tcl('set',paste('ng_data("',obj@env$ng_LinkedInstance,'.',dataName,'.','color','")',sep = ''), rep(value, n))
                     } else if(length(value) == n) {
                       tcl('set',paste('ng_data("',obj@env$ng_LinkedInstance,'.',dataName,'.','color','")',sep = ''), value)
                     } else {
                       stop('[ng_set_color]: length of specified color vector does not match with data dimension')
                     }
                     tcl('refresh_linked', obj@env$ng_LinkedInstance, dataName)
                   }
                   return(obj)
                 })



## ng_get_size
ng_get_size <- function(obj, dataName) {

  if(!is(obj, "NavGraph_handler")) {
    stop('[ng_get_color]: first arguments needs to be a NavGraph_handler object')
  }

  dnames <- sapply(obj@data, function(x)x@name)
  
  isMissing <- FALSE
  if(missing(dataName)) {
    if(length(obj@data) > 1) {
      cat("Specify data name. Choose from:\n")
      cat(paste('   ', paste(sapply(obj@data, function(x)x@name), collapse = ', '), '\n'))
      isMissing <- TRUE
    } else {
      dataName <- dnames[1]
    }
  } else {
    ## does dataName exist
    if(!(dataName %in% dnames)) {
      stop(paste('[ng_get_color]: data name "',dataName,'" does not exist in your NavGraph_handler', sep=''))
    }
  }

  if(!isMissing) {
    return(as.numeric(tcl('set',paste('ng_data("',obj@env$ng_LinkedInstance,'.',dataName,'.','size','")',sep = ''))))
  }
  
}


setReplaceMethod(
                 f = "ng_set_size",
                 signature = signature(obj = "NavGraph_handler"),
                 function(obj,dataName,value){
                   dnames <- sapply(obj@data, function(x)x@name)
                   
                   isMissing <- FALSE
                   if(missing(dataName)) {
                     if(length(obj@data) > 1) {
                       cat("Specify data name. Choose from:\n")
                       cat(paste('   ', paste(sapply(obj@data, function(x)x@name), collapse = ', '), '\n'))
                       isMissing <- TRUE
                     } else {
                       dataName <- dnames[1]
                     }
                   } else {
                     ## does dataName exist
                     if(!(dataName %in% dnames)) {
                       stop(paste('[ng_get_color]: data name "',dataName,'" does not exist in your NavGraph_handler', sep=''))
                     }
                   }

                   if(!isMissing){
                     ## where is the data stored
                     ind <- match(dataName,dnames)[1]
                     n <- dim(obj@data[[ind]]@data)[1]
                     
                     if(length(value) == 1) {
                       tcl('set',paste('ng_data("',obj@env$ng_LinkedInstance,'.',dataName,'.','size','")',sep = ''), rep(value, n))
                     } else if(length(value) == n) {
                       tcl('set',paste('ng_data("',obj@env$ng_LinkedInstance,'.',dataName,'.','size','")',sep = ''), value)
                     } else {
                       stop('[ng_set_color]: length of specified color vector does not match with data dimension')
                     }
                     tcl('refresh_linked', obj@env$ng_LinkedInstance, dataName)
                   }
                   return(obj)
                 })

