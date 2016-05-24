# S3 methods for datamap objects
`$.datamap` <- function(x,name) x[['.intern']][[name]]

`$<-.datamap` <- function(x,name,value) { x[['.intern']][[name]] <- value; x }

`print.datamap` <- function(x, ...){
	instIDX <- unlist(lapply(objects(x),bindingIsActive,x))
	print(
		list(
			installed=objects(x)[instIDX],
			extra=objects(x)[!instIDX],
			internal=objects(base::get('.intern',x),all.names=TRUE)
		)
	)
}

# Global list of mappers
Mappers <- new.env(hash=TRUE)

x <- new.env(hash=TRUE)

print.dataMappers <- function(x,...){
	type <- character()
	init <- get <- assign <- finalize <- logical()
	for(mapper in objects(Mappers,all.names=TRUE)){
		type <- append(type,mapper)
		init <- append(init,exists('init',envir=Mappers[[mapper]],inherits=FALSE))
		get <- append(get,exists('get',envir=Mappers[[mapper]],inherits=FALSE))
		assign <- append(assign,exists('assign',envir=Mappers[[mapper]],inherits=FALSE))
		finalize <- append(finalize,exists('finalize',envir=Mappers[[mapper]],inherits=FALSE))
	}
	print(data.frame(type=type,init=init,get=init,assign=assign, finalize=finalize))
}

# Inspecting mappers
print.dataMapper <- function(x,...){
	print(as.list(x))
}

newMapper <- function( type=NULL, init=NULL, get=NULL, assign=NULL, finalize=NULL, ...){

	if (is.null(init))
		stop("Need an init function at the least.")

	mapper <- new.env(hash=TRUE)

	base::assign('type',type,mapper);

	environment(init) <- mapper; 
	base::assign('init',init,mapper)

	if (!is.null(get)){
		environment(get) <- mapper
		base::assign('get',get,mapper)
	}
	if (!is.null(assign)){
		environment(assign) <- mapper
		base::assign('assign',assign,mapper)
	}
	if (!is.null(finalize)){
		environment(finalize) <- mapper
		base::assign('finalize',finalize,mapper)
	}

	elems <- list(...)
	if (length(elems)>0){
		lapply(
			names(elems),
			function(sym) {
				if (sym==''){
					warning('Arguments in ... must be named')
					return(NULL)
				}
				obj <- elems[[sym]]
				if (typeof(obj) == "closure")
					environment(obj) <- mapper
				base::assign(sym,obj,mapper)
			}
		)
	}

	base::assign(type,mapper,Mappers)
	class(Mappers[[type]]) <- c("dataMapper",class(Mappers[[type]]))
	invisible(Mappers[[type]])
}

# Make R CMD check shut up. Is there a better way?
type <- NULL;
binderGetOnly <- function(val,sym){
	if (!missing(val))
		warning(paste("Assignments to",sym,"are not possible with maps of type",type,"."))
	return(get(sym))
}

binderGetAssign <- function(val,sym){
	if (missing(val))
		return(get(sym))
	else
		return(assign(sym,val))
}

binderAssignOnly <- function(val,sym){
	if (missing(val))
		warning(paste("Cannot read variable",sym,"with maps of type",type,"."))
	else
		return(assign(sym,val))
}

newMap <- function(type=character(),...) {

	# Look up mapper type
	if (!exists(type,Mappers))
		stop(paste("Type",type,"not found!"))

	mapper <- base::get(type,Mappers)

	# The datamap object, an outward (map) facing environment
	# with a inward facing (intern) environment.
	intern <- new.env(hash=TRUE,globalenv())
	class(intern) <- c("datamapIntern",class(intern))
	map <- new.env(hash=TRUE,globalenv())
	class(map) <- c("datamap",class(map))
	base::assign('.intern',intern,map)
	base::assign('.map',map,intern)

	# Copy mapper into new map, fixing up environment for functions
	lapply(
		objects(mapper,all.names=TRUE),
		function(x) {
			obj <- base::get(x,mapper)
			if (typeof(obj) == "closure")
				environment(obj) <- intern
			base::assign(x,obj,intern)
		}
	)

	# init returns TRUE on success
	# Any other return type results in an error
	ret <- try(map$init(map,...))
	if (!is.logical(ret) || ret[1]!=TRUE)
		stop("init() must return TRUE to instantiate a datamap object.")

	if (!is.null(map$finalize))
		reg.finalizer(map,map$finalize)
	return(map)
}

install <- function(symbols,map){
	if (!inherits(map,"datamap"))
		stop("object not a datamap")

	if (!is.null(map$assign)){
		if (!is.null(map$get))
			binder <- binderGetAssign
		else
			binder <- binderAssignOnly
	} else if (!is.null(map$get)){
		binder <- binderGetOnly
	} else {
		stop("No get or assign methods for this map")
	}

	intern <- base::get('.intern',map)
	installed <- unlist(lapply(symbols,
		function(name){
			bind <- binder
			environment(bind) <- intern
			formals(bind)$sym=name
			makeActiveBinding(name,bind,map)
			name
		}
	))
	invisible(installed)
}

uninstall <- function(symbols,map){
	if (!inherits(map,"datamap"))
		stop("object not a datamap")
	rm(list=symbols,envir=map)
}

mapAttach <- function(map,pos=2,name=NULL,warn.conflicts=TRUE){
	# Create a UserDefinedDatabase object and attach to search path
	if (!inherits(map,"datamap"))
		stop("object not a datamap")
	if (missing(name))
		name <- paste('datamap',map$type,sep=':')
	attach(.Call('CreateUserDB',map,PACKAGE='datamap'),pos,name,warn.conflicts)
}
