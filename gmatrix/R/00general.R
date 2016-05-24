# General Functions
# 
# Author: nmorris
###############################################################################

setClass("gmatrix",
		representation(
				ptr="ANY",
				nrow = "integer",
				ncol = "integer",
				rownames="ANY",
				colnames="ANY",
				type="integer",
				device = "integer"),
		prototype = list(
				ptr=NULL,
				nrow=NULL,
				ncol=NULL,
				rownames=NULL,
				colnames=NULL,
				type=100L,
				device=100L)
) 



setClass("gvector",
		representation(
				ptr="ANY",
				length = "integer",
				names = "ANY",
				type="integer",
				device = "integer"),
		prototype = list(
				ptr=NULL,
				length=NULL,
				names=NULL,
				type=100L,
				device=0L)
) 

setClassUnion("index", members =  c("numeric", "logical", "character","gvector"))
#setClassUnion("numeric_matrix", members =  c("numeric", "matrix"))

.silent=FALSE

.onLoad <- function(lib, pkg) {
#	require(methods)
	.C("initialize_globals")
	computCap=.Call("get_device_info","major")+.Call("get_device_info","minor")/10
	if(length(computCap)<1)
		stop("No gpu devices detected.")
	#defualtdevice=which(computCap==max(computCap))[1]-1
	deviceOrder = order(computCap)-1 
	for(defualtdevice in  deviceOrder) {
		tmp = tryCatch(setDevice(defualtdevice), error = function(e) return(e))
		if("error" %in% class(tmp)) {
			warning("Unable to use device ", defualtdevice, ".\n")
			if(defualtdevice==max(deviceOrder))
				stop("No devices could be accesed.")
		} else {
			break
		}
	}
}

.onUnload <- function(libname) {
	#check_started()
	mem=.Call("get_device_info","totalGlobalMem")
	for(i in 0:(length(mem)-1)){
		tmp = .C("stopCublas", as.logical(.silent))
		tmp = .C("free_dev_states", as.logical(.silent))
		tmp = .C("deviceReset")
	}

}

listDevices = function() {
#	#(int* curdevice, int *memory, int *current, int *total, int *.silent)
#	current=as.integer(currentDevice)
#	total=as.integer(1)
#	memory=as.integer(20)
#	tmp=.C("RlistDevices", memory, current, total, silent)
#	total=tmp[[3]]
#	return(list(motmp[[1]][1:total]))
	tot=length(.Call("get_device_info","totalGlobalMem"))
	device=rep("", tot)
	device[.Call("get_device")+1]="Current Device"
	return(
			data.frame(
					row.names = make.names(.Call("get_device_info","name"), unique = TRUE),
					deviceNo = 0:(tot-1),
				#	globalMem=paste(round(-.Call("get_device_info","totalGlobalMem")/2^20), "Mb"),
					computeCape=.Call("get_device_info","major")+.Call("get_device_info","minor")/10,
					warpSize=.Call("get_device_info","warpSize"),
					clockRate=paste(.Call("get_device_info","clockRate")/1000, "MHz"),
					mpCount=.Call("get_device_info","multiProcessorCount"),
					currentDevice=device
			)
	)
}



setDevice = function(device,force=FALSE,silent=FALSE,...) {
#	WARNING("CHANGING GPU. CONSIDER DELETEING ANY VARIABLE CREATED ON A DIFFERENT GPU.
#					REFERENCING SUCH VARIABLES WILL HAVE UNDESIRED CONSEQUENCES.")
#	tmp=.C("stopCublas")
#	tmp=.C("free_dev_states")
	if(!is.numeric(device))
		stop("device must be numeric")
	if(length(device)!=1)
		stop("length(device) must be 1")
	if(!is.logical(silent))
		stop("silent must be numeric")
	if(length(silent)!=1)
		stop("length(silent) must be 1")
	
	tmp = .C("setDevice",as.integer(device), as.logical(silent))
#	tmp = .C("setFlagSpin")
    myset=integer(1)
	tmp=.C("startCublas", as.logical(silent), myset)
	setTuningPameters(force=force, ...)
	if(tmp[[2]]==1L) {#not sure why but sometimes cublas give error first time round... jerry riggin' it
        jerryrig=function() {
            x= g(matrix(0, 2, 2))
            y=  g(matrix(c(1, 1)))
            tmp=.Call("matrix_multiply", x, y, FALSE, FALSE, x@type)
        }
        tmp=tryCatch(jerryrig(), error=function(e) return("ERR"))
	}
}


getDevice = function() {
	return(.Call("get_device"))
}

checkDevice = function(x) {
	d=.Call("get_device")
	lapply(x, function(y) {
				#browser()
				if(y!=d)
					stop("You can only operate on gpu objects stored on the current gpu. Either change the current gpu device using setDevice(.) or transfer the object.")
				return(TRUE)
			})
	return(TRUE)
}

setTuningPameters=function(force=TRUE, threads_per_block=as.integer(2^8),
		total_states=as.integer(32*14*16), state=unclass(Sys.time())) {
	tmp=.C("set_threads_per_block", as.integer(threads_per_block))
	#SEXP setup_curand(SEXP in_total_states, SEXP in_seed, SEXP in_silent, SEXP in_force)
	tmp=.Call("setup_curand", as.integer(total_states), as.integer(state), as.logical(.silent), as.logical(force))
}
	
gset.seed=function(seed=unclass(Sys.time()), total_states=as.integer(32*14*16), silent=TRUE) {
	total_states=as.integer(total_states)
	if(total_states<as.integer(32*14*16))
		stop("Please make 'total_states' larger.")
	if(total_states>2^20)
		stop("Please make 'total_states' smaller.")
	state=as.integer(state)
	tmp=.Call("setup_curand", as.integer(total_states), as.integer(state), as.logical(silent), as.logical(TRUE))
	invisible(tmp)
}
#.resetDevice = function() {
#	warning("Any prviously created GPU variables will now be removed from the GPU.
#			Referencing such variables will have undesired consequences.")
#	temp = .C("stopCublas")#will result in fatal error if cublas not started
#	temp = .C("free_dev_states")
#	temp = .C("deviceReset")
#	.startDevice()
#}

ggc=function(silent=FALSE) {
	gc(verbose =FALSE)
	free=integer(1)
	tot=integer(1)
	tmp = .C("check_mem",free,tot,as.logical(silent))
}


.gpu_get =function(ptr, ln, tp) {
	if(tp==1L) { #float must be converted first
		ptr=.Call("gpu_convert", ptr,ln, 1L, 0L )
		tp=0L
	} 
	if(ln==0L) {
		if(tp==0 || tp==1)
			return(numeric(0))
		else if(tp==3)
			return(integer(0))
		else 
			return(logical(0))
	}
	
	ret=.Call("gpu_get", ptr, ln, tp)
	return(ret)
}



.type_num = function(type) {
	type=type[1]
	if(is.integer(type)) {
		if(type<4)
			return(type)
		else
			stop(paste("Invalid Type:", type),call. = FALSE)
	} else if(!is.character(type))
		stop("Invalid Type")
	else if(type=="d") 
		return(0L)
	else if(type=="s") 
		return(1L)
	else if(type=="i") 
		return(2L)
	else if(type=="l") 
		return(3L)
	else if(type=="double") 
		return(0L)
	else if(type=="single") 
		return(1L)
	else if(type=="integer") 
		return(2L)
	else if(type=="logical") 
		return(3L)
	else
		stop(paste("Invalid Type:", type),call. = FALSE)
}

.type_name= function(type) {
	if(type==0L) 
		return("double")
	else if(type==1L) 
		return("single")
	else if(type==2L) 
		return("integer")
	else if(type==3L) 
		return("logical")
	else
		stop("Invalid Type")
}


.Rclass_to_type=function(x) {
	if(is.integer(x))
		return(2L)
	else if(is.numeric(x))
		return(0L)
	else if(is.logical(x))
		return(3L)	
	else
		stop(paste("Unknown class:",class(x)))
}




.Rclass_to_type_int=function(x) {
	tmp=as.integer(x)
	if( as.integer(x)==x)
		return(2L)
	else if(is.numeric(x))
		return(0L)
	else if(is.logical(x))
		return(2L)	
	else
		stop(paste("Unknown class:",class(x)),call. = FALSE)
}


.convert_to_appropriate_class = function(x,type) {
	if(type==0L) {
		x=as.double(x)
	} else if(type==1L){
		x=as.double(x)
	} else if(type==2L){
		x=as.integer(x)
	} else if(type==3L){
		x=as.logical(x)
	} else
		stop("Invalid type")
	return(x)
}



#.check_input = function(IN1, IN2, IN3, allowed="DSIL") {#allowed is the types allowes D=doulbe, S=single, I=integer, L=logicla
#	if(is.missing(IN2))
#		stop("IN1 is missing")
#	if(allowed=="SFIL") {
#		if(is.missing(IN2))
#			return(TRUE)
#		else {
#			if(IN1@type!=IN2@type) {
#				stop("")
#			}
#		}
#		
#	}
#	
#}

convertType= function(x, to, dup=TRUE) {
	if(!(class(x) %in% c("gmatrix", "gvector")))
		stop("x is not a gpu object, so it's type cannot be converted")
	#cat(x@device, x@type, class(x),"\n")
	#print(x)
	checkDevice(x@device)
	totype=.type_num(to)
	fromtype=x@type
	if(fromtype!=totype) {
		if(totype==3L) {
			return(x!=0L)
		} else if(totype==2L && fromtype==3L){ #logical and int both stored as integers
			x@type=2L
			return(x)
		} else {
			ret=x
			ret@ptr=.Call("gpu_convert", x@ptr, length(x), fromtype, totype )
			ret@type=totype
			return(ret)
		}
	} else {
		if(dup)
			return(gdup(x))
		else
			return(x)
	}
}




g = function(x, type=NULL, dup=TRUE) {
	if(is.matrix(x))
		return(as.gmatrix(x,type=type,dup=dup))
	if(is.vector(x))
		return(as.gvector(x,type=type,dup=dup))
	stop("Input to 'g()' is allready a gpu object, or cannot be converted to a gpu object.")
}

h = function(x) {
	if(class(x)=="gmatrix")
		return(as.matrix(x))
	if(class(x)=="gvector") {
		ret=(as.vector(x))
		#if(!gnamestrip && length(names(x))>0)
		names(ret)=names(x)
		return(ret)
	}
	stop("Input to 'h' is not a gpu object.")
}





gdup=function(x, dev=getDevice()) {
	if(!(class(x) %in% c("gmatrix","gvector")))
		stop("not a gpu object.")
	if(x@device==dev) {
		oldDevice=getDevice()
		if(x@device!=oldDevice) {
			setDevice(x@device, silent=TRUE)
		}
		#checkDevice(x@device)
		mylength=length(x)
		ret=x
		ret@ptr=.Call("gpu_duplicate", x@ptr, mylength, x@type)
		if(x@device!=oldDevice) {
			setDevice(oldDevice, silent=TRUE)
		}
		return(ret)
	} else {
		ret=x
		device(ret)=dev
		return(ret)
	}
}

gnamestrip=function(x,dup=TRUE) {
	if(class(x) =="gvector")
		names(x)=NULL
	else if(class(x)=="gmatrix"){
		rownames(x)=NULL
		colnames(x)=NULL
	} else
		stop("not a gpu object.")
	if(dup)
		x=gdup(x)
	return(x)
} 
