# gvector definitions
# 
# Author: nmorris
###############################################################################

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

setMethod("initialize",
		"gvector",
		function(.Object,...) {
			#cat("test\n")
		#	browser()
			.Object <- callNextMethod()
			.Object@device=.Call("get_device")
			return(.Object)
		})


gvector = function(length, type="d") {
	return(g.rep(0, times=length, type=type))
}

g.rep =function(x, times=1L, each=1L, type=NULL) {
	if(is.null(type))
		tryCatch(typeno<- .Rclass_to_type(x),
				error=function(i) stop("Class of data cannot be converted to gpu type.", call. = FALSE))
	else
		typeno=.type_num(type)
	times=as.integer(times[1])
	each=as.integer(each[1])
	ret = new("gvector",length=as.integer(length(x)*times*each), type=typeno)
	if(length(ret)==0)
		return(ret)
	if(length(x)==1L) {
		tryCatch(x <- .convert_to_appropriate_class(x,typeno), error=function(e) stop("Invalid value for x, or invalid type.", call. = FALSE))
		ret@ptr = .Call("gpu_rep_1",x,length(ret) ,typeno)#gpu_rep_1(SEXP in_val, SEXP in_N)
	} else {
		
		if(class(x)!="gvector")
			x=as.gvector(x)
		x=convertType(x,typeno,dup=FALSE)
		x@names=NULL
		
		ret@ptr = .Call("gpu_rep_m", x@ptr, as.integer(length(x)), times, each , typeno)#gpu_rep_m(SEXP in_A,SEXP in_n, SEXP in_N, SEXP in_times_each);
	}
	
	return(ret)
}

.gcolon = function(a,b, type=NULL) {
	#SEXP gpu_seq( SEXP n_in, SEXP init_in, SEXP step_in, SEXP in_type  )
	a=a[1];b=b[1]
	if(is.null(type)) {
		typeno=.Rclass_to_type_int(a)
	} else {
		typeno=.type_num(type)
	}
	if(typeno==0L || typeno==1L) {
		a=as.double(a)
		b=as.double(b)
		n=as.integer(floor(abs(b-a)+1L))
		step=as.double(sign(b-a))
	} else if(typeno==2) {
		a=as.integer(a)
		b=as.integer(b)
		n=as.integer(abs(a-b)+1L)
		step=as.integer(sign(b-a))
	} else
		stop("Invalid type.")
	return(.Call("gpu_seq", n, a, step, typeno))
	
	
	#gpu_seq( SEXP n_in, SEXP init_in, SEXP step_in, SEXP in_type  )
}

setGeneric("%to%",
		function(from,to)
			standardGeneric("%to%")
)

setMethod("%to%", c("numeric","numeric"),
		function(from,to) {
			return(.gcolon(from,to))
		})

.seq = function(n,from,by, type=NULL) {
	
	if(is.null(type)) {
		typeno=min(.Rclass_to_type_int(from),.Rclass_to_type_int(by))
	} else {
		typeno=.type_num(type)
	}
	if(typeno==0L || typeno==1L) {
		from=as.double(from)
		by=as.double(by)
	} else if(typeno==2) {
		from=as.integer(from)
		by=as.integer(by)
	} else
		stop("Invalid type.",call. = FALSE)
	return(.Call("gpu_seq", as.integer(n), from, by, typeno))
	
	
	#gpu_seq( SEXP n_in, SEXP init_in, SEXP step_in, SEXP in_type  )
}

#modified from the R seq function
gseq = function (from = 1, to = 1, by = ((to - from)/(length.out - 1)), 
		length.out = NULL, along.with = NULL, type=NULL) 
{
	if ((One <- nargs() == 1L) && !missing(from)) {
		lf <- length(from)
		return(if (mode(from) == "numeric" && lf == 1L) .gcolon(1L,from,type) else if (lf) .gcolon(1L,lf,type) else gvector(0L,type="i"))
	}
	if (!missing(along.with)) {
		length.out <- length(along.with)
		if (One) 
			return(if (length.out) .gcolon(1L,length.out,type) else gvector(0L,type="i"))
	}
	else if (!missing(length.out)) 
		length.out <- ceiling(length.out)
	if (is.null(length.out)) 
		if (missing(by)) 
			.gcolon(from,to,type)
		else {
			del <- to - from
			if (del == 0 && to == 0) 
				return(as.gvector(to))
			n <- del/by
			if (!(length(n) && is.finite(n))) {
				if (length(by) && by == 0 && length(del) && del == 
						0) 
					return(as.gvector(from))
				stop("invalid (to - from)/by in seq(.)")
			}
			if (n < 0L) 
				stop("wrong sign in 'by' argument")
			if (n > .Machine$integer.max) 
				stop("'by' argument is much too small")
			dd <- abs(del)/max(abs(to), abs(from))
			if (dd < 100 * .Machine$double.eps) 
				return(as.gvector(from))
			if(abs(as.integer(n)-n)<100 * .Machine$double.eps){
				n=as.integer(round(n)+1L)
				tmp=.seq(n,from,by)
				tmp[n]=to
				return(tmp)
			} else {
				n=as.integer(n)
				return(.seq(n,from,by))
			}
			
		}
	else if (!is.finite(length.out) || length.out < 0L) 
		stop("length must be non-negative number")
	else if (length.out == 0L) 
		gvector(0L,type="i")
	else if (One) 
		.gcolon(1L,length.out)
	else if (missing(by)) {
		if (missing(to)) 
			to <- from + length.out - 1L
		if (missing(from)) 
			from <- to - length.out + 1L
		if (length.out > 2L) 
			.seq(length.out,from,by)
	}
	else if (missing(to)) 
		.seq(length.out,from,by)
		#from + (0L:(length.out - 1L)) * by
	else if (missing(from)) 
		.seq(length.out,to,-by)
		#to - ((length.out - 1L):0L) * by
	else stop("too many arguments")
}
		


setGeneric("as.gvector", useAsDefault=function(x, type=NULL, dup=TRUE) {
			#browser()
			tryCatch(fromtype<- .Rclass_to_type(x),
					error=function(i) stop("Class of data cannot be converted to gpu type.", call. = FALSE))
			
			if(is.null(type)) {
				typeno<-fromtype
			} else
				typeno=.type_num(type)
			
			ret = new("gvector",length=length(x), type=fromtype, names=names(x))
			if(length(ret)==0L)
				return(ret)
			
			ret@ptr = .Call("gpu_create", x, fromtype)
			
			if(typeno!=fromtype)
				ret=convertType(ret,typeno)
			
			if(length(names(x))!=0)
				if(length(tmp)==length(names(x)))
					ret@names=names(x)
			
			return(ret)
		})

setMethod("as.gvector", "gmatrix",
		function(x, type=x@type, dup=TRUE) {
			checkDevice(x@device)
			typeno=.type_num(type)
			mylength=as.integer(prod(dim(x)))
			ret = new("gvector", length=mylength, type=typeno)
			if(length(ret)==0L)
				return(ret)
			if(x@type==typeno) {
				if(dup)
					ret@ptr = .Call("gpu_duplicate", x@ptr, mylength, typeno)#duplicate x@ptr
				else
					ret@ptr = x@ptr
			} else {
				tmp=convertType(x,typeno,dup)
				ret@ptr=tmp@ptr
			}

			return(ret)
		})

setMethod("as.gvector", "gvector",
		function(x, type=x@type, dup=TRUE) {
			checkDevice(x@device)
			typeno=.type_num(type)
			if(x@type==typeno) {
				if(dup)
					ret = gdup(x)
				else
					ret = x
			} else {
				ret=convertType(x,typeno,dup)
			}
			return(ret)
		})


setMethod("as.numeric", "gvector",
		function(x) {
			checkDevice(x@device)
			ret=.gpu_get( x@ptr, x@length, x@type)
#			if(!namestrip && length(names(x))>0)
#				names(ret)=names(x)
			return(as.numeric(ret))
		})

setMethod("as.integer", "gvector",
		function(x) {
			checkDevice(x@device)
			ret=.gpu_get( x@ptr, x@length, x@type)
#			if(!namestrip && length(names(x))>0)
#				names(ret)=names(x)
			return(as.integer(ret))
		})

setMethod("as.logical", "gvector",
		function(x) {
			checkDevice(x@device)
			ret=.gpu_get( x@ptr, x@length, x@type)
#			if(length(names(x))>0)
#				names(ret)=names(x)
			return(as.logical(ret))
		})

#setGeneric("as.vector",
#		function(x, ...)
#			base::as.vector(x,...)
#)

as.vector.gvector =  function(x, mode=NULL) {
	checkDevice(x@device)
	ret=.gpu_get( x@ptr, x@length, x@type)
	if(!is.null(mode))
		mode(ret)=mode
	return(ret)
}


setMethod("as.vector", "gvector",as.vector.gvector)




setMethod("length", "gvector",
		function(x) {
			return(x@length)
		})

setMethod("names", "gvector",
		function(x) {
			return(x@names)
		})


setReplaceMethod("names", "gvector",
		function(x, value) {
			if(is.null(value))
				x@names <- NULL
			else if(length(x) != length(value))
				stop("'names' attribute must be the same length as the vector.")
			else
				x@names <- as.character(value)
			return(x)
		}
)
#setReplaceMethod("length", "gvector",
#		function(x, value) {
#			value=as.integer(value)
#			x@length <- value
#			return(x)
#		}
#)

#device = function(x) {
#	stop(paste("'device' not defined for class:",class(x)))
#}
setGeneric("device",
		function(x)
			standardGeneric("device")
)

setMethod("device", "gmatrix",
		function(x) {
			return(x@device)
		})
setMethod("device", "gvector",
		function(x) {
			return(x@device)
		})

setGeneric("device<-", function(x, value)
			standardGeneric("device<-"))

setReplaceMethod("device", "gmatrix",
		function(x, value) {
			curD=getDevice()
			value=as.integer(value)[1]
			if(x@device!=value) {
				if(x@device!=curD)
					setDevice(x@device, silent=TRUE)
				x=h(x)
				setDevice(value, silent=TRUE)
				x=g(x)
				if(getDevice()!=curD)
					setDevice(curD, silent=TRUE)
			}
			return(x)
		})

setReplaceMethod("device", "gvector",
		function(x, value) {
			curD=getDevice()
			value=as.integer(value)[1]
			if(x@device!=value) {
				if(x@device!=curD)
					setDevice(x@device, silent=TRUE)
				x=h(x)
				setDevice(value, silent=TRUE)
				x=g(x)
				if(getDevice()!=curD)
					setDevice(curD, silent=TRUE)
			}
			return(x)
		})


#type = function(x) {
#	stop(paste("'type' not defined for class:",class(x)))
#}
setGeneric("type", function(x)
			standardGeneric("type"))
setMethod("type", "gmatrix",
		function(x) {
			return(.type_name(x@type))
		})
setMethod("type", "gvector",
		function(x) {
			return(.type_name(x@type))
		})

setGeneric("type<-", function(x, value)
			standardGeneric("type<-"))
setReplaceMethod("type", "gmatrix",
		function(x, value) {
			checkDevice(x@device)
			typeno=.type_num(value)
			ret=convertType(x,typeno,dup=FALSE)
			return(ret)
		})

setReplaceMethod("type", "gvector",
		function(x, value) {
			checkDevice(x@device)
			typeno=.type_num(value)
			ret=convertType(x,typeno,dup=FALSE)
			return(ret)
		})

setMethod("show","gvector",
		function(object) {

			#checkDevice(object@device)
			cat(paste("gvector of length", object@length, "and type", sQuote(type(object)),"on gpu", device(object)), ":\n",sep="")
			#print(object@ptr)
			if(length(object)>0) {
				curD=getDevice()
				flg=FALSE
				if(object@device!=curD) {
					setDevice(object@device, silent=TRUE)
					flg=TRUE
				}
				tmp=h(object)
				if(flg)
					setDevice(curD, silent=TRUE)
				if(length(object)>20)
					cat("printing first 20 elements: \n")
				print(tmp[1:min(20,length(object))])
			}
		})



setMethod("t", "gvector", 
		function(x) {
			return(t(gmatrix(x)))
		}
)

setMethod("[", "gvector", 
                function(x, i, j,...,drop=TRUE) {
                        if(!missing(j)){
                                stop("incorrect number of dimenstions")
                        }
			if(missing(i))
				return(gdup(x))
                        checkDevice(x@device)
                        i=.check_make_valid(i, length(x),names(x))
                        ret=new("gvector", 
                                        ptr=.Call("gpu_numeric_index", x@ptr, length(x), i@ptr, length(i), x@type),
                                        length=length(i), type=x@type)

                        if(!is.null(names(x)))
                                names(ret)=names(x)[as.integer(i)]

                        return(ret)
                }
)



setReplaceMethod("[", "gvector", 
                function(x, i, j,..., value) {

                        if(!missing(j)){
                                stop("incorrect number of dimenstions")
                        }
                        if(!(class(value) %in% c("gvector","gmatrix")))
                                value=as.gvector(value, type=x@type, dup=FALSE)
                        checkDevice(c(x@device, value@device))
                        if(x@type!=value@type)
                                type(value)=x@type
			if(missing(i)){
				if(length(value)!=length(x))
					stop("Number of items to replace is not equal to the number of items")
				tmp=.Call("gpu_cpy",value@ptr, x@ptr, length(x),x@type)
				return(x)
			}
				
                        i=.check_make_valid(i, length(x),names(x),TRUE)

                        #gpu_numeric_index_set(SEXP A_in, SEXP n_A_in, SEXP val_in, SEXP n_val_in, SEXP index_in, SEXP n_index_in)
                        junk=.Call("gpu_numeric_index_set", x@ptr, length(x), value@ptr, length(value), i@ptr, length(i), x@type)


                        return(x)
                }
)

as.matrix.gvector =  function(x, ...) {
	checkDevice(x@device)
	ret=.gpu_get( x@ptr, x@length, x@type)
	if(length(names(x))>0)
		names(ret)=names(x)

	return(as.matrix(ret))
}

setMethod("as.matrix", signature(x = "gvector"), as.matrix.gvector 	)

