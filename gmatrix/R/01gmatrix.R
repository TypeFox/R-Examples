# TODO: gmatrix definitions
# 
# Author: nmorris
###############################################################################


setMethod("initialize",
		"gmatrix",
		function(.Object,...) {
			#cat("test\n")
			#browser()
			.Object <- callNextMethod()
			.Object@device=.Call("get_device")
			return(.Object)
		})

#.type_name = function(type) {
#	if(type==0L) 
#		return("double")
#	else if(type==1L) 
#		return("single")
#	else if(type==2L) 
#		return("integer")
#	else if(type==3L) 
#		return("logical")
#	else if(type=="double") 
#		return(0L)
#	else if(type=="single") 
#		return(1L)
#	else if(type=="integer") 
#		return(2L)
#	else if(type=="logical") 
#		return(3L)
#	else
#		stop("Invalid Type")
#}




gmatrix = function(data = NA, nrow = 1L, ncol = 1L, byrow = FALSE, dimnames = NULL, type=NULL, dup=TRUE) {
	if(class(data)=="gvector") {
		checkDevice(data@device)
		if(nrow<1 || ncol<1)
			stop("gmatrix must have strictly posative rows and columns.")
		if (missing(nrow)) 
			nrow <- ceiling(length(data)/ncol)
		else if (missing(ncol)) 
			ncol <- ceiling(length(data)/nrow)
		if(nrow==1L) {
			rn=NULL
			cn=names(data)
		} else if(ncol==1L) {
			rn=names(data)
			cn=NULL
		} else {
			cn=NULL
			rn=NULL
		}
		if(length(data) != (nrow*ncol) )
			stop("If data is of class gvector, then length(data) must be exactly (nrow*ncol).")
		checkDevice(data@device)
		
		ret = new("gmatrix",nrow=as.integer(nrow),ncol=as.integer(ncol),
				rownames=rn,colnames=cn, type=data@type)
		
		if(is.null(type)){
			if(dup)
				ret@ptr = gdup(data)@ptr
			else
				ret@ptr = data@ptr
		} else {
			ret@ptr = data@ptr
			ret = convertType(ret,to=type, dup=dup)
		}
	
	} else {
		if (missing(nrow)) 
			nrow <- ceiling(length(data)/ncol)
		else if (missing(ncol)) 
			ncol <- ceiling(length(data)/nrow)
		tmp=matrix(data = data, nrow = nrow, ncol = ncol, byrow = byrow, dimnames = dimnames)
		tryCatch(fromtype<- .Rclass_to_type(data),
				error=function(i) stop("Class of data cannot be converted to gpu type.", call. = FALSE))
		ret = new("gmatrix",nrow=nrow(tmp),ncol=ncol(tmp),
				rownames=rownames(tmp),colnames=colnames(tmp), type=fromtype)
		
		ret@ptr = .Call("gpu_create", tmp, fromtype)
		if(!is.null(type))
			ret=convertType(ret,to=type, dup=FALSE)
	}
	return(ret)
}

as.matrix.gmatrix =  function(x, ...) {
	checkDevice(x@device)
	ret=.gpu_get( x@ptr, x@nrow * x@ncol, x@type)
	dim(ret)=c(x@nrow , x@ncol)
	rownames(ret)=x@rownames
	colnames(ret)=x@colnames

	return(ret)
}

setMethod("as.matrix", signature(x = "gmatrix"), as.matrix.gmatrix)



setGeneric("as.gmatrix", useAsDefault=gmatrix)
setMethod("as.gmatrix", "gmatrix",
		function(data, type=data@type, dup=TRUE) {
			checkDevice(data@device)
			typeno=.type_num(type)
			if(data@type==typeno) {
				if(dup)
					ret = gdup(data)
				else
					ret = data
			} else {
				ret=convertType(data,typeno,dup)
			}
			return(ret)
		})

setMethod("as.gmatrix", "gvector",
		function(data, type=data@type, dup=TRUE) {
			checkDevice(data@device)
			if(length(data)<1L)
				stop("'gmatrix' must have strictly posative rows and columns.")
			typeno=.type_num(type)
			ret = new("gmatrix",nrow=length(data), ncol=1L,
					rownames=names(data), type=typeno)
			if(data@type==typeno) {
				if(dup)
					ret@ptr = gdup(data)@ptr
				else
					ret@ptr = data@ptr
			} else {
				ret@ptr=convertType(data,typeno)@ptr
			}
			return(ret)
		})
setMethod("as.gmatrix", "matrix",
		function(data,type=NULL, dup=TRUE) {
			if(nrow(data)<1 || ncol(data)<1)
				stop("'gmatrix' must have strictly posative rows and columns.")
			
			tryCatch(fromtype<- .Rclass_to_type(data),
					error=function(i) stop("Class of data cannot be converted to gpu type.", call. = FALSE))
			ret = new("gmatrix",nrow=nrow(data),ncol=ncol(data),
					rownames=rownames(data),colnames=colnames(data), type=fromtype)
			ret@ptr = .Call("gpu_create", data, fromtype)
			if(!is.null(type)) {
				typeno=.type_num(type)
				if(	.type_num(type)!=fromtype)
					ret=convertType(ret,typeno)
			}	
			return(ret)
		})

setMethod("as.numeric", "gmatrix",
		function(x) {
			checkDevice(x@device)
			#browser()
			if(x@type!=0)
				x=convertType(x,"double")
			ret=.gpu_get( x@ptr, x@nrow * x@ncol, x@type)
			#length(ret)=length(x)
			return(ret)
		})

as.vector.gmatrix =  function(x,mode=NULL) {
	checkDevice(x@device)

	ret=.gpu_get( x@ptr, x@nrow * x@ncol, x@type)
	if(!is.null(mode))
		mode(ret)=mode
	return(ret)
}


setMethod("as.vector", "gmatrix", as.vector.gmatrix)

setMethod("as.integer", "gmatrix",
		function(x) {
			checkDevice(x@device)
			ret=.gpu_get( x@ptr, x@nrow * x@ncol, x@type)
			return(as.integer(ret))
		})

setMethod("as.logical", "gmatrix",
		function(x) {
			checkDevice(x@device)
			ret=.gpu_get( x@ptr, x@nrow * x@ncol, x@type)
			return(as.logical(ret))
		})

setMethod("rownames", "gmatrix",
		function(x) {
			return(x@rownames)
		})

setMethod("colnames", "gmatrix",
		function(x) {
			return(x@colnames)
		})
setMethod("dimnames", "gmatrix",
		function(x) {
			return(list(x@rownames,x@colnames))
		})
setReplaceMethod("rownames", "gmatrix",
		function(x, value) {
			if(is.null(value)) 
				x@rownames=NULL
			else if(length(value)!=nrow(x))
				stop("length of rownames does not match the dimension of the array")
			else
				x@rownames=as.character(value)
			return(x)
		}
)
setReplaceMethod("colnames", "gmatrix",
		function(x, value) {
			if(is.null(value)) 
				x@colnames=NULL
			else if(length(value)!=ncol(x))
				stop("length of colnames does not match the dimension of the array")
			else
				x@colnames=as.character(value)
			return(x)
		}
)
setReplaceMethod("dimnames", "gmatrix",
		function(x, value) {
			if(is.null(value)) {
				rownames(x)=NULL
				colnames(x)=NULL
			} else if(length(value)!=2) {
				stop("dimnames must be a list of length2")
			} else {
				rownames(x)=value[[1]]
				colnames(x)=value[[2]]
			}
			return(x)
		}
)

#TODO: as.gmatrix for gpuvector
setMethod("dim",  "gmatrix",
		function(x) {
			return(c(x@nrow,x@ncol))
		})

setReplaceMethod("dim", "gmatrix",
		function(x, value) {
			value=as.integer(value)
			if(length(value)==1){
				if(value!=length(x))
					stop("Matrix dimensions result in different number of total elements.")
				x=new("gvector", ptr=x@ptr, length=length(x))
			} else if(length(value)==2) {
				if(value[1]*value[2] != x@nrow*x@ncol)
					stop("Matrix dimensions result in different number of total elements.")
				x@nrow <- as.integer(value[1])
				x@ncol <- as.integer(value[2])
			} else
				stop("'dim' must have a length of 1 or 2.")
			return(x)
		}
)

setMethod("length",  "gmatrix",
		function(x) {
			return(as.integer(x@nrow*x@ncol))
		})

setMethod("nrow",  "gmatrix",
		function(x) {
			return(x@nrow)
		})



setMethod("ncol", "gmatrix",
		function(x) {
			return(x@ncol)
		})

setMethod("show","gmatrix",
		function(object) {
			#checkDevice(object@device)
			cat(paste("gmatrix of size", object@nrow,"x", object@ncol, "and type", sQuote(type(object)),"on device", device(object)), ":\n",sep="")
			#print(object@ptr)
			curD=getDevice()
			flg=FALSE
			if(object@device!=curD) {
				setDevice(object@device, silent=TRUE)
				flg=TRUE
			}
			tmp=h(object)
			if(flg)
				setDevice(curD, silent=TRUE)
			if(nrow(object)>10 && ncol(object)>10)
				cat("printing upper left corner of the matrix: \n")
			else if(nrow(object)<=10 && ncol(object)>10)
				cat("printing first 10 columns of the matrix: \n")
			else if(nrow(object)>10 && ncol(object)<=10)
				cat("printing first 10 rows of the matrix: \n")
			print(tmp[1:min(10,nrow(object)),1:min(10,ncol(object)),drop=FALSE])
			
		})


setMethod("t", "gmatrix", 
		function(x) {
			checkDevice(x@device)
			return(
					new("gmatrix",
							ptr= .Call("gpu_naive_transpose", x@ptr, nrow(x), ncol(x), x@type),
							nrow=ncol(x), ncol=nrow(x),
							rownames=colnames(x), colnames=rownames(x),type=x@type)
			)
		}
)

setMethod("diag", "gmatrix", 
		function(x) {
			checkDevice(x@device)
			#gpu_diag_get(SEXP A_in, SEXP n_row_in, SEXP n_col_in)
			return(
					new("gvector",
							ptr= .Call("gpu_diag_get", x@ptr, nrow(x), ncol(x), x@type),
							length=min(ncol(x),nrow(x)), type=x@type)
			)
		}
)
setMethod("diag", "gvector", 
		function(x) {
			checkDevice(x@device)
			#gpu_diag_get(SEXP A_in, SEXP n_row_in, SEXP n_col_in)
			ret = gmatrix(0, length(x), length(x), type=x@type)
			diag(ret)=x
			return(ret)
		}
)

gident = function(n,val=1, type="d") {
	n=as.integer(n)[1]
	ret = gmatrix(0, n, n, type=type)
	diag(ret)=val
	return(ret)
}

setReplaceMethod("diag", "gmatrix", 
		function(x, value) {
			checkDevice(x@device)
			if(length(value)==1) {
				tryCatch(value <- .convert_to_appropriate_class(value,x@type), error=function(e) stop("Invalid value for diag.", call. = FALSE))
				tmp=.Call("gpu_diag_set_one", x@ptr, nrow(x), ncol(x),value, x@type)
				return(x)
			} else {
				if(class(value)!="gvector")
					value=as.gvector(value)
				if(x@type!=value@type)
					type(value)=x@type
				tmp=.Call("gpu_diag_set", x@ptr, nrow(x), ncol(x),value@ptr,length(value), x@type)
				return(x)
			}
		}
)




#setMethod("[", signature(x = "CsparseMatrix", i = "index", j = "missing",
#				drop = "logical"),

#.is.logical.gvector= function(x) {
#	if(class(x)="gvector"){
#		if(x@type=3L)
#			return(TRUE)
#		else
#			return(FALSE)
#	} else
#		return(FALSE)
#}
#
#.is.integer.gvector=function(x) {
#	if(class(x)="gvector"){
#		if(x@type=2L)
#			return(TRUE)
#		else
#			return(FALSE)
#	} else
#		return(FALSE)
#}
.check_make_valid =function(i, size, nmx, chk=FALSE) {
	if(is.vector(i) ) {
		if(is.logical(i) )
			i=as.gvector(which(i), type=2L)
		else if(is.character(i)) {
			i=as.gvector(match(i,nmx), type=2L)
		} else if(min(i,  na.rm = TRUE)<1)
			i=as.gvector((1:size)[i], type=2L)
		else
			i=as.gvector(i, type=2L)
	} else	if(class(i)=="gvector") {
		checkDevice(i@device)
		f=function(x) {if(is.na(x)) return(x)
			else return(sum(h(i),na.rm=TRUE))}
		if(i@type==3L) {
			i=which(i)
		} else if(f(min(i,retgpu=FALSE))<1) {
			i=as.gvector((1:size)[as.integer(i)], type=2L) #works but may be bottle neck in some cases
		} else
			i=as.gvector(i, type=2L, dup=FALSE) 
	} else
		stop("Gpu based indexing can only be performed with an index of class 'vector' or 'gvector'.")
	if(length(i)<1)
		stop("Index on gpu must select at least one thing.")
	if(chk) {
		maxi=max(i, retgpu=FALSE)
		mini=min(i, retgpu=FALSE)
		if(is.finite(maxi) && is.finite(mini)) {
			if(mini<1 || maxi>size )
				stop("The index cannot be outside the boundaries.")
		} else
			stop("The index must be finite and not missing.")
	}
	return(i)
}

.twoindex=	function(x, i, j, ..., drop=TRUE) {
	checkDevice(x@device)
	
	j=.check_make_valid(j,ncol(x),colnames(x),TRUE)
	i=.check_make_valid(i,nrow(x),rownames(x),TRUE)
	#browser()
	#SEXP gpu_gmatrix_index_both(SEXP A_in, SEXP n_row_A_in, SEXP n_col_A_in,
	#		SEXP index_row_in, SEXP n_index_row_in,SEXP index_col_in, SEXP n_index_col_in)
	ret=new("gmatrix", 
			ptr=.Call("gpu_gmatrix_index_both", x@ptr, nrow(x), ncol(x), i@ptr, length(i), j@ptr, length(j), x@type),
			nrow=length(i), ncol=length(j), type=x@type)
	
	
	
	if(!is.null(rownames(x)))
		rownames(ret)=rownames(x)[as.integer(i)]
	if(!is.null(colnames(x)))
		colnames(ret)=colnames(x)[as.integer(j)]
	if(drop) {
		if(nrow(ret)==1)
			ret=new("gvector",ptr=ret@ptr, length=ncol(ret), names=colnames(ret), type=ret@type)
		else if(ncol(ret)==1)
			ret=new("gvector",ptr=ret@ptr, length=nrow(ret), names=rownames(ret), type=ret@type)
	}
	
	return(ret)
}

.colindex=	function(x, i, j, ..., drop=TRUE) {
	checkDevice(x@device)
	
	j=.check_make_valid(j,ncol(x),colnames(x), TRUE)
	ret=new("gmatrix", 
			ptr=.Call("gpu_gmatrix_index_col", x@ptr, nrow(x), ncol(x), j@ptr, length(j),  x@type),
			nrow=nrow(x), ncol=length(j), type=x@type)
	#browser()
	if(!is.null(rownames(x)))
		rownames(ret)=rownames(x)
	if(!is.null(colnames(x)))
		colnames(ret)=colnames(x)[as.integer(j)]
	if(drop) {
		if(nrow(ret)==1)
			ret=new("gvector",ptr=ret@ptr, length=ncol(ret), names=colnames(ret), type=ret@type)
		else if(ncol(ret)==1)
			ret=new("gvector",ptr=ret@ptr, length=nrow(ret), names=rownames(ret), type=ret@type)
	}
	
	return(ret)
}



.rowindex=	function(x, i, j, ..., drop=TRUE) {
	checkDevice(x@device)
	if(nargs() == 2) { ## e.g. M[0] , M[TRUE],  M[1:2]
		return(.oneindex(x, i))
	} else {
		i=.check_make_valid(i,nrow(x),rownames(x),TRUE)
		ret=new("gmatrix", 
				ptr=.Call("gpu_gmatrix_index_row", x@ptr, nrow(x), ncol(x), i@ptr, length(i), x@type),
				nrow=length(i), ncol=ncol(x),type= x@type)
		
		if(!is.null(rownames(x)))
			rownames(ret)=rownames(x)[as.integer(i)]
		if(!is.null(colnames(x)))
			colnames(ret)=colnames(x)
		if(drop) {
			if(nrow(ret)==1)
				ret=new("gvector",ptr=ret@ptr, length=ncol(ret), names=colnames(ret), type=ret@type)
			else if(ncol(ret)==1)
				ret=new("gvector",ptr=ret@ptr, length=nrow(ret), names=rownames(ret), type=ret@type)
		}
	}
	return(ret)
}





.setrowindex = function(x, i, j,..., value) {	
	if(!(class(value) %in% c("gvector","gmatrix")))
		value=as.gvector(value, type=x@type)
	checkDevice(c(x@device,value@device))
	if(x@type!=value@type)
		type(value)=x@type
	if(nargs() == 3) { ## e.g. M[0] , M[TRUE],  M[1:2]
		return(.setoneindex(x, i, value))
	} else {
		i=.check_make_valid(i,nrow(x),names(x), TRUE)
		tmp=.Call("gpu_gmatrix_index_row_set", x@ptr, nrow(x),ncol(x), value@ptr, length(value), i@ptr, length(i), x@type)
	}
	return(x)
}

.setcolindex = function(x, i, j,..., value) {
	
	if(!(class(value) %in% c("gvector","gmatrix")))
		value=as.gvector(value, type=x@type)
	checkDevice(c(x@device,value@device))
	if(x@type!=value@type)
		type(value)=x@type
	
	j=.check_make_valid(j,ncol(x),names(x),TRUE)				
	#gpu_gmatrix_index_col_set(SEXP A_in, SEXP n_row_A_in, SEXP n_col_A_in,  SEXP val_in, SEXP n_val_in, SEXP index_in, SEXP n_index_in)
	tmp=.Call("gpu_gmatrix_index_col_set", x@ptr, nrow(x),  ncol(x),value@ptr, length(value), j@ptr, length(j), x@type)
	
	return(x)
}


.settwoindex = function(x, i, j,..., value) {
	
	if(!(class(value) %in% c("gvector","gmatrix")))
		value=as.gvector(value, type=x@type)
	checkDevice(c(x@device,value@device))
	if(x@type!=value@type)
		type(value)=x@type
	
	i=.check_make_valid(i,nrow(x),names(x),TRUE)
	j=.check_make_valid(j,ncol(x),names(x),TRUE)
	
	#SEXP gpu_gmatrix_index_both(SEXP A_in, SEXP n_row_A_in, SEXP n_col_A_in,
	#		SEXP index_row_in, SEXP n_index_row_in,SEXP index_col_in, SEXP n_index_col_in)
	tmp=.Call("gpu_gmatrix_index_both_set", x@ptr, nrow(x), ncol(x), value@ptr, length(value), i@ptr, length(i), j@ptr, length(j), x@type)
	return(x)
}


.oneindex=function(x, i) {
	checkDevice(x@device)
	if(is.character(i))
		stop("Index for this operation cannot be a character vector.")
	i=.check_make_valid(i,length(x),names(x))
	ret=new("gvector", 
			ptr=.Call("gpu_numeric_index", x@ptr, length(x), i@ptr, length(i), x@type),
			length=length(i),type= x@type)
	
	if(!is.null(names(x)))
		names(ret)=names(x)[as.integer(i)]
	
	return(ret)
}

.setoneindex=function(x, i, value) {
	
	if(is.character(i))
		stop("Index for this operation cannot be a character vector.")
	if(!(class(value) %in% c("gvector","gmatrix")))
		value=as.gvector(value, type=x@type)
	i=.check_make_valid(i,length(x),names(x), TRUE)
	checkDevice(c(x@device,value@device))
	if(x@type!=value@type)
		type(value)=x@type
	#gpu_numeric_index_set(SEXP A_in, SEXP n_A_in, SEXP val_in, SEXP n_val_in, SEXP index_in, SEXP n_index_in)
	junk=.Call("gpu_numeric_index_set", x@ptr, length(x), value@ptr, length(value), i@ptr, length(i), x@type)
	
	return(x)
}


setMethod("[", c("gmatrix","index","index"),    .twoindex)
setMethod("[",  c("gmatrix","index","missing"), .rowindex)
setMethod("[",  c("gmatrix","missing","index"), .colindex)
setMethod("[", c("gmatrix","missing","missing"),    function(x, i, j, ..., drop=TRUE) return(x))
setReplaceMethod("[", c("gmatrix","index","index"), .settwoindex)
setReplaceMethod("[", c("gmatrix","index","missing"), .setrowindex)
setReplaceMethod("[", c("gmatrix","missing","index"), .setcolindex)
setReplaceMethod("[", c("gmatrix","missing","missing"), 
	function(x, i, j,..., value) {
		if(!(class(value) %in% c("gvector","gmatrix")))
                	value=as.gvector(value, type=x@type, dup=FALSE)
		if(length(value)!=length(x))
			stop("Number of items to replace is not equal to the number of items")
		tmp=.Call("gpu_cpy",value@ptr, x@ptr, length(x),x@type)
		return(x)
	}
)







