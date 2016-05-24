# Various Operations
# 
# Author: nmorris
###############################################################################

.exprs_xy = bquote( {
			if(x@type!=y@type) {
				totype=min(c(x@type,y@type))
				if(x@type!=totype)
					x=convertType(x,totype)
				if(y@type!=totype)
					y=convertType(y,totype)
			} else if (x@type==3L) {
				x=convertType(x,2L)
				y=convertType(y,2L)
			}
			checkDevice(c(y@device,x@device))
		})

#			if(max(c(x@type,y@type))>1)
#				stop("Operation can only be performed on 'double' or 'single' type.", call. = FALSE)
.exprs_sf_xy = bquote({
			totype=min(c(x@type,y@type))
			if(x@type!=y@type || totype>1L) {
				if(totype>1L)
					totype=0L
				if(x@type!=totype)
					x=convertType(x,totype)
				if(y@type!=totype)
					y=convertType(y,totype)
			} 
			checkDevice(c(y@device,x@device))
		})

.exprs_compare_xy = bquote( {
			if(x@type!=y@type) {
				tmp=c(x@type,y@type)
				totype=min(tmp)
				if(totype==0 && max(tmp)==1)
					totype=1L #these comparisons lead to counterintuitive results otherwise
				if(x@type!=totype)
					x=convertType(x,totype)
				if(y@type!=totype)
					y=convertType(y,totype)
			}
			checkDevice(c(y@device,x@device))
		})

.exprs_l_xy = bquote({
				if(x@type!=3L)
					x=convertType(x,3L)
				if(y@type!=3L)
					y=convertType(y,3L)
				checkDevice(c(y@device,x@device))
			
		})
.exprs_gpu_cpu_xy = bquote({
			type2 =.Rclass_to_type(y)
			if(x@type!=type2) {
				if(x@type<2L) {
					y=as.double(y)					
				}else {
					totype=min(c(x@type,type2))
					if(x@type!=totype)
						x=convertType(x,totype)
					if(type2!=totype ) {
						if(totype==0L || totype==1L)
							y=as.double(y)
						else if( totype==2L)
							y=as.integer(y)
						else
							y=as.logical(y)
					}
				}
			} else if(type2==3L) {#both inputs are logical
				y=as.integer(y)
				x=convertType(x,2L)
			} 
			checkDevice(x@device)
		})

.exprs_sf_gpu_cpu_xy = bquote({ #the case where everything must be double or single
			y=as.double(y)
			if(x@type>1L) #if x is not single, then make it double
				x=convertType(x,0L)
			checkDevice(x@device)
		})
.exprs_l_gpu_cpu_xy = bquote({ #for logical operators
			if(x@type!=3L)
				x=convertType(x,3L)
			y=as.logical(y)
			checkDevice(x@device)
		})


###################################
#    matrix multiplication
###################################
setMethod("%*%", signature(x = "gmatrix", y = "gmatrix"), 
		function(x,y) {
			#browser()
			eval(.exprs_sf_xy)
			return(.Call("matrix_multiply", x, y, FALSE, FALSE, x@type))
		})

setMethod("%*%", signature(x = "gmatrix", y = "gvector"), 
		function(x,y) {
			eval(.exprs_sf_xy)
			ytmp=new("gmatrix",nrow=length(y),ncol=1L,ptr=y@ptr, type=y@type)
			return(.Call("matrix_multiply", x, ytmp, FALSE, FALSE, x@type))
		})
setMethod("%*%", signature(x = "gvector", y = "gmatrix"), 
		function(x,y) {
			eval(.exprs_sf_xy)
			xtmp=new("gmatrix",ncol=length(x),nrow=1L,ptr=x@ptr, type=x@type)
			return(.Call("matrix_multiply", xtmp, y, FALSE, FALSE, x@type))
		})
setMethod("%*%", signature(x = "gmatrix", y = "numeric"), 
		function(x,y) {
			if(x@type==1L) {
				y=as.gvector(y, type=1L)
			} else
				y=as.gvector(y)
			return(x%*%y)
		})
setMethod("%*%", signature(x = "numeric", y = "gmatrix"), 
		function(x,y) {
			if(y@type==1L) {
				x=as.gvector(x, type=1L)
			} else
				x=as.gvector(x)
			return(x%*%y)
		})
setMethod("%*%", signature(x = "gmatrix", y = "matrix"), 
		function(x,y) {
			if(x@type==1L) {
				y=as.gmatrix(y, type=1L)
			} else
				y=as.gmatrix(y)
			return(x%*%y)
		})
setMethod("%*%", signature(x = "matrix", y = "gmatrix"), 
		function(x,y) {
			if(y@type==1L) {
				x=as.gmatrix(x, type=1L)
			} else
				x=as.gmatrix(x)
			return(x%*%y)
		})

setMethod("%*%", signature(x = "gmatrix", y = "logical"), 
		function(x,y) {
			if(x@type==1L) {
				y=as.gvector(y, type=1L)
			} else
				y=as.gvector(y)
			return(x%*%y)
		})
setMethod("%*%", signature(x = "logical", y = "gmatrix"), 
		function(x,y) {
			if(y@type==1L) {
				x=as.gvector(x, type=1L)
			} else
				x=as.gvector(x)
			return(x%*%y)
		})

setMethod("%*%", signature(x = "gvector", y= "gvector"), 
		function(x,y) {
			#browser()
			return(gmatrix(x,nrow=1,dup=FALSE) %*% gmatrix(y,ncol=1,dup=FALSE))
		})
setMethod("%*%", signature(x = "numeric", y= "gvector"), 
		function(x,y) {
			#browser()
			if(y@type==1L) 
				return(gmatrix(x,nrow=1,dup=FALSE, type=1L) %*% gmatrix(y,ncol=1,dup=FALSE))
			else
				return(gmatrix(x,nrow=1,dup=FALSE) %*% gmatrix(y,ncol=1,dup=FALSE))
		})
setMethod("%*%", signature(x = "gvector", y= "numeric"), 
		function(x,y) {
			#browser()
			if(x@type==1L) 
				return(gmatrix(x,nrow=1,dup=FALSE) %*% gmatrix(y,ncol=1,dup=FALSE, type=1L))
			else
				return(gmatrix(x,nrow=1,dup=FALSE) %*% gmatrix(y,ncol=1,dup=FALSE))
		})
setMethod("%*%", signature(x = "logical", y= "gvector"), 
		function(x,y) {
			if(y@type==1L) 
				return(gmatrix(x,nrow=1,dup=FALSE, type=1L) %*% gmatrix(y,ncol=1,dup=FALSE))
			else
				return(gmatrix(x,nrow=1,dup=FALSE) %*% gmatrix(y,ncol=1,dup=FALSE))
		})
setMethod("%*%", signature(x = "gvector", y= "logical"), 
		function(x,y) {
			if(x@type==1L) 
				return(gmatrix(x,nrow=1,dup=FALSE) %*% gmatrix(y,ncol=1,dup=FALSE, type=1L))
			else
				return(gmatrix(x,nrow=1,dup=FALSE) %*% gmatrix(y,ncol=1,dup=FALSE))
		})

setMethod("%*%", signature(x = "gvector", y= "matrix"), 
		function(x,y) {
			if(x@type==1L) {
				y=as.gmatrix(y, type=1L)
			} else
				y=as.gmatrix(y)
			return(x%*%y)
		})
setMethod("%*%", signature(x = "matrix", y= "gvector"), 
		function(x,y) {
			if(y@type==1L) {
				x=as.gmatrix(x, type=1L)
			} else
				x=as.gmatrix(x)
			return(x%*%y)
		})

##############################################
#       GMM function to avoid the allocating C
##############################################
gmm = function(A, B, C, trA=FALSE, trB=FALSE, accum=FALSE) {
  if(class(A)=="gvector")
    A=gmatrix(A,nrow=1,dup=FALSE)
  else if(class(A)!="gmatrix")
    stop("A must be of class 'gmatrix' or 'gvector.'")
  if(class(B)=="gvector")
    B=gmatrix(B,ncol=1,dup=FALSE)
  else if(class(B)!="gmatrix")
    stop("B must be of class 'gmatrix' or 'gvector.'")
  if(class(C)=="gvector")
    C=gmatrix(C,ncol=1,dup=FALSE)
  else if(class(C)!="gmatrix")
    stop("C must be of class 'gmatrix' or 'gvector.'")
  accum=as.logical(accum)[1]
  
  if(A@type!=B@type || B@type !=C@type)
    stop("A, B and C must all have the same type.")
  if(A@type>1)
    stop("The 'gmm' function may only be used for type 'single' or 'double.' Please convert first.")
  # SEXP gpu_gmm(SEXP A_in, SEXP B_in, SEXP C_in, SEXP transa, SEXP transb,  SEXP accum, SEXP in_type)
  invisible(.Call("gpu_gmm", A, B, C, trA, trB, accum, A@type))
}


##############################################
#              Cross Product
##############################################
setMethod("crossprod", signature(x = "gmatrix", y = "missing"), 
		function(x,y){ 
			if(x@type>1L)
				type(x)=0
			checkDevice(x@device)
			return(.Call("matrix_multiply", x, x, TRUE, FALSE,x@type))
		})
setMethod("crossprod", signature(x = "gmatrix", y = "gmatrix"), 
		function(x,y){ 
			eval(.exprs_sf_xy)
			return(.Call("matrix_multiply", x, y, TRUE, FALSE,x@type))
		})

setMethod("crossprod", signature(x = "gmatrix", y = "matrix"), 
		function(x,y){
			if(x@type==1L)
				y=g(y, type=1L)
			else
				y=g(y)
			eval(.exprs_sf_xy)
			return(.Call("matrix_multiply", x, y, TRUE, FALSE, x@type))
		})

setMethod("crossprod", signature(x = "matrix", y = "gmatrix"), 
		function(x,y){ 
			if(y@type==1L)
				x=g(x, type=1L)
			else
				x=g(x)
			eval(.exprs_sf_xy)
			return(.Call("matrix_multiply",  x, y, TRUE, FALSE,x@type))
		})

setMethod("crossprod", signature(x = "matrix", y = "gmatrix"), 
		function(x,y){ 
			if(y@type==1L)
				x=g(x, type=1L)
			else
				x=g(x)
			eval(.exprs_sf_xy)
			return(.Call("matrix_multiply",  x, y, TRUE, FALSE,x@type))
		})


.cp = function(x,y) {
	flgx=!(class(x) %in% c("gmatrix", "gvector"))
	flgy=!(class(y) %in% c("gmatrix", "gvector"))
	if(flgx) {
		if(!flgy) {
			if(y@type==1)
				x=g(x,type=1L)
			else
				x=g(x)
		} else
			x=g(x)
	}
	if(flgy) {
		if(!flgx) {
			if(x@type==1)
				y=g(y,type=1L)
			else
				y=g(y)
		} else
			y=g(y)
	}
	eval(.exprs_sf_xy)
	if(class(x)=="gvector"){
		x=gmatrix(x,nrow=1,dup=FALSE)
		if(class(y)=="gvector") 
			y=gmatrix(y,ncol=1,dup=FALSE)
		return(.Call("matrix_multiply",  x, y, FALSE, FALSE,x@type))
	}
	#if here then y is a gvector (b/c one of the two must be a gvector) and x is a gmatrix
	y=gmatrix(y,ncol=1,dup=FALSE)
	return(.Call("matrix_multiply",  x, y, TRUE, FALSE,x@type))
}
setMethod("crossprod", signature(x = "numeric", y = "gmatrix"), .cp)
setMethod("crossprod", signature(x = "gmatrix", y = "numeric"),  .cp)
setMethod("crossprod", signature(x = "logical", y = "gmatrix"),  .cp)
setMethod("crossprod", signature(x = "gmatrix", y = "logical"),  .cp)
setMethod("crossprod", signature(x = "gvector", y = "gmatrix"),  .cp)
setMethod("crossprod", signature(x = "gmatrix", y = "gvector"),  .cp)
setMethod("crossprod", signature(x = "gvector", y = "gvector"),  .cp)
setMethod("crossprod", signature(x = "gvector", y = "missing"),  function(x,y) {.cp(x,x)})
setMethod("crossprod", signature(x = "gvector", y = "numeric"), .cp)
setMethod("crossprod", signature(x = "numeric", y = "gvector"),  .cp)
setMethod("crossprod", signature(x = "gvector", y = "logical"),  .cp)
setMethod("crossprod", signature(x = "logical", y = "gvector"),  .cp)
setMethod("crossprod", signature(x = "gvector", y = "matrix"), .cp)
setMethod("crossprod", signature(x = "matrix", y = "gvector"), .cp)

##############################################
#              Transpose Cross Product
##############################################

setMethod("tcrossprod", signature(x = "gmatrix", y = "missing"), 
		function(x,y){
			if(x@type>1L)
				type(x)=0
			checkDevice(x@device)
			return(.Call("matrix_multiply", x, x, FALSE, TRUE,x@type))
		})
setMethod("tcrossprod", signature(x = "gmatrix", y = "gmatrix"), 
		function(x,y){
			eval(.exprs_sf_xy)
			return(.Call("matrix_multiply", x, y, FALSE, TRUE,x@type))
		})
setMethod("tcrossprod", signature(x = "gmatrix", y = "matrix"), 
		function(x,y){
			if(x@type==1L)
				y=g(y, type=1L)
			else
				y=g(y)
			eval(.exprs_sf_xy)
			return(.Call("matrix_multiply",x, y, FALSE, TRUE, x@type))
		})
setMethod("tcrossprod", signature(x = "matrix", y = "gmatrix"), 
		function(x,y){
			if(y@type==1L)
				x=g(x, type=1L)
			else
				x=g(x)
			eval(.exprs_sf_xy)
			return(.Call("matrix_multiply",x ,y, FALSE, TRUE,  x@type))
		})

.tcp = function(x,y) {
	flgx=!(class(x) %in% c("gmatrix", "gvector"))
	flgy=!(class(y) %in% c("gmatrix", "gvector"))
	if(flgx) {
		if(!flgy) {
			if(y@type==1)
				x=g(x,type=1L)
			else
				x=g(x)
		} else
			x=g(x)
	}
	if(flgy) {
		if(!flgx) {
			if(x@type==1)
				y=g(y,type=1L)
			else
				y=g(y)
		} else
			y=g(y)
	}

	eval(.exprs_sf_xy)
	#browser()
	if(class(x)=="gvector"){
		if(class(y)=="gvector") {
			x=gmatrix(x,ncol=1,dup=FALSE)
			y=gmatrix(y,nrow=1,dup=FALSE)
			return(.Call("matrix_multiply",  x, y, FALSE, FALSE, x@type))
		} else {
		#	browser()
			x=gmatrix(x,nrow=1,dup=FALSE)
			return(.Call("matrix_multiply",  x, y, FALSE, TRUE, x@type))
		}
	}
	#if here then y is a gvector (b/c one of the two must be a gvector) and x is a gmatrix
	y=gmatrix(y,ncol=1,dup=FALSE)
	return(.Call("matrix_multiply",  x, y, FALSE, TRUE, x@type))
}

setMethod("tcrossprod", signature(x = "numeric", y = "gmatrix"), .tcp)
setMethod("tcrossprod", signature(x = "gmatrix", y = "numeric"), .tcp)
setMethod("tcrossprod", signature(x = "logical", y = "gmatrix"), .tcp)
setMethod("tcrossprod", signature(x = "gmatrix", y = "logical"), .tcp)
setMethod("tcrossprod", signature(x = "gvector", y = "gmatrix"), .tcp)
setMethod("tcrossprod", signature(x = "gmatrix", y = "gvector"), .tcp)
setMethod("tcrossprod", signature(x = "gvector", y = "gvector"), .tcp)
setMethod("tcrossprod", signature(x = "gvector", y = "missing"), function(x,y) {.tcp(x,x)})
setMethod("tcrossprod", signature(x = "gvector", y = "numeric"), .tcp)
setMethod("tcrossprod", signature(x = "numeric", y = "gvector"), .tcp)
setMethod("tcrossprod", signature(x = "gvector", y = "logical"), .tcp)
setMethod("tcrossprod", signature(x = "logical", y = "gvector"), .tcp)
setMethod("tcrossprod", signature(x = "gvector", y = "matrix"), .tcp)
setMethod("tcrossprod", signature(x = "matrix", y = "gvector"), .tcp)
#####################################################
#   Outer Product
#####################################################
gouter<-function(x,y,FUN="*", normalOrder=TRUE ){
	if(is.vector(x))
		if(class(y)=="gvector") {
			if(y@type==1L)
				x=as.gvector(x,type=1L)
			else
				x=as.gvector(x)
		} else
			x=as.gvector(x)
	if(is.vector(y))
		if(class(x)=="gvector") {
			if(x@type==1L)
				y=as.gvector(y,type=1L)
			else
				y=as.gvector(y)
		} else
			y=as.gvector(y)
	if(class(x)!="gvector" && class(y)!="gvector")
		stop("Expected input types of 'gvector' or 'vector'.")
	eval(.exprs_sf_xy)
	checkDevice(c(x@device,y@device))
	if(class(x)!="gvector")
		x=as.gvector(x)
	if(class(y)!="gvector")
		y=as.gvector(y)
#			if(op==1)
#				ret[i] = x[row] * y[col];
#			else if(op==2)
#				ret[i] = x[row] + y[col];
#			else if(op==3)
#				ret[i] =x[row] - y[col];
#			else if(op==4)
#				ret[i] =y[row] - x[col];
#			else if(op==5)
#				ret[i] =x[row] / y[col];
#			else if(op==6)
#				ret[i] =y[row] / x[col];
#			else if(op==7)
#				ret[i] = pow(x[row] , y[col]) ;
#			else if(op==8)
#				ret[i] = pow(y[row] , x[col]) ;
	if(FUN=="*")
		op=1L
	else if(FUN=="+")
		op=2L
	else if(FUN=="-"){
		if(normalOrder)
			op=3L
		else
			op=4L
	} else if(FUN=="/"){
		if(normalOrder)
			op=5L
		else
			op=6L
	}else if(FUN=="^"){
		if(normalOrder)
			op=7L
		else
			op=8L
	} else
		stop("Invalid function.")

	#gpu_outer(SEXP A_in, SEXP B_in,SEXP n_A_in, SEXP n_B_in, SEXP op_in)
	return(new("gmatrix",
					ptr=.Call("gpu_outer",x@ptr,y@ptr,length(x),length(y),op, x@type),
					nrow=length(x),ncol=length(y),
					rownames=names(x),colnames=names(y),
					type=x@type))
}


setMethod("%o%", signature(X = "gvector", Y= "gvector"), 
		function(X,Y) {
				return(gmatrix(X,ncol=1,dup=FALSE) %*% gmatrix(Y,nrow=1,dup=FALSE))
			
		})
setMethod("%o%", signature(X = "numeric", Y= "gvector"), 
		function(X,Y) {
			if(Y@type==1L) 
				return(gmatrix(X,ncol=1,dup=FALSE, type=1L) %*% gmatrix(Y,nrow=1,dup=FALSE))
			else
				return(gmatrix(X,ncol=1,dup=FALSE) %*% gmatrix(Y,nrow=1,dup=FALSE))
			
		})

setMethod("%o%", signature(X = "gvector", Y= "numeric"), 
		function(X,Y) {
	
			if(X@type==1L) 
				return(gmatrix(X,ncol=1,dup=FALSE) %*% gmatrix(Y,nrow=1,dup=FALSE, type=1L))
			else
				return(gmatrix(X,ncol=1,dup=FALSE) %*% gmatrix(Y,nrow=1,dup=FALSE))
		})

setMethod("%o%", signature(X = "logical", Y= "gvector"), 
		function(X,Y) {
			if(Y@type==1L) 
				return(gmatrix(X,ncol=1,dup=FALSE, type=1L) %*% gmatrix(Y,nrow=1,dup=FALSE))
			else
				return(gmatrix(X,ncol=1,dup=FALSE) %*% gmatrix(Y,nrow=1,dup=FALSE))
			
		})

setMethod("%o%", signature(X = "gvector", Y= "logical"), 
		function(X,Y) {
			if(X@type==1L) 
				return(gmatrix(X,ncol=1,dup=FALSE) %*% gmatrix(Y,nrow=1,dup=FALSE, type=1L))
			else
				return(gmatrix(X,ncol=1,dup=FALSE) %*% gmatrix(Y,nrow=1,dup=FALSE))
		})


################################################
#            Sums and Means
################################################
setMethod("colSums", "gmatrix",
		function(x, na.rm = FALSE, dims = 1) {
			if(na.rm)
				stop("na.rm not implemented on the GPU")
			if(dims!=1)
				stop("invalid 'dims'")
			return(as.gvector(  rep(1,nrow(x)) %*% x   , dup=FALSE))
		}
)



setMethod("rowSums",  "gmatrix",
		function(x, na.rm = FALSE, dims = 1) {
			if(na.rm)
				stop("na.rm not implemented on the GPU")
			if(dims!=1)
				stop("invalid 'dims'")
		#	browser()
			return(as.gvector(x %*% g.rep(1,ncol(x)), dup=FALSE))
		}
)

gRowLogSums = function(x, startCol=1L, endCol=ncol(x)) {
	
	if(class(x)!="gmatrix")
		stop("Object must be of class 'gmatrix.'")
	if(x@type>1L)
		stop("Invalid type")
	startCol=as.integer(startCol)
    endCol=as.integer(endCol)	
	if(startCol<1L || startCol>ncol(x) || startCol>endCol)
		stop("Invalid startCol.")
	if(endCol<1L || endCol>ncol(x) )
		stop("Invalid endCol.")
		
	checkDevice(x@device)
	ret = new("gvector", ptr=.Call("gpu_rowLogSums", ptr=x@ptr,
			nrow(x), endCol, startCol, x@type),
			length=nrow(x),type=x@type)
	return(ret)
}

setMethod("colMeans",  "gmatrix",
		function(x, na.rm = FALSE, dims = 1) {
		#	browser()
			if(na.rm)
				stop("na.rm not implemented on the GPU")
			if(dims!=1)
				stop("invalid 'dims'")
			return(g.rep(1,nrow(x)) %*% x/nrow(x))
		}
)
setMethod("rowMeans",  "gmatrix",
		function(x, na.rm = FALSE, dims = 1) {
			if(na.rm)
				stop("na.rm not implemented on the GPU")
			if(dims!=1)
				stop("invalid 'dims'")
			return(as.gvector(x %*% rep(1,ncol(x)), dup=FALSE)/ncol(x))
		}
)

setMethod("sum",  "gvector",
		function(x, na.rm = FALSE, retgpu=TRUE) {
			if(na.rm)
				stop("na.rm not implemented on the GPU")
			ret=.Call("gpu_sum",x@ptr, length(x), x@type)
			if(retgpu)
				ret=as.gvector(ret)
			return(ret)
		}
)
setMethod("mean",  "gvector",
		function(x, na.rm = FALSE, trim=0,retgpu=TRUE) {
			if(trim != 0)
				stop("trim not supported on GPU")
			if(na.rm)
				stop("na.rm not implemented on the GPU")
			checkDevice(x@device)
			ret=.Call("gpu_sum",x@ptr, length(x), x@type)/length(x)
			if(retgpu)
				ret=as.gvector(ret)
			return(ret)
		}
)
setMethod("sum",  "gmatrix", function(x,...)
			return(sum(as.gvector(x,dup=FALSE),...)))
setMethod("mean",  "gmatrix", function(x,...)
			return(mean(as.gvector(x,dup=FALSE),...)))


setMethod("max",  "gvector",
		function(x, na.rm = FALSE, retgpu=TRUE, pos=FALSE) {
			if(na.rm)
				stop("na.rm not implemented on the GPU")
			checkDevice(x@device)
			ret=.Call("gpu_max",x@ptr, length(x), x@type)
			if(retgpu)
				ret=as.gvector(ret)
			return(ret)
		}
)
setMethod("max",  "gmatrix", function(x,...)
			return(max(as.gvector(x,dup=FALSE),...)))

setMethod("min",  "gvector",
		function(x, na.rm = FALSE, retgpu=TRUE, pos=FALSE) {
			if(na.rm)
				stop("na.rm not implemented on the GPU")
			checkDevice(x@device)
			ret=.Call("gpu_min",x@ptr, length(x), x@type)
			if(retgpu)
				ret=as.gvector(ret)
			return(ret)

		}
)
setMethod("min",  "gmatrix", function(x,...)
			return(max(as.gvector(x,dup=FALSE),...)))


setMethod("sort",  "gvector",
		function(x, decreasing = FALSE, dup=TRUE) {
			if(length(names(x))>0) {
				warning("Names are being striped for sorting")
				x=gnamestrip(x)
			}
			#browser()
			if(dup)
				x=gdup(x)
			#if(x@type==0 && stable)
			#	warning("There may be problems with stable sorting for type 'double'. We suggest that the user set stable to 'FALSE' when sorting doubles. This is likely a problem with the thrust library which will be fixed in future releases.")
			#SEXP gpu_sort(SEXP A_in, SEXP n_in, SEXP stable_in, SEXP decreasing_in, SEXP in_type)
			checkDevice(x@device)
			
			if(x@type ==0L)
				ret=.Call("gpu_sort",x@ptr, length(x),FALSE, as.logical(decreasing), x@type)
			else if (x@type ==1L)
				ret=.Call("gpu_sort",x@ptr, length(x),TRUE, as.logical(decreasing), x@type)
			else {
				#ok... jerry rig this just so it works... this could use some attention so it isn't so strange
				outype=x@type
				type(x)=0L
				ret=.Call("gpu_sort",x@ptr, length(x),FALSE, as.logical(decreasing), x@type)
				type(x)=outype
			}
			return(x)
		}
)

#setMethod("order",  "gvector",

#suppressWarnings(setGeneric("order",
#				function(x,...)
#					base::order(x,...)
#		))

gorder = function(x,  decreasing = FALSE, stable=TRUE, sortx=FALSE) {
	if(class(x) != "gvector")
		x=as.gvector(x)
	else
		if(!sortx)
			x=gdup(x)
	#SEXP gpu_sort(SEXP A_in, SEXP n_in, SEXP stable_in, SEXP decreasing_in, SEXP in_type)
	#warning("figure out na last")
	checkDevice(x@device)
	ret=.Call("gpu_order",x@ptr, length(x),as.logical(stable), as.logical(decreasing), x@type)
	return(ret)
}


setMethod("which",  "gvector",
		function(x, arr.ind = FALSE, useNames = TRUE) {
			if(arr.ind)
				stop("Error in which: 'arr.ind' = TRUE has not been implemented on the GPU. Transfer to CPU first.")
			#SEXP gpu_sort(SEXP A_in, SEXP n_in, SEXP stable_in, SEXP decreasing_in, SEXP in_type)
			#gpu_which(SEXP A_in, SEXP n_in)
			if(x@type!=3L)
				stop("Error in which: argument to 'which' is not logical")
			checkDevice(x@device)
			ret=.Call("gpu_which",x@ptr, length(x))
			return(ret)
		}
)


setMethod("ifelse",  "gvector",
		function(test, yes, no) {
			if(!(class(yes) %in% c("gvector","gmatrix")))
				yes=g(yes)
			if(!(class(no) %in% c("gvector","gmatrix")))
				no=g(no)
			if(test@type!=3L)
				type(test)=3L
			if(yes@type!=no@type) {
				totype=min(c(yes@type,no@type))
				if(yes@type!=totype)
					yes=convertType(yes,totype)
				if(no@type!=totype)
					no=convertType(no,totype)
			} 

			#SEXP gpu_if(SEXP H_in, SEXP A_in, SEXP B_in,SEXP snh, SEXP sna, SEXP snb, SEXP in_type)
			checkDevice(c(test@device,yes@device,no@device))
			ret=new("gvector",
					length=length(test),
					names=names(test),
					type=yes@type,
					ptr=.Call("gpu_if",test@ptr,yes@ptr, no@ptr, length(test), length(yes),length(no), yes@type ))
			return(ret)
		}
)


setMethod("ifelse",  "gmatrix",
		function(test, yes, no) {
			if(!(class(yes) %in% c("gvector","gmatrix")))
				yes=g(yes)
			if(!(class(no) %in% c("gvector","gmatrix")))
				no=g(no)
			if(test@type!=3L)
				type(test)=3L
			if(yes@type!=no@type) {
				totype=min(c(yes@type,no@type))
				if(yes@type!=totype)
					yes=convertType(yes,totype)
				if(no@type!=totype)
					no=convertType(no,totype)
			} 
			
			#SEXP gpu_if(SEXP H_in, SEXP A_in, SEXP B_in,SEXP snh, SEXP sna, SEXP snb, SEXP in_type)
			checkDevice(c(test@device,yes@device,no@device))
			ret=new("gmatrix",
					nrow = test@nrow,
					ncol = test@ncol,
					rownames=test@rownames,
					colnames=test@rownames,
					type=yes@type,
					ptr=.Call("gpu_if",test@ptr,yes@ptr, no@ptr, length(test), length(yes),length(no), yes@type ))
			return(ret)
		}
)

#setMethod("min",  "gvector",
#		function(x, na.rm = FALSE, retgpu=TRUE, pos=FALSE) {
#			if(na.rm)
#				stop("na.rm not implemented on the GPU")
#			ret=.Call("gpu_min_pos",x@ptr, length(x))
#			if(!pos) {
#				ret=x[ret]
#				if(!retgpu)
#					ret=as.numeric(ret)
#			} else if(retgpu)
#				ret=as.gvector(ret)
#			
#			return(ret)
#		}
#)
#setMethod("min",  "gmatrix", function(x,...)
#			return(min(as.gvector(x,dup=FALSE),...)))


.exprs_Av=eval(substitute(substitute(e, list(x = bquote(A), y=bquote(v))), list(e = .exprs_xy)))
gmatTimesDiagVec=function(A,v) {
	if(class(A)!="gmatrix") {
		if(is.matrix(A))
			A=as.gmatrix(A, dup=FALSE)
		else
			stop("Input 'A' must be a 'matrix' or 'gmatrix'.")
	}
	if(class(v)!="gvector")
		v=as.gvector(v, dup=FALSE)
	if(ncol(A)!=length(v))
		stop("Mismatched dimensions.")
	eval(.exprs_Av)
	ret=A
	ret@ptr=.Call("gpu_mat_times_diag_vec",A@ptr, v@ptr, nrow(A), ncol(A), ret@type)
	return(ret)
}

.exprs_AB=eval(substitute(substitute(e, list(x = bquote(A), y=bquote(B))), list(e = .exprs_xy)))
gkroneckerProd=function(A,B) {
	if(class(A)=="gvector")
		A=gmatrix(A, dup=FALSE)
	if(class(B)=="gvector")
		B=gmatrix(B, dup=FALSE)
	if(class(A)!="gmatrix") {
		if(is.matrix(A)|| is.vector(B)) {
			if(class(B)%in%c("gmatrix","gvector")) {
				if(B@type==1L)
					A=as.gmatrix(A, type=1L)
				else
					A=as.gmatrix(A)
			} else
				A=as.gmatrix(A)
		}else
			stop("Input 'A' must be a 'vector', 'matrix', 'gmatrix' or 'gvector'.")
	}
	if(class(B)!="gmatrix") {
		if(is.matrix(B) || is.vector(B))
			if(class(A)%in%c("gmatrix","gvector")){
				if(A@type==1L) 
					B=as.gmatrix(B, type=1L)
				else
					B=as.gmatrix(B)
			} else
				B=as.gmatrix(B)
		else
			stop("Input 'B' must be a 'vector', 'matrix', 'gmatrix' or 'gvector'.")
	}
	.exprs_AB
#SEXP gpu_kronecker(SEXP A_in, SEXP B_in,SEXP n_A_row_in,SEXP n_A_col_in, SEXP n_B_row_in,SEXP n_B_col_in);
	ret = new("gmatrix",ptr=.Call("gpu_kronecker",A@ptr, B@ptr, nrow(A), ncol(A), nrow(B), ncol(B), A@type),
			nrow=nrow(A)*nrow(B),ncol=ncol(A)*ncol(B), type=A@type)
	return(ret)
}

gsumby=function(v,startPos,stopPos) {
	if(class(v)!="gvector")
		v=as.gvector(v, dup=FALSE)
	if(type(v)==3L)
		type(v)=2L
	if(!is.integer(startPos))
		startPos=as.integer(startPos)
	if(!is.integer(stopPos))
		stopPos=as.integer(stopPos)
	if(length(startPos)!=length(stopPos))
		stop("The length of startPos and stopPos must be the same.")
	#SEXP gpu_kernal_sumby(SEXP A_in, SEXP index1_in,SEXP index2_in,SEXP n_A_in,SEXP n_index_in);
	checkDevice(v@device) #Todo: make the C code utilize startPos of class gvector
	ret = new("gvector",
			ptr    = .Call("gpu_kernal_sumby", v@ptr, startPos, stopPos, length(v), length(stopPos), v@type),
			length = length(startPos),
			names  = names(startPos),
			type=v@type)
	return(ret)
}


setMethod("%x%", signature(X = "gvector", Y = "gvector"), function(X,Y)	return(as.gvector(gkroneckerProd(X,Y),dup=FALSE)))
setMethod("%x%", signature(X = "gvector", Y = "numeric"), function(X,Y)	return(as.gvector(gkroneckerProd(X,Y),dup=FALSE)) )
setMethod("%x%", signature(X = "numeric", Y = "gvector"),  function(X,Y)	return(as.gvector(gkroneckerProd(X,Y),dup=FALSE)))

setMethod("%x%", signature(X = "gmatrix", Y = "gmatrix"), function(X,Y)	return(gkroneckerProd(X,Y)))
setMethod("%x%", signature(X = "gmatrix", Y = "gvector"), function(X,Y)	return(gkroneckerProd(X,Y)))
setMethod("%x%", signature(X = "gvector", Y = "gmatrix"), function(X,Y)	return(gkroneckerProd(X,Y)))

setMethod("%x%", signature(X = "gmatrix", Y = "numeric"),function(X,Y)	return(gkroneckerProd(X,Y)))
setMethod("%x%", signature(X = "numeric", Y = "gmatrix"), function(X,Y)	return(gkroneckerProd(X,Y)))

setMethod("%x%", signature(X = "gmatrix", Y = "matrix"),function(X,Y)	return(gkroneckerProd(X,Y)))
setMethod("%x%", signature(X = "matrix", Y = "gmatrix"), function(X,Y)	return(gkroneckerProd(X,Y)))

setMethod("%x%", signature(X = "gvector", Y = "matrix"),function(X,Y)	return(gkroneckerProd(X,Y)))
setMethod("%x%", signature(X = "matrix", Y = "gvector"), function(X,Y)	return(gkroneckerProd(X,Y)))

####################################### 
# elementwise special functions
####################################### 
setMethod("!","gvector", function(x) return(x==0L) )
setMethod("!","gmatrix", function(x)  return(x==0L) )


setMethod("-", signature(e1 = "gmatrix", e2 = "missing"), function(e1) return(-1*e1))
setMethod("-", signature(e1 = "gvector", e2 = "missing"), function(e1) return(-1*e1))

setMethod("+", signature(e1 = "gmatrix", e2 = "missing"), function(e1) return(gdup(e1)))
setMethod("+", signature(e1 = "gvector", e2 = "missing"), function(e1) return(gdup(e1)))


#one_over=function(x)
#	return(1/x)
#setGeneric("one_over")

.str_unaryops ='
setMethod("MRNAME","gvector", function(x) {
			if(x@type>1L)
				x=convertType(x,0L)
			if(length(gvector)<1L)
				return(x)
			checkDevice(x@device)
			return(new("gvector",ptr=.Call("gpu_MNAME",x@ptr, length(x), x@type), 
							length=x@length, names=x@names, type=MTYPE))}) 
setMethod("MRNAME","gmatrix", function(x)  {
			if(x@type>1L)
				x=convertType(x,0L)
			checkDevice(x@device)
			return(new("gmatrix",ptr=.Call("gpu_MNAME",x@ptr, as.integer(prod(dim(x))),x@type), 
							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames, type=MTYPE))})'
		
					.unOps=		list(
						#	c('one_over','one_over', "x@type"), 
							c('sqrt','sqrt',  "x@type"),
							c('exp','exp',  "x@type"),
							c('expm1','expm1',  "x@type"),
							c('log','log',  "x@type"),
							c('log2','log2',  'x@type'),
							c('log10','log10',  'x@type'),
							c('log1p','log1p',  'x@type'),
							c('sin','sin',  'x@type'),
							c('cos','cos',  'x@type'),
							c('tan','tan',  'x@type'),
							c('asin','asin',  'x@type'),
							c('acos','acos',  'x@type'),
							c('atan','atan',  'x@type'),
							c('sinh','sinh',  'x@type'),
							c('cosh','cosh',  'x@type'),
							c('tanh','tanh',  'x@type'),
							c('asinh','asinh',' x@type'),
							c('acosh','acosh', 'x@type'),
							c('atanh','atanh',  'x@type'),
							c('abs','fabs',  'x@type'),
							c('lgamma','lgamma',  'x@type'),
							c('gamma','gamma',  'x@type'),
							c('sign','sign',  'x@type'),
							c('round','round', '2L'),
							c('ceiling','ceil',  '2L'),
							c('floor','floor',  '2L'),
							c('is.na','isna',  '3L'),
							c('is.nan','isnan',  '3L'),
							c('is.finite','isfinite',  '3L'),
							c('is.infinite','isinfinite',  '3L')					
					)

lapply(.unOps, function(op) {
								exprStr = gsub("MRNAME",op[1],.str_unaryops)
								exprStr = gsub("MNAME",op[2],exprStr)
								exprStr = gsub("MTYPE",op[3],exprStr)
								#cat(exprStr)
								eval(parse(text=exprStr))
							} )				

					
					
####################################### 
# elementwise binary operations
####################################### 


#Non-exchangeable binary operations
.strOpNonExchang = '
setMethod("MOP", signature(e1 = "gvector", e2 = "gvector"), 
		function(e1,e2) {
			MEXPRS
			if(length(e2)<1L)
				return(e2)
			else if(length(e1)<1L)
				return(e1)
			else if(length(e2)==1L){
				tryCatch(e2 <- .convert_to_appropriate_class(e2,e1@type), error=function(e) stop("Invalid value for operator: MOP.", call. = FALSE))
				#e2=as.double(e2)
				ret=new("gvector", length=length(e1), names=names(e1), type=MT1)
				ret@ptr=.Call("gpu_scaler_MNAME12", e1@ptr, e2, length(e1), e1@type)
			} else if(length(e1)==1L){
				tryCatch(e1 <- .convert_to_appropriate_class(e1,e2@type), error=function(e) stop("Invalid value for operator: MOP.", call. = FALSE))
				#e1=as.double(e1)
				ret=new("gvector", length=length(e2), names=names(e2), type=MT2)
				ret@ptr=.Call("gpu_scaler_MNAME21", e2@ptr, e1, length(e2), e2@type)	
			} else if(length(e1)==length(e2)) {
				ret=new("gvector",length=length(e1), names=names(e1), type=MT1)
				ret@ptr=.Call("gpu_same_size_MNAME12", e1@ptr, e2@ptr,  e2@length, e1@type)
			} else {
				if(length(e1)<length(e2)) {
					small=e1
					big=e2
					reverse=TRUE
				} else {
					small=e2
					big=e1
					reverse=FALSE
				}
				if((length(big) %% length(small)) !=0)
					stop("longer object length is not a multiple of shorter object length", call. = FALSE)
				
				ret=new("gvector",length=length(big), names=names(big), type=MT1)

				if(reverse)	
					ret@ptr=.Call("gpu_diff_size_MNAME21", big@ptr, small@ptr, big@length, small@length, e1@type)
				else
					ret@ptr=.Call("gpu_diff_size_MNAME12", big@ptr, small@ptr, big@length, small@length, e1@type)
			} 
			return(ret)
		}
)
.gvec.vec.MNAME = function(e1,e2) {
			if(length(e2)<1L)
				return(e2)
			else if(length(e1)<1L)
				return(e1)
			else if(length(e2)==1L){
				#tryCatch(e2 <- .convert_to_appropriate_class(e2,e1@type), error=function(e) stop("Invalid value for operator: MOP.", call. = FALSE))
				#e2=as.double(e2)
				MDIFEXPRS12
				ret=new("gvector", length=length(e1), names=names(e1), type=MT1)
				ret@ptr=.Call("gpu_scaler_MNAME12", e1@ptr, e2, length(e1), e1@type)
			} else {
				if(e1@type==1L)
					tmp=as.gvector(e2, type=1L)
				else
					tmp=as.gvector(e2)
				ret=e1 MOP tmp	
			}
			return(ret)
		}
.vec.gvec.MNAME = 	function(e1,e2) {		
			if(length(e2)<1L)
				return(e2)
			else if(length(e1)<1L)
				return(e1)
			else if(length(e1)==1L){
				#tryCatch(e1 <- .convert_to_appropriate_class(e1,e2@type), error=function(e) stop("Invalid value for operator: MOP.", call. = FALSE))
				#e1=as.double(e1)
				MDIFEXPRS21
				ret=new("gvector", length=length(e2), names=names(e2), type=MT2)
				ret@ptr=.Call("gpu_scaler_MNAME21", e2@ptr, e1, length(e2), e2@type)	
			} else {
				if(e2@type==1L)
					tmp=as.gvector(e1, type=1L)
				else
					tmp=as.gvector(e1)
				ret=tmp MOP e2	
			}
			return(ret)
		}

setMethod("MOP", signature(e1 = "gvector", e2 = "numeric"), .gvec.vec.MNAME )	
setMethod("MOP", signature(e1 = "gvector", e2 = "logical"), .gvec.vec.MNAME )	
setMethod("MOP", signature(e1 = "numeric", e2 = "gvector"), .vec.gvec.MNAME )		
setMethod("MOP", signature(e1 = "logical", e2 = "gvector"), .vec.gvec.MNAME )

setMethod("MOP", signature(e1 = "gmatrix", e2 = "gmatrix"), 
		function(e1,e2) {
			MEXPRS
			if(any(dim(e1)!=dim(e2)))
				stop("Dimensions of matrix do not match for elementwise operations.", call. = FALSE)
			myn=as.integer(prod(dim(e1)))
			ret=new("gmatrix", nrow = e1@nrow, ncol = e1@ncol,
					rownames=e1@rownames, colnames=e1@colnames, type=MT1)
			ret@ptr=.Call("gpu_same_size_MNAME12", e1@ptr, e2@ptr, myn, e1@type)
			return(ret)
		})
.gmat.vec.MNAME=function(e1,e2) {
			checkDevice(e1@device)
			tmp1=new("gvector",length=as.integer(prod(dim(e1))),ptr=e1@ptr, type=e1@type)
			if(length(tmp1)==1) {
				ret = tmp1 MOP e2
			}else if(length(tmp1) < length(e2)) {
				stop("Dimensions do not match for elementwise operations.", call. = FALSE)
			}else {
				ret=e1
				tmpret=tmp1 MOP e2
				ret@ptr=tmpret@ptr
				ret@type=tmpret@type
			}
			return(ret)
		}
.vec.gmat.MNAME =function(e1,e2) {
			checkDevice(e2@device)
			tmp2=new("gvector",length=as.integer(prod(dim(e2))),ptr=e2@ptr, type=e2@type)
			if(length(tmp2)==1) {
				ret = e1 MOP tmp2
			}else if(length(tmp2) < length(e1)) {
				stop("Dimensions do not match for elementwise operations.", call. = FALSE)
			}else {
				ret=e2
				tmpret=e1 MOP tmp2
				ret@ptr=tmpret@ptr
				ret@type=tmpret@type
			}
			return(ret)
		}

setMethod("MOP", signature(e1 = "gmatrix", e2 = "numeric"), .gmat.vec.MNAME)
setMethod("MOP", signature(e1 = "gmatrix", e2 = "logical"), .gmat.vec.MNAME)
setMethod("MOP", signature(e1 = "numeric", e2 = "gmatrix"), .vec.gmat.MNAME)
setMethod("MOP", signature(e1 = "logical", e2 = "gmatrix"), .vec.gmat.MNAME)

setMethod("MOP", signature(e1 = "gmatrix", e2 = "gvector"), 
		function(e1,e2) {
			checkDevice(c(e1@device,e2@device))
			tmp1=new("gvector",length=as.integer(prod(dim(e1))),ptr=e1@ptr, type=e1@type)
			if(length(tmp1)==1) {
				ret = tmp1 MOP e2
			}else if(length(tmp1) < length(e2)) {
				stop("Dimensions do not match for elementwise operations.", call. = FALSE)
			}else {
				ret=e1
				tmpret=tmp1 MOP e2
				ret@ptr=tmpret@ptr
				ret@type=tmpret@type
			}
			return(ret)
		})

setMethod("MOP", signature(e1 = "gvector", e2 = "gmatrix"), 
		function(e1,e2) {
			checkDevice(c(e1@device,e2@device))
			MEXPRS
			tmp2=new("gvector",length=as.integer(prod(dim(e2))),ptr=e2@ptr, type= e2@type)
			if(length(tmp2)==1) {
				ret = e1 MOP tmp2
			}else if(length(tmp2) < length(e1)) {
				stop("Dimensions do not match for elementwise operations.", call. = FALSE)
			}else {
				ret=e2
				tmpret=e1 MOP tmp2
				ret@ptr=tmpret@ptr
				ret@type=tmpret@type
			}
			return(ret)
		}
)

setMethod("MOP", signature(e1 = "gvector", e2 = "matrix"), 
		function(e1,e2) {
			if(e1@type==1L)
				tmp=as.gmatrix(e2,type=1L)
			else
				tmp=as.gmatrix(e2)
			return(e1 MOP tmp)
		}
)
setMethod("MOP", signature(e1 = "matrix", e2 = "gvector"), 
		function(e1,e2) {
			if(e2@type==1L)
				tmp=as.gmatrix(e1,type=1L)
			else
				tmp=as.gmatrix(e1)
			return(tmp MOP e2)
		}
)

setMethod("MOP", signature(e1 = "gmatrix", e2 = "matrix"), 
		function(e1,e2) {
			if(any(dim(e1)!=dim(e2)))
				stop("Dimensions of matrix do not match for elementwise operation.", call. = FALSE)
			if(e1@type==1L)
				return(e1 MOP as.gmatrix(e2, type=1L))
			else
				return(e1 MOP as.gmatrix(e2))
		})
setMethod("MOP", signature(e1 = "matrix", e2 = "gmatrix"), 
		function(e1,e2) {
			if(any(dim(e1)!=dim(e2)))
				stop("Dimensions of matrix do not match for elementwise operation.", call. = FALSE)
			if(e2@type==1L)
				return(as.gmatrix(e1, type=1L) MOP e2)
			else
				return(as.gmatrix(e1) MOP e2)
		})
'

#Exchangeable binary operations
.strOpExchang = '
setMethod("MOP", signature(e1 = "gvector", e2 = "gvector"), 
		function(e1,e2) {
			MEXPRS
			if(length(e2)<1L)
				return(e2)
			else if(length(e1)<1L)
				return(e1)
			else if(length(e1)<length(e2)) {
				small=e1
				big=e2
			} else {
				small=e2
				big=e1
			}
			if(length(small)==1L){
				tryCatch(small <- .convert_to_appropriate_class(small,big@type), error=function(e) stop("Invalid value for operator: MOP.", call. = FALSE))
				#small=as.double(small)
				ret=new("gvector", length=length(big), type = MT2)
				ret@ptr=.Call("gpu_scaler_MNAME", big@ptr, small, length(big), big@type)
			} else if(length(big)==length(small)) {
				ret=new("gvector",length=length(big), type=MT1)
				ret@ptr=.Call("gpu_same_size_MNAME", big@ptr, small@ptr,  small@length, big@type)
			} else if((length(big) %% length(small)) ==0) {
				ret=new("gvector",length=length(big), type=MT1)
				ret@ptr=.Call("gpu_diff_size_MNAME", big@ptr, small@ptr, big@length, small@length, big@type)
			} else
				stop("longer object length is not a multiple of shorter object length", call. = FALSE)
			return(ret)
		}
)
.gvec.vec.MNAME =
		function(e1,e2) {
			if(length(e2)<1L)
				return(e2)
			else if(length(e1)<1L)
				return(e1)
			else if(length(e2)==1L){
				MDIFEXPRS12
				ret=new("gvector", length=length(e1), type=MT1)
				#tryCatch(e2 <- .convert_to_appropriate_class(e2,e1@type), error=function(e) stop("Invalid value for operator: MOP.", call. = FALSE))
				ret@ptr=.Call("gpu_scaler_MNAME", e1@ptr, e2, length(e1), e1@type)
			} else {
				#if(length(e1)<length(e2)) 
				#	stop("for elementwise operations, the larger vector/matrix must be on the gpu", call. = FALSE)
				if(e1@type==1L)
					tmp=as.gvector(e2, type=1L)
				else
					tmp=as.gvector(e2)
				ret=e1 MOP tmp
			}
			return(ret)
		}
setMethod("MOP", signature(e1 = "gvector", e2 = "numeric"), .gvec.vec.MNAME )
setMethod("MOP", signature(e1 = "gvector", e2 = "logical"), .gvec.vec.MNAME )
setMethod("MOP", signature(e1 = "numeric", e2 = "gvector"), function(e1,e2)  .gvec.vec.MNAME(e2, e1))
setMethod("MOP", signature(e1 = "logical", e2 = "gvector"), function(e1,e2)  .gvec.vec.MNAME(e2, e1))


setMethod("MOP", signature(e1 = "gmatrix", e2 = "gmatrix"), 
		function(e1,e2) {
			#browser()
			if(any(dim(e1)!=dim(e2)))
				stop("Dimensions of matrix do not match for elementwise operation.", call. = FALSE)
			MEXPRS
			myn=as.integer(prod(dim(e1)))
			ret=new("gmatrix", nrow = e1@nrow, ncol = e1@ncol,rownames=e1@rownames,colnames=e1@colnames, type=MT1)
			ret@ptr=.Call("gpu_same_size_MNAME", e1@ptr, e2@ptr, myn,  e1@type)
			#browser()
			return(ret)
		})

.gmat.vec.MNAME = 	function(e1,e2) {
			checkDevice(e1@device)
			tmp1=new("gvector",length=as.integer(prod(dim(e1))),ptr=e1@ptr, type=e1@type)
			if(length(tmp1)==1) {
				ret=tmp1  MOP  e2
			}else if(length(tmp1)< length(e2)) {
				stop("Length of the matrix does not match for elementwise multiplication", call. = FALSE)
			}else {
				ret=e1
				tmpret=tmp1  MOP  e2
				ret@type=tmpret@type
				ret@ptr=tmpret@ptr
			}
			return(ret)
		}
setMethod("MOP", signature(e1 = "gmatrix", e2 = "numeric"), .gmat.vec.MNAME)
setMethod("MOP", signature(e1 = "gmatrix", e2 = "logical"), .gmat.vec.MNAME)
setMethod("MOP", signature(e1 = "numeric", e2 = "gmatrix"), function(e1,e2) .gmat.vec.MNAME(e2,e1))
setMethod("MOP", signature(e1 = "logical", e2 = "gmatrix"), function(e1,e2) .gmat.vec.MNAME(e2,e1))
setMethod("MOP", signature(e1 = "gmatrix", e2 = "gvector"), 
		function(e1,e2) {
			checkDevice(c(e1@device,e2@device))
			tmp1=new("gvector",length=as.integer(prod(dim(e1))),ptr=e1@ptr, type=e1@type)
			if(length(tmp1)==1) {
				ret=tmp1  MOP  e2
			}else if(length(tmp1)< length(e2)) {
				stop("Length of the matrix does not match for elementwise multiplication", call. = FALSE)
			}else {
				ret=e1
				tmpret=tmp1  MOP  e2
				ret@type=tmpret@type
				ret@ptr=tmpret@ptr
			}
			return(ret)
		})
setMethod("MOP", signature(e1 = "gvector", e2 = "gmatrix"), function(e1,e2) return(e2 MOP e1))
setMethod("MOP", signature(e1 = "gmatrix", e2 = "matrix"), 
		function(e1,e2) {
			if(any(dim(e1)!=dim(e2)))
				stop("Dimensions of matrix do not match for elementwise operation.", call. = FALSE)
			if(e1@type==1L)
				return(e1 MOP as.gmatrix(e2, type=1L))
			else
				return(e1 MOP as.gmatrix(e2))
			
		})
setMethod("MOP", signature(e1 = "matrix", e2 = "gmatrix"), 
		function(e1,e2) {
			if(any(dim(e1)!=dim(e2)))
				stop("Dimensions of matrix do not match for elementwise operation.", call. = FALSE)
			if(e2@type==1L)
				return(as.gmatrix(e1, type=1L) MOP e2)
			else
				return(as.gmatrix(e1) MOP e2)
		})
setMethod("MOP", signature(e1 = "gvector", e2 = "matrix"), 
		function(e1,e2) {
			if(e1@type==1L)
				tmp=as.gmatrix(e2,type=1L)
			else
				tmp=as.gmatrix(e2)
			return(e1 MOP tmp)
		}
)
setMethod("MOP", signature(e1 = "matrix", e2 = "gvector"), 
		function(e1,e2) {
			if(e2@type==1L)
				tmp=as.gmatrix(e1,type=1L)
			else
				tmp=as.gmatrix(e1)
			return(tmp MOP e2)
		}
)

'

.exprs_e1e2=eval(substitute(substitute(e, list(x = bquote(e1), y=bquote(e2))), list(e = .exprs_xy)))
.exprs_str = paste( deparse(.exprs_e1e2),collapse="\n")
.exprs_sf_e1e2=eval(substitute(substitute(e, list(x = bquote(e1), y=bquote(e2))), list(e = .exprs_sf_xy )))
.exprs_sf_str =  paste( deparse(.exprs_sf_e1e2),collapse="\n")
.exprs_l_e1e2=eval(substitute(substitute(e, list(x = bquote(e1), y=bquote(e2))), list(e = .exprs_l_xy  )))
.exprs_l_str =  paste( deparse(.exprs_l_e1e2),collapse="\n")
.exprs_gc_e1e2=eval(substitute(substitute(e, list(x = bquote(e1), y=bquote(e2))), list(e = .exprs_gpu_cpu_xy )))
.exprs_gc12_str =  paste( deparse(.exprs_gc_e1e2),collapse="\n")
.exprs_gc_e2e1=eval(substitute(substitute(e, list(x = bquote(e2), y=bquote(e1))), list(e = .exprs_gpu_cpu_xy )))
.exprs_gc21_str =  paste( deparse(.exprs_gc_e2e1),collapse="\n")

.exprs_l_gc_e1e2=eval(substitute(substitute(e, list(x = bquote(e1), y=bquote(e2))), list(e = .exprs_l_gpu_cpu_xy )))
.exprs_l_gc12_str =  paste( deparse(.exprs_l_gc_e1e2),collapse="\n")
.exprs_l_gc_e2e1=eval(substitute(substitute(e, list(x = bquote(e2), y=bquote(e1))), list(e = .exprs_l_gpu_cpu_xy )))
.exprs_l_gc21_str =  paste( deparse(.exprs_l_gc_e2e1),collapse="\n")

.exprs_compare_e1e2=eval(substitute(substitute(e, list(x = bquote(e1), y=bquote(e2))), list(e = .exprs_compare_xy )))
.exprs_compare_str =  paste( deparse(.exprs_compare_e1e2),collapse="\n")

.exprs_sf_gc_e1e2=eval(substitute(substitute(e, list(x = bquote(e1), y=bquote(e2))), list(e = .exprs_sf_gpu_cpu_xy )))
.exprs_sf_gc12_str =  paste( deparse(.exprs_sf_gc_e1e2),collapse="\n")
.exprs_sf_gc_e2e1=eval(substitute(substitute(e, list(x = bquote(e2), y=bquote(e1))), list(e = .exprs_sf_gpu_cpu_xy )))
.exprs_sf_gc21_str =  paste( deparse(.exprs_sf_gc_e2e1),collapse="\n")

setGeneric("%lgspadd%",
		function(e1,e2)
			standardGeneric("%lgspadd%")
)


.exOps=list(
		c("lgspadd","%lgspadd%", .exprs_str, "e1@type", "big@type", .exprs_gc12_str,.exprs_gc21_str  ),
		c("mult","*",      .exprs_str, "e1@type", "big@type", .exprs_gc12_str,.exprs_gc21_str  ),
		c("add" ,"+",      .exprs_str, "e1@type", "big@type", .exprs_gc12_str,.exprs_gc21_str  ),
		c("eq"  ,"==",     .exprs_compare_str, "3L", "3L", .exprs_gc12_str,    .exprs_gc21_str  ),
		c("ne"  ,"!=",     .exprs_compare_str, "3L", "3L" , .exprs_gc12_str,   .exprs_gc21_str  ),
		c("and"  ,"&",    .exprs_l_str, "3L", "3L", .exprs_l_gc12_str, .exprs_l_gc21_str  ),
		c("or"  ,"|",     .exprs_l_str, "3L", "3L", .exprs_l_gc12_str, .exprs_l_gc21_str  )
		)


		lapply(.exOps, function(op) {
					#browser()
					exprStr = 	gsub("MDIFEXPRS21",op[7],
							gsub("MDIFEXPRS12",op[6],
									gsub("MT2",op[5],
											gsub("MT1",op[4],
													gsub("MEXPRS",op[3],
															gsub("MNAME", op[1],
																	gsub("MOP",op[2],.strOpExchang)))))))
					#cat(exprStr)
					eval(parse(text=exprStr))
				} )


.nonExOps=list(	c("sub","-",.exprs_str, "e1@type", "e2@type",.exprs_gc12_str,.exprs_gc21_str ),
		c("div","/",        .exprs_sf_str, "e1@type", "e2@type",.exprs_sf_gc12_str,.exprs_sf_gc21_str ),
		c("pow","^",        .exprs_sf_str, "e1@type", "e2@type",.exprs_sf_gc12_str,.exprs_sf_gc21_str ),
		c("mod","%%",       .exprs_sf_str, "e1@type", "e2@type",.exprs_sf_gc12_str,.exprs_sf_gc21_str ),
		c("gt"  ,">",      .exprs_compare_str, "3L", "3L", .exprs_gc12_str,.exprs_gc21_str  ),
		c("lt"  ,"<",      .exprs_compare_str, "3L", "3L", .exprs_gc12_str,.exprs_gc21_str  ),
		c("gte" ,">=",     .exprs_compare_str, "3L", "3L", .exprs_gc12_str,.exprs_gc21_str  ),
		c("lte" ,"<=",     .exprs_compare_str, "3L", "3L", .exprs_gc12_str,.exprs_gc21_str  )
		)

lapply(.nonExOps, function(op) {
			exprStr = gsub("MDIFEXPRS21",op[7],
					gsub("MDIFEXPRS12",op[6],
							gsub("MT2",op[5],
									gsub("MT1",op[4],
											gsub("MEXPRS",op[3],
													gsub("MNAME", op[1],
															gsub("MOP",op[2],.strOpNonExchang)))))))
			eval(parse(text=exprStr))
		} )



#A= matrix(1:10,ncol=2)
#B= matrix(1:10,nrow=2)
#A=as.gmatrix(A)
#B=as.gmatrix(B)
#C=A%*%B
#D=A%*%B
#E=C*D
#
#A<-matrix(1:9000000,ncol=2)
#system.time(
#A<-A*(1:2)
#)
#
#A<-gmatrix(1:9000000,ncol=2)
#gc()
#system.time(B<-A*A)
#system.time(C<-A*1:2)
#system.time(D<-A*3)

#
#setMethod("one_over","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_one_over",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("sqrt","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_sqrt",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("sqrt","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_sqrt",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("exp","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_exp",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("exp","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_exp",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("expm1","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_expm1",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("expm1","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_expm1",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("log","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_log",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("log","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_log",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("log2","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_log2",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("log2","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_log2",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("log10","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_log10",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("log10","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_log10",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("log1p","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_log1p",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("log1p","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_log1p",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("sin","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_sin",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("sin","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_sin",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("cos","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_cos",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("cos","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_cos",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("tan","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_tan",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("tan","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_tan",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("asin","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_asin",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("asin","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_asin",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("acos","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_acos",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("acos","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_acos",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("atan","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_atan",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("atan","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_atan",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("sinh","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_sinh",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("sinh","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_sinh",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("cosh","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_cosh",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("cosh","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_cosh",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("tanh","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_tanh",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("tanh","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_tanh",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("asinh","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_asinh",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("asinh","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_asinh",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("acosh","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_acosh",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("acosh","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_acosh",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("atanh","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_atanh",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("atanh","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_atanh",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("abs","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_abs",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("abs","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_abs",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("lgamma","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_lgamma",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("lgamma","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_lgamma",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#
#setMethod("gamma","gmatrix", function(x)  
#			return(new("gmatrix",ptr=.Call("gpu_gamma",x@ptr, as.integer(prod(dim(x)))), 
#							nrow=x@nrow, ncol=x@ncol, rownames=x@rownames, colnames=x@colnames))) 
#setMethod("gamma","gvector", function(x)  
#			return(new("gvector",ptr=.Call("gpu_gamma",x@ptr, length(x)), 
#							length=x@length, names=x@names))) 
#qr(matrix(1,10,10))
#
#matrixRank=function(A) {
#	return(qr(A)$rank)
#}
#matrixRank(matrix(1,10,10))