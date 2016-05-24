###########################################################################
# Test Script
# 
# Author: nmorris
###############################################################################

setClass("gqr",
		representation(
				qr="gmatrix",
				qraux = "gvector"),
		prototype = list(
				qr=NULL,
				qraux=NULL)
) 




setGeneric("gqr.coef", function(qr,y) stop("'gqr.coef' can only act on an object of class 'gqr.'"))
setMethod("gqr.coef",  signature(qr = "gqr", y = "ANY"),
		function(qr, y) 
		{
			#browser()
			if(qr@qr@type>1L)
				type(qr@qr)=0L
			if(qr@qraux@type>1L)
				type(qr@qraux)=0L
			if(qr@qr@type!=qr@qraux@type) {
				totype=min(c(qr@qr@type,qr@qraux@type))
				type(qr@qr)=totype
				type(qr@qraux)=totype
			}

			if (class(y)!="gmatrix") 
				y <- as.gmatrix(y)
			else
				y=gdup(y)
			
			if(y@type>1L)
				type(y)=type(qr@qr)
			if(y@type!=qr@qr@type)
				stop("Type mismatch.")
			if(qr@qraux@type!=qr@qr@type)
				stop("Type mismatch.")
			
			checkDevice(c(y@device,qr@qr@device,qr@qraux@device))
			
			n <- as.integer(nrow(qr@qr))
			p <- as.integer(ncol(qr@qr))
			#k <- as.integer(qr$rank)
			ny <- as.integer(ncol(y))
			if(n!=nrow(y))
			   stop("Dimension mismatch")
			if (p == 0L) 
				stop("gmatrix in qr has a dimension of 0")

			if (p > n) 
				stop("Columns of qr matrix must be less than or equal to rows.")
				
			dummy=.Call("rcusolve_modqr_coef", qr@qr, qr@qraux, y)
			
			if (p < n)
				return(y[1:p,])
			else
				return(y)
		}
)

setMethod("solve", signature(a = "gmatrix", b = "ANY"),
		function (a, b, ...)
		{
			if(a@type>1L)
				type(a)=0L
			#if(!missing(b))
			#	if(b@type>1L)
			#		type(b)=type(a)
			#browser()
			a <- qr(a)
			nc <- ncol(a@qr)
			#if (a$rank != nc) 
			#	stop("singular matrix 'a' in 'solve'")
			if (missing(b)) {
				if (nc != nrow(a@qr)) 
					stop("only square matrices can be inverted")
				b <- gident(nc, type=a@qr@type)
				colnames(b) <- rownames(a@qr)
			}
			return(gqr.coef(a, b))
		}
)


setMethod("solve", signature(a = "gqr", b = "ANY"),
		function (a, b, ...) 
		{
			nc <- ncol(a@qr)
			nr <- nrow(a@qr)
#			if (a$rank != min(nc, nr)) 
#				if (a$rank != nc) 
#					stop("singular matrix 'a' in 'solve'")
			if (missing(b)) {
				if (nc != nr) 
					stop("only square matrices can be inverted")
				b <- gident(nc, type=a@qr@type)
				colnames(b) <- rownames(a@qr)
			}
			res <- gqr.coef(a, b)
			#res[is.na(res)] <- 0
			res
		}
)
setMethod("qr", signature(x = "gmatrix"),
			  function(x,...) 
		{	
			checkDevice(x@device)
			if(x@type>1L)
				type(x)=0L
			checkDevice(x@device)
			res=new("gqr",
					qr=gdup(x),
					qraux=gvector(ncol(x), type=x@type)
			)

			tmp <- .Call("rcusolve_qr", res@qr, res@qraux@ptr)
			#res@rank=as.integer(sum(res@qraux/max(res@qraux)>10^-6, retgpu=FALSE))
			#if (!is.null(cn <- colnames(x))) 
			#	colnames(res@qr) <- cn[res@pivot]
			if (!is.null(cn <- colnames(x))) 
				colnames(res@qr) <- cn
			if (!is.null(cn <- rownames(x))) 
				rownames(res@qr) <- cn
			return(res)
		}
)


setClass("gsvd",
		representation(
				U="gmatrix",
				S = "gvector",
				VT = "gmatrix"),
		prototype = list(
				U=NULL,
				S=NULL,
				VT=NULL)
) 

setGeneric("svd", function(x,...) base::svd(x,...))
setMethod("svd", signature(x = "gmatrix"),
		function (x) 
		{
			if(x@type>1L)
				type(x)=0L
			else
				x=gdup(x)
				
			if(nrow(x)<ncol(x))
				stop("SVD on the GPU only works for a matrix rows>=cols. Try transposing the 'x' matrix first.")
			res=new("gsvd",
					U=gmatrix(0,nrow=nrow(x), ncol=nrow(x), type=x@type),
					VT=gmatrix(0,nrow=ncol(x), ncol=ncol(x), type=x@type),
					S=gvector(min(ncol(x),nrow(x)), type=x@type)
			)
			tmp <- .Call("rcusolve_svd", x, res@S@ptr,res@U, res@VT)
			res
		}
)


setMethod("chol", signature(x = "gmatrix"),
		function (x, dup=TRUE) 
		{
			if(x@type>1L) {
				type(x)=0L
			} else {
				if(dup)
					x=gdup(x)
			}
			
				
				
			if(nrow(x)!=ncol(x))
				stop("chol on the GPU only works for a matrix rows=cols.")
				
			tmp <- .Call("rcusolve_chol", x)
			return(x)
		}
)


