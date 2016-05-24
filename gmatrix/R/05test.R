###########################################################################
# Test Script
# 
# Author: nmorris
###############################################################################
#library(gmatrix,lib.loc="/home/nmorris/.RLIBS")
gtest =function() {
	
	tn=100
	tsn=5
	scaler=10.13
	testscalergpu=list(
			as.gvector(scaler,type="d"),
			as.gvector(scaler,type="s"),
			as.gvector(scaler,type="i"),
			as.gvector(scaler,type="l")
	)
	testscalercpu=list(
			as.numeric(scaler),
			as.numeric(scaler),
			as.integer(scaler),
			as.logical(scaler))
	
	
	
	vec=((1:tn)-1)*(5/tn)
	testvecgpu=list(
			as.gvector(vec,type="d"),
			as.gvector(vec,type="s"),
			as.gvector(vec,type="i"),
			as.gvector(vec,type="l")
	)
	testveccpu=list(
			as.numeric(vec),
			as.numeric(vec),
			as.integer(vec),
			as.logical(vec))
	
	
	mat=matrix(0:(tn*tn-1),tn,tn)*(5/tn)
	testmatgpu=list(
			as.gmatrix(mat,type="d"),
			as.gmatrix(mat,type="s"),
			as.gmatrix(mat,type="i"),
			as.gmatrix(mat,type="l"))
	
	testmatcpu=list(
			matrix(as.numeric(mat),tn,tn),
			matrix(as.numeric(mat),tn,tn),
			matrix(as.integer(mat),tn,tn),
			matrix(as.logical(mat),tn,tn)
	)
	
	
	vec=5*((1:tsn)-1)/tsn
	testsvecgpu=list(
			as.gvector(vec,type="d"),
			as.gvector(vec,type="s"),
			as.gvector(vec,type="i"),
			as.gvector(vec,type="l")
	)
	testsveccpu=list(
			as.numeric(vec),
			as.numeric(vec),
			as.integer(vec),
			as.logical(vec))
	
	
	mat=matrix(0:(tsn*tsn-1),tsn,tsn)*(5/tn)
	testsmatgpu=list(
			as.gmatrix(mat,type="d"),
			as.gmatrix(mat,type="s"),
			as.gmatrix(mat,type="i"),
			as.gmatrix(mat,type="l"))
	
	testsmatcpu=list(
			matrix(as.numeric(mat),tn,tn),
			matrix(as.numeric(mat),tn,tn),
			matrix(as.integer(mat),tn,tn),
			matrix(as.logical(mat),tn,tn)
	)
	
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
	
	
	.type_name2= function(type) {
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
	
	warningslist=list()
	
	###################################
	# Test Matrix Multiplication
	###################################
	cputype=function(gpu,type) {
		ret=type
		if(!gpu)
			if(type==1)
				ret=0
		return(ret)
	}
	checkmm=function(x1,y1,mat1,mat2,gpu1,gpu2, type1=cputype(gpu1,i-1),type2=cputype(gpu2,j-1)) {
		yexpr=substitute(y1, list(i=i,j=j))
		xexpr=substitute(x1, list(i=i,j=j))
		tmpnm=deparse(expr=xexpr)
		mywarnings=character(0)
		x=tryCatch(eval(xexpr), error=function(e) {
					mywarnings<<-conditionMessage(e)
					return(-99999L)})
		if(as.logical(x[1]!=-99999L)) {
			y=eval(yexpr)
			
			if(class(x) %in% c("gmatrix","gvector")) {
				
				maxtype=max(type1,type2)
				mintype=min(type1,type2)
				if(maxtype>=1L & mintype==1L)
					out_type=1L
				else
					out_type=0L
				if(!gpu1 & type1==0L)
					if(type2==1L)
						out_type=1L
				if(!gpu2 & type2==0L)
					if(type1==1L)
						out_type=1L
				
				if(x@type!=out_type) {
					mywarnings=c(mywarnings,paste("Incorrect type returned:",.type_name(x@type)))
					#	browser()
				}
				x=h(x)
			}
			
			if(class(y) %in% c("gmatrix","gvector")) {
				y=h(y)
			}
			tmp=sum((y-x)^2)/sum(x^2)
			if(!is.numeric(tmp)||any(is.na(tmp)))
					mywarnings=c(mywarnings,paste(mywarnings, "result missing or not numeric."))
			else if(tmp>10^(-6))
					mywarnings=c(mywarnings,"not equal")
		}
		if(length(mywarnings)>0) {
			c1=ifelse(gpu1,paste(ifelse(mat1,"gmatrix","gvector"),"of type",.type_name(type1)),
					paste(ifelse(mat1,"matrix","vector"),"of type",.type_name2(type1)))
			c2=ifelse(gpu2,paste(ifelse(mat2,"gmatrix","gvector"),"of type",.type_name(type2)),
					paste(ifelse(mat2,"matrix","vector"),"of type",.type_name2(type2)))
			mywarnings=list(mywarnings)
			names(mywarnings)=paste("applying",opname," a",c1, "and a",c2,".")
			
			warningslist<<-c(warningslist,mywarnings)
		}
		return(mywarnings)
	}
	
	#warningslist=list()
	cat("Checking matrix multiplication, crossprod and tcrossprod... \n")
	multOps=list(
			list("matrix multiplication", function(a,b) a %*% b),
			list("crossprod", function(a,b) crossprod(a, b)),
			list("tcrossprod", function(a,b) tcrossprod(a , b))
	)
	for(op in multOps) {	
	#	browser()
		#suppressWarnings(setGeneric("%op%",op[[2]],where=globalenv()))
		#print(op)
		opname=op[[1]]
		opf=op[[2]]
		for(i in 1:4){
			for(j in 1:4){
				checkmm(opf(testmatgpu[[i]], testmatgpu[[j]]), opf(testmatcpu[[i]], testmatcpu[[j]]),mat1=TRUE,mat2=TRUE,gpu1=TRUE,gpu2=TRUE)
				checkmm(opf(testmatcpu[[i]], testmatgpu[[j]]), opf(testmatcpu[[i]], testmatcpu[[j]]),mat1=TRUE,mat2=TRUE,gpu1=FALSE,gpu2=TRUE)
				checkmm(opf(testmatgpu[[i]], testmatcpu[[j]]), opf(testmatcpu[[i]], testmatcpu[[j]]),mat1=TRUE,mat2=TRUE,gpu1=TRUE,gpu2=FALSE)
				
				if(opname!="tcrossprod"){
					checkmm(opf(testmatgpu[[i]], testvecgpu[[j]]), opf(testmatcpu[[i]], testveccpu[[j]]),mat1=TRUE,mat2=FALSE,gpu1=TRUE,gpu2=TRUE)
					checkmm(opf(testmatgpu[[i]], testveccpu[[j]]), opf(testmatcpu[[i]], testveccpu[[j]]),mat1=TRUE,mat2=FALSE,gpu1=TRUE,gpu2=FALSE)
					checkmm(opf(testmatcpu[[i]], testvecgpu[[j]]), opf(testmatcpu[[i]], testveccpu[[j]]),mat1=TRUE,mat2=FALSE,gpu1=FALSE,gpu2=TRUE)
				}
				checkmm(opf(testvecgpu[[i]], testmatgpu[[j]]), opf(testveccpu[[i]], testmatcpu[[j]]),mat1=FALSE,mat2=TRUE,gpu1=TRUE,gpu2=TRUE)
				checkmm(opf(testvecgpu[[i]], testmatcpu[[j]]), opf(testveccpu[[i]], testmatcpu[[j]]),mat1=FALSE,mat2=TRUE,gpu1=TRUE,gpu2=FALSE)
				checkmm(opf(testveccpu[[i]], testmatgpu[[j]]), opf(testveccpu[[i]], testmatcpu[[j]]),mat1=FALSE,mat2=TRUE,gpu1=FALSE,gpu2=TRUE)
				
				checkmm(opf(testvecgpu[[i]], testvecgpu[[j]]), opf(testveccpu[[i]], testveccpu[[j]]),mat1=FALSE,mat2=FALSE,gpu1=TRUE,gpu2=TRUE)
				checkmm(opf(testveccpu[[i]], testvecgpu[[j]]), opf(testveccpu[[i]], testveccpu[[j]]),mat1=FALSE,mat2=FALSE,gpu1=FALSE,gpu2=TRUE)
				checkmm(opf(testvecgpu[[i]], testveccpu[[j]]), opf(testveccpu[[i]], testveccpu[[j]]),mat1=FALSE,mat2=FALSE,gpu1=TRUE,gpu2=FALSE)		
			}
		}
	}
	#rm("%op%", envir=globalenv() )
	#warningslist
	cat("Checking outer product and kronecker product... \n")
	if(any(h(testvecgpu[[1]] %o% testvecgpu[[1]])!=testveccpu[[1]] %o% testveccpu[[1]] )) 
		warningslist<-c(warningslist,paste("%o% not working correctly"))
	if(any(h(testveccpu[[1]] %o% testvecgpu[[1]])!=testveccpu[[1]] %o% testveccpu[[1]] )) 
		warningslist<-c(warningslist,paste("%o% not working correctly"))
	if(any(h(testvecgpu[[1]] %o% testveccpu[[1]])!=testveccpu[[1]] %o% testveccpu[[1]] )) 
		warningslist<-c(warningslist,paste("%o% not working correctly"))
	if(any(h(testvecgpu[[1]] %o% testveccpu[[1]])!=testveccpu[[1]] %o% testveccpu[[1]] )) 
		warningslist<-c(warningslist,paste("%o% not working correctly"))
	if(type(testvecgpu[[2]] %o% testveccpu[[1]])!= "single") 
		warningslist<-c(warningslist,paste("%o% not working correctly (should return single)"))
	if(type(testveccpu[[1]] %o% testvecgpu[[2]])!= "single") 
		warningslist<-c(warningslist,paste("%o% not working correctly (should return single)"))
	
	if(any(h(gouter(testvecgpu[[1]] , testvecgpu[[1]]))!=testveccpu[[1]] %o% testveccpu[[1]] )) 
		warningslist<-c(warningslist,paste("%o% not working correctly"))
	if(any(h(gouter(testveccpu[[1]] , testvecgpu[[1]]))!=testveccpu[[1]] %o% testveccpu[[1]] )) 
		warningslist<-c(warningslist,paste("%o% not working correctly"))
	if(any(h(gouter(testvecgpu[[1]] , testveccpu[[1]]))!=testveccpu[[1]] %o% testveccpu[[1]] )) 
		warningslist<-c(warningslist,paste("%o% not working correctly"))
	
	if(any(h(gkroneckerProd(testmatgpu[[1]][1:10,1:10] , testmatgpu[[1]][1:10,1:10]))!=testmatcpu[[1]][1:10,1:10] %x% testmatcpu[[1]][1:10,1:10] )) 
		warningslist<-c(warningslist,paste("gkroneckerProd not working correctly"))
	if(any(h(gkroneckerProd(testmatcpu[[1]][1:10,1:10] , testmatgpu[[1]][1:10,1:10]))!=testmatcpu[[1]][1:10,1:10] %x% testmatcpu[[1]][1:10,1:10]  )) 
		warningslist<-c(warningslist,paste("gkroneckerProd not working correctly"))
	if(any(h(testmatgpu[[1]][1:10,1:10] %x% testmatcpu[[1]][1:10,1:10])!=testmatcpu[[1]][1:10,1:10] %x% testmatcpu[[1]][1:10,1:10] )) 
		warningslist<-c(warningslist,paste("%x% not working correctly"))
	if(any(h(testvecgpu[[1]] %x% testveccpu[[1]])!=testveccpu[[1]] %x% testveccpu[[1]])) 
		warningslist<-c(warningslist,paste("%x% not working correctly"))
	######################################################
	# Test Exchangeable + Nonexchangeable Binary elementwise operations    #
	######################################################
	exOps=list(
			c("mult","*",  "min", "2"), #name, operator, output type, maximum output type
			c("add" ,"+",  "min", "2"),
			c("eq"  ,"==", "3", "2"),
			c("ne"  ,"!=", "3", "2" ),
			c("and" ,"&",  "3", "2"),
			c("or"  ,"|",  "3", "2" ),
			c("sub","-",     "min", "2" ),
			c("div","/",     "min", "1"), #always returns at least a single
			c("pow","^",     "min", "1" ),#always returns at least a single
#		c("mod","%%",    "min", "2"), #cuda is not numericaly accurate on this...Todo: add note in
			c("gt"  ,">",    "3", "2"),
			c("lt"  ,"<",    "3", "2"),
			c("gte" ,">=",   "3", "2" ),
			c("lte" ,"<=",   "3", "2" )
	)
	checkbinobs=function(x1,y1,mat1,mat2,gpu1,gpu2, type1=i-1,type2=j-1, out_type=outt, op=binaryop, maxtype=as.integer(maxouttype)) {
		if(out_type=="min") {
			out_type=min(type1,type2)
			if(out_type>maxtype) {
				if(maxtype==2)
					out_type=maxtype
				else
					out_type=0L
			}
			if(gpu1 && !gpu2) {
				if(type1==1)
					out_type=1
				else if(type2==1)
					out_type=0
			}
			if(!gpu1 && gpu2) {
				if(type2==1)
					out_type=1
				else if(type1==1)
					out_type=0
			}
			
#		if(gpu1 && gpu2)
#			if(type1==1 || type2==1)
#				out_type=0
		} else
			out_type=as.numeric(out_type)
		
		yexpr=substitute(y1, list(i=i,j=j))
		xexpr=substitute(x1, list(i=i,j=j))
		mywarnings=character(0)
		x=tryCatch(eval(xexpr), error=function(e) {
					mywarnings<<-conditionMessage(e)
					return(-99999L)
				}
		)
		#browser()
		if(as.logical(x[1]!=-99999L)) {
			y=eval(yexpr)
			
			
			#	if(x@type<2)
			if(x@type!=out_type) {
				#	browser()
				mywarnings=c(mywarnings,paste("Incorrect type returned:",.type_name(x@type)))
			}
			
			x=h(x)
			#browser()
			delta=max(abs(x-y)/abs(x))
			#tmp=sum((y-x)^2)/sum(x^2)
			if(!is.finite(delta) || max(abs(x))<10^-6) {
				if(any(is.finite(y)!=is.finite(x)) & exOp[2]!="^") { #Todo: add note about "^"
					#browser()
					c(mywarnings,paste("not equal: rel error not defined"))
				}else {
					#		browser()
					good=is.finite(x) & is.finite(y)
					tmpx=x[good];tmpy=y[good]
					delta=(tmpx-tmpy)/tmpx
					delta=delta[is.finite(delta) & abs(tmpx)>10^-4]
					if(length(delta)>0)
						if(max(delta)>ifelse(max(type1,type2)==1| out_type==1,10^-4,10^-5)) {
							#browser()
							mywarnings=c(mywarnings,paste("not equal: rel error not defined"))
						}
				}
				
			} else if(delta>ifelse(max(type1,type2)==1 | out_type==1,10^-4,10^-9)) {
				mywarnings=c(mywarnings,paste("not equal: rel error = ", sum((y-x)^2)/sum(x^2)))
			}
		}
		tmpf=function(mat,gpu,type) {
			ap=ifelse(gpu,"g","")
			if(mat==0) 
				ret= paste("scaler " , ap, "vector",sep="")
			else if(mat==1)
				ret= paste(ap, "vector",sep="")
			else
				ret= paste(ap, "matrix",sep="")
			return(paste(ret,"of type",.type_name2(type)))
		}
		
		if(length(mywarnings)>0) {
			
			c1=tmpf(mat1,gpu1,type1)
			c2=tmpf(mat2,gpu2,type2)
			mywarnings=list(mywarnings)
			#browser()
			names(mywarnings)=paste("Operation " , op, " failed for: (" ,c1,") ", op," (" ,c2,").",sep="")
			warningslist<<-c(warningslist,mywarnings)
		}
		return(mywarnings)
	}
	

	cat("Checking Binary Operations... ")
	
	for(exOp in exOps) {
		#exOp=exOps[[8]]
		binaryop=exOp[2]
		maxouttype=exOp[4]
		outt=exOp[3]
		
		
#		suppressWarnings(
#				setGeneric("%op%",
#						function(a,b)
#							return(eval(parse(text=paste("a",binaryop,"b")))
#							),
#						where=globalenv()
#				)
#		)
		opf=function(a,b)
			return(eval(parse(text=paste("a",binaryop,"b"))))
		
		cat(binaryop," ")
		for(i in 1:4){
			for(j in 1:4){
				#	cat(i,j,"\n")
				
				checkbinobs(opf(testmatgpu[[i]], testmatgpu[[j]]), opf(testmatcpu[[i]], testmatcpu[[j]]),mat1=2,mat2=2,gpu1=TRUE,gpu2=TRUE)
				checkbinobs(opf(testmatcpu[[i]], testmatgpu[[j]]), opf(testmatcpu[[i]], testmatcpu[[j]]),mat1=2,mat2=2,gpu1=FALSE,gpu2=TRUE)
				checkbinobs(opf(testmatgpu[[i]], testmatcpu[[j]]), opf(testmatcpu[[i]], testmatcpu[[j]]),mat1=2,mat2=2,gpu1=TRUE,gpu2=FALSE)
				
				checkbinobs(opf(testmatgpu[[i]], testvecgpu[[j]]), opf(testmatcpu[[i]], testveccpu[[j]]),mat1=2,mat2=1,gpu1=TRUE,gpu2=TRUE)
				checkbinobs(opf(testmatgpu[[i]], testveccpu[[j]]), opf(testmatcpu[[i]], testveccpu[[j]]),mat1=2,mat2=1,gpu1=TRUE,gpu2=FALSE)
				checkbinobs(opf(testmatcpu[[i]], testvecgpu[[j]]), opf(testmatcpu[[i]], testveccpu[[j]]),mat1=2,mat2=1,gpu1=FALSE,gpu2=TRUE)
				
				checkbinobs(opf(testvecgpu[[i]], testmatgpu[[j]]), opf(testveccpu[[i]], testmatcpu[[j]]),mat1=1,mat2=2,gpu1=TRUE,gpu2=TRUE)
				checkbinobs(opf(testvecgpu[[i]], testmatcpu[[j]]), opf(testveccpu[[i]], testmatcpu[[j]]),mat1=1,mat2=2,gpu1=TRUE,gpu2=FALSE)
				checkbinobs(opf(testveccpu[[i]], testmatgpu[[j]]), opf(testveccpu[[i]], testmatcpu[[j]]),mat1=1,mat2=2,gpu1=FALSE,gpu2=TRUE)
				
				checkbinobs(opf(testvecgpu[[i]], testvecgpu[[j]]), opf(testveccpu[[i]], testveccpu[[j]]),mat1=1,mat2=1,gpu1=TRUE,gpu2=TRUE)
				checkbinobs(opf(testveccpu[[i]], testvecgpu[[j]]), opf(testveccpu[[i]], testveccpu[[j]]),mat1=1,mat2=1,gpu1=FALSE,gpu2=TRUE)
				checkbinobs(opf(testvecgpu[[i]], testveccpu[[j]]), opf(testveccpu[[i]], testveccpu[[j]]),mat1=1,mat2=1,gpu1=TRUE,gpu2=FALSE)
				
				checkbinobs(opf(testvecgpu[[i]], testscalergpu[[j]]), opf(testveccpu[[i]], testscalercpu[[j]]),mat1=1,mat2=0,gpu1=TRUE,gpu2=TRUE)
				checkbinobs(opf(testveccpu[[i]], testscalergpu[[j]]), opf(testveccpu[[i]], testscalercpu[[j]]),mat1=1,mat2=0,gpu1=FALSE,gpu2=TRUE)
				checkbinobs(opf(testvecgpu[[i]], testscalercpu[[j]]), opf(testveccpu[[i]], testscalercpu[[j]]),mat1=1,mat2=0,gpu1=TRUE,gpu2=FALSE)
				
				checkbinobs(opf(testscalergpu[[i]], testvecgpu[[j]]), opf(testscalercpu[[i]], testveccpu[[j]]),mat1=0,mat2=1,gpu1=TRUE,gpu2=TRUE)
				checkbinobs(opf(testscalergpu[[i]], testveccpu[[j]]), opf(testscalercpu[[i]], testveccpu[[j]]) ,mat1=0,mat2=1,gpu1=TRUE,gpu2=FALSE)
				checkbinobs(opf(testscalercpu[[i]], testvecgpu[[j]]), opf(testscalercpu[[i]], testveccpu[[j]]) ,mat1=0,mat2=1,gpu1=FALSE,gpu2=TRUE)
				
				checkbinobs(opf(testmatgpu[[i]], testscalergpu[[j]]), opf(testmatcpu[[i]], testscalercpu[[j]]),mat1=2,mat2=0,gpu1=TRUE,gpu2=TRUE)
				checkbinobs(opf(testmatcpu[[i]], testscalergpu[[j]]), opf(testmatcpu[[i]], testscalercpu[[j]]),mat1=2,mat2=0,gpu1=FALSE,gpu2=TRUE)
				checkbinobs(opf(testmatgpu[[i]], testscalercpu[[j]]), opf(testmatcpu[[i]], testscalercpu[[j]]),mat1=2,mat2=0,gpu1=TRUE,gpu2=FALSE)
				
				checkbinobs(opf(testscalergpu[[i]], testmatgpu[[j]]), opf(testscalercpu[[i]], testmatcpu[[j]]),mat1=0,mat2=2,gpu1=TRUE,gpu2=TRUE)
				checkbinobs(opf(testscalergpu[[i]], testmatcpu[[j]]), opf(testscalercpu[[i]], testmatcpu[[j]]) ,mat1=0,mat2=2,gpu1=TRUE,gpu2=FALSE)
				checkbinobs(opf(testscalercpu[[i]], testmatgpu[[j]]), opf(testscalercpu[[i]], testmatcpu[[j]]) ,mat1=0,mat2=2,gpu1=FALSE,gpu2=TRUE)
				
				
				checkbinobs(opf(testvecgpu[[i]], testsvecgpu[[j]]), opf(testveccpu[[i]], testsveccpu[[j]]),mat1=1,mat2=0,gpu1=TRUE,gpu2=TRUE)
				checkbinobs(opf(testveccpu[[i]], testsvecgpu[[j]]), opf(testveccpu[[i]], testsveccpu[[j]]),mat1=1,mat2=0,gpu1=FALSE,gpu2=TRUE)
				checkbinobs(opf(testvecgpu[[i]], testsveccpu[[j]]), opf(testveccpu[[i]], testsveccpu[[j]]),mat1=1,mat2=0,gpu1=TRUE,gpu2=FALSE)
				
				checkbinobs(opf(testsvecgpu[[i]], testvecgpu[[j]]), opf(testsveccpu[[i]], testveccpu[[j]]),mat1=0,mat2=1,gpu1=TRUE,gpu2=TRUE)
				checkbinobs(opf(testsvecgpu[[i]], testveccpu[[j]]), opf(testsveccpu[[i]], testveccpu[[j]]) ,mat1=0,mat2=1,gpu1=TRUE,gpu2=FALSE)
				checkbinobs(opf(testsveccpu[[i]], testvecgpu[[j]]), opf(testsveccpu[[i]], testveccpu[[j]]) ,mat1=0,mat2=1,gpu1=FALSE,gpu2=TRUE)
				
				
			}
		}
	}

	
	
	#########
	unOps=	list(
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
			c('abs','abs',  'x@type'),
			c('lgamma','lgamma',  'x@type'),
			c('gamma','gamma',  'x@type'),
			c('sign','sign',  'x@type'),
			c('round','round', '2L'),
			c('ceiling','ceil',  '2L'),
			c('floor','floor',  '2L'),
			c('is.na','isna',  '3L'),
			c('is.nan','isnan',  '3L'),
			c('is.finite','isfinite',  '3L'),
			c('is.infinite','isinfinite',  '3L'),
			c("!","!",""),
			c("-","-",""),
			c("+","+",""))
	
	checkunobs=function(x1,y, op=unaryop) {
		mywarnings=character()
		x=tryCatch(eval(x1), error=function(e) {
					mywarnings<<-conditionMessage(e)
					return(-99999L)
				}
		)
		
		if(as.logical(x[1]!=-99999L)) {
			
			x=h(x)
			y=suppressWarnings(y)
			tmp=sum((y-x)^2)/sum(y^2)
			if(!is.finite(tmp) || sum(x^2)<10^-3) {
				if(any(is.nan(y)!=is.nan(x)))
					c(mywarnings,paste("not equal: rel error not defined"))
				else {
					good=is.finite(x) & is.finite(y)
					if(sum((y[good]-x[good])^2)>10^-4)
						mywarnings=c(mywarnings,paste("not equal: rel error not defined"))
				}
				
			} else if(tmp>10^-4)
				mywarnings=c(mywarnings,paste("not equal: rel error = ", sum((y-x)^2)/sum(x^2)))
		}
		
		
		if(length(mywarnings)>0) {
			mywarnings=list(mywarnings)
			names(mywarnings)=paste("Operation " , unaryop, " failed.",sep="")
			warningslist<<-c(warningslist,mywarnings)
		}
		return(mywarnings)
	}
	cat("\nChecking Unary Operations/special functions... ")

	for(exOp in unOps) {
		unaryop=exOp[1]
		callexpr=parse(text=paste(unaryop,"(x)"))
		cat(unaryop," ")
		myop = function(x)
			return(eval(callexpr))
		
		for(i in 1:4) {
			checkunobs(myop(testvecgpu[[i]]), myop(testveccpu[[i]]))
			
		}
		
	}
	cat("\nChecking ifelse and which... \n")
	if(any(h(ifelse(gmatrix(c(T,F),4,4),c(1,2),c(1,2,3)))!=ifelse(matrix(c(T,F),4,4),c(1,2),c(1,2,3))))
		warningslist<-c(warningslist,"ifelse not working correctly")
	if(any(h(which(g.rep(c(T,F),100)))!=which(rep(c(T,F),100))))
		warningslist<-c(warningslist,"which or g.rep not working correctly")
	cat("Checking sort and order... \n")
	for(i in 0:2) {
		if(any(h(sort(as.gvector(c(1,4,3,2),type=i)))!=c(1,2,3,4)))
			warningslist<-c(warningslist,paste("sort not working correctly:",i))
		if(any(h(sort(as.gvector(c(1,4,3,2)),decreasing=TRUE))!=c(4,3,2,1)))
			warningslist<-c(warningslist,paste("sort decreasing=TRUE not working correctly:",i))
	}
	cat("Checking max, min, sum and col/row sums/means... \n")
	for(i in 0:2) {
		if(h(max(as.gvector(c(1,4,3,2),type=i)))!=4)
			warningslist<-c(warningslist,paste("max not working correctly for type:",i))
		if(h(min(as.gvector(c(1,4,3,2),type=i)))!=1)
			warningslist<-c(warningslist,paste("min not working correctly for type:",i))
		if(h(sum(as.gvector(c(1,4,3,2),type=i)))!=10)
			warningslist<-c(warningslist,paste("sum not working correctly for type:",i))
		if(h(mean(as.gvector(c(1,4,3,2),type=i)))!=2.5)
			warningslist<-c(warningslist,paste("mean not working correctly for type:",i))
		if(h(mean(as.gvector(c(1,4,3,2),type=i)))!=2.5)
			warningslist<-c(warningslist,paste("mean not working correctly for type:",i))
		checkunobs(rowMeans(testmatgpu[[i+1]]), rowMeans(testmatcpu[[i+1]]), "rowMeans" )
		checkunobs(rowSums( testmatgpu[[i+1]]), rowSums(testmatcpu[[i+1]]), "rowSums" )
		checkunobs(colMeans(testmatgpu[[i+1]]), colMeans(testmatcpu[[i+1]]), "colMeans" )
		checkunobs(colSums( testmatgpu[[i+1]]), colSums(testmatcpu[[i+1]]), "colSums" )
	}

	
	cat("Checking gsumby... \n")
	if(any(h(gsumby(1:10, c(1,5),c(4,10)))!=c(10,45))) 
		warningslist<-c(warningslist,paste("gsumby not working correctly"))
	
	cat("Checking transpose... \n")
	for(i in 1:4){
		if(any(h(t(testmatgpu[[i]])!=t(testmatcpu[[i]]))))
			warningslist<-c(warningslist,paste("Error in transpose type:",i-1))
		if(any(h(t(testvecgpu[[i]])!=t(testveccpu[[i]]))))
			warningslist<-c(warningslist,paste("Error in transpose type:",i-1))
	}
	
	cat("Checking logRowSums... \n")
	tmp=cbind(rep(7,10000),1,1,1)
	gtmp=g(log(tmp))
	if(any( abs(h(gRowLogSums(gtmp))-log(10))>10^-9))
		warningslist<-c(warningslist, "Error in 'gRowLogSums.'")
	
	cat("Checking indexing and diag functions... \n")
	for(i in 1:4) {
		for(j in 1:4)
			testvecgpu[[i]][2:5]=testvecgpu[[j]][7:10,]
		if(any(h(testvecgpu[[i]][2:5]!=testvecgpu[[j]][7:10]))) 
			warningslist<-c(warningslist,paste("Indexing problem for gvector"))
		
		testmatgpu[[i]][2:5,]=testmatgpu[[j]][7:10,]
		if(any(h(testmatgpu[[i]][2:5,]!=testmatgpu[[j]][7:10,]))) 
			warningslist<-c(warningslist,paste("Indexing problem for gvector"))
		
		testmatgpu[[i]][2:5,]=testmatgpu[[j]][7:10,]
		if(any(h(testmatgpu[[i]][2:5,]!=testmatgpu[[j]][7:10,]))) 
			warningslist<-c(warningslist,paste("Indexing problem for gvector"))
		
		testmatgpu[[i]][,2:5]=testmatgpu[[j]][,7:10]
		if(any(h(testmatgpu[[i]][,2:5]!=testmatgpu[[j]][,7:10]))) 
			warningslist<-c(warningslist,paste("Indexing problem for gvector"))
		
		testmatgpu[[i]][2:5,2:5]=testmatgpu[[j]][8:11,7:10]
		if(any(h(testmatgpu[[i]][2:5,2:5]!=testmatgpu[[j]][8:11,7:10]))) 
			warningslist<-c(warningslist,paste("Indexing problem for gvector"))
		
		diag(testmatgpu[[i]])=diag(testmatgpu[[j]])
		if(any(h(diag(testmatgpu[[i]])!=  convertType(diag(testmatgpu[[j]]) , as.integer(i-1))))) 
			warningslist<-c(warningslist,paste("Indexing problem for gvector"))
		
		diag(testmatgpu[[i]])=diag(testmatcpu[[j]])
		if(any(h(diag(testmatgpu[[i]])!=convertType(diag(testmatgpu[[j]]), as.integer(i-1)))))
			warningslist<-c(warningslist,paste("Indexing problem for gvector"))
		
		diag(testmatgpu[[i]])=(testscalercpu[[j]])
		if(any(h(diag(testmatgpu[[i]])!=convertType(testscalergpu[[j]] , as.integer(i-1))) ))
			warningslist<-c(warningslist,paste("Indexing problem for gvector"))
		
		diag(testmatgpu[[i]])=(testscalergpu[[j]])
		if(any(h(diag(testmatgpu[[i]])!=convertType(testscalergpu[[j]] , as.integer(i-1))) )) 
			warningslist<-c(warningslist,paste("Indexing problem for gvector"))
		
	}
	

	
	
	
	

	dosimcont =function(qdistfun=qnorm, rdistfun=rnorm, nsims=100000, cell_ct=100) {
		boundaries= qdistfun((0:cell_ct)/cell_ct)
		y=1:(cell_ct+1)
		myfun=approxfun(boundaries, y, method = "constant")
		mysims=as.numeric(rdistfun(nsims))
		mytbl=as.vector(table(myfun(mysims)))
		expected=nsims/cell_ct
		chisq = sum( (mytbl - expected)^2/expected)
		pval= pchisq(chisq,cell_ct-1, lower.tail=FALSE)
		return(pval)
	}
	
	dosimint =function(qdistfun, pdistfun, rdistfun, nsims=100000, max1=.99, max2=NULL) { #integer random vars
		if(is.null(max2)) {
			max2= qdistfun(max1)
		}
		boundaries=.5:(max2-.5)
		b_ct=length(boundaries)
		probs=pdistfun(boundaries)
		probs=c(probs[1],probs[2:b_ct]-probs[1:(b_ct-1)],1-probs[b_ct])
		mysims=as.numeric(rdistfun(nsims))
		mysims=ifelse(mysims>=max2,max2,mysims)
		#browser()
		mytbl=as.vector(table(factor(mysims, levels=0:max2)))
		expected=nsims*probs
		
		#deal with small expected cell counts so assymptotics works
		bad=expected<5
		if(sum(bad)>0) {
			prosb=c(probs[!bad], sum(probs[bad]))
			expected=c(expected[!bad], sum(expected[bad]))
			mytbl=c(mytbl[!bad], sum(mytbl[bad]))
			if(expected[length(expected)]<=.05) {
				bad=which(min(expected)==expected)[1]
				prosb=c(probs[-bad], sum(probs[bad]))
				expected=c(expected[!bad], sum(expected[bad]))
				mytbl=c(mytbl[-bad], sum(mytbl[bad]))
			}}
		
		#browser()
		#calculate statistic
		chisq = sum( (mytbl - expected)^2/expected)
		pval= pchisq(chisq,length(expected)-1, lower.tail=FALSE)
		return(pval)
	}
	
	
	
	
	mycheck=function(x, distfun) {
		mywarnings=character()
		ntests=length(x)
		x=if(any(x<.005/ntests))
			warningslist<<-c(warningslist,paste("Test of the goodness of fit of simulations from function",distfun, "yeilded the following p-values:",paste(as.vector(x), collapse=", ")))	
	}
	
	cat("Testing Random Generators...\n")
	#normal
	mycheck(sapply(c(.5,10), function(mean) {
						tmp = lapply(c(.01, .5, 1, 5, 10), function(sd) {
									#cat(mean,sd,"\n");
									dosimcont(
											qdistfun=function(q) qnorm(q,mean,sd=sd),
											rdistfun=function(n) grnorm(n,mean,sd=sd)
									)
								})
						return(unlist(tmp))
					}
			), "grnorm")
	
	
	#gamma
	mycheck(sapply(c(.5,10), function(rate) {
						tmp = lapply(c(.01, .5, 1, 5, 10), function(shape) {
									#cat(rate,shape,"\n");
									dosimcont(
											qdistfun=function(q) qgamma(q,shape,rate=rate),
											rdistfun=function(n) grgamma(n,shape,rate=rate)
									)
								})
						return(unlist(tmp))
					}
			), "grgamma")
	
	
	#beta
	mycheck(sapply(c(.25, .5, 1, 5, 10), function(shape1) {
						tmp = lapply(c(.5, 1, 5, 10), function(shape2) {
									#	cat(shape1,shape2,"\n");
									dosimcont(
											qdistfun=function(q) qbeta(q,shape1,shape2),
											rdistfun=function(n) grbeta(n,shape1,shape2),cell_ct=25
									)
								})
						return(unlist(tmp))
					}
			), "grbeta")
	
	#uniform
	mycheck(sapply(c(.01, .5, 1, 5, 10), function(min) {
						tmp = lapply(c(.01, .5, 1, 5, 10)+10, function(max) {
									dosimcont(
											qdistfun=function(q) qunif(q,min,max),
											rdistfun=function(n) grunif(n,min,max)
									)
								})
						return(unlist(tmp))
					}
			), "grunif")
	
	#binom
	mycheck(sapply(c(10),function(size) {
						tmp = lapply(c(.01, .25, .5, .8), function(p) {
									dosimint(
											qdistfun=function(q) qbinom(q,size,p),
											pdistfun=function(q) pbinom(q,size,p),
											rdistfun=function(n) grbinom(n,size,p),max2=size
									)})
						return(unlist(tmp))
					}
			), "grbinom")
	
	#pois
	mycheck(sapply(c(.01, .25, .5, .8), function(lamb) {
						
						dosimint(
								qdistfun=function(q) qpois(q,lamb),
								pdistfun=function(q) ppois(q,lamb),
								rdistfun=function(n) grpois(n,lamb),max1=.999999
						)
					}), "grpois")
	#todo: check rsample
	# tmp=cbind(rep(7,10000),1,1,1)
	# gtmp=g(log(tmp))
	# gRowLogSums(gtmp)
	# table(h(rsample(gtmp)))

	cat("Checking distribution functions...\n")
	#Distribution Checks
	checkd=function(a,b,distfun) {
		tmp=abs(as.numeric(a)-b)/b
		tmp=tmp[abs(b)>10^-20 & is.finite(b)]
		if(any(tmp > 10^-8))
			warningslist<<-c(warningslist,paste("Problem with function",distfun))	
			
	}
	checkd( gdnorm ((-10):10, mean=c(1,2), sd=c(1,2,3,4)), dnorm((-10):10, mean=c(1,2), sd=c(1,2,3,4)), "gdnorm")
	checkd( gdgamma((1:20)/5, shape=c(1,2), rate=c(1,2,3,4)) , dgamma((1:20)/5, shape=c(1,2), rate=c(1,2,3,4)),"gdgamma")
	checkd( gdunif((1:20)/5, min=c(1,2), max=3:6) , dunif((1:20)/5, min=c(1,2), max=3:6), "gdunif")
	checkd( gdbeta((1:20)/5, shape1=c(1,2), shape2=c(1,2,3,4)) , dbeta((1:20)/5, shape1=c(1,2), shape2=c(1,2,3,4)), "gdbeta")
	checkd( gdbinom(1:20, size=c(20,30), prob=c(1,2,3,4)/5) , dbinom(1:20, size=c(20,30), prob=c(1,2,3,4)/5), "gdbinom")
	checkd( gdpois(1:20, c(5,10)) , dpois(1:20, c(5,10)), "gdbinom")

	checkd( gdnorm ((-10):10, mean=c(1,2), sd=c(1,2,3,4), log=TRUE), dnorm((-10):10, mean=c(1,2), sd=c(1,2,3,4), log=TRUE), "gdnorm")
	checkd( gdgamma((1:20)/5, shape=c(1,2), rate=c(1,2,3,4), log=TRUE) , dgamma((1:20)/5, shape=c(1,2), rate=c(1,2,3,4), log=TRUE),"gdgamma")
	checkd( gdunif((1:20)/5, min=c(1,2), max=3:6, log=TRUE) , dunif((1:20)/5, min=c(1,2), max=3:6, log=TRUE), "gdunif")
	checkd( gdbeta((1:20)/5, shape1=c(1,2), shape2=c(1,2,3,4), log=TRUE) , dbeta((1:20)/5, shape1=c(1,2), shape2=c(1,2,3,4), log=TRUE), "gdbeta")
	checkd( gdbinom(1:20, size=c(20,30), prob=c(1,2,3,4)/5, log=TRUE) , dbinom(1:20, size=c(20,30), prob=c(1,2,3,4)/5, log=TRUE), "gdbinom")
	checkd( gdpois(1:20, c(5,10), log=TRUE) , dpois(1:20, c(5,10), log=TRUE), "gdbinom")
	
	checkd( gqnorm(log((0:10)/10), mean=c(1,2), sd=c(1,2,3,4), log.p=TRUE,warn=FALSE),
			qnorm(log((0:10)/10), mean=c(1,2), sd=c(1,2,3,4), log.p=TRUE), "gqnorm")
	checkd( gqnorm((0:10)/10, mean=c(1,2), sd=c(1,2,3,4),warn=FALSE),
			qnorm((0:10)/10, mean=c(1,2), sd=c(1,2,3,4)), "gqnorm")
	checkd( gpnorm(((-10):10)/2, mean=c(1,2), sd=c(1,2,3,4), log.p=TRUE,warn=FALSE),
			pnorm(((-10):10)/2, mean=c(1,2), sd=c(1,2,3,4), log.p=TRUE), "gpnorm")
	checkd( gpnorm(((-10):10)/2, mean=c(1,2), sd=c(1,2,3,4),warn=FALSE),
			pnorm(((-10):10)/2, mean=c(1,2), sd=c(1,2,3,4)), "gpnorm")
			
			

	if(.Call("cudaVersion")>=7000L) {
		cat("Checking solver functions...\n")
		checkD=function(x1,y1, funnm) {
			yexpr=substitute(y1)
			xexpr=substitute(x1)
			tmpnm=deparse(expr=xexpr)
			mywarnings=character(0)

			x=tryCatch(eval(xexpr), error=function(e) {
			mywarnings<<-conditionMessage(e)
			return(-99999L)})

			if(as.logical(x[1]!=-99999L)) {
				y=eval(yexpr)
				tmp=sum((y-x)^2)/sum(x^2)

				if(!is.numeric(tmp)||any(is.na(tmp)))
					mywarnings=c(mywarnings,paste(funnm, ":result missing or not numeric."))
				else if(tmp>10^(-6))
					mywarnings=c(mywarnings,paste(funnm, ":result not equal to host value."))
			}
			warningslist<<-c(warningslist,mywarnings)
		}

		#SVD
		hm=matrix(rnorm(10^2*2),20,10)
		gdm=g(hm)
		gsm=g(hm, type="s")
		checkD(svd(hm)$d, h(svd(gdm)@S), "svd double")
		checkD(svd(hm)$d, h(svd(gsm)@S), "svd single")

		#QR
		hb=matrix(rnorm(4),20,2)
		gdb=g(hb)
		gsb=g(hb, type="s")
		hcoefqr = qr.coef(qr(hm)   ,hb) 
		checkD(hcoefqr, h(gqr.coef(qr(gdm),gdb)) , "qr/gqr.coef double")
		checkD(hcoefqr, h(gqr.coef(qr(gsm),gsb)) , "qr/gqr.coef single")
		checkD(hcoefqr, h(gqr.coef(qr(gdm),hb))  , "qr/gqr.coef double / host")

		#QR
		hm =t(hm) %*% hm
		gdm=g(hm)
		gsm=g(hm, type="s")
		checkD(chol(hm), h(chol(gdm)), "chol double")
		checkD(chol(hm), h(chol(gsm)), "chol single")	
	}
	
	if(length(warningslist)==0)
		cat("No errors or warnings\n")
	else {
		cat("\n\nThe following errors/warnings were detected:\n")
		#browser()
		for(i in 1:length(warningslist)) {
			if(is.null(names(warningslist[i])))
				cat("Error:\n")
			else if(nchar(names(warningslist[i]))!=0){
				cat(names(warningslist[i]),":\n")
			} else{
				cat("Error:\n")
			}
			for(err in warningslist[i])
				cat("***",err,"\n")
		}
	}
	
	invisible(TRUE)
}
