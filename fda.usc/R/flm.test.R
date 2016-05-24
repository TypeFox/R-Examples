
##############################################
##          Functional Linear Model         ##
##             Goodness-of-fit              ##
##############################################

##############################################
## File created by Eduardo Garcia-Portugues ##
## using code from library fda.usc          ##
##############################################

# Generation of centred random variables with unit variance for the wild bootstrap
rwild=function(residuals,type="golden"){
  n=length(residuals)
  res=switch(type,
   golden=sample(c((1-sqrt(5))/2,(1+sqrt(5))/2),size=n,prob=c((5+sqrt(5))/10,(5-sqrt(5))/10),replace=TRUE), 
   Rademacher={	
  	sample(c(-1,1),size=n,prob=c(.5,.5),replace=TRUE)
  	},
  	normal=rnorm(n)
  	)
  return(residuals*res)
}

# Adot function to call Fortran code
Adot=function(X,inpr){
	
	# Check if the Fortran function is correctly loaded
	if(!is.loaded("adot")) stop("Fortran function adot is not loaded")
	
	# Compute the inproduct
	if(missing(inpr)){
		if(is.fd(X)){
			inpr=t(X$coefs)%*%inprod(X$basis,X$basis)%*%X$coefs
		}else if(is.fdata(X)){
			inpr=inprod.fdata(X,X)
		}else{
			stop("X is not a fd or fdata object")
		}
	}
	
	# Number of functions in X
	n=dim(inpr)[1]
	
	# Check for repeated interior products (will cause error in the adot subroutine)
	if(anyDuplicated(diag(inpr))){
		diag(inpr)=diag(inpr)*(1+rexp(n,rate=1000))
	}
	
	# Call Fortran function
	res=.Fortran("adot",n=as.integer(n),inprod=t(inpr)[upper.tri(inpr,diag=TRUE)],Adot_vec=as.double(rep(0,(n*n-n+2)/2)))$Adot_vec
	
	# Return result
	return(res)
	
}

# PCvM statistic
PCvM.statistic=function(X,residuals,p,Adot.vec){
	
	# Check if the Fortran function is correctly loaded
	if(!is.loaded("pcvm_statistic")) stop("Fortran function pcvm_statistic is not loaded")
	
	# Obtain the sample size
	n=length(residuals)
	
	# Check if the lengths of e and X are the same
	if(is.fd(X)){
		if(n!=dim(X$coefs)[2]) stop("Incompatible lenghts of X.fd and residuals")
	}else if(is.fdata(X)){
		if(n!=dim(X$data)[1]) stop("Incompatible lenghts of X.fdata and residuals")
	}
		
	# Compute Adot.vec if missing
	if(missing(Adot.vec)) Adot.vec=Adot(X)
	
	# Call Fortran function
	res=n^(-2)*pi^((p/2)-1)/gamma(p/2)*.Fortran("pcvm_statistic",n=as.integer(n),Adot_vec=Adot.vec,residuals=residuals,statistic=0)$statistic
	
	return(res)

}

# PCvM test for the composite hypothesis with bootstrap calibration
flm.test=function(X.fdata,Y,beta0.fdata=NULL,B=5000,est.method="pls",p=NULL,type.basis="bspline",verbose=TRUE,plot.it=TRUE,B.plot=100,G=200,...){
	
	# Check B.plot
	if(plot.it & B.plot>B) stop("B.plot must be less or equal than B")
	
	# Number of functions
	n=dim(X.fdata)[1]
	
	if(verbose) cat("Computing estimation of beta... ")	
	
	## COMPOSITE HYPOTHESIS ##
	if(is.null(beta0.fdata)){

		# Center the data first
		X.fdata=fdata.cen(X.fdata)$Xcen
		Y=Y-mean(Y)
		
		## 1. Optimal estimation of beta and the basis order ##
		
		if(est.method=="pc"){
			
			if(is.null(p)){

				# Method
				meth="PCvM test for the functional linear model using optimal PC basis representation"

				# Choose the number of basis elements: SICc is probably the best criteria
				mod.pc=fregre.pc.cv(fdataobj=X.fdata,y=Y,kmax=1:10,criteria="SICc")
				p.opt=length(mod.pc$pc.opt)
				ord.opt=mod.pc$pc.opt
				
				# PC components to be passed to the bootstrap
				pc.comp=mod.pc$fregre.pc$fdata.comp # pc.comp=mod.pc$fregre.pc$pc
				pc.comp$l=mod.pc$pc.opt
				
				# Express X.fdata and beta.est in the PC basis
				basis.pc=mod.pc$fregre.pc$fdata.comp$rotation
				if(length(pc.comp$l)!=1){
					X.est=fdata(mdata=mod.pc$fregre.pc$fdata.comp$x[,mod.pc$fregre.pc$l]%*%mod.pc$fregre.pc$fdata.comp$rotation$data[mod.pc$fregre.pc$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval) # X.est=fdata(mdata=mod.pc$fregre.pc$pc$x[,mod.pc$fregre.pc$l]%*%mod.pc$fregre.pc$pc$rotation$data[mod.pc$fregre.pc$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.pc$fregre.pc$fdata.comp$x[,mod.pc$fregre.pc$l]%*%t(mod.pc$fregre.pc$fdata.comp$rotation$data[mod.pc$fregre.pc$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval) # X.est=fdata(mdata=mod.pc$fregre.pc$pc$x[,mod.pc$fregre.pc$l]%*%t(mod.pc$fregre.pc$pc$rotation$data[mod.pc$fregre.pc$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.pc$fregre.pc$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.pc$fregre.pc$residuals
				
			}else{
			
				# Method
				meth=paste("PCvM test for the functional linear model using a representation in a PC basis of ",p,"elements") 
				
				# Estimation of beta on the given fixed basis
				mod.pc=fregre.pc(fdataobj=X.fdata,y=Y,l=1:p)
				p.opt=p
				ord.opt=mod.pc$l

				# PC components to be passed to the bootstrap
				pc.comp=mod.pc$pc
				pc.comp$l=mod.pc$l
				
				# Express X.fdata and beta.est in the basis
				if(p!=1){
					X.est=fdata(mdata=mod.pc$fdata.comp$x[,mod.pc$l]%*%mod.pc$fdata.comp$rotation$data[mod.pc$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.pc$fdata.comp$x[,mod.pc$l]%*%t(mod.pc$fdata.comp$rotation$data[mod.pc$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.pc$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.pc$residuals
			
			}
			
		}else if(est.method=="pls"){
			
			if(is.null(p)){

				# Method
				meth="PCvM test for the functional linear model using optimal PLS basis representation"
			
				# Choose the number of the basis: SICc is probably the best criteria
				mod.pls=fregre.pls.cv(fdataobj=X.fdata,y=Y,kmax=10,criteria="SICc") 
				p.opt=length(mod.pls$pls.opt)
				ord.opt=mod.pls$pls.opt
				
				# PLS components to be passed to the bootstrap
				pls.comp=mod.pls$fregre.pls$fdata.comp
				pls.comp$l=mod.pls$pls.opt
						
				# Express X.fdata and beta.est in the PLS basis
				basis.pls=mod.pls$fregre.pls$fdata.comp$rotation
				if(length(pls.comp$l)!=1){
					X.est=fdata(mdata=mod.pls$fregre.pls$fdata.comp$x[,mod.pls$fregre.pls$l]%*%mod.pls$fregre.pls$fdata.comp$rotation$data[mod.pls$fregre.pls$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.pls$fregre.pls$fdata.comp$x[,mod.pls$fregre.pls$l]%*%t(mod.pls$fregre.pls$fdata.comp$rotation$data[mod.pls$fregre.pls$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.pls$fregre.pls$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.pls$fregre.pls$residuals
				
			}else{
			
				# Method
				meth=paste("PCvM test for the functional linear model using a representation in a PLS basis of ",p,"elements") 
			
				# Estimation of beta on the given fixed basis
				mod.pls=fregre.pc(fdataobj=X.fdata,y=Y,l=1:p)
				p.opt=p
				ord.opt=mod.pls$l

				# PLS components to be passed to the bootstrap
				pls.comp=mod.pls$fdata.comp
				pls.comp$l=mod.pls$l
				
				# Express X.fdata and beta.est in the basis
				if(p!=1){
					X.est=fdata(mdata=mod.pls$fdata.comp$x[,mod.pls$l]%*%mod.pls$fdata.comp$rotation$data[mod.pls$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}else{
					X.est=fdata(mdata=mod.pls$fdata.comp$x[,mod.pls$l]%*%t(mod.pls$fdata.comp$rotation$data[mod.pls$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
				}
				beta.est=mod.pls$beta.est
				norm.beta.est=norm.fdata(beta.est)
				
				# Compute the residuals
				e=mod.pls$residuals
			
			}
			
		}else if(est.method=="basis"){
			
			if(is.null(p)){
			
				# Method
				meth=paste("PCvM test for the functional linear model using optimal",type.basis,"basis representation")
			
				# Choose the number of the bspline basis with GCV.S
				mod.basis=fregre.basis.cv(fdataobj=X.fdata,y=Y,basis.x=seq(5,30,by=1),basis.b=NULL,type.basis=type.basis,type.CV=GCV.S,verbose=FALSE,...)
				p.opt=mod.basis$basis.x.opt$nbasis
				ord.opt=1:p.opt
				
				# Express X.fdata and beta.est in the optimal basis
				basis.opt=mod.basis$basis.x.opt
				X.est=mod.basis$x.fd
				beta.est=mod.basis$beta.est
				norm.beta.est=norm.fd(beta.est)
				
				# Compute the residuals
				e=mod.basis$residuals
			
			}else{
				
				# Method
				meth=paste("PCvM test for the functional linear model using a representation in a",type.basis,"basis of ",p,"elements") 

				# Estimation of beta on the given fixed basis
				basis.opt=do.call(what=paste("create.",type.basis,".basis",sep=""),args=list(rangeval=X.fdata$rangeval,nbasis=p,...))
				mod.basis=fregre.basis(fdataobj=X.fdata,y=Y,basis.x=basis.opt,basis.b=basis.opt)
				p.opt=p
				ord.opt=1:p.opt

				# Express X.fdata and beta.est in the basis
				X.est=mod.basis$x.fd
				beta.est=mod.basis$beta.est
				norm.beta.est=norm.fd(beta.est)
				
				# Compute the residuals
				e=mod.basis$residuals
			
			}
			
		}else{
		
			stop(paste("Estimation method",est.method,"not implemented."))
			
		}
	
	## SIMPLE HYPOTHESIS ##
	}else{
	
		## 1. Optimal representation of X and beta0 ##
		
		# Choose the number of basis elements
		if(type.basis!="pc" & type.basis!="pls"){
			
			# Basis method: select the number of elements of the basis is it is not given
			if(is.null(p)){
				
				# Method
				meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in an optimal ",type.basis,"basis") 
				
				p.opt=min.basis(X.fdata,type.basis=type.basis,numbasis=seq(31,71,by=2),verbose=FALSE)$numbasis.opt
				if(p.opt==31){
					cat("Readapting interval...\n")
					p.opt=min.basis(X.fdata,type.basis=type.basis,numbasis=seq(5,31,by=2),verbose=FALSE)$numbasis.opt
				}else if(p.opt==71){
					cat("Readapting interval...\n")
					p.opt=min.basis(X.fdata,type.basis=type.basis,numbasis=seq(71,101,by=2),verbose=FALSE)$numbasis.opt
				}
				
				ord.opt=1:p.opt
				
			}else{
				
				# Method
				meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in a ",type.basis," basis of ",p,"elements") 
				
				p.opt=p
				ord.opt=1:p.opt
			}
			
			# Express X.fdata in the basis
			X.est=fdata2fd(X.fdata,type.basis=type.basis,nbasis=p)
			beta.est=fdata2fd(beta0.fdata,type.basis=type.basis,nbasis=p)
			
			# Compute the residuals
			e=Y-inprod(X.est,beta.est)

		}else if (type.basis=="pc"){
			
			if(is.null(p)) stop("Simple hypothesis with type.basis=\"pc\" need the number of components p") 
			# Method
			meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in a PC basis of ",p,"elements") 
			
			# Express X.fdata in a PC basis
			fd2pc=fdata2pc(X.fdata,ncomp=p)
			if(length(fd2pc$l)!=1){
				X.est=fdata(mdata=fd2pc$x[,fd2pc$l]%*%fd2pc$rotation$data[fd2pc$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}else{
				X.est=fdata(mdata=fd2pc$x[,fd2pc$l]%*%t(fd2pc$rotation$data[fd2pc$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}
			beta.est=beta0.fdata
			p.opt=p
			ord.opt=1:p.opt
			
			# Compute the residuals
			e=Y-inprod.fdata(X.est,beta.est)
			
		}else if(type.basis=="pls"){
		
			if(is.null(p)) stop("Simple hypothesis with type.basis=\"pls\" need the number of components p")

			# Method
			meth=paste("PCvM test for the simple hypothesis in a functional linear model, using a representation in a PLS basis of ",p,"elements") 
			
			# Express X.fdata in a PLS basis
			fd2pls=fdata2pls(X.fdata,Y,ncomp=p)
			if(length(fd2pls$l)!=1){
				X.est=fdata(mdata=fd2pls$x[,fd2pls$l]%*%fd2pls$rotation$data[fd2pls$l,],argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}else{
				X.est=fdata(mdata=fd2pls$x[,fd2pls$l]%*%t(fd2pls$rotation$data[fd2pls$l,]),argvals=X.fdata$argvals,rangeval=X.fdata$rangeval)
			}		
			beta.est=beta0.fdata
			p.opt=p
			ord.opt=1:p.opt
			
			# Compute the residuals
			e=Y-inprod.fdata(X.est,beta.est)
			
		}else{
		
			stop(paste("Type of basis method",type.basis,"not implemented."))

		}
	
	}
	
	## 2. Bootstrap calibration ##
	
	# Start up
	pcvm.star=numeric(B)
	e.hat.star=matrix(ncol=n,nrow=B)
	
	# Calculus of the Adot.vec
	Adot.vec=Adot(X.est)

	# REAL WORLD
	pcvm=PCvM.statistic(X=X.est,residuals=e,p=p.opt,Adot.vec=Adot.vec)
	
	# BOOTSTRAP WORLD 
	if(verbose) cat("Done.\nBootstrap calibration...\n ")
	if(verbose) pb=txtProgressBar(style=3)

	## COMPOSITE HYPOTHESIS ##
	if(is.null(beta0.fdata)){
					
		# Calculate the design matrix of the linear model
		# This allows to resample efficiently the residuals without estimating again the beta
		if(est.method=="pc"){
			
			if(is.null(p)){
				# Design matrix for the PC estimation
				X.matrix=mod.pc$fregre.pc$lm$x
			}else{
				# Design matrix for the PC estimation
				X.matrix=mod.pc$lm$x
			}
			
		}else if(est.method=="pls"){

			if(is.null(p)){	
				# Design matrix for the PLS estimation
				X.matrix=mod.pls$fregre.pls$lm$x
			}else{
				# Design matrix for the PLS estimation
				X.matrix=mod.pls$lm$x
			}		
		}else if(est.method=="basis"){

			if(is.null(p)){	
				# Design matrix for the basis estimation
				X.matrix=mod.basis$lm$x
			}else{
				# Design matrix for the basis estimation
				X.matrix=mod.basis$lm$x
			}		
		}
	  
	  # Projection matrix
	  P=(diag(rep(1,n))-X.matrix%*%solve(t(X.matrix)%*%X.matrix)%*%t(X.matrix))
	  
		# Bootstrap resampling
		for(i in 1:B){
		
			# Generate bootstrap errors
			e.hat=rwild(e,"golden")
	
			# Calculate Y.star
			Y.star=Y-e+e.hat
			
			# Residuals from the bootstrap estimated model
			e.hat.star[i,]=P%*%Y.star
		
			# Calculate PCVM.star
			pcvm.star[i]=PCvM.statistic(X=X.est,residuals=e.hat.star[i,],p=p.opt,Adot.vec=Adot.vec)
			
			if(verbose) setTxtProgressBar(pb,i/B)
		
		}
					
	## SIMPLE HYPOTHEIS ##
	}else{
		
		# Bootstrap resampling
		for(i in 1:B){
				
			# Generate bootstrap errors
			e.hat.star[i,]=rwild(e,"golden")
		
			# Calculate PCVM.star
			pcvm.star[i]=PCvM.statistic(X=X.est,residuals=e.hat.star[i,],p=p.opt,Adot.vec=Adot.vec)
				
			if(verbose) setTxtProgressBar(pb,i/B)
		
		}
					
	}

	## 3. MC estimation of the p-value and order the result ##
	
	# Compute the p-value
	pvalue=sum(pcvm.star>pcvm)/B
	
	## 4. Graphical representation of the integrated process ##
	
	if(verbose) cat("\nDone.\nComputing graphical representation... ")
	if(is.null(beta0.fdata) & plot.it){
	
		gamma=rproc2fdata(n=G,t=X.fdata$argvals,sigma="OU",par.list=list(theta=2/diff(range(X.fdata$argvals))))
		gamma=gamma/drop(norm.fdata(gamma))
		ind=drop(inprod.fdata(X.fdata,gamma))
		
		r=0.9*max(max(ind),-min(ind))
		u=seq(-r,r,l=200)
		mean.proc=numeric(length(u))
		mean.boot.proc=matrix(ncol=length(u),nrow=B.plot)
		
		res=apply(ind,2,function(iind){
			iind.sort=sort(iind,index.return=TRUE)
			stepfun(x=iind.sort$x,y=c(0,cumsum(e[iind.sort$ix])))(u)/sqrt(n)
		}
		)
		mean.proc=apply(res,1,mean)
	
		for(i in 1:B.plot){
		
			res=apply(ind,2,function(iind){
				iind.sort=sort(iind,index.return=TRUE)
				stepfun(x=iind.sort$x,y=c(0,cumsum(e.hat.star[i,iind.sort$ix])))(u)/sqrt(n)
			}
			)
			mean.boot.proc[i,]=apply(res,1,mean)
		
		}
		
		# Plot
		dev.new()
		plot(u,mean.proc,ylim=c(min(mean.proc,mean.boot.proc),max(mean.proc,mean.boot.proc))*1.05,type="l",xlab=expression(paste(symbol("\341"),list(X, gamma),symbol("\361"))),ylab=expression(R[n](u)),main="")
		for(i in 1:B.plot) lines(u,mean.boot.proc[i,],lty=2,col=gray(0.8))
		lines(u,mean.proc)
		text(x=0.75*u[1],y=0.75*min(mean.proc,mean.boot.proc),labels=sprintf("p-value=%.3f",pvalue))
		
	}
	if(verbose) cat("Done.\n")

	
	# Result: class htest
	names(pcvm)="PCvM statistic"
	result=structure(list(statistic=pcvm,boot.statistics=pcvm.star,p.value=pvalue,method=meth,B=B,type.basis=type.basis,beta.est=beta.est,p=p.opt,ord=ord.opt,data.name="Y=<X,b>+e"))
							
	class(result)="htest"
	return(result)
	
}


