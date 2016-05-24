
###################################
## Test using random projections ## 
###################################

# Test for the FLM using random projections
rp.flm.test=function(X.fdata,Y,beta0.fdata=NULL,est.method="pc",p=NULL,type.basis="bspline",B=5000,n.proj=50,verbose=TRUE,F.code=TRUE,sigma="vexponential",par.list=list(theta=diff(range(X.fdata$argvals))/20),same.rwild=FALSE,...){ 
	
	# Arguments:
	# X.fdata: Functional observations in the class fdata0
	# Y: scalar responses for the FLM
	# beta0: if NULL, composite hypothesis is tested. If !NULL, the simple hypothesis is contrasted (beta0 must be in the fdata class)
	# est.method: method for the estimation (composite hypothesis) or representation (simple hypothesis) of the beta0 coefficient either "pc", "pls" or "basis"
	# p: number of elements for the basis representation of beta0 and X with the est.method. If not supplied it is estimated from the data
	# type.basis: type of basis if est.method="basis"
	# B: number of bootstrap replicates
	# n.proj: number of projections (it can be a vector)
	# verbose: whether to show or not information about the testing progress
	
	# Sample size
	n=dim(X.fdata)[1]
	
	# Number of projections
	if(length(n.proj)>1){
    if(any(sort(n.proj)!=n.proj)) stop("vector n.proj must be sorted")
		vec.nproj=n.proj
		n.proj=max(n.proj)
	}else{
		vec.nproj=n.proj
	}
	
	if(verbose) cat("Computing estimation of beta... ")	

	## COMPOSITE HYPOTHESIS: optimal estimation of beta and the basis order ##
	if(is.null(beta0.fdata)){

		# Center the data first
		X.fdata=fdata.cen(X.fdata)$Xcen
		Y=Y-mean(Y)
		
		# PC
		if(est.method=="pc"){
			
			# Optimal p by SIC criterion
			if(is.null(p)){

				# Method
				meth="Random projection based test for the functional linear model using optimal PC basis representation"

				# Choose the number of basis elements: SIC is probably the best criteria
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
			
			# Fixed p
			}else{
			
				# Method
				meth=paste("Random projection based test for the functional linear model using a representation in a PC basis of ",p,"elements") 
				
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
		
		# PLS		
		}else if(est.method=="pls"){
			
			# Optimal p by SIC criterion
			if(is.null(p)){

				# Method
				meth="Random projection based test for the functional linear model using optimal PLS basis representation"
			
				# Choose the number of the basis: SIC is probably the best criteria
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
				
			# Fixed p
			}else{
			
				# Method
				meth=paste("Random projection based test for the functional linear model using a representation in a PLS basis of ",p,"elements") 
			
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
		
		# Deterministic basis
		}else if(est.method=="basis"){
			
			# Optimal p by GCV criterion
			if(is.null(p)){
			
				# Method
				meth=paste("Random projection based test for the functional linear model using optimal",type.basis,"basis representation")
			
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
			
			# Fixed p
			}else{
				
				# Method
				meth=paste("Random projection based test for the functional linear model using a representation in a",type.basis,"basis of ",p,"elements") 

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
	
	## REPRESENTATION: SIMPLE HYPOTHESIS: optimal representation of X and beta0  ##
	}else{
		
		# Method
		meth=paste("Random projection based test for the simple hypothesis in a functional linear model") 
		
		# Express X.fdata in a PC basis
		beta.est=beta0.fdata
		p.opt=NA
		
		# Compute the residuals
		e=drop(Y-inprod.fdata(X.fdata,beta.est))
			
	}
	
	## CALCULUS OF THE STATISTIC ##
	
	# Compute random projections for the statistic and bootstrap replicates
	if(verbose) cat("Done.\nComputing Projections... ")
  proj.gamma=rproc2fdata(n.proj,t=argvals(X.fdata),sigma=sigma,par.list=par.list,norm=TRUE)
  proj.X=inprod.fdata(X.fdata,proj.gamma) # matrix n x n.proj
	
	# Statistic
	rp.stat=rp.flm.statistic(proj.X=proj.X,residuals=e,verbose=FALSE,F.code=F.code)

	## BOOTSTRAP CALIBRATION ##
	
	# Start up
  rp.stat.star=array(NA,dim=c(n.proj,2,B))
  e.hat.star=array(NA,dim=c(B,n.proj,n))
	if(verbose) cat("Done.\nBootstrap calibration...\n ")
	if(verbose) pb=txtProgressBar(style=3)

	## COMPOSITE HYPOTHESIS ##
	if(is.null(beta0.fdata)){

		# Calculate the design matrix of the linear model depending on the chosen basis
		# This allows to resample efficiently the residuals without estimating again the beta

		# PC
		if(est.method=="pc"){
			
			# Design matrix for the PC estimation
			if(is.null(p)){
				X.matrix=mod.pc$fregre.pc$lm$x
			}else{
				X.matrix=mod.pc$lm$x
			}
			
		# PLS
		}else if(est.method=="pls"){

			# Design matrix for the PLS estimation
			if(is.null(p)){	
				X.matrix=mod.pls$fregre.pls$lm$x
			}else{
				X.matrix=mod.pls$lm$x
			}		
		
		# Deterministic basis
		}else if(est.method=="basis"){

			# Design matrix for the basis estimation
			if(is.null(p)){	
				X.matrix=mod.basis$lm$x
			}else{
				X.matrix=mod.basis$lm$x
			}	
			
		}
		
    # Projection matrix of the linear model
    lm.proj.X.matrix=X.matrix%*%solve(t(X.matrix)%*%X.matrix)%*%t(X.matrix)
    
		# Bootstrap resampling
		for(i in 1:B){
		
			# Generate bootstrap errors
			if(same.rwild){
				e.hat=matrix(rwild(e,"golden"),nrow=n.proj,ncol=n,byrow=TRUE)
			}else{
				e.hat=matrix(rwild(rep(e,n.proj),"golden"),nrow=n.proj,ncol=n,byrow=TRUE)
			}
			
			# Calculate Y.star
			Y.star=sweep(e.hat,2,Y-e,"+")

			# Residuals from the bootstrap estimated model
			e.hat.star[i,,]=Y.star%*%(diag(rep(1,n))-lm.proj.X.matrix)
		
			# Calculate the bootstrap statistics
			rp.stat.star[,,i]=rp.flm.statistic(proj.X=proj.X,residuals=e.hat.star[i,,],proj.X.ord=rp.stat$proj.X.ord,verbose=FALSE,F.code=F.code)$statistic
			
			if(verbose) setTxtProgressBar(pb,i/B)
		
		}

	## SIMPLE HYPOTHEIS ##
	}else{
		
		# Bootstrap resampling
		for(i in 1:B){

      # Generate bootstrap errors			
			if(same.rwild){
				e.hat.star[i,,]=matrix(rwild(e,"golden"),nrow=n.proj,ncol=n,byrow=TRUE)
			}else{
				e.hat.star[i,,]=matrix(rwild(rep(e,n.proj),"golden"),nrow=n.proj,ncol=n,byrow=TRUE)
			}
		
			# Calculate the bootstrap statistics
			rp.stat.star[,,i]=rp.flm.statistic(proj.X=proj.X,residuals=e.hat.star[i,,],proj.X.ord=rp.stat$proj.X.ord,verbose=FALSE)$statistic

			if(verbose) setTxtProgressBar(pb,i/B)
		
		}

	}
	
	# Compute the p-values
  pvalue=function(vec,value){mean(vec>value)}
	pval=matrix(nrow=n.proj,ncol=2)
	for (i in 1:n.proj){
		pval[i,1]=pvalue(rp.stat.star[i,1,],rp.stat$statistic[i,1]) # CvM
		pval[i,2]=pvalue(rp.stat.star[i,2,],rp.stat$statistic[i,2]) # KS
	}
	
	# Compute p-values depending for the vector of projections
	if(length(vec.nproj)==1){
	  if(vec.nproj==1){
	    rp.pvalue=pval[1,]
	  }else{
  		rp.pvalue=apply(pval,2,function(x){l=length(x);return(min(l/(1:l)*sort(x)))})
	  }
		names(rp.pvalue)=names(pval)=c("CvM","KS")
	}else{
		rp.pvalue=matrix(nrow=length(vec.nproj),ncol=2)
		for(k in seq_along(vec.nproj)){
      if(vec.nproj[k]==1){
        rp.pvalue[k,]=pval[1,]
      }else{
  			rp.pvalue[k,]=apply(pval[1:vec.nproj[k],],2,function(x){l=length(x);return(min(l/(1:l)*sort(x)))})
      }
		}
		colnames(rp.pvalue)=colnames(pval)=c("CvM","KS")
		rownames(rp.pvalue)=vec.nproj
	}
	
	# Return result
	if(verbose) cat("\nDone.\n")
  result=structure(list(statistics=rp.stat$statistic,boot.statistics=rp.stat.star,p.value=rp.pvalue,proj.p.values=pval,method=meth,B=B,
							n.proj=vec.nproj,type.basis=type.basis,beta.est=beta.est,p=p.opt,data.name="Y=<X,b>+e"))
	class(result)="htest"
  return(result)
	
}






# Statistic for testing the FLM using random projections
rp.flm.statistic=function(proj.X,residuals,proj.X.ord=NULL,verbose=FALSE,F.code=TRUE){
	
	# Sample size and number of projections
	n.proj=ncol(proj.X)
	if(is.null(dim(residuals))){
		residuals=matrix(residuals,nrow=n.proj,ncol=length(residuals),byrow=TRUE)
	}
	n=dim(residuals)[2]
	if(dim(residuals)[1]!=n.proj) stop("residuals should be a matrix of n.proj x n")
	
	# Matrix of statistics (CVM and KS, columns) projected in n.proj projections (rows)
  rp.stat=matrix(0,nrow=n.proj,ncol=2)
	
	# Order projections and compute their ranges
	if (verbose) cat("Computing Projections...\n")
	if (is.null(proj.X.ord)) proj.X.ord=apply(proj.X,2,order)
	if (verbose) cat("Computing Statistics...\n")
	
	if(F.code){
		
		# Check dll is loaded
		if(!is.loaded("rp_stat")) stop("rp_stat not loaded")
		
		# Statistic
		rp.stat=.Fortran("rp_stat",proj_X_ord=proj.X.ord,residuals=residuals,n_proj=n.proj,n=n,rp_stat_proj=rp.stat)$rp_stat_proj

	}else{
		
		# R implementation
		for (i in 1:n.proj){
			
			# Empirical process 
			y=cumsum(residuals[i,proj.X.ord[,i]])
			
			# Statistics (CVM and KS, rows)
			CVM=sum(y^2)/n^2
			KS=max(abs(y))/sqrt(n)
			rp.stat[i,]=c(CVM,KS)

		}

	}
	
	colnames(rp.stat)=c("CvM","KS")
	return(list(statistic=rp.stat,proj.X.ord=proj.X.ord))

}
