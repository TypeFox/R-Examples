## Lavaan to ram
lavaan2ram<-function(fitModel, digits=2, zero.print="0", ram.out=TRUE, fit=FALSE){
	parTable<-fitModel@ParTable
	parEst<-fitModel@Fit@est
	parSE<-fitModel@Fit@se
	fitInd<-fitMeasures(fitModel)
	ngroup<-fitModel@Data@ngroups
	if (ngroup>1){
		A<-S<-Ase<-Sse<-M<-Mse<-list()
		for(g in 1:ngroup) {
        	obsVar <- lavNames(parTable, "ov", group=g)
        	latVar <- lavNames(parTable, "lv", group=g)
        	varName<-c(obsVar, latVar)
			manifest<-length(obsVar)
			latent<-length(latVar)
			nrow<-length(varName)
			A[[g]]<-S[[g]]<-Ase[[g]]<-Sse[[g]]<-matrix(0,nrow,nrow,dimnames=list(varName, varName))
			M[[g]]<-Mse[[g]]<-matrix(0,nrow,1,dimnames=list(varName, 'M'))
			
			for (j in parTable$id){
				if (parTable$group[j]==g){
					if (parTable$op[j]=="~"){
						A[[g]][parTable$lhs[j], parTable$rhs[j]]<-parEst[j]
						Ase[[g]][parTable$lhs[j], parTable$rhs[j]]<-parSE[j]
					}
					if (parTable$op[j]=="=~"){
						A[[g]][parTable$rhs[j], parTable$lhs[j]]<-parEst[j]
						Ase[[g]][parTable$rhs[j], parTable$lhs[j]]<-parSE[j]
					}
					if (parTable$op[j]=="~~"){
						S[[g]][parTable$lhs[j], parTable$rhs[j]]<-parEst[j]
						Sse[[g]][parTable$lhs[j], parTable$rhs[j]]<-parSE[j]
						S[[g]][parTable$rhs[j], parTable$lhs[j]]<-parEst[j]
						Sse[[g]][parTable$rhs[j], parTable$lhs[j]]<-parSE[j]
					}
					if (parTable$op[j]=="~1"){
						M[[g]][parTable$lhs[j], 'M']<-parEst[j]
						Mse[[g]][parTable$lhs[j], 'M']<-parSE[j]
					}
				}
			}
			
		}
		
		## Print some results
		## Print the Model fit
		if (fit){
			cat("Model fit statistics and indices\n")
			print(fitInd)
		}
		if (ram.out){
			for (g in 1:ngroup){
				## Print the parameter estimates
				A.na<-A[[g]]
				A.na[A[[g]]==0]<-NA
				S.na<-S[[g]]
				S.na[S[[g]]==0]<-NA
				Ase.na<-Ase[[g]]
				Ase.na[Ase[[g]]==0]<-NA
				Sse.na<-Sse[[g]]
				Sse.na[Sse[[g]]==0]<-NA
				cat("\n--------------------\n")
	 			cat("Group ",g,"\n")
	  			cat("--------------------\n")
			
	  			cat("\n--------------------\n")
	  			cat("Parameter estimates:\n")
	  			cat("--------------------\n")
	  			cat("\nMatrix A\n\n")
	  			print(A.na, digits=digits,na.print = zero.print)
	  			cat("\nMatrix S\n\n")
	  			print(S.na,digits=digits,na.print = zero.print)
	  
	  			cat("\n----------------------------------------\n")
	  			cat("Standard errors for parameter estimates:\n")
	  			cat("----------------------------------------\n")
	  			cat("\nMatrix A\n\n")
	  			print(Ase.na,digits=digits,na.print = zero.print)
	  			cat("\nMatrix S\n\n")
	  			print(Sse.na,digits=digits,na.print = zero.print)
	  			cat("\n\n")
			}
		}
	  	lname<-NULL
  		if (latent>0) lname = varName[(manifest+1):nrow]
		invisible(list(A=A, S=S, Ase=Ase, Sse=Sse, M=M, Mse=Mse, fit=fitInd, lavaan=fitModel, nvar=nrow, manifest=manifest,latent=latent,lname=lname,varname=varName))        	
	}else{
		obsVar <- lavNames(parTable, "ov")
        latVar <- lavNames(parTable, "lv")
		varName<-c(obsVar, latVar)
		manifest<-length(obsVar)
		latent<-length(latVar)
	
		nrow<-length(varName)
		A<-S<-Ase<-Sse<-matrix(0,nrow,nrow,dimnames=list(varName, varName))
		M<-Mse<-matrix(0,nrow,1,dimnames=list(varName, 'M'))
		for (j in parTable$id){
			if (parTable$op[j]=="~"){
				A[parTable$lhs[j], parTable$rhs[j]]<-parEst[j]
				Ase[parTable$lhs[j], parTable$rhs[j]]<-parSE[j]
			}
			if (parTable$op[j]=="=~"){
				A[parTable$rhs[j], parTable$lhs[j]]<-parEst[j]
				Ase[parTable$rhs[j], parTable$lhs[j]]<-parSE[j]
			}
			if (parTable$op[j]=="~~"){
				S[parTable$lhs[j], parTable$rhs[j]]<-parEst[j]
				Sse[parTable$lhs[j], parTable$rhs[j]]<-parSE[j]
				S[parTable$rhs[j], parTable$lhs[j]]<-parEst[j]
				Sse[parTable$rhs[j], parTable$lhs[j]]<-parSE[j]
			}
			if (parTable$op[j]=="~1"){
				M[parTable$lhs[j], 'M']<-parEst[j]
				Mse[parTable$lhs[j], 'M']<-parSE[j]
			}
		}
	
	## Print some results
	## Print the Model fit
	if (fit){
	cat("Model fit statistics and indices\n")
	print(fitInd)
	}
	## Print the parameter estimates
	A.na<-A
	A.na[A==0]<-NA
	S.na<-S
	S.na[S==0]<-NA
	Ase.na<-Ase
	Ase.na[Ase==0]<-NA
	Sse.na<-Sse
	Sse.na[Sse==0]<-NA
	
	if (ram.out){
	  cat("\n--------------------\n")
	  cat("Parameter estimates:\n")
	  cat("--------------------\n")
	  cat("\nMatrix A\n\n")
	  print(A.na, digits=digits,na.print = zero.print)
	  cat("\nMatrix S\n\n")
	  print(S.na,digits=digits,na.print = zero.print)
	  
	  cat("\n----------------------------------------\n")
	  cat("Standard errors for parameter estimates:\n")
	  cat("----------------------------------------\n")
	  cat("\nMatrix A\n\n")
	  print(Ase.na,digits=digits,na.print = zero.print)
	  cat("\nMatrix S\n\n")
	  print(Sse.na,digits=digits,na.print = zero.print)
		cat("\n\n")
	}
  lname<-NULL
  if (latent>0) lname = varName[(manifest+1):nrow]
	invisible(list(A=A, S=S, Ase=Ase, Sse=Sse,  M=M, Mse=Mse, fit=fitInd, lavaan=fitModel, nvar=nrow, manifest=manifest,latent=latent,lname=lname,varname=varName))
	}
}

## Refit a model
ramReFit<-function(object, add, ram.out=FALSE, ...){   
   data<-object$data
   if (!is.list(object$model)){
		model<-paste(object$model, "\n", add)
		fitModel<-lavaan(model=model, data=data, ...)
		print(summary(fitModel, fit.measures=TRUE))
		invisible((list(model=model, lavaan=fitModel, data=data)))
	}else{
		if (object$type=='no'){
			model<-paste(object$model$no, "\n", add)
			results.no<-growth(model=model, data=data,...)
			print(summary(results.no, fit.measures=TRUE))
			ram.no<-lavaan2ram(results.no, ram.out=ram.out)
			invisible((list(lavaan=list(no=results.no), ram=list(no=ram.no), model=list(no=model), data=data, type=object$type)))
		}else if (object$type=="linear"){
			model<-paste(object$model$linear, "\n", add)
			results.linear<-growth(model=model, data=data,...)
			ram.linear<-lavaan2ram(results.linear, ram.out=ram.out)
			print(summary(results.linear,fit.measures=TRUE))
			invisible((list(lavaan=list(linear=results.linear), ram=list(linear=ram.linear), model=list(linear=model), data=data, type=object$type)))
		}else if (object$type=="latent"){
			model<-paste(object$model$latent, "\n", add)
			results.latent<-growth(model=model, data=data,...)
			ram.latent<-lavaan2ram(results.latent, ram.out=ram.out)
			print(summary(results.latent,fit.measures=TRUE))
			invisible((list(lavaan=list(latent=results.latent), ram=list(latent=ram.latent), model=list(latent=model), data=data, type=object$type)))
		}else if (object$type=="quadratic"){
			model<-paste(object$model$quadratic, "\n", add)
			results.quadratic<-growth(model=model, data=data,...)
			ram.quadratic<-lavaan2ram(results.quadratic, ram.out=ram.out)
			print(summary(results.quadratic,fit.measures=TRUE))
			invisible((list(lavaan=list(quadratic=results.quadratic), ram=list(quadratic=ram.quadratic), model=list(quadratic=model), data=data, type=object$type)))
		}else{
		stop("Refit is not allowed for this type of models!")
		}
	}
}

ramShowModel<-function(object){
	if (!is.list(object$model)){
		cat(object$model)
		cat("\n")
	}else{
		for (i in 1:length(object$model)){
			cat(object$model[[i]])
			cat("\n\n")
		}
	}
}