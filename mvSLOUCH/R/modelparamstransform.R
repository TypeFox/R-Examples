.par.transform<-function(params,EstimationParams,ModelParams=NULL){
## same parametrization for all models, since all are nested in the mvslouch
## if some parameters are not needed then one uses NA (NOT NULL)
    if (is.null(ModelParams)){
	ModelParams<-vector("list",12)
	names(ModelParams)<-c("A","B","mPsi","mPsi0","vY0","vX0","Syy","Syx","Sxy","Sxx","eigenSignsA","GivensQCsignsA")
    }
    if (!is.null(EstimationParams$Fixed$A)){ModelParams$A<-EstimationParams$Fixed$A}
    else{
	if (EstimationParams$kY==1){ModelParams$A<-matrix(params["A"],ncol=1,nrow=1)}
	else{
	    lAparams=switch(EstimationParams$Atype,
		SingleValueDiagonal={diag(params["A"],ncol=EstimationParams$kY,nrow=EstimationParams$kY)},
		Diagonal={diag(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],ncol=EstimationParams$kY,nrow=EstimationParams$kY)},
		SymmetricPositiveDefinite={.sym.par(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY)},
		Symmetric={.par.transform.symmetric.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY)},
		TwoByTwo={.par.transform.twobytwo.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))])},
		UpperTri={.par.transform.uppertri.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY)},
		LowerTri={.par.transform.lowertri.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY)},
		DecomposablePositive={.par.transform.decomp.pos.real.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY,NULL,ifelse(is.null(EstimationParams$minAeigen),0.01,EstimationParams$minAeigen))},
		DecomposableNegative={.par.transform.decomp.neg.real.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY,NULL,ifelse(is.null(EstimationParams$maxAeigen),0.01,EstimationParams$maxAeigen))},
		DecomposableReal={
		    .par.transform.decomp.real.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY)
		}, 
		Invertible={.par.transform.invert.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY,"qr",ifelse(is.null(EstimationParams$minRdiag),0.01,EstimationParams$minRdiag))}
	    )
	    if (is.element(EstimationParams$Atype,c("DecomposablePositive","DecomposableNegative","DecomposableReal","Invertible"))){
		ModelParams$A<-lAparams$A	    
		ModelParams$GivensQCsignsA<-lAparams$QcSigns
		if (is.element(EstimationParams$Atype,c("DecomposablePositive","DecomposableNegative","DecomposableReal"))){ 
		    ModelParams$eigenSignsA<-lAparams$eigenSigns
		}
	    }
	    else{ModelParams$A<-lAparams}
	}
	if (!is.null(EstimationParams$diagA)){
	    diag(ModelParams$A)=switch(EstimationParams$diagA,
		Positive={exp(diag(ModelParams$A))+ifelse(is.null(EstimationParams$minAdiag),0,EstimationParams$minAdiag)},
		Negative={(-1)*exp(diag(ModelParams$A))-ifelse(is.null(EstimationParams$maxAdiag),0,EstimationParams$maxAdiag)}
	    )
	}	  
	if (!is.null(EstimationParams$signsA)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$A[which(EstimationParams$signsA=="0")]<-0
	    ModelParams$A[which(EstimationParams$signsA==0)]<-0
	    ModelParams$A[which(EstimationParams$signsA=="-")]<- (-1)*exp(ModelParams$A[which(EstimationParams$signsA=="-")])
	    ModelParams$A[which(EstimationParams$signsA=="+")]<- exp(ModelParams$A[which(EstimationParams$signsA=="+")])
	}
	ModelParams$A[which(abs(ModelParams$A)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$B)){ModelParams$B<-EstimationParams$Fixed$B}
    else{

	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){ModelParams$B<-matrix(params["B"],ncol=1,nrow=1)}
	else{
	    ModelParams$B=switch(EstimationParams$Btype,
		MinusA={(-1)*ModelParams$A},
		SingleValue={matrix(params["B"],ncol=EstimationParams$kX,nrow=EstimationParams$kY)},
		SingleValueDiagonal={diag(params["B"],ncol=EstimationParams$kX,nrow=EstimationParams$kY)},
	        Diagonal={diag(params[(which(names(params)=="Bstart")):(which(names(params)=="Bend"))],ncol=EstimationParams$kX,nrow=EstimationParams$kY)},
		Any={matrix(params[(which(names(params)=="Bstart")):(which(names(params)=="Bend"))],ncol=EstimationParams$kX,nrow=EstimationParams$kY,byrow=TRUE)}
	    )
	}
	if (!is.null(EstimationParams$signsB)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$B[which(EstimationParams$signsB=="0")]<-0
	    ModelParams$B[which(EstimationParams$signsB==0)]<-0
	    ModelParams$B[which(EstimationParams$signsB=="-")]<- (-1)*exp(ModelParams$B[which(EstimationParams$signsB=="-")])
	    ModelParams$B[which(EstimationParams$signsB=="+")]<- exp(ModelParams$B[which(EstimationParams$signsB=="+")])
	}
    	ModelParams$B[which(abs(ModelParams$B)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$mPsi)){ModelParams$mPsi<-EstimationParams$Fixed$mPsi}
    else{
	if ((EstimationParams$kY==1)&&((EstimationParams$mPsitype=="Global")||((EstimationParams$mPsitype=="Regimes")&&(length(EstimationParams$RegimeTypes)==1)))){
	    ModelParams$mPsi<-matrix(params["Psi"],ncol=1,nrow=1)
	    if(EstimationParams$mPsitype=="Regimes"){colnames(ModelParams$mPsi)<-EstimationParams$RegimeTypes}
	}
	else{
	    ModelParams$mPsi=switch(EstimationParams$mPsitype,
		Global={matrix(params[(which(names(params)=="Psistart")):(which(names(params)=="Psiend"))],ncol=1)},
		Regimes={
		    mPsi<-matrix(params[(which(names(params)=="Psistart")):(which(names(params)=="Psiend"))],nrow=EstimationParams$kY,ncol=length(EstimationParams$RegimeTypes),byrow=FALSE)
		    names(mPsi)<-EstimationParams$RegimeTypes
		    mPsi
		}	    
	    )    
	}
    	if (!is.null(EstimationParams$signsmPsi)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$mPsi[which(EstimationParams$signsmPsi=="0")]<-0
	    ModelParams$mPsi[which(EstimationParams$signsmPsi==0)]<-0
	    ModelParams$mPsi[which(EstimationParams$signsmPsi=="-")]<- (-1)*exp(ModelParams$mPsi[which(EstimationParams$signsmPsi=="-")])
	    ModelParams$mPsi[which(EstimationParams$signsmPsi=="+")]<- exp(ModelParams$mPsi[which(EstimationParams$signsmPsi=="+")])
	}
    	ModelParams$mPsi[which(abs(ModelParams$mPsi)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$mPsi0)){ModelParams$mPsi0<-EstimationParams$Fixed$mPsi0}
    else{
	if (EstimationParams$kY==1){ModelParams$mPsi0<-params["Psi0"]}
	else{ModelParams$mPsi0<-params[(which(names(params)=="Psi0start")):(which(names(params)=="Psi0end"))]}
    	if (!is.null(EstimationParams$signsmPsi0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="0")]<-0
	    ModelParams$mPsi0[which(EstimationParams$signsmPsi0==0)]<-0
	    ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="-")]<- (-1)*exp(ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="-")])
	    ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="+")]<- exp(ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="+")])
	}
    	ModelParams$mPsi0[which(abs(ModelParams$mPsi0)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$vY0)){ModelParams$vY0<-EstimationParams$Fixed$vY0}
    else{
	if (EstimationParams$kY==1){ModelParams$vY0<-params["vY0"]}
	else{ModelParams$vY0<-params[(which(names(params)=="vY0start")):(which(names(params)=="vY0end"))]}
    	if (!is.null(EstimationParams$signsvY0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$vY0[which(EstimationParams$signsvY0=="0")]<-0
	    ModelParams$vY0[which(EstimationParams$signsvY0==0)]<-0
	    ModelParams$vY0[which(EstimationParams$signsvY0=="-")]<- (-1)*exp(ModelParams$vY0[which(EstimationParams$signsvY0=="-")])
	    ModelParams$vY0[which(EstimationParams$signsvY0=="+")]<- exp(ModelParams$vY0[which(EstimationParams$signsvY0=="+")])
	}
    	ModelParams$vY0[which(abs(ModelParams$vY0)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$vX0)){ModelParams$vX0<-EstimationParams$Fixed$vX0}
    else{
	if (EstimationParams$kX==1){ModelParams$vX0<-params["vX0"]}
	else{ModelParams$vX0<-params[(which(names(params)=="vX0start")):(which(names(params)=="vX0end"))]}
    	if (!is.null(EstimationParams$signsvX0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$vX0[which(EstimationParams$signsvX0=="0")]<-0
	    ModelParams$vX0[which(EstimationParams$signsvX0==0)]<-0
	    ModelParams$vX0[which(EstimationParams$signsvX0=="-")]<- (-1)*exp(ModelParams$vX0[which(EstimationParams$signsvX0=="-")])
	    ModelParams$vX0[which(EstimationParams$signsvX0=="+")]<- exp(ModelParams$vX0[which(EstimationParams$signsvX0=="+")])
	}
    	ModelParams$vX0[which(abs(ModelParams$vX0)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$Syy)){ModelParams$Syy<-EstimationParams$Fixed$Syy}
    else{ ## Generally symmetric structures for these matrices should suffice since we do Ct(C) anyway but just in case general 
	if (EstimationParams$kY==1){ModelParams$Syy<-matrix(params["Syy"],ncol=1,nrow=1)}
	else{
	    ModelParams$Syy=switch(EstimationParams$Syytype,
		SingleValueDiagonal={diag(params["Syy"],ncol=EstimationParams$kY,nrow=EstimationParams$kY)},
		Diagonal={diag(params[(which(names(params)=="Syystart")):(which(names(params)=="Syyend"))],ncol=EstimationParams$kY,nrow=EstimationParams$kY)},
		Symmetric={.sym.par(params[(which(names(params)=="Syystart")):(which(names(params)=="Syyend"))],EstimationParams$kY)},
		UpperTri={.par.transform.uppertri.matrix(params[(which(names(params)=="Syystart")):(which(names(params)=="Syyend"))],EstimationParams$kY)},
		LowerTri={.par.transform.lowertri.matrix(params[(which(names(params)=="Syystart")):(which(names(params)=="Syyend"))],EstimationParams$kY)},
		Any={matrix(params[(which(names(params)=="Syystart")):(which(names(params)=="Syyend"))],ncol=EstimationParams$kY,nrow=EstimationParams$kY,byrow=TRUE)}
	    )	
	}
	if (!is.null(EstimationParams$diagSyy)){
		diag(ModelParams$Syy)=switch(EstimationParams$diagSyy,
		    Positive={exp(diag(ModelParams$Syy))},
		    Negative={(-1)*exp(diag(ModelParams$Syy))}
		)
	}
	if (!is.null(EstimationParams$signsSyy)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$Syy[which(EstimationParams$signsSyy=="0")]<-0
	    ModelParams$Syy[which(EstimationParams$signsSyy==0)]<-0
	    ModelParams$Syy[which(EstimationParams$signsSyy=="-")]<- (-1)*exp(ModelParams$Syy[which(EstimationParams$signsSyy=="-")])
	    ModelParams$Syy[which(EstimationParams$signsSyy=="+")]<- exp(ModelParams$Syy[which(EstimationParams$signsSyy=="+")])
	}
    	ModelParams$Syy[which(abs(ModelParams$Syy)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$Syx)){ModelParams$Syx<-EstimationParams$Fixed$Syx} 
    else{
	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){ModelParams$Syx<-matrix(params["Syx"],ncol=1,nrow=1)}
	else{ModelParams$Syx<-matrix(params[(which(names(params)=="Syxstart")):(which(names(params)=="Syxend"))],ncol=EstimationParams$kX,nrow=EstimationParams$kY,byrow=TRUE)}
    	if (!is.null(EstimationParams$signsSyx)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$Syx[which(EstimationParams$signsSyx=="0")]<-0
	    ModelParams$Syx[which(EstimationParams$signsSyx==0)]<-0
	    ModelParams$Syx[which(EstimationParams$signsSyx=="-")]<- (-1)*exp(ModelParams$Syx[which(EstimationParams$signsSyx=="-")])
	    ModelParams$Syx[which(EstimationParams$signsSyx=="+")]<- exp(ModelParams$Syx[which(EstimationParams$signsSyx=="+")])
	}
    	ModelParams$Syx[which(abs(ModelParams$Syx)<1e-15)]<-0

    }
    if (!(is.null(EstimationParams$Fixed$Sxy))){ModelParams$Sxy<-EstimationParams$Fixed$Sxy} 
    else{
    	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){ModelParams$Sxy<-matrix(params["Sxy"],ncol=1,nrow=1)}
	else{ModelParams$Sxy<-matrix(params[(which(names(params)=="Sxystart")):(which(names(params)=="Sxyend"))],ncol=EstimationParams$kY,nrow=EstimationParams$kX,byrow=TRUE)}    
	if (!is.null(EstimationParams$signsSxy)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$Sxy[which(EstimationParams$signsSxy=="0")]<-0
	    ModelParams$Sxy[which(EstimationParams$signsSxy==0)]<-0
	    ModelParams$Sxy[which(EstimationParams$signsSxy=="-")]<- (-1)*exp(ModelParams$Sxy[which(EstimationParams$signsSxy=="-")])
	    ModelParams$Sxy[which(EstimationParams$signsSxy=="+")]<- exp(ModelParams$Sxy[which(EstimationParams$signsSxy=="+")])
	}
    	ModelParams$Sxy[which(abs(ModelParams$Sxy)<1e-15)]<-0

    }
    if (!is.null(EstimationParams$Fixed$Sxx)){ModelParams$Sxx<-EstimationParams$Fixed$Sxx}
    else{  ## Generally symmetric structures for these matrices should suffice since we do Ct(C) anyway but just in case general 
	if (EstimationParams$kX==1){ModelParams$Sxx<-matrix(params["Sxx"],ncol=1,nrow=1)}
	else{
	    ModelParams$Sxx=switch(EstimationParams$Sxxtype,
		SingleValueDiagonal={diag(params["Sxx"],ncol=EstimationParams$kX,nrow=EstimationParams$kX)},
		Diagonal={diag(params[(which(names(params)=="Sxxstart")):(which(names(params)=="Sxxend"))],ncol=EstimationParams$kX,nrow=EstimationParams$kX)},
		Symmetric={.sym.par(params[(which(names(params)=="Sxxstart")):(which(names(params)=="Sxxend"))],EstimationParams$kX)},
		Any={matrix(params[(which(names(params)=="Sxxstart")):(which(names(params)=="Sxxend"))],ncol=EstimationParams$kX,nrow=EstimationParams$kX,byrow=TRUE)}
	    )
	}
	if (!is.null(EstimationParams$signsSxx)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$Sxx[which(EstimationParams$signsSxx=="0")]<-0
	    ModelParams$Sxx[which(EstimationParams$signsSxx==0)]<-0
	    ModelParams$Sxx[which(EstimationParams$signsSxx=="-")]<- (-1)*exp(ModelParams$Sxx[which(EstimationParams$signsSxx=="-")])
	    ModelParams$Sxx[which(EstimationParams$signsSxx=="+")]<- exp(ModelParams$Sxx[which(EstimationParams$signsSxx=="+")])
	}
	ModelParams$Sxx[which(abs(ModelParams$Sxx)<1e-15)]<-0
    }
    ModelParams    
}

.par.inv.transform<-function(ModelParams,EstimationParams){
## same parametrization for all models, since all are nested in the mvslouch
## if some parameters are not needed then one uses NA (NOT NULL)
    params<-c() ## empty vector with which to start paramtrizing
    if (is.null(EstimationParams$Fixed$A)){
	if (!is.null(EstimationParams$signsA)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$A[which(EstimationParams$signsA=="-")]<- log((-1)*ModelParams$A[which(EstimationParams$signsA=="-")])
		ModelParams$A[which(EstimationParams$signsA=="+")]<- log(ModelParams$A[which(EstimationParams$signsA=="+")])
	}
	if (!is.null(EstimationParams$diagA)){
	    diag(ModelParams$A)=switch(EstimationParams$diagA,
		Positive={log(diag(ModelParams$A)-ifelse(is.null(EstimationParams$minAdiag),0,EstimationParams$minAdiag))},
		Negative={log((-1)*diag(ModelParams$A)-ifelse(is.null(EstimationParams$maxAdiag),0,EstimationParams$maxAdiag))}
	    )
	}
	if (EstimationParams$kY==1){Aparam<-c("A"=ModelParams$A[1,1])}
	else{
	    lAparams=switch(EstimationParams$Atype,
		SingleValueDiagonal<-{ModelParams$A[1,1]},
		SymmetricPositiveDefinite={.sym.unpar(ModelParams$A)},
		Symmetric={.par.inv.transform.symmetric.matrix(ModelParams$A)},
		TwoByTwo={.par.inv.transform.twobytwo.matrix(ModelParams$A)},
		UpperTri={.par.inv.transform.uppertri.matrix(ModelParams$A)},
		LowerTri={.par.inv.transform.lowertri.matrix(ModelParams$A)},
		DecomposablePositive={.par.inv.transform.decomp.pos.real.matrix(ModelParams$A,ModelParams$eigenSignsA,ModelParams$GivensQCsignsA,1e-15,ifelse(is.null(EstimationParams$minAeigen),0.01,EstimationParams$minAeigen))},
		DecomposableNegative={.par.inv.transform.decomp.neg.real.matrix(ModelParams$A,ModelParams$eigenSignsA,ModelParams$GivensQCsignsA,1e-15,ifelse(is.null(EstimationParams$maxAeigen),0.01,EstimationParams$maxAeigen))},
		DecomposableReal={.par.inv.transform.decomp.real.matrix(ModelParams$A,ModelParams$eigenSignsA,ModelParams$GivensQCsignsA)}, 
		Invertible={.par.inv.transform.invert.matrix(ModelParams$A,1e-15,"qr",NULL,ifelse(is.null(EstimationParams$minRdiag),0.01,EstimationParams$minRdiag))}
	    )
	    if (is.element(EstimationParams$Atype,c("DecomposablePositive","DecomposableNegative","DecomposableReal","Invertible"))){
		Aparam<-lAparams$vParams
	    }
	    else{Aparam<-lAparams}		    
	    
	    if((EstimationParams$Atype=="SingleValue")||(EstimationParams$Atype=="SingleValueDiagonal")){Aparam<-c("A"=Aparam)}
	    else{
		names(Aparam)<-sapply(1:length(Aparam),function(x){paste("A_",x,sep="")})
		names(Aparam)[1]<-"Astart"
		names(Aparam)[length(Aparam)]<-"Aend"	 
	    }
	}
	params<-c(params,Aparam)
    }
    
    if (is.null(EstimationParams$Fixed$B)){
	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){params<-c(params,"B"=ModelParams$B[1,1])}
	else{
	    if (!is.null(EstimationParams$signsB)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$B[which(EstimationParams$signsB=="-")]<- log((-1)*ModelParams$B[which(EstimationParams$signsB=="-")])
		ModelParams$B[which(EstimationParams$signsB=="+")]<- log(ModelParams$B[which(EstimationParams$signsB=="+")])
	    }
	    Bparam=switch(EstimationParams$Btype,
		MinusA={c(NA)},
		SingleValue={ModelParams$B[1,1]},
		SingleValueDiagonal<-{ModelParams$B[1,1]},
		Diagonal={diag(ModelParams$B)},
		Symmetric={.sym.unpar(ModelParams$B)},
		Any={c(t(ModelParams$B))}
	    )
	    if (!(is.na(Bparam))){
		if((EstimationParams$Btype=="SingleValue")||(EstimationParams$Btype=="SingleValueDiagonal")){Bparam<-c("B"=Bparam)}
		else{
		    names(Bparam)<-sapply(1:length(Bparam),function(x){paste("B_",x,sep="")})
		    names(Bparam)[1]<-"Bstart"
		    names(Bparam)[length(Bparam)]<-"Bend"	 
		}
		params<-c(params,Bparam)
	    }
	}
    }
    
    if (is.null(EstimationParams$Fixed$mPsi)){
    	if ((EstimationParams$kY==1)&&((EstimationParams$mPsitype=="Global")||((EstimationParams$mPsitype=="Regimes")&&(length(EstimationParams$RegimeTypes)==1)))){
	    Psiparam<-c("Psi"=ModelParams$mPsi[,1])
	}else{
	    if (!is.null(EstimationParams$signsmPsi)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$mPsi[which(EstimationParams$signsmPsi=="-")]<- log((-1)*ModelParams$mPsi[which(EstimationParams$signsmPsi=="-")])
		ModelParams$mPsi[which(EstimationParams$signsmPsi=="+")]<- log(ModelParams$mPsi[which(EstimationParams$signsmPsi=="+")])
	    }
	    Psiparam=switch(EstimationParams$mPsitype,
		Global={ModelParams$mPsi[,1]},
		Regimes={c(ModelParams$mPsi)}
	    )    
	    names(Psiparam)<-sapply(1:length(Psiparam),function(x){paste("Psi_",x,sep="")})
	    names(Psiparam)[1]<-"Psistart"
	    names(Psiparam)[length(Psiparam)]<-"Psiend"	 
	}
	params<-c(params,Psiparam)
    }
    if (is.null(EstimationParams$Fixed$mPsi0)){
    	if (EstimationParams$kY==1){Psi0param<-c("Psi0"=ModelParams$mPsi0)}
	else{
	    if (!is.null(EstimationParams$signsmPsi0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="-")]<- log((-1)*ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="-")])
		ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="+")]<- log(ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="+")])
	    }
    	    Psi0param<-ModelParams$mPsi0
	    names(Psi0param)<-sapply(1:length(Psi0param),function(x){paste("Psi0_",x,sep="")})
	    names(Psi0param)[1]<-"Psi0start"
	    names(Psi0param)[length(Psi0param)]<-"Psi0end"	 
	}
	params<-c(params,Psi0param)
    }
    if (is.null(EstimationParams$Fixed$vY0)){
    	if (EstimationParams$kY==1){vY0param<-c("vY0"=ModelParams$vY0)}
	else{
	    if (!is.null(EstimationParams$signsvY0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$vY0[which(EstimationParams$signsvY0=="-")]<- log((-1)*ModelParams$vY0[which(EstimationParams$signsvY0=="-")])
		ModelParams$vY0[which(EstimationParams$signsvY0=="+")]<- log(ModelParams$vY0[which(EstimationParams$signsvY0=="+")])
	    }
    	    vY0param<-ModelParams$vY0
	    names(vY0param)<-sapply(1:length(vY0param),function(x){paste("vY0_",x,sep="")})
	    names(vY0param)[1]<-"vY0start"
	    names(vY0param)[length(vY0param)]<-"vY0end"	 
	}
	params<-c(params,vY0param)
    }
    if (is.null(EstimationParams$Fixed$vX0)){
    	if (EstimationParams$kX==1){vX0param<-c("vX0"=ModelParams$vX0)}
	else{
	    if (!is.null(EstimationParams$signsvX0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$vX0[which(EstimationParams$signsvX0=="-")]<- log((-1)*ModelParams$vX0[which(EstimationParams$signsvX0=="-")])
		ModelParams$vX0[which(EstimationParams$signsvX0=="+")]<- log(ModelParams$vX0[which(EstimationParams$signsvX0=="+")])
	    }
	    vX0param<-ModelParams$vX0
	    names(vX0param)<-sapply(1:length(vX0param),function(x){paste("vX0_",x,sep="")})
	    names(vX0param)[1]<-"vX0start"
	    names(vX0param)[length(vX0param)]<-"vX0end"	 
	}
	params<-c(params,vX0param)
    }
    
    if (is.null(EstimationParams$Fixed$Syy)){ 
    ## Generally symmetric structures for these matrices should suffice since we do Ct(C) anyway but just in case general 
	if (!is.null(EstimationParams$signsSyy)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$Syy[which(EstimationParams$signsSyy=="-")]<- log((-1)*ModelParams$Syy[which(EstimationParams$signsSyy=="-")])
		ModelParams$Syy[which(EstimationParams$signsSyy=="+")]<- log(ModelParams$Syy[which(EstimationParams$signsSyy=="+")])
	}    
        if (!is.null(EstimationParams$diagSyy)){
		diag(ModelParams$Syy)=switch(EstimationParams$diagSyy,
		    Positive={log(diag(ModelParams$Syy))},
		    Negative={log((-1)*diag(ModelParams$Syy))}
		)
	}
	if (EstimationParams$kY==1){Syyparam<-c("Syy"=ModelParams$Syy[1,1])}
	else{	
	    Syyparam=switch(EstimationParams$Syytype,
		SingleValueDiagonal={ModelParams$Syy[1,1]},
		Diagonal={diag(ModelParams$Syy)},
		Symmetric={.sym.unpar(ModelParams$Syy)},
		UpperTri={.par.inv.transform.uppertri.matrix(ModelParams$Syy)},
		LowerTri={.par.inv.transform.lowertri.matrix(ModelParams$Syy)},
		Any={c(t(ModelParams$Syy))}
	    )
	    if (EstimationParams$Syytype!="OneSigmaDiagonal"){
		names(Syyparam)<-sapply(1:length(Syyparam),function(x){paste("Syy_",x,sep="")})
		names(Syyparam)[1]<-"Syystart"
		names(Syyparam)[length(Syyparam)]<-"Syyend"	 
	    }else{Syyparam<-c("Syy"=Syyparam)}
	}
	params<-c(params,Syyparam)	
    }
    
    if (is.null(EstimationParams$Fixed$Syx)){
    	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){Syxparam<-c("Syx"=ModelParams$Syx[1,1])}
	else{
	    if (!is.null(EstimationParams$signsSyx)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$Syx[which(EstimationParams$signsSyx=="-")]<- log((-1)*ModelParams$Syx[which(EstimationParams$signsSyx=="-")])
		ModelParams$Syx[which(EstimationParams$signsSyx=="+")]<- log(ModelParams$Syx[which(EstimationParams$signsSyx=="+")])
	    }
	    Syxparam<-c(t(ModelParams$Syx))
	    names(Syxparam)<-sapply(1:length(Syxparam),function(x){paste("Syx_",x,sep="")})
	    names(Syxparam)[1]<-"Syxstart"
    	    names(Syxparam)[length(Syxparam)]<-"Syxend"	 
	}
	params<-c(params,Syxparam)	
    }

    if (is.null(EstimationParams$Fixed$Sxy)){
	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){Sxytparam<-c("Sxy"=ModelParams$Sxy[1,1])}
	else{
	    if (!is.null(EstimationParams$signsSxy)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$Sxy[which(EstimationParams$signsSxy=="-")]<- log((-1)*ModelParams$Sxy[which(EstimationParams$signsSxy=="-")])
		ModelParams$Sxy[which(EstimationParams$signsSxy=="+")]<- log(ModelParams$Sxy[which(EstimationParams$signsSxy=="+")])
	    }
	    Sxyparam<-c(t(ModelParams$Sxy))
	    names(Sxyparam)<-sapply(1:length(Sxyparam),function(x){paste("Sxy_",x,sep="")})
	    names(Sxyparam)[1]<-"Sxystart"
    	    names(Sxyparam)[length(Sxyparam)]<-"Sxyend"	 
	}
	params<-c(params,Sxyparam)	
    }

    if (is.null(EstimationParams$Fixed$Sxx)){ 
    ## Generally symmetric structures for these matrices should suffice since we do Ct(C) anyway but just in case general 
	if (EstimationParams$kX==1){Sxxparam<-c("Sxx"=ModelParams$Sxx[1,1])}
	else{
	    if (!is.null(EstimationParams$signsSxx)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$Sxx[which(EstimationParams$signsSxx=="-")]<- log((-1)*ModelParams$Sxx[which(EstimationParams$signsSxx=="-")])
		ModelParams$Sxx[which(EstimationParams$signsSxx=="+")]<- log(ModelParams$Sxx[which(EstimationParams$signsSxx=="+")])
	    }
	    Sxxparam=switch(EstimationParams$Sxxtype,
		SingleValueDiagonal={ModelParams$Sxx[1,1]},
		Diagonal={diag(ModelParams$Sxx)},
		Symmetric={.sym.unpar(ModelParams$Sxx)},
		Any={c(t(ModelParams$Sxx))}
	    )
	    if (EstimationParams$Sxxtype!="OneSigmaDiagonal"){
		names(Sxxparam)<-sapply(1:length(Sxxparam),function(x){paste("Sxx_",x,sep="")})
		names(Sxxparam)[1]<-"Sxxstart"
		names(Sxxparam)[length(Sxxparam)]<-"Sxxend"	 
	    }else{Sxxparam<-c("Sxx"=Sxxparam)}
	}
	params<-c(params,Sxxparam)	
    }
    params    
}

