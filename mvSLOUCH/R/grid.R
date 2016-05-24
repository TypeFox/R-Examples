.GenerateModelGrid<-function(EstimationParams){
    MinMax<-matrix(0,ncol=2,nrow=1)
    i<-1
    step=c()
    if (is.null(EstimationParams$Fixed$A)){
	MinMax<-rbind(MinMax,EstimationParams$Aminmax)
	step<-c(step,EstimationParams$Astep)
	if((EstimationParams$kY==1)||(EstimationParams$Atype=="SingleValue")||(EstimationParams$Atype=="SingleValueDiagonal")){i<-i+1;rownames(MinMax)[i]<-"A"}
	else{
	    rownames(MinMax)[(i+1):(i+nrow(EstimationParams$Aminmax))]<-sapply(1:nrow(EstimationParams$Aminmax),function(x){paste("A_",x,sep="")})
	    rownames(MinMax)[i+1]<-"Astart"
    	    rownames(MinMax)[i+nrow(EstimationParams$Aminmax)]<-"Aend"	 
	    i<-i+nrow(EstimationParams$Aminmax)
	}
    }
    if (is.null(EstimationParams$Fixed$B)){
	MinMax<-rbind(MinMax,EstimationParams$Bminmax)
	step<-c(step,EstimationParams$Bstep)
    	if (((EstimationParams$kX==1)&&(EstimationParams$kY==1))||(EstimationParams$Btype=="SingleValue")||(EstimationParams$Btype=="SingleValueDiagonal")){i<-i+1;rownames(MinMax)[i]<-"B"}
	else{
	    rownames(MinMax)[(i+1):(i+nrow(EstimationParams$Bminmax))]<-sapply(1:nrow(EstimationParams$Bminmax),function(x){paste("B_",x,sep="")})
	    rownames(MinMax)[i+1]<-"Bstart"
    	    rownames(MinMax)[i+nrow(EstimationParams$Bminmax)]<-"Bend"	 
	    i<-i+nrow(EstimationParams$Bminmax)	
	}
    }
    if (is.null(EstimationParams$Fixed$mPsi)){
    	MinMax<-rbind(MinMax,EstimationParams$Psiminmax)
	step<-c(step,EstimationParams$Psistep)    
    	if ((EstimationParams$kY==1)&&((EstimationParams$mPsitype=="Global")||((EstimationParams$mPsitype=="Regimes")&&(length(EstimationParams$RegimeTypes)==1))))
	{i<-i+1;rownames(MinMax)[i]<-"Psi"}
	else{
	    rownames(MinMax)[(i+1):(i+nrow(EstimationParams$Psiminmax))]<-sapply(1:nrow(EstimationParams$Psiminmax),function(x){paste("Psi_",x,sep="")})
	    rownames(MinMax)[i+1]<-"Psistart"
    	    rownames(MinMax)[i+nrow(EstimationParams$Psiminmax)]<-"Psiend"	 
	    i<-i+nrow(EstimationParams$Psiminmax)	
	}
    }
    if (is.null(EstimationParams$Fixed$vY0)){
	MinMax<-rbind(MinMax,EstimationParams$vY0minmax)
	step<-c(step,EstimationParams$vY0step)    
    	if (EstimationParams$kY==1){i<-i+1;rownames(MinMax)[i]<-"vY0"}
	else{
	    rownames(MinMax)[(i+1):(i+nrow(EstimationParams$vY0minmax))]<-sapply(1:nrow(EstimationParams$vY0minmax),function(x){paste("vY0_",x,sep="")})
	    rownames(MinMax)[i+1]<-"vY0start"
    	    rownames(MinMax)[i+nrow(EstimationParams$vY0minmax)]<-"vY0end"	 
	    i<-i+nrow(EstimationParams$vY0minmax)	
	}
    }
    if (is.null(EstimationParams$Fixed$vX0)){
    	MinMax<-rbind(MinMax,EstimationParams$vX0minmax)
	step<-c(step,EstimationParams$vX0step)    
	if (EstimationParams$kX==1){i<-i+1;rownames(MinMax)[i]<-"vX0"}
	else{
	    rownames(MinMax)[(i+1):(i+nrow(EstimationParams$vX0minmax))]<-sapply(1:nrow(EstimationParams$vX0minmax),function(x){paste("vX0_",x,sep="")})
	    rownames(MinMax)[i+1]<-"vX0start"
    	    rownames(MinMax)[i+nrow(EstimationParams$vX0minmax)]<-"vX0end"	 
	    i<-i+nrow(EstimationParams$vX0minmax)	
	}
    }    
    if (is.null(EstimationParams$Fixed$Syy)){ 
    	MinMax<-rbind(MinMax,EstimationParams$Syyminmax)
	step<-c(step,EstimationParams$Syystep)    
    	if ((EstimationParams$kY==1)||(EstimationParams$Syytype=="OneSigmaDiagonal")){i<-i+1;rownames(MinMax)[i]<-"Syy"}
	else{
	    rownames(MinMax)[(i+1):(i+nrow(EstimationParams$Syyminmax))]<-sapply(1:nrow(EstimationParams$Syyminmax),function(x){paste("Syy_",x,sep="")})
	    rownames(MinMax)[i+1]<-"Syystart"
    	    rownames(MinMax)[i+nrow(EstimationParams$Syyminmax)]<-"Syyend"	 
	    i<-i+nrow(EstimationParams$Syyminmax)	
	}
    }
    if (is.null(EstimationParams$Fixed$Syx)){
    	MinMax<-rbind(MinMax,EstimationParams$Syxminmax)
	step<-c(step,EstimationParams$Syxstep)    
    	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){i<-i+1;rownames(MinMax)[i]<-"Syx"}
	else{
	    rownames(MinMax)[(i+1):(i+nrow(EstimationParams$Syxminmax))]<-sapply(1:nrow(EstimationParams$Syxminmax),function(x){paste("Syx_",x,sep="")})
	    rownames(MinMax)[i+1]<-"Syxstart"
    	    rownames(MinMax)[i+nrow(EstimationParams$Syxminmax)]<-"Syxend"	 
	    i<-i+nrow(EstimationParams$Syxminmax)		
	}
    }
    if (is.null(EstimationParams$Fixed$Sxy)){
    	MinMax<-rbind(MinMax,EstimationParams$Sxyminmax)
	step<-c(step,EstimationParams$Sxystep)    
    	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){i<-i+1;rownames(MinMax)[i]<-"Sxy"}
	else{
	    rownames(MinMax)[(i+1):(i+nrow(EstimationParams$Sxyminmax))]<-sapply(1:nrow(EstimationParams$Sxyminmax),function(x){paste("Sxy_",x,sep="")})
	    rownames(MinMax)[i+1]<-"Sxystart"
    	    rownames(MinMax)[i+nrow(EstimationParams$Sxyminmax)]<-"Sxyend"	 
	    i<-i+nrow(EstimationParams$Sxyminmax)		
	}
    }
    if (is.null(EstimationParams$Fixed$Sxx)){ 
    	MinMax<-rbind(MinMax,EstimationParams$Sxxminmax)
	step<-c(step,EstimationParams$Sxxstep)    
    	if ((EstimationParams$kX==1)||(EstimationParams$Syytype=="OneSigmaDiagonal")){i<-i+1;rownames(MinMax)[i]<-"Sxx"}
	else{
	    rownames(MinMax)[(i+1):(i+nrow(EstimationParams$Sxxminmax))]<-sapply(1:nrow(EstimationParams$Sxxminmax),function(x){paste("Sxx_",x,sep="")})
	    rownames(MinMax)[i+1]<-"Sxxstart"
    	    rownames(MinMax)[i+nrow(EstimationParams$Sxxminmax)]<-"Sxxend"	 
	    i<-i+nrow(EstimationParams$Sxxminmax)	
	}
    }    
    MinMax<-MinMax[-1,]
    .generategrid(MinMax,step)        
}


.generategrid<-function(MinMax,step){
## The rownames of MinMax have to be named this is crucial later for parameter identifiability
    Numvars<-nrow(MinMax)
    lValues<-sapply(1:Numvars,function(x){seq(MinMax[x,1],MinMax[x,2],by=step[x])},simplify=FALSE)
    vGridMargSize<-sapply(lValues,function(x){length(x)})
    mGrid<-matrix(NA,ncol=Numvars,nrow=prod(vGridMargSize))
    colnames(mGrid)<-rownames(MinMax)
    mGrid[,1]<-sapply(lValues[[1]],function(x){rep(x,prod(vGridMargSize[-1]))})
    if(Numvars>1){
        mGrid[,2:Numvars]<-sapply(2:Numvars,function(i){rep(sapply(lValues[[i]],function(x){rep(x,prod(vGridMargSize[-(1:i)]))}),prod(vGridMargSize[1:(i-1)]))})
    }
    mGrid
}
