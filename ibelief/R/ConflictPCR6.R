ConflictPCR6 <- function(MassIn,TabConflict){

	#MassIn: 2^n lines, nb experts columns (n is the number of classes)
	#[TabConflict] = ConflictTable(n,nbexperts);
	#Conflict: 2^n lines, 1 column


	n=nrow(MassIn);
	m=ncol(MassIn);

	if (missing(TabConflict)){     
		# check the number of input arguments
		TabConflict = ConflictTable(n,m);
	}

    Conflict=matrix(0,n,1);
	TabInd=matrix(1,m,1);
	condition=TRUE;


	while(condition){

		indCour=m;
		
		P=1;
		S=0;
		for (i in 1:m){
			P=P*MassIn[TabInd[i],i];
			S=S+MassIn[TabInd[i],i];    
		}
		
			
		if ((S!=0) && (P!=0)){
			cmp=1;
					
			OK=FALSE;
			while (!OK && cmp<ncol(TabConflict)+1){ 
				if(sum(is.element(unique(TabInd),TabConflict[,cmp]))==length(which(TabConflict[,cmp]!=0))){
					OK=TRUE;
				}
				cmp=cmp+1;
			}
			if(OK){
				for(j in 1:m){
					Conflict[TabInd[j]]=Conflict[TabInd[j]]+MassIn[TabInd[j],j]*P/S;
				}
			}
		}
		
		TabInd[indCour]=TabInd[indCour]+1;
		
		
		while(((TabInd[indCour]==n+1) && condition)){
			TabInd[indCour]=1;
			indCour=indCour-1;
			if(indCour==0){
				condition=FALSE;
				indCour=1;
		    }else{
				TabInd[indCour]=TabInd[indCour]+1;
		    }
		}
		
	}

	return(Conflict)
}
