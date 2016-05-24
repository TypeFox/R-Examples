Q.TradPerm <-
function(genotypeData,affectionData,split,sep,naString,model="allele",method="MH",repeatNum=1000)
{
    if (!is.element(model, c("allele", "dominant", "recessive"))){
	     stop("'model' should be 'allele', 'dominant' or 'recessive'.")
	 } 
	 if (!is.element(method, c("Inverse","MH","Peto"))){
	     stop("'method' should be 'Inverse','MH' or 'Peto'.")
	 }
	 if ( repeatNum<0 || repeatNum != round(repeatNum) ) {
         stop("'repeatNum' must be a positive integer.")
    }
	 #### 1. argument
	 ## genotypeData[n*1,string]
	 ## affectionData[n*1,string]
	 ## split=string split
	 ## sep=allele/allele
	 
	 #### 2. calculate frequency for real data
	 study_num=nrow(genotypeData)
	 stat=matrix(0,nrow=study_num,ncol=6)
	 sample=c()
	 ## colnames(stat)=c(case_11,case_12,case_22,control_11,control_12,control_22)
	 for(i in 1:study_num){
	    gen=as.matrix(unlist(strsplit(genotypeData[i,],split=split)),nrow=1)
	    aff=as.matrix(unlist(strsplit(affectionData[i,],split=split)),nrow=1)
	    temp=genotypeStat(gen,aff,1,naString,sep)$genotypeCount
		 stat[i,]=temp[-c(4,8)]
		 sample=c(sample,sum(stat[i,]))
	 }
	 
	 #### 3.make sure risk_allele(OR>1)
    case_1=2*stat[,1]+stat[,2]
	 case_2=2*stat[,3]+stat[,2]
	 control_1=2*stat[,4]+stat[,5]
	 control_2=2*stat[,6]+stat[,5]	
    OR=(case_1*control_2)/(case_2*control_1)
	 genotype=unique(gen)
	 allele12=c()
	 for(i in genotype){
	    allele12=c(allele12,unlist(strsplit(i,split=sep)))
	 }
	 risk_allele=sort(unique(allele12))[1]
	 risk_index=1
	 not_risk=length(which(OR<1))
	 if(not_risk>(study_num/2)){
	    risk_allele=sort(unique(allele12))[2]
		 risk_index=2
		 stat=stat[,c(3,2,1,6,5,4)]
	 }
	 
	 #### 4.calculate allele/dominant/recessive model
	 switch(model,
		 allele={case_1=2*stat[,1]+stat[,2]
		    case_2=2*stat[,3]+stat[,2]
		    control_1=2*stat[,4]+stat[,5]
		    control_2=2*stat[,6]+stat[,5]
		 },
		 dominant={case_1=stat[,1]+stat[,2]
		    case_2=stat[,3]
		    control_1=stat[,4]+stat[,5]
		    control_2=stat[,6]
		 },
		 recessive={case_1=stat[,1]
		    case_2=stat[,2]+stat[,3]
		    control_1=stat[,4]
		    control_2=stat[,5]+stat[,6]
		 }
	 )
	 
	 ### 5. calculate heterogeneity for real data 
	 switch(method,
		 Inverse={true_meta=rma(ai=case_1, bi=case_2, ci=control_1, di=control_2,
			 measure="OR", method="FE")},
		 MH={true_meta=rma.mh(ai=case_1,bi=case_2,ci=control_1,di=control_2,
			 measure="OR")},
		 Peto={true_meta=rma.peto(ai=case_1,bi=case_2,ci=control_1,di=control_2)}
    )
	 QE=true_meta$QE
	 QEp=true_meta$QEp
	 
	 ### 6. generate random data and calculate heterogeneity for random data
	 perm_case_11=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_case_12=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_case_22=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_control_11=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_control_12=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_control_22=matrix(0,nrow=study_num,ncol=repeatNum)
	 
	 perm_Qp=matrix(0,nrow=1,ncol=repeatNum)
	 
	 for(i in 1:repeatNum){
	    ## 6.1. generate random data and calculate frequency for genotypes
	    for(j in 1:study_num){  
		    gen=unlist(strsplit(genotypeData[j,],split=split))
			 aff=unlist(strsplit(affectionData[j,],split=split))
			 ## temp=genotypeStat(gen,aff,1,naString,sep)$genotypeCount;
			 naIndex=which(gen==naString)
			 if(length(naIndex)!=0){
			    gen=gen[-naIndex]
				 aff=aff[-naIndex]
			 }
			 temp=length(gen)
			 temp=sample(temp)
			 gen=gen[temp]
			 perm_stat=table(aff,gen)
		    rowNum=nrow(perm_stat)
		    if(rowNum != 2){
		       stop("no case/control samples\n")
		    }
		 
		    ## deal with genotype<3
		    colNum=ncol(perm_stat)
		    if(colNum != 3){
		       colName=colnames(perm_stat)
			    alleleName=c()
	          for(p in 1:colNum){
	             alleleName=c(alleleName,unlist(strsplit(colName[p], sep)))
	          }
	          temp=matrix(0,nrow=2,ncol=1)
			    
				 if(colNum==2){
				    if (alleleName[1] != alleleName[2]){
			             perm_stat=cbind(temp,perm_stat)
		          }else{
			          if (alleleName[3] != alleleName[4]){
				          perm_stat=cbind(perm_stat,temp)
			          }else{
                      perm_stat=cbind(perm_stat[,1,drop=FALSE],temp,perm_stat[,2,drop=FALSE])
			          }
		          }
				 }
				 
				 if(colNum==1){
				    if (alleleName[1] != alleleName[2]){
			          perm_stat=cbind(temp,perm_stat,temp)
		          }else{
		             perm_stat=cbind(perm_stat,temp,temp)
		          }
				 }
		    }
			 
			 ## risk allele
			 if(risk_index==2){
			    perm_stat=perm_stat[,c(3,2,1)]
			 }
			 perm_case_11[j,i]=perm_stat[2,1]
			 perm_case_12[j,i]=perm_stat[2,2]
			 perm_case_22[j,i]=perm_stat[2,3]
			 perm_control_11[j,i]=perm_stat[1,1]
			 perm_control_12[j,i]=perm_stat[1,2]
			 perm_control_22[j,i]=perm_stat[1,3]
		}
		 ## 6.2 calculate heterogeneity for random data
		 switch(model,
		    allele={perm_case_1=2*perm_case_11[,i]+perm_case_12[,i]
		       perm_case_2=2*perm_case_22[,i]+perm_case_12[,i]
		       perm_control_1=2*perm_control_11[,i]+perm_control_12[,i]
		       perm_control_2=2*perm_control_22[,i]+perm_control_12[,i]
			 },
			 dominant={perm_case_1=perm_case_11[,i]+perm_case_12[,i]
		       perm_case_2=perm_case_22[,i]
		       perm_control_1=perm_control_11[,i]+perm_control_12[,i]
		       perm_control_2=perm_control_22[,i]
			 },
			 recessive={perm_case_1=perm_case_11[,i]
		       perm_case_2=perm_case_12[,i]+perm_case_22[,i]
		       perm_control_1=perm_control_11[,i]
		       perm_control_2=perm_control_12[,i]+perm_control_22[,i]
			 }
		)
		 
		 switch(method,
		    Inverse={temp=rma(ai=perm_case_1, bi=perm_case_2, ci=perm_control_1, di=perm_control_2,
			    measure="OR", method="FE")},
			 MH={temp=rma.mh(ai=perm_case_1,bi=perm_case_2,ci=perm_control_1,di=perm_control_2,
			    measure="OR")},
			 Peto={temp=rma.peto(ai=perm_case_1,bi=perm_case_2,ci=perm_control_1,di=perm_control_2)}
		)
		 perm_Qp[1,i]=temp$QEp
	}
	 
	 ## 7 correct heterogeneity
	 corrected_Qp=length(which(perm_Qp<QEp))/repeatNum
	 
	 ## return value
	 result=list("risk_allele"=risk_allele,"Q"=QE,"Qp"=QEp,"corrected_Qp"=corrected_Qp)
	 return(result)
}
