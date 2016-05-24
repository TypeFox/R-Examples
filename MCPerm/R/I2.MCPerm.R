I2.MCPerm <-
function(case_11,case_12,case_22,control_11,control_12,control_22,
	 model="allele",method="MH",repeatNum=1000){
    if (!is.element(model, c("allele", "dominant", "recessive"))){
	     stop("'model' should be 'allele', 'dominant' or 'recessive'.")
	 } 
	 if (!is.element(method, c("Inverse","MH","Peto"))){
	     stop("'method' should be 'Inverse','MH' or 'Peto'.")
	 }
	 if ( repeatNum<0 || repeatNum != round(repeatNum) ) {
         stop("'repeatNum' must be a positive integer.")
    }
	 
	 study_num=length(case_11)
	 #### 1. make sure risk allele
	 case_1=2*case_11+case_12
	 case_2=2*case_22+case_12
	 control_1=2*control_11+control_12
	 control_2=2*control_22+control_12
	 OR=(case_1*control_2)/(case_2*control_1)
	 risk_allele=1
	 not_risk=length(which(OR<1))
	 if(not_risk>(study_num/2)){
	    risk_allele=2
		 temp=case_11
		 case_11=case_22
		 case_22=temp
		 temp=control_11
		 control_11=control_22
		 control_22=temp
	 }
	 
	 #### 2.calculate allele/dominant/recessive
	 switch(model,
	    allele={case_1=case_11*2+case_12
		    case_2=case_22*2+case_12
			 control_1=control_11*2+control_12
			 control_2=control_22*2+control_12
		},
		dominant={case_1=case_11+case_12
			 case_2=case_22
		    control_1=control_11+control_12
			 control_2=control_22
      },
		recessive={case_1=case_11
			 case_2=case_12+case_22
			 control_1=control_11
		    control_2=control_12+control_22
		}
	 )
     
	 #### 3. heterogeneity of real data
	 switch(method,
	    Inverse={ true_meta=rma(ai=case_1, bi=case_2, ci=control_1, di=control_2,
		    measure="OR", method="FE")
		 },
		 MH={ true_meta=rma.mh(ai=case_1, bi=case_2, ci=control_1, di=control_2,
			 measure="OR")
		 },
	    Peto={ true_meta=rma.peto(ai=case_1,bi=case_2,ci=control_1,di=control_2)
		 }
	 )
	 QE=true_meta$QE
	 I2=max(100*(QE-(study_num-1))/QE,0)
	 
	 #### 4. generate random data
	 perm_case_11=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_case_12=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_case_22=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_control_11=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_control_12=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_control_22=matrix(0,nrow=study_num,ncol=repeatNum)
	 
	 count_case=case_11+case_12+case_22
	 count_control=control_11+control_12+control_22
	 sample=count_case+count_control
	 count_11=case_11+control_11
	 count_12=case_12+control_12
	 
	 for(i in 1:study_num){
	    perm_case_11[i,]=rhyper(repeatNum,m=count_case[i],n=count_control[i],k=count_11[i])
		 perm_control_11[i,]=count_11[i]-perm_case_11[i,]
		 temp_count_case=count_case[i]-perm_case_11[i,]
		 temp_count_control=count_control[i]-perm_control_11[i,]
		 perm_case_12[i,]=rhyper(repeatNum,m=temp_count_case,n=temp_count_control,k=count_12[i])
		 perm_control_12[i,]=count_12[i]-perm_case_12[i,]
		 perm_case_22[i,]=temp_count_case-perm_case_12[i,]
		 perm_control_22[i,]=temp_count_control-perm_control_12[i,]
	 }
	 
	 #### 5. heterogeneity for random data
	 switch(model,
		 allele={perm_case_1=2*perm_case_11+perm_case_12
		    perm_case_2=2*perm_case_22+perm_case_12
			 perm_control_1=2*perm_control_11+perm_control_12
			 perm_control_2=2*perm_control_22+perm_control_12
		 },
		 dominant={perm_case_1=perm_case_11+perm_case_12
		    perm_case_2=perm_case_22
			 perm_control_1=perm_control_11+perm_control_12
			 perm_control_2=perm_control_22
		 },
		 recessive={perm_case_1=perm_case_11
		    perm_case_2=perm_case_12+perm_case_22
			 perm_control_1=perm_control_11
			 perm_control_2=perm_control_12+perm_control_22
		 }
	 )
	 perm_I2=matrix(0,nrow=1,ncol=repeatNum)
	 
	 for(i in 1:repeatNum){
	    switch(method,
		    Inverse={temp=rma(ai=perm_case_1[,i], bi=perm_case_2[,i], 
			    ci=perm_control_1[,i], di=perm_control_2[,i], measure="OR", method="FE")
			 },
			 MH={temp=rma.mh(ai=perm_case_1[,i],bi=perm_case_2[,i],
			    ci=perm_control_1[,i],di=perm_control_2[,i],measure="OR")
			 },
			 Peto={temp=rma.peto(ai=perm_case_1[,i],bi=perm_case_2[,i],
			    ci=perm_control_1[,i],di=perm_control_2[,i])
			 }
		 )
		 perm_I2[1,i]=max(100*(temp$QE-(study_num-1))/temp$QE,0)
	 }
	 corrected_I2p=length(which(perm_I2>I2))/repeatNum
	 
	 ### return value
	 result=list("risk_allele"=risk_allele,"I2"=I2,"corrected_I2p"=corrected_I2p)
	 return(result)
}
