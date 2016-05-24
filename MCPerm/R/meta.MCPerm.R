meta.MCPerm <-
function(case_11,case_12,case_22,control_11,control_12,control_22,
	 model="allele",fixed_method="MH",random_method="DL",Qp_alpha=0.01,repeatNum=1000){
    if (!is.element(model, c("allele", "dominant", "recessive"))){
	     stop("'model' should be 'allele', 'dominant' or 'recessive'.")
	 } 
	 if (!is.element(fixed_method, c("Inverse","MH","Peto"))){
	     stop("'fixed_method' should be 'Inverse','MH' or 'Peto'.")
	 }
	 if (!is.element(random_method, c("DL", "HE", "SJ", "ML", "REML", "EB"))){
	     stop("'random_method' should be 'DL', 'HE', 'SJ', 'ML', 'REML' or 'EB'.")
	 }
	 
	 if(length(Qp_alpha)==1){
	     if(Qp_alpha<0 || Qp_alpha>1){
	         stop("'Qp_alpha' should be in 0~1.")
	     }
	 }else{
	     stop("'Qp_alpha' should be in 0~1.")
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
	 switch(fixed_method,
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
	 QEp=true_meta$QEp
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
	 
	 perm_lnOR=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_VARlnOR=matrix(0,nrow=study_num,ncol=repeatNum)
	 perm_Qp=matrix(0,nrow=1,ncol=repeatNum)
	 perm_I2=matrix(0,nrow=1,ncol=repeatNum)
	 
	 perm_merged_LnOR=matrix(0,nrow=1,ncol=repeatNum)
	 perm_merged_VARLnOR=matrix(0,nrow=1,ncol=repeatNum)
	 perm_p=matrix(0,nrow=1,ncol=repeatNum)
	 for(i in 1:repeatNum){
	    switch(fixed_method,
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
		 perm_lnOR[,i]=temp$yi
		 perm_VARlnOR[,i]=temp$vi
		 perm_Qp[1,i]=temp$QEp
		 perm_I2[1,i]=max(100*(temp$QE-(study_num-1))/temp$QE,0)
		 
		 perm_merged_LnOR[1,i]=temp$b
		 perm_merged_VARLnOR[1,i]=(temp$se)*(temp$se)
		 perm_p[1,i]=temp$pval
	 }
	 corrected_Qp=length(which(perm_Qp<QEp))/repeatNum
	 corrected_I2p=length(which(perm_I2>I2))/repeatNum
	 
	 #### 6. based on corrected heterogeneity,
	 #### calculate p value of merged lnOR for true data and random data
	 if(corrected_Qp<Qp_alpha || corrected_Qp<Qp_alpha){
	    true_meta=rma(ai=case_1,bi=case_2,ci=control_1,di=control_2,measure="OR",method=random_method)
	    for(i in 1:repeatNum){
		    temp=rma(yi=perm_lnOR[,i],vi=perm_VARlnOR[,i],method=random_method) 
			 perm_merged_LnOR[1,i]=temp$b
			 perm_merged_VARLnOR[1,i]=(temp$se)*(temp$se)
		    perm_p[1,i]=temp$pval
		 }
	 }
	 corrected_p=length(which(perm_p<true_meta$pval))/repeatNum
	 
	 ### return value
	 corrected_result=matrix(c(QE,I2,true_meta$b,QEp,NA,true_meta$pval,
	    corrected_Qp,corrected_I2p,corrected_p),nrow=3,byrow=TRUE)
	 colnames(corrected_result)=c("Q","I2","merged_lnOR")
	 rownames(corrected_result)=c("true_stat","p.value","p.corrected")
	 result=list("corrected_result"=corrected_result,"risk_allele"=risk_allele,
	    "true_merged_LnOR"=true_meta$b,"true_merged_LnOR_VAR"=(true_meta$se)*(true_meta$se),
		 "true_merged_LnOR_p"=true_meta$pval,"true_merged_LnOR_ci.lb"=true_meta$ci.lb,
		 "true_merged_LnOR_ci.ub"=true_meta$ci.ub,"study_num"=study_num,
	    "true_LnOR"=true_meta$yi[1:study_num],"true_VARLnOR"=true_meta$vi,
	    "perm_case_11"=perm_case_11,"perm_case_12"=perm_case_12,"perm_case_22"=perm_case_22,
		 "perm_control_11"=perm_control_11,"perm_control_12"=perm_control_12,
		 "perm_control_22"=perm_control_22,"perm_LnOR"=perm_lnOR,"perm_VARLnOR"=perm_VARlnOR,
		 "perm_Qp"=perm_Qp,"perm_I2"=perm_I2,"perm_merged_LnOR"=perm_merged_LnOR,
		 "perm_merged_VARLnOR"=perm_merged_VARLnOR,"perm_p"=perm_p,"sample"=sample,
		 "model"=model,"fixed_method"=fixed_method,"random_method"=random_method,
		 "Qp_alpha"=Qp_alpha,"repeatNum"=repeatNum)
	 class(result)='PermMeta'
	 return(result)
}
