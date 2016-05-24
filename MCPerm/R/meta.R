meta <-
function( model,fixed_method,random_method,Qp_alpha,
     case_11,case_12,case_22,control_11,control_12,control_22,label=NULL,dataset=NULL){
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
	 
	 if(is.null(dataset)){
	     if(!is.vector(case_11) || !is.vector(case_22) || !is.vector(case_12)){
		     stop("'case_11' or 'case_12' or 'case_22' is not a vector.")
		 }
	     if(!is.vector(control_11) || !is.vector(control_22) || !is.vector(control_12)){
		     stop("'control_11' or 'control_12' or 'control_22' is not a vector.")
		 }
		 if(is.null(label)){
	         label=paste("study",1:length(case_11),sep=" ")
	     }
		 if(length(label)!=length(case_11)){
			 stop("Study labels not of same length as data.")
		 }
	 }else{
	     case_11=dataset[,case_11]
	     case_12=dataset[,case_12]
		 case_22=dataset[,case_22]
		 control_11=dataset[,control_11]
		 control_12=dataset[,control_12]
		 control_22=dataset[,control_22]
		 
		 if(!is.null(label)){
		     len=length(label)
			 if(len==1){
			     label=as.character(dataset[,label])
			 }else{
			     if(len==2){
			         label=paste(dataset[,label[1]], dataset[,label[2]], sep=",")
			     }else{
			         temp=paste(dataset[,label[1]],dataset[,label[2]], sep=",")
				     for(i in 3:len){
				         temp=paste(temp,dataset[,label[i]],sep=",")
				     }
				     label=temp
			     }
			 } 
		 }else{
	         label=paste("study",1:length(case_11),sep=" ")
	     }    		 
	 } 
	 
	 if(model=="allele"){
	     case_1=2*case_11+case_12
		 case_2=2*case_22+case_12
		 control_1=2*control_11+control_12
		 control_2=2*control_22+control_12
	 }
	 if(model=="dominant"){
	     case_1=case_11+case_12
		 case_2=case_22
		 control_1=control_11+control_12
		 control_2=control_22
	 }
	 if(model=="recessive"){
	     case_1=case_11
		 case_2=case_12+case_22
		 control_1=control_11
		 control_2=control_12+control_22
	 }

     ## heterogeneity
	 if(fixed_method=="Inverse"){
	     fixed_result=rma(ai=case_1, bi=case_2, ci=control_1, di=control_2,measure="OR", method="FE",slab=label)
	 }else{
	     if(fixed_method=="MH"){
		     fixed_result=rma.mh(ai=case_1,bi=case_2,ci=control_1,di=control_2,measure="OR",slab=label)
		 }else{
		     fixed_result=rma.peto(ai=case_1,bi=case_2,ci=control_1,di=control_2,slab=label)
		 }
	 }
	 
	 Qp_value=fixed_result$QEp
	 
	 random_result=NULL
	 if(Qp_value<Qp_alpha || Qp_value==Qp_alpha){
		     random_result=rma(ai=case_1,bi=case_2,ci=control_1,di=control_2,measure="OR",
			     method=random_method,slab=label)
			 return(random_result)
	 }else{
	     return(fixed_result)
	 }
}
