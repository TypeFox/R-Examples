compare<-function(x,y,adjusted=FALSE,abseps = 0.000001){
  if(x$iid==FALSE){  
    stop(paste("iid=FALSE (default) have been chosen in the input of timeROC() for estimation of the first object. \n",sep=""))
    }
  if(y$iid==FALSE){  
    stop(paste("iid=FALSE (default) have been chosen in the input of timeROC() for estimation of the second object. \n",sep=""))
  }  
  if(!identical(x$times,y$times)){
    stop(paste("The two objects you want to compare should have been computed for the same vector of times\n",sep=""))
  }
  if(!x$n==y$n){
    stop(paste("The two objects you want to compare were not fitted with the same subjects.\n
               This function does not deal with such a case.",sep=""))
  } 
  # functions to compute adjusted p-values  
  func1<-function(a,nb_stat,MatCorr){
    p<-1-pmvnorm(lower=rep(-abs(a),nb_stat),upper=rep(abs(a),nb_stat),rep(0,nb_stat),sigma=MatCorr,abseps =abseps)
    result<-p[1]
    return(as.numeric(result))
  }
  func2<-function(MatCorr,stat){
    nb_stat<-length(stat)
    p_val<-rep(NA,nb_stat)
    for (i in 1:nb_stat){ 
      p_val[i]<-func1(stat[i],nb_stat=nb_stat,MatCorr=MatCorr)
    }  
    return(p_val)
  }
  func3<-function(mat_contrast,AUCvector,VarcovMat){
    diff<-mat_contrast%*%AUCvector
    MatVarcovdiff<-mat_contrast%*%VarcovMat%*%t(mat_contrast)
    MarCorrDiff<-diag(sqrt(diag(MatVarcovdiff)^(-1))) %*%MatVarcovdiff%*%t(diag(sqrt(diag(MatVarcovdiff)^(-1))))
    stat_test<-diff/sqrt(diag(MatVarcovdiff))
    return(list(diff=diff,
                MatVarcovdiff=MatVarcovdiff,
                MarCorrDiff=MarCorrDiff,
                p_a=func2(MarCorrDiff,stat_test),
                p=func2(matrix(1,length(x$times),length(x$times)),stat_test)  
    )
    )
  }
  if(class(x)=="ipcwcompetingrisksROC"){
    #compute the difference of AUC estimates
    DELTA_1<-x$AUC_1-y$AUC_1
    DELTA_2<-x$AUC_2-y$AUC_2  
    #compute sd of the difference of AUC estimates
    sd_DELTA_1<-apply(x$inference$mat_iid_rep_1-y$inference$mat_iid_rep_1,2,sd)/sqrt(x$n)
    sd_DELTA_2<-apply(x$inference$mat_iid_rep_2-y$inference$mat_iid_rep_2,2,sd)/sqrt(x$n)
    #test statistics
    Test_Stat_1<-abs(DELTA_1)/sd_DELTA_1
    Test_Stat_2<-abs(DELTA_2)/sd_DELTA_2  
    # compute p-values
    p_1<-2*pnorm(Test_Stat_1,lower.tail=FALSE)
    p_2<-2*pnorm(Test_Stat_2,lower.tail=FALSE)
    # compute p-values and adjusted p-values, to account for the multiplicity of tests (as many test statistics as there are times)
    if(adjusted==TRUE){
      if( sum(is.na(Test_Stat_1))>=1 | sum(is.na(Test_Stat_2))>=1){
        stop(paste("You cannot compute adjusted p-values if some AUC estimates are NA.\n
               Please choose the times for which you compute AUCs such as there are not NA values for AUC estimates",sep=""))
       }
      mat_cont<-cbind(diag(1,length(x$times)),diag(-1,length(x$times)))
      var_vect_AUC_1<-var(cbind(x$inference$mat_iid_rep_1,y$inference$mat_iid_rep_1))/(x$n)
      var_vect_AUC_2<-var(cbind(x$inference$mat_iid_rep_2,y$inference$mat_iid_rep_2))/(x$n)
      all_1<-func3(mat_cont,c(x$AUC_1,y$AUC_1),var_vect_AUC_1)
      all_2<-func3(mat_cont,c(x$AUC_2,y$AUC_2),var_vect_AUC_2)
      p_all_1<-rbind(p_1,all_1$p_a)
      rownames(p_all_1)<-c("Non-adjusted","Adjusted")
      p_all_2<-rbind(p_2,all_2$p_a)
      rownames(p_all_2)<-c("Non-adjusted","Adjusted")
      }
    if(adjusted==FALSE){
      return(list(p_values_AUC_1=p_1, p_values_AUC_2=p_2))
      }else{
      return(list(p_values_AUC_1=p_all_1, p_values_AUC_2=p_all_2,Cor=list(Cor_1=all_1$MarCorrDiff,Cor_2=all_2$MarCorrDiff)))
    }
  }
  if(class(x)=="ipcwsurvivalROC"){
    #compute the difference of AUC estimates
    DELTA_1<-x$AUC-y$AUC 
    #compute sd of the difference of AUC estimates
    sd_DELTA_1<-apply(x$inference$mat_iid_rep_1-y$inference$mat_iid_rep_1,2,sd)/sqrt(x$n)
    #test statistics
    Test_Stat_1<-abs(DELTA_1)/sd_DELTA_1
    # compute p-values
    p_1<-2*pnorm(Test_Stat_1,lower.tail=FALSE)
    if(adjusted==TRUE){
      if( sum(is.na(Test_Stat_1))>=1 ){
        stop(paste("You cannot compute adjusted p-values if some AUC estimates are NA.\n
               Please choose the times for which you compute AUCs such as there are not NA values for AUC estimates",sep=""))
      }
      mat_cont<-cbind(diag(1,length(x$times)),diag(-1,length(x$times)))
      var_vect_AUC_1<-var(cbind(x$inference$mat_iid_rep_1,y$inference$mat_iid_rep_1))/(x$n)   
      all_1<-func3(mat_cont,c(x$AUC,y$AUC),var_vect_AUC_1)
      p_all_1<-rbind(p_1,all_1$p_a)
      rownames(p_all_1)<-c("Non-adjusted","Adjusted")
    }
    if(adjusted==FALSE){
      return(list(p_values_AUC=p_1))
    }else{
      return(list(p_values_AUC=p_all_1,Cor=all_1$MarCorrDiff))
    }
  }
}
