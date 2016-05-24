best.reg <-
function(x){
 
  est.reg1<-function(x,col_indep){
    ##creating a fill matrix
    cols<-ncol(x);rows<-nrow(x);estimated_matrix<-matrix(ncol=cols,nrow=rows)
    ##defining dependent and independent variables
    if (col_indep==1){deps<-2:cols} else 
      if (col_indep==cols){deps<-1:(cols-1)} else {
        deps<-c(1:(col_indep-1),(col_indep+1):cols)} 
    ndeps<-length(deps)
    indep<-x[,col_indep]
    rcoefs<-numeric()
    #deciding which variable to regress indep on
    for (i in 1:ndeps){
      vari<-deps[i]
      lm_fit<-lm(log10(x[,vari])~log10(indep)); lm_sum<-summary.lm(lm_fit)
      rcoefs[i]<-sqrt(lm_sum$r.squared)
    }
    ranks<-rank(rcoefs)
    newindep<-indep
    ## run lms and fill in missing values for independent variable
    for (m in 1:length(deps)){
      a<-m-1
      whichvari<-ifelse(ranks==(length(rcoefs)-a),1,0); strongest<-sum(deps*whichvari)
      indeps_lm<-lm(log10(newindep)~log10(x[,strongest]))
      indeps_coef<-indeps_lm$coefficients
      logestimate_indep<-indeps_coef[1]+indeps_coef[2]*log10(x[,strongest])
      estimate_indep<-10^logestimate_indep
      missings<-ifelse(is.na(newindep),1,0)
      nonmissings<-ifelse(is.na(newindep),0,1)
      fillnonmissing<-ifelse(is.na(newindep),0,newindep)
      fillmissing<-ifelse(is.na(newindep),estimate_indep,0)
      newindep<-fillmissing+fillnonmissing
    }
    return(newindep)
  }
  
  
  estimated.matrix<-matrix(ncol=ncol(x),nrow=nrow(x))
  for (i in 1:ncol(x)){
    estimated.matrix[,i]<-est.reg1(x,i)
  }
  return(estimated.matrix)
}
