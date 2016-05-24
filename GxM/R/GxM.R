GxM <-
function(dataset,rGvalue=0.995,modelname,zeroset=character(),closedform=FALSE,K=8,coreNumber=1,manualinitial=NULL,priority=1,gradientlevel=2) {

  ## replace rG=1 with rGvalue
  dataset$rG[abs(dataset$rG-1) < 1e-4] = rGvalue;
  listdata = split(dataset, dataset$id);

  ## initialize model parameters
  if(!(modelname %in% c('Chol','CholGxM','NLMainGxM','CorrGxM','CholNonLin','CorrNonLin'))) { 
    cat('Model name should be Chol, CholGxM, NLMainGxM, CorrGxM, CholNonLin or CorrNonLin.');  
    stop(); 
  }

  ## model 3
  if(modelname == 'Chol') { 
    fullparam = c('muM','aM','cM','eM','muP','aC','cC','eC','aU','cU','eU');
    misspell = (!(zeroset %in% fullparam));
    if(sum(misspell)>0) {
      cat('Unexpected input in zeroset:', zeroset[misspell],'\n'); 
      cat('Parameters of model Chol should be member(s) of\n', fullparam,'\n');
      stop();
    }
    usedfullparam = setdiff(fullparam,zeroset);
    
    misspell = (!(names(manualinitial) %in% usedfullparam));
    if(sum(misspell)>0) { cat(names(manualinitial)[misspell], 'can not be manually initialized.\n'); }
    usedmanualinitial = manualinitial[!misspell];
    
    starting = initialization_chol(dataset=dataset,rGvalue=rGvalue,fullparam=fullparam,zeroset=zeroset,manualinitial=usedmanualinitial);
    if(closedform)  { fMP = fMP_Chol_closedform; } else { fMP = fMP_Chol; } 
    GxMmle=GxMestimator(listdata=listdata,fMP=fMP,zeroset=zeroset,closedform=closedform,starting=starting,K=K,coreNumber=coreNumber,usedfullparam=usedfullparam,usedmanualinitial=usedmanualinitial,priority=priority,gradientlevel=gradientlevel);
    GxMmle$par = refineGxM(GxMmle$par, modelname='Chol');
  }

  ## model 4
  if(modelname == 'CholGxM') { 
    fullparam = c('muM','aM','cM','eM','muP','aC','cC','eC','aU','cU','eU','alphaC','alphaU','kappaC','kappaU','epsilonC','epsilonU');
    misspell = (!(zeroset %in% fullparam));
    if(sum(misspell)>0) {
      cat('Unexpected input in zeroset:', zeroset[misspell],'\n'); 
      cat('Parameters of model CholGxM should be member(s) of\n', fullparam,'\n');
      stop();
    }
    usedfullparam = setdiff(fullparam,zeroset);
    
    misspell = (!(names(manualinitial) %in% usedfullparam));
    if(sum(misspell)>0) { cat(names(manualinitial)[misspell], 'can not be manually initialized.\n'); }
    usedmanualinitial = manualinitial[!misspell];
    
    starting = initialization_chol(dataset=dataset,rGvalue=rGvalue,fullparam=fullparam,zeroset=zeroset,manualinitial=usedmanualinitial);
    if(closedform)  { fMP = fMP_CholGxM_closedform; } else { fMP = fMP_CholGxM; } 
    GxMmle=GxMestimator(listdata=listdata,fMP=fMP,zeroset=zeroset,closedform=closedform,starting=starting,K=K,coreNumber=coreNumber,usedfullparam=usedfullparam,usedmanualinitial=usedmanualinitial,priority=priority,gradientlevel=gradientlevel);
    GxMmle$par = refineGxM(GxMmle$par, modelname='CholGxM');
  }

  ## model 5
  if(modelname == 'NLMainGxM') { 
    fullparam = c('muM','aM','cM','eM','muP','beta1','beta2','aU','cU','eU','alphaU','kappaU','epsilonU');
    misspell = (!(zeroset %in% fullparam));
    if(sum(misspell)>0) {
      cat('Unexpected input in zeroset:', zeroset[misspell],'\n'); 
      cat('Parameters of model NLMainGxM should be member(s) of\n', fullparam,'\n');
      stop();
    }
    usedfullparam = setdiff(fullparam,zeroset);
    
    misspell = (!(names(manualinitial) %in% usedfullparam));
    if(sum(misspell)>0) { cat(names(manualinitial)[misspell], 'can not be manually initialized.\n'); }
    usedmanualinitial = manualinitial[!misspell];
    
    starting = initialization_nlmain(dataset=dataset,rGvalue=rGvalue,fullparam=fullparam,zeroset=zeroset,manualinitial=usedmanualinitial);
    if(closedform)  { fMP = fMP_NLMainGxM_closedform; } else { fMP = fMP_NLMainGxM; } 
    GxMmle=GxMestimator(listdata=listdata,fMP=fMP,zeroset=zeroset,closedform=closedform,starting=starting,K=K,coreNumber=coreNumber,usedfullparam=usedfullparam,usedmanualinitial=usedmanualinitial,priority=priority,gradientlevel=gradientlevel);
    GxMmle$par = refineGxM(GxMmle$par, modelname='NLMainGxM');
  }

  ## model 6
  if(modelname == 'CorrGxM') { 
    fullparam = c('muM','aM','cM','eM','muP','aP','cP','eP','rA','rC','rE','alphaP','kappaP','epsilonP');
    misspell = (!(zeroset %in% fullparam));
    if(sum(misspell)>0) {
      cat('Unexpected input in zeroset:', zeroset[misspell],'\n'); 
      cat('Parameters of model CorrGxM should be member(s) of\n', fullparam,'\n');
      stop();
    }
    usedfullparam = setdiff(fullparam,zeroset);
    
    misspell = (!(names(manualinitial) %in% usedfullparam));
    if(sum(misspell)>0) { cat(names(manualinitial)[misspell], 'can not be manually initialized.\n'); }
    usedmanualinitial = manualinitial[!misspell];
    
    starting = initialization_corr(dataset=dataset,rGvalue=rGvalue,fullparam=fullparam,zeroset=zeroset,manualinitial=usedmanualinitial);
    if(closedform)  { fMP = fMP_CorrGxM_closedform; } else { fMP = fMP_CorrGxM; } 
    GxMmle=GxMestimator(listdata=listdata,fMP=fMP,zeroset=zeroset,closedform=closedform,starting=starting,K=K,coreNumber=coreNumber,usedfullparam=usedfullparam,usedmanualinitial=usedmanualinitial,priority=priority,gradientlevel=gradientlevel);
    GxMmle$par = refineGxM(GxMmle$par, modelname='CorrGxM');
  }

  ## model 7
  if(modelname == 'CholNonLin') { 
    if(closedform) { cat('Computation using closed-form likelihood for model CholNonLin is not provided.\n'); stop(); }
    fullparam = c('muM','aM','cM','eM','muP','aC','cC','eC','aU','cU','eU','gamma1','gamma2','gamma3','delta1','delta2','delta3');
    misspell = (!(zeroset %in% fullparam));
    if(sum(misspell)>0) {
      cat('Unexpected input in zeroset:', zeroset[misspell],'\n'); 
      cat('Parameters of model CholNonLin should be member(s) of\n', fullparam,'\n');
      stop();
    }
    usedfullparam = setdiff(fullparam,zeroset);
    
    misspell = (!(names(manualinitial) %in% usedfullparam));
    if(sum(misspell)>0) { cat(names(manualinitial)[misspell], 'can not be manually initialized.\n'); }
    usedmanualinitial = manualinitial[!misspell];
    
    starting = initialization_chol(dataset=dataset,rGvalue=rGvalue,fullparam=fullparam,zeroset=zeroset,manualinitial=usedmanualinitial);
    fMP=fMP_CholNonLin;
    GxMmle=GxMestimator(listdata=listdata,fMP=fMP,zeroset=zeroset,closedform=closedform,starting=starting,K=K,coreNumber=coreNumber,usedfullparam=usedfullparam,usedmanualinitial=usedmanualinitial,priority=priority,gradientlevel=gradientlevel);
    GxMmle$par = refineGxM(GxMmle$par, modelname='CholNonLin');
  }

  ## model 8
  if(modelname == 'CorrNonLin') {
    if(closedform) { cat('Computation using closed-form likelihood for model CorrNonLin is not provided.\n'); stop(); }
    fullparam = c('muM','aM','cM','eM','muP','aP','cP','eP','rA','rC','rE','lambda1','lambda2','lambda3');
    misspell = (!(zeroset %in% fullparam));
    if(sum(misspell)>0) {
      cat('Unexpected input in zeroset:', zeroset[misspell],'\n'); 
      cat('Parameters of model CorrNonLin should be member(s) of\n', fullparam,'\n');
      stop();
    }
    usedfullparam = setdiff(fullparam,zeroset);
    
    misspell = (!(names(manualinitial) %in% usedfullparam));
    if(sum(misspell)>0) { cat(names(manualinitial)[misspell], 'can not be manually initialized.\n'); }
    usedmanualinitial = manualinitial[!misspell];
    
    starting = initialization_corr(dataset=dataset,rGvalue=rGvalue,fullparam=fullparam,zeroset=zeroset,manualinitial=usedmanualinitial);
    fMP = fMP_CorrNonLin;
    GxMmle=GxMestimator(listdata=listdata,fMP=fMP,zeroset=zeroset,closedform=closedform,starting=starting,K=K,coreNumber=coreNumber,usedfullparam=usedfullparam,usedmanualinitial=usedmanualinitial,priority=priority,gradientlevel=gradientlevel);
    GxMmle$par = refineGxM(GxMmle$par, modelname='CorrNonLin');
  }

  return(GxMclass(loglikelihood=GxMmle$loglikelihood, BIC=GxMmle$BIC, par=GxMmle$par, hess=GxMmle$hess, gradient=GxMmle$gradient,
                  modelname=modelname, zeroset=zeroset, closedform=closedform, K=K, coreNumber=coreNumber));
}
