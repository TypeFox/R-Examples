checkGxM <-
function(GxMmle,dataset,rGvalue=0.995,localpercentage=2,localrange=20,outlyingextent=2) {

  ## replace rG=1 with rGvalue
  dataset$rG[abs(dataset$rG-1) < 1e-4] = rGvalue;
  listdata = split(dataset, dataset$id);

  zeroset = GxMmle@zeroset;
  closedform = GxMmle@closedform;
  K = GxMmle@K; 
  coreNumber = GxMmle@coreNumber;

  GxMmle@gradient[is.na(GxMmle@gradient)] = 0;

  ## model 3
  if(GxMmle@modelname == 'Chol') { 
    fullparam = c('muM','muP','aM','cM','eM','aC','cC','eC','aU','cU','eU');
    checkparam = list(muM=0,muP=0,aM=0,cM=0,eM=0,aC=0,cC=0,eC=0,aU=0,cU=0,eU=0); 
    if(closedform)  { fMP = fMP_Chol_closedform; } else { fMP = fMP_Chol; } 
  }
  ## model 4
  if(GxMmle@modelname == 'CholGxM') { 
    fullparam = c('muM','muP','aM','cM','eM','aC','cC','eC','aU','cU','eU','alphaC','alphaU','kappaC','kappaU','epsilonC','epsilonU');
    checkparam = list(muM=0,muP=0,aM=0,cM=0,eM=0,aC=0,cC=0,eC=0,aU=0,cU=0,eU=0,alphaC=0,alphaU=0,kappaC=0,kappaU=0,epsilonC=0,epsilonU=0); 
    if(closedform)  { fMP = fMP_CholGxM_closedform; } else { fMP = fMP_CholGxM; } 
  }
  ## model 5
  if(GxMmle@modelname == 'NLMainGxM') { 
    fullparam = c('muM','muP','aM','cM','eM','aU','cU','eU','beta1','beta2','alphaU','kappaU','epsilonU');
    checkparam = list(muM=0,muP=0,aM=0,cM=0,eM=0,aU=0,cU=0,eU=0,beta1=0,beta2=0,alphaU=0,kappaU=0,epsilonU=0); 
    if(closedform)  { fMP = fMP_NLMainGxM_closedform; } else { fMP = fMP_NLMainGxM; } 
  }
  ## model 6
  if(GxMmle@modelname == 'CorrGxM') { 
    fullparam = c('muM','muP','aM','cM','eM','aP','cP','eP','rA','rC','rE','alphaP','kappaP','epsilonP');
    checkparam = list(muM=0,muP=0,aM=0,cM=0,eM=0,aP=0,cP=0,eP=0,rA=0,rC=0,rE=0,alphaP=0,kappaP=0,epsilonP=0); 
    if(closedform)  { fMP = fMP_CorrGxM_closedform; } else { fMP = fMP_CorrGxM; } 
  }
  ## model 7
  if(GxMmle@modelname == 'CholNonLin') { 
    fullparam = c('muM','muP','aM','cM','eM','aC','cC','eC','aU','cU','eU','gamma1','gamma2','gamma3','delta1','delta2','delta3');
    checkparam = list(muM=0,muP=0,aM=0,cM=0,eM=0,aC=0,cC=0,eC=0,aU=0,cU=0,eU=0,gamma1=0,gamma2=0,gamma3=0,delta1=0,delta2=0,delta3=0); 
    fMP = fMP_CholNonLin;
  }
  ## model 8
  if(GxMmle@modelname == 'CorrNonLin') { 
    fullparam = c('muM','muP','aM','cM','eM','aP','cP','eP','rA','rC','rE','lambda1','lambda2','lambda3');
    checkparam = list(muM=0,muP=0,aM=0,cM=0,eM=0,aP=0,cP=0,eP=0,rA=0,rC=0,rE=0,lambda1=0,lambda2=0,lambda3=0); 
    fMP = fMP_CorrNonLin;
  }

  for (i in 1:length(GxMmle@par)) { 
    ind = which(fullparam %in% (names(GxMmle@par)[i])); 
    checkparam[[ind]] = GxMmle@par[i]; 
  }

  graddirect = max(which(max(abs(GxMmle@gradient)) == abs(GxMmle@gradient)));

  localparam = unlist(checkparam);
  localcomp = GxMmle@par[graddirect];
  localleftbound = localcomp - abs(localcomp)*localpercentage/100;
  localrightbound = localcomp + abs(localcomp)*localpercentage/100;

  if ( (names(localcomp) %in% c('rA','rC','rE')) & (localcomp > 0) )  localrightbound = localcomp;
  if ( (names(localcomp) %in% c('rA','rC','rE')) & (localcomp < 0) )  localleftbound = localcomp;

  localseq = seq(from=localleftbound, to=localrightbound, length.out=localrange);
  localresponse = rep(0, localrange);

  ind = max(which(fullparam %in% (names(GxMmle@par)[graddirect]))); 

  if(closedform) {
    if(coreNumber == 1) {
      for (k in 1:localrange) { localparam[[ind]] = localseq[k]; localresponse[k] = likelihood_closedform(param=localparam,paramnames=fullparam,zeroset=zeroset,listdata=listdata,fMP=fMP); }
    }
    if(coreNumber > 1) {
      cluster = makeCluster(coreNumber);
      for (k in 1:localrange) { localparam[[ind]] = localseq[k]; localresponse[k] = likelihoodcluster_closedform(param=localparam,paramnames=fullparam,zeroset=zeroset,listdata=listdata,fMP=fMP,cl=cluster); }
      stopCluster(cluster);
    }
  } else {
    PointsWK = ghpoints3(K=K);
    if(coreNumber == 1) {
      for (k in 1:localrange) { localparam[[ind]] = localseq[k]; localresponse[k] = likelihood(param=localparam,paramnames=fullparam,zeroset=zeroset,listdata=listdata,fMP=fMP,PointsW=PointsWK); }
    }
    if(coreNumber > 1) {
      cluster = makeCluster(coreNumber);
      for (k in 1:localrange) { localparam[[ind]] = localseq[k]; localresponse[k] = likelihoodcluster(param=localparam,paramnames=fullparam,zeroset=zeroset,listdata=listdata,fMP=fMP,PointsW=PointsWK,cl=cluster); }
      stopCluster(cluster);
    }
  }  # end of  "if(closedform)"

  gapindex = max(which(max(abs(diff(localresponse))) == abs(diff(localresponse))));
  localparam1 = checkparam; localparam1[[ind]] = localseq[gapindex];
  localparam2 = checkparam; localparam2[[ind]] = localseq[gapindex + 1];
  if(closedform) {
    localresponse1 = unlist(lapply(X=listdata, FUN=fMP, param=localparam1)); 
    localresponse2 = unlist(lapply(X=listdata, FUN=fMP, param=localparam2));
  } else {
    localresponse1 = unlist(lapply(X=listdata, FUN=fMP, param=localparam1, PointsW=PointsWK)); 
    localresponse2 = unlist(lapply(X=listdata, FUN=fMP, param=localparam2, PointsW=PointsWK));
  }  # end of  "if(closedform)"
  diffresponse = abs(localresponse1-localresponse2); 
  outlierID = 0;
  if(max(diffresponse) > (outlyingextent)) { 
    outlierID = names(listdata)[max(which(max(diffresponse)==diffresponse))]; 
    cat('The pair of observations with ID number', outlierID, 'may cause singular issue in numerical optimization.\n');
    cat('Fitting the data without this pair may produce a more reasonable result.\n');
  }

  return(list(locallikelihood=localresponse, localfMP=data.frame(localresponse1,localresponse2,diffresponse), outlierID=outlierID));
}
