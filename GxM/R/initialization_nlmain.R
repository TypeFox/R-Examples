initialization_nlmain <-
function(dataset, rGvalue, fullparam, zeroset, manualinitial) {

  beta1out = is.element('beta1', zeroset); 
  beta2out = is.element('beta2', zeroset);

  if(length(manualinitial) > 0) {
    muPout = is.element('muP', zeroset); 
    muPoutinitial = is.element('muP', names(manualinitial));
    if (muPout && muPoutinitial)  { cat('muP can not be manually initialized.\n'); }

    beta1outinitial = is.element('beta1', names(manualinitial));
    if (beta1out && beta1outinitial)  { cat('beta1 can not be manually initialized.\n'); }

    beta2outinitial = is.element('beta2', names(manualinitial));
    if (beta2out && beta2outinitial)  { cat('beta2 can not be manually initialized.\n'); }
  } 

  M = c(dataset$M1,dataset$M2);
  P = c(dataset$P1,dataset$P2);
  if(beta1out) { 
    Mlinear = NULL; 
  } else { Mlinear = M; }
  if(beta2out) { 
    Msquare = NULL; 
  } else { Msquare = M^2; }

  X = cbind(rep(1,length(M)),Mlinear,Msquare);
  Lmfit = lm.fit(x=X,y=P);
  coeff = Lmfit$coef;

  muP = coeff[1];
  paramnames = 'muP'; 
  paramvalues = muP;
  if(beta1out) { 
    beta1 = 0; 
    if(beta2out) { beta2 = 0;
    } else { beta2 = coeff[2]; paramnames = c(paramnames,'beta2'); paramvalues = c(paramvalues, beta2); }
  } else { 
    beta1 = coeff[2]; paramnames = c(paramnames,'beta1'); paramvalues = c(paramvalues, beta1);
    if(beta2out)  { beta2 = 0;
    } else { beta2 = coeff[3]; paramnames = c(paramnames,'beta2'); paramvalues = c(paramvalues, beta2); }
  }  
  names(paramvalues) = paramnames;
  lower = paramvalues - 10;
  upper = paramvalues + 10;

  dataset$P1 = Lmfit$resid[1:(length(P)/2)];
  dataset$P2 = Lmfit$resid[-c(1:(length(P)/2))];

  refineind = which(names(manualinitial) %in% c('muP','beta1','beta2'));
  if (sum(refineind) > 0) {
    manualinitial = manualinitial[-refineind];
  }
  param = initialization_chol(dataset=dataset, rGvalue=rGvalue, fullparam=fullparam, zeroset=union(c('muP','beta1','beta2'),zeroset),manualinitial=manualinitial);

  initial = c(paramvalues, param$initial);
  lower = c(lower, param$lower);
  upper = c(upper, param$upper);
  param = rbind(initial, lower, upper);

  paramnames = colnames(param);
  initial = param[1,];

  return(list(paramnames=paramnames,initial=param[1,],lower = param[2,],upper=param[3,]));
}
