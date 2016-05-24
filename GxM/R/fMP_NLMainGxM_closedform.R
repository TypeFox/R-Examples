fMP_NLMainGxM_closedform <-
function(listdata, param) {

  ## model 5 closedform

  m = c(listdata$M1,listdata$M2); 
  p = c(listdata$P1,listdata$P2); 
  rG = listdata$rG;

  aM = param$aM;  cM = param$cM;  eM = param$eM;
  aU = param$aU;  cU = param$cU;  eU = param$eU;

  muM = param$muM;  muP = param$muP;
  beta1 = param$beta1;  beta2 = param$beta2; 
  alphaU = param$alphaU;  kappaU = param$kappaU;  epsilonU = param$epsilonU;

  entry11 = aM^2 + cM^2 + eM^2   
  entry12 = aM^2*rG + cM^2 
  SigmaM = matrix(c(entry11, entry12, 
                    entry12, entry11),2,2)
  invSigmaM = solve(SigmaM);
  fm = log(2*pi) + 0.5*log(det(SigmaM)) + 0.5*(m-muM) %*% invSigmaM %*%(m-muM); 

  muP_M = muP + beta1*m + beta2*(m^2);
  entry11 = (aU+alphaU*m[1])^2 + (cU+kappaU*m[1])^2 + (eU+epsilonU*m[1])^2;
  entry12 = (aU+alphaU*m[1])*(aU+alphaU*m[2])*rG + (cU+kappaU*m[1])*(cU+kappaU*m[2]);
  entry22 = (aU+alphaU*m[2])^2 + (cU+kappaU*m[2])^2 + (eU+epsilonU*m[2])^2;
  covP_M = matrix(c(entry11, entry12, 
                    entry12, entry22),2,2)

  fp_M = log(2*pi) + 0.5*log(det(covP_M)) + 0.5*(p-muP_M) %*% solve(covP_M) %*% (p-muP_M);
  return(fm+fp_M);
}
