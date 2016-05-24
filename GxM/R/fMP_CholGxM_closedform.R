fMP_CholGxM_closedform <-
function(listdata, param) {

  ## model4 closedform

  m = c(listdata$M1,listdata$M2); 
  p = c(listdata$P1,listdata$P2); 
  rG = listdata$rG;

  aM = param$aM;  cM = param$cM;  eM = param$eM;
  aU = param$aU;  cU = param$cU;  eU = param$eU;

  muM = param$muM;  muP = param$muP;
  aC = param$aC;  cC = param$cC;  eC = param$eC;
  alphaC = param$alphaC;  alphaU = param$alphaU;
  kappaC = param$kappaC;  kappaU = param$kappaU;
  epsilonC = param$epsilonC;  epsilonU = param$epsilonU;

  entry11 = aM^2 + cM^2 + eM^2   
  entry12 = aM^2*rG + cM^2 
  SigmaM = matrix(c(entry11, entry12, 
                    entry12, entry11),2,2)
  invSigmaM = solve(SigmaM);
  fm = log(2*pi) + 0.5*log(det(SigmaM)) + 0.5*(m-muM) %*% invSigmaM %*%(m-muM); 

  SigmaAC = matrix(c(1, rG, 0,
                     rG, 1, 0,
                     0,  0, 1),3,3);
  SigmaACM = matrix(c(aM, aM*rG, cM, aM*rG, aM, cM), 3,2)

  muAC_M = SigmaACM %*% invSigmaM %*% (m-muM);
  covAC_M = SigmaAC - SigmaACM %*% invSigmaM %*% t(SigmaACM);
  muE_M = (m - muM - aM*muAC_M[1:2] - cM*muAC_M[3]) / eM; 
  muP_M = muP + (aC+alphaC*m)*muAC_M[1:2] + (cC+kappaC*m)*muAC_M[3] + (eC+epsilonC*m)*muE_M;

    am1 = aC+alphaC*m[1] - (eC+epsilonC*m[1])*aM / eM;
    am2 = aC+alphaC*m[2] - (eC+epsilonC*m[2])*aM / eM;
    cm1 = cC+kappaC*m[1] - (eC+epsilonC*m[1])*cM / eM;
    cm2 = cC+kappaC*m[2] - (eC+epsilonC*m[2])*cM / eM;
    trans = matrix(c(am1, 0, cm1,
                     0, am2, cm2),2,3,byrow=T);
    cov_Cpart = trans%*%covAC_M%*%t(trans);
    entry11 = (aU+alphaU*m[1])^2 + (cU+kappaU*m[1])^2 + (eU+epsilonU*m[1])^2; 
    entry12 = (aU+alphaU*m[1])*(aU+alphaU*m[2])*rG + (cU+kappaU*m[1])*(cU+kappaU*m[2]);
    entry22 = (aU+alphaU*m[2])^2 + (cU+kappaU*m[2])^2 + (eU+epsilonU*m[2])^2;
    cov_Upart = matrix(c(entry11, entry12, 
                         entry12, entry22),2,2)
  covP_M = cov_Cpart + cov_Upart;  # cov(P|M)=E(cov(P|ACM))+cov(E(P|ACM));

  fp_M = log(2*pi) + 0.5*log(det(covP_M)) + 0.5*(p-muP_M) %*% solve(covP_M) %*% (p-muP_M);
  return(fm+fp_M);
}
