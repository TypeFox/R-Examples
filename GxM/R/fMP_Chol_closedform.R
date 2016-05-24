fMP_Chol_closedform <-
function(listdata, param) {

  ## model 3 closedform

  m = c(listdata$M1,listdata$M2); 
  p = c(listdata$P1,listdata$P2); 
  rG = listdata$rG;

  aM = param$aM;  cM = param$cM;  eM = param$eM;
  aU = param$aU;  cU = param$cU;  eU = param$eU;

  muM = param$muM;  muP = param$muP;
  aC = param$aC;  cC = param$cC;  eC = param$eC;

  entry11 = aM^2 + cM^2 + eM^2;   
  entry12 = aM^2*rG + cM^2;
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
  muP_M = muP + aC*muAC_M[1:2] + cC*muAC_M[3] + eC*muE_M;

    am = aC - eC*aM / eM;
    cm = cC - eC*cM / eM;
    trans = matrix(c(am, 0, cm,
                     0, am, cm),2,3,byrow=T);
    cov_Cpart = trans%*%covAC_M%*%t(trans);
    entry11 = aU^2 + cU^2 + eU^2; 
    entry12 = aU^2*rG + cU^2;
    cov_Upart = matrix(c(entry11, entry12, 
                         entry12, entry11),2,2)
  covP_M = cov_Cpart + cov_Upart;  # cov(P|M)=E(cov(P|ACM))+cov(E(P|ACM));

  fp_M = log(2*pi) + 0.5*log(det(covP_M)) + 0.5*(p-muP_M) %*% solve(covP_M) %*% (p-muP_M);
  return(fm+fp_M);
}
