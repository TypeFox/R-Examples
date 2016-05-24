fMP_CholNonLin <-
function(listdata, param, PointsW) {

  ## model 7

  m = c(listdata$M1,listdata$M2); 
  p1 = listdata$P1;  p2 = listdata$P2; 
  rG = listdata$rG;

  aM = param$aM;  cM = param$cM;  eM = param$eM;
  aU = param$aU;  cU = param$cU;  eU = param$eU;
  muM = param$muM;
  muP = param$muP;  
  aC = param$aC;  cC = param$cC;  eC = param$eC;
  gamma1 = param$gamma1;  gamma2 = param$gamma2;  gamma3 = param$gamma3;
  delta1 = param$delta1;  delta2 = param$delta2;  delta3 = param$delta3;


  Points = PointsW$Points;
  W=PointsW$W;
  K3 = length(W);

  entry1 = aM^2 + cM^2 + eM^2;   
  entry2 = aM^2*rG + cM^2; 
  SigmaM = matrix(c(entry1, entry2, 
                    entry2, entry1),2,2);
  invSigmaM = solve(SigmaM);
  fm = log(2*pi) + 0.5*log(det(SigmaM)) + 0.5*(m-muM)%*%invSigmaM%*%(m-muM); 

  SigmaAC = matrix(c( 1, rG, 0,
                     rG,  1, 0,
                      0,  0, 1), 3,3)
  SigmaACM = matrix(c( aM, aM*rG, cM, aM*rG, aM, cM), 3,2)
  covAC_M = SigmaAC - SigmaACM %*%invSigmaM%*%t(SigmaACM);
  trans = sqrt(2)*t(chol(covAC_M));
  tmpA = SigmaACM %*% invSigmaM; 
  tmpB = tmpA %*% c(muM,muM);
  tmpC = tmpA %*% m - tmpB;  
  tmpD = matrix(rep(tmpC,K3),3,K3)
  AC = tmpD + trans %*% Points;
  CC = t(matrix(rep(AC[3,],2),K3,2))

  tmpA = matrix(rep(m-muM,K3),2,K3)    
  E = (tmpA - aM*AC[1:2,] - cM*CC) / eM;    
  muP_ACM = matrix(muP,2,K3)+aC*AC[1:2,]+cC*CC+eC*E+gamma1*AC[1:2,]^2+gamma2*CC^2+gamma3*E^2; 

  entry11 = (aU+delta1*AC[1,])^2 + (cU+delta2*AC[3,])^2 + (eU+delta3*E[1,])^2; 
  entry12 = (aU+delta1*AC[1,])*(aU+delta1*AC[2,])*rG + (cU+delta2*AC[3,])^2;
  entry22 = (aU+delta1*AC[2,])^2 + (cU+delta2*AC[3,])^2 + (eU+delta3*E[2,])^2;

  det_covP_ACM = entry11*entry22 - entry12^2;
  infind = det_covP_ACM < 1e-6;
  entry11[infind] = 1e-3;
  entry12[infind] = 0;
  entry22[infind] = 1e-3;

  tmpA = (p1-muP_ACM[1,])*entry22 - (p2-muP_ACM[2,])*entry12;  
  tmpB = -(p1-muP_ACM[1,])*entry12 + (p2-muP_ACM[2,])*entry11;  
  tmpC = (tmpA*(p1-muP_ACM[1,]) + tmpB*(p2-muP_ACM[2,])) / det_covP_ACM;
  inversevalue = 2*pi*sqrt(det_covP_ACM)*exp(0.5*tmpC);

  fp_M = -log(sum(W/inversevalue)) + 1.5*log(pi);

  return(fm+fp_M);
}
