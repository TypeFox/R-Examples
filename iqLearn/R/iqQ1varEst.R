iqQ1varEst <-
function (object, H1CVar, s1sInts, method, ...){
  cmeanResids = object$cmeanResids;
  A1 = object$A1;
  
  if (method=="homo"){
    homo = T;
  }
  else if (method=="hetero"){
    homo = F;
  }
  else{
    stop ("method must be one of {homo, hetero}");
  }

  if (homo){
    stdDev = sd (as.vector (cmeanResids));
    stdResids = as.vector (cmeanResids)/stdDev;
    gammaHat0 = NULL;
    gammaHat1 = NULL;
    s1VarFit = NULL;
    sigPos = stdDev;
    sigNeg = stdDev;
  }
  else{
    lRes2 = log (cmeanResids^2);
    s1sVars = as.data.frame (H1CVar);
    s1names = names (s1sVars);
    s1sInts = as.vector (s1sInts);

    if (length (s1sInts) > 0){  
      s1v. = as.matrix (cbind (1, s1sVars, A1,
        A1*s1sVars[,s1sInts])); 
      colnames (s1v.) = c ("intercept", s1names, "A1",
                 paste(s1names[s1sInts], "A1", sep=":")); 
      p10 = ncol (s1sVars);
      p11 = ncol (s1sVars[,s1sInts]);
      p1 = ncol (s1v.);
      H10. = s1v.[, 1:(p10+1)];
      H11. = s1v.[, (p10+2):p1]*A1;
      
      ## first-stage log linear regression for var function
      s1VarFit = lm (lRes2 ~ s1v. - 1, ...);
      gammaHat = s1VarFit$coefficients;
      sig = exp (s1v. %*% gammaHat/2);
      stdResids = as.vector (cmeanResids/sig);
      gammaHat[1] = gammaHat[1] + 2*log (sd (stdResids));
      gammaHat0 = gammaHat[1:(p10+1)];
      gammaHat1 = gammaHat[(p10+2):p1];
      stdDev = NULL;
      sig = exp (s1v. %*% gammaHat/2);    
      stdResids = as.vector (cmeanResids)/sig;
      sigPos = exp (as.matrix (cbind (1, s1sVars, 1,
        s1sVars[,s1sInts]))%*%gammaHat/2); 
      sigNeg = exp (as.matrix (cbind (1, s1sVars, -1,
        -s1sVars[,s1sInts])) %*% gammaHat/2);
    }
    else{
      s1v. = as.matrix (cbind (1, s1sVars, A1)); 
      colnames (s1v.) = c ("intercept", s1names, "A1"); 
      p10 = ncol (s1sVars);
      p1 = ncol (s1v.);
      H10. = s1v.[, 1:(p10+1)];
      H11. = 1;
      
      ## first-stage log linear regression for var function
      s1VarFit = lm (lRes2 ~ s1v. - 1, ...);
      gammaHat = s1VarFit$coefficients;
      sig = exp (s1v. %*% gammaHat/2);
      stdResids = as.vector (cmeanResids/sig);
      gammaHat[1] = gammaHat[1] + 2*log (sd (stdResids));
      gammaHat0 = gammaHat[1:(p10+1)];
      gammaHat1 = gammaHat[p1];
      stdDev = NULL;
      sig = exp (s1v. %*% gammaHat/2);    
      stdResids = as.vector (cmeanResids)/sig;
      sigPos = exp (as.matrix (cbind (1, s1sVars, 1))%*%gammaHat/2); 
      sigNeg = exp (as.matrix (cbind (1, s1sVars, -1)) %*% gammaHat/2);
    }      
  }

  list ("stdDev"=stdDev, "stdResids"=stdResids, "gammaHat0"=gammaHat0,
  "gammaHat1"=gammaHat1, "s1VarFit"=s1VarFit, "homo"=homo,
  "sigPos"=sigPos, "sigNeg"=sigNeg, "s1sInts"=s1sInts);
}
