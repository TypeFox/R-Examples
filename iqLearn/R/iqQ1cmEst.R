iqQ1cmEst <-
function (object, H1CMean, A1, s1cmInts, ...){
  cmResp = object$contrast;
  s1cmVars = as.data.frame (H1CMean);
  s1names = names (s1cmVars);
  s1cmInts = as.vector (s1cmInts);

  if (length (s1cmInts) > 0){
    s1cm. = as.matrix (cbind (1, s1cmVars, A1,
      A1*s1cmVars[,s1cmInts])); 
    colnames (s1cm.) = c ("intercept", s1names, "A1",
               paste(s1names[s1cmInts], "A1", sep=":")); 
    p10 = ncol (s1cmVars);
    p11 = ncol (s1cmVars[,s1cmInts]);
    p1 = ncol (s1cm.);
    H10. = s1cm.[, 1:(p10+1)];
    H11. = s1cm.[, (p10+2):p1]*A1;
    
    ## first-stage regression for contrast function mean
    s1cmFit = lm (cmResp ~ s1cm. - 1, ...);
    betaHat1 = s1cmFit$coefficients;
    betaHat10 = betaHat1[1:(p10+1)];
    betaHat11 = betaHat1[(p10+2):p1];
    
    cmeanResids = cmResp - s1cm.%*%betaHat1;
    cmPos = as.matrix (cbind (1, s1cmVars, 1,
      s1cmVars[,s1cmInts]))%*%betaHat1; 
    cmNeg = as.matrix (cbind (1, s1cmVars, -1,
      -s1cmVars[,s1cmInts]))%*%betaHat1; 
  }
  else{
    s1cm. = as.matrix (cbind (1, s1cmVars, A1)); 
    colnames (s1cm.) = c ("intercept", s1names, "A1"); 
    p10 = ncol (s1cmVars);
    p1 = ncol (s1cm.);
    H10. = s1cm.[, 1:(p10+1)];
    H11. = 1;
    
    ## first-stage regression for contrast function mean
    s1cmFit = lm (cmResp ~ s1cm. - 1, ...);
    betaHat1 = s1cmFit$coefficients;
    betaHat10 = betaHat1[1:(p10+1)];
    betaHat11 = betaHat1[p1];
    
    cmeanResids = cmResp - s1cm.%*%betaHat1;
    cmPos = as.matrix (cbind (1, s1cmVars, 1))%*%betaHat1; 
    cmNeg = as.matrix (cbind (1, s1cmVars, -1))%*%betaHat1; 
  }

  
  list ("betaHat10"=betaHat10, "betaHat11"=betaHat11,
        "s1cmFit"=s1cmFit, "cmeanResids"=cmeanResids, "cmPos"=cmPos,
        "cmNeg"=cmNeg, "s1cmInts"=s1cmInts, "A1"=A1);     
}
