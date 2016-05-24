iqQ1MainEst <-
function (object, H1Main, A1, s1mainInts, ...){
  mainResp = object$main;
  s1mainVars = as.data.frame (H1Main);
  s1names = names (s1mainVars);
  s1mainInts = as.vector (s1mainInts);

  if (length (s1mainInts) > 0){  
    s1m. = as.matrix (cbind (1, s1mainVars, A1,
      A1*s1mainVars[,s1mainInts])); 
    colnames (s1m.) = c ("intercept", s1names, "A1",
               paste(s1names[s1mainInts], "A1", sep=":")); 
    p10 = ncol (s1mainVars);
    p11 = ncol (s1mainVars[,s1mainInts]);
    p1 = ncol (s1m.);
    H10. = s1m.[, 1:(p10+1)];
    H11. = s1m.[, (p10+2):p1]*A1;
    
    ## first-stage regression for the main effect term
    s1MainFit = lm (mainResp ~ s1m. - 1, ...);
    alphaHat = s1MainFit$coefficients;
    alphaHat0 = alphaHat[1:(p10+1)];
    alphaHat1 = alphaHat[(p10+2):p1];

    mainPos = as.matrix (cbind (1, s1mainVars, 1,
      s1mainVars[,s1mainInts])) %*% alphaHat;
    mainNeg = as.matrix (cbind (1, s1mainVars, -1,
      -s1mainVars[,s1mainInts])) %*% alphaHat;
  }
  else{
    s1m. = as.matrix (cbind (1, s1mainVars, A1)); 
    colnames (s1m.) = c ("intercept", s1names, "A1"); 
    p10 = ncol (s1mainVars);
    p1 = ncol (s1m.);
    H10. = s1m.[, 1:(p10+1)];
    H11. = 1;
    
    ## first-stage regression for the main effect term
    s1MainFit = lm (mainResp ~ s1m. - 1, ...);
    alphaHat = s1MainFit$coefficients;
    alphaHat0 = alphaHat[1:(p10+1)];
    alphaHat1 = alphaHat[p1];

    mainPos = as.matrix (cbind (1, s1mainVars, 1)) %*% alphaHat;
    mainNeg = as.matrix (cbind (1, s1mainVars, -1)) %*% alphaHat;
  }
        
  list ("alphaHat0"=alphaHat0, "alphaHat1"=alphaHat1,
        "s1MainFit"=s1MainFit, "mainPos"=mainPos, "mainNeg"=mainNeg,
        "s1mainInts"=s1mainInts, "A1"=A1);        
}
