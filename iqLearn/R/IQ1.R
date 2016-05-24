IQ1 <-
function (mainObj, cmObj, sigObj, dens, h1main, h1cm, h1sig){
  A1 = cmObj$A1;
  if (dens=="norm"){
    normal = T;
  }
  else if (dens=="nonpar"){
    normal = F;
  }
  else{
    stop ("dens must be one of {norm, nonpar}");
  }

  if (length (mainObj$s1mainInts) > 0){
    mainPos = sum (c (1, h1main, 1,
      h1main[mainObj$s1mainInts])*c (mainObj$alphaHat0,
                                     mainObj$alphaHat1));
    mainNeg = sum (c (1, h1main, -1,
      -h1main[mainObj$s1mainInts])*c (mainObj$alphaHat0,
                                      mainObj$alphaHat1));
  }
  else{
    mainPos = sum (c (1, h1main, 1)*c (mainObj$alphaHat0,
      mainObj$alphaHat1));
    mainNeg = sum (c (1, h1main, -1)*c (mainObj$alphaHat0,
      mainObj$alphaHat1));
  }    


  if (length (cmObj$s1cmInt) > 0){
    cmPos = sum (c (1, h1cm, 1,
      h1cm[cmObj$s1cmInts])*c(cmObj$betaHat10, cmObj$betaHat11)); 
    cmNeg = sum (c (1, h1cm, -1,
      -h1cm[cmObj$s1cmInts])*c(cmObj$betaHat10, cmObj$betaHat11));
  }
  else{
    cmPos = sum (c (1, h1cm, 1)*c(cmObj$betaHat10, cmObj$betaHat11));  
    cmNeg = sum (c (1, h1cm, -1)*c(cmObj$betaHat10, cmObj$betaHat11)); 
  }
  
  if (sigObj$homo){
    sigPos = sigObj$sigPos;
    sigNeg = sigObj$sigNeg;
  }
  else{
    if (length (sigObj$s1sInts) > 0){
      sigPos = exp (sum (c (1, h1sig, 1,
        h1sig[sigObj$s1sInts])*c(sigObj$gammaHat0, sigObj$gammaHat1))/2);
      sigNeg = exp (sum (c (1, h1sig, -1,
        -h1sig[sigObj$s1sInts])*c(sigObj$gammaHat0,
    sigObj$gammaHat1))/2);
    }
    else{
      sigPos = exp (sum (c (1, h1sig, 1)*c(sigObj$gammaHat0,
        sigObj$gammaHat1))/2); 
      sigNeg = exp (sum (c (1, h1sig, -1)*c(sigObj$gammaHat0,
        sigObj$gammaHat1))/2);      
    }
  }
  
  ## Estimate Q1 for both A1=1 and A1=-1 using either normal density
  ## or empirical estimate
  if (normal){
    q1Pos = mainPos + cmPos*(1-2*pnorm (-cmPos/sigPos)) +
      sqrt(2/pi)*sigPos*exp (-cmPos^2/(2*sigPos^2)); 
    q1Neg = mainNeg + cmNeg*(1-2*pnorm (-cmNeg/sigNeg)) +
      sqrt(2/pi)*sigNeg*exp (-cmNeg^2/(2*sigNeg^2)); 
  }
  else{
    opPos = sum (abs (cmPos + sigPos*sigObj$stdResids)*(A1==1)) /
      sum(A1==1);  
    opNeg = sum (abs (cmNeg + sigNeg*sigObj$stdResids)*(A1==-1)) /
      sum(A1==-1); 
    q1Pos = mainPos + opPos;
    q1Neg = mainNeg + opNeg;
  }
  
  q1opt = sign (q1Pos - q1Neg);
  return (list ("q1Pos"=q1Pos, "q1Neg"=q1Neg, "q1opt"=q1opt));
}
