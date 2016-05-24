qLearnS1Est <-
function (object, H1q, A1, s1ints, ...){
  
  s1vars = as.data.frame (H1q);
  s1names = names (s1vars);
  s1ints = as.vector (s1ints);
  Ytilde = object$Ytilde;

  if (length (s1ints) > 0){
    s1. = as.matrix (cbind (1, s1vars, A1, A1*s1vars[,s1ints]));
    colnames (s1.) = c ("intercept", s1names, "A1",
               paste(s1names[s1ints], "A1", sep=":")); 
    p10 = ncol (s1vars);
    p11 = ncol (s1vars[,s1ints]);
    p1 = ncol (s1.);
    H10 = s1.[, 1:(p10+1)];
    H11 = s1.[, (p10+2):p1]*A1;
    
    ## second-stage regression
    s1Fit = lm (Ytilde ~ s1. - 1, ...);
    betaHat1 = s1Fit$coefficients;
    betaHat10 = betaHat1[1:(p10+1)];
    betaHat11 = betaHat1[(p10+2):p1];
    
    ## vector of optimal second-stage txts
    optA1 = sign (H11 %*% betaHat11);
  }
  else{
    s1. = as.matrix (cbind (1, s1vars, A1));
    colnames (s1.) = c ("intercept", s1names, "A1"); 
    p10 = ncol (s1vars);
    p1 = ncol (s1.);
    H10 = s1.[, 1:(p10+1)];
    H11 = 1;
    
    ## second-stage regression
    s1Fit = lm (Ytilde ~ s1. - 1, ...);
    betaHat1 = s1Fit$coefficients;
    betaHat10 = betaHat1[1:(p10+1)];
    betaHat11 = betaHat1[p1];
    
    ## vector of optimal second-stage txts
    optA1 = sign (betaHat11);
  }
  
  list ("betaHat10"=betaHat10, "betaHat11"=betaHat11, "optA1"=optA1,
        "s1Fit"=s1Fit, "s1ints"=s1ints);  
}
