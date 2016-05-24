qLearnS2Est <-
function (H2, Y, A2, s2ints, ...){
  
  s2vars = as.data.frame (H2);
  s2names = names (s2vars);
  s2ints = as.vector (s2ints);

  if (length (s2ints) > 0){
    s2. = as.matrix (cbind (1, s2vars, A2, A2*s2vars[,s2ints]));
    colnames (s2.) = c ("intercept", s2names, "A2",
               paste(s2names[s2ints], "A2", sep=":")); 
    p20 = ncol (s2vars);
    p21 = ncol (s2vars[,s2ints]);
    p2 = ncol (s2.);
    H20 = s2.[, 1:(p20+1)];
    H21 = s2.[, (p20+2):p2]*A2;
    
    ## second-stage regression
    s2Fit = lm (Y ~ s2. - 1, ...);
    betaHat2 = s2Fit$coefficients;
    betaHat20 = betaHat2[1:(p20+1)];
    betaHat21 = betaHat2[(p20+2):p2];

    ## vector of optimal second-stage txts
    optA2 = sign (H21 %*% betaHat21);
    
    ## maximization
    Ytilde = H20 %*% betaHat20 + abs (H21 %*% betaHat21);
  }
  else{
    s2. = as.matrix (cbind (1, s2vars, A2));
    colnames (s2.) = c ("intercept", s2names, "A2"); 
    p20 = ncol (s2vars);
    p2 = ncol (s2.);
    H20 = s2.[, 1:(p20+1)];
    H21 = 1;
    
    ## second-stage regression
    s2Fit = lm (Y ~ s2. - 1, ...);
    betaHat2 = s2Fit$coefficients;
    betaHat20 = betaHat2[1:(p20+1)];
    betaHat21 = betaHat2[p2];

    ## vector of optimal second-stage txts
    optA2 = sign (betaHat21);
    
    ## maximization
    Ytilde = H20 %*% betaHat20 + abs (betaHat21);
  }
    
  list ("betaHat20"=betaHat20, "betaHat21"=betaHat21,
  "Ytilde"=Ytilde, "optA2"=optA2, "s2Fit"=s2Fit, "s2ints"=s2ints); 
}
