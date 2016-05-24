iqQ2Est <-
function (H2, Y, A2, s2ints, ...){
  s2vars = as.data.frame (H2);
  s2names = names (s2vars);
  s2ints = as.vector (s2ints);

  s2. = as.matrix (cbind (1, s2vars, A2, A2*s2vars[,s2ints]));
  colnames (s2.) = c ("intercept", s2names, "A2",
             paste(s2names[s2ints], "A2", sep=":")); 
  p20 = ncol (s2vars);
  p2 = ncol (s2.);
  H20. = s2.[, 1:(p20+1)];
  H21. = s2.[, (p20+2):p2]*A2;
  
  ## second-stage regression
  s2Fit = lm (Y ~ s2. - 1., ...);
  betaHat2 = s2Fit$coefficients;
  betaHat20 = betaHat2[1:(p20+1)];
  betaHat21 = betaHat2[(p20+2):p2];
  
  ## vector of optimal second-stage txts
  optA2 = sign (H21. %*% betaHat21);
  
  ## main effect term and contrast function for model building
  main = H20. %*% betaHat20;
  contrast = H21. %*% betaHat21; 
  
  list ("betaHat20"=betaHat20, "betaHat21"=betaHat21, "s2Fit"=s2Fit,
        "optA2"=optA2, "main"=main, "contrast"=contrast,
        "s2ints"=s2ints, "A2"=A2);     
}
