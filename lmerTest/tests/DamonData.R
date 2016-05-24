testSat <- TRUE
testKen <- TRUE

if(testSat){
  ## checking on big data sets
  
  require(lmerTest)
  
  load(system.file("testdata", "ttDamon.RData", package="lmerTest"))
  
  testFirstR02.test=lmer(RT2LogR~condition*v_matri_plus_n_indef_freq_sq+
                           (1+condition|Subject),ttDamon, REML=TRUE)
  
  an.sat <- anova(testFirstR02.test)
  
  TOL <- 1e-4 # for the check
  
  #with 4 decimals should agree with SAS output
  #numbers before decimals should agree with SAS output
  stopifnot(
    all.equal(an.sat[,"DenDF"], c(76.206, 6074.6, 6146.1), tol = TOL), 
    all.equal(an.sat[, "F.value"], c(9.2182, 2.8360, 6.6544), tol = TOL)
    , TRUE)
  
  if(testKen){
    if(require(pbkrtest))
      an.kenw <- anova(testFirstR02.test, ddf="kenw")
    stopifnot(ncol(an.kenw) == 6)
  }  
    
}



