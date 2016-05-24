ARsim <-
function(margSimFn, acceptFn, N, keepTestInfo=FALSE){
  margSims <- replicate(N, eval(margSimFn()), simplify=FALSE);
  testList <- lapply(margSims, acceptFn);
  keepers <- sapply(testList, function(x){x[[1]]}, simplify=TRUE); 
  if (keepTestInfo){
    testInfo <- sapply(testList, function(x){x[[2]]}, simplify=FALSE); # testList has the T/F and testInfo
    return(list(acceptSims=margSims[keepers], testVals=testInfo[keepers]));
  } #one more list layer
  else{ return(list(acceptSims=margSims[keepers])); } #one less list layer (listr of Xs)
#  sum(keepers);
}

