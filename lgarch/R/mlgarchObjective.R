mlgarchObjective <-
function(pars, aux)
{
  #check parameters:
  if( any(is.na(pars)) || any(pars<=aux$lower) || any(pars>=aux$upper) ){
    chk.conds <- FALSE
  }else{
    chk.conds <- TRUE
  } #end if..else..

  if(chk.conds){
    #recursion:
    uadj <- mlgarchRecursion1(pars, aux)
    if(aux$yanyrowiszeron > 0){
      uadj <- uadj[-aux$yzerowhichrows,]
    }
    #compute objective value:
    mS <- matrix(NA,aux$m,aux$m) #m x m covariance matrix
    mS[lower.tri(mS)] <- pars[aux$cov.indx]
    mS[upper.tri(mS)] <- t(mS)[upper.tri(t(mS))]
    diag(mS) <- pars[aux$sigma2u.indx]
    mSlnDet <- as.numeric(determinant(mS, logarithm=TRUE)$modulus)
    mSinv <- solve(mS)
    tuadjmSuadj <- rowSums( (uadj%*%mSinv)*uadj )
    term1 <- -aux$ynonzerorowsn*aux$m*log(2*pi)/2
    term2 <- -aux$ynonzerorowsn*mSlnDet/2
    term3 <- -sum(tuadjmSuadj)/2
    objective.value <- term1 + term2 + term3
    #check objective value:
    if(is.na(objective.value) || abs(objective.value) == Inf){
      objective.value <- aux$objective.penalty
    }
  }else{
    objective.value <- aux$objective.penalty
  } #end if(chk.conds)
  return(objective.value)
}
