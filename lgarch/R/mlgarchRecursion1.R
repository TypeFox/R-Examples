mlgarchRecursion1 <-
function(pars, aux)
{
  #matrices:
  lny2adj <- aux$lny2adj
  uadj <- matrix(0,aux$nmaxpq,aux$m)
  mInnov <- t(matrix(rep(pars[aux$const.indx],aux$n),aux$m,aux$n))
  if(aux$xreg.k > 0){
    xpars <- matrix(pars[aux$xreg.indx],aux$m,aux$xreg.m)
    mInnov <- mInnov + aux$xreg %*% t(xpars)
    if(aux$maxpq > 0){
      innovMeans <- colMeans(mInnov)
      mInnov <- rbind(t(matrix(rep(innovMeans,aux$maxpq),aux$m,aux$maxpq)),
        mInnov)
    } #end if(..maxpq > 0)
  }else{
    if(aux$maxpq > 0){
      innovMeans <- pars[aux$const.indx]
      mInnov <- rbind(t(matrix(rep(innovMeans,aux$maxpq),aux$m,aux$maxpq)),
        mInnov)
    } #end if(..maxpq > 0)
  } #end if(..$xreg.k > 0)
  if(aux$ar > 0){
    phi1 <- matrix(pars[aux$ar.indx],aux$m,aux$m)
  }else{ phi1 <- matrix(0,aux$m,aux$m) }
  if(aux$ma > 0){
    theta1 <- matrix(pars[aux$ma.indx],aux$m,aux$m)
  }else{ theta1 <- matrix(0,aux$m,aux$m) }

  #recursion:
  if(aux$c.code){
    tmp <- VARMARECURSION1(as.numeric(aux$maxpq),
      as.numeric(aux$nmaxpq), as.numeric(aux$m),
      as.numeric(uadj), as.numeric(lny2adj), as.numeric(mInnov),
      as.numeric(phi1), as.numeric(theta1),
      as.numeric(aux$yiszeroadj))
    uadj <- matrix(tmp$mU,aux$nmaxpq,aux$m)
  }else{
    for(i in c(aux$maxpq+1):aux$nmaxpq){
      Elny2adj <- mInnov[i,] + phi1%*%lny2adj[c(i-1),] + theta1%*%uadj[c(i-1),]
      if(aux$yanyrowiszeroadj[i]==1){
        for(j in 1:aux$m){
          if(aux$yiszeroadj[i,j]==1){
            lny2adj[i,j] <- Elny2adj[j]
          } #end if(is zero)
        } #end for(j..)
      }
      uadj[i,] <- lny2adj[i,] - Elny2adj
    } #end for(i..) loop
  } #end if(aux$c.code)

  #output:
  if(aux$maxpq > 0){
    uadj <- uadj[-c(1:aux$maxpq),]
  }
  if(aux$verboseRecursion){
    if(aux$c.code){
      lny2adj <- matrix(tmp$mY,aux$nmaxpq,aux$m)
    } #end if(aux$c.code)
    if(aux$maxpq > 0){
      lny2adj <- lny2adj[-c(1:aux$maxpq),]
    } #end if(aux$maxpq > 0)
    uadj <- cbind(uadj, lny2adj)
    colnames(uadj) <- c(paste("u",1:aux$m,sep=""),
      paste("lny2no",1:aux$m,sep=""))
  } #end if(aux$verboseRecursion)
  return(uadj)
}
