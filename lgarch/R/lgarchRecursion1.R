lgarchRecursion1 <-
function(pars, aux)
{
  #make vectors:
  if(aux$xreg.k == 0){
    innov <- pars[1]
    innovMean <- mean(innov)
    innov <- rep(innov, max(aux$n+1,aux$nmaxpq))
  }else{
    innov <- pars[1] + aux$xreg %*% pars[aux$xreg.indx]
    innovMean <- mean(innov)
    innov <- c(rep(innovMean, max(1,aux$maxpq)), innov)
  }
  #Version 0.4: Elny2 <- innovMean/(1-sum(pars[aux$ar.indx]))
  if(aux$mean.correction){
    lny2adj <- c(rep(aux$Elny2,1:max(1,aux$maxpq)), aux$lny2mc)
  }else{
    lny2adj <- c(rep(aux$Elny2,1:max(1,aux$maxpq)), aux$lny2)
  }
  uadj <- aux$zerosaux

  #recursion:
  if(aux$maxpq > 0){ phi1 <- pars[aux$ar.indx] }else{ phi1 <- 0 }
  if(aux$ma > 0){ theta1 <- pars[aux$ma.indx] }else{ theta1 <- 0 }
  if(aux$c.code){
    iStart <- max(2,aux$maxpq+1) - 1
    iEnd <- max(aux$nmaxpq, aux$n+1)
    tmp <- ARMARECURSION1(as.integer(iStart), as.integer(iEnd),
      as.numeric(phi1), as.numeric(theta1), as.numeric(aux$yzeroadj),
      as.numeric(innov), as.numeric(lny2adj), as.numeric(uadj))
    uadj <- tmp$uadj
  }else{
    for(i in max(2,aux$maxpq+1):max(aux$nmaxpq, aux$n+1)){
      imin1 <- i - 1
      if(aux$yzeroadj[i]==0){
        lny2adj[i] <- innov[i] + phi1*lny2adj[imin1] + theta1*uadj[imin1]
        uadj[i] <- 0
      }else{
        uadj[i] <- lny2adj[i] - innov[i] - phi1*lny2adj[imin1] - theta1*uadj[imin1]
      } #end if(not-0)
    } #end for loop
  } #end if(c.code)

  #output:
  if(aux$verboseRecursion){
    if(aux$c.code){
      result <- cbind(tmp$uadj, tmp$lny2adj)
    }else{
      result <- cbind(uadj, lny2adj)
    } #end if(c.code)
    result <- result[-c(1:max(1,aux$maxpq)),]
  }else{
    result <- uadj[-c(1:max(1,aux$maxpq))]
  } #end if(verboseRecursion)
  return(result)
}
