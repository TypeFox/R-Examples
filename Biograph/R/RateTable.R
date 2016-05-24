RateTable <-
function (Bdata,occup,trans)
{ if (!exists(as.character(substitute(occup)))) 
	  print ("Message from RateTable: 'occup' does not exists")
  if (!exists(as.character(substitute(trans))))
      print ("Message from RateTable: 'trans' does not exists")
  namstates <- attr (Bdata,"param")$namstates
  numstates <- length (namstates)
  iagelow <- attr(Bdata,"param")$iagelow
  iagehigh <- attr(Bdata,"param")$iagehigh
  nage <- attr(Bdata,"param")$nage
  tstate <-occup$state_occup
  tsjt <- occup$tsjt
  meanage <- trans$meanage
  trans <- trans$trans


  agelist2 <- c(iagelow:iagehigh)
  ist2 <-agelist2
  # trans <-trans$trans
  Stable <- array(0,c(length(agelist2)+2,numstates,numstates+4))
  dimnames (Stable) <- list(Age=c(iagelow:iagehigh,"Total","MeanAge"),State=namstates,Case=c("Occup","PY","Leaving",namstates,"Censored"))
  for (itrans in 1:numstates) {
    Stable[1:length(agelist2),itrans,] <- cbind(tstate[,itrans],round(tsjt[,itrans],2),apply(trans[,,itrans],1,sum),trans[,,itrans])
    }
    zz2 <- length(agelist2) + 1
  Stable[zz2,,] <- apply(Stable,c(2,3),sum)
  zz <- length(agelist2)
  zz3 <- zz + 2
  for (kk in 1:numstates) Stable[zz3,kk,4:(4+numstates)]<- meanage[kk,]

  # save state occupancy at censoring for microsimulation
 censored_by_age <- Stable[1:nage,,numstates+4]
 censored_by_age <- cbind(rep(0,nrow(censored_by_age)),censored_by_age)
 censored_by_age[,1] <- apply(censored_by_age,1,sum)  # first column: total censored at age ix
 colnames(censored_by_age) <- c("Total",colnames(censored_by_age)[2:(numstates+1)])

 return (list( Stable = aperm(Stable,c(1,3,2)),
               censored_by_age = censored_by_age))
}
