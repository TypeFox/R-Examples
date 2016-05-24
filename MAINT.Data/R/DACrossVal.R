DACrossVal <- function(data,grouping,TrainAlg,EvalAlg=EvalClrule,Strfolds=TRUE,kfold=10,CVrep=20,prior="proportions",...)
{
  fold <- function(n,kfold,fi) return((round((fi-1)*n/kfold)+1):round(fi*n/kfold))

  codes <- levels(grouping)
  nk <- table(grouping)
  k <- nrow(nk)
  nk <- as.numeric(nk)
  trep <- kfold*CVrep
  EvalRes <- array(dim=c(trep,k,2))
  dimnames(EvalRes)[[2]] <- codes
  dimnames(EvalRes)[[3]] <- c("Nk","Clerr")
  n <- sum(nk)
  if (prior[1]=="proportions") prior <- nk/n
  if (Strfolds) permut <- vector("list",k)
  for (i in 1:CVrep) {
	if (Strfolds) for (grp in 1:k) permut[[grp]] <- sort.int(runif(nk[grp]),index.return=TRUE)$ix
	else permut <- sort.int(runif(n),index.return=TRUE)$ix
	for (j in 1:kfold)  {
  		if (Strfolds) {
			out <- which(grouping==codes[1])[permut[[1]]][fold(nk[1],kfold,j)]
			for (grp in 2:k) out <- c(out,which(grouping==codes[grp])[permut[[grp]]][fold(nk[grp],kfold,j)])
		}
		else out <- permut[fold(n,kfold,j)]
		tres <- TrainAlg(data[-out,],grouping[-out],prior=prior,...)
		EvalRes[(i-1)*kfold+j,,] <- EvalAlg(tres,data[out,],grouping[out],k=k)
	}
  }
  EvalRes  # return(EvalRes)
}

EvalClrule <- function(darule,VData,Vgrp,k)
{
	grpcodes <- levels(Vgrp)
	errates <- array(dim=k)
	Nk <- array(dim=k)
        clres <- predict(darule,VData)$class
	for (grpInd in 1:k)  {
		Nk[grpInd] <- length(Vgrp[Vgrp==grpcodes[grpInd]])
		thisgrpclres <- clres[Vgrp==grpcodes[grpInd]]
		levels(thisgrpclres) <- grpcodes
		thisgrperr <- thisgrpclres[grpcodes[grpInd]!=thisgrpclres]
		if (Nk[grpInd]>0) errates[grpInd] <- length(thisgrperr)/Nk[grpInd]
		else errates[grpInd] <- 0
	}
      	list(err=errates,Nk=Nk)  #  return(list(err=errates,Nk=Nk))
      	cbind(Nk,errates)  #  return(cbind(Nk,errates))
}

