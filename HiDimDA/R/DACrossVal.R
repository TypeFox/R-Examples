### DACrossVal.R  (2012-07-07)
###    
###
### Copyright 2011 A. Pedro Duarte Silva
###
### This file is part of the `HiDimDA' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

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
  dimnames(EvalRes)[[3]] <- c("Ng","Clerr")
  n <- sum(nk)
  if (prior[1]=="proportions") prior <- nk/n
  if (Strfolds) permut <- vector("list",k)
  out <- vector("list",kfold)
  for (i in 1:CVrep) {
	if (Strfolds) for (grp in 1:k) permut[[grp]] <- sort.int(runif(nk[grp]),index.return=TRUE)$ix
	else permut <- sort.int(runif(n),index.return=TRUE)$ix
	for (j in 1:kfold)  {
  		if (Strfolds) {
			out[[j]] <- which(grouping==codes[1])[permut[[1]]][fold(nk[1],kfold,j)]
			for (grp in 2:k) out[[j]] <- c(out[[j]],which(grouping==codes[grp])[permut[[grp]]][fold(nk[grp],kfold,j)])
		}
		else out[[j]] <- permut[fold(n,kfold,j)]
		if (!all(is.finite(out[[j]]))) 
			stop(paste("DACrossVal is not able to create",kfold,"folds, for group sizes equal to:\n",
				paste(nk,collapse=" ")))
	}
	for (j in 1:kfold)  {
		tres <- TrainAlg(data[-out[[j]],],grouping[-out[[j]]],k=k,grpcodes=codes,prior=prior,...)
		EvalRes[(i-1)*kfold+j,,] <- EvalAlg(tres,data[out[[j]],],grouping[out[[j]]],k=k,grpcodes=codes)	
	}
  }
  EvalRes  # return(EvalRes)
}

EvalClrule <- function(darule,VData,Vgrp,k,grpcodes)
{
	errates <- array(dim=k)
	Ng <- array(dim=k)
	if (is.null(dim(VData))) dim(VData) <- c(length(Vgrp),1)
        clres <- predict(darule,VData,grpcodes=grpcodes)$class
	for (grpInd in 1:k)  {
		Ng[grpInd] <- length(Vgrp[Vgrp==grpcodes[grpInd]])
		thisgrpclres <- clres[Vgrp==grpcodes[grpInd]]
		levels(thisgrpclres) <- grpcodes
		thisgrperr <- thisgrpclres[grpcodes[grpInd]!=thisgrpclres]
		if (Ng[grpInd]>0) errates[grpInd] <- length(thisgrperr)/Ng[grpInd]
		else errates[grpInd] <- 0
	}
      	list(err=errates,Ng=Ng)  #  return(list(err=errates,Ng=Ng))
      	cbind(Ng,errates)  #  return(cbind(Ng,errates))
}

