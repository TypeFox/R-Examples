# The EMbC Package for R
# 
# Copyright 2013, 2014, 2015 Joan Garriga <jgarriga@ceab.csic.es>, Aitana Oltra <aoltra@ceab.csic.es>, John R.B. Palmer <johnrbpalmer@gmail.com>, Frederic Bartumeus <fbartu@ceab.csic.es>
#   
#   EMbC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or  (at your option) any later version.
# 
# EMbC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses.


# prior smoothing function
# -------------------------

priorSmth <- function(bCP,smth){
	tSmth <- smth*3600/2
	s <- ncol(bCP@X)-1
	t <- ncol(bCP@X)
	X <- bCP@X
	U <- bCP@U
	T <- c(bCP@spn,0)
	for (j in seq(nrow(X))){
		bCP@X[j,c(s,t)] <- c(0,0)
		bCP@U[j,c(s,t)] <- c(0,0)
		lLoc <- j
		dTime <- 0
		while (lLoc > 1 & dTime <= tSmth){
			lLoc <- lLoc -1
			dTime <- dTime + T[lLoc]
			if (dTime <= tSmth){
				bCP@X[j,c(s,t)] <- bCP@X[j,c(s,t)] + X[lLoc,c(s,t)] * U[lLoc,c(s,t)] * rep((tSmth-dTime)/tSmth,2)
				bCP@U[j,c(s,t)] <- bCP@U[j,c(s,t)] + U[lLoc,c(s,t)] * rep((tSmth-dTime)/tSmth,2)
				}
			}
		rLoc <- j
		dTime <- 0
		while (rLoc < nrow(X) & dTime <=tSmth){
			bCP@X[j,c(s,t)] <- bCP@X[j,c(s,t)] + X[rLoc,c(s,t)] * U[rLoc,c(s,t)] * rep((tSmth-dTime)/tSmth,2)
			bCP@U[j,c(s,t)] <- bCP@U[j,c(s,t)] + U[rLoc,c(s,t)] * rep((tSmth-dTime)/tSmth,2)
			dTime <- dTime + T[rLoc]
			rLoc <- rLoc +1
			}
		if (bCP@U[j,s] != 0) bCP@X[j,s] <- bCP@X[j,s]/bCP@U[j,s]
		if (bCP@U[j,t] != 0) bCP@X[j,t] <- bCP@X[j,t]/bCP@U[j,t]
		bCP@U[j,c(s,t)] <- bCP@U[j,c(s,t)]/(rLoc-lLoc)
		}
	return(bCP)
	}


# posterior (singles) smoothing function
# ------------------------------------

postSmth <- function(pth,dlta){
  singls <- which(apply(array(seq(pth@n-2)),1,function(j){(pth@A[j]!=pth@A[j+1] & pth@A[j+1]!=pth@A[j+2] & pth@A[j]==pth@A[j+2])}))
  singls <- array(singls) +1
  grpLst <- list()
  sngGrp <- singls[1]
  for (s in seq(2,dim(singls))){
    if (singls[s]==singls[s-1]+1){
      sngGrp <- c(sngGrp,singls[s])
    } else {
      grpLst <- append(grpLst,list(sngGrp))
      sngGrp <- singls[s]
    }
  }
  grpLst <- append(grpLst,list(sngGrp))
  for (g in seq(length(grpLst))){
    if (!any(pth@A[grpLst[[g]]]==pth@k+1)){
		if (length(grpLst[[g]])>1){
		  dltGrp <- apply(array(grpLst[[g]]),1,function(loc){pth@W[loc,pth@A[loc]]-pth@W[loc,pth@A[loc-1]]})
		  while (min(dltGrp)<dlta){
			chg <- which.min(dltGrp)
			loc <- grpLst[[g]][chg]
			pth@A[loc] <- pth@A[loc-1]
			dltGrp[max(1,(chg-1)):min(length(dltGrp),(chg+1))] <- 1.0
		  }
		} else {
		  loc <- grpLst[[g]]
		  if ((pth@W[loc,pth@A[loc]]-pth@W[loc,pth@A[loc-1]])<dlta){
			pth@A[loc] <- pth@A[loc-1]
		  }
		}
	}
  }
  return(pth)}
