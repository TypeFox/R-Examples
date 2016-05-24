cov.sel.np<-function(T, Y, X, alg, scope, thru, thro, thrc, dat, data.0, data.1,covar, ...){


if(thru<0 || thru>1 || thro<0 || thro>1 || thrc<0){stop("Threshold value(s) out of range")}


if(any(sapply(X,is.ordered))){
		fac<-names(which(sapply(X,is.factor)))[-which(sapply(X,is.ordered))]
	}else{
		fac<-names(which(sapply(X,is.factor)))
		}	
## computing thresholds for factor variables depending on number of levels
maxbw<-1-1/as.numeric(lapply(mapply(levels,dat),length)[fac])
thr<-thru*maxbw 
covars<-covar
thrs<-thr

if(alg==1){
	
datlist<-list(dat,data.0,data.1)
reslist<-list("T ~","Y~","Y~")	
covlist<-vector("list", 3)
bwlist<-vector("list", 3)
## Algoritm A
## Step 1
## creating npregbw object
for(i in 1:3){
	f1 <- as.formula(paste(reslist[[i]], paste(paste(covars), collapse= "+")))
	bw<-npregbw(f1,data=datlist[[i]], ...)
	Xfac<-names(which(bw$iuno))[which((thrs-bw$bw[which(bw$iuno)])>0)]
	Xord<-names(which(bw$iord))[which(bw$bw[which(bw$iord)]<thro)]
	Xcon<-names(which(bw$icon))[which(bw$bw[which(bw$icon)]<thrc)]
	Xs<-sort(unique(c(Xfac,Xord,Xcon,scope)))
	Xfac<-sort(unique(c(Xfac,fac[match(scope,fac)[is.na(match(scope,fac))==FALSE]])))
	bws<-bw$bw[match(Xs,bw$xnames)]
	covlist[[i]]<-Xs
	bwlist[[i]]<-bws
	if(i==1){if(length(Xs)==0){ covlist[[3]]<-covlist[[2]]<-covlist[[1]]
						bwlist[[3]]<-bwlist[[2]]<-bwlist[[1]]		
								break}
		
				covars<-Xs
				maxbwT<-1-1/as.numeric(lapply(mapply(levels,dat),length)[Xfac])
thrs<-thru*maxbwT			}
}

## Return Values
l <- list(X.T = covlist[[1]], Q.0 = covlist[[2]], Q.1 = covlist[[3]],
  bandwidthsQ.0 = bwlist[[2]], bandwidthsQ.1 = bwlist[[3]], regtype=bw$pregtype,bwtype=bw$type, covar = covar)
class(l) <- "cov.sel"
invisible(return(l))

}else if(alg==2){
	
datlist<-list(data.0,dat,data.1,dat)
reslist<-list("Y~","T~","Y~","T~")	
covlist<-vector("list", 4)
bwlist<-vector("list", 4)
## Algoritm B
## Step 1
for(i in 1:4){
	f1 <- as.formula(paste(reslist[[i]], paste(paste(covars), collapse= "+")))
	## creating npregbw object for control
	bw<-npregbw(f1,data=datlist[[i]], ...)
	Xfac<-names(which(bw$iuno))[which((thrs-bw$bw[which(bw$iuno)])>0)]
	Xord<-names(which(bw$iord))[which(bw$bw[which(bw$iord)]<thro)]
	Xcon<-names(which(bw$icon))[which(bw$bw[which(bw$icon)]<thrc)]
	Xs<-sort(unique(c(Xfac,Xord,Xcon,scope)))
	Xfac<-sort(unique(c(Xfac,fac[match(scope,fac)[is.na(match(scope,fac))==FALSE]])))
	bws<-bw$bw[match(Xs,bw$xnames)]
  covlist[[i]]<-Xs
	bwlist[[i]]<-bws
	if(i==2 || i==4){
		covars<-covar
		thrs<-thr
	}
	if(i==1){if(length(Xs)==0){ covlist[[2]]<-covlist[[1]]
						bwlist[[2]]<-bwlist[[1]]		
								i<-2}else{covars<-Xs
		maxbw<-1-1/as.numeric(lapply(mapply(levels,dat),length)[Xfac])
		thrs<-thru*maxbw}
	}
	if(i==3){if(length(Xs)==0){ covlist[[4]]<-covlist[[3]]
						bwlist[[4]]<-bwlist[[3]]		
								break}else{		
		covars<-Xs
		maxbw<-1-1/as.numeric(lapply(mapply(levels,dat),length)[Xfac])
		thrs<-thru*maxbw}
	}	
	
}	
		
## Return Values
l <- list(X.0 = covlist[[1]], X.1 = covlist[[3]], Z.0 = covlist[[2]], Z.1 = covlist[[4]],
 bandwidthsZ.0 = bwlist[[2]], bandwidthsZ.1 = bwlist[[4]], regtype=bw$pregtype,bwtype=bw$type, covar = covar)
class(l) <- "cov.sel"
invisible(return(l))
		
}else if(alg==3){
			
datlist1<-list(dat,data.0,data.1)
reslist1<-list("T ~","Y~","Y~")	
covlist1<-vector("list", 3)
bwlist1<-vector("list", 3)
## Algoritm A
## Step 1
## creating npregbw object
for(i in 1:3){
	f1 <- as.formula(paste(reslist1[[i]], paste(paste(covars), collapse= "+")))
	bw<-npregbw(f1,data=datlist1[[i]], ...)
	Xfac<-names(which(bw$iuno))[which((thrs-bw$bw[which(bw$iuno)])>0)]
	Xord<-names(which(bw$iord))[which(bw$bw[which(bw$iord)]<thro)]
	Xcon<-names(which(bw$icon))[which(bw$bw[which(bw$icon)]<thrc)]
	Xs<-sort(unique(c(Xfac,Xord,Xcon,scope)))
	Xfac<-sort(unique(c(Xfac,fac[match(scope,fac)[is.na(match(scope,fac))==FALSE]])))
	bws<-bw$bw[match(Xs,bw$xnames)]
  covlist1[[i]]<-Xs
	bwlist1[[i]]<-bws
	if(i==1){if(length(Xs)==0){ covlist1[[3]]<-covlist1[[2]]<-covlist1[[1]]
						bwlist1[[3]]<-bwlist1[[2]]<-bwlist1[[1]]		
								break}
			covars<-Xs
			maxbwT<-1-1/as.numeric(lapply(mapply(levels,dat),length)[Xfac])
			thrs<-thru*maxbwT		
				}

}
		
covars<-covar
thrs<-thr
datlist2<-list(data.0,dat,data.1,dat)
reslist2<-list("Y~","T~","Y~","T~")	
covlist2<-vector("list", 4)
bwlist2<-vector("list", 4)
## Algoritm B
## Step 1
for(i in 1:4){
	f1 <- as.formula(paste(reslist2[[i]], paste(paste(covars), collapse= "+")))
	## creating npregbw object for control
	bw<-npregbw(f1,data=datlist2[[i]], ...)
	Xfac<-names(which(bw$iuno))[which((thrs-bw$bw[which(bw$iuno)])>0)]
	Xord<-names(which(bw$iord))[which(bw$bw[which(bw$iord)]<thro)]
	Xcon<-names(which(bw$icon))[which(bw$bw[which(bw$icon)]<thrc)]
	Xs<-sort(unique(c(Xfac,Xord,Xcon,scope)))
	Xfac<-sort(unique(c(Xfac,fac[match(scope,fac)[is.na(match(scope,fac))==FALSE]])))
	bws<-bw$bw[match(Xs,bw$xnames)]
  covlist2[[i]]<-Xs
	bwlist2[[i]]<-bws
	if(i==2 || i==4){
		covars<-covar
		thrs<-thr
	}
	if(i==1){if(length(Xs)==0){ covlist2[[2]]<-covlist2[[1]]
						bwlist2[[2]]<-bwlist2[[1]]		
								i<-2}else{covars<-Xs
		maxbw<-1-1/as.numeric(lapply(mapply(levels,dat),length)[Xfac])
		thrs<-thru*maxbw}
	}
	if(i==3){if(length(Xs)==0){ covlist2[[4]]<-covlist2[[3]]
						bwlist2[[4]]<-bwlist2[[3]]		
								break}else{		
		covars<-Xs
		maxbw<-1-1/as.numeric(lapply(mapply(levels,dat),length)[Xfac])
		thrs<-thru*maxbw}
	}	

}		
			
## Return Values
l <- list(X.T = covlist1[[1]], Q.0 = covlist1[[2]], Q.1 = covlist1[[3]],
  X.0 = covlist2[[1]], X.1 = covlist2[[3]], Z.0 = covlist2[[2]], Z.1 = covlist2[[4]],
  bandwidthsQ.0 = bwlist1[[2]], bandwidthsQ.1 = bwlist1[[3]],
  bandwidthsZ.0 = bwlist2[[2]], bandwidthsZ.1 = bwlist2[[4]], regtype=bw$pregtype,bwtype=bw$type, covar = covar)
class(l) <- "cov.sel"
invisible(return(l))
	
			
}else{
stop("Wrong Selected Algorithm")
}
}