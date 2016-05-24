######################################
#	Internal function called by dtRates(...)
#
#
segMap = function(ephy, tau) {
	tol <- 0.0001;
#	tH <- max(branching.times(as.phylo.bammdata(ephy)));
	tH <- max(ephy$end);
	remainder <- numeric(max(ephy$edge[,1]));
	dtsegs <- vector("list",nrow(ephy$edge));
	for (i in 1:nrow(ephy$edge)) {
		if (remainder[ephy$edge[i,1]] > 0) {
			if (ephy$begin[i]/tH+remainder[ephy$edge[i,1]] > ephy$end[i]/tH) {
				remainder[ephy$edge[i,2]] <- ephy$begin[i]/tH + remainder[ephy$edge[i,1]] - ephy$end[i]/tH;
				segs <- ephy$begin[i]/tH;
			}
			else {
				segs <- seq(ephy$begin[i]/tH+remainder[ephy$edge[i,1]], ephy$end[i]/tH, tau);
				segs <- c(ephy$begin[i]/tH,segs);
			}
		}
		else {
			segs <- seq(ephy$begin[i]/tH, ephy$end[i]/tH, tau);
		}
		if (length(segs) > 1) {
			if (ephy$end[i]/tH - tail(segs,1) > tol) {
				remainder[ephy$edge[i,2]] <- tau - (ephy$end[i]/tH-tail(segs,1));
				segs <- c(segs,ephy$end[i]/tH);
			}
			segs <- rep(segs,each=2);
			segs <- segs[-c(1,length(segs))];
			segs <- matrix(segs,ncol=2,byrow=TRUE);
			segs <- cbind(rep(ephy$edge[i,2],nrow(segs)), segs);
		}
		else {
			if (remainder[ephy$edge[i,1]] == 0) {
				remainder[ephy$edge[i,2]] <- tau - (ephy$end[i]/tH-tail(segs,1));
			}
			segs <- matrix(c(ephy$edge[i,2], ephy$begin[i]/tH, ephy$end[i]/tH),nrow=1,ncol=3); 
		}
		dtsegs[[i]] <- segs;
	}
	dtsegs <- do.call(rbind,dtsegs);
	dtsegs[,2] <- dtsegs[,2]*tH;
	dtsegs[,3] <- dtsegs[,3]*tH;
	return(dtsegs);	
}



# segMap = function(nodes,begin,end,tau)
# {
	# foo = function(x,tau)
	# {
		# len = (x[3] - x[2])/tau; if (len%%1 == 0) len = len+1;
		# ret = seq(x[2],x[3],length.out=len);
		# if(length(ret) == 1) return(matrix(x,nrow=1));
		# #ret = seq(x[2],x[3],length.out=length(ret));
		# ret = rep(ret,each=2); ret=ret[-c(1,length(ret))];
		# ret = matrix(ret,ncol=2,byrow=TRUE);
		# return(cbind(matrix(rep(as.integer(x[1]),nrow(ret)),ncol=1), ret));
	# }
	# times = cbind(nodes,begin,end);
	# ret = apply(times,1,foo,tau);
	# return(do.call(rbind,ret));	
# }
