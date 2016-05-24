mkdtsegsPhylo <- function(x, tau, edge) {
	tol <- 0.0001;
	xremainder <- numeric(max(edge[,1]));
	dtsegs <- vector("list",nrow(edge));
	for (i in 1:nrow(x)) {
		if (xremainder[edge[i,1]] > 0) {
			if (x[i,1]+xremainder[edge[i,1]] > x[i,3]) {
				xremainder[edge[i,2]] = x[i,1]+xremainder[edge[i,1]]-x[i,3];
				xx <- x[i,1];
				yy <- x[i,2];
			}
			else {
				xx <- seq(x[i,1]+xremainder[edge[i,1]],x[i,3],tau);
				xx <- c(x[i,1],xx);
				yy <- rep(x[i,2],length(xx));
			}
		}
		else {
			xx <- seq(x[i,1],x[i,3],tau);
			yy <- rep(x[i,2],length(xx));
		}
		if (length(xx) > 1) {
			if (x[i,3] - tail(xx,1) > tol) {
				xremainder[edge[i,2]] <- tau - (x[i,3]-tail(xx,1));
				xx <- c(xx,x[i,3]);
				yy <- c(yy,x[i,2]);
			}
			xx <- rep(xx,each=2);
			xx <- xx[-c(1,length(xx))];
			xx <- matrix(xx,ncol=2,byrow=TRUE);
			yy <- rep(yy,each=2);
			yy <- yy[-c(1,length(yy))];
			yy <- matrix(yy,ncol=2,byrow=TRUE);
			segs <- cbind(xx,yy,rep(edge[i,2],nrow(xx)));	
		}
		else {
			if (xremainder[edge[i,1]] == 0) {
				xremainder[edge[i,2]] <- tau - (x[i,3]-tail(xx,1));
			}
			segs <- matrix(c(x[i,1],x[i,3],x[i,2],x[i,4],edge[i,2]),nrow=1,ncol=5);
		}
		dtsegs[[i]] <- segs;
	}
	return(do.call(rbind,dtsegs));
}

# OLD VERSION
# mkdtsegs = function(x,tau,phy,tH)
# {
	# #bn = sqrt((x[3]-x[1])^2 + (x[4]-x[2])^2);
	# #len = bn/tau; if (len %% 1 == 0) len = len + 1;
	# len = (phy$end[match(x[5],phy$edge[,2])]/tH-phy$begin[match(x[5],phy$edge[,2])]/tH)/tau; if (len %% 1 == 0) len = len + 1;
	
	# j = seq(x[1],x[3],length.out=len);
	# if(length(j) == 1) return(matrix(x[c(1,3,2,4,5)],nrow=1));
		
	# k = seq(x[2],x[4],length.out = len);
	
	# j = rep(j,each=2); j = j[-c(1,length(j))];
	# j = matrix(j,ncol=2,byrow=TRUE);	
	# k = rep(k,each=2); k = k[-c(1,length(k))];
	# k = matrix(k,ncol=2,byrow=TRUE);	
	# l = matrix(rep(x[5],nrow(j)),ncol=1);
	# return(cbind(j,k,l));	
# }
