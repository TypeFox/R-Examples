mkdtsegsPolar <- function(x, tau, edge) {
	tol <- 0.0001;
	xremainder <- numeric(max(edge[,1]));
	dtsegs <- vector("list", nrow(edge));
	for (i in 1:nrow(x)) {
		xtol <- abs(tol*cos(x[i,5]));
		d <- sqrt((x[i,3]-x[i,1])^2+(x[i,4]-x[i,2])^2);
		if (xremainder[edge[i,1]] > 0) {
			if ((x[i,5] > pi/2 && x[i,5] < 3*pi/2) && x[i,5] != pi) {
				if (x[i,1]+xremainder[edge[i,1]]*cos(x[i,5]) < x[i,3]) {
					xremainder[edge[i,2]] <- abs(xremainder[edge[i,1]]-d);
					xx <- x[i,1];
					yy <- x[i,2];
				}
				else {
					xx <- seq(x[i,1]+xremainder[edge[i,1]]*cos(x[i,5]),x[i,3],tau*cos(x[i,5]));
					yy <- xx*tan(x[i,5]);
					xx <- c(x[i,1],xx);
					yy <- c(x[i,2],yy);
				}
			}
			else if ((x[i,5] < pi/2 || x[i,5] > 3*pi/2) && x[i,5] != 0) {
				if (x[i,1]+xremainder[edge[i,1]]*cos(x[i,5]) > x[i,3]) {
					xremainder[edge[i,2]] <- abs(xremainder[edge[i,1]]-d);
					xx <- x[i,1];
					yy <- x[i,2];
				}
				else {
					xx <- seq(x[i,1]+xremainder[edge[i,1]]*cos(x[i,5]),x[i,3],tau*cos(x[i,5]));
					yy <- xx*tan(x[i,5]);
					xx <- c(x[i,1],xx);
					yy <- c(x[i,2],yy);
				}
			}
			else if (x[i,5] == pi/2) {
				if (x[i,2]+xremainder[edge[i,1]] > x[i,4]) {
					xremainder[edge[i,2]] <- abs(xremainder[edge[i,1]]-d);
					xx <- x[i,1];
					yy <- x[i,2];
				}
				else {
					yy <- seq(x[i,2]+xremainder[edge[i,1]],x[i,4],tau);
					yy <- c(x[i,2],yy);
					xx <- rep(x[i,1],length(yy));
				}
			}		
			else if (x[i,5] == 3*pi/2) {
				if (x[i,2]-xremainder[edge[i,1]] < x[i,4]) {
					xremainder[edge[i,2]] <- abs(xremainder[edge[i,1]]-d);
					xx <- x[i,1];
					yy <- x[i,2];
				}
				else {
					yy <- seq(x[i,2]-xremainder[edge[i,1]],x[i,4],-tau);
					yy <- c(x[i,2],yy);
					xx <- rep(x[i,1],length(yy));
				}
			}		
			else if (x[i,5] == pi) {
				if (x[i,1]-xremainder[edge[i,1]] < x[i,3]) {
					xremainder[edge[i,2]] <- abs(xremainder[edge[i,1]]-d);
					xx <- x[i,1];
					yy <- x[i,2];
				}
				else {	
					xx <- seq(x[i,1]-xremainder[edge[i,1]],x[i,3],-tau);
					xx <- c(x[i,1],xx);
					yy <- rep(x[i,2],length(xx));
				}
			}		
			else {
				if (x[i,1]+xremainder[edge[i,1]] > x[i,3]) {
					xremainder[edge[i,2]] <- abs(xremainder[edge[i,1]]-d);
					xx <- x[i,1];
					yy <- x[i,2];
				}
				else {	
					xx <- seq(x[i,1]+xremainder[edge[i,1]],x[i,3],tau);
					xx <- c(x[i,1],xx);
					yy <- rep(x[i,2],length(xx));
				}
			}
		}		
		else {
			if (x[i,5] == 0 || x[i,5] == pi) {
				xx <- seq(x[i,1],x[i,3],tau*cos(x[i,5]));
				yy <- rep(x[i,2],length(xx));
			}
			else if (x[i,5] == pi/2 || x[i,5] == 3*pi/2) {
				yy <- seq(x[i,2],x[i,4],tau*sin(x[i,5]));
				xx <- rep(x[i,1],length(yy));
			}
			else {
				xx <- seq(x[i,1],x[i,3],tau*cos(x[i,5]));
				yy <- xx*tan(x[i,5]);	
			}
		}		
		if (length(xx) > 1) {
			d <- sqrt((x[i,3]-tail(xx,1))^2+(x[i,4]-tail(yy,1))^2);
			if (abs(x[i,3] - tail(xx,1)) > xtol) {
				xremainder[edge[i,2]] <- abs(tau-d);
				xx <- c(xx,x[i,3]);
				yy <- c(yy,x[i,4]);
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
			d <- sqrt((x[i,3]-tail(xx,1))^2+(x[i,4]-tail(yy,1))^2);
			if (xremainder[edge[i,1]] == 0) {
				xremainder[edge[i,2]] <- abs(tau-d);
			}
			segs <- matrix(c(x[i,1],x[i,3],x[i,2],x[i,4],edge[i,2]),nrow=1,ncol=5);
		}
		dtsegs[[i]] <- segs;
	}
	return(do.call(rbind, dtsegs));
}
