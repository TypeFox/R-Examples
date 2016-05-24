
getBestShiftConfiguration <- function(x, expectedNumberOfShifts , threshold = 5){
	
	if (class(x) == 'bammdata') {
		x <- credibleShiftSet(x, expectedNumberOfShifts, threshold, set.limit = 0.95);	
	} else if (class(x) == 'credibleshiftset') {

	} else {
		stop("Argument x must be of class bammdata or credibleshiftset\n");
	}
	
	class(x) <- 'bammdata';	
	subb <- subsetEventData(x, index=x$indices[[1]]);
	
	# Drop all non-core shifts after adding root:
	coreshifts <- c((length(x$tip.label) + 1), x$coreshifts);
	coreshifts <- intersect(subb$eventData[[1]]$node, coreshifts); 
	 
		
 	for (i in 1:length(subb$eventData)){
 		#subb$eventData[[i]] <- subb$eventData[[i]][subb$eventData[[i]]$node %in% coreshifts,];
 		if (i == 1){
 			ff <- subb$eventData[[i]];
 		}
 		ff <- rbind(ff, subb$eventData[[i]]);
 	}

	xn <- numeric(length(coreshifts));
	xc <- character(length(coreshifts));
	
	if (x$type == 'diversification'){
		dff <- data.frame(generation = xn, leftchild=xc, rightchild=xc, abstime=xn, lambdainit=xn, lambdashift=xn, muinit = xn, mushift = xn, stringsAsFactors=F);	
		for (i in 1:length(coreshifts)){
			if (coreshifts[i] <= length(x$tip.label)){
			# Node is terminal:	
				dset <- c(x$tip.label[coreshifts[i]], NA)
				
			}else{
				# node is internal.
				tmp <- extract.clade(as.phylo(x), node= coreshifts[i]);
				dset <- tmp$tip.label[c(1, length(tmp$tip.label))];		
			}
			
			tmp2 <- ff[ff$node == coreshifts[i], ];
			
			dff$leftchild[i] <- dset[1];
			dff$rightchild[i] <- dset[2];
			dff$abstime[i] <- mean(tmp2$time);
			dff$lambdainit[i] <- mean(tmp2$lam1);
			dff$lambdashift[i] <- mean(tmp2$lam2);
			dff$muinit[i] <- mean(tmp2$mu1);
			dff$mushift[i] <- mean(tmp2$mu2);
		}	
		best_ed <- getEventData(as.phylo(x), eventdata=dff);	
	}else if (x$type == 'trait'){
		dff <- data.frame(generation = xn, leftchild=xc, rightchild=xc, abstime=xn, betainit=xn, betashift=xn, stringsAsFactors=F);					
		for (i in 1:length(coreshifts)){
			if (coreshifts[i] <= length(x$tip.label)){
			# Node is terminal:	
				dset <- c(x$tip.label[coreshifts[i]], NA)
				
			}else{
				# node is internal.
				tmp <- extract.clade(as.phylo(x), node= coreshifts[i]);
				dset <- tmp$tip.label[c(1, length(tmp$tip.label))];		
			}
			
			tmp2 <- ff[ff$node == coreshifts[i], ];
			
			dff$leftchild[i] <- dset[1];
			dff$rightchild[i] <- dset[2];
			dff$abstime[i] <- mean(tmp2$time);
			dff$betainit[i] <- mean(tmp2$lam1);
			dff$betashift[i] <- mean(tmp2$lam2);

		}			
		best_ed <- getEventData(as.phylo(x), eventdata=dff, type = 'trait');	
	}else{
		stop("error in getBestShiftConfiguration; invalid type");
	}

	return(best_ed);
 
}

