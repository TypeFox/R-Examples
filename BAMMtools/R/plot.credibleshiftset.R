plot.credibleshiftset <- function(x, plotmax=9, method='phylogram', pal = 'RdYlBu', shiftColor = 'black', spex = "s", add.freq.text = TRUE, use.plot.bammdata = TRUE, border = TRUE, legend = FALSE, send2pdf = FALSE, logcolor=FALSE, breaksmethod='linear', color.interval=NULL, JenksSubset=20000, ...)
{
	if (class(x) != "credibleshiftset") {
		stop('arg x must be of class "credibleshiftset"');
	}
	if (!spex %in% c('s', 'e', 'netdiv')) {
		stop("arg spex must be 's', 'e' or 'netdiv'. ")
	}
	if ((spex == "e" || spex == "se") && x$type == "trait") {
		warning("arg spex not meaningful for BAMMtrait");
		spex <- "s";
	}
	cset.bamm <- as.bammdata(x);
	if (plotmax > 9 && send2pdf == FALSE) {
	    plotmax <- 9;
	    cat("arg plotmax coerced to 9\n");
	}
	mm <- min(x$number.distinct, plotmax);
	if (send2pdf) {
	    pdf("credibleshiftset.pdf");
	}
	if (mm == 1) {
		if (legend) {
			m <- matrix(c(1,2),byrow=TRUE,nrow=1,ncol=2);
			layout(m,respect=TRUE);
		} else {
			par(mfrow=c(1,1));
		}
	} else if (mm <= 2) {
		if (legend) {
			m <- matrix(c(1,2,3),byrow=TRUE,nrow=1,ncol=3);
			layout(m,respect=TRUE);	
		} else {
			par(mfrow=c(1,2));
		}		
	} else if (mm <= 4) {
		if (legend) {
			m <- matrix(c(1,2,0,3,4,5), byrow=TRUE, nrow=2,ncol=3);
			layout(m,respect=TRUE);	
		} else {
			par(mfrow = c(2,2));
		}
	} else if (mm <= 6) {
		if (legend) {
			m <- matrix(c(1,2,0,3,4,7,5,6,0),byrow=TRUE,nrow=3,ncol=3);
			layout(m,respect=TRUE);
		} else {
			par(mfrow=c(2,3));
		}
	} else {
		if (legend) {
			m <- matrix(c(1,2,3,0,4,5,6,10,7,8,9,0),byrow=TRUE,nrow=3,ncol=4);
			layout(m,respect=TRUE);
		} else {
			par(mfrow=c(3,3));	
		}
	}
	cat("Omitted", max(x$number.distinct,mm) - min(x$number.distinct,mm), "plots\n");
	if (use.plot.bammdata) {
    	cset.bamm <- dtRates(cset.bamm, 0.01);
	    colorbreaks <- assignColorBreaks(cset.bamm$dtrates$rates,spex=spex, logcolor=logcolor, method=breaksmethod, JenksSubset=JenksSubset);
	}
	for (i in 1:mm) {
	    sed <- subsetEventData(cset.bamm, index=x$indices[[i]]);
		par(mar = c(2,2,2,2));
		if (use.plot.bammdata) {
            plot.bammdata(sed, method=method, pal=pal, spex=spex, colorbreaks=colorbreaks, par.reset=FALSE, logcolor=logcolor, ...);
		}
		else {
		    if (method=="polar") method = "fan";
		    plot.phylo(as.phylo.bammdata(cset.bamm),type=method,show.tip.label=FALSE);
		}
		if (add.freq.text) mtext(sprintf("f = %.2g",x$frequency[i]),3);	
		if (border) box();
		#shiftnodes <- getShiftNodesFromIndex(cset.bamm, i);
		shiftnodes <- x$shiftnodes[[i]];
		# shiftnode_parents <- cset.bamm$edge[match(shiftnodes, cset.bamm$edge[,2],nomatch=0), 1];
	 #    root <- (shiftnode_parents == (cset.bamm$Nnode + 2));
  #   	if (sum(root) > 0) {
  #       	isShiftNodeParent <- integer(length(shiftnodes));
  #       	isShiftNodeParent[root] <- 1;
  #       	isShiftNodeParent[!root] <- sed$eventVectors[[1]][match(shiftnode_parents[!root], cset.bamm$edge[,2])];
  #   	} 
  #   	else {
  #   	    isShiftNodeParent <- sed$eventVectors[[1]][match(shiftnode_parents, cset.bamm$edge[,2])];
  #   	}
  #   	isShiftNode <- match(shiftnodes, sed$eventData[[1]]$node);
  #   	time <- sed$eventData[[1]][isShiftNode, 2] - sed$eventData[[1]][isShiftNodeParent, 2];
  #   	if (spex == "s") {
  #   		lam1 <- sed$eventData[[1]][isShiftNodeParent, 3]; 
  #   		lam2 <- sed$eventData[[1]][isShiftNodeParent, 4];
  #   	    AcDc <- exponentialRate(time, lam1, lam2) > sed$eventData[[1]][isShiftNode, 3];	
  #   	}
  #   	else if (spex == "e") {
  #   		mu1 <- sed$eventData[[1]][isShiftNodeParent, 5]; 
  #   		mu2 <- sed$eventData[[1]][isShiftNodeParent, 6];
  #   		AcDc <- exponentialRate(time, mu1, mu2) > sed$eventData[[1]][isShiftNode, 5];
  #   	}
  #   	else {
  #   		lam1 <- sed$eventData[[1]][isShiftNodeParent, 3]; 
  #   		lam2 <- sed$eventData[[1]][isShiftNodeParent, 4];
  #   		mu1 <- sed$eventData[[1]][isShiftNodeParent, 5]; 
  #   		mu2 <- sed$eventData[[1]][isShiftNodeParent, 6];
  #   		AcDc <- (exponentialRate(time, lam1, lam2)-exponentialRate(time, mu1, mu2)) > (sed$eventData[[1]][isShiftNode, 3]-sed$eventData[[1]][isShiftNode, 5]);
  #  		}
  #   	bg <- rep("blue", length(AcDc));
  #   	bg[which(AcDc == FALSE)] <- "red";

  		bg <- rep(shiftColor, length(shiftnodes));
		cex <- 0.75 + 5 * x$marg.probs[as.character(shiftnodes)];
		if (use.plot.bammdata) {
			cex <- cex[match(sed$eventData[[1]]$node, shiftnodes, nomatch=0)];
			bg <- bg[match(sed$eventData[[1]]$node, shiftnodes, nomatch=0)];
			shiftnodes <- shiftnodes[match(sed$eventData[[1]]$node, shiftnodes, nomatch=0)];
			addBAMMshifts(sed, 1, method, cex=cex, bg=transparentColor(bg, 0.5), shiftnodes = shiftnodes, par.reset=FALSE);	
		}
		else {
			r <- cset.bamm$edge.length[match(shiftnodes,cset.bamm$edge[,2])];
			lastPP <- get("last_plot.phylo",envir=.PlotPhyloEnv);
			XX <- lastPP$xx[shiftnodes];
			YY <- lastPP$yy[shiftnodes];
			if (method == "phylogram") {
				XX <- XX - 0.5*r;
			}
			else {
				theta <- atan2(YY, XX);
				XX <- XX - 0.5*r*cos(theta);
				YY <- YY - 0.5*r*sin(theta);
			}
			points(XX, YY, pch = 21, cex = cex, col = 1, bg = transparentColor(bg,0.5));
		}
	}
	if (legend) {
		par(mar=c(2.1,1.1,2.1,2.1));
		plot.new();			
		plot.window(xlim=c(0,1),ylim=c(0,1.25));
		points(rep(0.5,8),rev(seq(0,1.125,length.out=8)),pch=21,bg="white",cex = rev(0.75 + 5*seq(0.125,1,0.125)));
		mtext("marginal shift probability",cex=0.75,line=0);
		text(x=rep(0.8,8),y=rev(seq(0,1.125,length.out=8)),labels=sprintf("%10.3f",rev(seq(0.125,1,0.125))));
	}
	if (send2pdf) dev.off();
}
