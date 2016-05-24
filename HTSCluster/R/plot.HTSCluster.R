`plot.HTSCluster` <-
function (x,file.name = FALSE, 
graphs = c("map", "map.bycluster", "lambda"), data=NA, ...) 
{
	
    	if (class(x) != "HTSCluster") {
        	stop(paste(sQuote("x"), sep = ""), " must be of class ", 
            paste(dQuote("HTSCluster"), sep = ""), sep = "")
    	}
#	if(cluster.choice != "ICL" & cluster.choice != "BIC" & length(cluster.choice) > 1) {
#		stop(paste(sQuote("cluster.choice"), sep = ""), " must be one of ",
#			paste(dQuote("ICL"), sep = ""), " or ", paste(dQuote("BIC"), sep = ""), 
#			" or \n one of the cluster sizes run for object ", paste(sQuote("x"), sep = ""))
#	}

	conds <- x$conds;
	BIC <- x$BIC
	ICL <- x$ICL
	labels <- x$labels
	probaPost <- x$probaPost
	lambda <- x$lambda
	pi <- x$pi
	g <- length(pi)

	## Plot dimensions
#	if(g > 36) {
#		stop("Visualization of individual clusters not yet possible for a 
#			large number of clusters (> 36).")
#	}

	if(g == 1) plot.dim <- c(1,1)
	if(g == 2) plot.dim <- c(1,2)
	if(g == 3 | g == 4) plot.dim <- c(2,2)
	if(g == 5 | g == 6) plot.dim <- c(2,3)
	if(g == 7 | g == 8) plot.dim <- c(2,4)
	if(g == 9 | g == 10) plot.dim <- c(2,5)
	if(g == 11 | g == 12) plot.dim <- c(3,4)
	if(g == 13 | g == 14 | g == 15 | g == 16) plot.dim <- c(4,4)
	if(g == 17 | g == 18 | g == 19 | g == 20) plot.dim <- c(4,5)
	if(g > 20 & g <= 25) plot.dim <- c(5,5)
	if(g > 25 & g <= 30) plot.dim <- c(5,6)
	if(g > 30 & g <= 36) plot.dim <- c(6,6)

	if(file.name != FALSE) pdf(paste(file.name));

	map <- apply(matrix(probaPost, ncol=g), 1, max)
	if("map" %in% graphs) {
	## MAP histogram for all clusters
	if(file.name == FALSE) par(ask = TRUE);
	par(mfrow = c(1,1))
	hist(map, col = "grey", breaks = 50, xlab = "Max conditional probability",
		main = "Max conditional probability: All clusters", cex.axis = 1.25, cex.lab = 1.25,
		cex.main = 1.5, font = 2)
	}

	## MAP histogram for each individual cluster
	if("map.bycluster" %in% graphs) {
	if(file.name == FALSE) par(ask = TRUE);
	par(mfrow = c(1,1), mar = c(4,4,2,2), oma = c(2,1,2,0))
	for(i in 1:g) {
		if(sum(labels == i) > 1) {
			map.clust <- apply(matrix(probaPost[labels == i,], ncol=g), 1, max)
		}
		if(sum(labels == i) == 1) {
			map.clust <- max(probaPost[labels == i,])
		}
		if(sum(labels == i) == 0) {
			map.clust <- NA
		}
		if(i == 1) {
			boxplot(map.clust, at = 1, xlim = c(0.5,g+.5), col = "grey",
				ylim = c(0,1), width = .35)
		}
		if(i > 1 & is.na(map.clust[1]) == FALSE) {
			boxplot(map.clust, at = i, col = "grey", add = TRUE, width = .35)
		}
	}
	mtext(side = 1, at = 1:g, 1:g, line = 1)
	mtext(side = 3, outer = TRUE, "Max conditional probability: By cluster", line = 0, font = 2, cex = 1.25)
	mtext(side = 2, outer = TRUE, "Max conditional probability", line = 0)
	mtext(side = 1, outer = TRUE, "Cluster", line = -1)
	}

	## Weighted histograms
	if("weighted.histograms" %in% graphs) {
	if(is.na(data) == TRUE) stop("Must provide original data for weighted histograms to be plotted.")
	n <- nrow(data)
	data.w <- data/rowSums(data)
	if(file.name == FALSE) par(ask = TRUE);
	if(g > 36) {
		stop("Weighted histograms not yet possible for a 
			large number of clusters (> 36).")
	}
	for(var in 1:dim(data)[2]) {
		par(mfrow = plot.dim, mar = c(4,4,2,2), oma = c(1,1,2,0)) 
		for(clus in 1:g) {
	
			if(sum(labels == clus) == 0) {
				plot(c(0,1),c(0,1), col = "white", bty = "n", xaxt = "n",
					yaxt = "n", xlab = "", ylab = "")
				mtext(side = 3, paste("Cluster ", clus, sep = ""), line = 0, 
					cex = 1)
				next;
			}
			wh <- weighted.hist(data.w[,var], w = probaPost[,clus], 
				breaks = 50, freq = F, plot = F)
			breaks <- wh$breaks
			mids <- wh$mids
			data.w.sub <- data.w[labels == clus,var]
			o <- order(data.w.sub)
			ref <- rep(0, length(breaks)-1)
			for(br in 2:length(breaks)) {
				ref[br-1] <- length(which(data.w.sub >= breaks[br-1] & 
					data.w.sub < breaks[br]))
			}
			ref <- ref/length(data.w.sub)
			ymax <- max(wh$density, ref)
			plot(mids-min(breaks), ref, col = "white", ylim = c(0,ymax), ylab = "Density",
				xlab = expression(paste(y[i*j*l], "/", w[i], sep = "")), xaxt = "n",
				cex.axis = 1.1, cex.lab = 1.2, yaxt = "n",
				xlim = c(0-min(breaks), max(breaks)-min(breaks)))
			weighted.hist(data.w[,var], w = probaPost[,clus], ylab = "", 
				breaks = 50, freq = F, col = "grey", add = T)
			points(mids-min(breaks), ref, col = "red", pch = 19, cex = .8)
			lines(mids-min(breaks), ref, col = "red")
			mtext(side = 3, paste("Cluster ", clus, sep = ""), line = 0, 
			cex = 1)
		}
		mtext(side = 3, outer = T, line = 0, paste("Weighted histograms: Variable", var), cex = 1,
			font = 2)
	}
	}

	## Lambda and pi plot
	if("lambda" %in% graphs) {
	if(file.name == FALSE) par(ask = TRUE);
	if(length(unique(conds)) == 1) {
		par(oma = c(0,1,0,0),mfrow = c(1,1))
		bar <- barplot(unlist(lambda[1,]), ylim = c(0, max(lambda)+.25), 
			col = "lightgrey",
			width = unlist(pi), names.arg = 1:g, xlab = "Cluster", 
			cex.lab = 1.25, cex.names = 1.25, yaxt = "n", main = "", cex.main = 1.5)
		axis(side = 2, at = seq(from = 0, to = 2.5, by = 0.5), 
			labels = seq(from = 0, to = 2.5, by = 0.5))
		mtext(side = 2, line = 3, at = 1, expression(hat(lambda)[1*k]), cex = 1.5)
		abline(h = 0, lwd = 2)
		abline(h = 1, lty = 2)
		pi.lab <- paste(round(unlist(pi)*100,1), "%", sep = "")
		text(bar, unlist(lambda[1,]) + .15, pi.lab, xpd = TRUE, col = "blue")
		mtext(side = 3, outer = F, line = 0,
			expression(paste("Relative values of ", (lambda)[k], 
			", by cluster", sep = "")), cex = 1.25,
			font = 2)
	}
#	if(length(unique(conds)) == 2) {
#		par(oma = c(0,1,0,0), mfrow = c(1,1))
#		bar <- barplot(unlist(lambda[1,]), ylim = c(-(max(lambda)+.25), max(lambda)+.25), 
#			col = "lightgrey",
#			width = unlist(pi), names.arg = 1:g, xlab = "Cluster", 
#			cex.lab = 1.25, cex.names = 1.25, yaxt = "n", main = "", cex.main = 1.5)
#		axis(side = 2, at = seq(from = -2.5, to = 2.5, by = 0.5), 
#			labels = c(seq(from = 2.5, to = 0, by = -0.5), seq(from = 0.5, to = 2.5, 
#			by = 0.5)))
#		barplot(unlist(-lambda[2,]), add = T, col = "darkgrey", width = unlist(pi), 
#			names.arg = 1:g, cex.names = 1.25, yaxt = "n")
#		mtext(side = 2, line = 3, at = 1, expression(hat(lambda)[1*k]), cex = 1.5)
#		mtext(side = 2, line = 3, at = -1, expression(hat(lambda)[2*k]), cex = 1.5)
#		abline(h = 0, lwd = 2)
#		abline(h = c(-1,1), lty = 2)
#		pi.lab <- paste(round(unlist(pi)*100,1), "%", sep = "")
#		text(bar, unlist(lambda[1,]) + .15, pi.lab, xpd = TRUE, col = "blue")
#		mtext(side = 3, outer = F, line = 0,
#			expression(paste("Relative values of ", (lambda)[k], 
#			", by cluster", sep = "")), cex = 1.25,
#			font = 2)
#	}
#	if(length(unique(conds)) > 2 & lambda.plot == "stars") {
#
#		par(oma = c(3*dim(lambda)[1],0,0,0), mfrow = c(1,1))
#		lambda.scale <- t(lambda) / rowSums(t(lambda))
#		rownames(lambda.scale) <- paste(rownames(lambda.scale),
#			" (", round(pi*100,2), "%", ")", sep = "")
#		st <- stars(lambda.scale, draw.segments = TRUE, scale = FALSE,
#			mar = c(0,0,3,0), full = FALSE,
#			col.segments = 1:dim(lambda)[1], 
#			main = expression(paste("Relative values of ", (lambda)[k], 
#			", by cluster", sep = "")), cex = 1.25, len = pi/max(pi))
##		text(st[,1], st[,2]+.55, paste("(", round(pi*100,2), "%", ")", sep = ""))
#		par(xpd=NA)
#		legend(mean(par("usr")[1:2]-.5),mean(par("usr")[3]-.1), col = 1:dim(lambda)[1], 
#			paste("Condition", unique(conds)), pch = 19)
#
#	}
#	if(length(unique(conds)) > 2 & lambda.plot == "pie") {
#		df <- data.frame(x = as.vector(lambda), lambda=rep(rownames(lambda),g),
#                 cluster = rep( paste("Cluster ", 1:g, "\npi = ", round(pi,2),  
#			"\n(", table(labels), " obs)", sep = ""),
#                   	each = nrow(lambda)))
#		df$cluster <- factor(df$cluster, 
#			levels = rep( paste("Cluster ", 1:g, "\npi = ", round(pi,2), 
#			"\n(", table(labels), " obs)", sep = "")))
#		levels(df$lambda) <- rownames(lambda)
#		gg <- ggplot(df, aes(x = lambda, fill = factor(lambda))) + 
#			geom_bar(aes(weight=x)) + theme(legend.position="none") +
# 			scale_fill_manual(name = "Group", values = rev(brewer.pal(nrow(lambda), "Set3"))) + 
#			facet_wrap(~cluster) + coord_polar() +
#  			labs(x = "", y = "") + theme_bw() + theme(axis.text.x = element_blank())
#		print(gg)
#	}
	if(length(unique(conds)) >= 2) {
		s <- x$s
		ss <- rep(NA, length(unique(conds))) 
		for(c in 1:length(unique(conds))) {
			ss[c] <- sum(s[which(conds == unique(conds)[c])])
		}
		b <- barplot(lambda*ss, beside=FALSE, col=1:length(unique(conds)), width = pi, legend.text = rownames(lambda),
        		args.legend = list(y=-0.1, horiz=TRUE, bty="n"),
				names.arg=rep("", ncol(lambda)), space=0)
		mtext(side=2, outer=F, expression(paste(lambda[jk], s[j.])), 
			line=2)
		mtext(side=1, outer=F, "Cluster", line=2.5)
		mtext(side=1, at=b, line = 0, 1:ncol(lambda), cex = 0.75)
		if(length(unique(conds)) == 2) abline(h=0.5, lty=2, col="grey")
	}
      }
	if(file.name != FALSE)  dev.off();
      }
