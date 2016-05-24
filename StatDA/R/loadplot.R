"loadplot" <-
function(fa.object, titlepl = "Factor Analysis", crit = 0.3, length.varnames=2)
{
#Reimann-plot of loadings matrix
# crit...plot all loadings > abs(crit)
# length.varnames ... number of letters for abbreviating the variable names in the plot
tplace=-0.5 # where to place text labels
	l <- fa.object$loa
	lnam <- dimnames(l)[[1]]
	k <- dim(l)[2]
	p <- dim(l)[1]
	for(j in 1:p) {
		lnam[j] <- substring(lnam[j], 1, length.varnames)
	}
	plot(cbind(c(0, 1, 1, 0, 0), c(-1, -1, 1, 1, -1)), type = "l", axes = FALSE, 
		xlab = "", ylab = "")
	segments(0, 0, 1, 0)
	segments(0, 0.5, 1, 0.5, lty = 2)
	segments(0, -0.5, 1, -0.5, lty = 2)
	bb <- apply(l^2, 2, sum)/sum(l^2)
	bb1 <- cumsum(bb)
	aa <- eigen(fa.object$loa %*% cor(fa.object$sco) %*% t(fa.object$loa) + 
		diag(fa.object$uni))$val
	aa1 <- cumsum(aa/sum(aa))
	#text(0, 1.1, "0%")
	mtext("0%",side=3,at=0,line=tplace)
	for(i in 1:k) {
		segments(bb1[i], -1, bb1[i], 1)
		lcrit <- abs(l[, i]) > crit
		lsel <- l[lcrit, i]
		names(lsel) <- lnam[lcrit]
		#text(bb1[i], 1.1, paste(round(100 * aa1[i]), "%", sep = ""))
		mtext(paste(round(100 * aa1[i]), "%", sep = ""),side=3,at=bb1[i],line=tplace)
		if(i == 1) {
#text(runif(length(lsel), 0, bb1[1]), lsel, 
#	names(lsel))
			chardist <- bb[1]/(length(lsel) + 1)
			text(seq(from = chardist, by = chardist, length = length(
				lsel)), lsel, names(lsel))
			#text(bb1[i]/2, -1.1, paste("F", round(i), sep = ""))
			mtext(paste("F", round(i), sep = ""),side=1,at=bb1[i]/2,line=tplace)
		}
		else {
#text(runif(length(lsel), bb1[i - 1], bb1[i]), 
#	lsel, names(lsel))
			chardist <- (bb1[i] - bb1[i - 1])/(length(lsel) + 1)
			text(seq(from = bb1[i - 1] + chardist, by = chardist, 
				length = length(lsel)), lsel, names(lsel))
			#text(bb[i]/2 + bb1[i - 1], -1.1, paste("F", round(i), sep = ""))
			mtext(paste("F", round(i), sep = ""),side=1,at=bb[i]/2 + bb1[i - 1],line=tplace)
		}
	}
#axis(2,at=a<-c(-1,-0.5,0,0.5,1),labels=a)
par(las=1)
	mtext("-1",side=2,at=-1,line=tplace)
	mtext("-0.5",side=2,at=-0.5,line=tplace)
	mtext("0",side=2,at=0,line=tplace)
	mtext("+0.5",side=2,at=0.5,line=tplace)
	mtext("+1",side=2,at=1,line=tplace)
	#text(-0.05, -1, "-1")
	#text(-0.05, -0.5, "-0.5")
	#text(-0.05, 0, "0")
	#text(-0.05, 0.5, "+0.5")
	#text(-0.05, 1, "+1")
	title(titlepl)
}
