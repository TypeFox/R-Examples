## Probability distribution of a node and its parents
## and illustration of the pruning process 

setMethod("ppplot", signature="PSTf", 
	def=function(object, path, gain, C,  cex.plot=1, nsize=0.3, nlab=TRUE, psize=nsize/2, 
		pruned.col="red", div.col="green", ...) {

		A <- object@alphabet
		cpal <- c(object@cpal)

		if (has.cdata(object)) {
			c.A <-  alphabet(object@cdata)
			c.cpal <- cpal(object@cdata)
		} else { 
			c.A <- A 
			c.cpal <- cpal
		}
	
		nbstate <- length(A)
		oolist <- list(...)

		if (length(path)==1) { path <- seqdecomp(path) }

		sl <- length(path)

		## Retrieving probabilities
		prob <- matrix(nrow=nbstate, ncol=sl+1)
		N <- matrix(nrow=1, ncol=sl+1)

		node <- query(object, "e", output="all")
		prob[,(sl+1)] <- as.numeric(node@prob)
		N[,sl+1] <- node@n
	
		lsuff <- lsuffix(object, path)
		lsp <- sl-length(lsuff)+1	

		for (j in lsp:sl) {
			node <- query(object, paste(path[j:sl], collapse="-"), output="all")
			prob[,j] <- as.numeric(node@prob)
			N[,j] <- node@n

		}

		## Plotting path
		barw <- 1
		gsep <- 0.1
		ppsep <- psize/4
		poff <- 0

		nC <- if (!missing(C)) { length(C) } else { 0 }

		plot(NULL, 
			xlim=c(1-nsize, sl+2),
			ylim=c(0,(nsize+gsep+1+(nC*(psize+ppsep))+ ((nC>0)*gsep))),
			axes=FALSE,
			xlab="L (memory)", ylab="",
			...)

		## Tag as div 
		if (!missing(C)) {
			div <- matrix(nrow=nC, ncol=sl)
			pruned <- matrix(nrow=nC, ncol=sl)

			for (j in lsp:sl) {
				idpar <- 1
				if (gain=="G1") {
					for (i in 1:nC) { 
						div[idpar, j] <- G1(prob[,j], prob[,(j+1)], C=C[i])
						idpar <- idpar+1
					}
				} else if (gain=="G2") {
					for (i in 1:nC) { 
						div[idpar, j] <- G2(prob[,j], prob[,(j+1)], C=C[i], N=N[,j])
						idpar <- idpar+1
					}
				}
			}

			pruned[,lsp] <- !div[,lsp]
			for (j in (lsp+1):sl) {
				for (pp in 1:nC) {		
					pruned[pp, j] <- !div[pp, j] & pruned[pp, j-1]
				}
			}

			ppar.lab.pos <- NULL
			poff <- poff+(psize/2)

			for (pp in 1:nC) {
				segments(1, poff, sl+1, poff, col="grey", lwd=3)
	
				for (i in 1:sl) {
					pcol <- if (pruned[pp, i]) {pruned.col} else if ( div[pp,i] ) {div.col} else {"grey"}
					symbols(x=i, y=poff, circles=psize, bg=pcol, add=TRUE, inches=FALSE)
				}
	
				symbols(x=sl+1, y=poff,	circles=psize, bg="grey", add=TRUE, inches=FALSE)

				ppar.lab.pos <- c(ppar.lab.pos, poff)
				poff <- poff+psize+ppsep
			}

			ppar.lab <- paste("C", 1:nC, sep="")
			
			mtext(paste("C=", round(C,2)), side=2, at=ppar.lab.pos, 
				las=2, cex=cex.plot, adj=1)

		}

		## Plotting path and next symbol probability distributions
		poff <- poff+(nsize/2)
		prob.yBottom <- poff+(nsize/2)+gsep
	
		segments(1, poff, sl+1, poff, col="grey", lwd=3)

		for (i in 1:sl) {
			segments(i, poff, i, poff+(nsize/2)+gsep, col="grey", lwd=3)
	
			symbols(x=i, y=poff, circles=nsize, bg=c.cpal[which(path[i]==c.A)], add=TRUE, inches=FALSE)
			if (nlab) { text(x=i, y=poff, labels=path[i], cex=cex.plot) }
			plotProb(i-nsize, prob.yBottom , i+nsize, prob.yBottom+1, prob=t(prob[,i]), 
				cpal=cpal)
		}

		## ROOT NODE
		segments(sl+1, poff, sl+1, poff+(nsize/2)+gsep, col="grey", lwd=3)

		symbols(sl+1, y=poff, circles=nsize, bg="grey", add=TRUE, inches=FALSE)
		if (nlab) { text(x=sl+1, y=poff, labels="e") }

		## Plotting next symbol probability distributions
		plotProb((sl+1)-nsize, prob.yBottom, (sl+1)+nsize, prob.yBottom+1, 
			prob=t(prob[,sl+1]), cpal=cpal)

		axis(1, at=(1:(sl+1)), labels=sl:0, pos=-0.04, cex.axis=cex.plot)

		plabpos <- seq(from=prob.yBottom, to=(prob.yBottom+1), by=0.2)
		plab <- plabpos-prob.yBottom

		axis(2, at=plabpos, 
			labels=plab, 
			## las=2, 
			cex.axis=cex.plot)

		mtext("Prob", side=2, at=(max(plabpos)+min(plabpos))/2, line=3)
	}
)

