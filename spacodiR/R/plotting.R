#####
spacodi.treeplot <-function(spacodi.permutations, phy, cex=list(pch=1.5, tip=0.5, legend=0.75), transp=0.8, sig.plot=TRUE, cut.off=0.05, cols=list("white", "gray", "black"), main=TRUE, outfile=NULL, add.id=FALSE, ...){
#	op=par(no.readonly=TRUE)
#	on.exit(par(op))
	
	if(class(spacodi.permutations)=="list") base=as.data.frame(spacodi.permutations[[1]]) else base=as.data.frame(spacodi.permutations)
	if(nrow(base)!=phy$Nnode)stop("SPACoDi.plot requires data (whether NA not) for every node in the spacodi.permutations.") 
	
	base=base[order(base$node.ID),]
	
	
	sp.parm=strsplit(names(spacodi.permutations)[1],".",fixed=TRUE)[[1]][2]
	
	
	if(sig.plot==TRUE) {
		rand=as.data.frame(spacodi.permutations$randomization.test)
		if(length(rand)!=0) {
			rand=rand 
		} else {
			stop("Cannot perform significance plot without there being a randomization.test within the spacodi.permutations")
		}
	}
	
	if(any(names(base)==sp.parm)) {

		# SP.PARM PLOT: Ist, Pst, Bst, PIst, Tst, Ust, TAUst
		B=as.matrix(base[[sp.parm]])
		if(!is.null(outfile))pdf(outfile)
#		par(mfrow = c(1, ncol(B)))
		for(i in 1:ncol(B)){
			bst.i <- B[,i]
			if(sig.plot==TRUE) {
				col.i=array(dim=c(nrow(base), 2))
				col.i[,1]=as.numeric(base$node.ID)
				for(x in 1:nrow(base)){
					if(!is.na(match(base$node.ID[x],rand$node.ID)->foo.t)) {
						r.row=rand[which(rand$node.ID==base$node.ID[x]),]
						diff=(r.row[[paste("obs",sp.parm,sep=".")]]-r.row[[paste("m.exp",sp.parm,sep=".")]])
						if(!is.na(diff)) {
							if((diff > 0) && r.row$p.value <= cut.off) {
								col.i[x,2]=cols[[2]]
							} else if((diff < 0) && r.row$p.value <= cut.off) {
								col.i[x,2]=cols[[3]]
							}
						} else {
							col.i[x,2]=NA
						}
					}
					else {col.i[x,2]=NA}
				}
				col.i=as.vector(col.i[,2])
			} else {
				col.i <- 1+(bst.i-min(bst.i, na.rm = TRUE))*99/(max(bst.i, na.rm = TRUE)-min(bst.i, na.rm = TRUE))
				col.i <- round(col.i)
				cols  <- diverge_hcl(100)
				col.i <- cols[c(col.i, recursive = TRUE)]
				col.i[is.na(col.i)] <- "gray"
			}
			if(i != ncol(B)){
				plot(phy, show.tip.label = FALSE, no.margin = TRUE, cex=cex$tip, ...)
			} else {
				plot(phy, no.margin = TRUE, cex=cex$tip, ...)
			}
		    if(sig.plot) {
				nodelabels(bg=add.transparency(col.i,transp), col="gray", frame="circle", pch=21, cex=cex$pch)
				legend("bottomleft", legend=c(paste(sp.parm, ": as expected",sep=""), paste(sp.parm, ": larger than expected",sep=""), paste(sp.parm, ": smaller than expected",sep="")), fill=add.transparency(c(NA, cols[[2]],cols[[3]]),transp), bg = add.transparency("white", 0.5), box.lty="blank", border="gray", cex=cex$legend)
		    } else {
				nodelabels(bg=add.transparency(col.i,transp), col="gray", frame="circle", pch=21, cex=cex$pch)
				legend("bottomleft", legend = sprintf("%6.2f", rev(seq(min(bst.i, na.rm = TRUE), max(bst.i, na.rm = TRUE), length = 9))),  fill = add.transparency(rev(diverge_hcl(9)),transp), bg = add.transparency("white", 0.5), box.lty="blank", border="gray", cex=cex$legend)
		    }
			if(add.id==TRUE) {
				if(!is.null(phy$node.label)) phy$node.label=NULL
				nodelabels(adj=c(-0.5,-0.5),cex=0.6, frame="none")
			}
		}
		if(main)title(main = paste("Values of ",sp.parm,sep=""), outer = TRUE, line = -1)
		if(!is.null(outfile))dev.off()
	} else {stop("Cannot interpret spacodi.permutations for plotting.")}
#	invisible()
}	



#####
phy.dotplot<-function(sp.plot, phy, edge.width=0.2, lab.adj=c(0, 0), tips.adj=c(0.5, 3), tips.cex=1.0, pch=21, print.labs=FALSE, outfile=NULL,...) {
	op=par(no.readonly=TRUE)
	on.exit(par(op))
	if(!is.null(outfile))pdf(outfile)
	
	obj=sp.plot
	tips=phy$tip.label
	
	obj=as.phylocom(obj, picante=FALSE)
	obj=obj[which(obj$samples!=0),]
	plots=split(obj, obj$plot)
	n=ceiling(sqrt(length(plots)->l.plot))
	layout(matrix(1:n^2,n,n))
	if(n^2>16)warning("Tip labels are unlikely to plot properly with this many trees.")
	for(i in 1:l.plot) {
		tip.labs=vector()
		for(j in 1:length(tips)){
			if(!is.na(match(as.vector(tips[j]), as.vector(plots[[i]]$species)))) tip.labs[j]=1 else tip.labs[j]=0
		}
		plot(phy, no.margin = TRUE, show.tip.label = FALSE, direction="u", edge.width=edge.width, ...)
		tiplabels(pch=pch, bg=ifelse(tip.labs == 1, "black", NA), col=NA, adj=tips.adj, cex=tips.cex)
		legend("bottomleft", legend = names(plots)[i], fill = NULL, bg = NULL, box.lty="blank", border="white", adj=lab.adj)
	}
	if(!is.null(outfile))dev.off()
	if(print.labs) {res=list(plots, phy$tip.label); names(res)=c("groups", "tip.labels"); print(res); dev.new(); plot(phy,direction="u")}
	invisible()
}



#####
spacodi.permutplot<-function (spacodi.permutations, cex=list(pch=1.5, rand=0.1, node=0.5, legend=0.75), transp=0.8, col = list("black", 
"lightgray"), bg = list("white", "lightgray", "black"), all.points = TRUE, 
add.id = TRUE, sig.plot = TRUE, cut.off = 0.05, envelope = FALSE, 
outfile = NULL, ...) 
{
	
	sp.parm=strsplit(names(spacodi.permutations)[1],".",fixed=TRUE)[[1]][2]

#   op = par(no.readonly = TRUE)
#   on.exit(par(op))
    if (is.atomic(cex)) 
	cex = list(cex, cex/10, cex/2)
    o.orig = spacodi.permutations[[paste("observed",sp.parm,sep=".")]][which(!is.na(spacodi.permutations[[paste("observed",sp.parm,sep=".")]][[sp.parm]])), ]
    o = as.vector(o.orig[[sp.parm]])
    e = spacodi.permutations[[paste("expected",sp.parm,sep=".")]]
    lim.set <- function(obj, scl) {
        out = lapply(obj, function(x) {
					 foo = x + (x * scl)
					 return(foo)
					 })
        return(unlist(out))
    }
    x.lim = lim.set(c(min(-o.orig$node.time), max(-o.orig$node.time)), 
					0.05)
    if (all.points == TRUE) {
        range.array = array(dim = c(nrow(e), 2))
        for (rr in 1:nrow(e)) {
            r.foo = e[rr, which(!is.na(e[rr, ]))]
            if (length(r.foo) > 0) {
                range.array[rr, 1] = min(r.foo)
                range.array[rr, 2] = max(r.foo)
            }
            else {
                lapply(range.array[rr, ], function(obj) obj = NA)
            }
        }
        Bst.lim.y = c(l.b <- min(as.vector(o), as.vector(unlist(e)), na.rm = TRUE), u.b <- max(as.vector(o), as.vector(unlist(e)), na.rm = TRUE))
        y.scl = 0.35
        y.lim.exp = lim.set(Bst.lim.y, y.scl)
        y.lim = lim.set(c(min(o, na.rm = TRUE), max(o, na.rm = TRUE)), y.scl)
        if (diff(y.lim.exp) > diff(y.lim)) 
		warning(paste("Some expected ",sp.parm," were excluded from the plot.",sep=""))
    } else {
        y.lim = lim.set(c(min(o, na.rm = TRUE), max(o, na.rm = TRUE)), y.scl)
    }
    if (!is.null(outfile)) pdf(file = outfile)
    plot(-o.orig$node.time, o.orig[[sp.parm]], main = paste(sp.parm, " permutations through time", sep=""), 
		 xlab = "branching times", ylab = sp.parm, ylim = y.lim, 
		 xlim = c(x.lim[1], 0), type = "n", ...)
    if (all.points == TRUE) {
        e = as.data.frame(e)
        for (pp in 1:nrow(e)) {
            b.time = o.orig$node.time[which(o.orig$node.ID == row.names(e)[pp])]
            points(rep(-b.time, length(e[pp, ])), e[pp, ], cex = cex$rand, col = col[[2]], bg = bg[[2]], pch = 21)
        }
    }
	
# optional envelope plotting [DEPRECATE]
    if (envelope == TRUE) {
        poor = which(is.na(range.array[, 1]))
        if (length(poor) > 0) 
		min.array = range.array[-poor, 1]
        else min.array = range.array[, 1]
        if (length(poor) > 0) 
		max.array = range.array[-poor, 2]
        else max.array = range.array[, 2]
        if (length(poor) > 0) 
		o.nodes = -o.orig$node.time[-poor]
        else o.nodes = -o.orig$node.time
        lines(smooth.spline(o.nodes, min.array), lty = 2)
        lines(smooth.spline(o.nodes, max.array), lty = 2)
    }
	
# color nodes by significance of diversity partitioning
    if (sig.plot == TRUE) {
        if (any(names(spacodi.permutations) == "randomization.test")) {
            rr = spacodi.permutations$randomization.test
        } else {
            stop("Must have randomization.test results supplied with spacodi.permutations.")
        }
        bg.sig = vector()
        if (length(bg) == 3) 
		bg.cols = bg
        else bg.cols = list("white", "black", "gray")
        for (node in 1:nrow(rr)) {
            if (rr[node, "p.value"] > cut.off || is.na(rr[node, "p.value"])) {
                bg.sig[node] = bg.cols[[1]]
            } else if (rr[node, "p.value"] <= cut.off && ((rr[node, paste("obs",sp.parm,sep=".")] - rr[node, paste("m.exp",sp.parm,sep=".")]) > 0)) {
                bg.sig[node] = bg.cols[[2]]
            } else {
                bg.sig[node] = bg.cols[[3]]
            }
        }
        points(-o.orig$node.time, o.orig[[sp.parm]], cex = cex$pch, col = col[[1]], pch = 21, bg = add.transparency(bg.sig, transp), ...)
        legend("topleft", legend = c(paste(sp.parm,": as expected",sep=""), paste(sp.parm, ": larger than expected", sep=""), 
									paste(sp.parm, ": smaller than expected", sep="")), fill = add.transparency(c(bg.cols[[1]], 
									bg.cols[[2]], bg.cols[[3]]),transp), bg = add.transparency("white", 0.5), box.lty = "blank", 
									border = "gray", cex = cex$legend)
    } else if(sig.plot == FALSE) {
        points(-o.orig$node.time, o.orig[[sp.parm]], cex = cex$pch,  col = col[[1]], pch = 21, bg = add.transparency(bg[[1]],transp), ...)
    }
	
# add node labels
    if (add.id == TRUE) {
        textxy(-o.orig$node.time, o.orig[[sp.parm]], o.orig$node.ID, 
			   m = c(mean(-o.orig$node.time), mean(o.orig[[sp.parm]])), 
			   cx = cex$node)
    }
    if (!is.null(outfile)) 
	dev.off()
#   invisible()
}

#####
add.transparency <-
function (col, alpha) {
    tmp <- col2rgb(col)/255
    rgb(tmp[1, ], tmp[2, ], tmp[3, ], alpha = alpha)
}

#####
textxy <-
function (X, Y, labs, cx = 0.5, dcol = "black", m = c(0, 0)) 
{
    posXposY <- ((X >= m[1]) & ((Y >= m[2])))
    posXnegY <- ((X >= m[1]) & ((Y < m[2])))
    negXposY <- ((X < m[1]) & ((Y >= m[2])))
    negXnegY <- ((X < m[1]) & ((Y < m[2])))
    if (sum(posXposY) > 0) 
	text(X[posXposY], Y[posXposY], labs[posXposY], adj = c(-0.3,  -0.3), cex = cx, col = dcol)
    if (sum(posXnegY) > 0) 
	text(X[posXnegY], Y[posXnegY], labs[posXnegY], adj = c(-0.3,  1.3), cex = cx, col = dcol)
    if (sum(negXposY) > 0) 
	text(X[negXposY], Y[negXposY], labs[negXposY], adj = c(1.3,  -0.3), cex = cx, col = dcol)
    if (sum(negXnegY) > 0) 
	text(X[negXnegY], Y[negXnegY], labs[negXnegY], adj = c(1.3,  1.3), cex = cx, col = dcol)
}


