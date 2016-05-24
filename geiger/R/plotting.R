plot.rjmcmc <- function (x, trace = TRUE, density = TRUE, smooth = FALSE, bwf,
	auto.layout = TRUE, ask = dev.interactive(), ...)
{
	.detachDiversitree()

    oldpar <- NULL
    on.exit(par(oldpar))
    if (auto.layout) {
        mfrow <- .set.mfrow(Nchains = nchain(x), Nparms = nvar(x),
						   nplots = trace + density)
        oldpar <- par(mfrow = mfrow)
    }
    for (i in 1:nvar(x)) {
        y <- mcmc(as.matrix(x)[, i, drop = FALSE], start(x),
				  end(x), thin(x))
        if (trace)
			if(all(y>0)) log="y" else log=""
			traceplot(y, smooth = smooth, log=log, ...)
        if (density) {
            if (missing(bwf))
			densplot(y, ...)
            else densplot(y, bwf = bwf, ...)
        }
        if (i == 1)
		oldpar <- c(oldpar, par(ask = ask))
    }
}

plot.rjmcmc.list<-
function (x, trace = TRUE, density = TRUE, smooth = TRUE, bwf,
auto.layout = TRUE, ask = dev.interactive(), ...)
{
    oldpar <- NULL
    on.exit(par(oldpar))
    if (auto.layout) {
        mfrow <- .set.mfrow(Nchains = nchain(x), Nparms = nvar(x),
						   nplots = trace + density)
        oldpar <- par(mfrow = mfrow)
    }
	nc=nchain(x)
    for (i in 1:nvar(x)) {
        if (trace)
			y=sapply(x, function(z) z[,i])
			if(all(y>0)) log="y" else log=""
			traceplot(x[, i, drop = FALSE], smooth = smooth, log=log, lwd=1/nc,
				  ...)
        if (density) {
            if (missing(bwf))
			densplot(x[, i, drop = FALSE], ...)
            else densplot(x[, i, drop = FALSE], bwf = bwf, ...)
        }
        if (i == 1)
		oldpar <- c(oldpar, par(ask = ask))
    }
}



.samples.plot=function(x, par=c("jumps","shifts"), ...){
	par=match.arg(par, c("shifts","jumps"))

	ff=list(...)
	if("burnin"%in%names(ff)){
		if(!.withinrange(ff$burnin, 0, 1)) stop("Supply 'burnin' as a fraction (between 0 and 1)")
	}
	if("level"%in%names(ff)){
		if(!.withinrange(ff$level, 0, 1)) stop("Supply 'level' as a fraction (between 0 and 1)")
	}
	switch(par,
		   jumps=.jumps.plot(x, ...),
		   shifts=.shifts.plot(x, ...)
		   )
}

plot.auteurMCMCMC=function(x, par=c("jumps","shifts"), ...){
	.samples.plot(x, par, ...)
}

plot.auteurMCMC=function(x, par=c("jumps","shifts"), ...){
	.samples.plot(x, par, ...)
}


#plotting function for comparing posterior densities of estimates
#author: JM EASTMAN 2010

.traceplot <-
function(obj, col, alpha, lwd=1, hpd=0.95, bars=TRUE, legend.control=list(plot=TRUE, pos=NA, cex=1, pt.cex=1, pch=22, title=""), truncate=list(min=NA, max=NA), xlim=list(min=NA, max=NA), ylim=NULL, ...){

	.infer.y <-
	function(xx, x, y) {
		f=lm(y~x)
		p=unname(coef(f))
		yy=xx*p[2]+p[1]
		return(yy)
	}

	.nearest.pair <-
	function(val, x) {
		dev=abs(x-val)
		names(dev)=1:length(dev)
		return(sort(as.numeric(names(sort(dev))[1:2])))
	}


	if(!is.data.frame(obj)) {
		if(is.vector(obj)) {
			obj=as.data.frame(obj)
		} else if(is.matrix(obj)) {
			obj=data.frame(obj, stringsAsFactors=FALSE)
		} else {
			stop("Object must be supplied as a vector, matrix, or dataframe.")
		}
	}

	# prepare (or assign) necessary arguments
	control=list(plot=TRUE,pos=NA,cex=1,pt.cex=1,pch=22,title="")
	control[names(legend.control)]<-legend.control
	legend.control=control
	if(missing(col)) col=gray.colors(ncol(obj))
	if(length(col)!=ncol(obj)) col=gray.colors(ncol(obj))
	if(missing(alpha)) alpha=0.8
	line.col=.transparency(col, min(c(0.95, alpha+alpha*0.25)))
	col=.transparency(col, alpha)

	# truncate data
	if(any(!is.na(truncate))) {
		trunc.todo=c("min","max")[match(names(!is.na(truncate)),c("min","max"))]
		for(t in 1:length(trunc.todo)) {
			if(length(truncate[[t]])!=ncol(obj) & length(truncate[[t]]==1)) {
				truncate[[t]]=rep(truncate[[t]],ncol(obj))
				warning(paste("assuming truncation value for ", sQuote(trunc.todo[t]), " is the same for all columns in 'obj'.", sep=""))
			}
			for(c in 1:ncol(obj)) {
				if(trunc.todo[t]=="min") {
					if(!is.na(truncate$min[c])) {
						mi=sapply(obj[,c], function(x) if(is.na(x) | x<truncate$min[c]) return(TRUE) else return(FALSE))
						obj[mi,c]=NA
					}
				} else if(trunc.todo[t]=="max") {
					if(!is.na(truncate$max[c])) {
						ma=sapply(obj[,c], function(x) if(is.na(x) | x>truncate$max[c]) return(TRUE) else return(FALSE))
						obj[ma,c]=NA
					}
				}
			}
		}
	}

	# prepare density objects and compute HDRs
	dd=lapply(1:ncol(obj), function(x) {y=range(obj[,x],na.rm=TRUE); density(obj[,x], from=min(y), to=max(y), na.rm=TRUE, n=1024, kernel="cosine")})
	if(!is.null(hpd)) {
		hh=lapply(1:ncol(obj), function(z) {zz=obj[!is.na(obj[,z]),z]; hdr(zz, hpd=hpd)})
	} else {
		hh=lapply(1:ncol(obj), function(z) {zz=obj[!is.na(obj[,z]),z]; c(min(zz), max(zz))})
	}
	xx=lapply(1:ncol(obj),function(z) {zz=obj[!is.na(obj[,z]),z]; return(c(min(zz[which(zz>=hh[[z]][1])]), max(zz[which(zz<=hh[[z]][2])])))})
	density.yy=lapply(1:length(dd), function(z) {sapply(xx[[z]], function(w) {xs=.nearest.pair(w, dd[[z]]$x); .infer.y(w, dd[[z]]$x[xs], dd[[z]]$y[xs])})})

	# include HDR points in density objects (necessary for precision in plotting)
	dd.new=lapply(1:length(dd), function(z) {
				  xxx=c(xx[[z]],dd[[z]]$x)
				  yyy=c(density.yy[[z]],dd[[z]]$y)
				  dx=sort(c(xx[[z]],dd[[z]]$x),index=TRUE)
				  return(dd.out=list(x=dx$x, y=yyy[dx$ix]))
				  })
	for(z in 1:length(dd)) {
		dd[[z]]$x=dd.new[[z]]$x
		dd[[z]]$y=dd.new[[z]]$y
	}

	ii=lapply(1:length(dd), function(z) sapply(dd[[z]]$x, function(w) .withinrange(w, min(hh[[z]]), max(hh[[z]]))))

	# find xlims for plotting
	lims.x.tmp=unlist(lapply(1:length(dd), function(z) return(c(max(dd[[z]]$x), min(dd[[z]]$x)))))
	lims.x=c(min(lims.x.tmp), max(lims.x.tmp))
	lims.x=c(min(lims.x)-0.1*min(lims.x), max(lims.x)+0.1*max(lims.x))
	if(any(!is.na(xlim))) {
		xlim.todo=c("min","max")[match(names(!is.na(xlim)),c("min","max"))]
		for(t in 1:length(xlim.todo)) {
			if(xlim.todo[t]=="min") lims.x[1]=xlim$min
			if(xlim.todo[t]=="max") lims.x[2]=xlim$max
		}
	}

	lims.y.tmp=unlist(lapply(1:length(dd), function(z) return(c(max(dd[[z]]$y), min(dd[[z]]$y)))))
	lims.y=c(-0.05, 1.05)*max(lims.y.tmp)
	if(is.null(ylim)) yylim=lims.y else yylim=ylim

	# PLOTTING of densities and legend
    plot(x = NULL, xlim = lims.x, ylim = yylim, ylab = "density", bty = "n", type = "n", ...)
	q=seq(0,min(lims.y)+0.4*min(lims.y),length=ncol(obj)+2)
	q=q[-c(1:2)]
	for(i in 1:length(dd)) {
		dat=dd[[i]]
		index=ii[[i]]
		xs=xx[[i]]
		ys=density.yy[[i]]
		polygon(c(dat$x[index], rev(dat$x[index])), c(dat$y[index],rep(0, length(dat$y[index]))), col=col[i], border=NA)
		dat$x=c(min(dat$x, na.rm=TRUE),dat$x,max(dat$x, na.rm=TRUE))
		dat$y=c(0,dat$y,0)
		lines(dat, col=line.col[i], lwd=lwd)
		if(bars) {
			arrows(xx[[i]][1], q[i], xx[[i]][2], q[i], code = 1, length = 0.0, col = line.col[i], lwd=lwd)
			points(data.frame(list(c(xx[[i]][1], xx[[i]][2]), c(q[i], q[i]))), pch=21, col=line.col[i], bg=line.col[i], cex=0.5*lwd)
		}
	}
	if(all(!is.null(names(obj))) & legend.control$plot) {
		if(is.na(legend.control$pos)) {
			mm=c(unlist(obj))
			mm=mm[!is.na(mm)]
			if(abs(mean(mm)-min(mm))>abs(mean(mm)-max(mm))) pos="topleft" else pos="topright"
		} else {
			pos=legend.control$pos
		}
		legend(x=pos, names(obj), pch=legend.control$pch, col="gray", pt.bg=col, bty="n", cex=legend.control$cex, title=legend.control$title, pt.cex=legend.control$pt.cex)
	}
}



#plotting function for comparing posterior densities of evolutionary process 'shifts' along a phylogeny
#author: JM EASTMAN 2010

.process.shifts<-function(phy, shifts, level) {
	if(!"hphylo"%in%class(phy)) stop("Supply 'phy' as an 'hphylo' object")
	hits.tmp<-apply(shifts,1,sum,na.rm=TRUE)
	hits=length(hits.tmp[hits.tmp>0])
	branches.tmp=apply(shifts, 2, function(x) sum(x, na.rm=TRUE))
	if(hits>0) branches.tmp=branches.tmp/hits
	branches=branches.tmp[branches.tmp>=level]
	branches=branches[order(branches, decreasing=TRUE)]
	if(!length(branches)) branches=NULL
	return(list(branch.shift.probability=branches))
}


.shifts.plot <-
function(samples, burnin=0, level=0.01, paint.branches=TRUE, colors=256, legend=TRUE, ...) {

  ## require colorspace
	phy=samples$phy

	color.length=17

	posterior.samples=list(rates=samples$rates, shifts=samples$shifts)

	tt=sapply(posterior.samples, function(x) "mcmc"%in%class(x))
	if(!all(tt)){
		stop("'samples' must contain an object of class 'rjmcmc' -- see to.auteur().")
	}

	phy=reorder(samples$phy)
	if(!"hphylo"%in%class(phy)) stop("Supply 'phy' as an 'hphylo' object")

	ps.tmp=lapply(posterior.samples, function(x) {mm=match(phy$hash[phy$edge[,2]], colnames(x)); return(x[,mm])})
	names(ps.tmp)=names(posterior.samples)
	posterior.samples=ps.tmp

	# collect data
	shifts=posterior.samples$shifts
	burnin=ceiling(burnin*nrow(shifts))
	if(burnin>0) shifts=shifts[-c(1:(burnin)),]
	shifts.res=.process.shifts(phy, shifts, level)
	ests=posterior.samples[[(which(names(posterior.samples)%in%c("rates"))->target)]]
	if(burnin>0) ests=ests[-c(1:(burnin)),]

	# determine whether to use logspace for plotting
	if(any(ests<=0,na.rm=TRUE)) logspace=FALSE else logspace=TRUE
	if(logspace) median.ests<-exp(apply(log(ests),2,median,na.rm=TRUE)) else median.ests<-apply(ests,2,median,na.rm=TRUE)

	param=names(posterior.samples)[target]

	# collect edge colors (for rates)
	if(paint.branches) {

		colors.branches.tmp=.branchcol.plot(phy, as.data.frame(ests), plot=FALSE, colors=list(branches=colors, legend=color.length, missing=1), log=logspace)
		colors.branches=colors.branches.tmp$col
	} else {
		colors.branches=1
	}

	# collect node colors (for shifts)
	ccx=colorspace::diverge_hcl(color.length, power = 0.5)
	c.seq=round(seq(-1,1,length=color.length),digits=1)
	all.nodes=seq(1:(Ntip(phy)+Nnode(phy)))[-(Ntip(phy)+1)]

	if(!is.null(shifts.res$branch.shift.probability)) {
		hh=phy$hash
		nodes=match(names(shifts.res$branch.shift.probability), phy$hash)
		NN=match(names(ests),phy$hash)
		xee=phy$edge[,2]
		shift.direction=sapply(all.nodes, function(x) {
					a=.get.ancestor.of.node(x, phy)
					if(hh[a]%in%colnames(ests)) comp=ests[,hh[a]] else comp=NULL
					this=ests[,hh[x]]
					if(!length(comp)) { # dealing with first descendant of root
						d=.get.desc.of.node(a, phy)
						d=d[which(d!=x)]
						d.shifts=shifts[, hh[d]]
						comp=ests[,hh[d]]
						x.res=sapply(1:length(d.shifts), function(y) {
								 if(is.na(d.shifts[y])) return(0)
								 if(d.shifts[y]==0) {
									zz=this[y]-comp[y]
									if(zz>0) {
										return(1)
									} else {
										if(zz<0) return(-1) else return(0)
									}
								 } else {
									return(0)
								 }
							}
						)
						x.res=mean(x.res[x.res!=0])
						if(is.na(x.res)) x.res=0
					} else {
						yy=this-comp
						zz=yy
						zz[yy>0]=1
						zz[yy<0]=-1
						zz[yy==0]=0
						x.res=mean(zz[zz!=0])
						if(is.na(x.res)) x.res=0
					}
					if(hh[x]%in%names(shifts.res$branch.shift.probability)) return(x.res) else return(0)
				}
		)
		colors.nodes=ccx[match(round(shift.direction,1),c.seq)]
		names(colors.nodes)=hh[all.nodes]
		colors.nodes=colors.nodes[match(phy$hash[phy$edge[,2]], names(colors.nodes) )]
	} else {
		colors.nodes=NULL
		shift.direction=rep(0,nrow(phy$edge))
	}

	## PLOTTING OF TREE ##
	if(legend) {
		def.par <- par(no.readonly = TRUE)
		on.exit(par(def.par))
		layout(matrix(c(1,2,1,3,1,4), 3, 2, byrow=TRUE), widths=c(20,5), respect=FALSE)
	}

	if(paint.branches){
		plot(phy, edge.color=colors.branches, no.margin=TRUE, ...)
	} else {
		plot(phy, no.margin=TRUE, ...)
	}

	NN=phy$hash[phy$edge[,2]]
	ll<-cc<-rr<-rep(0,length(NN))
	if(!is.null(shifts.res$branch.shift.probability)) {
		branches=names(shifts.res$branch.shift.probability)
		ll[match(branches, NN)]=1
		cc[match(branches, NN)]=shifts.res$branch.shift.probability
		rr=colors.nodes
	}

	edgelabels.auteur(text=NULL, pch=ifelse(ll==1, 21, NA), cex=4*cc, col=.transparency(rr,0.95), bg=.transparency(rr,0.5), lwd=0.5)
	## END PLOTTING of TREE ##

	if(legend) {
		legend.seq=seq(1,color.length,by=2)
		point.seq=rev(seq(0.2,0.8,length=length(legend.seq)))

		# shift direction
		plot(rep(-0.5, length(point.seq)), point.seq, xlim=c(-1,2), ylim=c(0,1), cex=2, pch=21, col = .transparency(rev(ccx[legend.seq]),0.95), bg = .transparency(rev(ccx[legend.seq]),0.5), bty="n", xaxt="n", yaxt="n")
		mtext("shift direction",side=3,line=-3,cex=0.75)
		text(rep(1, length(point.seq)), point.seq, adj=1, labels=sprintf("%9.2f", rev(c.seq[legend.seq])))

		# posterior estimates
		if(any(colors.branches!=1)) {
			cbt=colors.branches.tmp$legend.seq
			lchars=sapply(cbt, floor)
			if(all(lchars>0)) {
				ldec=1
			} else if(all(lchars==0)) {
				ldec=min(5, max(nchar(range(cbt))))
			} else {
				ldec=3
			}

			lnchar=max(nchar(lchars))+ldec
			plot(rep(-0.5, length(point.seq)), point.seq, xlim=c(-1,2), ylim=c(0,1), cex=2, pch=22, col = "darkgray", bg = rev(colors.branches.tmp$legend.col[legend.seq]), bty="n", xaxt="n", yaxt="n")
			mtext(paste("posterior ",param,sep=""),side=3,line=-3,cex=0.75)
			text(rep(1, length(point.seq)), point.seq, adj=1, labels=sprintf(paste("%",max(10,lnchar),".",ldec,"f",sep=""), cbt[legend.seq]))

		}

		# posterior probabilities of shift
		plot(rep(-0.5, length(point.seq)), point.seq, xlim=c(-1,2), ylim=c(0,1), cex=4*(seq(1, 0, length=9)), pch=21, col = .transparency("darkgray",0.8), bg = .transparency("white",0.8), bty="n", xaxt="n", yaxt="n")
		mtext("shift probability",side=3,line=-3,cex=0.75)
		text(rep(1, length(point.seq)), point.seq, adj=1, labels=sprintf("%10.3f", seq(1, 0, length=9)))

		# reset plotting device
		invisible()
	}

	# GENERATE TABULAR OUTPUT of RESULTS

	allres=data.frame(matrix(NA, nrow=nrow(phy$edge), ncol=3))
	shift.direction=shift.direction[match(phy$edge[,2],all.nodes)]
	shift.probability=cc
	allres=data.frame(branch=phy$edge[,2],shift.direction,shift.probability,median.ests)
	names(allres)[ncol(allres)]=paste("median",param,sep=".")
	rownames(allres)=phy$hash[phy$edge[,2]]
	return(res=allres)
}

.process.jumps<-function(phy, jumps) {
	if(!"hphylo"%in%class(phy)) stop("Supply 'phy' as an 'hphylo' object")
	m=match(phy$hash[phy$edge[,2]], colnames(jumps))
	jumps=jumps[,m]
	jump.prob <- apply(jumps,2,function(x) sum(x!=0,na.rm=TRUE)/length(x[!is.na(x)]))
	jump.counts <- apply(jumps,2,function(x) if(any(x!=0,na.rm=TRUE)) return(mean(x[x!=0],na.rm=TRUE)) else return(0))
	return(list(jump.probability=jump.prob, jump.counts=jump.counts))
}

.jumps.plot <-
function(samples, burnin=0, level=0.01, paint.branches=TRUE, colors=256, legend=TRUE, ...) {

  ## require colorspace
	phy=samples$phy

	if(!"hphylo"%in%class(phy)) stop("Supply 'phy' as an 'hphylo' object")
	if("mcmc"%in%class(samples$jumps)){
		posterior.samples=samples$jumps
	} else {
		stop("'samples' must contain an object of class 'rjmcmc' -- see load.rjmcmc().")
	}
	phy=reorder(phy)
	zz=match(phy$hash[phy$edge[,2]], colnames(posterior.samples))
	ps.tmp=posterior.samples[,zz]
	jumps=ps.tmp
	rates=samples$rates[,zz]

	# collect data

	jumpvar=samples$log[,"jumpvar"]
	if(burnin!=0){
		burnin=ceiling(burnin*nrow(jumps))
		jumps=jumps[-c(1:(burnin)),]
		rates=rates[-c(1:(burnin)),]
		jumpvar=jumpvar[-c(1:burnin)]
	}
	bm.var=rates
	j.var=jumps*c(jumpvar)


	jumps.res=.process.jumps(phy, jumps)

	j.prob=jumps.res$jump.probability
	j.count=jumps.res$jump.counts

	cols.x=pretty(c(1,max(5,max(j.count))),5)
	color.length=length(cols.x)
	rrx=colorspace::diverge_hcl(1+2*color.length, power = 1.5)
	rrx=rrx[(color.length+2):length(rrx)]

	NN=phy$edge[,2]
	hh=phy$hash[NN]
	ll<-cc<-rr<-rep(0,length(NN))
	if(any(j.prob>level)) {
		branches=names(j.prob[j.prob>level])
		ll[match(branches, hh)]=1
	}
	cc=j.prob
	cc[is.na(cc)]=0
	rr=sapply(j.count, function(x) {tmp=rrx[min(which(abs(x-cols.x)==min(abs(x-cols.x))))]; tmp})

	## PLOTTING OF TREE ##

	if(legend) {
		def.par <- par(no.readonly = TRUE)
		on.exit(par(def.par))
		layout(matrix(c(1,2,1,3), 2, 2, byrow=TRUE), widths=c(20,5), respect=FALSE)
	}

	if(paint.branches){
		# scalars=j.var+bm.var
		scalars=bm.var
		colors.branches=.branchcol.plot(phy, as.data.frame(scalars), plot=FALSE, colors=list(branches=colors, legend=17, missing=1))
		plot(phy, edge.color=colors.branches$col, no.margin=TRUE, ...)
	} else {
		plot(phy, no.margin=TRUE, ...)
	}


	edgelabels.auteur(text=NULL, pch=ifelse(ll==1, 21, NA), cex=4*cc, col=.transparency("red",0.95), bg=.transparency(rr,0.5), lwd=0.5)
	## END PLOTTING of TREE ##

	if(legend){
		legend.seq=seq(1,color.length,by=2)
		point.seq=rev(seq(0.2,0.8,length=color.length))
		prob.seq=rev(seq(0.2,0.8,length=9))

		# jumps count
		plot(rep(-0.5, length(point.seq)), point.seq, xlim=c(-1,2), ylim=c(0,1), cex=2, pch=21, col = .transparency("red",0.95), bg = .transparency(rev(rrx),0.5), bty="n", xaxt="n", yaxt="n")
		mtext("inferred jumps",side=3,line=-3,cex=0.75)
		text(rep(1, length(point.seq)), point.seq, adj=1, labels=sprintf("%i", rev(cols.x)))

		# jump probability
		plot(rep(-0.5, length(prob.seq)), prob.seq, xlim=c(-1,2), ylim=c(0,1), cex=4*(seq(1, 0, length=9)), pch=21, col = .transparency("darkgray",0.8), bg = .transparency("white",0.8), bty="n", xaxt="n", yaxt="n")
		mtext("jump probability",side=3,line=-3,cex=0.75)
		text(rep(1, length(prob.seq)), prob.seq, adj=1, labels=sprintf("%10.3f", seq(1, 0, length=9)))

		# reset plotting device
		invisible()
	}

	jumps.res=data.frame(jumps.res)
	jumps.res=cbind(branch=NN, jumps.res)

	return(jumps.res)
}

# general phylogenetic plotting utility, given a named vector or data.frame of values that can be associated with phy$edge[,2]
# author: JM EASTMAN 2010
# note: small values are given bluish hues, large values reddish hues; median values are given gray hues
# added some stuff; numeric vector input didn't seem to be supported.
.branchcol.plot <- function (phy, cur.rates, colors = list(branches = 256, legend = 17, missing = 1),
	digits = 3, plot = TRUE, legend = TRUE, legend.title = "", log = FALSE, ...) {


    ## require colorspace

	if (!"hphylo" %in% class(phy)) stop("Supply 'phy' as an 'hphylo' object");

	if ("data.frame" %in% class(cur.rates)) {
		if (!is.null(colnames(cur.rates))) {
			cur.rates <- cur.rates[,match(phy$hash[phy$edge[,2]], colnames(cur.rates))];
		} else {
			names(cur.rates) <- phy$hash[phy$edge[,2]];
			warning("Rates assumed to be ordered as in 'phy$edge'");
		}
	} else if ("numeric" %in% class(cur.rates)) {
		if (!is.null(names(cur.rates))) {
			cur.rates <- cur.rates[match(phy$hash[phy$edge[,2]], names(cur.rates))];
		} else {
			names(cur.rates) <- phy$hash[phy$edge[,2]];
			warning("Rates assumed to be ordered as in 'phy$edge'");
		}
		cur.rates <- as.data.frame(t(cur.rates));
	} else {
		stop("Expecting either a data.frame or named numeric vector");
	}

	cur.rates <- apply(cur.rates, 2, function(x) {
		if (any(!is.na(x))) {
			return(median(x, na.rm=TRUE));
		} else {
			return(NA);
		}
	})
	if (log) {
		ests <- log(cur.rates);
	} else {
		ests <- cur.rates;
	}

	ms <- median(ests, na.rm=TRUE);
	mm <- sapply(ests, function(x) x - ms);
	cce <- colorspace::diverge_hcl(2 * colors$branches + 1, power = 0.5);
	lcce <- cce[round(seq(1, length(cce), length=colors$legend))];
	e.seq <- seq(-max(abs(mm + 0.05 * ms), na.rm=TRUE), max(abs(mm + 0.05 * ms), na.rm=TRUE), length = 2 * colors$branches + 1);
	lseq <- e.seq+ms;
	lseq <- seq(min(lseq), max(lseq), length=colors$legend);
	lcce <- cce[round(seq(1, length(cce), length=colors$legend))];
	if (log) lseq <- exp(rev(lseq)) else lseq <- rev(lseq)

	ucr <- unique(cur.rates);
	ucr <- ucr[!is.na(ucr)];
	if (length(ucr) == 1) {
		mp <- cce[round(length(cce)/2)];
		colors.branches <- rep(mp, length(mm));
		colors.branches[is.na(mm)] <- colors$missing;
	} else {
		colors.branches <- sapply(mm, function(x) {
			if (is.na(x)) {
				return(colors$missing);
			} else {
				cce[which(min(abs(e.seq - x)) == abs(e.seq - x))];
			}
		})
	}

	if (plot) {
		plot.phylo(phy, cex=0.1, edge.color=colors.branches, ...)
		if (legend) {
			legend("topright", title=legend.title, cex=0.5, pt.cex=1, text.col="darkgray",
				   legend = sprintf(paste("%", 2*digits, paste(digits, "f", sep=""), sep="."), lseq),
				   pch=21, ncol=1, col = "darkgray", pt.bg = rev(lcce), box.lty="blank", border="white");
		}
	} else {
		return(list(col=colors.branches,legend.seq=lseq,legend.col=lcce));
	}
}




#general phylogenetic plotting utility, which is a modification of ape::edgelabels, plotting edge symbols at the temporal beginning of the branch
#author: E PARADIS 2009 and JM EASTMAN 2010
#note: may not be trustworthy where lastPP$type is not "phylogram"

edgelabels.auteur <-
function (text, edge, adj = c(0.5, 0.5), frame = "rect", pch = NULL,
thermo = NULL, pie = NULL, piecol = NULL, col = "black",
bg = "lightgreen", horiz = FALSE, width = NULL, height = NULL,
...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(edge)) {
        sel <- 1:dim(lastPP$edge)[1]
        subedge <- lastPP$edge
    }
    else {
        sel <- edge
        subedge <- lastPP$edge[sel, , drop = FALSE]
    }
    if (lastPP$type == "phylogram") {
        if (lastPP$direction %in% c("rightwards", "leftwards")) {
            XX <- (lastPP$xx[subedge[, 1]])
            YY <- lastPP$yy[subedge[, 2]]
        }
        else {
            XX <- lastPP$xx[subedge[, 2]]
            YY <- (lastPP$yy[subedge[, 1]])
        }
    }
    else {
        XX <- (lastPP$xx[subedge[, 2]])
        YY <- (lastPP$yy[subedge[, 2]])
    }
	if(missing(text)) text=lastPP$edge[,2]

    BOTHlabels(text, sel, XX, YY, adj, frame, pch, thermo, pie,
					 piecol, col, bg, horiz, width, height, ...)
}





.acegram=function(phy, dat, cex.node=2, cex.tip=2, labs=TRUE, ...){
	root=Ntip(phy)+1
	dd=dat
	names(dd)=match(names(dat),phy$tip.label)
	xx=c(ace(dat, phy, CI=FALSE, method="pic")$ace,dd)
	alpha=xx[names(xx)==root]
	xx=xx[names(xx)!=root]
	mm=match(phy$edge[,2],names(xx))
	hist=data.frame(phy$edge, phy$edge.length, xx[mm])
	names(hist)=c("ancestor","descendant","edge","phenotype")
	hist$time=sapply(hist$descendant, function(x) {anc=c(x,.get.ancestors.of.node(x,phy)); anc=anc[anc!=root]; sum(phy$edge.length[match(anc, phy$edge[,2])])})
	mm=max(abs(alpha-unlist(hist$phenotype)))
	root=Ntip(phy)+1
	pp=pretty(c(0,max(hist$time)))
	plot(x=NULL, y=NULL, xlim=rev(range(pp)), ylim=range(pretty(c(-mm+alpha, mm+alpha))), bty="n", xlab="time", ylab="phenotypic value")
	hist$ptime=abs(hist$time-max(hist$time))
	mbt=max(branching.times(phy))
	for(i in 1:nrow(hist)) {
		start=ifelse(hist$ancestor[i]==root, alpha, hist$phenotype[which(hist$descendant==hist$ancestor[i])])
		stime=ifelse(hist$ancestor[i]==root, mbt, hist$ptime[which(hist$descendant==hist$ancestor[i])])

		end=hist$phenotype[i]
		etime=hist$ptime[i]
		lines(c(stime,etime),c(start,end),col=.transparency("gray25",0.75),...)
	}
	nn=hist$descendant<=Ntip(phy)
	if(labs) {
		ll=phy$tip.label[hist$descendant[nn]]
		tt=hist$ptime[nn]
		yy=hist$phenotype[nn]
		text(tt+0.01*max(tt),yy,ll, cex=cex.tip, pos=4)
		points(hist$ptime[!nn],hist$phenotype[!nn],bg=.transparency("white",0.75),pch=21,cex=cex.node)
	} else {
		points(hist$ptime,hist$phenotype,bg=.transparency("white",0.75),pch=21,cex=ifelse(hist$descendant<=Ntip(phy), cex.tip, cex.node))
	}
	points(0,alpha,bg=.transparency("white",0.75),pch=21,cex=cex.node)
}


# Plotting utility from coda
# Author: Martyn Plummer
.set.mfrow <- function (Nchains = 1, Nparms = 1, nplots = 1, sepplot = FALSE)
{
  mfrow <- if (sepplot && Nchains > 1 && nplots == 1) {
    if (Nchains == 2) {
      switch(min(Nparms, 5), c(1, 2), c(2, 2), c(3, 2),
             c(4, 2), c(3, 2))
    }
    else if (Nchains == 3) {
      switch(min(Nparms, 5), c(2, 2), c(2, 3), c(3, 3),
             c(2, 3), c(3, 3))
    }
    else if (Nchains == 4) {
      if (Nparms == 1)
        c(2, 2)
      else c(4, 2)
    }
    else if (any(Nchains == c(5, 6, 10, 11, 12)))
      c(3, 2)
    else if (any(Nchains == c(7, 8, 9)) || Nchains >= 13)
      c(3, 3)
  }
  else {
    if (nplots == 1) {
      mfrow <- switch(min(Nparms, 13), c(1, 1), c(1, 2),
                      c(2, 2), c(2, 2), c(3, 2), c(3, 2), c(3, 3),
                      c(3, 3), c(3, 3), c(3, 2), c(3, 2), c(3, 2),
                      c(3, 3))
    }
    else {
      mfrow <- switch(min(Nparms, 13), c(1, 2), c(2, 2),
                      c(3, 2), c(4, 2), c(3, 2), c(3, 2), c(4, 2),
                      c(4, 2), c(4, 2), c(3, 2), c(3, 2), c(3, 2),
                      c(4, 2))
    }
  }
  return(mfrow)
}
