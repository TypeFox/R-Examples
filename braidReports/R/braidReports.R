potentiationPlot <- function(bfit,levs,dnum=1,concs=NULL,control=list()) {
	if (is.null(control$xlog)) { control$xlog <- TRUE }
	if (is.null(control$ylog)) { control$ylog <- FALSE }
	if (is.null(control$addscales)) { control$addscales <- TRUE }
	if (is.null(control$ciplot)) { control$ciplot <- TRUE }
	else if (control$ciplot && (is.null(bfit$ciPass) || !bfit$ciPass)) {
		warning("Confidence intervals cannot be plotted if 'bfit' does not have valid bootstrapped coefficients")
	}
	if (is.null(bfit$ciPass) || !bfit$ciPass) { control$ciplot <- FALSE }
	if (is.null(control$cialph)) { control$cialph <- TRUE }
	if (is.null(control$pts)) { control$pts <- FALSE }
	if (is.null(concs)) {
		cpos <- bfit$conc1>0 & bfit$conc2>0
		if (dnum==1) {
			mnconc <- min(bfit$conc1[cpos])
			mxconc <- max(bfit$conc1[cpos])
		} else {
			mnconc <- min(bfit$conc2[cpos])
			mxconc <- max(bfit$conc2[cpos])
		}
		if (!is.null(control$pconcs)) {
			mnconc <- min(mnconc,min(control$pconcs))
			mxconc <- max(mxconc,max(control$pconcs))
		}
		if (control$xlog) { concs <- exp(seq(log(mnconc),log(mxconc),length=100)) }
		else { concs <- seq(mnconc,mxconc,length=100) }
	}
	concv <- rep(concs,times=length(levs))
	levv <- rep(levs,each=length(concs))
	if (control$ciplot) {
		if (dnum==1) { actmat <- calcBRAIDconfint(bfit,function(parv) evalBRAIDrsm(concv,levv,parv)) }
		else { actmat <- calcBRAIDconfint(bfit,function(parv) evalBRAIDrsm(levv,concv,parv)) }
		if (!is.null(control$ytrans)) { actmat <- control$ytrans(actmat) }
		pldf <- data.frame(conc=concv,level=factor(levv),act=actmat[,2],actmin=actmat[,1],actmax=actmat[,3])
	} else {
		if (dnum==1) { pldf <- data.frame(conc=concv,level=factor(levv),act=evalBRAIDrsm(concv,levv,bfit$fullpar)) }
		else { pldf <- data.frame(conc=concv,level=factor(levv),act=evalBRAIDrsm(levv,concv,bfit$fullpar)) }
		if (!is.null(control$ytrans)) { pldf$act <- control$ytrans(pldf$act) }
	}
	if (!is.null(control$xtrans)) { pldf$conc <- control$xtrans(pldf$conc) }
	pplt <- ggplot(pldf,aes_string(x="conc",y="act"))
	if (control$ciplot) {
		if (control$cialph) { pplt <- pplt+geom_ribbon(aes_string(ymin="actmin",ymax="actmax",fill="level"),alpha=0.3) }
		else { pplt <- pplt+geom_ribbon(aes_string(ymin="actmin",ymax="actmax",fill="level")) }
	}
	if (!is.null(control$lnsize)) { pplt <- pplt+geom_line(aes_string(colour="level"),size=control$lnsize) }
	else { pplt <- pplt+geom_line(aes_string(colour="level")) }
	if (control$pts) {
		if (is.null(control$pconcs)) {
			cpos <- bfit$conc1>0 & bfit$conc2>0
			if (dnum==1) { ptconcs <- unique(bfit$conc1[cpos]) }
			else { ptconcs <- unique(bfit$conc2[cpos]) }
			ptconcs <- ptconcs[ptconcs>=min(concs) & ptconcs<=max(concs)]
		} else { ptconcs <- control$pconcs }
		ptconcv <- rep(ptconcs,times=length(levs))
		ptlevv <- rep(levs,each=length(ptconcs))
		if (dnum==1) { ptsdf <- data.frame(conc=ptconcv,level=factor(ptlevv),act=evalBRAIDrsm(ptconcv,ptlevv,bfit$fullpar)) }
		else { ptsdf <- data.frame(conc=ptconcv,level=factor(ptlevv),act=evalBRAIDrsm(ptlevv,ptconcv,bfit$fullpar)) }
		if (!is.null(control$xtrans)) { ptsdf$conc <- control$xtrans(ptsdf$conc) }
		if (!is.null(control$ytrans)) { ptsdf$act <- control$ytrans(ptsdf$act) }
		if (!is.null(control$ptsize)) { pplt <- pplt+geom_point(aes_string(colour="level"),data=ptsdf,size=control$ptsize) }
		else { pplt <- pplt+geom_point(aes_string(colour="level"),data=ptsdf) }
	}
	if (!is.null(control$thresh)) {
		if (!is.null(control$ytrans)) { control$thresh <- control$ytrans(control$thresh) }
		if (!is.null(control$lnsize)) { pplt <- pplt+geom_hline(yintercept=control$thresh,linetype=2,size=control$lnsize) }
		else { pplt <- pplt+geom_hline(yintercept=control$thresh,linetype=2) }
	}
	if (control$addscales) {
		if (control$xlog) { pplt <- pplt+scale_x_log10() }
		if (control$ylog) { pplt <- pplt+scale_y_log10() }
	}
	return(pplt)
}

calculateIAE <- function(parv,lev,macs,same=FALSE) {
	nconcs <- 100
	if (length(macs)==1) { macs <- c(macs,macs) }
	if (same && macs[1]!=macs[2]) { stop("If 'same' is TRUE, maximum achievable concentrations must be equal.") }
	isIncr <- sign(parv[10]-parv[7])
	iaev <- rep(0,length(lev))
	for (i in 1:length(lev)) {
		E <- lev[i]
		if (sign(E-parv[7])==0 || sign(E-parv[7])!=sign(parv[10]-parv[7])) { iaev[i] <- Inf; next }
		if (sign(parv[10]-E)==0 || sign(parv[10]-E)!=sign(parv[10]-parv[7])) { iaev[i] <- 1; next }
		if (parv[6]<0) {
			eps <- (1-(parv[6]^2/4))^(parv[5]*sqrt(parv[3]*parv[4]))
			bs1 <- eps*(parv[8]-E)/(E-parv[7]) - (1-eps)*(parv[10]-parv[8])/(parv[10]-parv[7])
			if (bs1<=0) { ctop1 <- Inf }
			else { ctop1 <- parv[1]*(1/bs1)^(1/parv[3]) }
			bs2 <- eps*(parv[9]-E)/(E-parv[7]) - (1-eps)*(parv[10]-parv[9])/(parv[10]-parv[7])
			if (bs2<=0) { ctop2 <- Inf }
			else { ctop2 <- parv[2]*(1/bs2)^(1/parv[4]) }
		} else {
			if ((E-parv[7])/(parv[8]-parv[7]) >= 1) { ctop1 <- Inf }
			else { ctop1 <- parv[1]*((E-parv[7])/(parv[8]-E))^(1/parv[3]) }
			if ((E-parv[7])/(parv[9]-parv[7]) >= 1) { ctop2 <- Inf }
			else { ctop2 <- parv[2]*((E-parv[7])/(parv[9]-E))^(1/parv[4]) }
		}
		ctop1 <- min(ctop1,macs[1])
		ctop2 <- min(ctop2,macs[2])
		uc1 <- ctop1*seq(1,2*nconcs-1,by=2)/(2*nconcs)
		uc2 <- ctop2*seq(1,2*nconcs-1,by=2)/(2*nconcs)
		sc1 <- rep(uc1,each=length(uc2))
		sc2 <- rep(uc2,times=length(uc1))
		act <- evalBRAIDrsm(sc1,sc2,parv)
		if (same) { rel <- sum((isIncr*act)<(isIncr*E) & (sc1+sc2)<macs[1]) }
		else { rel <- sum((isIncr*act)<(isIncr*E)) }
		rae <- ((ctop1*ctop2)/(macs[1]*macs[2]))*(rel/length(sc1))
		if (same) { rae <- 2*rae }
		iaev[i] <- sqrt(1/rae)
	}
	return(iaev)
}

makeBRAIDreport <- function(brdAnalysis,compounds,iaelevs,macs,control=list()) {
	# Process input parameters and perform sanity checks
	if (is.null(control$actlabel)) { control$actlabel <- "Effect" }
	if (is.null(control$abbs)) {
		control$abbs <- c(gsub("[^[:alnum:]]","",substr(compounds[1],1,3)),gsub("[^[:alnum:]]","",substr(compounds[2],1,3)))
	}
	if (is.null(control$irreg)) { control$irreg <- FALSE }
	if (is.null(control$clog)) { control$clog <- TRUE }
	if (is.null(control$zlog)) { control$zlog <- FALSE }
	if (is.null(control$palette) || control$palette=="gjet" || control$palette=="revgjet") {
		gjet <- c("#007F00", "green", "#7FFF00", "yellow", "#FF7F00", "red", "#7F0000")
		if (!is.null(control$palette) && control$palette=="revgjet") { control$palette <- rev(gjet) }
		else  { control$palette <- gjet }
	}
	if (is.null(control$plot) || control$plot) {
		if (requireNamespace("gridExtra",quietly=TRUE) && packageVersion("gridExtra")>=package_version("2.0.0") &&
				requireNamespace("gtable",quietly=TRUE)) {
			if (is.null(control$plot)) { control$plot <- TRUE }
		} else {
			warning("Plotting using 'makeBRAIDreport' requires the package 'gridExtra' (v2.0.0 or greater).")
			control$plot <- FALSE
		}
	}
	if (is.null(control$return)) { control$return <- FALSE }
	if (length(macs)==1) { macs <- c(macs,macs) }
	if (compounds[1]==compounds[2] && macs[1]!=macs[2]) { stop("If compounds are the same, maximum achievable concentrations must be equal.") }
	
	# Construct data frame for plotting
	syn_f <- data.frame(conc1=brdAnalysis$concs[,1],conc2=brdAnalysis$concs[,2],act=brdAnalysis$act,
					predc=fitted(brdAnalysis$braidFit))
	# Derive additive surface from single-agent fits if available, from RSM fit if not
	if (!is.null(brdAnalysis$hfit1)) { h1par <- coef(brdAnalysis$hfit1$allfits[[brdAnalysis$hfit1$bestModIdx]]) }
	else { h1par <- NULL }
	if (!is.null(brdAnalysis$hfit2)) { h2par <- coef(brdAnalysis$hfit2$allfits[[brdAnalysis$hfit2$bestModIdx]]) }
	else { h2par <- NULL }
	addpar <- getAdditiveParameters(brdAnalysis$braidFit$fullpar,h1par,h2par)
	
	# Calculate interpolated activities for "Interpolation error" map
	syn_f$intact <- syn_f$act
	if (!control$irreg) {
		posind <- syn_f$conc1>0 & syn_f$conc2>0
		margx <- syn_f$conc1>0 & syn_f$conc2==0
		margy <- syn_f$conc1==0 & syn_f$conc2>0
		if (control$clog) {
			syn_f$intact[posind] <- gaussInterp2d(log(syn_f$conc1[posind]),log(syn_f$conc2[posind]),
										syn_f$act[posind],log(syn_f$conc1[posind]),log(syn_f$conc2[posind]))
			if (any(margx)) { syn_f$intact[margx] <- gaussInterp1d(log(syn_f$conc1[margx]),syn_f$act[margx],log(syn_f$conc1[margx])) }
			if (any(margy)) { syn_f$intact[margy] <- gaussInterp1d(log(syn_f$conc2[margy]),syn_f$act[margy],log(syn_f$conc2[margy])) }
		} else {
			syn_f$intact[posind] <- gaussInterp2d(syn_f$conc1[posind],syn_f$conc2[posind],
										syn_f$act[posind],syn_f$conc1[posind],syn_f$conc2[posind])
			if (any(margx)) { syn_f$intact[margx] <- gaussInterp1d(syn_f$conc1[margx],syn_f$act[margx],syn_f$conc1[margx]) }
			if (any(margy)) { syn_f$intact[margy] <- gaussInterp1d(syn_f$conc2[margy],syn_f$act[margy],syn_f$conc2[margy]) }
		}
	}
	syn_f$pred <- evalBRAIDrsm(syn_f$conc1,syn_f$conc2,brdAnalysis$braidFit$fullpar)
	syn_f$lpred <- evalBRAIDrsm(syn_f$conc1,syn_f$conc2,addpar)
	syn_f$lpredc <- evalBRAIDrsm(brdAnalysis$braidFit$conc1,brdAnalysis$braidFit$conc2,addpar)
	
	# Transform x- and y-coordinates if desired (e.g. molar to micromolar)
	pmacs <- macs
	if (!is.null(control$xtrans)) {
		syn_f$conc1 <- control$xtrans(syn_f$conc1)
		pmacs[1] <- control$xtrans(macs[1])
	}
	if (!is.null(control$ytrans)) {
		syn_f$conc2 <- control$ytrans(syn_f$conc2)
		pmacs[2] <- control$ytrans(macs[2])
	}
	
	# Build activity and error map plots
	plts <- list()
	plts[[1]] <- responseMap(act~conc1+conc2,syn_f,interpolate=FALSE,margins=TRUE,irreg=control$irreg,logscale=control$clog)
	plts[[2]] <- responseMap(act~conc1+conc2,syn_f,interpolate=TRUE,margins=TRUE,irreg=control$irreg,logscale=control$clog)
	plts[[3]] <- responseMap(lpred~conc1+conc2,syn_f,interpolate=TRUE,margins=TRUE,irreg=control$irreg,logscale=control$clog)
	plts[[4]] <- responseMap(pred~conc1+conc2,syn_f,interpolate=TRUE,margins=TRUE,irreg=control$irreg,logscale=control$clog)
	plts[[5]] <- responseMap((intact-act)~conc1+conc2,syn_f,interpolate=FALSE,margins=TRUE,irreg=control$irreg,logscale=control$clog)
	plts[[6]] <- responseMap((lpredc-act)~conc1+conc2,syn_f,interpolate=FALSE,margins=TRUE,irreg=control$irreg,logscale=control$clog)
	plts[[7]] <- responseMap((predc-act)~conc1+conc2,syn_f,interpolate=FALSE,margins=TRUE,irreg=control$irreg,logscale=control$clog)
	
	# Format plots with matching color-scales, appropriate axis labels and titles
	if (is.null(control$xunit)) {
		xlab <- bquote(Conc[.(control$abbs[1])])
		xlab2 <- control$abbs[1]
	} else {
		if (class(control$xunit)=="expression")  {
			xlab <- bquote(Conc[.(control$abbs[1])]~(.(control$xunit[[1]])))
			xlab2 <- bquote(.(control$abbs[1])~(.(control$xunit[[1]])))
		} else  {
			xlab <- bquote(Conc[.(control$abbs[1])]~(.(control$xunit)))
			xlab2 <- paste(control$abbs[1]," (",control$xunit,")",sep="")
		}
	}
	if (is.null(control$yunit)) {
		ylab <- bquote(Conc[.(control$abbs[2])])
		ylab2 <- control$abbs[2]
	}
	else {
		if (class(control$yunit)=="expression")  {
			ylab <- bquote(Conc[.(control$abbs[2])]~(.(control$yunit[[1]])))
			ylab2 <- bquote(.(control$abbs[2])~(.(control$yunit[[1]])))
		} else  {
			ylab <- bquote(Conc[.(control$abbs[2])]~(.(control$yunit)))
			ylab2 <- paste(control$abbs[2]," (",control$yunit,")",sep="")
		}
	}
	totmin <- Inf
	totmax <- -Inf
	for (i in 1:4) {
		totmin <- min(totmin,min(plts[[i]]$data$z))
		totmax <- max(totmax,max(plts[[i]]$data$z))
		for (j in 1:length(plts[[i]]$layers)) {
			if (!is.null(plts[[i]]$layers[[j]]$data$z)) {
				totmin <- min(totmin,min(plts[[i]]$layers[[j]]$data$z))
				totmax <- max(totmax,max(plts[[i]]$layers[[j]]$data$z))
			}
		}
	}
	formatThisResponseMap <- function(plt,interp=TRUE,xl=NULL,yl=NULL,ttl=NULL) {
		nplt <- plt
		if (is.null(control$ztrans)) { nplt <- plt+scale_fill_gradientn(control$actlabel,colours=control$palette,limits=c(totmin,totmax)) }
		else { nplt <- plt+scale_fill_gradientn(control$actlabel,colours=control$palette,limits=c(totmin,totmax),
						labels=function(z) signif(control$ztrans(z),2)) }
		if (interp) {
			if (iaelevs[1]>min(plt$data$z) && iaelevs[1]<max(plt$data$z)) { 
				if (is.null(control$lnsize)) { nplt <- nplt+geom_contour(aes_string(z="z"),colour="black",breaks=iaelevs[1]) }
				else { nplt <- nplt+geom_contour(aes_string(z="z"),colour="black",breaks=iaelevs[1],size=control$lnsize) }
			}
			if (length(iaelevs)>1 && iaelevs[2]>min(plt$data$z) && iaelevs[2]<max(plt$data$z)) { 
				if (is.null(control$lnsize)) { nplt <- nplt+geom_contour(aes_string(z="z"),colour="white",breaks=iaelevs[2]) }
				else { nplt <- nplt+geom_contour(aes_string(z="z"),colour="white",breaks=iaelevs[2],size=control$lnsize) }
			}
			if (compounds[1]==compounds[2]) {
				if (min(plt$data$x+plt$data$y)>pmacs[1] || max(plt$data$x+plt$data$y)<pmacs[1]) { tplt <- NULL }
				else {
					tplt <- data.frame(x=c(pmacs[1],seq(pmacs[1]-min(plt$data$y),min(plt$data$x),length=100),0))
					tplt$y <- pmacs[1]-tplt$x
					tplt$z <- 0
				}
			} else {
				if (pmacs[1]<min(plt$data$x) | pmacs[2]<min(plt$data$y)) { tplt <- NULL }
				else if (pmacs[1]>max(plt$data$x)) {
					if (pmacs[2]>max(plt$data$y)) { tplt <- NULL }
					else { tplt <- data.frame(x=c(0,Inf),y=pmacs[2],z=0) }
				} else if (pmacs[2]>max(plt$data$y)) { tplt <- data.frame(x=pmacs[1],y=c(0,Inf),z=0) }
				else { tplt <- data.frame(x=c(pmacs[1],pmacs[1],0),y=c(0,pmacs[2],pmacs[2]),z=0) }
			}
			if (!is.null(tplt)) {
				if (is.null(control$lnsize)) { nplt <- nplt+geom_path(data=tplt,colour="black",linetype=2) }
				else { nplt <- nplt+geom_path(data=tplt,colour="black",linetype=2,size=control$lnsize) }
			}
		}
		nplt <- nplt+labs(x=xl,y=yl,title=ttl)
		return(nplt)
	}
	plts[[1]] <- formatThisResponseMap(plts[[1]],interp=FALSE,xl=xlab,yl=ylab,ttl="Observed Effect (Raw)")+coord_fixed()
	plts[[2]] <- formatThisResponseMap(plts[[2]],interp=TRUE,xl=xlab,yl=ylab,ttl="Observed Effect")+coord_fixed()
	plts[[3]] <- formatThisResponseMap(plts[[3]],interp=TRUE,xl=xlab,yl=ylab,ttl="Additive Effect")+coord_fixed()
	plts[[4]] <- formatThisResponseMap(plts[[4]],interp=TRUE,xl=xlab,yl=ylab,ttl="BRAID Effect")+coord_fixed()
	errmin <- Inf
	errmax <- -Inf
	for (i in 5:7) {
		errmin <- min(errmin,min(plts[[i]]$data$z))
		errmax <- max(errmax,max(plts[[i]]$data$z))
		for (j in 1:length(plts[[i]]$layers)) {
			if (!is.null(plts[[i]]$layers[[j]]$data$z)) {
				errmin <- min(errmin,min(plts[[i]]$layers[[j]]$data$z))
				errmax <- max(errmax,max(plts[[i]]$layers[[j]]$data$z))
			}
		}
	}
	plts[[5]] <- plts[[5]]+scale_fill_gradient2("Error",limits=c(errmin,errmax))+coord_fixed()
	plts[[5]] <- plts[[5]]+labs(x=xlab,y=ylab,title="Interpolation Error")
	plts[[6]] <- plts[[6]]+scale_fill_gradient2("Error",limits=c(errmin,errmax))+coord_fixed()
	plts[[6]] <- plts[[6]]+labs(x=xlab,y=ylab,title=bquote(Additive~Error~(R^2==.(round(cor(syn_f$act,syn_f$lpredc),3)))))
	plts[[7]] <- plts[[7]]+scale_fill_gradient2("Error",limits=c(errmin,errmax))+coord_fixed()
	plts[[7]] <- plts[[7]]+labs(x=xlab,y=ylab,title=bquote(BRAID~Error~(R^2==.(round(cor(syn_f$act,syn_f$predc),3)))))
	
	# Make Potentiation Plots
	if (control$clog) {
		pxlevs <- c(0,10^(ceiling(log10(min(brdAnalysis$concs[brdAnalysis$concs[,1]>0 & brdAnalysis$concs[,2]>0,1])))+(0:2)))
		pylevs <- c(0,10^(ceiling(log10(min(brdAnalysis$concs[brdAnalysis$concs[,1]>0 & brdAnalysis$concs[,2]>0,2])))+(0:2)))
	} else {
		pxmax <- max(brdAnalysis$concs[brdAnalysis$concs[,1]>0 & brdAnalysis$concs[,2]>0,1])
		pxscl <- 10^floor(log10(pxmax))
		pxnum <- floor(pxmax/(3*pxscl))
		pxlevs <- (0:3)*pxnum*pxscl
		pymax <- max(brdAnalysis$concs[brdAnalysis$concs[,1]>0 & brdAnalysis$concs[,2]>0,2])
		pyscl <- 10^floor(log10(pymax))
		pynum <- floor(pymax/(3*pyscl))
		pylevs <- (0:3)*pynum*pyscl
	}
	pxcontrol <- list(xlog=control$clog,ylog=control$zlog)
	if (!is.null(control$ztrans)) { pxcontrol$ytrans <- control$ztrans }
	pxcontrol$thresh <- iaelevs[1]
	if (!is.null(control$lnsize)) { pxcontrol$lnsize <- control$lnsize }
	if (!control$irreg) { pxcontrol$pts <- TRUE }
	else { pxcontrol$pts <- FALSE }
	if (!is.null(control$ptsize)) { pxcontrol$ptsize <- control$ptsize }
	pycontrol <- pxcontrol
	if (!is.null(control$xtrans)) { pxcontrol$xtrans <- control$xtrans }
	if (!is.null(control$ytrans)) { pycontrol$xtrans <- control$ytrans }
	if (pxcontrol$pts && brdAnalysis$corrconc) {
		posind <- brdAnalysis$concs[,1]>0 & brdAnalysis$concs[,2]>0
		pxcontrol$pconcs <- sort(unique(brdAnalysis$concs[posind,1]))
		pycontrol$pconcs <- sort(unique(brdAnalysis$concs[posind,2]))
	}
	plts[[8]] <- potentiationPlot(brdAnalysis$braidFit,pylevs,dnum=1,control=pxcontrol)
	plts[[8]] <- plts[[8]]+labs(x=xlab,y=control$actlabel,colour=ylab2,fill=ylab2,title=paste("Potentiation of",compounds[1]))
	if (!is.null(control$ytrans)) {
		plts[[8]] <- plts[[8]]+scale_colour_discrete(label=function(s) as.character(control$ytrans(as.numeric(s))))+
								scale_fill_discrete(label=function(s) as.character(control$ytrans(as.numeric(s))))
	}
	plts[[9]] <- potentiationPlot(brdAnalysis$braidFit,pxlevs,dnum=2,control=pycontrol)
	plts[[9]] <- plts[[9]]+labs(x=ylab,y=control$actlabel,colour=xlab2,fill=xlab2,title=paste("Potentiation of",compounds[2]))
	if (!is.null(control$xtrans)) {
		plts[[9]] <- plts[[9]]+scale_colour_discrete(label=function(s) as.character(control$xtrans(as.numeric(s))))+
								scale_fill_discrete(label=function(s) as.character(control$xtrans(as.numeric(s))))
	}
	if (control$clog || control$zlog) {
		thm <- theme_get()
		tkcol <- thm$axis.ticks$colour
		tklen <- thm$axis.ticks.length
	}
	for (i in 1:9) {
		plts[[i]] <- plts[[i]]+theme(plot.margin=unit(c(0,0,0,0),"npc"),legend.margin=unit(-0.01,"npc"),
														legend.key.width=unit(0.01,"npc"),legend.key.height=unit(0.0125,"npc"),
														text=element_text(size=7),line=element_line(size=0.5))
		if (i<8) {
			if (control$clog) {
				if (is.null(control$lnsize)) {
					plts[[i]] <- plts[[i]]+annotation_logticks(sides="bl",short=0.333*tklen,mid=0.666*tklen,
													long=tklen,colour=tkcol)
				} else {
					plts[[i]] <- plts[[i]]+annotation_logticks(sides="bl",short=0.333*tklen,mid=0.666*tklen,
													long=tklen,colour=tkcol,size=control$lnsize)
				}
			}
		} else {
			if (control$clog) {
				if (control$zlog) { csides <- "bl" }
				else { csides <- "b" }
			} else if (control$zlog) { csides <- "l" }
			else { next }
			if (is.null(control$lnsize)) {
				plts[[i]] <- plts[[i]]+annotation_logticks(sides=csides,short=0.333*tklen,mid=0.666*tklen,
												long=tklen,colour=tkcol)
			} else {
				plts[[i]] <- plts[[i]]+annotation_logticks(sides=csides,short=0.333*tklen,mid=0.666*tklen,
												long=tklen,colour=tkcol,size=control$lnsize)
			}
		}
	}
	
	bsumm <- coef(summary(brdAnalysis$braidFit))
	fpar <- brdAnalysis$braidFit$fullpar
	if (!is.null(control$xtrans)) {
		bsumm[1,] <- control$xtrans(bsumm[1,])
		fpar[1] <- control$xtrans(fpar[1])
	}
	if (!is.null(control$ytrans)) {
		bsumm[2,] <- control$ytrans(bsumm[2,])
		fpar[2] <- control$ytrans(fpar[2])
	}
	if (!is.null(control$ztrans)) {
		erows <- intersect(c("E0","EfA","EfB","EfAB"),row.names(bsumm))
		bsumm[erows,] <- control$ztrans(bsumm[erows,])
		fpar[7:10] <- control$ztrans(fpar[7:10])
	}
	ptab1 <- makeParamTable(bsumm,fpar,parse=FALSE)
	ptab2 <- makeParamTable(bsumm,fpar,parse=TRUE)
	
	itab1 <- getAdditiveIndexTable(brdAnalysis,iaelevs,macs,same=(compounds[1]==compounds[2]))
	itab2 <- calcBRAIDconfint(brdAnalysis$braidFit,function(parv) calculateIAE(parv,iaelevs,macs,same=(compounds[1]==compounds[2])))
	if (!is.null(control$levtext)) { irows <- paste("IAE[",control$levtext,"]",sep="") }
	else if (!is.null(control$ztrans)) { irows <- paste("IAE[",signif(control$ztrans(iaelevs),2),"]",sep="") }
	else { irows <- paste("IAE[",signif(iaelevs,2),"]",sep="") }
	
	otab1 <- calcBRAIDconfint(brdAnalysis$braidFit,function(parv) invertBRAIDrsm(iaelevs[1],DB=pylevs,parv=parv))
	otab2 <- calcBRAIDconfint(brdAnalysis$braidFit,function(parv) invertBRAIDrsm(iaelevs[1],DA=pxlevs,parv=parv))
	if (!is.null(control$ytrans)) {
		otab2 <- array(control$ytrans(otab2),dim=dim(otab2))
		orows1 <- signif(control$ytrans(pylevs),3)
	} else { orows1 <- signif(pylevs,3) }
	if (!is.null(control$xtrans)) {
		otab1 <- array(control$xtrans(otab1),dim=dim(otab1))
		orows2 <- signif(control$xtrans(pxlevs),3)
	} else { orows2 <- signif(pxlevs,3) }
	if (!is.character(ylab2)) { ocols1 <- c(gsub("(\\s)|(\")","",capture.output(print(ylab2))),"Est") }
	else { ocols1 <- c(ylab2,"Est") }
	if (!is.character(xlab2)) { ocols2 <- c(gsub("(\\s)|(\")","",capture.output(print(xlab2))),"Est") }
	else { ocols2 <- c(xlab2,"Est") }
	
	if (control$return) {
		tblist <- list()
		tblist[[1]] <- ptab1
		tblist[[2]] <- itab1
		rownames(tblist[[2]]) <- irows
		tblist[[3]] <- itab2
		rownames(tblist[[3]]) <- irows
		tblist[[4]] <- cbind(orows1,otab1)
		tblist[[5]] <- cbind(orows2,otab2)
	}
	if (control$plot) {
		grid.newpage()
		tbgb <- list()
		tt <- gridExtra::ttheme_default(core=list(fg_params=list(fontsize=6,parse=TRUE),padding=unit(c(2,2),"mm")),
										colhead=list(fg_params=list(fontsize=6,parse=TRUE),padding=unit(c(2,2),"mm")),
										rowhead=list(fg_params=list(fontsize=6,parse=TRUE),padding=unit(c(2,2),"mm")))
		tbgb[[1]] <- addTableTitle(gridExtra::tableGrob(as.data.frame(ptab2),rows=rownames(ptab2),cols=NULL,theme=tt),"Best BRAID Fit")
		rg <- rectGrob(gp=gpar(col="white"))
		itab1 <- array(paste(signif(itab1[,2],3),"~(",signif(itab1[,1],3),"-",signif(itab1[,3],3),")",sep=""),dim=c(nrow(itab1),1))
		itab2 <- array(paste(signif(itab2[,2],3),"~(",signif(itab2[,1],3),"-",signif(itab2[,3],3),")",sep=""),dim=c(nrow(itab2),1))
		tbgb[[2]] <- gridExtra::arrangeGrob(addTableTitle(gridExtra::tableGrob(as.data.frame(itab1),
						rows=irows,cols=NULL,theme=tt),"Additive Indices"),addTableTitle(gridExtra::tableGrob(as.data.frame(itab2),
						rows=irows,cols=NULL,theme=tt),"BRAID Indices"),nrow=2)
		otab1 <- array(paste(signif(otab1[,2],3),"~(",signif(otab1[,1],3),"-",signif(otab1[,3],3),")",sep=""),dim=c(nrow(otab1),1))
		otab2 <- array(paste(signif(otab2[,2],3),"~(",signif(otab2[,1],3),"-",signif(otab2[,3],3),")",sep=""),dim=c(nrow(otab2),1))
		otab1 <- cbind(orows1,otab1)
		otab2 <- cbind(orows2,otab2)
		if (!is.null(control$levtext)) { raethr <- control$levtext[1] }
		else if (!is.null(control$ztrans)) { raethr <- round(control$ztrans(iaelevs[1])) }
		else { raethr <- round(iaelevs[1]) }
		if (!is.null(control$xunit)) {
			if (class(control$xunit)=="expression")  {
				otitle1 <- bquote(EC[.(raethr)]~of~.(compounds[1])~(.(control$xunit[[1]])))
			} else { otitle1 <- bquote(EC[.(raethr)]~of~.(compounds[1])~(.(control$xunit))) }
		} else { otitle1 <- bquote(EC[.(raethr)]~of~.(compounds[1])) }
		if (!is.null(control$yunit)) {
			if (class(control$yunit)=="expression")  {
				otitle2 <- bquote(EC[.(raethr)]~of~.(compounds[2])~(.(control$yunit[[1]])))
			} else { otitle2 <- bquote(EC[.(raethr)]~of~.(compounds[2])~(.(control$yunit))) }
		} else { otitle2 <- bquote(EC[.(raethr)]~of~.(compounds[2])) }
		tbgb[[3]] <- gridExtra::arrangeGrob(addTableTitle(gridExtra::tableGrob(as.data.frame(otab1),rows=NULL,cols=ocols1,theme=tt),
						otitle1),addTableTitle(gridExtra::tableGrob(as.data.frame(otab2),rows=NULL,cols=ocols2,theme=tt),otitle2),nrow=2)
		
		grid.draw(gridExtra::arrangeGrob(plts[[1]],tbgb[[1]],tbgb[[2]],plts[[2]],plts[[3]],plts[[4]],plts[[5]],
								plts[[6]],plts[[7]],rg,rg,rg,tbgb[[3]],plts[[8]],plts[[9]],heights=c(13,13,13,1,14),
								top=paste(compounds[1],"vs.",compounds[2]),ncol=3))
	}
	if (control$return) { invisible(c(plts,tblist)) }
	else { invisible(NULL) }
}

getAdditiveParameters <- function(fpar,h1par=NULL,h2par=NULL) {
	addpar <- fpar
	addpar[c(5,6)] <- c(1,0)
	isIncr <- sign(addpar[10]-addpar[7])
	if (!is.null(h1par)) {
		addpar[1] <- exp(h1par[4])
		addpar[c(3,7,8)] <- h1par[1:3]
	}
	if (!is.null(h2par)) {
		# par2 <- coef(brdAnalysis$hfit2$allfits[[brdAnalysis$hfit2$bestModIdx]])
		addpar[2] <- exp(h2par[4])
		addpar[c(4,7,9)] <- h2par[1:3]
		if (!is.null(h1par)) { addpar[7] <- isIncr*min(isIncr*h1par[2],isIncr*h2par[2]) }
	}
	if (addpar[10]==fpar[8] || addpar[10]==fpar[9] ||
			isIncr*addpar[10]<max(isIncr*addpar[8],isIncr*addpar[9])) {
		addpar[10] <- isIncr*max(isIncr*addpar[8],isIncr*addpar[9])
	}
	return(addpar)
}

getAdditiveIndexTable <- function(brdAnalysis,iaelevs,macs,same) {
	ostart <- brdAnalysis$braidFit$ostart
	fixed <- brdAnalysis$braidFit$fixed
	fullpar <- brdAnalysis$braidFit$fullpar
	isIncr <- sign(fullpar[10]-fullpar[7])
	parv2fullpar <- function(parv) {
		tfullpar <- fullpar
		tfullpar[which(is.na(fixed))] <- parv
		if (is.na(fixed[10])) {
			if (is.na(fixed[8])||is.na(fixed[9])) {
				if (!is.na(fixed[8]) && ostart[8]==ostart[10]) { tfullpar[8] <- tfullpar[10] }
				else if (!is.na(fixed[9]) && ostart[9]==ostart[10]) { tfullpar[9] <- tfullpar[10] }
				else { tfullpar[10] <- tfullpar[10]+isIncr*max(isIncr*tfullpar[8:9]) }
			}
			else if (ostart[8]==ostart[10] && ostart[9]==ostart[10]) { tfullpar[8:9] <- tfullpar[10] }
		} else if ((is.na(fixed[8]) || is.na(fixed[9])) && ostart[10]==isIncr*max(isIncr*ostart[8:9])) {
			tfullpar[10] <- isIncr*max(isIncr*tfullpar[8:9])
		}
		return(tfullpar)
	}
	
	nboot <- nrow(brdAnalysis$braidFit$bCoefs)
	if (!is.null(brdAnalysis$hfit1)) { nboot <- min(nboot,nrow(brdAnalysis$hfit1$bCoefs)) }
	else { h1par <- NULL }
	if (!is.null(brdAnalysis$hfit2)) { nboot <- min(nboot,nrow(brdAnalysis$hfit2$bCoefs)) }
	else { h2par <- NULL }
	fpar <- brdAnalysis$braidFit$fullpar
	if (!is.null(brdAnalysis$hfit1)) { h1par <- coef(brdAnalysis$hfit1$allfits[[brdAnalysis$hfit1$bestModIdx]]) }
	if (!is.null(brdAnalysis$hfit2)) { h2par <- coef(brdAnalysis$hfit2$allfits[[brdAnalysis$hfit2$bestModIdx]]) }
	av <- calculateIAE(getAdditiveParameters(fpar,h1par,h2par),iaelevs,macs,same)
	avmat <- av
	for (i in 1:nboot) {
		fpar <- parv2fullpar(brdAnalysis$braidFit$bCoefs[i,])
		if (!is.null(brdAnalysis$hfit1)) { h1par <- brdAnalysis$hfit1$bCoefs[i,] }
		if (!is.null(brdAnalysis$hfit2)) { h2par <- brdAnalysis$hfit2$bCoefs[i,] }
		addpar <- getAdditiveParameters(fpar,h1par,h2par)
		if (length(iaelevs)>1) { avmat <- cbind(avmat,calculateIAE(getAdditiveParameters(fpar,h1par,h2par),iaelevs,macs,same)) }
		else { avmat <- c(avmat,calculateIAE(getAdditiveParameters(fpar,h1par,h2par),iaelevs,macs,same)) }
	}
	if (length(iaelevs)==1) { avmat <- array(avmat,dim=c(1,length(avmat))) }
	cimat <- apply(avmat,1,quantile,brdAnalysis$braidFit$ciLevs)
	if (length(iaelevs)==1) { return(array(c(cimat[1],av,cimat[2]),dim=c(1,3))) }
	else { return(cbind(cimat[1,],av,cimat[2,])) }
}

makeParamTable <- function(bsumm,fpar,parse=TRUE) {
	if (parse) { cistr <- function(v) paste(v[2],"~(",v[1],"-",v[3],")",sep="") }
	else { cistr <- function(v) paste(v[2]," (",v[1]," - ",v[3],")",sep="") }
	rnms <- row.names(bsumm)
	bsumm[,1] <- signif(bsumm[,1],3)
	bsumm[,2] <- signif(bsumm[,2],3)
	bsumm[,3] <- signif(bsumm[,3],3)
	tnms <- c("ID['M,A']","ID['M,B']","n[a]","n[b]","E[0]")
	tab <- rep("",4)
	for (i in 1:4) { tab[i] <- cistr(bsumm[i,]) }
	if ("E0" %in% rnms) { tab <- c(tab,cistr(bsumm["E0",])) }
	else { tab <- c(tab,as.character(fpar[7])) }
	if (!("EfAB" %in% rnms)) {
		if ("EfA" %in% rnms || "EfB" %in% rnms) {
			tnms <- c(tnms,"E['f,A']","E['f,B']")
			if ("EfA" %in% rnms) { tab <- c(tab,cistr(bsumm["EfA",])) }
			else { tab <- c(tab,as.character(fpar[8])) }
			if ("EfB" %in% rnms) { tab <- c(tab,cistr(bsumm["EfB",])) }
			else { tab <- c(tab,as.character(fpar[9])) }
			if (fpar[10]!=fpar[8] && fpar[10]!=fpar[9]) {
				tnms <- c(tnms,"E[f]")
				tab <- c(tab,as.character(fpar[10]))
			}
		} else {
			tnms <- c(tnms,"E[f]")
			tab <- c(tab,as.character(fpar[10]))
		}
	} else {
		if ("EfA" %in% rnms) {
			tnms <- c(tnms,"E['f,A']")
			tab <- c(tab,cistr(bsumm["EfA",]))
		}
		if ("EfB" %in% rnms) {
			tnms <- c(tnms,"E['f,B']")
			tab <- c(tab,cistr(bsumm["EfB",]))
		}
		tnms <- c(tnms,"E[f]")
		tab <- c(tab,cistr(bsumm["EfAB",]))
	}
	if ("kappa" %in% rnms) {
		tnms <- c(tnms,"kappa")
		tab <- c(tab,cistr(bsumm["kappa",]))
	}
	if ("delta" %in% rnms) {
		tnms <- c(tnms,"delta")
		tab <- c(tab,cistr(bsumm["delta",]))
	}
	tab <- array(tab,dim=c(length(tab),1))
	row.names(tab) <- tnms
	return(tab)
}

addTableTitle <- function(tg,titleStr="Title",fs=7) {
	ttgrob <- textGrob(titleStr,gp=gpar(fontsize=fs))
	wt <- max(convertWidth(grobWidth(ttgrob),"mm",FALSE),convertWidth(sum(tg$widths),"mm",FALSE))
	ht1 <- convertHeight(grobHeight(ttgrob),"mm",FALSE)+unit(3,"mm")
	ht2 <- convertHeight(sum(tg$heights),"mm",FALSE)
	tabgrob <- gtable::gtable(wt,unit.c(ht1,ht2))
	tabgrob <- gtable::gtable_add_grob(tabgrob,ttgrob,1,1)
	tabgrob <- gtable::gtable_add_grob(tabgrob,tg,2,1)
	return(tabgrob)
}
