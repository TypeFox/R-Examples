# plot=TRUE; auto.layout=TRUE; plot.se.profile=TRUE; force.fit=FALSE; test.LOD=FALSE; find.LOD=FALSE; find.best.dilution=FALSE; return.fits=FALSE; grid.len=NULL; unk.replicate=NULL; verbose=FALSE; robust="mean"; fit.4pl=FALSE; lod.ci=95; plot.log="x"; control.jags=list(n.iter=1e5, jags.seed=1, n.thin=1e1, keep.jags.samples=FALSE); xlab=NULL; ylab=NULL; var.model="constant"
ncal <- function(...) UseMethod("ncal") 
ncal.formula = function (formula, data, 
    bcrm.fit=FALSE, bcrm.model="t4", robust="mean", bcrm.prior="default",
    force.fit=TRUE, fit.4pl=FALSE, return.fits=FALSE, 
    plot=TRUE, auto.layout=TRUE, plot.se.profile=TRUE, plot.log="x", plot.unknown=TRUE,
    test.LOD=FALSE, find.LOD=FALSE, find.LOQ=FALSE, grid.len=50, lod.ci=95,
    unk.replicate=NULL, find.best.dilution=FALSE, unk.median=FALSE,
    control.jags=list(n.iter=1e5, jags.seed=1, n.thin=NULL, 
        keep.jags.samples=FALSE, n.adapt=1e3), 
    cex=.5, additional.plot.func=NULL, check.out.of.range=1, xlab=NULL, ylab=NULL, main=NULL,
    var.model=c("constant","power"), log.both.sides=FALSE,
    control.crm.fit=list(max.iter=20),
    verbose=FALSE, 
...)
{

    var.model=match.arg(var.model)
    
    # populate default values for control.jags
    if (is.null(control.jags$n.iter)) control.jags$n.iter=1e5
    if (is.null(control.jags$jags.seed)) control.jags$jags.seed=1
    if (is.null(control.jags$n.thin)) control.jags$n.thin=NULL
    if (is.null(control.jags$keep.jags.samples)) control.jags$keep.jags.samples=T
    if (is.null(control.jags$n.adapt)) control.jags$n.adapt=1e3
    if (is.null(control.crm.fit$max.iter)) control.crm.fit$max.iter=20
    
    # for backward compatibility with the rumi function call
    if (is.null(formula)) {
        formula=log(fi) ~ expected_conc
    }
    
#    if (is.null(data$starting_conc) & is.null(data$expected_conc)) stop("ERROR: starting_conc and expected_conc missing from data\n\n")
#    if (is.null(data$expected_conc)) data$expected_conc=data$starting_conc/data$dilution
    
   # if (any(is.na(model.matrix(formula, data[data$well_role=="Standard",])[,2]))) {
   #     stop("some predictors are NA")
   #}        ## this wont work for WLS formula 
    
    if (all.names(formula)[2]=="log") log.transform=T else log.transform=F
    
    outcome.coln=all.vars(formula)[1]
    predictor.coln=all.vars(formula)[2]
   
    xlab=ifelse(is.null(xlab),"Concentration",xlab)
    ylab=ifelse(is.null(ylab),ifelse(log.transform,"log(","")%+%outcome.coln%+%ifelse(log.transform,")",""), ylab)
    
    ## somehow not working, replace with some dummy code     
    #ylab=ifelse(is.null(ylab), "dummy label", ylab)
    
    if (plot.se.profile) find.LOQ=T
    
    if (log.transform) {
        if (any(data[[outcome.coln]]<=0)) warning(outcome.coln%+%" var cannot be 0 or less, will set to 1")
        data[[outcome.coln]][data[[outcome.coln]]<=0]=1
    }
                    
    # convert bead_type to analyte, backward compatibility
    if (is.null(data$analyte) & !is.null(data$bead_type)) data$analyte=data$bead_type
    
    # checking the data
    if (is.null(data$assay_id)) stop("assay_id missing from data")
    if (is.null(data$analyte)) stop("analyte missing from data")
    if (is.null(data$well_role)) {
        if (!is.null(data$sample_id)) stop("well_role missing from data")
        # otherwise, assume all are standard wells
        #warning("no well_role column, assume all samples are standard samples")
        data$well_role="Standard"
    } else {
        if (all(data$well_role!="Standard")) stop("There are not Standard wells.")
    }
    if (is.null(data$sample_id) & any(data$well_role!="Standard")) {
        stop("sample_id missing from data with not just standard samples\n\n")
    }
    
    if (is.null(data[[outcome.coln]])) {cat("ERROR: y. missing from data\n\n"); stop()}
    if (!all(is.numeric(data[[outcome.coln]]))) stop("ERROR: y. column not numeric\n\n")
    
    if (!plot) plot.se.profile=F
    
    # filter out rows with fi being NA
    data=data[!is.na(data[[outcome.coln]]),]
    
    ana=sort(unique(data$analyte))
    #ana=sort(setdiff(unique(data$analyte), "blank"))
    assays=sort(setdiff(unique(data$assay_id), "blank"))
    high.low=NULL
    LOQ=NULL
    LOQ.30=NULL
    out=data.frame()
    fits=list()
    for (a in ana) {
        
        # get fits
        if (bcrm.fit) {
            if (verbose) {cat("\n");myprint (a)}        
            dat.std.bcrm=data[data$well_role=="Standard" & data$analyte==a,] 
            if (nrow(dat.std.bcrm)==0) next    
            fits = bcrm (formula, data=dat.std.bcrm, parameterization="gh", error.model=bcrm.model, prior=bcrm.prior, mean.model=ifelse(fit.4pl,"4PL","5PL"),
                n.iter=control.jags$n.iter, jags.seed=control.jags$jags.seed, n.thin=control.jags$n.thin, keep.jags.samples=control.jags$keep.jags.samples, 
                T.init=NULL, prior.sensitivity=NULL, n.adapt=control.jags$n.adapt, verbose=verbose)
        } 
        
        for (p in assays) {
            if (verbose) myprint(a)    
            dat.a.p=data[data$assay_id==p & data$analyte==a,]
            if (nrow(dat.a.p)==0) next    
            dat.std= dat.a.p[dat.a.p$well_role=="Standard",] 
            if (nrow(dat.std)==0) next
            
            std.low =log(min (dat.std[[predictor.coln]]))
            std.high=log(max (dat.std[[predictor.coln]]))
            # y.low and y.high are setting range in plotting
            y.low=min(data[[outcome.coln]]); if (log.transform) y.low=log(y.low)
            y.high=max(data[[outcome.coln]]); if (log.transform) y.high=log(y.high)
            
            # add columns like fi.avg
            if (!outcome.coln%+%".avg" %in% names(dat.std)) {
                stand.avg=aggregate(dat.std[[outcome.coln]], by = list(dat.std$analyte, dat.std[[predictor.coln]], dat.std$assay_id), mean)
                names(stand.avg) = c("analyte", predictor.coln,"assay_id", outcome.coln%+%".avg")
                dat.std = merge(dat.std, stand.avg, by=c("analyte", predictor.coln,"assay_id"), all.x=FALSE, all.y=TRUE)
            }
            
            # fit standard curves
            if (bcrm.fit) {
                fit = get.single.fit(fits, assay_id=p) 
            } else {
                fit = crm.fit(formula, data = dat.std, var.model=var.model, robust=robust, fit.4pl=fit.4pl, verbose=verbose, max.iter=control.crm.fit$max.iter, log.both.sides=log.both.sides)
            }
            if (return.fits & !bcrm.fit) fits[[p%+%a]]=fit
            if (verbose>2) print("debug 100")
            
            out.ana=data.frame()
            if (auto.layout & plot) {
                if (plot.se.profile) par(mfrow=c(2,2)) else par(mfrow=c(1,1))
            }
    
            if (!is.null(fit)) {
    
                df.=nrow(dat.std)
                dd.95=qt(1-(1-95/100)/2, df=df.-ifelse(fit.4pl,4,5))
                dd=qt(lod.ci/100, df=df.-ifelse(fit.4pl,4,5))
        
                # sometimes, even though fit is not NULL, some standard errors of the parameter estimates are negative
                if (verbose>2) print("debug 300")
                bad.se=fit[["bad.se"]]
                if (verbose & bad.se) print("bad.se")
                #if (bad.se & plot.se.profile) empty.plot()
                
                coef.tmp=coef(fit,type="classical")
                if("g" %in% names(coef.tmp)) {
                    #gnls.fit returns g-h parameterization
                    coef.tmp=gh2cla(coef.tmp)
                }
                incr. = coef.tmp[startsWith(names(coef.tmp),"b")]<0 
                if (verbose>2) print("debug 500")
    
                # estimate unknown concentrations
                dat.unk=dat.a.p[dat.a.p$well_role!="Standard",]
                if (nrow(dat.unk)>0 & is.null(dat.unk$dilution)) dat.unk$dilution=1
                sam = unique(dat.unk$sample_id)
                #if (!is.null(sam)) sam=sort(sam) # do not think it is a good idea anymore. v15-16
                for (s in sam) {
                    
                    if (verbose) myprint(s)                    
                    dil=sort(unlist(unique(dat.unk[dat.unk$sample_id==s, "dilution"])))
                    
                    out.s=data.frame()
                    for (dil. in dil) {
                        if (verbose) myprint(dil.)
                        # for each sample/dilution, there may be replicates
                        dat.unk.s.d = dat.unk[dat.unk$dilution == dil. & dat.unk$sample_id == s,]
                        if (is.null(unk.replicate)) unk.replicate=nrow(dat.unk.s.d)
                        y.=dat.unk.s.d[[outcome.coln]]
                        if (log.transform) y.=log(y.)
                        
                        estimated = getConc(fit, ifelse(unk.median, median(y.), mean(y.)), n.replicate=length(y.), x.range=exp(c(std.low,std.high)), y.range=c(y.low,y.high), 
                            verbose=verbose, check.out.of.range=check.out.of.range)
                        estimated = unname(estimated)
                        if (verbose) print(estimated)
    
                        if (!bad.se) {
    
                            # test limits of detection
                            if (test.LOD) {
                                test.1.not.rejected=(estimated[1] - std.low)/estimated[2] < dd | estimated[2]==Inf
                                test.2.not.rejected=(std.high - estimated[1])/estimated[2] < dd | estimated[2]==Inf
                                if(test.1.not.rejected & !test.2.not.rejected) estimated[1] = std.low
                                if(test.2.not.rejected & !test.1.not.rejected) estimated[1] = std.high
                                # if both rejected, set to the closer one
                                if(test.1.not.rejected & test.2.not.rejected) estimated[1] = ifelse (estimated[1] > log((exp(std.low)+exp(std.high))/2), std.high, std.low)
                                # set std err to Inf for those touched
                                if(test.1.not.rejected | test.2.not.rejected) estimated[2] = Inf
                            }
                            
                            out.s=rbind (out.s, data.frame (c(dat.unk.s.d[1,], "est.log.conc"=estimated[1]+log(dil.), "se"=estimated[2], 
                                "est.conc"=exp(estimated[1]+log(dil.)), "lb.conc"=exp(estimated[1]+log(dil.)-dd.95*estimated[2]), "ub.conc"=exp(estimated[1]+log(dil.)+dd.95*estimated[2])  ),
                                stringsAsFactors =FALSE))
                        
                        } else {
                            out.s=rbind (out.s, data.frame (c(dat.unk.s.d[1,], "est.log.conc"=estimated[1]+log(dil.), "se"=NA, 
                                "est.conc"=exp(estimated[1]+log(dil.)), "lb.conc"=NA, "ub.conc"=NA  ), stringsAsFactors =FALSE))
                        }
                        # let outcome.coln be the average
                        out.s[nrow(out.s),outcome.coln]=ifelse(log.transform, exp(mean(y.)), mean(y.))
                    }
                    if (!find.best.dilution) {
                        if(length(dil)>1) warning ("There are most than one dilutions for this sample, and we are returning all. sample_id: "%+%s)
                        out.ana=rbind(out.ana, out.s)
                    } else {
                        out.ana=rbind(out.ana, out.s[which.min(out.s$se), ])
                    }
                }
                
                if (plot) {
                    #if (!is.null(dat.std$replicate)) pch=c(1,19:15)[dat.std$replicate] else 
                    pch=rep(1,nrow(dat.std))
                    suppressWarnings(
#                      if (!WLS) {
                        plot(fit, type=ifelse(bcrm.fit,"p","all"), main=ifelse(is.null(main),p%+%", "%+%a,main), cex=cex, log=plot.log, xlab=xlab, pch=pch, ylab=ylab, ...)      
#                      } else {
#                        conc.est.grid=myplot.gnls (fit, myDat=dat.std, main=ifelse(is.null(main),p%+%", "%+%a,main), x1=NULL, xn=NULL, y1=NULL, yn=NULL, 
#                            check.out.of.range=check.out.of.range, unk.replicate=1, verbose=verbose)
#                      }      
                    )
                    
                    if (!is.null(additional.plot.func)) additional.plot.func()
                    
                    if (plot.unknown) {
#                        if (!WLS) {
                            suppressWarnings( plot(fit, type=ifelse(bcrm.fit,"p","all"), main=ifelse(is.null(main),p%+%", "%+%a,main), cex=cex, pch=pch,log=plot.log,xlab=xlab,ylab=ylab, ...) )
#                        } else {
#                            # do nothing b/c plot.gnls makes two plots and compute conc.est.grid
#                        }
                        # add unknown samples
                        if (nrow(out.ana)>0) {
                            if(log.transform) y.tmp=log(out.ana[[outcome.coln]]) else y.tmp=out.ana[[outcome.coln]]
                            points(out.ana$est.conc/out.ana$dilution, y.tmp, pch="*", col=2)
                        }
                    }
                }
            
                # if there is no sample, unk.replicate may still be null
                if (is.null(unk.replicate)) unk.replicate=1
                    
                # plot se profile file, compute lod, loq
                if (!bad.se & (find.LOD | plot.se.profile | find.LOQ)) {
                
                    if(verbose) print("nCal 600")
                    # LoDi
                    fit.low=FivePL.t(std.low, coef(fit))
                    fit.high=FivePL.t(std.high, coef(fit))
                    tmp.flag = fit.low>fit.high
                    if (!incr.) tmp.flag=!tmp.flag
                    if (tmp.flag) {
                        if (verbose) print("out of here")
                        next
                    }
                    mid.log.fi = unname((fit.low+fit.high)/2)
                    if(verbose) print("nCal 620")
                    # low end   
#                    if (!WLS) {
                      
                        x.hat.low=getConc(fit, seq(fit.low, mid.log.fi, length=grid.len+1)[-1], n.replicate=rep(unk.replicate,grid.len), x.range=exp(c(std.low,std.high)), 
                                            y.range=c(y.low,y.high), verbose=verbose, check.out.of.range=check.out.of.range)[,c(1,2,6,7,8,4,5)]
                          
                        same.as.low = x.hat.low[,1] - dd*x.hat.low[,2] < std.low
                        #same.as.low = exp(x.hat.low[,1]) - dd * exp(x.hat.low[,1]) * x.hat.low[,2] < exp(std.low)
                        same.as.low.rle=rle(same.as.low)
                        tmp=nrow(x.hat.low)-last(same.as.low.rle$lengths)
                        x.low=x.hat.low[ifelse(tmp>0,tmp,1), 1]
                        if(verbose) print("nCal 640")
                        # high end
                        x.hat.high=getConc(fit, seq(fit.high, mid.log.fi, length=grid.len+1)[-1], n.replicate=rep(unk.replicate,grid.len), x.range=exp(c(std.low,std.high)), 
                            y.range=c(y.low,y.high), check.out.of.range=check.out.of.range)[,c(1,2,6,7,8,4,5)]
                        same.as.high = x.hat.high[,1] + dd * x.hat.high[,2] > std.high
                        #same.as.high = exp(x.hat.high[,1]) + dd * exp(x.hat.high[,1]) * x.hat.high[,2] > exp(std.high)
                        same.as.high.rle=rle(same.as.high)
                        tmp=nrow(x.hat.high)-last(same.as.high.rle$lengths)
                        x.high=x.hat.high[ifelse(tmp>0,tmp,1), 1]
                        high.low=rbind(high.low,c(p,a,x.high,x.low))
                        conc.est.grid=rbind(x.hat.low, x.hat.high[nrow(x.hat.high):1,])
                    
#                    } else {
#                        conc.est.grid = conc.est.grid[,c(1,2,6,7,8,4,5)]                    
#                    }
                    
                    # plot vertical lines at x.high and x.low
                    if (find.LOD & plot.se.profile) {
                        abline(v=exp(x.low))
                        abline(v=exp(x.high))
                    }
                    
#                    print(sum(same.as.low.rle[[2]]))
#                    print(sum(same.as.high.rle[[2]]))
#                    if (same.as.low.rle[[2]]!=1 | same.as.high.rle[[2]]!=1) {
#                        stop("rle not 1")
#                    }
    
                    if (find.LOQ) {
                    
                        if(verbose>2) print("debug 700")
                        
                         conc.est.grid=conc.est.grid[conc.est.grid[,2]!=Inf & !is.na(conc.est.grid[,2]),]
                        
                        #conc.est.grid=subset(conc.est.grid, log.conc>=0.05)
                        # define se profile as se/estimated conc on the conc scale
                        se.profile=conc.est.grid[,5]/exp(conc.est.grid[,1])*100
#                        # define se profile as percent cv vs estimated conc, using length of CI, similar to last
#                        se.profile=(conc.est.grid[,7]-conc.est.grid[,6])/4/exp(conc.est.grid[,1])*100
#                        # define se profile as se/estimated conc on log scale
#                        se.profile=conc.est.grid[,2]/conc.est.grid[,1]*100
                                                
                        # LOQ is defined on linear scale
                        rle1 = rle(abs(unname(se.profile))<20)
                        if (sum(rle1[[2]])==1) {
                            loq.ind=which(rle1[[2]])
                            LOQ=rbind(LOQ, c(p,a,exp(conc.est.grid[c(ifelse(loq.ind==1, 1, 1+cumsum(rle1[[1]])[loq.ind-1]), cumsum (rle1[[1]])[loq.ind]),1])))
                        } else if (sum(rle1[[2]])==0) {
                            LOQ=rbind(LOQ, c(p,a,NA,NA))
                        } else {
                            LOQ=rbind(LOQ, c(p,a,NA,NA))
                            #stop("rle1 not quite right, "%+%p%+%a)
                        }
                        
                        rle30 = rle(abs(unname(se.profile))<30)
                        if (sum(rle30[[2]])==1) {
                            loq.ind=which(rle30[[2]])
                            LOQ.30=rbind(LOQ.30, c(p,a,exp(conc.est.grid[c(ifelse(loq.ind==1, 1, 1+cumsum(rle30[[1]])[loq.ind-1]), cumsum (rle30[[1]])[loq.ind]),1])))
                        } else if (sum(rle30[[2]])==0) {
                            LOQ.30=rbind(LOQ.30, c(p,a,NA,NA))
                        } else {
                            LOQ.30=rbind(LOQ.30, c(p,a,NA,NA))
                            #stop("rle30 not quite right, "%+%p%+%a)
                        }
                        
                    }
                    
                    # plot error profile
                    if (plot.se.profile) {
                        
                        #abline(v=c(exp(x.low), exp(x.high)), col=1)
                        abline(v=LOQ[nrow(LOQ),3:4], col="gray")
                        #abline(v=LOQ.30[nrow(LOQ.30),3:4], col=3)
                        mylegend(legend=c("LOQ"), col=c("gray",2,1), lty=c(1,0), pch=c("","*"), x=2)
                        #mylegend(legend=c("LOQ","estimate of unknowns"), col=c("gray",2,1), lty=c(1,0), pch=c("","*"), x=2)                                            
    
                        # se vs conc
                        plot(exp(conc.est.grid[,1]), conc.est.grid[,2],type="n", xlab="estimate", ylab="variance ("%+%unk.replicate%+%" replicate"%+%ifelse(unk.replicate>1,"s","")%+%")", log="x", 
                            #ylim=c(0,max(conc.est.grid[,2])^2), 
                            ylim = c(0, min(c(3, max(conc.est.grid[,2])^2))),
                            xlim=range(dat.std[[predictor.coln]]))
                        lines(exp(conc.est.grid[,1]), (conc.est.grid[,4]),type="l", col="blue")
                        lines(exp(conc.est.grid[,1]), (conc.est.grid[,3]),type="l", col="red")
                        lines(exp(conc.est.grid[,1]), conc.est.grid[,2]^2,type="l", col="black")
                        mylegend(legend=c("total","if curve is perfectly known","if y is perfectly known"),col=c("black","red","blue"), lty=1, x=2)
                        
                        # %CV vs conc
                        plot(exp(conc.est.grid[,1]), se.profile,type="l", xlab="estimate", ylab="% CV ("%+%unk.replicate%+%" replicate"%+%ifelse(unk.replicate>1,"s","")%+%")", log="x", 
                            ylim=c(0,100), xlim=range(dat.std[[predictor.coln]]))
                        abline(h=20, col="gray")
                        #abline(h=30, col=3)
                        #abline(h=40, col=4)
                        #abline(h=50, col=6)
                        #abline(v=c(exp(x.low), exp(x.high)))
                        mylegend(legend=c("LOQ"), col=c("gray",3,4,6,1), lty=1, x=2)
                        
                    }    
                    
                } else {
                    if (find.LOD) high.low=rbind(high.low, c(p,a,NA,NA))
                    if (find.LOQ) {
                        LOQ=rbind(LOQ, c(p,a,NA,NA))
                        LOQ.30=rbind(LOQ.30, c(p,a,NA,NA))
                    }
    
                }
    
            } else {
    
                if (plot) {
                    plot (formula, dat.std, log=plot.log, main="FAILED: "%+%p%+%", "%+%a, cex=.1)
                    if (plot.se.profile) empty.plot()
                }                
                bad.se=T # needed after this if statement
                
            }
            
            out=rbind(out, out.ana)
                
        } # end for loop plate
    
    } # end for loop analyte
    
    if (find.LOD) {
        if (!is.null(LOD)) {
            LOD=as.data.frame(high.low, stringsAsFactors=FALSE) # without stringsAsFactors=FALSE, sometimes x.high/x.low are set to factor
            names(LOD)=c("assay","analyte","rlodi","llodi")
            LOD$llodi=exp(as.double(LOD$llodi))
            LOD$rlodi=exp(as.double(LOD$rlodi))
            LOD=LOD[,c("assay", "analyte", "llodi", "rlodi")]
        }
        attr(out, "LOD")=LOD
    }
    if (find.LOQ) {
        if (!is.null(LOQ)) {
            LOQ=as.data.frame(LOQ, stringsAsFactors=FALSE) # without stringsAsFactors=FALSE, sometimes x.high/x.low are set to factor
            names(LOQ)=c("assay","analyte","lloq","rloq")
            LOQ$lloq=as.double(LOQ$lloq)
            LOQ$rloq=as.double(LOQ$rloq)
        }
        attr(out, "LOQ")=LOQ
        
#        LOQ.30=as.data.frame(LOQ.30, stringsAsFactors=FALSE) # without stringsAsFactors=FALSE, sometimes x.high/x.low are set to factor
#        names(LOQ.30)=c("assay","analyte","lloq30","rloq30")
#        LOQ.30$lloq30=as.double(LOQ.30$lloq30)
#        LOQ.30$rloq30=as.double(LOQ.30$rloq30)
#        attr(out, "LOQ30")=LOQ.30
    }
    if (return.fits) attr(out, "fits")=fits
    return (out)
}

ncal.character = function (file, is.luminex.xls, formula, bcrm.fit, verbose=FALSE, ...) {
    if (is.luminex.xls) dat=read.luminex.xls(file, verbose) else dat=read.csv(file, as.is=TRUE)
    #mywrite.csv(dat, file=file)
    ncal.formula(formula, dat, bcrm.fit=bcrm.fit, ...)
}

rumi = function(data, ...) ncal.formula(formula=log(fi)~expected_conc, data, ...)
    


# n.replicate=1; check.out.of.range=1; x.range=NULL; y.range=NULL; verbose=FALSE
getConc=function(fit, y, n.replicate=1, check.out.of.range=1, x.range=NULL, y.range=NULL, verbose=FALSE) {
    
    if (verbose>2) print("getConc 100")
    param=coef(fit,type="classical")
    names(param)=substr(names(param),1,1)
    b=param["b"]; c=param["c"]; d=param["d"]; e=param["e"]
    if (is.na(param["f"])) {
        fit.4pl=TRUE
        f=1
    } else {
        fit.4pl=FALSE
        f=param["f"]
    }    
        
    if (class(fit)[1]=="drc" | class(fit)[1]=="crm") {
                    
        if (verbose>2) print("getConc 200")
        
        # point estimate
        x = FivePL.x.inv(y, param) 
                    
        # variance part 1
        A=((d-c)/(y-c))^(1/f)-1
        B=(d-c)/(y-c)
        dx.dy = -e/(b*f)*A^{1/b-1}*B^{1/f+1}/(d-c)
        sigma2 = getVarComponent(fit)
        s1=dx.dy^2*sigma2/n.replicate /x^2
        
        if (verbose>2) print("getConc 300")
        # variance part 2
        dx.db = suppressWarnings(-log(A)*A^(1/b)*b^(-2)*e) # suppressWarnings to avoid NaN warning
        dx.dc = e/b*A^(1/b-1)*B^(1/f-1) *(1/f)*(d-y)/(y-c)^2
        dx.dd = e/b*A^(1/b-1)*B^(1/f-1) *(1/f)/(y-c)
        dx.de = A^(1/b)
        dx.df = e/b*A^(1/b-1)*B^(1/f) *suppressWarnings(log(B))*(-1/f^2)
        if(fit.4pl) {
            D.beta= cbind(dx.db, dx.dc, dx.dd, dx.de)    
        } else {
            D.beta= cbind(dx.db, dx.dc, dx.dd, dx.de, dx.df)
        }        
        V.beta = vcov(fit,type="classical")
        s2=diag(D.beta%*%V.beta%*%t(D.beta)) / x^2    
#        print(D.beta%*%V.beta) # from jags output, the correlation between e and f is very different from drm output. 
#        print(s2)
    
        if (verbose>2) print("getConc 400")
        se.t = sqrt(s1 + s2)
        lower.bound=suppressWarnings(exp(log(x)-2*se.t))
        upper.bound=suppressWarnings(exp(log(x)+2*se.t))
                
    } else if (class(fit)[1]=="bcrm") {
    
        if (verbose>2) print("getConc 500")
        params=fit$coef.samples
        x.samples = FivePL.x.inv(y=rep(y,each=nrow(params)), param=rep.matrix(params, times=length(y), by.row=TRUE)) 
        x.samples=matrix(x.samples,nrow=nrow(params))
        
        # need to do this here before computing s2, b/c if there is Inf in x.samples, s2 will be NaN
        
        x=apply(x.samples,2,median)        
        
        if (verbose>2) print("getConc 600")
        # variance part 1, assuming curve is perfectly known
        A=((d-c)/(y-c))^(1/f)-1
        B=(d-c)/(y-c)
        dx.dy = -e/(b*f)*A^{1/b-1}*B^{1/f+1}/(d-c)
        sigma2 = getVarComponent(fit)
        s1=dx.dy^2*sigma2/n.replicate /x^2                      
        
        if (verbose>2) print("getConc 700")
        # variance part 2, assuming y is perfectly known, based on MCMC samples of parameters
        t.samples=suppressWarnings(log(x.samples))
        s2=apply(t.samples,2,var,na.rm=TRUE)
    
        se.t = sqrt(s1 + s2)
        if (verbose) myprint(s1, s2, se.t)
        lower.bound=apply(t.samples,2,quantile,.025,na.rm=TRUE)
        lower.bound=lower.bound-1.96*sqrt(s1)
        lower.bound=lower.bound+1.96*sqrt(s1)
        lower.bound=exp(lower.bound)
        upper.bound=apply(t.samples,2,quantile,.975,na.rm=TRUE)
        upper.bound=upper.bound-1.96*sqrt(s1)
        upper.bound=upper.bound+1.96*sqrt(s1)
        upper.bound=exp(upper.bound)
    
        x=apply(x.samples,2,median)        
    } else if (class(fit)[1]=="gnls") {  ## for WLS
      if (verbose>2) print("getConc 200")
      
      # point estimate
      x = FivePL.x.inv(y, param) 
      
      # variance part 1
      A=((d-c)/(y-c))^(1/f)-1
      B=(d-c)/(y-c)
      dx.dy = -e/(b*f)*A^{1/b-1}*B^{1/f+1}/(d-c)
      
      sigma2 = summary(fit)$sigma^2
      g= abs(y)^coef(fit$modelStruct$varStruct)
      s1 = dx.dy^2 * sigma2 *g^2/n.replicate/x^2
           
      if (verbose>2) print("getConc 300")
      # variance part 2
      dx.db = suppressWarnings(-log(A)*A^(1/b)*b^(-2)*e) # suppressWarnings to avoid NaN warning
      dx.dc = e/b*A^(1/b-1)*B^(1/f-1) *(1/f)*(d-y)/(y-c)^2
      dx.dd = e/b*A^(1/b-1)*B^(1/f-1) *(1/f)/(y-c)
      dx.de = A^(1/b)
      dx.df = e/b*A^(1/b-1)*B^(1/f) *suppressWarnings(log(B))*(-1/f^2)
      if(fit.4pl) {
        D.beta= cbind(dx.db, dx.dc, dx.dd, dx.de)    
      } else {
        D.beta= cbind(dx.db, dx.dc, dx.dd, dx.de, dx.df)
      }        
      V.beta = vcov(fit,type="classical")
      s2 = diag(sigma2*D.beta %*% V.beta %*% t(D.beta))/x^2    
      #        print(D.beta%*%V.beta) # from jags output, the correlation between e and f is very different from drm output. 
      #        print(s2)
      
      if (verbose>2) print("getConc 400")
      se.t = sqrt(s1 + s2)
      #se.t = se.t[1,1]
      lower.bound=suppressWarnings(exp(log(x)-2*se.t))
      upper.bound=suppressWarnings(exp(log(x)+2*se.t))     
      
    }
    
    if (verbose>2) print("getConc 800")
    if (check.out.of.range==1) {
        # by construction, x.range has two elements, and the first is small and the second is large
        se.t[x>x.range[2] | x<x.range[1]]=Inf
        x[x>x.range[2]]=x.range[2]
        x[x<x.range[1]]=x.range[1]/2
    } else if (check.out.of.range==2) {
        se.t[x==Inf | x==0]=Inf
        x[x==Inf]=x.range[2]
        x[x==0]=x.range[1]/2
    } else stop ("check out of range mode not supported")
    
    # replace by se.t[x>x.range[2] | x<x.range[1]]=Inf
#    # if estimated conc outside estimated asymptote, set se to Inf
#    se1 = Inf
#    se.t=ifelse(y<c | y>d, se1, se.t)
    
    # may not be properly vectorized            
#    if(check.out.of.range) {
#        # if MFI outside standards MFI, set of Inf
#        if ( y < y.range[1] ) {
#            if (verbose) print("y < y.range[1]")
#            return (c("log.conc"=log(x.ninf), "s.e."=se1, "concentration"=NaN, "lower.bound"=NaN, "upper.bound"=NaN, "s1"=NaN, "s2"=NaN, "se.x"=NaN))
#        } else if ( y > y.range[2] ) {
#            if (verbose) print("y > y.range[2]")
#            return (c("log.conc"=log(x.inf), "s.e."=se1, "concentration"=NaN, "lower.bound"=NaN, "upper.bound"=NaN, "s1"=NaN, "s2"=NaN, "se.x"=NaN))
#        }
#    
#        # if estimated conc outside expected conc, set to min and max of standard conc, cannot use x.inf and x.ninf
#        if ( x < x.range[1] ) {
#            if (verbose) print("x < x.range[1]")
#            return (c("log.conc"=log(x.range[1]), "s.e."=se1, "concentration"=NaN, "lower.bound"=NaN, "upper.bound"=NaN, "s1"=NaN, "s2"=NaN, "se.x"=NaN))
#        } else if ( x > x.range[2] ) {
#            if (verbose) print("x > x.range[2]")
#            return (c("log.conc"=log(x.range[2]), "s.e."=se1, "concentration"=NaN, "lower.bound"=NaN, "upper.bound"=NaN, "s1"=NaN, "s2"=NaN, "se.x"=NaN))
#        }
#    }
    
    
    #need to fix log x warning
    if (verbose>2) print("getConc 900")
    res=drop(cbind(
        "log.conc"=unname(log(x)), 
        "s.e."=se.t, 
        "concentration"=unname(x), 
        "lower.bound"=lower.bound, 
        "upper.bound"=upper.bound, 
        "s1"=unname(s1), 
        "s2"=unname(s2), 
        "se.x"=unname(x*se.t) 
    ))
    return (res)
}
