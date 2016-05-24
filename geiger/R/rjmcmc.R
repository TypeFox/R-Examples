.rjmcmc.bm.multi=function(phy, dat, SE, ngen, samp, type=c("jump-rbm", "rbm", "jump-bm", "bm"),...){
	con=list(...)
	trees=phy
	fb=con$filebase
	nm=names(phy)
	if(is.null(nm)) nm=1:length(phy)
	if(!is.null(fb)) files=paste(fb, nm, sep="_") else files=nm
	FUN=lapply
	FUN(1:length(trees), function(idx){
		rjmcmc.bm(phy=trees[[idx]], dat=dat, SE=SE, ngen=ngen, samp=samp, type=type, filebase=files[idx], ...)
	})
}

.link.root=function(rootd, rootv, nd, d, v){
    open=nd[which(d==0)]
    rooto=rootd[rootd%in%open]
    if(length(rooto)){
        uu=unique(v[nd%in%rooto])
        if(!length(uu)==1) stop("encountered unexpected error")
    } else {
        uu=rootv
    }
    uu
}


rjmcmc.bm <- function ( phy, dat, SE=NA, ngen=50000, samp=100, type=c("jump-rbm", "rbm", "jump-bm", "bm"), ... )
{
    
#   SE=NA; ngen=50000; samp=100; type=c("jump-rbm")

    # method
	typs=c("jump-rbm", "rbm", "jump-bm", "bm")
	type=match.arg(type, typs)
	
	if(samp>ngen) stop("increase 'ngen' or decrease 'samp'")
	if("multiPhylo"%in%class(phy)) return(.rjmcmc.bm.multi(phy, dat, SE, ngen, samp, type, ...))
	
	# controller objects 
	tmp=make.gbm(phy, dat, SE=SE, type=type, ...)
	ct=tmp$control
	ct$thin=samp
	ct$ngen=ngen
	cache=tmp$cache
	sp=tmp$start
    
	lik=ct$lik
    rootd=cache$desc$fdesc[[cache$n.tip+1]]
	nd=argn(lik)$rates
	max_j=max(attr(attr(ct$dlnJUMP,"density"),"count"))
	jumpsedgewise=.jumps.edgewise(cache$phy)
	
	# runtime objects
	cur.root=sp$root
    cur.rootrate=sp$rootrate
	cur.rates=sp$rates
	cur.delta=sp$delta
	cur.nR=sum(cur.delta)
	cur.jumps=sp$jumps
	cur.nJ=sum(cur.jumps)
	cur.jumpvar=sp$jumpvar
	cur.jumprates=cur.jumps*cur.jumpvar
	cur.scalars=cur.jumprates+cur.rates
	cur.SE=sp$se
	if(is.na(cur.SE)) cur.mod=ct$beta*lik(cur.scalars, cur.root) else cur.mod=ct$beta*lik(cur.scalars, cur.root, cur.SE)
		
	tickerFreq=ceiling(ngen/30)
	
	### begin rjMCMC
    for (i in 1:ngen) {
		
		## GENERALIZED IMPLEMENTATION ##
		#		-- multiple rate categories
		#		-- point process (jumps)
		#		-- measurement error (IN PROGRESS)
		
		# initialize updates
		new.root=cur.root
        new.rootrate=cur.rootrate
		new.rates=cur.rates
		new.delta=cur.delta
		new.nR=cur.nR
		new.jumps=cur.jumps
		new.nJ=cur.nJ
		new.jumpvar=cur.jumpvar
		new.jumprates=cur.jumprates
		new.scalars=cur.scalars
		new.SE=cur.SE
		
		while(1) {
			lnLikelihoodRatio <- lnHastingsRatio <- lnPriorRatio <- 0
			cur.proposal=min(which(runif(1)<ct$prop.cs))
			
			if (cur.proposal==1) {												## update BM process
				if(runif(1) < ct$bm.jump) {                                         # adjust rate categories
					if(is.null(ct$constrainSHIFT) & runif(1) < ct$mergesplit.shift) {
						nr=.splitormerge(x=cur.rates, delta=cur.delta, control=ct, cache=cache)
						new.rates=nr$x
						new.delta=nr$delta
                        new.rootrate=.link.root(rootd, cur.rootrate, nd, new.delta, new.rates) ## JME
						new.nR=sum(new.delta)
						new.scalars=new.rates+cur.jumprates
						lnHastingsRatio=nr$lnpriorproposalRatio 
						lnPriorRatio=0
						subprop="mergesplit"
						break()	
					} else if(cur.nR>0){											# shift local rate
						nr=.adjustshift(x=cur.rates, delta=cur.delta, rootv=cur.rootrate, control=ct, cache=cache)
						new.rates=nr$new.values
						new.delta=nr$new.delta 
						new.scalars=new.rates+cur.jumprates
						lnHastingsRatio=nr$lnHastingsRatio 
						lnPriorRatio=.dlnratio(cur.scalars, new.scalars, ct$dlnRATE)
						subprop="moveshift"
						break()						
					} else {
						next()
					}
				} else {														## update JUMP process
					if(runif(1) < ct$mergesplit.shift & is.null(ct$constrainJUMP)){					
						if(cur.nJ==0) {												# adjust N jumps
							if(runif(1)<0.5 & cur.nJ<max_j) {								# -- 0 -> 1
								new.nJ=1
								tmp=.adjustjump(cur.jumps, add=TRUE, drop=FALSE, swap=FALSE, control=ct, cache=cache)
								subprop="incjump"
							} else {														# -- 0 -> 0
								tmp=list(jumps=cur.jumps, lnHastingsRatio=0)
								subprop="movejump"
							}  
						} else {
							if(runif(1)<0.5 & cur.nJ<max_j){								# -- j -> j + 1
								new.nJ=cur.nJ+1
								tmp=.adjustjump(cur.jumps, add=TRUE, drop=FALSE, swap=FALSE, control=ct, cache=cache)
								subprop="incjump"
							} else {															
								if(cur.nJ==max_j & runif(1)<0.5){							# -- max_j -> max_j
									tmp=list(jumps=cur.jumps, lnHastingsRatio=0)
									subprop="movejump"
								} else {													# -- j -> j - 1
									new.nJ=cur.nJ-1
									tmp=.adjustjump(cur.jumps, add=FALSE, drop=TRUE, swap=FALSE, control=ct, cache=cache)
									subprop="decjump"	
								}
							}
						}
					} else if(cur.nJ>0){											# adjust location of a jump
						if(cur.nJ==max_j & runif(1)<0.5){									# -- max_j -> max_j
							tmp=list(jumps=cur.jumps, lnHastingsRatio=0)
							subprop="movejump"
						} else {															# swap a jump
							tmp=.adjustjump(cur.jumps, add=FALSE, drop=FALSE, swap=TRUE, control=ct, cache=cache)
							subprop="movejump"
						}
					} else {
						next()
					}
					new.jumps=tmp$jumps
					new.jumprates=new.jumps*cur.jumpvar
					new.scalars=new.jumprates+cur.rates
					lnHastingsRatio=tmp$lnHastingsRatio
					lnPriorRatio=.dlnratio(cur.nJ, new.nJ, ct$dlnJUMP)
					lnPriorRatio=lnPriorRatio + .dlnratio(cur.scalars, new.scalars, ct$dlnRATE)
					break()
				}
			} else if(cur.proposal==2){													
				if(runif(1) < ct$tune.scale & cur.nR>0) {						## tune local rate
					nr=.tune.rate(rates=cur.rates, ct)
					new.rates=nr$values
                    new.rootrate=.link.root(rootd, cur.rootrate, nd, new.delta, new.rates)
					new.scalars=new.rates+cur.jumprates
					lnHastingsRatio=nr$lnHastingsRatio
					lnPriorRatio=.dlnratio(cur.scalars, new.scalars, ct$dlnRATE)
					subprop="ratetune"
					break()						
				} else {														## scale rates
					if(runif(1) < ct$bm.jump){										# BM rates
						sc=.proposal.slidingwindow(1, ct$prop.width, ct$rate.lim)
						new.rates=sc$v*cur.rates
                        new.rootrate=sc$v*cur.rootrate
						new.scalars=new.rates+cur.jumprates
						lnHastingsRatio=sc$lnHastingsRatio
						lnPriorRatio=.dlnratio(cur.scalars, new.scalars, ct$dlnRATE)
						subprop="ratescale"
						break()						
					} else if(cur.nJ>0){											# JUMP effectsize
						nr=.tune.rate(cur.jumpvar, ct)
						new.jumpvar=nr$values
						new.jumprates=new.jumpvar*cur.jumps
						new.scalars=cur.rates+new.jumprates
						lnHastingsRatio=nr$lnHastingsRatio
						lnPriorRatio=.dlnratio(cur.scalars, new.scalars, ct$dlnRATE)
						lnPriorRatio=lnPriorRatio + .dlnratio(cur.jumpvar, new.jumpvar, ct$dlnPULS)
						subprop="jumpvar"
						break()						
					} else {
						next()
					}
				} 
			} else if(cur.proposal==3) {										## adjust root state
				nr=.tune.value(cur.root, ct)
				new.root=nr$v
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=.dlnratio(cur.root, new.root, ct$dlnROOT)
				subprop="rootstate"
				break()
			} else if(cur.proposal==4) {										## adjust SE
				nr=.tune.SE(cur.SE, ct)
				new.SE=nr$values
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=.dlnratio(cur.SE, new.SE, ct$dlnSE)
				subprop="SE"
				break()										
			}
		}
		
        if(is.na(lnPriorRatio)) {
            stop(print(nr))
        }
        .check.root(rootd, new.rootrate, nd, new.delta, new.rates, c(subprop, nr))
        
		ct$n.subprop[subprop]=ct$n.subprop[subprop]+1				
		
		# compute fit of proposed model
		if(is.na(new.SE)) new.mod=ct$beta*lik(new.scalars, new.root) else new.mod=ct$beta*lik(new.scalars, new.root, new.SE)
		r=.proc.lnR(i, subprop, cur.mod, new.mod, lnPriorRatio, lnHastingsRatio, heat=1, control=ct)
		
		if (runif(1) <= r$r) {			## adopt proposal ##
			new.root->cur.root
            new.rootrate->cur.rootrate
			new.rates->cur.rates
			new.delta->cur.delta
			new.nR->cur.nR
			new.jumps->cur.jumps
			new.nJ->cur.nJ
			new.jumpvar->cur.jumpvar
			new.jumprates->cur.jumprates
			new.scalars->cur.scalars
			new.SE->cur.SE
			
			ct$n.subaccept[subprop] = ct$n.subaccept[subprop]+1
			
			new.mod->cur.mod
		} 
		
		# iteration-specific functions
		if(i%%tickerFreq==0 & ct$summary) {
			if(i==tickerFreq) cat("|",rep(" ",9),toupper("generations complete"),rep(" ",9),"|","\n")
			cat(". ")
		}
		if(i%%samp==0 & ct$summary) {
			pr=.gbm.prior.lik(cur.scalars, cur.jumpvar, cur.nJ, cur.nR, cur.root, ct)

			parms=list(principal=list(shifts=list(delta=cur.delta, shifts=cur.rates), jumps=list(delta=cur.jumps)), gen=i, lnL=cur.mod, 
					   lnLp=pr, qlnL.p=lnPriorRatio, qlnL.h=lnHastingsRatio, jumpvar=cur.jumpvar, SE=cur.SE, root=cur.root)
			.parlog.rjmcmc(init=FALSE, end=FALSE, parameters=parms, control=ct)
		}
	}
	
	# end rjMCMC
	.cleanup.rjmcmc(ct, cache)
}

.check.root=function(rootd, rootv, nd, d, v, ...){
    open=nd[which(d==0)]
    rooto=rootd[rootd%in%open]
    if(length(rooto)){
        uu=unique(v[nd%in%rooto])
        if(!all(uu==rootv)) {
            xx=as.list(...)
            for(i in 1:length(xx)) print(xx[[i]])
            stop("encountered unexpected error in link between 'rates' and 'root'")
        }
    }
}


