require(histogram)

##### Fonction calcNL extracted from histogram package ################
calcNL = function (y, BL, right) {
    n <- length(y)
    bn <- length(BL)
    NL <- rep(0, bn)
    if (right) {
        k <- sum(y <= BL[1])
        NL[1] <- k
        ind <- 2
        while (k < n) {
            k <- k + 1
            rbd <- BL[ind]
            while (y[k] > rbd) {
                NL[ind] <- k - 1
                ind <- ind + 1
                rbd <- BL[ind]
            }
            NL[ind] <- k
        }
    }
    else {
        k <- sum(y < BL[1])
        NL[1] <- k
        ind <- 2
        while (k < n) {
            k <- k + 1
            rbd <- BL[ind]
            while (y[k] >= rbd) {
                NL[ind] <- k - 1
                ind <- ind + 1
                rbd <- BL[ind]
            }
            NL[ind] <- k
        }
    }
    if (ind < bn) 
        NL[ind:bn] <- n
    NL
}
#############################################################################

##### Fonction histgreedy extracted from histogram package ################
histgreedy = function (BL, NL, n, binmax, verbose = verbose) {
    nB <- length(BL)
    if (nB == 2) 
        return(BL)

    inclikelihood <- function(i, j) {
        if ((i + 1) == j) 
            return(NULL)
        else if (BL[i] == BL[j]) 
            return(rep(-Inf, j - i - 1))
        else {
            k <- (i + 1):(j - 1)
            old <- (NL[j] - NL[i]) * log((NL[j] - NL[i])/(BL[j] - 
                BL[i])/n)
            new <- rep(0, length(k))
            indv <- (NL[j] - NL[k] > 0)
            new[indv] <- ((NL[j] - NL[k]) * log((NL[j] - NL[k])/(BL[j] - 
                BL[k])/n))[indv]
            indv <- (NL[k] - NL[i] > 0)
            new[indv] <- new[indv] + ((NL[k] - NL[i]) * log((NL[k] - 
                NL[i])/(BL[k] - BL[i])/n))[indv] - old
            new
        }
    }

    breaks <- c(1, nB)
    increment <- c(-Inf, inclikelihood(breaks[1], breaks[2]), 
        -Inf)
    while ((max(increment) > 0) && (length(breaks) < (binmax + 
        1))) {
        maxi <- which.max(increment)
        here <- sum(breaks < maxi)
        i <- breaks[here]
        j <- breaks[here + 1]
        debut <- increment[1:i]
        debut[i] = -Inf
        fin <- increment[j:nB]
        fin[1] = -Inf
        gche <- inclikelihood(i, maxi)
        drte <- inclikelihood(maxi, j)
        increment <- c(debut, gche, -Inf, drte, fin)
        breaks <- c(breaks[1:here], maxi, breaks[(here + 1):length(breaks)])
    }
    breaks <- BL[breaks]
    if (verbose) 
        message(paste("- Pre-selected finest partition with", 
            length(breaks) - 1, "bins."))
    return(breaks)
}
#############################################################################

##### Fonction DynamicExtreme extracted from histogram package ################
DynamicExtreme = function (FctWeight, n, Dmax, mini = TRUE, msg = TRUE) {
	
    OptimizePath <- function(D) {
        ancestor0 = matrix(nrow = n - D + 1)
        cumWeight0 = matrix(nrow = n - D + 1)
        for (i in D:n) {
            tmp <- cumWeight[(D - 1):(i - 1), D - 1] + weight[D:i, 
                i + 1]
            if (mini) 
                ancestor0[i - D + 1] <- which.min(tmp)
            else ancestor0[i - D + 1] <- which.max(tmp)
            cumWeight0[i - D + 1] <- tmp[ancestor0[i - D + 1]]
        }
        ancestor[D:n, D] <<- ancestor0 + (D - 2)
        cumWeight[D:n, D] <<- cumWeight0
    }
    
    if (n == 1) 
        return(list(extreme = FctWeight(1, 2), ancestor = cbind(0, 
            0)))
    cumWeight <- matrix(nrow = n, ncol = Dmax)
    ancestor <- matrix(0, nrow = n, ncol = Dmax)
    weight <- matrix(nrow = n + 1, ncol = n + 1)
    if (msg) 
        message("- Computing weights for dynamic programming algorithm.")
    for (i in 1:n) weight[i, (i + 1):(n + 1)] = mapply(FctWeight, 
        rep(i, n - i + 1), (i + 1):(n + 1))
    cumWeight[, 1] <- weight[1, 2:(n + 1)]
    if (msg) 
        message("- Now performing dynamic optimization.")
    mapply(OptimizePath, 2:Dmax)
    extreme <- cumWeight[n, ]
    list(extreme = extreme, ancestor = ancestor)
}
#############################################################################

##### Fonction DynamicList adapted from histogram package ################
DynamicList = function (C, B, D) {
	
	PathList = function (C, D) {
	    L <- nrow(C)
	    for (i in D:1) {
	        L <- c(C[L[1], i], L)
	    }
	    return(L)
	}
    
    L <- PathList(C, D) + 1
    bounds <- B[L]
    return(bounds)
}
#############################################################################



################ Utilities for Testimation simulations ##################
################ Utilities for Testimation simulations ##################
################ Utilities for Testimation simulations ##################
################ Utilities for Testimation simulations ##################
TBuildEstim = function(descript,X) {
	f = switch(descript[[1]],					
		'kernel' = Tkernel(X=X,bw=descript$bw,kernel=descript$kernel),
		'histo' = Thisto(X=X,D=descript$D,breaks=descript$breaks),
		'parametric' = Tparametric(X=X,name=descript$name)
		)
	
}

Tkernel = function(X,bw=NULL,kernel) {
	# X is a sample
	# bw bandwidth
	# kernel specifies the kernel type as in the procedure density
	# kernel can be any kernel of kde R-function, see density
	if (is.null(bw)) stop('Bandwidth should be provided')
	R = range(X)
	dX = diff(R)
	# cuts = c(R[1]-dX/20,R[2]+dX/20)
	# f=density(X,bw=bw,kernel=kernel,n=1024,from=cuts[1],to=cuts[2])
	f=density(X,bw=bw,kernel=kernel,n=1024)
	cuts = c(f$x[1],f$x[1024])
	# I = sum((f$y[-1]+f$y[-length(f$y)])/2*diff(f$x)); f$y = f$y/I
	list(f=approxfun(f$x,f$y,yleft=0,yright=0),cuts=cuts,descript=list('kernel',kernel=kernel,bw=bw))
}

Thisto = function(H,X=NULL,D=NULL,breaks=NULL){
	# H is the result of hist procedure
	# X is a sample
	# D is the number of bin of a regular histogram
	# breaks specifies the bin breaks of the histogram
	if (!is.null(X)) 
		if (is.null(D)&&is.null(breaks)) 
			stop('Bin number or breaks should be provided')
		else {
			if (is.null(breaks)) { # regular case defined by the bin number
				R = range(X)
				breaks=seq(from=R[1],to=R[2],l=D+1)
				breaks[D+1] = breaks[D+1]+1e-9
			} else {
				D = length(breaks)-1
				if (min(X)<breaks[1]) breaks[1]=min(X)-1e-9
				if (breaks[D+1]<max(X)) breaks[D+1]=max(X)+1e-9
			}
			H=hist(X,breaks=breaks,plot=F)
			}
	D = length(H$breaks)-1		
	histo = function(z) {
		dz = sapply(z,function(y)  d=sum(H$breaks<=y)) # compute bin index
		y = rep(0,l=length(z))
		ind = which((1<=dz)&(dz<=D))
		y[ind] = H$density[dz[ind]]
		y
	}
	list(f=histo,cuts=H$breaks,descript=list('histo',D=D,breaks=breaks))
}

TBuildKernel=function(X,n=length(X),bwtab=NULL,kerneltab=NULL) {
	# kernel can be any kernel of density R-function, see density
	if (is.null(kerneltab)) kerneltab="epanechnikov"
	if (is.null(bwtab))  bwtab=diff(range(X))/2/(ceiling(n/log(n)):1)
	flist = list()
	i = 0
	for(kernel in kerneltab){
		for (bw in bwtab) {
			i = i+1
			flist[[i]] = Tkernel(X=X,bw=bw,kernel=kernel)
			}
		}
	flist
}

TBuildRegularHisto=function(X,n=length(X),Dmax=NULL,Dtab=NULL) {
	
	if (is.null(Dmax)) Dmax=ceiling(n/log(n))
	if (is.null(Dtab)) Dtab=1:Dmax
	flist = list()
	i = 0
	for (D in Dtab) {
		i = i+1
		flist[[i]] = Thisto(X=X,D=D)
		flist[[i]]$descript = list('histo',D=D,breaks=NULL)
		}
	flist
}

TBuildIrregularHisto=function(X,n=length(X),Dmax=NULL,
	greedyfirst=TRUE,grid=c("data","regular","quantiles"),breaks=NULL,verbose=FALSE){
		
	grid = grid[1]
	epsilon = 1e-8
	a = min(X)
  	b = max(X)+epsilon
  	y = sort((X-a)/(b-a)) # normalize range to [0,1]
	flist = list()
  	
	if ((grid=="regular")|(grid=="quantiles")) {
		if (is.null(Dmax)|(Dmax<0)) {
			grid='data'; greedyfirst=TRUE; 
			warning('Dmax not defined, use greedy procedure first')
			}
		BL = (0:Dmax)/Dmax
		if (verbose) 
			message(paste("- Using a regular grid with ",length(BL)-1," bins as finest grid.",sep=""))
  		}
	if (grid=="quantiles") {
		BL = quantile(y,probs=BL)
		BL[Dmax+1] = BL[Dmax+1]+epsilon
		if (verbose) 
			message("- Using a regular quantile grid with ",length(BL)-1," bins as finest grid.",sep="")
		}
	if (grid=="data") {
		BL = c(y[1],y[-n]+diff(y)/2,y[n]+epsilon)
  		BL = unique(sort(BL))
		if (verbose) 
			message("- Using finest grid based on observations.")
  		}
  		
	NL=calcNL(y,BL,right=F)
	
	# Greedy procedure
	if (greedyfirst==TRUE) {
		D = length(BL)-1 
		binmax<-ceiling(min(max(100,D^(1/3)),D))
		if ((D>=2)&(D>binmax)){
			if (verbose)
				message(paste("- Using greedy procedure to recursively build",
					" a finest partition with at most ",binmax," bins.",sep=""))
    		BL=histgreedy(BL,NL,n,binmax,verbose=verbose)
    		NL=calcNL(y,BL,right=F)
			}
 		}
	Dmax = length(BL)-1
	
	weightfunction<-function(i,j) {
		dN = NL[j]-NL[i]
		dBL = BL[j]-BL[i]
		# compute the loglikelihood part of interval [B(i),B(j)[
		if (dN==0)
			W<-0
		else
			if ((dBL<log(n)^1.5/n)&(grid=="data"))
				return(-Inf) # remove small bins
			else
				W <- dN*log(dN/dBL)
		if (dN<log(n)) W<- -Inf
		return(W)
		}
	tmp = DynamicExtreme(weightfunction,n=length(BL)-1,Dmax=Dmax,mini=FALSE,msg=verbose)
	i = 0
	for (D in 1:Dmax) {
		i = i+1
		breaks = DynamicList(tmp$ancestor,BL,D)
		flist[[i]] = Thisto(X=X,breaks=a+(b-a)*breaks)
		flist[[i]]$descript = list('histo',D=D,breaks=a+(b-a)*breaks)
	}
	flist
}

Tparametric = function(X,name){
	#c("uniform","normal","lognormal","beta","exp","gamma","chisq")) {
	# X is a sample
	# name is the R-name of the distribution: 'norm', 'exp', 'gamma', 'chisq', ...
	# returns the function related to the mle or moment estimate of the parameter
	eps = 1e-6
	
	n = length(X)
	R = range(X); mini = R[1]; maxi = R[2]
	
	if (name =='unif') {
		par = list(min=mini-eps, max=maxi+eps)
		f=function(x){dunif(x,min=par[[1]],max=par[[2]])}
		cuts = qunif(p=c(eps,1-eps),min=par[[1]],max=par[[2]])
		# moment method: list(f=function(x){dunif(x,min=0,max=2*mean(X))},cuts=cuts)
		}
	if (name =='norm') {
		par = list(mean=mean(X),sd=sd(X))
		f=function(x){dnorm(x,mean=par[[1]],sd=par[[2]])}
		cuts = qnorm(p=c(eps,1-eps),mean=par[[1]],sd=par[[2]])
		}
	if (name =='lnorm') {
		if (mini<=0) return(NA)
		# méthode des moments
		m1 = mean(X)
		s2 = var(X)
		sdlog = log(1+s2/m1^2)
		meanlog = log(m1)-sdlog/2
		par = list(meanlog=meanlog,sdlog=sdlog)
		if (par$sdlog<=0) return(NA)
		f=function(x){dlnorm(x,meanlog=par[[1]],sdlog=par[[2]])}
		cuts = qlnorm(p=c(eps,1-eps),meanlog=par[[1]],sdlog=par[[2]])
		}
	if (name =='exp') {
		if (mini<=0) return(NA)
		m = mean(X)
		par = list(rate=1/m)
		if (par$rate<=0) return(NA)
		f=function(x){dexp(x,rate=par[[1]])}
		cuts = qexp(p=c(eps,1-eps),rate=par[[1]])
		}
	if (name =='gamma') {
		if (mini<=0) return(NA)
		# méthode des moments
		m1 = mean(X)
		s2 = var(X)
		par = list(shape=m1^2/s2,scale=s2/m1)
		if ((par$shape<=0)|(par$scale<=0)) return(NA)
		f=function(x){dgamma(x,shape=par[[1]],scale=par[[2]])}
		cuts = qgamma(p=c(eps,1-eps),shape=par[[1]],scale=par[[2]])
		}
	if (name =='chisq') {
		if (mini<=0) return(NA)
		df=max(1,round(mean(X)))
		par = list(df=df)
		f=function(x){dchisq(x,df=par[[1]])}
		cuts = qchisq(p=c(eps,1-eps),df=par[[1]])
		}
	if (name =='beta') {
		# méthode des moments
		if ((mini<=0)|(maxi>=1)) return(NA)
		m1=mean(X)
		s2=var(X)
		tmp = (m1*(1-m1)/s2-1)
		par = list(shape1=m1*tmp,shape2=(1-m1)*tmp)
		f=function(x){dbeta(x,shape1=par[[1]],shape2=par[[2]])}
		cuts = qbeta(p=c(eps,1-eps),shape1=par[[1]],shape2=par[[2]])
		}
	list(f=f,cuts=cuts,descript=list('parametric',name=name,par=par))
}

TBuildParametric=function(X,namelist=NULL) {
	
	if (is.null(namelist)) namelist=c('unif','norm','exp','lnorm','gamma','chisq','beta') # ,'exp-tr','gamma-tr','chisq-tr')
	
	flist = list()
	i = 0
	for (name in namelist) {
		i = i+1
		flist[[i]] = Tparametric(X=X,name=name)
		}
	flist = flist[!is.na(flist)]
}

TBuildList = function(X,family=c('Kernel','RegularHisto','IrregularHisto','Parametric'),
	kerneltab=NULL,bwtab=NULL,Dmax=NULL,Dtab=NULL) {
		
	flist = NULL;
	for (i in 1:length(family)) {
		tmplist = switch(family[i],
			'Kernel' = TBuildKernel(X,kerneltab=kerneltab,bwtab=bwtab),
			'IrregularHisto' = TBuildIrregularHisto(X,Dmax=Dmax),
			'RegularHisto' = TBuildRegularHisto(X,Dmax=Dmax,Dtab=Dtab),
			'Parametric' = TBuildParametric(X)
			)
		flist = c(flist,tmplist)
	}
	flist
}

Tnorm = function(g1,p) {
	g2 = g1
	g2$f = function(x) {0*x}
	g2$cuts = NULL
	Tdistance(g1,g2,p)
}

Tdistance = function(g1,g2,p,npts=100) {
	# compute Lp distance between g1 and g2
	# if p==0 uses Hellinger distance
	
	p[tolower(p)=='hellinger'] = 0
	p = unlist(p)
			
	phi = function(a,b,p,nd) {
		
		if (is.histo) { # two histograms 
			I = rep(0,length(p))
			for (i in 1:length(p))
				if (p[i]==0) I[i] = (sqrt(g2$f(a))-sqrt(g1$f(a)))^2
				else I[i] = (abs(g2$f(a)-g1$f(a)))^p[i]
			return((b-a)*I)
		}
		
		######### phi INTERNAL FUNCTION #######
		phi.internal2 = function(a,b,p,nd){				
			I = rep(0,length(p))
			for (i in 1:length(p)) 
				if (p[i]==0)
					I[i] = integrate(function(x) (sqrt(g1$f(x))-sqrt(g2$f(x)))^2,a,b,stop.on.error=F,subdivisions=nd)$value
				else 
					I[i] = integrate(function(x) abs(g1$f(x)-g2$f(x))^p[i],a,b,stop.on.error=F,subdivisions=nd)$value
			I
		}
		######### END phi.rec INTERNAL FUNCTION #######		
		return(phi.internal2(a,b,p,nd))
	}
		
	cuts = sort(unique(c(g1$cuts,g2$cuts)))
	
	is.histo = ((exists('descript',g1)&&(g1$descript[[1]]=='histo'))
		&(exists('descript',g2)&&(g2$descript[[1]]=='histo'))) 
	
	L0 = diff(range(cuts))
	nbins = length(cuts)-1

	dist = rep(0,length(p))
	for (d in 1:nbins) {
		if (cuts[d]+1e-9>=cuts[d+1]-1e-9) next
		nd = ceiling(npts*(cuts[d+1]-cuts[d])/L0)
		I = phi(cuts[d]+1e-9,cuts[d+1]-1e-9,p,nd)
		dist = dist+I
		}	
	
	for (i in 1:length(p))
		if (p[i]==0) dist[i] = sqrt(dist[i]/2)
		else dist[i] = dist[i]^(1/p[i])

	dist[p==0] = min(1,dist[p==0]) # ensures that Hellinger distance is less than 1
	dist
	}

###################### MAIN FUNCTION #########################
###################### MAIN FUNCTION #########################
###################### MAIN FUNCTION #########################
###################### MAIN FUNCTION #########################

DensityTestim = function(X,p=1/2,family=NULL,test=c('birge','baraud'),theta=1/4,last=c('full','training'),
	plot=TRUE,verbose=TRUE,wlegend=TRUE,kerneltab=NULL,Dmax=NULL,bwtab=NULL,
	do.MLHO=FALSE,do.LSHO=FALSE,start=c('LSHO','MLHO'),csqrt=1,
	H2dist=NULL,allImageX2=NULL,flist=NULL,
	...) {
					
	test = match.arg(test)		
	last = match.arg(last)
	if (verbose) print(paste('Uses test',test))		
					
	if (is.null(family)) 
		family=c('Kernel','RegularHisto','IrregularHisto','Parametric')
		
		
	################# INTERNAL PROCEDURE #################		
	################# INTERNAL PROCEDURE #################		
	################# INTERNAL PROCEDURE #################		
	################# INTERNAL PROCEDURE #################		
	HO = function() {
				
		# ML stands for Maximum Likelihood
		# MLE for Maximum Likehood Estimate
		
		Thellinger2 = function(i,j) Tdistance(flist[[i]],flist[[j]],'Hellinger')^2
		
		BuildAllImage = function() {
			# Compute the image of X2 for each estimate of flist
			lik = sapply(1:length(flist),
				function(m) allImageX2[m,] <<- flist[[m]]$f(X2)
				)
		}
	
		RobustTest <- function(i,j,affinity=NULL){ 
			# internal function
			# Compute the robust test between estimate i and j in flist
			# Return TRUE if j is chosen and FALSE otherwise
			# This procedure use balls defined by an affinity and a radius depending on theta 
			# Default value of theta 1/4
			# Default affinity (=NULL) is Hellinger affinity
			
			if (i==j) return(FALSE)
			
			gi=flist[[i]]
			gj=flist[[j]]
			
			if (is.na(allImageX2[i,1])) allImageX2[i,] <<- gi(X2)
			if (is.na(allImageX2[j,1])) allImageX2[j,] <<- gj(X2)
			
			Yi = allImageX2[i,]
			Yj = allImageX2[j,]
			
			ind = which(Yi+Yj>0)
			Yi = Yi[ind]
			Yj = Yj[ind]
			sYi = sqrt(Yi)
			sYj = sqrt(Yj)
			
			###### Birgé (2006)'s test
			if (test=='birge') {
				if (is.na(affinity)||(affinity==1)) return(FALSE)
				omega<-acos(affinity)
				P = prod(
					(sin(omega*theta)*sYi+sin(omega*(1-theta))*sYj) / 
					(sin(omega*theta)*sYj+sin(omega*(1-theta))*sYi)
					) - 1
			}			
			
			###### Baraud (2011)'s test
			if (test=='baraud') {
				sij = gi ##### to be the mean of i and j
				sij$f = function(z) (gi$f(z)+gj$f(z))/2 ### now it is the mean
				sij$cuts = sort(unique(c(gi$cuts,gj$cuts)))
				if (!((exists('descript',gi)&&(gi$descript[[1]]=='histo'))
					&(exists('descript',gj)&&(gj$descript[[1]]=='histo')))) sij$descript = NULL
				P = Tdistance(gi,sij,'Hellinger')^2-Tdistance(gj,sij,'Hellinger')^2 + sqrt(2)*sum((sYj-sYi)/sqrt((Yj+Yi)))/n2
			}
			
			# if (P>=0) j else i
			H2test[i,j] <<- (P>=0)
			H2test[j,i] <<- (!H2test[i,j])
			(P>=0)
		}
	
		M = length(flist); 
		thr2 = csqrt^2/n2
		
		# Matrix of squared Hellinger distances
		if (is.null(H2dist)) H2dist<-matrix(NA,nrow=M,ncol=M)
		
		# Matrix of the Hellinger test result
		H2test<-matrix(NA,nrow=M,ncol=M)
	
		# Value at X of the flist estimates	
		if (is.null(allImageX2)) {
			allImageX2 <<- matrix(NA,nrow=M,ncol=n2)
			BuildAllImage()
			} 
		
		if ((start[1]=='MLHO')|(do.MLHO)) {
			ML = apply(allImageX2,1,function(x) sum(log(x[x>0])))/n2
			i0 = i0.ML = which.max(ML)
			MLHO = flist[[i0.ML]]
			do.MLHO = T
			} else i0.ML = MLHO = NULL
	
		if ((start[1]=='LSHO')|(do.LSHO)) {
			LS=-2*rowMeans(allImageX2)
			for(m in 1:M) LS[m] = LS[m]+Tnorm(flist[[m]],2)^2			
			i0 = i0.LS = which.min(LS)
			LSHO=flist[[i0.LS]]
			do.LSHO = T
			} else i0.LS = LSHO = NULL
			
		if (is.numeric(start[1])) i0=start[1]
			
		# Compute the distance between i0-select estimate and the others
		# Provide the radius G2 of i0-neighborhood
		# Here, all distances to i0 need to be computed as no distance have been computed yet
		G2 = 0;
		for (i in (1:M)[-i0]){
			if (is.na(H2dist[i0,i])) H2dist[i0,i]=H2dist[i,i0]=Thellinger2(i0,i)
			if (RobustTest(i0,i,affinity=1-H2dist[i,i0])) G2<-max(G2,H2dist[i0,i])
		}
		
		# Potential (further) good points : J<-J[-1] are in the previous neighborhood
		# As H2dist[i0,i0] is NA, i0 should not be selected
		# We consider only points in the ]corona 2/sqrt(n), G2]
		J = which((thr2<H2dist[i0,])&(H2dist[i0,]<=G2)) # we use only points far enough from i0
		
		# Try to find a better point (= with smaller neighborhood of better guys)
		while(length(J)>0){
			# order the potential points with their distance to the actual "best" guy
			J = J[order(H2dist[i0,J],decreasing=T)]
		
			# we start from the closest points to the actual "best" guy
			j = J[1] # new current point for computing G2
			J = J[-1] # the other active points
			
			# the set of tested points around location j so-called j-tested points
			jSet = j
				
			# compute the radius of the neighborhood of j "the competitor"	
			G2tmp = 0
			for(i in (1:M)[-j]){
				# check whether i is at distance not larger than csqrt/sqrt(n) from a j-tested point
				if (any(H2dist[i,jSet]<=thr2,na.rm=T)) next 
				# compute distance from i to j
				if (is.na(H2dist[i,j])) {
					H2dist[j,i]=H2dist[i,j]=Thellinger2(j,i)
					# check whether i is at distance not larger than csqrt/sqrt(n) from j
					if (H2dist[i,j]<=thr2) next
					}
				jSet = c(jSet,i)
				# do the robust test
				if (RobustTest(j,i,affinity=1-H2dist[i,j])) 
					if (H2dist[i,j]>G2tmp) {
						G2tmp = H2dist[i,j]
						if (G2tmp>=G2) break # too far ... bad j ... keep the actual "best" guy 
					}
				}	
		
			if (G2tmp<G2){
				# we found a new "best" guy
				G2 = G2tmp; i0 = j;
				J = J[which((thr2<H2dist[i0,J])&(H2dist[i0,J]<=G2))]  		### OK <- new version
				}	
			}
			
		THO = flist[[i0]]
		
		if (last[1]=='full') THO = TBuildEstim(THO$descript,X)

		if ((last[1]=='full')&(do.MLHO)) {
			if (i0==i0.ML) MLHO = THO
			else MLHO = TBuildEstim(MLHO$descript,X)
			}
		
		if ((last[1]=='full')&(do.LSHO)) {
			if (i0==i0.LS) LSHO = THO
			else LSHO = TBuildEstim(LSHO$descript,X)
			}	
	
		if (plot) {
			x = seq(from=min(X),to=max(X),l=500)
			y = THO$f(x)
			lwd = 2; lty = 1; col = 'red'; leg = 'T'
			if (do.MLHO) {
				y = cbind(y,MLHO$f(x))
				lwd = c(lwd,1); lty = c(lty,1); col = c(col,'black'); leg = c(leg,'ML')
			}
			if (do.LSHO) {
				y = cbind(y,LSHO$f(x))
				lwd = c(lwd,1); lty = c(lty,2); col = c(col,'black'); leg = c(leg,'LS')
			}
			matplot(x,y,col=col,lwd=lwd,lty=lty,type='l',ylab='',...)
			#title('T estimate')
			if (wlegend) legend("topright",legend=leg,col=col,lwd=lwd,lty=lty)
			}
			
		if (verbose) print(THO$descript)
		if (verbose&do.MLHO) print(MLHO$descript)
		if (verbose&do.LSHO) print(LSHO$descript)
	
		list(THO=THO,MLHO=MLHO,LSHO=LSHO,M=M,comput=sum(!is.na(H2test))/2,total=M*(M-1)/2,
			iTHO=i0,iMLHO=i0.ML,iLSHO=i0.LS,H2dist=H2dist,allImageX2=allImageX2,flist=flist)
	
	} 
	################# END INTERNAL PROCEDURE #################		
	################# END INTERNAL PROCEDURE #################		
	################# END INTERNAL PROCEDURE #################		
	################# END INTERNAL PROCEDURE #################		
	
	n = length(X)
	n1 = ceiling(p*n); n2 = n-n1;
	# X1 training sample
	# X2 testing sample used to select the T-estimate from flist (see below)
	X1 = X[1:n1]; X2 = X[(n1+1):n]
	
	if (is.function(csqrt)) csqrt=csqrt(n2)

	# list of density estimates defined as functions	
	if (is.null(flist)) flist = TBuildList(X1,family=family,kerneltab=kerneltab,bwtab=bwtab,Dmax=Dmax)
	
	######## call the main procedure #############	
  	HO()
}

