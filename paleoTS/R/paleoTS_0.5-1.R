as.paleoTS<- function (mm, vv, nn, tt, MM = NULL, genpars = NULL, label = NULL, 
    start.age = NULL, oldest = c("first", "last"), reset.time = TRUE) 
{
    oldest<- match.arg(oldest)
    y <- list(mm = mm, vv = vv, nn = nn, tt = tt, MM = MM, genpars = genpars, 
        label = label, start.age = start.age)
    if (oldest == "last") {
        oo <- length(y$mm):1
        y$mm <- y$mm[oo]
        y$vv <- y$vv[oo]
        y$nn <- y$nn[oo]
        y$tt <- y$tt[oo]
    }
    if(y$tt[1] > y$tt[2])	timeDir<- "decreasing"
    else 					timeDir<- "increasing"
    
    if (reset.time) {
        if (y$tt[1] != 0) {
            sa <- y$tt[1]
            if (!is.null(y$start.age) && sa != y$start.age) 
                stop("Age of first sample does not match start.age argument")
            if(timeDir == "decreasing")	y$tt <- sa - y$tt   # decreasing ages (e.g., Ma)  
            else					    y$tt <- y$tt - sa   # increasing ages (elapsed time)
            y$start.age<- sa	
        }
    }
    y$timeDir<- timeDir
    class(y) <- "paleoTS"
    return(y)
}



as.paleoTSfit<- function(logL, parameters, modelName, method, K, n, se)
{
  ic<- IC(logL=logL, K=K, n=n, method='AICc')
  y<- list(logL=logL, AICc=ic, parameters=parameters, modelName=modelName, method=method, se=se, K=K, n=n)
  class(y)<- "paleoTSfit"
  return(y)	
}



read.paleoTS<- function (file = NULL, oldest = "first", reset.time = TRUE, ...) 
{
    if (is.null(file)) {
        ff <- file.choose()
        x <- read.table(ff, ...)
        lab1 <- basename(ff)
    }
    else {
        x <- read.table(file = file, ...)
        #lab1 <- paste(getwd(), file)
        lab1<- file
    }
    xr <- as.paleoTS(mm = x[, 2], vv = x[, 3], nn = x[, 1], tt = x[, 
        4], label = lab1, oldest = oldest, reset.time = reset.time)
    return(xr)
}


modelCurves<- function(x, w, np=500)
# returns list of model means, upper and lower 95% probability envelopes
{
  ee<- ii<- array(dim=np)  # set up arrays
  mn<- w$modelName
  mp<- w$par
  x0<- ifelse(w$method=="AD", x$mm[1], w$par["anc"])
  ttp<- seq(x$tt[1], x$tt[length(x$tt)], length.out=np)
  #ttp<- x$tt
  #ttp<- sort(c(seq(x$tt[1], x$tt[length(x$tt)], length.out=np), x$tt))
  #print(ttp)
  
  # list of models implemented
  okModels<- c("URW", "GRW", "Stasis", "OU", paste("Punc-", 1:100, sep=""))
  if (mn %in% okModels){  
	  if(mn=="URW"){ ee<- rep(x0, np); vv<- mp["vstep"]*ttp}  
	  if(mn=="GRW"){ ee<- x0+mp["mstep"]*ttp; vv<- mp["vstep"]*ttp}
	  if(mn=="Stasis"){ ee<- rep(mp["theta"], np); vv<- rep(mp["omega"],np) }
	  if(mn=="OU"){ if(!is.null(x$start.age)) tto<- x$start.age-ttp
	  				else	tto<- ttp
	  				ee<- mp["theta"] * (1-exp(-mp["alpha"]*tto)) + mp["anc"]*exp(-mp["alpha"]*tto)
	  				vv<- (mp["vstep"]/(2*mp["alpha"])) * (1-exp(-2*mp["alpha"]*tto)) }
	  if(grepl("Punc", mn)){
	  		# handle time shifts
	  		sp<- w$par[grep("shift", names(w$par))]  # shift samples
	  		st<- x$tt[sp]  # shift times
	  		ng<- length(sp)+1
	  		ggt<- cut(ttp, breaks=c(min(x$tt)-1, st, max(x$tt)+1), right=TRUE, labels=1:ng)
	
	
	  		# extract needed paramaters
	  		th<- w$par[grep("theta", names(w$par))]  # theta estimates
	  		om<- w$par[grep("omega", names(w$par))]  # omega estimates
	  		if(length(om)==1)	om<- rep(om, ng)	 # make into vector if needed
			
			# get ee and vv
			ee<- th[ggt]
			vv<- om[ggt]  		
		}  
	} else {ee<-NA; vv<- NA; warning(paste("modelFit argument not implemented for model", mn, "\n"))}
	
	
  if (!is.null(x$start.age))	tto<- x$start.age - ttp  else tto<- ttp 
  res<- list(tt=ttp, ee=ee, ll=ee-1.96*sqrt(vv), uu=ee+1.96*sqrt(vv))	
  return(res)	
}


plot.paleoTS<- function (x, nse = 1, pool = FALSE, add = FALSE, modelFit = NULL, 
    pch = 19, lwd = 1.5, ylim=NULL, ...) 
{
    if (pool) 
        x <- pool.var(x, ret.paleoTS = TRUE)
    se <- sqrt(x$vv/x$nn)
    lci <- x$mm - (nse * se)
    uci <- x$mm + (nse * se)
    xx <- x
    if (!is.null(x$start.age)) {
        if(x$timeDir=="decreasing") { x$tt <- x$start.age - x$tt; xl <- rev(range(x$tt))}
        else						{ x$tt<- x$tt + x$start.age; xl<- range(x$tt)}
    }
    else xl <- range(x$tt)
    if (!is.null(modelFit)) {
        mlab <- paste(modelFit$modelName, "expectation [95% prob. interval]")
        mc <- modelCurves(xx, w = modelFit)
        if (is.na(mc$ee[1])) 
            modelFit <- NULL
    }
    if (is.null(modelFit)) 
        yl <- c(uci, lci)
    else yl <- c(uci, lci, mc$ll, mc$uu)
    if(is.null(ylim)) ylim<- range(yl)
    if (!add) 
        plot(range(x$tt), ylim=ylim, typ = "n", pch = 19, 
            xlab = "Time", ylab = "Trait Mean", xlim = xl, ...)
    if (!is.null(modelFit)) {
        if (!is.null(x$start.age)) 
            mc$tt <- x$start.age - mc$tt
        polygon(c(mc$tt, rev(mc$tt)), c(mc$uu, rev(mc$ll)), col = "wheat2", 
            border = "white")
        lines(mc$tt, mc$ee, col = "tan", lwd = 2)
    }
    lines(x$tt, x$mm, lwd = lwd, ...)
    segments(x$tt, lci, x$tt, uci, lty = 1, lwd = lwd, ...)
    points(x$tt, x$mm, pch = pch, cex = 1.2, ...)
    mtext(x$label, cex = 0.7, col = "grey", font = 3)
    if (!is.null(modelFit)) 
        mtext(mlab, side = 4, cex = 0.8, col = "tan", font = 2)
}




akaike.wts<- function(aa)
## aa is a vector of AIC or AICc values
{
  	# get subset, dropping NAs
  	okset<- !is.na(aa)
    aas<- aa[okset]

    ma<- min(aas)
    delt<- aas - ma
    denom<- sum(exp(-delt/2))
    ww<- exp(-delt/2)/denom
    names(ww)<- names(aa)
    
    aw<- ww[okset]
    return(aw)	
}



IC<- function(logL, K, n=NULL, method=c("AICc", "AIC", "BIC"))
# compute Information Criteria from log-likelihood, # parameters (K), and 
# sample size (n), if needed.
{
 method<-match.arg(method)	
 if ((method=="AICc" || method=="BIC") && is.null(n))  stop('AICc requires n.')
 if (method=="AIC")	ic<- -2*logL + 2*K
 if (method=="AICc")	ic<- -2*logL + 2*K + (2*K*(K+1))/(n-K-1)
 if (method=="BIC")	ic<- -2*logL + K*log(n)
 
 return(ic)
}



pool.var<- function (y, nn = NULL, minN = NULL, ret.paleoTS = FALSE) 
{
    if (class(y) == "paleoTS") {
        if (all(y$nn == 1)) 	vp <- mean(y$vv)
        else vp <- sum(y$vv * (y$nn - 1))/sum(y$nn - 1)
        #else vp <- sum(y$vv * (y$nn - 1))/sum(y$nn - 1
    }
    else vp <- sum(y * (nn - 1)/sum(nn - 1))

	# replace all vv's, or just those with nn<minN
    if (ret.paleoTS) {
        yn <- y
        if(is.null(minN))  yn$vv <- rep(vp, length(y$mm)) else yn$vv[yn$nn<minN]<- vp
        return(yn)
    }
    else {
        return(vp)
    }
}

test.var.het<- function (y, method="Bartlett")
# test for variance heterogeneity among samples in a paleoTS object
{

 vp<- pool.var(y)
 NN<- sum(y$nn)
 k<- length(y$vv)
 top<- (NN-k)*log(vp)-(sum((y$nn-1)*log(y$vv)))
 bot<- 1+ (1/(3*(k-1)))*( (sum(1/(y$nn-1))) - (1/(NN-k)))
 TT<- top/bot
 p.val<- pchisq(TT, df=k-1, lower.tail=FALSE)

 w<-list(stat=TT, p.value=p.val, df=k-1)
 return (w)
}

ln.paleoTS <- function (y)
# returns paleoTS, with data approx ln-transformed
# mean(ln[y])= ln(mean[y]); var(ln[y])=(sd[y]/mean[y])^2
{
 logx<- y
 logx$mm<- log(y$mm)
 logx$vv<- (sqrt(y$vv)/y$mm)^2
 
 return (logx)
}

std.paleoTS <- function (y, zero="start")
# returns paleoTS, converted in phenotypic SD units
# mm -> (mm - mean(mm) )/sqrt(vp); vv -> vv/vp
# optionally set starting mean value to zero
{
 vp<- pool.var(y)
 sp<- sqrt(vp)
 
 ys <- y 
 ys$mm<- (y$mm - mean(y$mm)) /sp
 ys$vv<- y$vv/vp
 
 if (zero=="start")
 	ys$mm <- ys$mm - ys$mm[1]
 
 return(ys)	
}


sub.paleoTS <- function (y, ok=NULL, k=0.1)
# subsample a paleoTS, either from steps given by T/F vector 'ok'
#  proportion 'k' of samples, chosen randomly
{
 ys<- y
 ns<- length(y$mm)
 take<- array(FALSE, dim=ns)
 if (!is.null(ok) )
   take<- ok
 else
   take[sample(1:ns, size=round(k*ns))]<- TRUE

 ys$mm<- y$mm[take]
 ys$vv<- y$vv[take]
 ys$nn<- y$nn[take]
 ys$tt<- y$tt[take]
 ys$MM<- y$MM[take]
 ys$label<- paste ("Subsetted from--", y$label)

 return(ys)
}




compareModels<- function(..., silent=FALSE)
{
  modelList<- list(...)
  
  # make sure all are paleoTSfit objects, and all use same method (AD or Joint)
  classv<- sapply(modelList, FUN=class)
  methv<- sapply(modelList, FUN=function(x){x$method})
  nv<- sapply(modelList, FUN=function(x){x$n})
  nm<- length(modelList)
  
  if(!all(classv=='paleoTSfit'))  	stop("All objects must be of class 'paleoTSfit'")
  if(!all(methv==methv[1]))			stop(paste("All model fits must use the same method (AD or Joint)", sep='\n'))
  else method<- methv[1]
  if(!all(nv==nv[1]))				stop("Objects have differing n.")
  else nn<- nv[1]
 
  
  # construct data frame and parameter list
  logL<- sapply(modelList, FUN=function(x){x$logL})
  K<- sapply(modelList, FUN=function(x){x$K})
  AICc<- sapply(modelList, FUN=function(x){x$AICc})
  Akaike.wt<- round(akaike.wts(AICc),3)
  df<- data.frame(logL, K, AICc, Akaike.wt)
  rn<- sapply(modelList, FUN=function(x){x$modelName})  # extract model names
  rn<- make.unique(rn)  # handles situation if there are non-unique model names
  row.names(df)<- rn
  
  pl<- lapply(modelList, FUN=function(x){x$parameters})
  names(pl)<- row.names(df)
 
  # print information
    if(!silent)
  	{
  	  	cat ('\nComparing ', nm, ' models [n = ', nn, ',', ' method = ', methv[1], ']\n\n', sep='')
  	  	print (df)
  	}
 
 if(silent)		return(list(modelFits=df, parameters=pl))
 else 			invisible(df)
}


fit3models<- function (y, silent = FALSE, method = c("Joint", "AD"), ...) 
{
    args<- list(...)
    check.var<- TRUE
	if(length(args)>0) if(args$pool==FALSE)	check.var<- FALSE
	if (check.var){
			tv<- test.var.het(y)
			pv<- round(tv$p.value, 0)
			wm<- paste("Sample variances not equal (P = ", pv,  "); consider using argument pool=FALSE", collapse="")
			if(pv <= 0.05)	warning(wm)
	}
    
    method <- match.arg(method)
    if (method == "AD") {
        m1 <- opt.GRW(y, ...)
        m2 <- opt.URW(y, ...)
        m3 <- opt.Stasis(y, ...)
    }
    else if (method == "Joint") {
        m1 <- opt.joint.GRW(y, ...)
        m2 <- opt.joint.URW(y, ...)
        m3 <- opt.joint.Stasis(y, ...)
    }
    mc <- compareModels(m1, m2, m3, silent = silent)
    invisible(mc)
}

fit4models<- function(y, silent=FALSE, method=c("Joint", "AD"), ...)
{
    args<- list(...)
    check.var<- TRUE
	if(length(args)>0) if(args$pool==FALSE)	check.var<- FALSE
	if (check.var){
			tv<- test.var.het(y)
			pv<- round(tv$p.value, 0)
			wm<- paste("Sample variances not equal (P = ", pv,  "); consider using argument pool=FALSE", collapse="")
			if(pv <= 0.05)	warning(wm)
	}

  method <- match.arg(method)
  if (method == "AD") {
        m1 <- opt.GRW(y, ...)
        m2 <- opt.URW(y, ...)
        m3 <- opt.Stasis(y, ...)
        m4 <- opt.StrictStasis(y, ...)
    }
  else if (method == "Joint") {
        m1 <- opt.joint.GRW(y, ...)
        m2 <- opt.joint.URW(y, ...)
        m3 <- opt.joint.Stasis(y, ...)
        m4 <- opt.joint.StrictStasis(y, ...)
    }
 
  mc <- compareModels(m1, m2, m3, m4, silent = silent)
  invisible(mc)
}

fit9models<- function(y, silent=FALSE, method=c("Joint", "AD"), ...)
{
    args<- list(...)
    check.var<- TRUE
	if(length(args)>0) if(args$pool==FALSE)	check.var<- FALSE
	if (check.var){
			tv<- test.var.het(y)
			pv<- round(tv$p.value, 0)
			wm<- paste("Sample variances not equal (P = ", pv,  "); consider using argument pool=FALSE", collapse="")
			if(pv <= 0.05)	warning(wm)
	}

  method <- match.arg(method)
  if (method == "AD") {
  		if(!silent) cat("Fitting simple models...\n")
        m4 <- opt.GRW(y, ...)
        m3 <- opt.URW(y, ...)
        m2 <- opt.Stasis(y, ...)
        m1 <- opt.StrictStasis(y, ...)
    }
  else if (method == "Joint") {      
        if(!silent)	cat("Fitting simple models...\n")
        m4 <- opt.joint.GRW(y, ...)
        m3 <- opt.joint.URW(y, ...)
        m2 <- opt.joint.Stasis(y, ...)
        m1 <- opt.joint.StrictStasis(y, ...)
    }
    
  if(!silent)	cat("Fitting punctuational model...\n")
  m5 <- fitGpunc(y, ng=2, method=method, silent=silent, ...)
  if(!silent)	cat("Fitting Stasis-URW model...\n")
  m6 <- fitModeShift(y, order="Stasis-RW", rw.model="URW", method=method, silent=silent, ...)
  if(!silent)	cat("Fitting Stasis-GRW model...\n")
  m7 <- fitModeShift(y, order="Stasis-RW", rw.model="GRW", method=method, silent=silent, ...)
  if(!silent)  cat("Fitting URW-Stasis model...\n")  
  m8 <- fitModeShift(y, order="RW-Stasis", rw.model="URW", method=method, silent=silent, ...)  
  if(!silent)  cat("Fitting GRW-Stasis model...\n")  
  m9 <- fitModeShift(y, order="RW-Stasis", rw.model="GRW", method=method, silent=silent,...) 

 
  mc <- compareModels(m1, m2, m3, m4, m5, m6, m7, m8, m9, silent = silent)
  invisible(mc)
}



logL.joint.StrictStasis<- function (p, y)
{
  theta<- p[1]
  n<- length(y$mm)
  VV<- diag(y$vv/y$nn)
  detV<- det(VV)
  invV<- solve(VV)
  M<- rep(theta, n)	
  #S<- dmvnorm(x$mm, mean = M, sigma = VV, log = TRUE)
  S<- dmnorm(y$mm, mean = M, varcov = VV, log = TRUE)
  
  return(S)
}

opt.joint.StrictStasis<- function (y, pool = TRUE, cl = list(fnscale = -1), hess = FALSE) 
{
    if (pool) 
        y <- pool.var(y, ret.paleoTS = TRUE)
    p0 <- mean(y$mm)
    names(p0)<- "theta"
    w <- optim(p0, fn = logL.joint.StrictStasis, control = cl, method = "Brent", lower=min(y$mm), upper=max(y$mm), 
    	hessian = hess, y = y)
    names(w$par)<- "theta"
    if (hess) 
        w$se <- sqrt(diag(-1 * solve(w$hessian)))
    wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName="StrictStasis", method="Joint", K=1, n=length(y$mm), se=w$se)
    return(wc)
}

logL.StrictStasis<- function(p,y)
{
  M <- p[1]
  dy <- diff(y$mm)
  nd <- length(dy)
  sv <- y$vv/y$nn
  svD <- sv[2:(nd + 1)]
  anc <- y$mm[1:nd]
  S <- dnorm(x = dy, mean = M - anc, sd = sqrt(svD), log = TRUE)
  return(sum(S))
}

opt.StrictStasis<- function (y, pool = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE) 
{
    p0 <- mean(y$mm)
    names(p0) <- c("theta")
    if (is.null(cl$ndeps)) 
        cl$ndeps <- min(abs(p0/10000), 1e-7)
    if (pool) 
        y <- pool.var(y, ret.paleoTS = TRUE)
    if (meth == "L-BFGS-B") 
        w <- try(optim(p0, fn = logL.StrictStasis, method = meth, lower = c(NA, 
            0), control = cl, hessian = hess, y = y), silent = TRUE)
    else w <- try(optim(p0, fn = logL.StrictStasis, method = meth, 
        control = cl, hessian = hess, y = y), silent = TRUE)
    if (class(w) == "try-error") {
        cl$ndeps <- rep(1e-09, length(p0))
        if (meth == "L-BFGS-B") 
            w <- try(optim(p0, fn = logL.StrictStasis, method = meth, 
                lower = c(NA, 0), control = cl, hessian = hess, 
                y = y), silent = TRUE)
        else w <- try(optim(p0, fn = logL.StrictStasis, method = meth, 
            control = cl, hessian = hess, y = y), silent = TRUE)
        if (class(w) == "try-error") {
            warning("opt.StrictStasis failed ", immediate. = TRUE)
            w$par <- c(NA, NA)
            w$value <- NA
        }
    }
    if (hess) 
        w$se <- sqrt(diag(-1 * solve(w$hessian)))
    else w$se <- NULL
    wc <- as.paleoTSfit(logL = w$value, parameters = w$par, modelName = "StrictStasis", 
        method = "AD", K = 1, n = length(y$mm) - 1, se = w$se)
    return(wc)
}




logL.GRW<- function(p,y)
# function to return log-likelihood of step mean and variance (M,V)= p,
# given a paleoTS object
{
  # get parameter values: M is Mstep, V is Vstep
  M<- p[1]
  V<- p[2]
  
  dy<- diff(y$mm)
  dt<- diff(y$tt)
  nd<- length(dy)  # number of differences
  sv<- y$vv/y$nn
  svA<- sv[1:nd]
  svD<- sv[2:(nd+1)]
  svAD<- svA + svD

  S<- dnorm(x=dy, mean=M*dt, sd=sqrt(V*dt + svAD), log=TRUE)
  return(sum(S))
}


logL.URW<- function(p,y)
# function to return log-likelihood of step mean and variance (M,V)= p,
# given a paleoTS object
{
  # get parameter values: V is Vstep
  V<- p[1]
    
  dy<- diff(y$mm)
  dt<- diff(y$tt)
  nd<- length(dy)  # number of differences
  sv<- y$vv/y$nn
  svA<- sv[1:nd]
  svD<- sv[2:(nd+1)]
  svAD<- svA + svD

  S<- dnorm(x=dy, mean=rep(0,nd), sd=sqrt(V*dt + svAD), log=TRUE)
  return(sum(S))
}




logL.Stasis <- function(p, y)
## logL of stasis model
{
 # get parameter estimates
 M<-p[1]	# M is theta
 V<-p[2]	# V is omega
 dy<- diff(y$mm)
 nd<- length(dy)
 
 sv<- y$vv/y$nn
 svD<- sv[2:(nd+1)]  # only need sampling variance of descendant
 anc<- y$mm[1:nd]

 S<- dnorm(x=dy, mean=M-anc, sd=sqrt(V + svD), log=TRUE)
 return(sum(S))
}


logL.Mult<- function (p, yl, model=c("GRW", "URW", "Stasis"))
# calculate logL over multiple sequences
#  here, yl is a list of paleoTS sequences
{
  model<- match.arg(model)
  Smult<-0
  nseq<- length(yl)
  for (i in 1:nseq)
    {
     if (model=="URW")
          Smult<- Smult + logL.URW(p,yl[[i]])
     else if (model=="GRW")
          Smult<- Smult + logL.GRW (p,yl[[i]])
     else if (model=="Stasis")
     	  Smult<- Smult + logL.Stasis(c(p[i], p[nseq]), yl[[i]])
    }
  return (Smult)
}


logL.SameMs <- function (p, yl)
# computes logL over >1 sequence, of model in which all sequences have the 
# same directionality (Mstep), with different step variances
# y is list of nseq paleoTS objects, p is array of K+1 parameters {m, v-1,..v-nseq}
{
  Sm<-0
  nseq<- length(yl)
  for (i in 1:nseq)
   	 Sm<- Sm + logL.GRW(p=c(p[1],p[i+1]), yl[[i]])
   	
  return(Sm) 	
}


logL.SameVs <- function (p, yl)
# computes logL over >1 sequence, of model in which all sequences have the 
# same step variance, with different steo means
# y is list of nseq paleoTS objects, p is array of K+1 parameters {m-1,..m-nseq, vs}
{
  Sm<-0
  nseq<- length(yl)
  for (i in 1:nseq)
   	 Sm<- Sm + logL.GRW(p=c(p[i],p[nseq+1]), yl[[i]])
   	
  return(Sm) 	
}


mle.GRW<- function(y)
# Gives analytical parameter estimates (GRW), assuming:
#   evenly spaced samples (constant dt)
#	same sampling variance for each dx (=2*Vp/n)
# Will try with reasonable values even if assumptions are violated
{
 nn<- length(y$mm)-1 
 tt<- (y$tt[nn+1]-y$tt[1])/nn
 eps<- 2*pool.var(y)/round(median(y$nn))  # sampling variance
 dy<- diff(y$mm)
 my<- mean(dy)
 
 mhat<- my/tt
 vhat<- (1/tt)*( (1/nn)*sum(dy^2) - my^2 - eps)
 
 w<- c(mhat, vhat)
 names(w)<- c("mstep", "vstep")
 return(w)
}


mle.URW<- function(y)
# Gives analytical parameter estimates (URW), assuming:
#   evenly spaced samples (constant dt)
#	same sampling variance for each dx (=2*Vp/n)
# Will try with reasonable values even if assumptions are violated
{
 nn<- length(y$mm)-1 
 tt<- (y$tt[nn+1]-y$tt[1])/nn
 eps<- 2*pool.var(y)/round(median(y$nn))  # sampling variance
 dy<- diff(y$mm)
 my<- mean(dy)
 
 vhat<- (1/tt)*( (1/nn)*sum(dy^2) - eps)
 
 w<- vhat
 names(w)<- "vstep"
 return(w)
}


mle.Stasis <- function (y)
# analytical solution to stasis model
{
  ns<- length(y$mm)
  vp<- pool.var(y)
  th<- mean(y$mm[2:ns])
  om<- var(y$mm[2:ns]) - vp/median(y$nn)
  
  w<- c(th, om)
  names(w)<- c("theta", "omega")
  return(w)
}



opt.GRW<- function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize estimates of step mean and variance for GRW model
# y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.GRW(y)
  if (p0[2] <= 0)	p0[2]<- 1e-7
  names(p0)<- c("mstep", "vstep")
  if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
  
  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
   w<- try(optim(p0, fn=logL.GRW, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
  else 
   w<- try(optim(p0, fn=logL.GRW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if (class(w)=="try-error")
    {
	cl$ndeps<- rep(1e-9, length(p0)) 	#p0/10000
    if (meth=="L-BFGS-B")
     w<- try(optim(p0, fn=logL.GRW, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
    else 
     w<- try(optim(p0, fn=logL.GRW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

    if (class(w)=="try-error")	# if still fails
	  {
		  warning("opt.GRW failed ", immediate.=TRUE)
		  w$par<- c(NA,NA)
		  w$value<- NA
	  }
    }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='GRW', method='AD', K=2, n=length(y$mm)-1, se=w$se)
  
  return (wc)
}


opt.URW<- function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize estimates of step mean and variance for GRW model
# y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.URW(y)
  if (p0 <= 0)	p0<- 1e-7
  names(p0)<- "vstep"
  if (is.null(cl$ndeps))		cl$ndeps<- p0/1e4
  
  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
   w<- try(optim(p0, fn=logL.URW, method=meth, lower=0, control=cl, hessian=hess, y=y), silent=TRUE)
  else 
   w<- try(optim(p0, fn=logL.URW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if (class(w)=="try-error")
    {
	cl$ndeps<- rep(1e-9, length(p0)) 	#p0/10000
	if (meth=="L-BFGS-B")
     w<- try(optim(p0, fn=logL.URW, method=meth, lower=0, control=cl, hessian=hess, y=y), silent=TRUE)
    else 
     w<- try(optim(p0, fn=logL.URW, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

    if (class(w)=="try-error")	# if still fails
	  {
		  warning("opt.URW failed ", immediate.=TRUE)
		  w$par<- NA
		  w$value<- NA
	  }
    }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='URW', method='AD', K=1, n=length(y$mm)-1, se=w$se)
  
  return (wc)
}



opt.Stasis<- function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize estimates of step mean and variance for GRW model
# y is an paleoTS object
{
  # get initial parameter estimates
  p0<- mle.Stasis(y)
  if (p0[2] <= 0 || is.na(p0[2]))	p0[2]<- 1e-7
  names(p0)<- c("theta", "omega")
  if (is.null(cl$ndeps))		cl$ndeps<- abs(p0/1e4)
  
  # pool variance, if needed
  if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)

  if (meth=="L-BFGS-B")
   w<- try(optim(p0, fn=logL.Stasis, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
  else 
   w<- try(optim(p0, fn=logL.Stasis, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

  # if optim failed, set ndeps based on p0
  if (class(w)=="try-error")
    {
	cl$ndeps<- rep(1e-9, length(p0)) 	#p0/10000
    if (meth=="L-BFGS-B")
     w<- try(optim(p0, fn=logL.Stasis, method=meth, lower=c(NA,0), control=cl, hessian=hess, y=y), silent=TRUE)
    else 
     w<- try(optim(p0, fn=logL.Stasis, method=meth, control=cl, hessian=hess, y=y), silent=TRUE)

   if (class(w)=="try-error")	# if still fails
	  {
		  warning("opt.Stasis failed ", immediate.=TRUE)
		  w$par<- c(NA,NA)
		  w$value<- NA
	  }
    }

  # add more information to results (p0, SE, K, n, IC scores)
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='Stasis', method='AD', K=2, n=length(y$mm)-1, se=w$se)
  
  return (wc)
}



opt.Mult<- function (yl, cl=list(fnscale=-1), model=c("GRW", "URW", "Stasis"), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates single model across multiple sequences
# pool=TRUE will pool variances _within_ sequences
{
  if (class(yl)=="paleoTS")
     stop("opt.RW.mult is only for multiple paleoTS sequences\n")
  nseq<- length(yl)
     
  # pool variances if needed
  if (pool){
  	for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }

  model<- match.arg(model)
  if (model=="GRW"){
     	ll<- c(NA,0)
     	pp<- sapply(yl, mle.GRW)
     	p0<- apply(pp, 1, median)	}
  else if (model=="URW")
     { ll<- 0
     	pp<- sapply(yl, mle.URW)
     	p0<- median(pp)	
     	names(p0)<- "vstep"}
  else if(model=="Stasis")
  	 {  ll<- c(rep(NA, nseq), 0)
  	 	pp<- sapply(yl, mle.Stasis)
  	 	p0<- c(pp[1,], median(pp[2,]))
  	 	names(p0)<- c(paste("theta", 1:nseq, sep=""), "omega")  		
  	 }
  K<- length(p0)
  if (is.null(cl$ndeps))	cl$ndeps<- rep(1e-8, length(p0))

  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.Mult, method=meth, lower=ll, control=cl, hessian=hess, yl=yl, model=model), silent=TRUE)
  else if (meth=="BFGS")
    w<- try(optim(p0, fn=logL.Mult, method=meth, control=cl, hessian=hess, yl=yl, model=model), silent=TRUE)
  
  # add more information to results (p0, SE, K, n, IC scores)
  nn<- sapply(yl, FUN=function(x) length(x$mm)) 
  n<- sum(nn) - nseq
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste(model,'.Mult', sep=''), method="AD", K=K, n=n, se=w$se)
  
  return (wc)
}


opt.joint.Mult<- function (yl, cl=list(fnscale=-1), model=c("GRW", "URW", "Stasis"), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates single model across multiple sequences
# pool=TRUE will pool variances _within_ sequences
{
  if (class(yl)=="paleoTS" || class(yl[[1]]) != "paleoTS")
     stop("opt.joint.Mult() is only for a list of multiple paleoTS sequences\n")
  nseq<- length(yl)
     
  # pool variances if needed
  if (pool){
  	for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }

  model<- match.arg(model)
  if (model=="GRW"){
     	ll<- c(rep(NA, nseq), NA,0)
     	anc0<- sapply(yl, FUN=function(x) x$mm[1])
     	P0<- sapply(yl, FUN=mle.GRW)
     	p0<- c(anc0, apply(P0, 1, median))
     	if(p0[nseq+2]<=0)	p0[nseq+2]<- 1e-7
     	names(p0)<- c(paste("anc", 1:nseq, sep=""), "mstep", "vstep")	
     	K<- nseq+2	 	}
  else if (model=="URW"){
     	ll<- c(rep(NA, nseq),0)
     	anc0<- sapply(yl, FUN=function(x) x$mm[1])
     	P0<- sapply(yl, FUN=mle.URW)
     	p0<- c(anc0, median(P0))	
     	if(p0[nseq+1]<=0)	p0[nseq+1]<- 1e-7
     	     	names(p0)<- c(paste("anc", 1:nseq, sep=""), "vstep")
     	K<- nseq+1	 	}
  else if (model=="Stasis"){
     	ll<- c(rep(NA, nseq),0)
     	P0<- sapply(yl, FUN=mle.Stasis)
     	p0<- c(P0[1,], median(P0[2,]))	
     	if(p0[nseq+1]<=0)	p0[nseq+1]<- 1e-7
     	names(p0)<- c(paste("theta", 1:nseq, sep=""), "omega")
     	K<- nseq+1	 	}

  if(is.null(cl$ndeps)) cl$ndeps<- rep(1e-8, length(p0))
  
  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.joint.Mult, method=meth, lower=ll, control=cl, hessian=hess, yl=yl, model=model), silent=TRUE) 
  else     w<- try(optim(p0, fn=logL.joint.Mult, method=meth, control=cl, hessian=hess, yl=yl, model=model), silent=TRUE)

  # add more information to results (p0, SE, K, n, IC scores)
  n<- sum(sapply(yl, FUN=function(x)length(x$mm)))

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))	
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste(model,'.Mult', sep=''), method='Joint', K=K, n=n, se=w$se)
  
  return (wc)
}

logL.joint.Mult<- function (p, yl, model=c("GRW", "URW", "Stasis"))
# calculate logL over multiple sequences
#  here, yl is a list of paleoTS sequences
{
  model<- match.arg(model)
  Smult<-0
  nseq<- length(yl)
  for (i in 1:nseq)
    {
     if (model=="URW")
          Smult<- Smult + logL.joint.URW(p=c(p[i], p[nseq+1]), yl[[i]])	# first nseq are anc_i, then others
     else if (model=="GRW")
          Smult<- Smult + logL.joint.GRW (p=c(p[i], p[nseq+1], p[nseq+2]), yl[[i]])
     else if (model=="Stasis")
          Smult<- Smult + logL.joint.Stasis (p=c(p[i], p[nseq+1]), yl[[i]])
    }
  return (Smult)
}

  
fitMult<- function(yl, model=c("GRW", "URW", "Stasis", "covTrack"), method=c("Joint", "AD"), pool=TRUE, zl=NULL, hess=FALSE)
{
  model<- match.arg(model)
  method<- match.arg(method)
	
  if(model=="covTrack"){
  	if(method=="Joint")	w<- opt.joint.covTrack.Mult(yl, zl, pool=pool, hess=hess)
  	if(method=="AD")	w<- opt.covTrack.Mult(yl, zl, pool=pool, hess=hess)	
  }	else{
  if(method=="AD")		w<- opt.Mult(yl, model=model, pool=pool)
  if(method=="Joint")	w<- opt.joint.Mult(yl, model=model, pool=pool)
  }
  	
  return(w)
}

fitSimple<- function(y, model=c("GRW", "URW", "Stasis", "StrictStasis", "OU", "covTrack"), method=c("Joint", "AD"), pool=TRUE, z=NULL, hess=FALSE)
{
  model<- match.arg(model)
  method<- match.arg(method)
  
  if(pool){
  		tv<- test.var.het(y)
		pv<- round(tv$p.value, 0)
		wm<- paste("Sample variances not equal (P = ", pv,  "); consider using argument pool=FALSE", collapse="")
		if(pv <= 0.05)	warning(wm)
  	}

  if(method=="Joint"){
  	if(model=="GRW")			w<- opt.joint.GRW(y, pool=pool, hess=hess)
  	if(model=="URW")			w<- opt.joint.URW(y, pool=pool, hess=hess)
  	if(model=="Stasis")			w<- opt.joint.Stasis(y, pool=pool, hess=hess)
  	if(model=="StrictStasis")	w<- opt.joint.StrictStasis(y, pool=pool, hess=hess)  	
  	if(model=="OU")				w<- opt.joint.OU(y, pool=pool, hess=hess)
  	if(model=="covTrack")	{
  								if(is.null(z))	stop("Covariate [z] needed for covTrack model.")
  								w<- opt.joint.covTrack(y, z, pool=pool, hess=hess)}
  }
  if(method=="AD"){
  	if(model=="GRW")			w<- opt.GRW(y, pool=pool, hess=hess)
  	if(model=="URW")			w<- opt.URW(y, pool=pool, hess=hess)
  	if(model=="Stasis")			w<- opt.Stasis(y, pool=pool, hess=hess)
  	if(model=="StrictStasis")	w<- opt.StrictStasis(y, pool=pool, hess=hess)  	
  	if(model=="OU")				stop("AD method not available for OU model.  Consider using Joint method.")
  	if(model=="covTrack")	{
  								if(is.null(z))	stop("Covariate [z] needed for covTrack model.")
	  							w<- opt.joint.covTrack(y, z, pool=pool, hess=hess) }
  }
	
  return(w)			
}

opt.RW.SameMs<- function (yl, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates shared Ms model across multiple sequences
{
  if (class(yl)=="paleoTS")
  	stop("Function opt.SameMs() is only meaningful for multiple sequences.\n")
  nseq<- length(yl)	
   
  # pool variances if needed
  if (pool){
  	for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }
  
  # generate initial parameter estimates {ms, v1,..vk}
  p0m<- array(dim=nseq)
  p0v<- array(dim=nseq)
  for (i in 1:nseq)
   {
   	gg<- mle.GRW(yl[[i]])	
  	p0m[i]<- gg[1]
  	if (gg[2] <= 0)	gg[2]<- 1e-7	# handle negative initial estimates
  	p0v[i]<- gg[2]
   }
  p0<- c(median(p0m), p0v)  # shared Ms, followed by separate Vs for each sequence
  names(p0)<- c("mstep", paste("vstep", 1:nseq, sep=""))
if (is.null(cl$ndeps))	cl$ndeps<- rep(1e-8, length(p0))  

  # optimize logL
  ll<- c(NA, rep(0,nseq))
  if (meth=="L-BFGS-B")
    w<- try(optim(p0, fn=logL.SameMs, method=meth, lower=ll, control=cl, hessian=hess, yl=yl), silent=TRUE)
  else
    w<- try(optim(p0, fn=logL.SameMs, method=meth, control=cl, hessian=hess, yl=yl), silent=TRUE)

  # add more information to results (p0, SE, K, n, IC scores)
  n<-0
  for (i in 1:nseq)	n<- n + (length(yl[[i]]$mm)-1)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName="sameMs.Mult", method='AD', K=nseq+1, n=n, se=w$se)
  
  return (wc)
}


opt.RW.SameVs<- function (yl, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
# estimates shared Vs model across multiple sequences
{
  if (class(yl)=="paleoTS")
  	stop("Function opt.SameVs() is only meaningful for multiple sequences.\n")
  nseq<- length(yl)	
   
  # pool variances if needed
  if (pool){
  	for (i in 1:nseq)	yl[[i]]<- pool.var(yl[[i]], ret.paleoTS=TRUE)  }
  
  # generate initial parameter estimates {ms, v1,..vk}
  p0m<- array(dim=nseq)
  p0v<- array(dim=nseq)
  for (i in 1:nseq)
   {
   	gg<- mle.GRW(yl[[i]])	
  	p0m[i]<- gg[1]
  	if (gg[2] <= 0)	gg[2]<- 1e-7	# handle negative initial estimates
  	p0v[i]<- gg[2]
   }
  p0<- c(p0m, median(p0v))  # separate Ms, followed by shared Vs for each sequence
  names(p0)<- c(paste("mstep", 1:nseq, sep=""), "vstep")
  if (is.null(cl$ndeps))	cl$ndeps<- rep(1e-8, length(p0))  
  
  # optimize logL
  ll<- c(rep(NA,nseq), 0)
  if (meth=="L-BFGS-B")
    w<- optim(p0, fn=logL.SameVs, method=meth, lower=ll, control=cl, hessian=hess, yl=yl)
  else
    w<- optim(p0, fn=logL.SameVs, method=meth, control=cl, hessian=hess, yl=yl)

  # add more information to results (p0, SE, K, n, IC scores)
  n<-0
  for (i in 1:nseq)	n<- n + (length(yl[[i]]$mm)-1)

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName="sameVs.Mult", method='AD', K=nseq+1, n=n, se=w$se)

  return (wc)
}



sim.GRW <- function (ns=20, ms=0, vs=0.1, vp=1, nn=rep(20,ns), tt=0:(ns-1))
# simulates GRW; ns= number of samples, ms=mean and vs=variance of the step distribution,
#  vp=population variance, tt=ages of the samples
{
 MM<- array(dim=ns)
 mm<- array(dim=ns)
 vv<- array(dim=ns)
 dt<- diff(tt)
 
 inc<- rnorm(ns-1, ms*dt, sqrt(vs*dt))	# evolutionary increments
  
 MM<- cumsum(c(0,inc))	# true means
 mm<- MM + rnorm(ns, 0, sqrt(vp/nn))	# true means plus sampling error
 vv<- rep(vp, ns)
 
 gp<- c(ms, vs)
 names(gp)<- c("mstep", "vstep")
 
 res<- as.paleoTS(mm=mm, vv=vv, nn=nn, tt=tt, MM=MM, genpars=gp, label="Created by sim.GRW", reset.time=FALSE) 
 return(res)
}


sim.Stasis <- function(ns=20, theta=0, omega=0, vp=1, nn=rep(20,ns), tt=0:(ns-1))
# simulate stasis
{
 xmu<- rnorm(ns, mean=theta, sd=sqrt(omega))
 xobs<- xmu + rnorm(ns, 0, sqrt(vp/nn))
 gp<- c(theta, omega)
 names(gp)<- c("theta", "omega")
   	
 x <- as.paleoTS(mm=xobs,vv=rep(vp,ns),nn=nn,tt=tt,MM=xmu,genpars=gp,label="Created by sim.Stasis", reset.time=FALSE) 
 return(x)	 	 	
}

lynchD<- function (y, gen.per.t=1e6, pool=TRUE, method=c('Joint', 'AD'), ...)
# compute Lynch's rate metric, Delta
# gen.per.t is the number of generations per unit tt
# ... are further arguments to opt.URW()
{
  method<- match.arg(method)
  vp<- pool.var(y)
  if (method == 'AD')
    {
      wu<- opt.URW(y, pool=pool)
      vs<- unname(wu$par)	
    }
  if (method == 'Joint')
    {
      wu<- opt.joint.URW(y, pool=pool)
      vs<- unname(wu$par[2])	
    }
    
  D<- 0.5*(vs/vp)/gen.per.t
  
  drift.min<- 5e-5
  drift.max<- 5e-3
  
  if (D < drift.min)	res<- "Slower than drift range"
  else if (D > drift.max)	res<- "Faster than drift range"
  else					res<- "Within range of drift"
  
  w<- list(D=D, pooled.var=vp, gen.per.t=gen.per.t, vstep=vs, drift.range=c(drift.min, drift.max), result=res) 
  return(w)	
}

## functions added from paleoTSalt


logL.joint.GRW<- function (p, y)
# returns logL of GRW model for paleoTS object x
# p is vector of parameters: alpha, ms, vs
{
 # prepare calculations
 anc<- p[1]
 ms<- p[2]
 vs<- p[3]
 n<- length(y$mm)
 
 # compute covariance matrix
 VV<- vs*outer(y$tt, y$tt, FUN=pmin)
 diag(VV)<- diag(VV) + y$vv/y$nn 
 
 # compute logL based on multivariate normal
 M<- rep(anc, n) + ms*y$tt
 S<- dmnorm(y$mm, mean=M, varcov=VV, log=TRUE)
 
 return(S)		
}

logL.joint.URW<- function(p, y)
# returns logL of URW model for paleoTS object x
# p is vector of parameters: alpha, vs
{
 # prepare calculations
 anc<- p[1]
 vs<- p[2]
 n<- length(y$mm)
 
 # compute covariance matrix
 VV<- vs*outer(y$tt, y$tt, FUN=pmin)
 diag(VV)<- diag(VV) + y$vv/y$nn 	# add sampling variance
 
 # compute logL based on multivariate normal
 M<- rep(anc, n)
 #S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
 #S<- dmvnorm(x$mm, mean=M, sigma=VV, log=TRUE)
 S<- dmnorm(y$mm, mean=M, varcov=VV, log=TRUE)

 return(S)
}

logL.joint.Stasis<- function (p, y)
# returns logL of Stasis model for paleoTS object x
# p is vector of parameters: theta, omega
{
 # prepare calculations
 theta<- p[1]
 omega<- p[2]
 n<- length(y$mm)
 
 # compute covariance matrix
 VV<- diag(omega + y$vv/y$nn)  # omega + sampling variance
 
 # compute logL based on multivariate normal
 M<- rep(theta, n)
 #S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
 #S<- dmvnorm(x$mm, mean=M, sigma=VV, log=TRUE)
 S<- dmnorm(y$mm, mean=M, varcov=VV, log=TRUE)
 
 return(S)	 	
}


opt.joint.GRW<- function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize GRW model using alternate formulation
{
 ## check if pooled, make start at tt=0
 if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
 if(y$tt[1] != 0)	stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  
 ## get initial estimates
 p0<- array(dim=3)
 p0[1]<- y$mm[1]	
 p0[2:3]<- mle.GRW(y)
 if (p0[3]<=0)	p0[3]<- 1e-7
 names(p0)<- c("anc", "mstep", "vstep")
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-8
 
 if (meth=="L-BFGS-B")	w<- optim(p0, fn=logL.joint.GRW, control=cl, method=meth, lower=c(NA,NA,0), hessian=hess, y=y)
 else 					w<- optim(p0, fn=logL.joint.GRW, control=cl, method=meth, hessian=hess, y=y)

 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
 else			w$se<- NULL
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='GRW', method='Joint', K=3, n=length(y$mm), se=w$se)
  
 return (wc)
}

opt.joint.URW<- function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize URW model using alternate formulation
{
 ## check if pooled
 if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
 if(y$tt[1] != 0)	stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")
 
 ## get initial estimates
 p0<- array(dim=2)
 p0[1]<- y$mm[1]
 p0[2]<- min(c(mle.URW(y), 1e-7)) ## handles negative vstep estimates
 names(p0)<- c("anc","vstep")
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-8
 
 if (meth=="L-BFGS-B")	w<- optim(p0, fn=logL.joint.URW, control=cl, method=meth, lower=c(NA,0), hessian=hess, y=y)
 else					w<- optim(p0, fn=logL.joint.URW, control=cl, method=meth, hessian=hess, y=y)

 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
 else			w$se<- NULL
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='URW', method='Joint', K=2, n=length(y$mm), se=w$se)
  
 return (wc)
}

opt.joint.Stasis<- function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize Stasis model using alternate formulation
{
 ## check if pooled
 if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
 if(y$tt[1] != 0)	stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")
  
 ## get initial estimates
 p0<- mle.Stasis(y)
 if(p0[2]<=0 || is.na(p0[2]))	p0[2]<- 1e-7
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-9
  
 if(meth=="L-BFGS-B")  w<- optim(p0, fn=logL.joint.Stasis, control=cl, method=meth, lower=c(NA,0), hessian=hess, y=y)
 else				  w<- optim(p0, fn=logL.joint.Stasis, control=cl, method=meth, hessian=hess, y=y)
 
 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
 else			w$se<- NULL
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='Stasis', method='Joint', K=2, n=length(y$mm), se=w$se)
  
 return (wc)
	
}


 # functions to compute mean and variance of OU process
 ou.M<- function(anc, theta, aa, tt) theta*(1 - exp(-aa*tt)) + anc*exp(-aa*tt)
 ou.V<- function(vs, aa, tt)        (vs/(2*aa))*(1 - exp(-2*aa*tt))


sim.OU<- function (ns=20, anc=0, theta=10, alpha=0.3, vs=0.1, vp=1, nn=rep(20, ns), tt=0:(ns-1))
## generate a paleoTS sequence according to an OU model
{
    
    MM <- array(dim = ns)
    mm <- array(dim = ns)
    vv <- array(dim = ns)
    dt <- diff(tt)
    MM[1] <- anc
    x <- rnorm(nn[1], mean = MM[1], sd = sqrt(vp))
    mm[1] <- mean(x)
    vv[1] <- var(x)
    for (i in 2:ns) {
        ex<- ou.M(MM[i-1], theta, alpha, dt[i-1])
        vx<- ou.V(vs, alpha, dt[i-1])
        MM[i]<- rnorm(1, ex, sqrt(vx))    
        x <- rnorm(nn[i], mean = MM[i], sd = sqrt(vp))
        mm[i] <- mean(x)
        vv[i] <- var(x)
    }
    
    gp <- c(anc, theta, alpha, vs)
    names(gp) <- c("anc", "theta", "alpha", "vs")
    res <- as.paleoTS(mm = mm, vv = vv, nn = nn, tt = tt, MM = MM, 
        genpars = gp, label = "Created by sim.OU()", reset.time=FALSE)

    return(res)
}

logL.joint.OU<- function(p, y)
# returns logL of OU model for paleoTS object x
# 
{
 # prepare calculations
 anc<- p[1]
 vs<- p[2]
 theta<- p[3]
 aa<- p[4]
 n<- length(y$mm)
 
 # compute covariance matrix
 ff<- function (a,b) abs(a-b)
 VV<- outer(y$tt, y$tt, FUN=ff)
 VV<- exp(-aa*VV)
 VVd<- ou.V(vs,aa,y$tt)
 VV2<- outer(VVd,VVd,pmin)
 VV<- VV*VV2
 diag(VV)<- VVd + y$vv/y$nn 	# add sampling variance
 #detV<- det(VV)
 #invV<- solve(VV)
 
 # compute logL based on multivariate normal
 M<- ou.M(anc, theta, aa, y$tt)
 #S<- -0.5*log(detV) - 0.5*n*log(2*pi) - ( t(x$mm-M)%*%invV%*%(x$mm-M) ) / (2)
 #S<- dmvnorm(t(x$mm), mean=M, sigma=VV, log=TRUE)
 S<- dmnorm(t(y$mm), mean=M, varcov=VV, log=TRUE)
 return(S)		
}


opt.joint.OU<- function (y, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE)
# optimize OU model using tree methods
{
 ## check if pooled
 if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)
 if(y$tt[1] != 0)	stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")
 
 ## get initial estimates
 w0<- mle.GRW(y)
 halft<- (y$tt[length(y$tt)]-y$tt[1])/4			# set half life to 1/4 of length of sequence
 p0<- c(y$mm[1], w0[2]/10, y$mm[length(y$mm)], log(2)/halft)
 names(p0)<- c("anc","vstep","theta","alpha")
 #print(p0)
 if (is.null(cl$ndeps))	cl$ndeps<- abs(p0/1e4)
 cl$ndeps[cl$ndeps==0]<- 1e-8
 
  
 if(meth=="L-BFGS-B")  w<- optim(p0, fn=logL.joint.OU, control=cl, method=meth, lower=c(NA,1e-10,NA,1e-8), hessian=hess, y=y)
 else 				  w<- optim(p0, fn=logL.joint.OU, control=cl, method=meth, hessian=hess, y=y) 

 if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
 else			w$se<- NULL
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='OU', method='Joint', K=4, n=length(y$mm), se=w$se)
  
 return (wc)
}


##### start of punctuate functions

cat.paleoTS<- function (y)
# concatenates multiple paleoTS objects, with y a list of paleoTS objects
{
 x<- y[[1]]
 for (i in 2:length(y))
   {
   	x$mm<- append(x$mm, y[[i]]$mm)
   	x$vv<- append(x$vv, y[[i]]$vv)
   	x$tt<- append(x$tt, y[[i]]$tt)
   	x$MM<- append(x$MM, y[[i]]$MM)
   	x$nn<- append(x$nn, y[[i]]$nn)
   }
  
  return (x)   	
}


sim.punc<- function (ns=c(10,10), theta=c(0,1), omega=rep(0,length(theta)), nn=rep(30,sum(ns)), tt=0:(sum(ns)-1), vp=1)
# simulate punctuated sequence; theta and omega are vectors of paramters
# ns is vector of ns in each sub-sequence
{
  nr<- length(theta)
  xl<- list()
  for (i in 1:nr)
   { 
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else 		
   	 {
   	  start.i<- sum(ns[1:(i-1)])+1
   	  end.i<- start.i + ns[i] -1		
   	 }
   	 
   	 xl[[i]]<- sim.Stasis(ns=ns[i], theta=theta[i], omega=omega[i], vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   }
  
  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.punc()"
  shft<- cumsum(ns[-nr])+1
  y$genpars<- c(theta, omega, shft)
  names(y$genpars)<- c(paste("theta",1:nr,sep=""), paste("omega",1:nr,sep=""), paste("shift",1:(nr-1),sep=""))
  return (y)  	
}


split4punc<- function (y, gg, overlap=TRUE)
# divides a paleoTS object (y) into a several paleoTS objects, according to vector 'gg'
# gg is a vectors of 1,2,3.. indicating groupings
# overlap=TRUE means that the adjacent samples are included 
{
  yl<- list()
  ng<- max(gg)
  for (i in 1:ng)
   {
   	 ok<- gg==i
   	 if(i>1 & overlap==TRUE)	ok[max(which(gg==i-1))]<- TRUE   # this is now right!
   	
   	 yl[[i]]<- as.paleoTS(y$mm[ok],y$vv[ok],y$nn[ok],y$tt[ok],y$MM[ok])
   }
  return (yl)	
}


logL.punc<- function (p, y, gg)
# logL of punctuation, with shifts 
{
  ng<- length(p)/2
  th<- p[1:ng]
  om<- p[(ng+1):(2*ng)]
  
  xl<- split4punc(y,gg)
  S<-0
  for (i in 1:ng)
    S<- S+ logL.Stasis(p=c(th[i], om[i]), xl[[i]])
  
  return (S)	
}

logL.punc.omega<- function(p,y,gg)
# logL of punctuation, assuming omega is same over all sections
{
  ng<- length(p)-1
  th<- p[1:ng]
  om<- p[ng+1]
  
  xl<- split4punc(y,gg)
  S<-0
  for (i in 1:ng)
    S<- S + logL.Stasis(p=c(th[i], om), xl[[i]])
  
  return (S)	  	
}

opt.punc<- function(y, gg, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE, oshare) 
{
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
 
  ng<- max(gg)
  mg<- tapply(y$mm, gg, mean)
  mv<- tapply(y$mm, gg, var)
  pn<- paste("theta", 1:ng, sep="")
  if(oshare)		{ p0<- c(mg, mean(mv)); K<- 2*ng; pn<- c(pn, "omega")}
  else				{ p0 <- c(mg, mv); K <- 3*ng-1; pn<- c(pn, paste("omega", 1:ng, sep=""))}
  names(p0)<- pn
  
  cl$ndeps <- p0/100
  if (oshare)
  {
   if (meth == "L-BFGS-B") w <- optim(p0, fn=logL.punc.omega, gg=gg, method = meth, lower = c(rep(NA,ng), 0), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn=logL.punc.omega, gg=gg, method = meth, control = cl, hessian=hess, y=y)
  }
  else
  {
   if (meth == "L-BFGS-B")  w <- optim(p0, fn=logL.punc, gg=gg, method = meth, lower = c(rep(NA,ng), rep(0,ng)), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn = logL.punc, gg=gg, method = meth, control = cl, hessian = hess, y = y)
  }
  
  # add more information to results
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste('Punc', ng-1, sep='-'), method='AD', K=K, n=length(y$mm)-1, se=w$se)

  return(wc)
}

logL.joint.punc<- function (p, y, gg)
# logL of punctuation, with shifts 
{
  ng<- length(p)/2
  th<- p[1:ng]
  om<- p[(ng+1):(2*ng)]
  
  M<- th[gg]  # vector of MVN means
  VV<- diag(om[gg] + y$vv/y$nn)  # vcv matrix
  #S<- dmvnorm(y$mm, mean=M, sigma=VV, log=TRUE)
  S<- dmnorm(y$mm, mean=M, varcov=VV, log=TRUE)
  
  return (S)	
}

logL.joint.punc.omega<- function(p,y,gg)
# logL of punctuation, assuming omega is same over all sections
{
  ng<- length(p)-1
  th<- p[1:ng]
  om<- p[ng+1]
  
  M<- th[gg]  # vector of MVN means
  omv<- rep(om, max(gg))
  VV<- diag(omv[gg] + y$vv/y$nn)  # vcv matrix
  #S<- dmvnorm(y$mm, mean=M, sigma=VV, log=TRUE)
  S<- dmnorm(y$mm, mean=M, varcov=VV, log=TRUE)
  
  return (S)	  	
}


opt.joint.punc<- function(y, gg, cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE, oshare) 
{
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
  if(y$tt[1] != 0)	stop("Initial time must be 0.  Use as.paleoTS() or read.paleoTS() to correctly process ages.")
 
  ng<- max(gg)
  mg<- tapply(y$mm, gg, mean)
  mv<- tapply(y$mm, gg, var)
  pn<- paste("theta", 1:ng, sep="")
  if(oshare)		{ p0<- c(mg, mean(mv)); K<- 2*ng; pn<- c(pn, "omega")}
  else				{ p0 <- c(mg, mv); K <- 3*ng-1; pn<- c(pn, paste("omega", 1:ng, sep=""))}
  names(p0)<- pn
  
  cl$ndeps <- p0/100
  if (oshare)
  {
   if (meth == "L-BFGS-B") w <- optim(p0, fn=logL.joint.punc.omega, gg=gg, method = meth, lower = c(rep(NA,ng), 0), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn=logL.joint.punc.omega, gg=gg, method = meth, control = cl, hessian=hess, y=y)
  }
  else
  {
   if (meth == "L-BFGS-B")  w <- optim(p0, fn=logL.joint.punc, gg=gg, method = meth, lower = c(rep(NA,ng), rep(0,ng)), control=cl, hessian=hess, y=y)
   else w <- optim(p0, fn = logL.joint.punc, gg=gg, method = meth, control = cl, hessian = hess, y = y)
  }
  
  # add more information to results
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste('Punc', ng-1, sep='-'), method='Joint', K=K, n=length(y$mm), se=w$se)

  return(wc)
}



shifts<- function(ns, ng, minb=7) 
{
	aaL<- combn(ns, ng-1, simplify=FALSE)
	aa<- matrix(unlist(aaL), nrow=ng-1)
	top<- rep(0, ncol(aa))
	bot<- rep(ns, ncol(aa))
	aaM<- unname(rbind(top, aa, bot))
	
	daa<- apply(aaM,2,diff)
	ok<- apply(daa, 2, function(x) all(x >= minb))		
	
	if(ng>2)	return(aa[,ok])	
	else		return(matrix(aa[,ok], nrow=1))
}



shift2gg<- function (ss, ns)
# ss is vector of shift points, ns is # samples
{
  z<- c(0,ss,ns+1)
  cc<- cut(1:ns, breaks=z, right=TRUE)
  gg<- as.numeric(cc)
  return(gg)	
}



fitGpunc<-function (y, ng = 2, minb = 7, pool = TRUE, oshare = TRUE, method = c("Joint", 
   "AD"), silent = FALSE, hess = FALSE, parallel = FALSE, ...) 
{
   method <- match.arg(method)
	if (pool){
			tv<- test.var.het(y)
			pv<- round(tv$p.value, 0)
			wm<- paste("Sample variances not equal (P = ", pv,  "); consider using argument pool=FALSE", collapse="")
			if(pv <= 0.05)	warning(wm)
	}


   if (ng == 1) {
      warning("Fitting stasis model (because ng=1)")
      if (method == "AD") 
         ww <- opt.Stasis(y, pool = pool, hess = hess, ...)
      else if (method == "Joint") 
         ww <- opt.joint.Stasis(y, pool = pool, hess = hess, ...)
      return(ww)
   }
   ns <- length(y$mm)
   GG <- shifts(ns, ng, minb = minb)
   nc <- ncol(GG)
   if (!silent) 
      cat("Total # hypotheses: ", nc, "\n")
   i<- 1  # done to avoid R Check error in response to foreach iterators
   
   
   if (parallel==TRUE) {
      # set up cluster
      cores<-parallel::detectCores()
      cl<-parallel::makeCluster(cores-1)
      doParallel::registerDoParallel(cl)
      # run loop
      wl<- foreach (i = 1:nc, .packages=c('paleoTS')) %dopar% {
         #if (!silent)  cat(i, " ")
         gg <- shift2gg(GG[, i], ns)
         if (method == "AD") 
            w <- opt.punc(y, gg, oshare = oshare, pool = pool, 
               hess = hess, ...)
         else if (method == "Joint") 
            w <- opt.joint.punc(y, gg, oshare = oshare, pool = pool, 
               hess = hess, ...)
         w  #return w to list wl
      }
      parallel::stopCluster(cl)	      # kill cluster
   
   # non-parallel version
   }else if (parallel==FALSE) {
      wl<- foreach (i = 1:nc, .packages=c('paleoTS')) %do% {
         if (!silent) cat(i, " ")
         gg <- shift2gg(GG[, i], ns)
         if (method == "AD") 
            w <- opt.punc(y, gg, oshare = oshare, pool = pool, 
               hess = hess, ...)
         else if (method == "Joint") 
            w <- opt.joint.punc(y, gg, oshare = oshare, pool = pool, 
               hess = hess, ...)
         w  #return w to list wl
      }
   }
   
   logl<-sapply(wl, function(x) x$logL) #extract logl values from wl
   if(!silent) cat("\n")
   winner <- which.max(logl)
   ww <- wl[[winner]]
   ss <- GG[, winner]
   names(ss) <- paste("shift", 1:(ng - 1), sep = "")
   ww$parameters <- append(ww$parameters, ss)
   ww$all.logl <- logl
   ww$GG <- GG
   return(ww)
}




sim.GRW.shift <- function (ns=c(10,10), ms=c(0,1), vs=c(0.5,0.5), nn=rep(30,sum(ns)), tt=0:(sum(ns)-1), vp=1)
# simulate sequence of GRW with different parameter values in different segments
# ns is vector of ns in each sub-sequence, similar for other parameters
{
  nr<- length(ms)
  cns<- cumsum(ns)
  xl<- list()
  for (i in 1:nr)
   { 
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else 		
   	 {
   	  start.i<- cns[i-1]+1
   	  end.i<- cns[i-1]+ns[i]	
   	 }
   	 
   	 xl[[i]]<- sim.GRW(ns=ns[i], ms=ms[i], vs=vs[i], vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   	 if (i>1)
   	 	xl[[i]]$mm<- xl[[i]]$mm + xl[[i-1]]$mm[length(xl[[i-1]]$mm)]
   }
  
  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.GRW.shift()"
  shft<- cumsum(ns[-nr])+1
  y$genpars<- c(ms, vs, shft)
  names(y$genpars)<- c(paste("mstep",1:nr,sep=""), paste("vstep",1:nr,sep=""), paste("shift",1:(nr-1),sep=""))  
  return (y)  	
}

opt.GRW.shift<- function(y, ng=2, minb=7, model=1, pool=TRUE, silent=FALSE)
## optimize for shifted GRW dynamics (with some min n per section)
## models:	1  grw (same Vs, diff Ms)
#			2  grw (same Ms, diff Vs)
#			3  urw (diff Vs)
#			4  grw (diff Ms, diff Vs)
{
 ns<- length(y$mm)
 GG<- shifts(ns, ng, minb=minb)

 nc<- ncol(GG)
 if (!silent)	cat ("\n\nTotal # hypotheses: ", nc, "\n")
 wl<- list()
 logl<- array(-Inf, dim=nc)
 
 for (i in 1:nc)
  {
    if (!silent) 	cat (i, " ")
    gg<- shift2gg(GG[,i], ns)
    yl<- split4punc(y, gg)
    #print(yl[[1]])
    if (model==1)	{ w<- opt.RW.SameVs(yl, pool=pool); w$modelName<- paste('GRWsameVs-shift', ng-1, sep='-'); w$K<- w$K + ng-1}
    if (model==2)	{ w<- opt.RW.SameMs(yl, pool=pool); w$modelName<- paste('GRWsameMs-shift', ng-1, sep='-'); w$K<- w$K + ng-1}
    if (model>=3)
     {
      wli<- list()
      totS<-0
      totpar<- numeric()
      for (j in 1:ng) 
       {
      	if (model==3){	
      		wli[[j]]<- opt.URW(yl[[j]], pool=pool)
      		names(wli[[j]]$parameters)<- paste("vstep", j, sep="")   }
      	if (model==4){	
      		wli[[j]]<- opt.GRW (yl[[j]], pool=pool)
      		names(wli[[j]]$parameters)<- paste(c("mstep","vstep"), j, sep="")   }
      	totS<- totS + wli[[j]]$logL
      	totpar<- c(totpar, wli[[j]]$parameters)
      	kk<- length(totpar)+ng-1
       }
    
      ifelse(model==3, mn<- 'URW-shift', mn<- 'GRW-shift')
      w<- as.paleoTSfit(logL=totS, parameters=totpar, modelName=paste(mn, ng-1, sep='-'), method='AD', K=kk, n=ns-1, se=NULL)	
     }
         
    logl[i]<- w$logL
    wl[[i]]<- w
  }

 # add more information to results
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ss<- GG[,winner]
 names(ss)<- paste("shift", 1:(ng-1), sep='')
 ww$parameters<- append(ww$parameters, ss)
 ww$all.logl<- logl
 ww$GG<- GG
 return(ww)
 
}

sim.sgs <- function (ns=c(20,20,20), theta=0, omega=1, ms=1, vs=0.1, nn=rep(30, sum(ns)), tt=0:(sum(ns)-1), vp=1)
# simulate stasis-grw-stasis sequence, take theta2 to be final value after grw part
{
  xl<- list()
  for (i in 1:3)
   { 
   	 if (i==1)
   	  {  start.i<- 1
   	  	 end.i<- ns[i]  }
   	 else 		
   	 {
   	  start.i<- sum(ns[1:(i-1)])+1
   	  end.i<- start.i + ns[i] -1	
   	 }
   	 
   	 if (i==2)
   	 	xl[[i]]<- sim.GRW(ns[2], ms, vs, nn=nn[start.i:end.i], tt=tt[start.i:end.i], vp=vp)
   	 
   	 else
   	 	xl[[i]]<- sim.Stasis(ns=ns[i], theta=theta, omega=omega, vp=vp, nn=nn[start.i:end.i], tt[start.i:end.i])
   }

  ## add offsets
  xl[[2]]$mm<- xl[[2]]$mm + xl[[1]]$MM[ns[1]]
  xl[[3]]$mm<- xl[[3]]$mm + xl[[2]]$MM[ns[2]]
	
  y<- cat.paleoTS(xl)
  y$label<- "Created by sim.sgs()"
  y$genpars <- c(theta, omega, ms, vs)
  names(y$genpars)<- c("theta","omega", "ms","vs") 
  return(y)
}


logL.sgs<- function(p, y, gg, model="GRW")
## log-likelihood of sgs model
{
  yl<- split4punc(y, gg)
  if (model=="GRW") {
  	ms<- p[1] 
  	vs<- p[2]
  	th<- p[3:4]
  	om<- p[5:6] }
  else {
  	vs<- p[1]
  	th<- p[2:3]
  	om<- p[4:5]	}
  
  l1<- logL.Stasis(p=c(th[1], om[1]), yl[[1]])
  if (model=="URW")	l2<- logL.URW(p=vs, yl[[2]])  
  else				l2<- logL.GRW(p=c(ms, vs), yl[[2]])
  l3<- logL.Stasis(p=c(th[2], om[2]), yl[[3]])
  
  logl<- l1+l2+l3
  return(logl)  	
}

logL.sgs.omega<- function(p, y, gg, model="GRW")
## log-likelihood of sgs model, omega shared over stasis segments
{
  yl<- split4punc(y, gg)
  if (model=="GRW") {
  	ms<- p[1]
  	vs<- p[2]
  	th<- p[3:4]
  	om<- p[5]	}
  else {
  	vs<- p[1]
  	th<- p[2:3]
  	om<- p[4]	}
  
  l1<- logL.Stasis(p=c(th[1], om), yl[[1]])
  if (model=="URW")	l2<- logL.URW(p=vs, yl[[2]])
  else				l2<- logL.GRW(p=c(ms, vs), yl[[2]])
  l3<- logL.Stasis(p=c(th[2], om), yl[[3]])
  
  logl<- l1+l2+l3
  return(logl)  		
}


opt.sgs<- function(y,gg,cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE, oshare=TRUE, model="GRW")
# do optimization for sgs model
{
 yl<- split4punc(y, gg)
 hat<- mle.GRW(yl[[2]])
 if (hat[1] == 0)	hat[1]<- 1e-3
 if (hat[2] < 1e-4)	hat[2]<- 1e-4
 if (model=="URW")	p0<- c(hat[2], mean(yl[[1]]$mm), mean(yl[[3]]$mm)) 
 else 				p0<- c(hat, mean(yl[[1]]$mm), mean(yl[[3]]$mm)) 
 
 if (oshare)	
 	{ p0<- append(p0, mean(var(yl[[1]]$mm), var(yl[[3]]$mm)))
 	  if (model=="GRW")		{ K<-7; pn<- c("ms","vs","theta1","theta2","omega"); lw<- c(NA,0,NA,NA,0) }
 	  else					{ K<-6; pn<- c("vs","theta1","theta2","omega"); lw<- c(0,NA,NA,0) }	}  
 else
 	{ p0<- append(p0, c(var(yl[[1]]$mm), var(yl[[3]]$mm)) )
 	  if (model=="GRW")	{ K<- 8; pn<- c("ms","vs","theta1","theta2","omega1","omega2"); lw<- c(NA,0,NA,NA,0,0) }
 	  else 				{ K<- 7; pn<- c("vs","theta1","theta2","omega1","omega2"); lw<- c(0,NA,NA,0,0) }
 	}
 
 cl$ndeps <- p0/100
 names(p0)<- pn

 if (oshare)
  {
   if (meth == "L-BFGS-B") w <- optim(p0, fn=logL.sgs.omega, gg=gg, method=meth, lower=lw, control=cl, hessian=hess, y=y, model=model)
   else w <- optim(p0, fn=logL.sgs.omega, gg=gg, method = meth, control = cl, hessian=hess, y=y)
  }
  else
  {
   if (meth == "L-BFGS-B")  w <- optim(p0, fn=logL.sgs, gg=gg, method= meth, lower=lw, control=cl, hessian=hess, y=y, model=model)
   else w <- optim(p0, fn = logL.sgs, gg=gg, method = meth, control = cl, hessian = hess, y = y)
  }

  # add more information to results
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else 			w$se<- NULL

  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste('SGS', model, sep='-'), method='AD', K=K, n=length(y$mm)-1, se=w$se)  
  return(wc) 
}

fit.sgs<- function(y, minb=7, oshare=TRUE, pool=TRUE, silent=FALSE, hess=FALSE, meth="L-BFGS-B", model="GRW")
## optimize for stasis-GRW-stasis dynamics (with some min n per section)
{
 if (pool)	y<- pool.var(y, ret.paleoTS=TRUE)  # pool variances
 ns<- length(y$mm)
 ng<-3
 GG<- shifts(ns, ng, minb=minb)

 nc<- ncol(GG)
 if (!silent)	cat ("\n\nTotal # hypotheses: ", nc, "\n\n", "i\tshifts\tlogL\n")
 wl<- list()
 logl<- array(-Inf, dim=nc)
 
 for (i in 1:nc) 
  {
    gg<- shift2gg(GG[,i], ns)
	w<- opt.sgs(y, gg, oshare=oshare, hess=hess, meth=meth, model=model)              
    if (!silent) 	cat (i, "\t", GG[,i], "\t", round(w$logL, 3), "\n")
    logl[i]<- w$logL
    wl[[i]]<- w
  }

 # add more information to results
 if (!silent) cat("\n")
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ss<- GG[,winner]
 names(ss)<- paste("shift", 1:(ng-1), sep='')
 ww$parameters<- append(ww$parameters, ss)
 ww$all.logl<- logl
 ww$GG<- GG
 return(ww)
}


sim.covTrack<- function (ns=20, b=1, evar=0.1, z, nn=rep(20, times=ns), tt=0:(ns-1), vp=1)
# simulates tracking optimum, with covariate z
{
  # check on length of covariate
  if(length(z)==ns)	
  	{ 
  	 z<- diff(z)
	 warning("Covariate z is same length of sequence (ns); using first difference of z as the covariate.")
	 }
  
  MM <- array(dim = ns)
  mm <- array(dim = ns)
  vv <- array(dim = ns)
  dt <- diff(tt)
  inc <- rnorm(ns-1, b*z, sqrt(evar))
  MM <- cumsum(c(0, inc))
  mm <- MM + rnorm(ns, 0, sqrt(vp/nn))  # add sampling error
  vv <- rep(vp, ns)
  gp <- c(b, evar)
  names(gp) <- c("slope.b", "evar")
  res <- as.paleoTS(mm = mm, vv = vv, nn = nn, tt = tt, MM = MM, 
        genpars = gp, label = "Created by sim.covTrack()")
  return(res)	
}




logL.covTrack<- function(p, y, z)
# z is covariate; pars = b (slope), evar (variance)
# IMPT: z is of length ns-1; one for each AD transition IN ORDER
{
  b<- p[1]
  evar<- p[2]
  dy <- diff(y$mm)
  dt <- diff(y$tt)
  nd <- length(dy)
  sv <- y$vv/y$nn
  svA <- sv[1:nd]
  svD <- sv[2:(nd + 1)]
  svAD <- svA + svD
  #S <- -0.5 * log(2 * pi * (evar + svAD)) - ((dy - (b*z))^2)/(2 * (evar + svAD))
  S<- dnorm(x=dy, mean=b*z, sd=sqrt(evar+svAD), log=TRUE)
  return(sum(S))	
}


opt.covTrack<- function (y, z, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE) 
{
    # check if z is of proper length; first difference if necessary
    ns<- length(y$mm)
    if(length(z)==length(y$mm))	
  	 { 
  	  z<- diff(z)
	  warning("Covariate z is same length of sequence (ns); using first difference of z as the covariate.")
	 }
	if (length(z) != ns-1)  stop("Covariate length [", length(z), "] does not match the sequence length [", ns, "]\n" )

    
    # get initial estimates by regression
    reg<- lm(diff(y$mm) ~ z-1)
    p0<- c(coef(reg), var(resid(reg)))
    names(p0) <- c("b", "evar")
    
    # pool variances if needed and do optimization
    if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
    if (is.null(cl$ndeps)) 
        cl$ndeps <- abs(p0/10000)
    cl$ndeps[cl$ndeps==0]<- 1e-8  ## will fail o.w. if any p0=0
    if (meth == "L-BFGS-B") 
        w <- optim(p0, fn=logL.covTrack, method = meth, lower = c(NA, 0), control = cl, hessian = hess, y=y, z=z)
    else  w<- optim(p0, fn=logL.covTrack, method=meth, control=cl, hessian=hess, y=y, z=z)
   

    if (hess)  w$se <- sqrt(diag(-1 * solve(w$hessian)))
    else w$se <- NULL
    
	wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='TrackCovariate', method='AD', K=2, n=length(y$mm)-1, se=w$se)
	return(wc)
}

logL.Mult.covTrack<- function (p, yl, zl)
# y is a _list_ of paleoTS objects, and z is a _list_ of covariates
{
  Smult<- 0
  nseq<- length(yl)
  for (i in 1:nseq)
    {	Smult<- Smult + logL.covTrack(p, yl[[i]], zl[[i]]) }
  return (Smult)	
}


opt.covTrack.Mult<- function (yl, zl, cl=list(fnscale=-1), pool=TRUE, hess=FALSE)
# y and z are lists of paleoTS, and covariates, respectively
{
 if (class(yl) == "paleoTS") 
      { stop("opt.covTrack.Mult is only for multiple paleoTS sequences\n") }
 nseq <- length(yl)
 if (pool) {
        for (i in 1:nseq) yl[[i]] <- pool.var(yl[[i]], ret.paleoTS = TRUE) }
        
 # check lengths of z
 for (i in 1:nseq)
 {
    ns<- length(yl[[i]]$mm)
    if(length(zl[[i]])==length(yl[[i]]$mm))	
  	 { 
  	  zl[[i]]<- diff(zl[[i]])
	  warning("Covariate z is same length of sequence (ns); using first difference of z as the covariate.")
	 }
	if (length(zl[[i]]) != ns-1)  stop("Covariate length [", length(zl[[i]]), "] does not match the sequence length [", ns, "]\n" )
 	
 }

        
 regMat<- array(dim=c(nseq, 3))  # hold results of regressions for sequences separately
 colnames(regMat)<- c("b0", "b1", "evar")
 for(i in 1:nseq)
 	{
 		w<- lm(diff(yl[[i]]$mm) ~ zl[[i]])
 		aa<- anova(w)
 		regMat[i,]<- c(coef(w), aa$Sum[2]/aa$Df[2])  # estimates of b0, b1, and evar
 	}
 p0<- apply(regMat, 2, median)  # initial estimates are median of separate regressions 
 w<- optim(p0, fn=logL.Mult.covTrack, control=cl, method="L-BFGS-B", lower=c(NA,0), hessian=hess, yl=yl, zl=zl)
  
 if (hess) 	w$se<- sqrt(diag(-1 * solve(w$hessian)))
 else		w$se<- NULL
 
 ff<- function(x) length(x$mm)-1
 n<- sum(sapply(yl, ff))
 
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='trackCovariate.Mult', method='AD', K=2, n=n, se=w$se)
 return(wc)	
}

logL.Mult.joint.covTrack<- function (p, yl, zl)
# y is a _list_ of paleoTS objects, and z is a *list* of covariates
{
  Smult<- 0
  nseq<- length(yl)
  for (i in 1:nseq)
    {	Smult<- Smult + logL.joint.covTrack(p, yl[[i]], zl[[i]]) }
  return (Smult)	
}


opt.joint.covTrack.Mult<- function (yl, zl, cl=list(fnscale=-1), pool=TRUE, hess=FALSE)
# y and z are lists of paleoTS, and covariates, respectively
{
 if (class(yl) == "paleoTS") 
      { stop("opt.covTrack.Mult is only for multiple paleoTS sequences\n") }
 nseq <- length(yl)
 if (pool) {
        for (i in 1:nseq) yl[[i]] <- pool.var(yl[[i]], ret.paleoTS = TRUE) }
        
 # check lengths of z
 for (i in 1:nseq)
 {
    ns<- length(yl[[i]]$mm)
	if (length(zl[[i]]) != ns)  stop("Covariate length [", length(zl[[i]]), "] does not match the sequence length [", ns, "]\n" )
 	
 }
        
 regMat<- array(dim=c(nseq, 3))  # hold results of regressions for sequences separately
 colnames(regMat)<- c("b0", "b1", "evar")
 for(i in 1:nseq)
 	{
 		w<- lm(yl[[i]]$mm ~ zl[[i]])
 		aa<- anova(w)
 		regMat[i,]<- c(coef(w), aa$Sum[2]/aa$Df[2])  # estimates of b0, b1, and evar
 	}
 p0<- apply(regMat, 2, median)  # initial estimates are median of separate regressions
 w<- optim(p0, fn=logL.Mult.joint.covTrack, control=cl, method="L-BFGS-B", lower=c(NA,NA,0), hessian=hess, yl=yl, zl=zl)
 
 if (hess) 	w$se<- sqrt(diag(-1 * solve(w$hessian)))
 else		w$se<- NULL
 
 ff<- function(x) length(x$mm)
 n<- sum(sapply(yl, ff))
 
 wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='trackCovariate.Mult', method='Joint', K=3, n=n, se=w$se)
 return(wc)	
}



logL.joint.covTrack<- function(p, y, z)
# z is covariate; pars = b0 (intercept), b1 (slope), evar (variance)
# z is of length ns; one for each sample in the time-series
{
  b0<- p[1]
  b1<- p[2]
  evar<- p[3]
  sv <- y$vv/y$nn
  S<- dnorm(x=y$mm, mean=b0 + b1*z, sd=sqrt(evar+sv), log=TRUE)
  return(sum(S))	
}

opt.joint.covTrack<- function (y, z, pool=TRUE, cl=list(fnscale=-1), meth="L-BFGS-B", hess=FALSE) 
{
    # check if z is of proper length; first difference if necessary
    ns<- length(y$mm)
	if (length(z) != ns)  stop("Covariate length [", length(z), "] does not match the sequence length [", ns, "]\n" )
    
    # check if both covariate and trait are trended; if so, warn that this approach is probably unreliable
    r.trait<- cor(y$mm, y$tt)
    r.cov<- cor(z, y$tt)
    mess<- c("Both the trait and covariate are strongly trended.  The joint approach is probably unreliable in this situation; consider using the AD parameterization instead.")
    if(abs(r.trait) > 0.6 && abs(r.cov > 0.6))	warning(mess)
    
    
    # get initial estimates by regression
    reg<- lm(y$mm ~ z)
    p0<- c(coef(reg), var(resid(reg)))
    names(p0) <- c("b0", "b1", "evar")
    
    # pool variances if needed and do optimization
    if (pool) y <- pool.var(y, ret.paleoTS = TRUE)
    if (is.null(cl$ndeps)) 
        cl$ndeps <- rep(1e-8, length(p0))
    if (meth == "L-BFGS-B") 
        w <- optim(p0, fn=logL.joint.covTrack, method = meth, lower = c(NA,NA,0), control = cl, hessian = hess, y=y, z=z)
    else  w<- optim(p0, fn=logL.joint.covTrack, method=meth, control=cl, hessian=hess, y=y, z=z)
   

    if (hess)  w$se <- sqrt(diag(-1 * solve(w$hessian)))
    else w$se <- NULL
    
	wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName='TrackCovariate', method='Joint', K=3, n=length(y$mm), se=w$se)
	return(wc)
}



LRI<- function(x, gen.per.t=1e6, draw=TRUE)
{
  n<- length(x$mm)
  x$tt<- x$tt*gen.per.t  # convert to generational timescale (e.g., 1e6 when time is measured in Myr and generation time = 1 yr)
  sp<- sqrt(pool.var(x))
  dy<- outer(x$mm, x$mm, FUN='-')  # outer() much faster than looping
  dt<- outer(x$tt, x$tt, FUN='-')
  dy<- abs(dy)
  dt<- abs(dt)
   
  dy<- dy/sp
  ut<- upper.tri(dy)
  logI<- log10(dt[ut])
  logR<- log10(dy[ut]) - logI
   
  # subset only non-zero rates
  ok<- is.finite(logR)
  if (sum(ok) != length(logI))	warning("Some zero rates were ignored.")
  
  # function to fit line in LRI plot robustly, using min abs deviations as per Gingerich 1993
  opt.lad<- function (x, y)
	{
	  ok<- is.finite(x) & is.finite(y)
	  xok<- x[ok]
	  yok<- y[ok]
	  w.ls<- lm(yok ~ xok)  # use LS coef as starting point in optimization
	
	  # function minimized for least abs deviation
	  fad<- function(p, x, y)  sum( abs(y - p[2]*x - p[1]) )
	  w.lad<- optim(w.ls$coef, fn=fad, x=xok, y=yok)
	  return(w.lad$par)	
	}
    
  # call robust regression function
  lf<- opt.lad(logI[ok], logR[ok])
  lf[3]<- 10^lf[1]   # generational rate
  names(lf)<- c("Intercept", "slope", "GenerationalRate")
  
  # do LRI plot, if desired
  if (draw)
   {
   	xl<- c(0, max(logI[ok]))
   	yl<- c(min(logR[ok]), max(lf[1], max(logR[ok])))
   	plot(logI[ok], logR[ok], xlim=xl, ylim=yl, xlab='log10 Interval [generations]', ylab='log10 Rate', cex=0.6)
   	abline(lf[1], lf[2], col='black', lwd=2)
   	title("LRI plot")
   	mtext(paste('data label: ', x$lab), cex=0.6, font=3, col='darkgrey')
   	restext<- paste('Slope = ', round(lf[2],3), '\n', 'Intercept = ', round(lf[1],3), '\n', 'Generational Rate = ', round(lf[3],5), '\n', sep='')
   	text(0, min(logR[ok]), restext, adj=c(0,0), cex=0.7, font=2)
   }
  
  #w<- list(b0=b0, b1=b1, h0=h0, logR=logR, logI=logI, dy=dy)  
  return(lf)
}


## Complex mode functions from PNAS ##
logL.joint.URW.Stasis<- function(p, y, gg)
 {
 	# preliminaries
 	anc<- p[1]
 	vs<- p[2]
 	theta<- p[3]
 	omega<- p[4]
 	n<- length(y$mm)

	# get vector of means
	st.seg<- which(gg==2)
 	rw.seg<- which(gg==1)
 	tt.rw<- y$tt[rw.seg]
	M<- c(rep(anc, length(rw.seg)), rep(theta, length(st.seg)))
	M<- unname(M)	

 	# compute covariance matrix
 	VVst<- diag(omega, nrow=length(st.seg))
 	VVrw<- vs * outer(tt.rw, tt.rw, FUN=pmin)
 	VVtot<- array(0, dim=c(n,n))
 	VVtot[st.seg, st.seg]<-  VVst
 	VVtot[rw.seg, rw.seg]<-  VVrw
 	diag(VVtot)<- diag(VVtot) + y$vv/y$nn  # add sampling error
 	#cat(round(p,3), "\t[")
 	
	# logL from mvn
	S<- dmnorm(y$mm, mean=M, varcov=VVtot, log=TRUE)
	#S<- dmvnorm(y$mm, mean=M, sigma=VVtot, log=TRUE)
	#cat(S, "]\n")
 	return(S)
 	#return(VVtot)
 }
 
 

 logL.joint.GRW.Stasis<- function(p, y, gg)
 {
 	# preliminaries
 	anc<- p[1]
 	ms<- p[2]
 	vs<- p[3]
 	theta<- p[4]
 	omega<- p[5]
 	n<- length(y$mm)

	# get vector of means
 	#ts<- max(y$tt[gg==1])  # time of shift point
	st.seg<- which(gg==2)
 	rw.seg<- which(gg==1)
 	tt.rw<- y$tt[rw.seg]
	M<- c(rep(anc, sum(gg==1)) + ms*y$tt[gg==1], rep(theta, length(st.seg)))
	M<- unname(M)	

 	# compute covariance matrix
 	VVst<- diag(omega, nrow=length(st.seg))
 	VVrw<- vs * outer(tt.rw, tt.rw, FUN=pmin)
 	VVtot<- array(0, dim=c(n,n))
 	VVtot[st.seg, st.seg]<-  VVst
 	VVtot[rw.seg, rw.seg]<-  VVrw
 	diag(VVtot)<- diag(VVtot) + y$vv/y$nn  # add sampling error
 	#cat(round(p,3), "\t[")

 	
	# logL from mvn
	S<- dmnorm(y$mm, mean=M, varcov=VVtot, log=TRUE)
	#S<- dmvnorm(y$mm, mean=M, sigma=VVtot, log=TRUE)
	#cat(S, "]\n")
 	return(S)
 	#return(VVtot)
 }
 


 
 opt.joint.RW.Stasis<- function(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
 {
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
  rw.model<- match.arg(rw.model)
  
  # get initial parameter estimates
  small<- 1e-8
  if(rw.model=="URW")  {p0rw<- mle.URW(sub.paleoTS(y, ok=gg==1)); K<- 5} # assumes shift point is free parameter
  else				   {p0rw<- mle.GRW(sub.paleoTS(y, ok=gg==1)); K<- 6}
  p0st<- mle.Stasis(sub.paleoTS(y, ok=gg==2))
  if(p0rw["vstep"] <= small)	p0rw["vstep"]<- 100*small
  if(p0st["omega"] <= small) 	p0st["omega"]<- 100*small
  p0anc<- y$mm[1]
  names(p0anc)<- "anc"
  p0<- c(p0anc, p0rw, p0st)
  #print(p0)
    
  #cl$parscale <-
  ll.urw<- c(NA,small,NA,small)
  ll.grw<- c(NA,NA,small,NA,small)
  if(rw.model=="URW")	{ ll<- ll.urw}
  else					{ ll<- ll.grw}
  if(meth!="L-BFGS-B")	ll<- NULL  # sets so to determine meth

  if(rw.model=="URW")	{
  		w <- try(optim(p0, fn=logL.joint.URW.Stasis, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)  
  		if(class(w)=="try-error"){
  				cl<- list(fnscale=-1, parscale=c(1,100,1,10))
  				w <- try(optim(p0, fn=logL.joint.URW.Stasis, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)  
  				}			
  } else if(rw.model=="GRW"){	
  		w <- try(optim(p0, fn=logL.joint.GRW.Stasis, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y) , silent=TRUE)
  		if(class(w)=="try-error"){
  				cl<- list(fnscale=-1, parscale=c(1,10,100,1,10))
 				w <- try(optim(p0, fn=logL.joint.GRW.Stasis, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y) , silent=TRUE)
 				}}
  # add more information to results
  if(class(w)=="try-error")	{
  		wc<- as.paleoTSfit(logL=NA, parameters=NA, modelName=paste(rw.model, "Stasis", sep='-'), method='Joint', K=K, n=length(y$mm), se=NULL)
  		return(wc)
  		}
  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste(rw.model, "Stasis", sep='-'), method='Joint', K=K, n=length(y$mm), se=w$se)

  return(wc)	 	
 }
 
 opt.AD.RW.Stasis<- function(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
 {
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
  rw.model<- match.arg(rw.model)
  
  # split and optimize AD separately for segments
  yl<- split4punc(y, gg)
  if(rw.model=="URW")	{ w1<- opt.URW(yl[[1]], pool=pool, meth=meth, hess=hess); K<- 4}
  else 					{ w1<- opt.GRW(yl[[1]], pool=pool, meth=meth, hess=hess); K<- 5}
  w2<- opt.Stasis(yl[[2]], pool=pool, meth=meth, hess=hess)
  
  
   # add more information to results
  if (hess)	{	se1<- sqrt(diag(-1*solve(w1$hessian)))
  				se2<- sqrt(diag(-1*solve(w2$hessian)))
  				se<- c(se1, se2)  } 					
  else			se<- NULL
  wc<- as.paleoTSfit(logL=w1$logL+w2$logL, parameters=c(w1$par, w2$par), modelName=paste(rw.model, "Stasis", sep='-'), method='AD', K=K, n=length(y$mm)-1, se=se)

  return(wc)	 	
 }
 
 
  logL.joint.Stasis.URW<- function(p, y, gg)
 {
 	# preliminaries
 	#cat(round(p,6), "\t[")
 	theta<- p[1]
 	omega<- p[2]
 	vs<- p[3]
 	n<- length(y$mm)

	# get vector of means
	M<- rep(theta, n)
	M<- unname(M)	

 	# compute covariance matrix
 	st.seg<- which(gg==1)
 	VVst<- diag(omega, nrow=length(st.seg))
 	rw.seg<- which(gg==2)
 	tt.rw<- y$tt[rw.seg] - y$tt[rw.seg[1]-1]
 	VVrw<- vs * outer(tt.rw, tt.rw, FUN=pmin)
 	VVtot<- array(0, dim=c(n,n))
 	VVtot[st.seg, st.seg]<- VVst
 	VVtot[rw.seg, rw.seg]<- VVrw
 	diag(VVtot)<- diag(VVtot) + y$vv/y$nn  # add sampling error
 	
	# logL from mvn
	S<- dmnorm(y$mm, mean=M, varcov=VVtot, log=TRUE)
	#S<- dmvnorm(y$mm, mean=M, sigma=VVtot, log=TRUE)
	#cat(S, "]\n")
 	return(S)
 	#return(VVtot)
 }
 
   logL.joint.Stasis.GRW<- function(p, y, gg)
 {
 	# preliminaries
 	theta<- p[1]
 	omega<- p[2]
 	ms<- p[3]
 	vs<- p[4]
 	n<- length(y$mm)

	# get vector of means
 	st.seg<- which(gg==1)
 	rw.seg<- which(gg==2)
 	#ts<- max(y$tt[gg==1])  # time of shift point
 	tt.rw<- y$tt[rw.seg] - y$tt[rw.seg[1]-1]
	M<- c(rep(theta, length(st.seg)), theta + ms*tt.rw)
	M<- unname(M)	

 	# compute covariance matrix
 	VVst<- diag(omega, nrow=length(st.seg))
 	VVrw<- vs * outer(tt.rw, tt.rw, FUN=pmin)
 	VVtot<- array(0, dim=c(n,n))
 	VVtot[st.seg, st.seg]<-  VVst
 	VVtot[rw.seg, rw.seg]<- VVrw
 	diag(VVtot)<- diag(VVtot) + y$vv/y$nn  # add sampling error
 	#cat(round(p,3), "\t[")
 	
	# logL from mvn
	S<- dmnorm(y$mm, mean=M, varcov=VVtot, log=TRUE)
	#S<- dmvnorm(y$mm, mean=M, sigma=VVtot, log=TRUE)
	#cat(S, "]\n")
 	return(S)
 	#return(VVtot)
 }

 
 opt.joint.Stasis.RW<- function(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
 {
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
  rw.model<- match.arg(rw.model)
  
  # get initial parameter estimates
  small<- 1e-8
  if(rw.model=="URW")  {p0rw<- mle.URW(sub.paleoTS(y, ok=gg==2)); K<- 4} # assumes shift point is free parameter
  else				   {p0rw<- mle.GRW(sub.paleoTS(y, ok=gg==2)); K<- 5}
  p0st<- mle.Stasis(sub.paleoTS(y, ok=gg==1))
  if(p0rw["vstep"] <= small)	p0rw["vstep"]<- 100*small
  if(p0st["omega"] <= small) 	p0st["omega"]<- 100*small
  p0<- c(p0st, p0rw)
  #cat(p0, "\n\n")  
    
  
  ll.urw<- c(NA,small,small)
  ll.grw<- c(NA,small,NA,small)
  if(rw.model=="URW")	{ ll<- ll.urw}
  else					{ ll<- ll.grw}
  if(meth!="L-BFGS-B")	ll<- -Inf  # sets so to determine meth


  if(rw.model=="URW")	{	
  		w <- try(optim(p0, fn=logL.joint.Stasis.URW, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)  
  		if(class(w)=="try-error"){
  				cl<- list(fnscale=-1, parscale=c(1,10,100))
		  		w <- try(optim(p0, fn=logL.joint.Stasis.URW, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)    				
  				}
  }else if(rw.model=="GRW")	{
  		w <- try(optim(p0, fn=logL.joint.Stasis.GRW, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)  
  		if(class(w)=="try-error"){
  				cl<- list(fnscale=-1, parscale=c(1,10,10,100))
		  		w <- try(optim(p0, fn=logL.joint.Stasis.GRW, gg=gg, method = meth, lower = ll, control=cl, hessian=hess, y=y), silent=TRUE)    				
 				}}
 				
  # add more information to results
  if(class(w)=="try-error")	{
  		wc<- as.paleoTSfit(logL=NA, parameters=NA, modelName=paste("Stasis", rw.model, sep='-'), method="Joint", K=K, n=length(y$mm), se=NULL)
  		return(wc)
  		}

  if (hess)		w$se<- sqrt(diag(-1*solve(w$hessian)))
  else			w$se<- NULL
  wc<- as.paleoTSfit(logL=w$value, parameters=w$par, modelName=paste("Stasis", rw.model, sep='-'), method='Joint', K=K, n=length(y$mm), se=w$se)

  return(wc)	 	
 }
 
  opt.AD.Stasis.RW<- function(y, gg, rw.model=c("URW", "GRW"), cl=list(fnscale=-1), pool=TRUE, meth="L-BFGS-B", hess=FALSE)
 {
  if (pool) {y<- pool.var(y, ret.paleoTS=TRUE) } # pool variances in sequences
  rw.model<- match.arg(rw.model)
  
  # split and optimize AD separately for segments
  yl<- split4punc(y, gg)
  w1<- opt.Stasis(yl[[1]], pool=pool, meth=meth, hess=hess)
  if(rw.model=="URW")	{ w2<- opt.URW(yl[[2]], pool=pool, meth=meth, hess=hess); K<- 4}
  else 					{ w2<- opt.GRW(yl[[2]], pool=pool, meth=meth, hess=hess); K<- 5}
 
  
  
   # add more information to results
  if (hess)	{	se1<- sqrt(diag(-1*solve(w1$hessian)))
  				se2<- sqrt(diag(-1*solve(w2$hessian)))
  				se<- c(se1, se2)  } 					
  else			se<- NULL
  wc<- as.paleoTSfit(logL=w1$logL+w2$logL, parameters=c(w1$par, w2$par), modelName=paste("Stasis", rw.model, sep='-'), method='AD', K=K, n=length(y$mm)-1, se=se)

  return(wc)	 	
 }

 

 fitModeShift<- function(y, minb=7, pool=TRUE, order=c("Stasis-RW", "RW-Stasis"), rw.model=c("URW","GRW"), method=c('Joint', 'AD'), silent=FALSE, hess=FALSE,...)
## optimize models (with some min n per section) that shift from GRW/URW to Stasis, or vice versa
{
 method<- match.arg(method) 
 order<- match.arg(order)
 rw.model<- match.arg(rw.model)
  
 ns<- length(y$mm)
 ng<- 2
 GG<- shifts(ns, ng, minb=minb)

 nc<- ncol(GG)
 if (!silent)	cat ("Total # hypotheses: ", nc, "\n")
 wl<- list()
 logl<- array(-Inf, dim=nc)
 
 for (i in 1:nc)
  {
    if (!silent) 	cat (i, " ")
    # the different gg for AD and J is required for the interpretation of the "shift" parameters to be the same across parameterizations
    gg<- shift2gg(GG[,i], ns)

    if(method=='AD') {
    		if(order=="Stasis-RW")  		w<- opt.AD.Stasis.RW(y, gg, rw.model=rw.model, pool=pool, hess=hess, ...)
    		else if (order=="RW-Stasis")	w<- opt.AD.RW.Stasis(y, gg, rw.model=rw.model, pool=pool, hess=hess, ...)
    	}
    if (method=='Joint'){
    		if(order=="Stasis-RW")  		w<- opt.joint.Stasis.RW(y, gg, rw.model=rw.model, pool=pool, hess=hess, ...)
    		else if (order=="RW-Stasis")	w<- opt.joint.RW.Stasis(y, gg, rw.model=rw.model, pool=pool, hess=hess, ...)
		}    			

    logl[i]<- w$logL
    wl[[i]]<- w
  }
 if(!silent)	cat("\n")
 winner<- which.max(logl)
 ww<- wl[[winner]]
 ss<- GG[,winner]
 names(ss)<- paste("shift", 1:(ng-1), sep='')
 ww$parameters<- append(ww$parameters, ss)
 ww$all.logl<- logl
 ww$GG<- GG
 return(ww)
} 


bootSimpleComplex<- function(y, simpleFit, complexFit, nboot=99, minb=7, ret.full.distribution=FALSE, parallel=FALSE, ...)
# function to to parametric bootstrapping to test if complex model is better than a simple model
# AICc seems often too liberal in favoring complex models
# consider comparing simpleFit to fit of all complex models?
{
	# make sure models to be tested are available
	sName<- simpleFit$modelName
	cName<- complexFit$modelName
	sAvailModels<- c("StrictStasis", "Stasis", "URW", "GRW")
	cAvailModels<- c("Punc-1", "URW-Stasis", "GRW-Stasis", "Stasis-URW", "Stasis-GRW")
	if (! (sName %in% sAvailModels) )	stop(paste("Simple model [", sName, "]", "not among available models [", sAvailModels, "]"))
	if (! (cName %in% cAvailModels) )	stop(paste("Complex model [", cName, "]", "not among available models [", cAvailModels, "]"))
	
	
	# compute observed LR
	LR<- -2*(simpleFit$logL - complexFit$logL)

	# generate bootstrap samples, fit simple and complex
	xboot<- list()
	pp<- simpleFit$par
	ns<- length(y$mm)
	vp<- pool.var(y)

   if (parallel==TRUE) {
      cores<-parallel::detectCores()
      cl<-parallel::makeCluster(cores-1)
      doParallel::registerDoParallel(cl)

      wl<- foreach (i = 1:nboot, .packages=c('paleoTS')) %dopar% {
		
		# create simulated dataset, fit simple model
		if (sName == "GRW"){
				xboot<- sim.GRW(ns=ns, ms=pp['mstep'], vs=pp['vstep'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.GRW(xboot, ...)	}
		if (sName == "URW") {
				xboot<- sim.GRW(ns=ns, ms=0, vs=pp['vstep'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.URW(xboot, ...)	}
		if (sName == "Stasis"){
				xboot<- sim.Stasis(ns=ns, theta=pp['theta'], omega=pp['omega'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.Stasis(xboot, ...)	}				
		if (sName == "StrictStasis") {
				xboot<- sim.Stasis(ns=ns, theta=pp['theta'], omega=0, vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.StrictStasis(xboot, ...)	}				

		# fit complex model
		if(cName=="Punc-1")	compBFit<- fitGpunc(xboot, ng=2, method="Joint", silent=TRUE, ...)
		if(cName=="URW-Stasis") compBFit<- fitModeShift(xboot, order="RW-Stasis", rw.model="URW", method="Joint", silent=TRUE, ...)
		if(cName=="GRW-Stasis") compBFit<- fitModeShift(xboot, order="RW-Stasis", rw.model="GRW", method="Joint", silent=TRUE, ...)
		if(cName=="Stasis-URW") compBFit<- fitModeShift(xboot, order="Stasis-RW", rw.model="URW", method="Joint", silent=TRUE, ...)
		if(cName=="Stasis-GRW") compBFit<- fitModeShift(xboot, order="Stasis-RW", rw.model="GRW", method="Joint", silent=TRUE, ...)

		bres<- list(simpBFit=simpBFit, compBFit=compBFit, xboot=xboot)
		bres
      }
      parallel::stopCluster(cl)	      # kill cluster
	}
	if(parallel==FALSE){
	      wl<- foreach (i = 1:nboot, .packages=c('paleoTS')) %do% {
		
		# create simulated dataset, fit simple model
		if (sName == "GRW"){
				xboot<- sim.GRW(ns=ns, ms=pp['mstep'], vs=pp['vstep'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.GRW(xboot, ...)	}
		if (sName == "URW") {
				xboot<- sim.GRW(ns=ns, ms=0, vs=pp['vstep'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.URW(xboot, ...)	}
		if (sName == "Stasis"){
				xboot<- sim.Stasis(ns=ns, theta=pp['theta'], omega=pp['omega'], vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.Stasis(xboot, ...)	}				
		if (sName == "StrictStasis") {
				xboot<- sim.Stasis(ns=ns, theta=pp['theta'], omega=0, vp=vp, nn=y$nn, tt=y$tt)
				xboot$vv<- y$vv
				simpBFit<- opt.joint.StrictStasis(xboot, ...)	}				

		# fit complex model
		if(cName=="Punc-1")	compBFit<- fitGpunc(xboot, ng=2, method="Joint", silent=TRUE, ...)
		if(cName=="URW-Stasis") compBFit<- fitModeShift(xboot, order="RW-Stasis", rw.model="URW", method="Joint", silent=TRUE, ...)
		if(cName=="GRW-Stasis") compBFit<- fitModeShift(xboot, order="RW-Stasis", rw.model="GRW", method="Joint", silent=TRUE, ...)
		if(cName=="Stasis-URW") compBFit<- fitModeShift(xboot, order="Stasis-RW", rw.model="URW", method="Joint", silent=TRUE, ...)
		if(cName=="Stasis-GRW") compBFit<- fitModeShift(xboot, order="Stasis-RW", rw.model="GRW", method="Joint", silent=TRUE, ...)

		bres<- list(simpBFit=simpBFit, compBFit=compBFit, xboot=xboot)
		bres
      }	
	}
	
	#return(wl)
	
	# compile information
	logL.simp<- sapply(wl, FUN=function(x) x$simpBFit$logL)
	logL.comp<- sapply(wl, FUN=function(x) x$compBFit$logL)
	LRnull<- -2* (logL.simp - logL.comp)
	cc<- is.finite(LRnull)	
	p.value<- (sum(LRnull[cc] > LR) +1)  / (length(LRnull[cc])+1)
	res<- list(LRobs=LR, p.value=p.value)
	
	if(ret.full.distribution){	
			res$nullLR<- LRnull[cc]	}
			#res$boot<- 	lapply(wl, FUN=function(x) x$xboot)	}
			
	return(res)	
}

ESD<- function(y, dt, model=c("GRW", "URW", "Stasis", "allThree"), method=c("Joint", "AD"), pool=TRUE, ...)
{
	model<- match.arg(model)
	method<- match.arg(method)
	
	if(model=="GRW"){
		if(method=="Joint")		w<- opt.joint.GRW(y, pool=pool, ...)
		else if (method=="AD")	w<- opt.GRW(y, pool=pool, ...)
		pp<- w$par
		esd<- pp["mstep"]^2 * dt*2 + pp["vstep"] * dt
	}
	if(model=="URW"){
		if(method=="Joint")		w<- opt.joint.URW(y, pool=pool, ...)
		else if (method=="AD")	w<- opt.URW(y, pool=pool, ...)
		pp<- w$par
		esd<- pp["vstep"] * dt
		
	}
	if(model=="Stasis"){
		if(method=="Joint")		w<- opt.joint.Stasis(y, pool=pool, ...)
		else if (method=="AD")	w<- opt.Stasis(y, pool=pool, ...)
		pp<- w$par
		esd<- 2*pp["omega"]
		
	}
	if(model=="allThree"){
		w<- fit3models(y, silent=TRUE, method=method, pool=pool)
		aw<- w$modelFits$Akaike.wt
		mn<- rownames(w$modelFits)
		
		esd.grw<- ESD(y, dt, model="GRW", method=method, pool=pool)
		esd.urw<- ESD(y, dt, model="URW", method=method, pool=pool)
		esd.sta<- ESD(y, dt, model="Stasis", method=method, pool=pool)
		
		esdV<- c(esd.grw, esd.urw, esd.sta)
		esd<- sum(esdV * aw)
	}
  
  return(unname(esd))
	
}








