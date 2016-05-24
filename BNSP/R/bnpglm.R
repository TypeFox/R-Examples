bnpglm <- function(formula,family,data,offset,sampler="slice",StorageDir,ncomp,sweeps,burn,thin=1,seed,prec,
                   V,Vdf,Mu.nu,Sigma.nu,Mu.mu,Sigma.mu,Alpha.xi,Beta.xi,Alpha.alpha,Beta.alpha,Turnc.alpha,
                   Xpred,offsetPred,...){
    # Match call
    call <- match.call()
    # Family
    family.indicator <- match(c(family),c("poisson","binomial","negative binomial","beta binomial","generalized poisson","com-poisson","hyper-poisson","ctpd"))
    if (is.na(family.indicator)){
        stop('family must be character, and one of "poisson", "binomial", "negative binomial", "beta binomial", "generalized poisson", "com-poisson", "hyper-poisson", "ctpd"')
    }
    # Data environment & design matrix
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    X<-model.matrix(formula,data=data)[,-1]
    # Dimensions
    n<-NROW(X)
    p<-NCOL(X)
    # Prior parameters
    if (missing(seed)) seed<-as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
    if (missing(V)) V<-diag(1,p)/p
    if (missing(Vdf)) Vdf<-p
    if (missing(Mu.nu)) Mu.nu<-rep(0,p)
    if (missing(Sigma.nu)) Sigma.nu<-diag(1,p)*100
    if (missing(Mu.mu) & p > 1) Mu.mu<-apply(X,2,mean)
    if (missing(Mu.mu) & p == 1) Mu.mu<-mean(X)
    if (missing(Sigma.mu) & p > 1) Sigma.mu<-diag((apply(X,2,max)-apply(X,2,min))^2)
    if (missing(Sigma.mu) & p == 1) Sigma.mu<-(max(X)-min(X))^2
    if (missing(Alpha.alpha)) Alpha.alpha<-2
    if (missing(Beta.alpha)) Beta.alpha<-4
    if (missing(Turnc.alpha)) Turnc.alpha<-0.25
    if (p > 1){
        xbar<-apply(X,2,mean)
        xsd<-apply(X,2,sd)
    } else if (p == 1){
        xbar<-mean(X)
        xsd<-sd(X)
    }
    # Family specific responses and offset terms
    if (family.indicator==1 | family.indicator==3 | family.indicator > 4){
        Y <- model.response(mf, "any")
        offset <- as.vector(model.offset(mf))
        if (!is.null(offset)){
           if (length(offset) != NROW(Y))
               stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                   length(offset), NROW(Y)), domain = NA)
        }
        if (missing(offset)) offset<-rep(1,n)
    } else if (family.indicator==2 | family.indicator==4){
        Y1 <- model.response(mf, "any")
        if (NCOL(Y1)==1){
	        if (any(Y1 < 0 | Y1 > 1)) stop("y values must be 0 <= y <= 1")
	        offset <- array(1,n)
	        Y<-Y1
	    } else if (NCOL(Y1) == 2){
            offset <- Y1[, 1] + Y1[, 2]
	        Y<-Y1[,1]
	      } else
	         stop(paste("For the binomial family, y must be",
			             "a vector of 0 and 1's or a 2 column",
			             "matrix where col 1 is no. successes",
			             "and col 2 is no. failures"))
    }
    # Family specific prior parameters
    if (family.indicator==1){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==1) stop(paste("For Poisson mixtures, argument Alpha.xi must be of length 1"))
        if (!missing(Alpha.xi)) if (Alpha.xi < 0) stop(paste("For Poisson mixtures, argument Alpha.xi must be positive"))
        if (missing(Alpha.xi)) Alpha.xi<-1.0
        if (!missing(Beta.xi)) if (!length(Beta.xi)==1) stop(paste("For Poisson mixtures, argument Beta.xi must be of length 1"))
        if (!missing(Beta.xi)) if (Beta.xi < 0) stop(paste("For Poisson mixtures, argument Beta.xi must be positive"))
        if (missing(Beta.xi)) Beta.xi<-0.1
    } else if (family.indicator==2){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==1) stop(paste("For Binomial mixtures, argument Alpha.xi must be of length 1"))
        if (!missing(Alpha.xi)) if (Alpha.xi < 0) stop(paste("For Binomial mixtures, argument Alpha.xi must be positive"))
        if (missing(Alpha.xi)) Alpha.xi<-1.0
        if (!missing(Beta.xi)) if (!length(Beta.xi)==1) stop(paste("For Binomial mixtures, argument Beta.xi must be of length 1"))
        if (!missing(Beta.xi)) if (Beta.xi < 0) stop(paste("For Binomial mixtures, argument Beta.xi must be positive"))
        if (missing(Beta.xi)) Beta.xi<-1.0
    } else if (family.indicator==3){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==2) stop(paste("For Negative Binimial mixtures, argument Alpha.xi must be of length 2"))
        if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >0) stop(paste("For Negative Binimial mixtures, vector Alpha.xi must have positive elements"))
        if (missing(Alpha.xi)) Alpha.xi<-c(1.0,1.0)
        if (!missing(Beta.xi)) if (!length(Beta.xi)==2) stop(paste("For Negative Binimial mixtures, argument Beta.xi must be of length 2"))
        if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For Negative Binimial mixtures, vector Beta.xi must have positive elements"))
        if (missing(Beta.xi)) Beta.xi<-c(0.1,0.1)
    } else if (family.indicator==4){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==2) stop(paste("For Beta Binimial mixtures, argument Alpha.xi must be of length 2"))
        if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >0) stop(paste("For Beta Binimial mixtures, vector Alpha.xi must have positive elements"))
        if (missing(Alpha.xi)) Alpha.xi<-c(1.0,1.0)
        if (!missing(Beta.xi)) if (!length(Beta.xi)==2) stop(paste("For Beta Binimial mixtures, argument Beta.xi must be of length 2"))
        if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For Beta Binimial mixtures, vector Beta.xi must have positive elements"))
        if (missing(Beta.xi)) Beta.xi<-c(0.1,0.1)
    } else if (family.indicator==5){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==2) stop(paste("For Generalized Poisson mixtures, argument Alpha.xi must be of length 2"))
        if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >0) stop(paste("For Generalized Poisson mixtures, vector Alpha.xi must have positive elements"))
        if (missing(Alpha.xi)) Alpha.xi<-c(1.0,1.0)
        if (!missing(Beta.xi)) if (!length(Beta.xi)==2) stop(paste("For Generalized Poisson mixtures, argument Beta.xi must be of length 2"))
        if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For Generalized Poisson mixtures, vector Beta.xi must have positive elements"))
        if (missing(Beta.xi)) Beta.xi<-c(0.1,1.0)
    } else if (family.indicator==6){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==2) stop(paste("For COM-Poisson mixtures, argument Alpha.xi must be of length 2"))
        if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >0) stop(paste("For COM-Poisson mixtures, vector Alpha.xi must have positive elements"))
        #if (missing(Alpha.xi)) Alpha.xi<-c(1.0,1.0)
        if (missing(Alpha.xi)) Alpha.xi<-c(1.0,0.5)
        if (!missing(Beta.xi)) if (!length(Beta.xi)==2) stop(paste("For COM-Poisson mixtures, argument Beta.xi must be of length 2"))
        if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For COM-Poisson mixtures, vector Beta.xi must have positive elements"))
        #if (missing(Beta.xi)) Beta.xi<-c(1,1)
        if (missing(Beta.xi)) Beta.xi<-c(0.1,0.5)
    } else if (family.indicator==7){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==2) stop(paste("For Hyper-Poisson mixtures, argument Alpha.xi must be of length 2"))
        if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >0) stop(paste("For Hyper-Poisson mixtures, vector Alpha.xi must have positive elements"))
        if (missing(Alpha.xi)) Alpha.xi<-c(1.0,0.5)
        if (!missing(Beta.xi)) if (!length(Beta.xi)==2) stop(paste("For Hyper-Poisson mixtures, argument Beta.xi must be of length 2"))
        if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For Hyper-Poisson mixtures, vector Beta.xi must have positive elements"))
        if (missing(Beta.xi)) Beta.xi<-c(0.1,0.5)
    } else if (family.indicator==8){
        if (!missing(Alpha.xi)) if (!length(Alpha.xi)==3) stop(paste("For CTPD mixtures, argument Alpha.xi must be of length 3"))
        if (!missing(Alpha.xi)) if (sum(Alpha.xi < 0) >1) stop(paste("For CTPD mixtures, vector Alpha.xi must have 2 positive elements"))
        if (missing(Alpha.xi)) Alpha.xi<-c(1.0,1.0,0.0)
        if (!missing(Beta.xi)) if (!length(Beta.xi)==3) stop(paste("For CTPD mixtures, argument Beta.xi must be of length 3"))
        if (!missing(Beta.xi)) if (sum(Beta.xi < 0) >0) stop(paste("For CTPD mixtures, vector Beta.xi must have positive elements"))
        if (missing(Beta.xi)) Beta.xi<-c(0.1,0.1,100)
    }
    #Precision
    if (!missing(prec)) if (length(prec) == 1) prec<-c(prec,prec,prec)
    if (missing(prec)) if (family.indicator < 5 | family.indicator == 6 | family.indicator == 7) prec <- c(10,10)
    if (missing(prec)) if (family.indicator == 5) prec <- c(10,0.02)
    if (missing(prec)) if (family.indicator == 8) prec <- c(10,10,10)
    # Predictions
    if (missing(Xpred)){
        npred <- 0
        Xpred <- 1
    } else{
        npred <- NROW(Xpred)
    }
    if (npred > 0){
           if (NCOL(Xpred) != p)
               stop(gettextf("Xpred matrix includes %d covariates, but it should include %d: the number covariates in the model",
                   NCOL(Xpred), p), domain = NA)
    }
    if ((!missing(offsetPred)) & (npred > 0)){
           if (length(offsetPred) != npred)
               stop(gettextf("number of prediction offsets is %d, but should equal %d: the number of prediction scenarios",
                   length(offsetPred), npred), domain = NA)
    }
    if (missing(offsetPred) & (npred > 0)) offsetPred <- rep(round(mean(offset)),npred)
    if (missing(offsetPred) & (npred == 0)) offsetPred <- 1.0
    allEqlI <- 0; if (length(unique(offsetPred)) == 1) allEqlI <- 1
    meanReg <- array(0,npred)
    modeReg <- array(0,npred)
    Q05Reg <- array(0,npred)
    Q10Reg <- array(0,npred)
    Q15Reg <- array(0,npred)
    Q20Reg <- array(0,npred)
    Q25Reg <- array(0,npred)
    Q50Reg <- array(0,npred)
    Q75Reg <- array(0,npred)
    Q80Reg <- array(0,npred)
    Q85Reg <- array(0,npred)
    Q90Reg <- array(0,npred)
    Q95Reg <- array(0,npred)
    # Family specific predictions
    if (family.indicator==1 | family.indicator==3 | family.indicator>4){
        maxy <- max(Y)+1000 # maxy <- max(Y)+max(10,floor(0.1*max(Y)))
    } else if (family.indicator==2 | family.indicator==4){
        maxy <- max(offsetPred)+1
    }
    denReg <- array(0,npred*maxy)
    denVar <- array(0,npred*maxy)
    #Sampler
    sampler.indicator <- match(sampler,c("slice","truncated"))
    if (is.na(sampler.indicator)){
        stop(c(sampler," 'sampler' not recognized"))
    }
    # Storage directory & files
    WF <- 1
    if (!missing(StorageDir)){
        StorageDir <- path.expand(StorageDir)
        ncharwd <- nchar(StorageDir)}
    if (!missing(StorageDir)) if (!(substr(StorageDir,ncharwd,ncharwd)=="/")) StorageDir <- paste(StorageDir,"/",sep="")
    if (!missing(StorageDir)) if (!file.exists(StorageDir)) dir.create(StorageDir)
    if (missing(StorageDir)) {WF <- 0; StorageDir <- paste(getwd(),"/",sep="")}
    on.exit(if (WF==0) file.remove(
    paste(StorageDir,"BNSP.Th.txt",sep=""),
    paste(StorageDir,"BNSP.Sigmah.txt",sep=""),
    paste(StorageDir,"BNSP.SigmahI.txt",sep=""),
    paste(StorageDir,"BNSP.nuh.txt",sep=""),
    paste(StorageDir,"BNSP.muh.txt",sep=""),
    paste(StorageDir,"BNSP.xih.txt",sep=""),
    paste(StorageDir,"BNSP.alpha.txt",sep=""),
    paste(StorageDir,"BNSP.compAlloc.txt",sep=""),
    paste(StorageDir,"BNSP.nmembers.txt",sep=""),
    paste(StorageDir,"BNSP.Updated.txt",sep=""),
    paste(StorageDir,"BNSP.MeanReg.txt",sep=""),
    paste(StorageDir,"BNSP.Q05Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q10Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q15Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q20Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q25Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q50Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q75Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q80Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q85Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q90Reg.txt",sep=""),
    paste(StorageDir,"BNSP.Q95Reg.txt",sep="")))
    #H-clustering
    dis<-dist(cbind(Y,X), method = "euclidean")
    hclst<-hclust(dis, method="ward.D")
    compA<-cutree(hclst, k=ncomp)-1
    #Call C
    out<-.C("OneResLtnt", as.integer(seed), as.double(unlist(c(X))), as.integer(Y), as.double(offset),
            as.integer(sweeps), as.integer(burn), as.integer(thin), as.integer(ncomp), as.integer(n), as.integer(p),
            as.double(V), as.double(Vdf),
            as.double(Mu.nu), as.double(Sigma.nu),
            as.double(Mu.mu), as.double(Sigma.mu),
            as.double(Alpha.xi),as.double(Beta.xi),
            as.double(Alpha.alpha),as.double(Beta.alpha),as.double(Turnc.alpha),
            as.double(xbar), as.double(xsd), as.double(sum(Y)/sum(offset)), as.double(prec),
            as.integer(family.indicator), as.integer(sampler.indicator),
            as.integer(npred), as.double(Xpred), as.double(offsetPred), as.integer(allEqlI), as.integer(maxy),
            as.double(meanReg), as.double(modeReg), as.double(Q05Reg), as.double(Q10Reg), as.double(Q15Reg),
            as.double(Q20Reg),as.double(Q25Reg),as.double(Q50Reg), as.double(Q75Reg), as.double(Q80Reg),
            as.double(Q85Reg),as.double(Q90Reg),as.double(Q95Reg), as.double(denReg),as.double(denVar),
            as.character(StorageDir),as.integer(WF),as.integer(compA))
    #Output
    location<-33
    meanReg <- out[[location+0]][1:npred]
    modeReg <- out[[location+1]][1:npred]
    Q05Reg <- out[[location+2]][1:npred]
    Q10Reg <- out[[location+3]][1:npred]
    Q15Reg <- out[[location+4]][1:npred]
    Q20Reg <- out[[location+5]][1:npred]
    Q25Reg <- out[[location+6]][1:npred]
    Q50Reg <- out[[location+7]][1:npred]
    Q75Reg <- out[[location+8]][1:npred]
    Q80Reg <- out[[location+9]][1:npred]
    Q85Reg <- out[[location+10]][1:npred]
    Q90Reg <- out[[location+11]][1:npred]
    Q95Reg <- out[[location+12]][1:npred]
    denReg <- matrix(out[[location+13]][1:(maxy*npred)],nrow=npred,ncol=maxy,byrow=TRUE)
    denVar <- matrix(out[[location+14]][1:(maxy*npred)],nrow=npred,ncol=maxy,byrow=TRUE)
    denVar <- denVar - denReg^2
    fit <- list(call=call,seed=seed,meanReg=meanReg,modeReg=modeReg,Q05Reg=Q05Reg,Q10Reg=Q10Reg,
    Q15Reg=Q15Reg,Q20Reg=Q20Reg,Q25Reg=Q25Reg,Q50Reg=Q50Reg,Q75Reg=Q75Reg,Q80Reg=Q80Reg,Q85Reg=Q85Reg,
    Q90Reg=Q90Reg,Q95Reg=Q95Reg,denReg=denReg,denVar=denVar)
    class(fit) <- 'bnp'
    return(fit)
}
