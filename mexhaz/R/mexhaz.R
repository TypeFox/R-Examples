mexhaz <- function(formula,data,expected=NULL,base=c("weibull","exp.bs","pw.cst"),degree=3,knots=NULL,bo.max=NULL,n.gleg=20,init=NULL,random=NULL,n.aghq=10,fnoptim=c("nlm","optim"),verbose=100,method="Nelder-Mead",iterlim=10000,print.level=1,...){

    time0 <- as.numeric(proc.time()[3])
    FALCenv <- environment()

    call <- match.call()
    base <- match.arg(base)
    fnoptim <- match.arg(fnoptim)

    mf <- match.call(expand.dots=FALSE)
    m <- match(c("formula","data"),names(mf),0L)
    if (m[1]==0){
        stop("The 'formula' argument is required...")
    }
    if (m[2]==0){
        stop("The 'data' argument is required...")
    }
    name.data <- paste(substitute(data),sep="")
    if (length(name.data)>1){
        name.data <- name.data[2]
    }

        if (base=="exp.bs" & !degree%in%c(1:3)){
        stop("This function can only be used to estimate log-hazards described by B-splines of degree 1 to 3...")
    }

    Frailty.Adapt <- function(nodes, nodessquare, logweights, clust, clustd, delta, expect, betal, betaL, A, var, muhatcond){
        obj <- .Call("FrailtyAdapt", nodes, nodessquare, logweights, clust, clustd, expect, betal, betaL, A, var, muhatcond)
    }
    IntBs1 <- function(x,nph,param,leint,whint,knots){
        .Call("IntBs1",x,nph,param,leint,whint,knots)
    }
    IntBs23 <- function(x,nph,timecat,param,deg,n,lw,matk,totk){
        .Call("IntBs23",x,nph,timecat,param,deg,n,lw,as.double(matk),as.double(totk))
    }
    IntPwCst <- function(nph,param,leint,lerem,whint){
        .Call("IntPwCst",nph,param,leint,lerem,whint)
    }

    # Functions for computing the part of the B-spline bases that depends only on the knots and can therefore be calculated only once at the beginning of the function (used in combination with the IntBs23Stat function to estimate the hazard and the cumulative hazard)

    if (base=="exp.bs" & degree%in%c(2:3)){

        # For cubic B-splines
        Transf3 <- function(vec.knots){
            Dim <- (length(vec.knots)-5)
            Res <- matrix(NA,4,Dim)
            for (i in 1:Dim){
                TempK <- vec.knots[i:(i+5)]
                Res[1,i] <- 1/((TempK[6]-TempK[3])*(TempK[5]-TempK[3])*(TempK[4]-TempK[3]))
                Res[2,i] <- 1/((TempK[5]-TempK[2])*(TempK[4]-TempK[2])*(TempK[4]-TempK[3]))
                Res[3,i] <- 1/((TempK[5]-TempK[2])*(TempK[5]-TempK[3])*(TempK[4]-TempK[3]))
                Res[4,i] <- 1/((TempK[4]-TempK[1])*(TempK[4]-TempK[2])*(TempK[4]-TempK[3]))
            }
            return(Res)
        }

        # For quadratic B-splines
        Transf2 <- function(vec.knots){
            Dim <- (length(vec.knots)-3)
            Res <- matrix(NA,2,Dim)
            for (i in 1:Dim){
                TempK <- vec.knots[i:(i+3)]
                Res[1,i] <- 1/((TempK[4]-TempK[2])*(TempK[3]-TempK[2]))
                Res[2,i] <- 1/((TempK[3]-TempK[1])*(TempK[3]-TempK[2]))
            }
            return(Res)
        }

        # Creating the points and weights of Gauss-Legendre quadrature (for computing the cumulative hazard when base="exp.bs" and degree in c(2:3))
        if (n.gleg<=0 | round(n.gleg,0)!=n.gleg){
            stop("The 'n.gleg' argument must be a strictly positive integer.")
        }
        gl <- gauss.quad(n=n.gleg,kind="legendre")
        gln <- gl$nodes
        lglw <- log(gl$weights)
    }

    mf <- mf[c(1L,m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- "na.pass"
    mf[[1L]] <- quote(stats::model.frame)
    tot.formula <- terms(formula,data=data,specials="nph")
    indNph <- attr(tot.formula,"specials")$nph
    if (!is.null(indNph)){
        nphterm <- attr(tot.formula,"variables")[[1+indNph]]
        nTerm <- deparse(nphterm[[2L]],width.cutoff=500L,backtick=TRUE)
        nTerm2 <- deparse(nphterm,width.cutoff=500L,backtick=TRUE)
        FormulaN <- as.formula(paste("~",nTerm))
    }
    else {
        nTerm2 <- 0
        FormulaN <- as.formula("~1")
    }
    mf$formula <- update(tot.formula,paste(".~.",nTerm2,sep="-"))
    FormulaF <- update(tot.formula,paste("NULL~.",nTerm2,sep="-"))
    Xlevel.formula <- terms(FormulaF,data=data)
    m <- eval(mf,parent.frame())
    if (nrow(m)==0){
        stop("No non-missing observations...")
    }

    Y <- model.extract(m,"response")
    if (!inherits(Y,"Surv")){
        stop("Response must be a Surv() object...")
    }
    Survtype <- attr(Y, "type")
    if ((ncol(Y)==2) && (Survtype!="right")){
        stop(paste("mexhaz does not support \"", Survtype, "\" type of censoring with (0, end] survival data ",
               sep = ""))
    }
    else if ((ncol(Y)!=2)){
        stop(paste("mexhaz does not support survival data in counting process format...",
               sep = ""))
    }
    time.obs <- Y[,1]   # Follow-up time
    status.obs <- Y[,2]   # Status variable

    Test.Status <- unique(status.obs[!is.na(status.obs)])
    if (!identical(Test.Status[order(Test.Status)],c(0,1))){
        stop("The event variable should take on the values 0 or 1...")
    }

    data.fix <- model.frame(FormulaF,data=data,na.action=na.pass)
    data.nph <- model.frame(FormulaN,data=data,na.action=na.pass)
    n.obs.tot <- dim(data)[1] # Total number of observations

    # Expected hazard (a vector of 0 is created if NULL)
    if (!is.null(expected)){
        lambda.pop <- data[,expected]
        if (min(lambda.pop,na.rm=TRUE)<0){
            stop("The expected hazard for some observations is negative!")
        }
        BoolExp <- 1
    }
    else {
        lambda.pop <- rep(0,n.obs.tot)
        BoolExp <- 0
    }

    # Remove observations containing missing values and perform several formatting operations on the data
    if (!is.null(random)){
        random.obs <- data[,random]
        Idx.NA.Rdm <- which(is.na(random.obs))
        if (length(Idx.NA.Rdm)>0){
            warning("Cluster information was missing for some observations. These observations were consequently removed...")
            random.obs <- random.obs[-Idx.NA.Rdm]
            time.obs <- time.obs[-Idx.NA.Rdm]
            status.obs <- status.obs[-Idx.NA.Rdm]
            lambda.pop <- lambda.pop[-Idx.NA.Rdm]
            data.fix <- data.fix[-Idx.NA.Rdm,,drop=FALSE]
            data.nph <- data.nph[-Idx.NA.Rdm,,drop=FALSE]
        }
        random.obs <- as.factor(random.obs)
        clust <- levels(random.obs)
        IdxRnd <- order(random.obs)
        # The dataset MUST be ordered by the levels of the clustering variable
        random.obs <- random.obs[IdxRnd]
        time.obs <- time.obs[IdxRnd]
        status.obs <- status.obs[IdxRnd]
        lambda.pop <- lambda.pop[IdxRnd]
        data.fix <- data.fix[IdxRnd,,drop=FALSE]
        data.nph <- data.nph[IdxRnd,,drop=FALSE]
    }
    TotData <- cbind(time.obs,status.obs,lambda.pop,data.fix,data.nph)
    Xlevels <- .getXlevels(Xlevel.formula,TotData)
    Idx.Non.NA <- which(apply(TotData,1,function(vec){sum(is.na(vec))})==0)
    if (length(Idx.Non.NA)<n.obs.tot){
        warning("Covariables information was missing for some observations. These observations were consequently removed...")
    }
    time.obs <- time.obs[Idx.Non.NA]
    n.obs <- length(time.obs)
    if (n.obs==0){
        stop("No non-missing values for some covariables...")
    }
    Idx.Time.Neg <- which(time.obs<0)
    if (length(Idx.Time.Neg)>0){
        stop("Some observations have a negative follow-up time...")
    }
    Idx.Time.0 <- which(time.obs==0)
    if (length(Idx.Time.0)>0){
        warning("Some observations had a follow-up time of 0. A value of 1/730.5 (half a day) was substituted. Please check to see if it is appropriate or deal with 0 follow-up time values before using the mexhaz function.")
    }
    time.obs[Idx.Time.0] <- 1/730.5
    max.time <- max(time.obs)
    Bo <- ifelse(is.null(bo.max),max.time,max(max.time,bo.max))
    status.obs <- status.obs[Idx.Non.NA]
    n.events <- sum(status.obs)
    data.fix <- data.fix[Idx.Non.NA,,drop=FALSE]
    fix.obs <- model.matrix(FormulaF,data=data.fix,drop.unused.levels=TRUE)
    names.fix <- colnames(fix.obs)[-1]
    data.nph <- data.nph[Idx.Non.NA,,drop=FALSE]
    nph.obs <- model.matrix(FormulaN,data=data.nph,drop.unused.levels=TRUE)
    nbtd <- dim(nph.obs)[2]
    names.nph <- colnames(nph.obs)[-1]
    if (sum(!names.nph%in%names.fix)){
        stop("Some variables in the 'nph()' sub-formula are not included in the formula...")
    }
    lambda.pop <- lambda.pop[Idx.Non.NA]

    # Remove unnecessary objects
    rm(TotData,data.fix,data.nph)

    # Creation of the different objects necessary to compute the hazards

    # Weibull hazard
    if (base=="weibull"){
        log.time.obs <- log(time.obs)
        n.td.base <- 1
        n.ntd <- dim(fix.obs)[2]
        if (n.ntd>1){
            which.ntd <- c(1,2+1:(n.ntd-1))
        }
        else {
            which.ntd <- 1
        }
        n.td.nph <- length(names.nph)
        if (n.td.nph>0){
            names.nph <- paste("Rho",names.nph,sep="*")
            which.td <- c(2,(2+(n.ntd-1))+1:n.td.nph)
        }
        else {
            which.td <- 2
        }
        fix.obs <- fix.obs[,-1,drop=FALSE]
        nph.obs <- nph.obs[,-1,drop=FALSE]
        param.names <- c("Lambda","Rho",names.fix,names.nph)
        n.par.fix <- n.td.base+n.ntd+n.td.nph
        param.init <- rep(0,n.par.fix)
        param.init[1:2] <- 0.1
    }

    # Hazard modelled by the exponential of a B-spline
    else if (base=="exp.bs"){

        # Baseline hazard-related objects
        if (!is.null(knots)){
            if (min(knots,na.rm=TRUE)<=0 | max(knots,na.rm=TRUE)>=Bo | is.na(sum(knots))){
                stop("The 'knots' argument should be a vector of values strictly between 0 and max(bo.max,time.max) where time.max is the maximum follow-up time in the dataset")
            }
            else {
                if (sum(abs(knots-knots[order(knots)]))>0){
                    warning("The B-spline interior knots were not given in increasing order. They were consequently re-ordered.")
                    knots <- knots[order(knots)]
                }
                if (length(unique(knots))<length(knots)){
                    warning("There were duplicate values in the vector of B-spline interior knots. Duplicate values were removed.")
                    knots <- unique(knots)
                }
            }
        }
        cuts <- c(0,knots,Bo)
        interv <- c(knots,Bo)-c(0,knots)
        time.cat <- cut(time.obs,breaks=cuts)
        time.cat <- as.numeric(time.cat)-1
        vec.knots <- c(rep(0,degree),knots,rep(Bo,degree))
        if (degree==3){
            MatK <- Transf3(vec.knots)
        }
        else if (degree==2){
            MatK <- Transf2(vec.knots)
        }
        n.td.base <- degree + length(knots)
        names.base <- paste("BS",degree,".",1:n.td.base,sep = "")

        # For non time-dependent effects
        n.ntd <- dim(fix.obs)[2]
        if (n.ntd>1){
            which.ntd <- c(1,(n.td.base+1)+1:(n.ntd-1))
        }
        else {
            which.ntd <- 1
        }

        # For time-dependent effects
        n.td.nph <- length(names.nph)*n.td.base
        if (n.td.nph>0){
            which.td <- c(1+1:n.td.base,(n.td.base+n.ntd)+1:n.td.nph)
        }
        else {
            which.td <- c(1+1:n.td.base)
        }
        names.nph <- unlist(sapply(names.nph,function(x){paste(x,names.base,sep="*")}))

        nph.obs <- t(nph.obs) # Matrix of time-dependent effects has to be transposed for use by the IntBs'xx'Stat functions
        param.names <- c("Intercept",names.base,names.fix,names.nph)
        n.par.fix <- n.td.base+n.ntd+n.td.nph
        param.init <- rep(0,n.par.fix)
        param.init[1:(n.td.base+1)] <- -1
    }

    # Hazard modelled by the exponential of a piecewise constant function
    else if (base=="pw.cst") {

        # Baseline hazard-related objects
        if (!is.null(knots)){
            if (min(knots,na.rm=TRUE)<=0 | max(knots,na.rm=TRUE)>=Bo | is.na(sum(knots))){
                stop("The 'knots' argument should be a vector of values strictly between 0 and max(bo.max,time.max) where time.max is the maximum follow-up time in the dataset")
            }
            else {
                if (sum(abs(knots-knots[order(knots)]))>0){
                    warning("The knots defining the intervals on which the hazard is constant were not given in increasing order. They were consequently re-ordered.")
                    knots <- knots[order(knots)]
                }
                if (length(unique(knots))<length(knots)){
                    warning("There were duplicate values in the vector of knots defining the intervals on which the hazard is constant. Duplicate values were removed.")
                    knots <- unique(knots)
                }
            }
        }
        cuts <- c(0,knots,Bo)
        interv <- c(knots,Bo)-c(0,knots)
        time.cat <- cut(time.obs,breaks=cuts)
        names.base <- levels(time.cat)
        time.cat <- as.numeric(time.cat)-1
        time.remain <- time.obs-cuts[time.cat+1]
        n.td.base <- length(knots)+1

        # For non time-dependent effects
        fix.obs <- fix.obs[,-which(colnames(fix.obs)%in%colnames(nph.obs)),drop=FALSE]
        names.fix <- colnames(fix.obs)
        if (!is.null(names.fix)){
            n.ntd <- dim(fix.obs)[2]
        }
        else {
            n.ntd <- 0
        }
        if (n.ntd>0){
            which.ntd <- c(n.td.base+1:n.ntd)
        }
        else {
            which.ntd <- NULL
        }

        # For time-dependent effects
        n.td.nph <- length(names.nph)*n.td.base
        if (n.td.nph>0){
            which.td <- c(1:n.td.base,(n.td.base+n.ntd)+1:n.td.nph)
        }
        else {
            which.td <- c(1:n.td.base)
        }
        if (length(names.nph)>0){
            names.nph <- as.vector(sapply(names.nph,function(x){paste(x,names.base,sep="*")}))
        }

        nph.obs <- t(nph.obs) # Matrix of time-dependent effects has to be transposed for use by the IntPwCstStat function

        param.names <- c(names.base,names.fix,names.nph)
        n.par.fix <- n.td.base+n.ntd+n.td.nph
        param.init <- rep(0,n.par.fix)
        param.init[1:n.td.base] <- -1
    }

    # Creation of objects related to the random effect
    n.clust <- 1
    n.rand <- 0
    if (!is.null(random)){
        n.rand <- 1
        random.obs <- random.obs[Idx.Non.NA]
        status.one <- which(status.obs==1)
        lambda.pop.delta <- lambda.pop[status.one]
        n.clust <- length(clust) # Number of clusters
        n.by.clust <- as.vector(table(random.obs)) # Size of each cluster
        n.by.clust.delta <- as.vector(table(random.obs[status.one]))
        Idx.clust.no.obs <- which(n.by.clust==0)
        if (length(Idx.clust.no.obs)>0){
            warning("Some clusters had no non-missing values for some covariables. These observations were consequently removed...")
            cat("The following clusters had no non-missing values for some covariables:\n")
            print(names(table(random.obs))[Idx.clust.no.obs])
            n.by.clust <- n.by.clust[-Idx.clust.no.obs]
            n.by.clust.delta <- n.by.clust.delta[-Idx.clust.no.obs]
            clust <- clust[-Idx.clust.no.obs]
            n.clust <- length(clust)
        }
        parent.cst.adj <- rep(0,n.clust)
        param.names <- c(param.names,paste(random," (sd)",sep=""))

        # Creation of the points and weights of Gauss-Hermite quadrature
        if (n.aghq<=0 | round(n.aghq,0)!=n.aghq){
            stop("The 'n.aghq' argument must be a strictly positive integer.")
        }
        gq <- gauss.quad(n=n.aghq,kind="hermite")
        x.H <- gq$nodes
        x.H.2 <- x.H^2
        log.rho.H <- log(gq$weights)
        log.rho.H[log.rho.H==-Inf] <- -.Machine$double.xmax
    }

    n.par <- n.td.base+n.ntd+n.td.nph+n.rand

    # Initial values
    if (!is.null(init)){
        if (length(init)!=n.par){
            print(data.frame(Time.Dep.baseline=n.td.base,Non.TD.Effects=n.ntd,TD.Effects=n.td.nph,Random.Effect=n.rand))
            stop("Wrong number of init parameters",call.=FALSE)
        }
    }
    else {
        init <- c(param.init,rep(0.1,n.rand))
    }

    # Creation of variables that will be modified by the Hazard and LL.Tot functions
    parent.neval <- 0
    parent.param <- init
    parent.ptd <- rep(99,(n.td.base+n.td.nph))
    parent.pntd <- rep(99,n.ntd)
    if (base=="exp.bs"){
        parent.logHaz <- Inf
        parent.logCum <- Inf
    }
    parent.logHazBeta <- rep(0,n.obs)
    parent.logCumBeta <- rep(0,n.obs)
    parent.logLik <- -.Machine$double.xmax

    # Functions that compute the hazards (the computed individual-specific hazards are stored in the main function environment)

    if (base=="weibull"){
        Hazard <- function(p.td.H,p.ntd.H,env=FALCenv){
            if (p.td.H[1]>0 & p.ntd.H[1]>0){
                log.p.LT.1 <- log(p.ntd.H[1]) + as.vector(fix.obs%*%p.ntd.H[-1])
                log.p.LT.2 <- log(p.td.H[1]) + as.vector(nph.obs%*%p.td.H[-1])
                l.lambda.beta <- log.p.LT.2 +log.time.obs*(exp(log.p.LT.2)-1)+log.p.LT.1
                l.Lambda.beta <- log.time.obs*exp(log.p.LT.2)+log.p.LT.1
                valtot <- sum(l.lambda.beta) + sum(l.Lambda.beta)
                Test <- sum((is.nan(valtot)) | (valtot==Inf))
            }
            else {
                Test <- 1
            }
            if (!Test){
                env$parent.logHazBeta <- l.lambda.beta
                env$parent.logCumBeta <- l.Lambda.beta
            }
            return(Test)
        }
    }

    else if (base=="exp.bs"){
        Hazard <- function(p.td.H,p.ntd.H,env=FALCenv){
            if (degree==1) {
                temp.H <- IntBs1(time.obs,nph=nph.obs,param=p.td.H,leint=interv,whint=time.cat,knots=vec.knots)
                l.lambda <- temp.H[1:n.obs]
                l.Lambda <- temp.H[(n.obs)+1:n.obs]
            }
            else {
                test1 <- (sum(abs(p.td.H-env$parent.ptd))>0)
                if (test1) {
                    env$parent.logHaz <- Inf
                    env$parent.logCum <- Inf
                    temp.H <- IntBs23(time.obs,nph=nph.obs,timecat=time.cat,param=p.td.H,deg=degree,gln,lglw,matk=MatK,totk=vec.knots)
                    valtot <- temp.H[2*n.obs+1]
                    test2 <- ((is.nan(valtot)) | (valtot==Inf))
                    if (!test2){
                        env$parent.logHaz <- temp.H[1:n.obs]
                        env$parent.logCum <- temp.H[(n.obs)+1:n.obs]
                    }
                }
                l.lambda <- env$parent.logHaz
                l.Lambda <- env$parent.logCum
            }
            if (n.ntd!=0){
                beta.x <- as.vector(fix.obs%*%p.ntd.H)
            }
            else {
                beta.x <- 0
            }
            valtot <- sum(l.lambda + l.Lambda + beta.x)
            Test <- ((is.nan(valtot)) | (valtot==Inf))
            if (!Test){
                env$parent.logHazBeta <- l.lambda + beta.x
                env$parent.logCumBeta <- l.Lambda + beta.x
            }
            return(Test)
        }
    }

    else if (base=="pw.cst"){
        Hazard <- function(p.td.H,p.ntd.H,env=FALCenv){
            temp.H <- IntPwCst(nph=nph.obs,param=p.td.H,leint=interv,lerem=time.remain,whint=time.cat)
            l.lambda <- temp.H[1:n.obs]
            l.Lambda <- temp.H[(n.obs)+1:n.obs]
            if (n.ntd!=0){
                beta.x <- as.vector(fix.obs%*%p.ntd.H)
            }
            else {
                beta.x <- 0
            }
            valtot <- sum(l.lambda + l.Lambda + beta.x)
            Test <- ((is.nan(valtot)) | (valtot==Inf))
            if (!Test){
                env$parent.logHazBeta <- l.lambda + beta.x
                env$parent.logCumBeta <- l.Lambda + beta.x
            }
            return(Test)
        }
    }

    # Function that controls what is printed during the optimisation procedure
    if (verbose>0){
        verbose.ll <- function(){
            if (!((iv <- parent.neval/verbose)-floor(iv))){
                time1 <- as.numeric(proc.time()[3])
                print(data.frame(Eval=parent.neval,LogLik=-parent.logLik,Time=time1-time0,row.names=""))
                cat("Param\n")
                print(round(parent.param,4))
                cat("\n")
            }
        }
    }
    else {
        verbose.ll <- function(){}
    }

    # Function that actually computes the log-likelihood

    LL.Tot <- function(p.LT,mu.hat.LT=0,env=FALCenv){

        p.td <- p.LT[which.td]
        p.ntd <- p.LT[which.ntd]

        if (Hazard(p.td,p.ntd)){
            res.LT <- .Machine$double.xmax
        }
        else {
            if (!is.null(random)){
                var.w <- max(p.LT[n.par]^2,.Machine$double.xmin)
                logHazBeta <- env$parent.logHazBeta[status.one]
                res.temp.LT <- Frailty.Adapt(nodes=x.H, nodessquare=x.H.2, logweights=log.rho.H, clust=n.by.clust, clustd=n.by.clust.delta, expect=lambda.pop.delta, betal=logHazBeta, betaL=parent.logCumBeta, A=parent.cst.adj, var=var.w, muhatcond=mu.hat.LT)
                if (mu.hat.LT==1)
                    res.LT <- res.temp.LT[1:n.clust]
                else if (mu.hat.LT==2)
                    res.LT <- res.temp.LT[n.clust+(1:n.clust)]
                else {
                    env$parent.cst.adj <- res.temp.LT[2*n.clust+(1:n.clust)]
                    res.LT <- res.temp.LT[3*n.clust+1]
                }
            }
            else {
                if (BoolExp==1){
                    temp.l <- exp(parent.logHazBeta)+lambda.pop
                    if (sum(temp.l-lambda.pop)==0){
                        res.LT <- .Machine$double.xmax
                    } # Prevents a kind of catastrophic cancellation
                    else {
                        log.lambda <- log(temp.l)
                        log.lambda[log.lambda==Inf] <- .Machine$double.xmax
                        res.LT <- -sum(-exp(parent.logCumBeta) + status.obs*log.lambda)
                        res.LT <- min(res.LT,.Machine$double.xmax)
                    }
                }
                else {
                    log.lambda <- parent.logHazBeta
                    log.lambda[log.lambda==Inf] <- .Machine$double.xmax
                    res.LT <- -sum(-exp(parent.logCumBeta) + status.obs*log.lambda)
                    res.LT <- min(res.LT,.Machine$double.xmax)
                }
            }
        }
        if (mu.hat.LT==0){
            res.LT[is.nan(res.LT) | abs(res.LT)==Inf] <- .Machine$double.xmax
            if (sum(p.LT-parent.param)){
                verbose.ll()
                env$parent.neval <- parent.neval + 1
            }
            env$parent.param <- p.LT+0
            env$parent.ptd <- p.td
            env$parent.pntd <- p.ntd
            env$parent.logLik <- res.LT
        }
        return(res.LT)

    }

    # Launching the optimisation procedure
    if (fnoptim=="nlm"){
        time0 <- as.numeric(proc.time()[3])
        mod.lik <- nlm(LL.Tot,init,hessian=TRUE,iterlim=iterlim,print.level=print.level,...)
        param.fin <- mod.lik$estimate
        code.fin <- mod.lik$code
        loglik <- -mod.lik$minimum
        iterations <- mod.lik$iterations
    }
    else if (fnoptim=="optim"){
        time0 <- as.numeric(proc.time()[3])
        mod.lik <- optim(init,LL.Tot,hessian=TRUE,method=method,...)
        param.fin <- mod.lik$par
        code.fin <- mod.lik$convergence
        loglik <- -mod.lik$value
        iterations <- "---"
    }
    names(param.fin) <- param.names

    # Hessian matrix
    hessian.fin <- mod.lik$hessian
    vcov <- try(solve(hessian.fin),silent=TRUE)
    if (is.character(vcov)){
        warning("Unable to invert the Hessian matrix...")
    }
    colnames(vcov) <- rownames(vcov) <- param.names

    # Empirical Bayes Estimates
    if (!is.null(random)) {
        mu.hat.fin <- LL.Tot(param.fin,mu.hat.LT=1)
        mu.hat.df <- data.frame(Cluster=clust,Mu.Hat=mu.hat.fin)
        param.fin[n.par] <- abs(param.fin[n.par])
    }
    else {
        mu.hat.df <- 0
    }

    time1 <- as.numeric(proc.time()[3])
    PrintData <- data.frame(Name=name.data,N.Obs.Tot=n.obs.tot,N.Obs=n.obs,N.Events=n.events,N.Clust=n.clust,row.names="")
    PrintDetails <- data.frame(Iter=iterations,Eval=parent.neval+1,Base=base,Nb.Leg=n.gleg,
               Nb.Aghq=n.aghq,Optim=fnoptim,Method=ifelse(fnoptim=="optim",method,"---"),
               Code=code.fin,LogLik=loglik,Total.Time=(time1-time0),row.names="")

    # Part of the results printed on screen
    cat("\nData\n")
    print(PrintData)
    cat("\nDetails\n")
    print(PrintDetails)

    res.FAR <- list(dataset=name.data,
         call=call,
         formula=formula,
         xlevels=Xlevels,
         n.obs.tot=n.obs.tot,
         n.obs=n.obs,
         n.events=n.events,
         n.clust=ifelse(!is.null(random),n.clust,1),
         n.time.0=length(Idx.Time.0),
         base=base,
         max.time=max.time,
         bounds=c(0,Bo),
         degree=ifelse(base=="exp.bs",degree,NA),
         knots=knots,
         names.ph=names.fix,
         random=ifelse(!is.null(random),random,NA),
         coefficients=param.fin,
         std.errors=sqrt(diag(vcov)),
         vcov=vcov,
         mu.hat=mu.hat.df,
         n.par=n.par,
         n.gleg=n.gleg,
         n.aghq=n.aghq,
         fnoptim=fnoptim,
         method=ifelse(fnoptim=="optim",method,NA),
         code=code.fin,
         loglik=loglik,
         iter=iterations,
         eval=parent.neval+1,
         time.elapsed=time1-time0)

    class(res.FAR) <- "mexhaz"
    res.FAR

}
