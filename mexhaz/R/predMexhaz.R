predMexhaz <- function(model,time.pts,data.val=data.frame(.NotUsed=NA),conf.int=c("none","delta","simul"),nb.sim=10000){

    time0 <- as.numeric(proc.time()[3])
    conf.int <- match.arg(conf.int)
    call <- model$call
    formula <- model$formula
    Xlev <- model$xlevels
    base <- model$base
    degree <- model$degree
    knots <- model$knots
    n.gleg <- model$n.gleg
    Bo <- model$bounds[2]
    max.time <- model$max.time
    coef <- model$coefficients
    vcov <- model$vcov
    df <- model$n.obs-model$n.par

    Idx.T.NA <- which(is.na(time.pts) | time.pts<0)
    if (length(Idx.T.NA)>0){
        stop("The 'time.pts' argument contains NA or negative values...")
    }
    Idx.T.Max <- which(time.pts>max.time)
    if (length(Idx.T.Max)>0){
        warning(paste("The model cannot be used to predict survival for times greater than ",max.time," (maximum follow-up time on which the model estimation was based). Consequently, these time values were removed from the 'time.pts' vector.",sep=""))
        time.pts <- time.pts[-Idx.T.Max]
    }
    time.pts <- time.pts[time.pts!=0]
    if (length(time.pts)==0){
        stop("The 'time.pts' argument contains no values for which predictions can be made...")
    }
    time.pts <- time.pts[order(time.pts)]
    nb.time.pts <- length(time.pts)
    if (nb.time.pts>1 & dim(data.val)[1]>1) {
        stop("Predictions can be made for n individuals at 1 time point or for 1 individual at m time points but not for several individuals at several time points...")
    }

    if (nb.time.pts>1){
        typepred <- "multitime"
    }
    else {
        typepred <- "multiobs"
    }

    if (conf.int=="delta"){
        if (base=="weibull"){
            warning("Computation of confidence intervals by the Delta Method is not available for base='weibull'. Consider using simulations instead.")
            delta <- FALSE
        }
        else {
            delta <- TRUE
        }
    }
    else {
        delta <- FALSE
    }

    if (conf.int=="simul" & (nb.sim<=0 | (round(nb.sim,0)!=nb.sim))){
        warning("The 'nb.sim' argument must be a strictly positive integer...")
    }

    tot.formula <- terms(formula,data=data.val,specials="nph")
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
    data.obs.temp <- cbind(time.pts, data.val, row.names = NULL)
    FormulaF <- update(tot.formula,paste("NULL~.",nTerm2,sep="-"))
    Test <- model.frame(FormulaF,data=data.obs.temp,xlev=Xlev,na.action=NULL) # Will raise an error if one of the variables necessary to the model is missing from data.val
    if (dim(Test)[1]==0){
        stop("The formula is incompatible with the data provided. If bs() was used to model the effect of some of the variables included in the formula, make sure that the spline boundaries were specified...")
    }

    Idx.Dv.NA <- which(apply(Test,1,function(x){sum(is.na(x))})>0)
    if (length(Idx.Dv.NA)>0){
        warning("Some rows in the 'data.val' data.frame had NA values. Consequently, these rows were deleted.")
        Test <- Test[-Idx.Dv.NA,,drop=FALSE]
    }
    if (dim(Test)[1]==0){
        stop("The 'data.val' dataframe contains no valid row on which predictions can be based.")
    }

    data.obs <- cbind(time.pts,Test,row.names=NULL)
    time.pts <- data.obs[,1]
    nb.time.pts <- length(time.pts)

    FALCenv <- environment()

    IntBs1 <- function(x,nph,param,leint,whint,knots){
        .Call("IntBs1",x,nph,param,leint,whint,knots)
    }
    IntBs23 <- function(x,nph,timecat,param,deg,n,lw,matk,totk){
        .Call("IntBs23",x,nph,timecat,param,deg,n,lw,as.double(matk),as.double(totk))
    }
    IntPwCst <- function(nph,param,leint,lerem,whint){
        .Call("IntPwCst",nph,param,leint,lerem,whint)
    }
    DeltaBs1 <- function(x,nph,param,others,leint,whint,knots,varcov){
        .Call("DeltaBs1",x,nph,param,others,leint,whint,knots,varcov)
    }
    DeltaBs23 <- function(x,nph,timecat,param,others,varcov,k,n,lw,matk,totk){
        .Call("DeltaBs23",x,nph,timecat,param,others,varcov,k,n,lw,matk,totk)
    }
    DeltaPwCst <- function(nph,param,others,leint,lerem,whint,varcov){
        .Call("DeltaPwCst",nph,param,others,leint,lerem,whint,varcov)
    }

    Surv <- NULL
    lambda <- NULL
    BInf1D <- NULL
    BInf2D <- NULL
    BSup1D <- NULL
    BSup2D <- NULL
    Var.log.Haz <- NULL
    Var.log.Cum <- NULL

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

        # Creating the points and weights of Gauss-Legendre quadrature (for computing the cumulative hazard when base="exp.bs")
        gl <- gauss.quad(n=n.gleg,kind="legendre")
        gln <- gl$nodes
        lglw <- log(gl$weights)
    }

    # Formatting the data
    fix.obs <- model.matrix(FormulaF,data=Test)
    names.fix <- colnames(fix.obs)[-1]
    nph.obs <- model.matrix(FormulaN,data=Test)
    nbtd <- dim(nph.obs)[2]
    names.nph <- colnames(nph.obs)[-1]

    # Creation of the different objects necessary to compute the hazards

    # Weibull hazard
    if (base=="weibull"){
        log.time.pts <- log(time.pts)
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
            which.td <- c(2,(2+(n.ntd-1))+1:n.td.nph)
        }
        else {
            which.td <- 2
        }
        fix.obs <- fix.obs[,-1,drop=FALSE]
        nph.obs <- nph.obs[,-1,drop=FALSE]
    }

    # Hazard modelled by the exponential of a B-spline
    else if (base=="exp.bs"){

        # Baseline hazard-related objects
        cuts <- c(0,knots,Bo)
        interv <- c(knots,Bo)-c(0,knots)
        time.cat <- cut(time.pts,breaks=cuts)
        time.cat <- as.numeric(time.cat)-1
        vec.knots <- c(rep(0,degree),knots,rep(Bo,degree))
        if (degree==3){
            MatK <- Transf3(vec.knots)
        }
        else if (degree==2){
            MatK <- Transf2(vec.knots)
        }
        n.td.base <- degree + length(knots)

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
        nph.obs <- t(nph.obs) # Matrix of time-dependent effects has to be transposed for use by the IntBs'xx'Stat functions
    }

    # Hazard modelled by the exponential of a piecewise constant function
    else if (base=="pw.cst") {

        # Baseline hazard-related objects
        cuts <- c(0,knots,Bo)
        interv <- c(knots,Bo)-c(0,knots)
        time.cat <- cut(time.pts,breaks=cuts)
        names.base <- levels(time.cat)
        time.cat <- as.numeric(time.cat)-1
        time.remain <- time.pts-cuts[time.cat+1]
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
        nph.obs <- t(nph.obs) # Matrix of time-dependent effects has to be transposed for use by the IntPwCstStat function
    }

    # Functions that compute the hazards (the computed individual-specific hazards are stored in the main function environment)

    if (base=="weibull"){
        HazardP <- function(p.td.H,p.ntd.H,Delta=delta,env=FALCenv){
            log.p.LT.1 <- log(p.ntd.H[1]) + as.vector(fix.obs%*%p.ntd.H[-1])
            log.p.LT.2 <- log(p.td.H[1]) + as.vector(nph.obs%*%p.td.H[-1])
            l.lambda.beta <- log.p.LT.2 +log.time.pts*(exp(log.p.LT.2)-1)+log.p.LT.1
            l.Lambda.beta <- log.time.pts*exp(log.p.LT.2)+log.p.LT.1
            env$lambda <- exp(l.lambda.beta)
            env$Surv <- exp(-exp(l.Lambda.beta))
        }
    }

    else if (base=="exp.bs"){
        HazardP <- function(p.td.H,p.ntd.H,Delta=delta,env=FALCenv){
            if (degree==1) {
                temp.H <- IntBs1(time.pts,nph=nph.obs,param=p.td.H,leint=interv,whint=time.cat,knots=vec.knots)
            }
            else {
                temp.H <- IntBs23(time.pts,nph=nph.obs,timecat=time.cat,param=p.td.H,deg=degree,gln,lglw,matk=MatK,totk=vec.knots)
            }
            l.lambda <- temp.H[1:nb.time.pts]
            l.Lambda <- temp.H[(nb.time.pts)+1:nb.time.pts]
            if (n.ntd!=0){
                beta.x <- as.vector(fix.obs%*%p.ntd.H)
            }
            else {
                beta.x <- 0
            }
            env$lambda <- exp(l.lambda + beta.x)
            env$Surv <- exp(-exp(l.Lambda + beta.x))
            if (Delta){
                if (degree==1){
                    Var <- DeltaBs1(time.pts, nph=nph.obs, param=c(p.ntd.H,p.td.H), others=t(fix.obs), leint=interv, whint=time.cat, knots=vec.knots, varcov=vcov[c(which.ntd,which.td),c(which.ntd,which.td)])
                }
                else {
                    Var <- DeltaBs23(time.pts, nph=nph.obs, timecat=time.cat, param=c(p.ntd.H,p.td.H), others=t(fix.obs), varcov=vcov[c(which.ntd,which.td),c(which.ntd,which.td)], k=degree, n=gln, lw=lglw, matk=MatK, totk=vec.knots)
                }
                env$Var.log.Haz <- Var[1:nb.time.pts]
                env$Var.log.Cum <- Var[(nb.time.pts) + 1:nb.time.pts]
                env$BInf1D <- exp(l.lambda + beta.x + qnorm(0.025)*sqrt(env$Var.log.Haz))
                env$BSup1D <- exp(l.lambda + beta.x + qnorm(0.975)*sqrt(env$Var.log.Haz))
                env$BSup2D <- exp(-exp(l.Lambda + beta.x + qnorm(0.025)*sqrt(env$Var.log.Cum)))
                env$BInf2D <- exp(-exp(l.Lambda + beta.x + qnorm(0.975)*sqrt(env$Var.log.Cum)))
            }
        }
    }

    else if (base=="pw.cst"){
        HazardP <- function(p.td.H,p.ntd.H,Delta=delta,env=FALCenv){
            temp.H <- IntPwCst(nph=nph.obs,param=p.td.H,leint=interv,lerem=time.remain,whint=time.cat)
            l.lambda <- temp.H[1:nb.time.pts]
            l.Lambda <- temp.H[(nb.time.pts)+1:nb.time.pts]
            if (n.ntd!=0){
                beta.x <- as.vector(fix.obs%*%p.ntd.H)
            }
            else {
                beta.x <- 0
            }
            env$lambda <- exp(l.lambda + beta.x)
            env$Surv <- exp(-exp(l.Lambda + beta.x))
            if (Delta){
                Var <- DeltaPwCst(nph=nph.obs, param=c(p.ntd.H,p.td.H), others=t(fix.obs), leint=interv, lerem=time.remain, whint=time.cat, varcov=vcov[c(which.ntd,which.td),c(which.ntd,which.td)])
                env$Var.log.Haz <- Var[1:nb.time.pts]
                env$Var.log.Cum <- Var[(nb.time.pts) + 1:nb.time.pts]
                env$BInf1D <- exp(l.lambda + beta.x + qnorm(0.025)*sqrt(env$Var.log.Haz))
                env$BSup1D <- exp(l.lambda + beta.x + qnorm(0.975)*sqrt(env$Var.log.Haz))
                env$BSup2D <- exp(-exp(l.Lambda + beta.x + qnorm(0.025)*sqrt(env$Var.log.Cum)))
                env$BInf2D <- exp(-exp(l.Lambda + beta.x + qnorm(0.975)*sqrt(env$Var.log.Cum)))
            }
        }
    }

    if (base!="weibull" & delta){
        HazardP(coef[which.td],coef[which.ntd])
        res1 <- data.frame(hazard=lambda,hazard.inf=BInf1D,hazard.sup=BSup1D,
                           surv=Surv,surv.inf=BInf2D,surv.sup=BSup2D)
    }
    else if (conf.int=="simul" & nb.sim>1 & round(nb.sim,0)==nb.sim){
        Res1 <- matrix(0,nb.time.pts,nb.sim)
        Res2 <- matrix(0,nb.time.pts,nb.sim)
        Coef <- mvrnorm(nb.sim,mu=coef,Sigma=vcov)
        for (i in 1:nb.sim){
            p.td <- Coef[i,which.td]
            p.ntd <- Coef[i,which.ntd]
            HazardP(p.td,p.ntd)
            Res1[,i] <- lambda
            Res2[,i] <- Surv
        }
        BInf1 <- apply(Res1,1,quantile,prob=0.025)
        BSup1 <- apply(Res1,1,quantile,prob=0.975)
        BInf2 <- apply(Res2,1,quantile,prob=0.025)
        BSup2 <- apply(Res2,1,quantile,prob=0.975)
        HazardP(coef[which.td],coef[which.ntd])
        res1 <- data.frame(hazard=lambda,hazard.inf=BInf1,hazard.sup=BSup1,
                           surv=Surv,surv.inf=BInf2,surv.sup=BSup2)
    }
    else {
        HazardP(coef[which.td],coef[which.ntd])
        res1 <- data.frame(hazard=lambda,hazard.inf=NA,hazard.sup=NA,
                       surv=Surv,surv.inf=NA,surv.sup=NA)
    }
    res.PS <- list(call=call,
                   results=cbind(time.pts,Test,res1),
                   variances=data.frame(var.log.haz=Var.log.Haz,var.log.cum=Var.log.Cum),
                   type=typepred,ci.method=conf.int,
                   nb.sim=ifelse(conf.int=="simul",nb.sim,NA))

    class(res.PS) <- "predMexhaz"
    res.PS
}
