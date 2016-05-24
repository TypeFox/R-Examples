geeUOmega<-function(geeOutput){
    fit<-geeOutput
    # I think this gives the same as m<-match.call() from within gee 
    m<-fit$call
    # next few lines from various places in gee function
    m$R <- m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <- m$corstr <- m$Mv <- m$silent <- m$contrasts <- m$family <- m$scale.fix <- m$scale.value <- m$v4.4compat <- NULL
    if (is.null(m$id)) 
        m$id <- as.name("id")
    if (!is.null(m$na.action) && m$na.action != "na.omit") {
        warning("Only 'na.omit' is implemented for gee\ncontinuing with 'na.action=na.omit'")
        m$na.action <- as.name("na.omit")
    }
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    y <- as.matrix(model.extract(m, "response"))
    x <- model.matrix(Terms, m, contrasts)
    Q <- qr(x)
    if (Q$rank < ncol(x)) 
        stop("rank-deficient model matrix")
    N <- rep(1, length(y))
    if (dim(y)[2] == 2) {
        N <- as.vector(y %*% c(1, 1))
        y <- y[, 1]
    } else {
        if (dim(y)[2] > 2) 
            stop("Only binomial response matrices (2 columns)")
    }
    offset <- model.extract(m, offset)
    if (length(offset) <= 1) 
        offset <- rep(0, length(y))
    if (length(offset) != length(y)) 
        stop("offset and y not same length")

    x <- model.matrix(Terms, m, fit$contrasts)
    eta <- as.vector(x %*% fit$coefficients)
    ### end of lines from gee
    # correction to gee: forgot to add back offset
    # residuals and fitted.values were wrong if offset was present
    # calculate dmu/deta without offset added 
    DmuDeta<-fit$family$mu.eta(eta)
    # calculate mu with offset added
    mu <- as.vector(fit$family$linkinv(eta+offset))
    fit$fitted.values <- mu
    # correction to gee: residuals are wrong if original dim(y)==2
    # i.e., wrong if N!=rep(1,length(y))
    # change residuals from 'y-mu' to 'y/N - mu'
    fit$residuals <- y/N - mu

     calc.U.Omega<-function(
        S=fit$residuals,
        ID=fit$id,
        dmu.deta=DmuDeta,
        A=fit$family$variance(mu),
        R=fit$working.correlation,
        X=x,p=length(fit$coefficients),
        phi=fit$scale){
        
        
        uid<-unique(ID)
        K<-length(uid)
        U<-matrix(NA,K,p)
        Omega<-array(NA,dim=c(K,p,p))
        sqrtAoverN<-sqrt(A/N)
        for (i in 1:K){
            pick<- ID==uid[i]
            dmu.detai<-dmu.deta[pick]
            ni<-length(dmu.detai)
            if (ni==1){
                Di<- dmu.detai*t(as.matrix(X[pick,]))
                Vi<- phi * sqrtAoverN[pick]*sqrtAoverN[pick]
                ViInv<- 1/Vi
                U[i,]<-  t(Di) * ViInv * S[pick]
                Omega[i,,]<-  ViInv* t(Di) %*% Di
            } else {
                Di<- diag(dmu.detai) %*% X[pick,]
                Vi<- phi * diag(sqrtAoverN[pick]) %*% R[1:ni,1:ni] %*%  diag(sqrtAoverN[pick])
                ViInv<-solve(Vi)
                U[i,]<- t(Di) %*% ViInv %*% S[pick]
                Omega[i,,]<-  t(Di) %*% ViInv %*% Di
            }
       }
        list(U=U,Omega=Omega)
    }
    out<-calc.U.Omega()
    fit$u<-out$U
    fit$omega<-out$Omega
    fit
}