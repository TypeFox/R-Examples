#se n=1000: value out of range in 'gammafn'
#warning se "lm" con "glm"???

`davies.test` <-
function (obj, seg.Z, k = 10, alternative = c("two.sided", "less", "greater"),
    type=c("lrt","wald"), values=NULL, dispersion=NULL) {
#    extract.t.value.U<-function(x){
#        #estrae il t-value dell'ultimo coeff in un oggetto restituito da lm.fit
#        #non serve... in realta' viene usata extract.t.value.U.glm()
#            #x<-x$obj
#            R<-qr.R(x$qr)
#            p<-ncol(R)
#            n<-length(x$fitted.values)
#            invR<-backsolve(R,diag(p))
#            hat.sigma2<-sum(x$residuals^2)/(n-p)
#            #solve(crossprod(qr.X(x$qr)))
#            V<-tcrossprod(invR)*hat.sigma2
#            tt<-x$coefficients[p]/sqrt(V[p,p])
#            tt}
#-------------------------------------------------------------------------------
daviesLM<-function(y, z, xreg, weights, offs, values, k, alternative){
#Davies test with sigma unknown
#--------------
#> gammaA<-function(x){
#  x^(x-.5)*exp(-x)*sqrt(2*pi)*(1+1/(12*x)+1/(288*x^2)-139/(51840*x^3) -571/(2488320*x^4))
#   }
#exp(lgamma())
       fn="pmax(x-p,0)"
       y<-y-offs
       n<-length(y)
       n1<-length(values)
       RIS<-matrix(,n1,2)

       X.psi<-matrix(,n,length(fn))
       df.res<- n - ncol(xreg) - length(fn)
       for(i in 1:n1){
          for(j in 1:length(fn)) X.psi[,j]<-eval(parse(text=fn[[j]]), list(x=z, p=values[i]))
          xx1.new<-cbind(X.psi,xreg)
          #lrt
          #mu1.new<-xx1.new%*%solve(crossprod(xx1.new), crossprod(xx1.new,y))
          #rss1<-sum((y-mu1.new)^2)
          #sigma2<-if(missing(sigma)) rss1/(n-ncol(xx1.new)) else sigma^2
          #RIS[i]<-((rss0-rss1)/ncol(X.psi))/sigma2
          #Wald
          invXtX1<-try(solve(crossprod(sqrt(weights)*xx1.new)), silent=TRUE)
          if(class(invXtX1)!="try-error"){
            hat.b<-drop(invXtX1%*%crossprod(weights*xx1.new,y))
            mu1.new<-xx1.new%*%hat.b
            devE<-sum((weights*(y-mu1.new)^2))
            hat.sigma<- sqrt(devE/df.res)
            RIS[i,1]<-hat.b[1]/(hat.sigma*sqrt(invXtX1[1, 1]))
            Z<-hat.b[1]/(sqrt(invXtX1[1, 1]))
            D2<- Z^2 + devE
            RIS[i,2]<-Z^2/D2 #beta
            }
       }
       valori<-values[!is.na(RIS[,1])]
       RIS<- RIS[!is.na(RIS[,1]),]
       V<-sum(abs(diff(asin(RIS[,2]^.5))))

       onesided <- TRUE
       if (alternative == "less") {
        M <- min(RIS[,1])
        best<-valori[which.min(RIS[,1])]
        p.naiv <- pt(M, df=df.res, lower.tail = TRUE)
          } else if (alternative == "greater") {
        M <- max(RIS[,1])
        best<-valori[which.max(RIS[,1])]
        p.naiv <- pt(M, df=df.res, lower.tail = FALSE)
          } else {
        M <- max(abs(RIS[,1]))
        best<-valori[which.max(abs(RIS[,1]))]
        p.naiv <- pt(M, df=df.res, lower.tail = FALSE)
        onesided <- FALSE
        }
        u<-M^2/((n-ncol(xx1.new))+ M^2)
        approxx<-V*(((1-u)^((df.res-1)/2))*gamma(df.res/2+.5))/(2*gamma(df.res/2)*pi^.5)
        p.adj <- p.naiv + approxx
        p.adj <- ifelse(onesided, 1, 2) * p.adj
        p.adj<-list(p.adj=p.adj, valori=valori, ris.valori=RIS[,1], best=best)
        return(p.adj)
#       M<-max(abs(RIS[,1]))
#       u<-M^2/((n-ncol(xx1.new))+ M^2)
#       approxx<-V*(((1-u)^((df.res-1)/2))*gamma(df.res/2+.5))/(2*gamma(df.res/2)*pi^.5)
#       p.naiv<-pt(-abs(M), df=df.res) #naive p-value
#       p.adj<-2*(p.naiv+approxx) #adjusted p-value (upper bound)
#       p.adj<-min(p.adj, 1)
#       p.adj<-list(p.adj=p.adj, valori=values, ris.valori=RIS[,1], approxx=approxx, p.naiv=p.naiv)
#       return(p.adj)
       }
#--------------------------------
daviesGLM<-function(y, z, xreg, weights, offs, values=NULL, k, list.glm, alternative){
#Davies test for GLM (via LRT or Wald)
    est.dispGLM<-function(object){
        df.r <- object$df.residual
        dispersion <- if(object$family$family%in%c("poisson","binomial")) 1 else object$dev/df.r
        dispersion
        }        
    extract.t.value.U.glm<-function(object,dispersion,isGLM=TRUE){
        #estrae il t-value dell'ultimo coeff in un oggetto restituito da lm.wfit/glm.fit
        est.disp <- FALSE
        df.r <- object$df.residual
        if (is.null(dispersion))
          dispersion <- if(isGLM&&(object$family$family%in%c("poisson","binomial"))) 1
        else if (df.r > 0) {
            est.disp <- TRUE
            if (any(object$weights == 0))
                warning("observations with zero weight not used for calculating dispersion")
            sum((object$weights * object$residuals^2)[object$weights > 0])/df.r
        }  else {
            est.disp <- TRUE
            NaN
        }
        dispersion<-max(c(dispersion, 1e-10))
        p <- object$rank
        p1 <- 1L:p
        Qr <- object$qr
        coef.p <- object$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        tvalue <- coef.p[1]/sqrt(covmat[1,1]) #<0.4.0-0 era coef.p[p]/sqrt(covmat[p,p])
        tvalue
        }#end extract.t.value.U.glm
#--------------
       fn<-"pmax(x-p,0)"
       dev0<-list.glm$dev0
       eta0<-list.glm$eta0
       family=list.glm$family
       type<-list.glm$type
       dispersion<-list.glm$dispersion
       n<-length(y)
       r<-length(fn)
       n1<-length(values)
       RIS<-rep(NA, n1)
       X.psi<-matrix(,n,length(fn))
       for(i in 1:n1){
          for(j in 1:length(fn)) X.psi[,j]<-eval(parse(text=fn[[j]]), list(x=z, p=values[i]))
          xreg1<-cbind(X.psi,xreg)
          o1<-glm.fit(x = xreg1, y = y, weights = weights, offset = offs,
                  family=family, etastart=eta0)
          dev<-o1$dev
          if (is.list(o1) && ncol(xreg1)==o1$rank) {
            RIS[i]<- if(type=="lrt") sqrt((dev0-dev)/est.dispGLM(o1))*sign(o1$coef[1]) else extract.t.value.U.glm(o1,dispersion)
            }
          }
      valori<-values[!is.na(RIS)]
      ris.valori<-RIS[!is.na(RIS)]
      V<-sum(abs(diff(ris.valori)))
      #-----Questo e' se il test di riferimento e' una \chi^2_r. (Dovresti considerare il LRT non segnato)
      #V<-sum(abs(diff(sqrt(RIS))))#nota sqrt
      #M<- max(RIS)
      #approxx<-(V*(M^((r-1)/2))*exp(-M/2)*2^(-r/2))/gamma(r/2)
      #p.naiv<-1-pchisq(M,df=r) #naive p-value
      #p.adj<-min(p.naiv+approxx,1) #adjusted p-value (upper bound)
    onesided <- TRUE
    if (alternative == "less") {
        M <- min(ris.valori)
        best<-valori[which.min(ris.valori)]
        p.naiv <- pnorm(M, lower.tail = TRUE)
    }
    else if (alternative == "greater") {
        M <- max(ris.valori)
        best<-valori[which.max(ris.valori)]
        p.naiv <- pnorm(M, lower.tail = FALSE)
    }
    else {
        M <- max(abs(ris.valori))
        best<-valori[which.max(abs(ris.valori))]
        p.naiv <- pnorm(M, lower.tail = FALSE)
        onesided <- FALSE
    }
    approxx<-V*exp(-(M^2)/2)/sqrt(8*pi)
    p.adj <- p.naiv + approxx
    p.adj <- ifelse(onesided, 1, 2) * p.adj
    p.adj<-list(p.adj=p.adj, valori=valori, ris.valori=ris.valori, best=best)
    return(p.adj)
    }
#-------------------------------------------------------------------------------
    if(!inherits(obj, "lm")) stop("A 'lm', 'glm', or 'segmented' model is requested")
    if(class(seg.Z)!="formula") stop("'seg.Z' should be an one-sided formula")
    if(k<=1) stop("k>1 requested! k>=10 is recommended")
    if(k<10) warnings("k>=10 is recommended")
    alternative <- match.arg(alternative)
    type        <- match.arg(type)
    if(length(all.vars(seg.Z))>1) warning("multiple segmented variables ignored in 'seg.Z'",call.=FALSE)
    isGLM<-"glm"%in%class(obj)
    Call<-mf<-obj$call
    mf$formula<-formula(obj)
    m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
    formulaOrig<-formula(obj)
    if(class(obj)[1]=="segmented"){
        if(!is.null(eval(obj$call$obj)$call$data)) mf$data <- eval(obj$call$obj)$call$data
        mf$formula<-update.formula(mf$formula,paste("~.-",paste(obj$nameUV$V, collapse="-")))
        for(i in 1:length(obj$nameUV$U)) assign(obj$nameUV$U[i], obj$model[,obj$nameUV$U[i]], envir=parent.frame())
        formulaOrig<-update.formula(formulaOrig, paste("~.-",paste(obj$nameUV$V, collapse="-")))
        }
    mf <- eval(mf, parent.frame())
    weights <- as.vector(model.weights(mf))
    offs <- as.vector(model.offset(mf))
    if(!is.null(Call$weights)){ #"(weights)"%in%names(mf)
      names(mf)[which(names(mf)=="(weights)")]<-all.vars(Call$weights) #as.character(Call$weights)
      #aggiungere???
      # mf["(weights)"]<-weights
      }
    mt <- attr(mf, "terms")
    interc<-attr(mt,"intercept")
    y <- model.response(mf, "any")
    XREG <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    n <- nrow(XREG)
    if (is.null(weights)) weights <- rep(1, n)
    if (is.null(offs)) offs <- rep(0, n)
    name.Z <- all.vars(seg.Z)
    Z<-XREG[,match(name.Z, colnames(XREG))]
    if(!name.Z %in% names(coef(obj))) XREG<-XREG[,-match(name.Z, colnames(XREG)),drop=FALSE]
    list.glm<-list(dev0=obj$dev, eta0=obj$linear.predictor, family=family(obj),
      type=type, dispersion=dispersion)
    
    if(is.null(values)) values<-seq(sort(Z)[2], sort(Z)[(n - 1)], length = k)
       #values<-seq(min(z), max(z), length=k+2)
       #values<-values[-c(1,length(values))]
    if(class(obj)=="lm" || identical(class(obj),c("segmented","lm")) ) {
        if(n<=300) { 
          rr<-daviesLM(y=y, z=Z, xreg=XREG, weights=weights, 
        offs=offs, values=values, k=k, alternative=alternative)
           } else { 
           list.glm$family<-gaussian()
           list.glm$type<-"wald"
           rr<-daviesGLM(y=y, z=Z, xreg=XREG, weights=weights, 
              offs=offs, values=values, k=k, list.glm=list.glm, alternative=alternative)
          }
        }
    if(identical(class(obj),c("glm","lm")) || identical(class(obj),c("segmented","glm","lm"))) rr<-daviesGLM(y=y, z=Z, xreg=XREG, weights=weights, 
        offs=offs, values=values, k=k, list.glm=list.glm, alternative=alternative)
      best<-rr$best
      p.adj<-rr$p.adj
      valori<-rr$valori
      ris.valori<-rr$ris.valori
    if(is.null(obj$family$family)) {
          famiglia<-"gaussian"
          legame<-"identity"} else {
               famiglia<-obj$family$family
               legame<-obj$family$link
          }
    out <- list(method = "Davies' test for a change in the slope",
#        data.name=paste("Model = ",famiglia,", link =", legame,
#        "\nformula =", as.expression(formulaOrig),
#        "\nsegmented variable =", name.Z),
        data.name=paste("formula =", as.expression(formulaOrig), ",   method =", obj$call[[1]] ,
        "\nmodel =",famiglia,", link =", legame, if(isGLM) paste(", statist =", type) else NULL ,
        "\nsegmented variable =", name.Z),
        statistic = c("'best' at" = best),
        parameter = c(n.points = length(valori)), p.value = min(p.adj,1),
        alternative = alternative, process=cbind(psi.values=valori, stat.values=ris.valori))
    class(out) <- "htest"
    return(out)
}


