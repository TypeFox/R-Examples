gcrq <-
function(formula, tau=c(.1,.25,.5,.75,.9), data, subset, weights, na.action, transf=NULL,
    y=TRUE, interc=FALSE, foldid=NULL, nfolds=10, cv=FALSE, n.boot=0, eps=0.0001, ...){
#Growth Charts via QR
#**************weights??
#se cv=TRUE restituisce anche una componente 'cv' che e' una matrice di n.righe=n.valori di lambda e colonne nfolds
#eps in control??
#foldid, nfold usati se lambda in ps() e' un vettore
#... 
#... passati all'interno di ps()? vedi
bspline <- function(x, ndx, xlr = NULL, knots=NULL, deg = 3, deriv = 0, outer.ok=FALSE) {
    # x: vettore di dati
    # xlr: il vettore di c(xl,xr)
    # ndx: n.intervalli in cui dividere il range
    # deg: il grado della spline
    # Restituisce ndx+deg basis functions per ndx-1 inner nodi
    #ci sono "ndx+1" nodi interni + "2*deg" nodi esterni
#    require(splines)
  if(is.null(knots)) {
    if (is.null(xlr)) {
        xl <- min(x) - 0.01 * diff(range(x))
        xr <- max(x) + 0.01 * diff(range(x))
    }
    else {
        if (length(xlr) != 2)
            stop("quando fornito, xlr deve avere due componenti")
        xl <- xlr[1]
        xr <- xlr[2]
    }
    dx <- (xr - xl)/ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
      }
      #else {
      #if(length(knots)!=(ndx+1+2*deg)) stop("errore nel numero di nodi fornito")
      #}
    B <- splineDesign(knots, x, ord = deg + 1, derivs = rep(deriv, length(x)), outer.ok=outer.ok)
    B
    }
    call <- match.call()
    if (missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    intercMt<-attr(mt,"intercept")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) names(Y) <- nm
        }
    if(!is.null(transf)) Y <- eval(parse(text=transf), list(y=Y)) 
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else stop("error in the design matrix")#matrix(, NROW(Y), 0L)
    attrContr<-attr(X, "contrasts")
    weights <- as.vector(model.weights(mf))
    n<-nrow(X)
    if(!is.null(weights) && !is.numeric(weights)) stop("'weights' must be a numeric vector")
    if(!is.null(weights) && any(weights < 0)) stop("negative weights not allowed")
    tf <- terms(formula, specials = c("ps","ridge"))
    id.ps<-attr(tf,"specials")$ps #puo' essere un vettore se ci sono piu' termini..non dipende se c'e' intercetta o meno..
    id.ridge<-attr(tf,"specials")$ridge
    #NB: (id.ps e id.ridge) posizione della variabile ps(x) nel modelframe (include y ma non dipende da l'intercetta)
    #---- se vuoi controllare il numero di termini ps
    #Il problema e' che devi modificare le funzioni fitter ncross.rq.fitB() et al.. in cui devi modificare le matrici R e r per consentire 
    #monot, noncrossing ecc.. (in realta' noncrossing credo che sia gia' garantito perche' e' soltanto sui coeff.)
    #if(length(id.ps)>1) stop("A single smooth term is allowed")
    #------------------------------------------------------
    nomiCoefUNPEN<-names(mf)[-c(1,id.ps,id.ridge)]
    testo.ps<-names(mf)[id.ps] #colnames(X)[id.ps-1+interc] #questa e' la stringa del tipo csda(..)
    testo.ridge<-names(mf)[id.ridge]
    #se NON ci sono termini di penalizz.
    if(length(testo.ps)<=0 && length(testo.ridge)<=0){
#        fitter<- if(length(tau)>1) get("arq.fitMulti") else get("arq.fit")
#        fit<-fitter(y=Y,X=X,tau=tau,g=g,beta0=b.start,control=control)
      fit<-ncross.rq.fitX(y=Y, X=X, taus=tau, eps=eps) #X gia' include l'intercetta
      if(n.boot>0){
          coef.boot<-array(, dim=c(nrow(as.matrix(fit$coef)), length(tau), n.boot))
          for(i in 1:n.boot) {
            id<-sample(1:n, size=n, replace=TRUE)
            coef.boot[,,i]<-ncross.rq.fitX(y=Y[id], X=X[id,,drop=FALSE], taus=tau, eps=eps)$coef
            
                  } #end i=1,..,n.boot
          } #end if n.boot


      #class(fit)<-"cgrq"
    } else {
    #se ci sono..
    drop.id<-lambda<-S<-B<-BB<-nomiCoefPEN<-nomiVariabPEN<-NULL
    if(length(testo.ps)>0){
#        z<-with(mf,eval(parse(text=testo.ps))) #questo va bene se ce ne solo uno di testo.ps
        l<-lapply(testo.ps,function(xx)with(mf,eval(parse(text=xx),envir=data)))
        vNdx<-sapply(l,function(xx){if(is.null(attr(xx,"ndx"))) round(min(40,n/4)) else attr(xx,"ndx")})
        #vNdx<-sapply(l,function(xx)attr(xx,"ndx"))
        #  f.Ndx<-if(is.null(vNdx[[1]])) vNdx<- round(min(40,n/4))
        vDeg<-sapply(l,function(xx)attr(xx,"deg"))
        vMonot<-sapply(l,function(xx)attr(xx,"monot"))
        vDiff<-sapply(l,function(xx)attr(xx,"pdiff"))
        lambda<-unlist(sapply(l,function(xx)attr(xx,"lambda")))
        nomeX<-unlist(sapply(l,function(xx)attr(xx,"nomeX")))
        var.pen<-unlist(sapply(l,function(xx)attr(xx,"var.pen"))) #c'e' bisogno di unlist()?
        #nomeBy<-unlist(sapply(l,function(xx)attr(xx,"nomeBy")))
        origName<-names(mf)
        nomiPS<-all.vars(formula)[c(1,match(nomeX,all.vars(formula)))] #estrae i nomi delle variabili nella formula in comune con nomeX
        #nomiPS<-all.vars(formula)[c(1,match(nomeBy,all.vars(formula)))] #estrae i nomi delle variabili nella formula in comune con nomeX
        names(mf)[c(1,id.ps)]<-nomiPS
        nomiPS<-nomiPS[-1]
        mVariabili<-mf[,colnames(mf)%in%nomeX,drop=FALSE] #uguale a mf[,id.ps,drop=FALSE]???
        if(ncol(mVariabili)!=length(l)) stop("A linear variable in ps()?")
        #byVariabili<-if(nomeBy!="NA") mf[,colnames(mf)%in%nomeBy,drop=FALSE] else matrix(1,nrow(mVariabili),ncol(mVariabili))
        #S<-B<-vector("list",ncol(mVariabili))
        nomiCoefPEN<-B<-BB<-Bderiv<- vector(length=ncol(mVariabili) , "list")
        for(j in 1:ncol(mVariabili)) {
          B[[j]]<- bspline(mVariabili[,j],ndx=vNdx[j],deg=vDeg[j]) 
          #---------?????????????????????????
          #S[[length(S)+1]]<- if(!vCyc[j]) crossprod(diff(diag(ncol(B[[j]])),diff=vDiff[j]))
          #per i disegni..
          xdisegno<-seq(min(mVariabili[,j]),max(mVariabili[,j]), l=100)
          BB1<-bspline(xdisegno,ndx=vNdx[j], deg=vDeg[j])
          Bderiv[[j]]<-bspline(xdisegno,ndx=vNdx[j], deg=vDeg[j], deriv=1)
          attr(BB1,"covariate.n")<- mVariabili[,j] #NB mVariabili[,j] (che viene assegnato a attr(,"covariate.n")) contiene altri attributi "ndx", "deg", "pdiff", "monot", "lambda","nomeX"
          attr(BB1,"covariate.35")<- xdisegno
          BB[[j]]<-BB1 
          nomiCoefPEN[[j]]<-paste(nomiPS[j],"ps",1:ncol(B[[j]]),sep=".")
              }
        nomiVariabPEN<-nomiPS
        drop.id<-id.ps-1+intercMt
        } #end length(testo.ps)>0
#    if(length(testo.ridge)>0){
#.......[SNIP].......
#        }
    X<- X[,-match(testo.ps,colnames(X)), drop=FALSE]
    if(!interc && ("(Intercept)"%in%colnames(X))) X<- X[,-match("(Intercept)",colnames(X)), drop=FALSE]
    if(interc && !("(Intercept)"%in%colnames(X))) X<-cbind("(Intercept)"=1,X)
    #id.interc<-match("(Intercept)",colnames(X))
    #X<-if(is.na(id.interc)) X[,-drop.id, drop=FALSE] else X[,-c(id.interc,drop.id), drop=FALSE]
    #if(interc) X<-cbind("(Intercept)"=1,X)
#    if(!is.null(lambda)&&(length(lambda)!=length(S))) stop("lambda tutti noti o ignoti")
    if(length(id.ps)>1 && (length(lambda)!=length(id.ps))) stop("a unique 'lambda' should be set in ps() with several smooth terms")
    if(length(lambda)>1 && length(id.ps)==1){
      lambdas<-lambda
      r.cv<-gcrq.rq.cv(Y, B[[1]], X, tau, interc, vMonot, vNdx, lambda, vDeg, vDiff, var.pen, cv, nfolds, foldid, eps=eps)
      lambda<-r.cv[[1]]
    } else {
      cv<-FALSE
    }
    names(lambda)<-nomiVariabPEN
#    if(ncol(X)<=0){
#      #NB: B is a list now! vMonot, and vNdx and vDeg vDiff are vectors!
#      fit<-ncross.rq.fitB(y=Y, B=B, taus=tau, monotone=vMonot, ndx=vNdx, lambda=lambda, deg=vDeg,
#        dif=vDiff, var.pen=var.pen, eps=eps)
#      if(n.boot>0){
#          coef.boot<-array(, dim=c(nrow(fit$coef), ncol(fit$coef), n.boot))
#          for(i in 1:n.boot) {
#            id<-sample(1:n, size=n, replace=TRUE)
#            coef.boot[,,i]<-ncross.rq.fitB(y=Y[id], B=NULL, x=mVariabili[id,,drop=FALSE], taus=tau, monotone=vMonot, ndx=vNdx, lambda=lambda, 
#                  deg=vDeg, dif=vDiff, var.pen=var.pen, eps=eps)$coef
#                  } #end i=1,..,n.boot
#          } #end if n.boot
#        } else {}
      fit<-ncross.rq.fitXB(y=Y, B=B, X=X, taus=tau, interc=interc, monotone=vMonot, ndx=vNdx, lambda=lambda, deg=vDeg,
        dif=vDiff, var.pen=var.pen, eps=eps)
      #boot
      if(n.boot>0){
          coef.boot<-array(, dim=c(nrow(as.matrix(fit$coef)), length(tau), n.boot))
          for(i in 1:n.boot) {
            id<-sample(1:n, size=n, replace=TRUE)
            coef.boot[,,i]<-ncross.rq.fitXB(y=Y[id], B=NULL, x=mVariabili[id,,drop=FALSE], X=X[id,, drop=FALSE], taus=tau, 
                interc=interc, monotone=vMonot, ndx=vNdx, lambda=lambda, deg=vDeg, dif=vDiff, var.pen=var.pen, eps=eps)$coef
                  } #end i=1,..,n.boot
          } #end if n.boot
      nomiCoefUNPEN<-rownames(as.matrix(fit$coefficients))[1:ncol(X)]
#      }
    nn<-c(nomiCoefUNPEN,unlist(nomiCoefPEN))
    #if(ncol(X)>0) nn<-c("(Intercept)",nn)
    if(length(tau)>1) rownames(fit$coefficients)<-nn else names(fit$coefficients)<-nn
    names(BB)<-nomiVariabPEN
    fit$info.smooth<-list(monotone=vMonot, ndx=vNdx, lambda=lambda, deg=vDeg, dif=vDiff)
    #names(BB)<-names(fit$lambda)<-names(fit$edf.smooth)<-nomiVariabPEN
    #fit$B<-B  #n righe..
    fit$BB<-BB #35 righe e attr "covariate.n" "covariate.35". #NB attr(,"covariate.n") contiene altri attributi "ndx", "deg", "pdiff", "monot", "lambda","nomeX"
    fit$Bderiv<-Bderiv
    } #fine del "se ci sono termini smooth"
    #--------------------------------------
    if(y) fit$y<-Y
    fit$contrasts <- attrContr  
    colnames(mf)[id.ps]<-testo.ps #devi sostituire i nome altrimenti .getXlevels() non funziona
    #non capisco (o non ricordo) perche' lo avevo messo..
    fit$xlevels <- .getXlevels(mt, mf) 
    fit$taus<-tau
    fit$call<-call
    if(n.boot){ 
      fit$boot.coef<- coef.boot
      }
    if(cv) {
        fit$cv <- cbind(lambdas,r.cv[[2]])
        colnames(fit$cv)<-c("lambdas",paste("fold",1:ncol(r.cv[[2]]), sep=""))
        #fit$foldid<-attr(r.cv[[2]], "foldid")
        attr(fit$cv, "foldid")<-attr(r.cv[[2]], "foldid")
        }
    class(fit)<-c("gcrq")
    fit
    }
