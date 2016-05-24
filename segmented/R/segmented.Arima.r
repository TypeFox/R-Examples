#`segmented.default` <-
#o1<-segmented(out.lm, seg.Z=~x +z,psi=list(x=c(30,60),z=.3), control=seg.control(display=FALSE, n.boot=20, seed=1515))
#o2<-ss(out.lm, seg.Z=~x +z,psi=list(x=c(30,60),z=.3), control=seg.control(display=FALSE, n.boot=20, seed=1515))

#o2<-ss(out.lm, seg.Z=~x +z,psi=list(x=c(30,60),z=.3), control=seg.control(display=FALSE, n.boot=0))
#o2<-ss(o, seg.Z=~age, psi=41, control=seg.control(display=FALSE, n.boot=0))

segmented.Arima<-
function(obj, seg.Z, psi, control = seg.control(), model = TRUE, ...) {
#Richiede control$f.obj that should be a string like "sum(x$residuals^2)" or "x$dev"
#-----------------
dpmax<-function(x,y,pow=1){
#deriv pmax
        if(pow==1) ifelse(x>y, -1, 0)
         else -pow*pmax(x-y,0)^(pow-1)
         }
#-----------
    n.Seg<-1
    if(missing(psi)){if(length(all.vars(seg.Z))>1) stop("provide psi") else psi<-Inf}
    if(length(all.vars(seg.Z))>1 & !is.list(psi)) stop("`psi' should be a list with more than one covariate in `seg.Z'")
    if(is.list(psi)){
      if(length(all.vars(seg.Z))!=length(psi)) stop("A wrong number of terms in `seg.Z' or `psi'")
      if(any(is.na(match(all.vars(seg.Z),names(psi), nomatch = NA)))) stop("Variables in `seg.Z' and `psi' do not match")
      n.Seg <- length(psi)
      }
    if(length(all.vars(seg.Z))!=n.Seg) stop("A wrong number of terms in `seg.Z' or `psi'")
    it.max <- old.it.max<- control$it.max
    toll <- control$toll
    visual <- control$visual
    stop.if.error<-control$stop.if.error
    n.boot<-control$n.boot
#    n.boot<-0
    size.boot<-control$size.boot
    gap<-control$gap
    random<-control$random
    pow<-control$pow
    visualBoot<-FALSE
    if(n.boot>0){
        if(!is.null(control$seed)) {
            set.seed(control$seed)
            employed.Random.seed<-control$seed
              } else {
            employed.Random.seed<-eval(parse(text=paste(sample(0:9, size=6), collapse="")))
            set.seed(employed.Random.seed)
              }
        if(visual) {visual<-FALSE; visualBoot<-TRUE}# warning("`display' set to FALSE with bootstrap restart", call.=FALSE)}
        if(!stop.if.error) stop("Bootstrap restart only with a fixed number of breakpoints")
     }
    last <- control$last
    K<-control$K
    h<-min(abs(control$h),1)
    if(h<1) it.max<-it.max+round(it.max/2)
    name.Z <-all.vars(seg.Z)
    if(length(name.Z)!=n.Seg) stop("errore strano 1")
    Z<-sapply(name.Z, function(xx) eval(parse(text=xx))) #e' sempre una matrice
    if(length(name.Z)!=ncol(Z)) stop("errore strano 2")
    n<-nrow(Z)
    n.psi<- length(unlist(psi))

    #################
    if(ncol(Z)==1 && length(psi)==1 && n.psi==1 && !any(is.na(psi))) { if(psi==Inf) psi<-median(Z)}
    #################

    if(ncol(Z)==1 && is.vector(psi) && (is.numeric(psi)||is.na(psi))){
        psi <- list(as.numeric(psi))
        names(psi)<-name.Z
        }
    if (!is.list(psi) || is.null(names(psi))) stop("psi should be a *named* list")
    id.nomiZpsi <- match(colnames(Z), names(psi))
    if ((ncol(Z)!=length(psi)) || any(is.na(id.nomiZpsi))) stop("Length or names of Z and psi do not match")
    nome <- names(psi)[id.nomiZpsi]
    psi <- psi[nome]
    initial.psi<-psi
    for(i in 1:length(psi)) {
        if(any(is.na(psi[[i]]))) psi[[i]]<-if(control$quant) {quantile(Z[,i], prob= seq(0,1,l=K+2)[-c(1,K+2)], names=FALSE)} else {(min(Z[,i])+ diff(range(Z[,i]))*(1:K)/(K+1))}
        }
    a <- sapply(psi, length)
    #per evitare che durante il processo iterativo i psi non siano ordinati
    id.psi.group <- rep(1:length(a), times = a) #identificativo di apparteneza alla variabile

    Z<-matrix(unlist(mapply(function(x,y)rep(x,y),Z,a,SIMPLIFY = TRUE)),nrow=n,byrow = TRUE)
    #negli altri metodi Z e' una lista per cui la linea di sopra diventa
    #Z<-matrix(unlist(mapply(function(x,y)rep(x,y),Z,a,SIMPLIFY = TRUE)),nrow=n)
    colnames(Z) <- nomiZ.vett <- rep(nome, times = a) #SERVE??? si perche' Z e' senza colnames
    
    psi <- unlist(psi)
    #se psi e' numerico, la seguente linea restituisce i valori ordinati all'interno della variabile..
    psi<-unlist(tapply(psi,id.psi.group,sort))
    k <- ncol(Z)
    PSI <- matrix(rep(psi, rep(n, k)), ncol = k)
    
    #controllo se psi e' ammissibile..
    c1 <- apply((Z <= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo <)
    c2 <- apply((Z >= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo >)
    if(sum(c1 + c2) != 0 || is.na(sum(c1 + c2)) ) stop("starting psi out of the admissible range")

    #ripetizioni <- as.numeric(unlist(sapply(table(nomiZ)[order(unique(nomiZ))], function(xxx) {1:xxx})))
    ripetizioni <- as.vector(unlist(tapply(id.psi.group, id.psi.group, function(x) 1:length(x) )))
    nomiU <- paste("U", ripetizioni, sep = "")
    nomiU <- paste(nomiU, nomiZ.vett, sep = ".")
    nomiV <- paste("V", ripetizioni, sep = "")
    nomiV <- paste(nomiV, nomiZ.vett, sep = ".")
    nnomi <- c(nomiU, nomiV)

#    U <- pmax((Z - PSI), 0)^pow[1]#U <- pmax((Z - PSI), 0)
#    #V <- dpmax(Z,PSI,pow=pow[2])# ifelse((Z > PSI), -1, 0)
#    V<-ifelse((Z > PSI), -1, 0)
#
#    for(i in 1:k) {
#        mfExt[nomiU[i]] <- U[,i]
#        mfExt[nomiV[i]] <- V[,i]
#        }
#    Fo <- update.formula(formula(obj), as.formula(paste(".~.+", paste(nnomi, collapse = "+"))))
#    Fo.noV <- update.formula(formula(obj), as.formula(paste(".~.+", paste(nomiU, collapse = "+"))))
#    call.ok <- update(obj, formula = Fo,  evaluate=FALSE, data = mfExt) #objF <- update(obj0, formula = Fo, data = KK)
#    call.noV <- update(obj, formula = Fo.noV,  evaluate=FALSE, data = mfExt) #objF <- update(obj0, formula = Fo, data = KK)
    XREG<-eval(obj$call$xreg)
    nomiXREG<-setdiff(names(obj$coef),c("intercept", paste("ar",1:100,sep=""), paste("ma",1:100,sep="")))
    XREG<-matrix(XREG, ncol=length(nomiXREG))
    colnames(XREG)<-nomiXREG
    mio.init<-mio.init.noV<-NULL
    call.ok <- update(obj,  xreg = cbind(XREG,U,V), init=mio.init, evaluate=FALSE) #, data = mfExt) #objF <- update(obj0, formula = Fo, data = KK)
    call.noV <- update(obj, xreg = cbind(XREG,U), init=mio.init.noV,  evaluate=FALSE) #, data = mfExt) #objF <- update(obj0, formula = Fo, data = KK)
#    call.noV <- update(obj, formula = Fo.noV,  evaluate=FALSE, data = mfExt) #objF <- update(obj0, formula = Fo, data = KK)

    if (it.max == 0) {
      obj1 <- eval(call.noV) #, envir=mfExt)
      return(obj1)
    }

    #obj1 <- eval(call.ok, envir=mfExt)
    initial <- psi
    obj0 <- obj
    
#    browser()
    
    dev0<- -obj$loglik
    if(is.na(dev0)) dev0<-10
    
    list.obj <- list(obj)
    nomiOK<-nomiU

    opz<-list(toll=toll,h=h,stop.if.error=stop.if.error,dev0=dev0,visual=visual,it.max=it.max,
        nomiOK=nomiOK, id.psi.group=id.psi.group, gap=gap, visualBoot=visualBoot, pow=pow)

    opz$call.ok<-call.ok
    opz$call.noV<-call.noV
    opz$nomiU<-nomiU
    opz$nomiV<-nomiV
#    opz$fn.obj <- fn.obj    

    if(n.boot<=0){
      obj<-seg.Ar.fit(obj, XREG, Z, PSI, opz)
      } else {
      obj<-seg.Ar.fit.boot(obj, XREG, Z, PSI, opz, n.boot=n.boot, size.boot=size.boot, random=random) #jt, nonParam
      }

    if(!is.list(obj)){
        warning("No breakpoint estimated", call. = FALSE)
        return(obj0)
        }
    if(!is.null(obj$obj$df.residual)){
      if(obj$obj$df.residual==0) warning("no residual degrees of freedom (other warnings expected)", call.=FALSE)
      }
    id.psi.group<-obj$id.psi.group
    nomiOK<-obj$nomiOK #sarebbe nomiU
#--
    nomiVxb<-paste("psi",sapply(strsplit(nomiOK,"U"), function(x){x[2]}), sep="")
    #nomiFINALI<-unique(sapply(strsplit(nomiOK, split="[.]"), function(x)x[2])) #nomi delle variabili con breakpoint stimati!
    nomiFINALI<-unique(sub("U[1-9]*[0-9].", "", nomiOK))
    #se e' stata usata una proc automatica "nomiFINALI" sara' differente da "name.Z"
    nomiSenzaPSI<-setdiff(name.Z,nomiFINALI)
    if(length(nomiSenzaPSI)>=1) warning("no breakpoints found for: ", paste(nomiSenzaPSI," "), call. = FALSE)

#--
    it<-obj$it
    psi<-obj$psi
    psi.values<-if(n.boot<=0) obj$psi.values else obj$boot.restart
    U<-obj$U
    V<-obj$V
#    return(obj)

    #if(any(table(rowSums(V))<=1)) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close")
    for(jj in colnames(V)) {
        VV<-V[, which(colnames(V)==jj), drop=FALSE]
        sumV<-abs(rowSums(VV))
        if( #(any(diff(sumV)>=2)|| #se ci sono due breakpoints uguali
            any(table(sumV)<=1) && stop.if.error) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close each other")
        }
    rangeZ<-obj$rangeZ
    obj<-obj$obj
    k<-length(psi)
#    beta.c<-if(k == 1) coef(obj)["U"] else coef(obj)[paste("U", 1:ncol(U), sep = "")]

    beta.c<- coef(obj)[nomiU]
    psi.values[[length(psi.values) + 1]] <- psi
    id.warn <- FALSE
    if (n.boot<=0 && it > it.max) { #it >= (it.max+1)
        warning("max number of iterations attained", call. = FALSE)
        id.warn <- TRUE
    }
    Vxb <- V %*% diag(beta.c, ncol = length(beta.c))

#    #se usi una procedura automatica devi cambiare ripetizioni, nomiU e nomiV, e quindi:
#    length.psi<-tapply(as.numeric(as.character(names(psi))), as.numeric(as.character(names(psi))), length)
#    forma.nomiU<-function(xx,yy)paste("U",1:xx, ".", yy, sep="")
#    forma.nomiVxb<-function(xx,yy)paste("psi",1:xx, ".", yy, sep="")
#    nomiU   <- unlist(mapply(forma.nomiU, length.psi, name.Z)) #in realta' non serve, c'era gia'!
#    nomiVxb <- unlist(mapply(forma.nomiVxb, length.psi, name.Z))
#    nnomi <- c(nomiU, nomiVxb)
#    
#    colnames(U)<-nomiU
#    colnames(Vxb)<-nomiVxb
#
##    for(i in 1:ncol(U)) {
#        mfExt[nomiU[i]]<-mf[nomiU[i]]<-U[,i]
#        mfExt[nomiVxb[i]]<-mf[nomiVxb[i]]<-Vxb[,i]
#        }


    nnomi <- c(nomiU, nomiVxb)
#    browser()
#    Fo <- update.formula(formula(obj0), as.formula(paste(".~.+", paste(nnomi, collapse = "+"))))
#    objF <- update(obj0, formula = Fo,  evaluate=FALSE, data = mfExt)
#    objF<- eval(objF, envir=mfExt)
#browser()

#    objF <- update(obj0,  xreg = cbind(XREG,U,Vxb), evaluate=TRUE) #, data = mfExt) #objF <- update(obj0, formula = Fo, data = KK)
    XREG.ok<-cbind(XREG, U, Vxb)
    colnames(XREG.ok)[((ncol(XREG.ok)-length(nnomi)+1):ncol(XREG.ok))]<- nnomi
    objF <- update(obj0,  xreg = XREG.ok, evaluate=TRUE) #, data = mfExt) #objF <- update(obj0, formula = Fo, data = KK)

    #Puo' capitare che psi sia ai margini e ci sono 1 o 2 osservazioni in qualche intervallo. Oppure ce ne
    #sono di piu' ma hanno gli stessi valori di x
    #objF$coef puo' avere mancanti.. names(which(is.na(coef(objF))))

    if(any(is.na(objF$coef)) && stop.if.error){
     stop("at least one coef estimate is NA: breakpoint(s) at the boundary? (possibly with many x-values replicated)", call. = FALSE)
    }
# CONTROLLARE!!!!
#    objF$offset<- obj0$offset
    #sostituire i valori: objF include le U e V, obj solo le U
    if(!gap){
     #names.coef <- names(objF$coefficients)
     #names(obj$coefficients)[match(nomiV, names(coef(obj)))]<-nomiVxb
     #objF$coefficients[names.coef]<-obj$coefficients[names.coef]
     names.coef <- names(obj$coef)
     objF$coef[names.coef]<-obj$coef[names.coef]
     objF$coef[nomiVxb]<-rep(0, k)
     #if(!is.null(objF$fitted.values)) objF$fitted.values<-obj$fitted.values
     if(!is.null(objF$residuals)) objF$residuals<-obj$residuals
     if(!is.null(objF$weights)) objF$weights<-obj$weights
     if(!is.null(objF$aic)) objF$aic<-obj$aic + 2*k
     }
    if(any(is.na(objF$coef))){ #Se gap==FALSE qui non ci possono essere NA (sono sostituiti dagli 0)
    stop("some estimate is NA: premature stopping with a large number of breakpoints?",
      call. = FALSE)
      }


    Cov <-  try(vcov(objF), silent=TRUE)
    idd <- match(nomiVxb, names(coef(objF)))
    if(class(Cov)!="try-error") {
        vv <- if (length(idd) == 1) Cov[idd, idd] else diag(Cov[idd, idd])
        } else {
      vv<-NA
        }
    a<-tapply(id.psi.group, id.psi.group, length) #ho sovrascritto "a" di sopra, ma non dovrebbe servire..
    ris.psi<-matrix(,length(psi),3)
    colnames(ris.psi) <- c("Initial", "Est.", "St.Err")
    rownames(ris.psi) <- nomiVxb
    ris.psi[,2]<-psi
    ris.psi[,3]<-sqrt(vv)
    a.ok<-NULL
    for(j in name.Z){
        if(j %in% nomiFINALI) {
          a.ok[length(a.ok)+1]<-a[1]
          a<-a[-1]
          } else {
          a.ok[length(a.ok)+1]<-0
          } #ifelse(name.Z %in% nomiFINALI,1,0)
        }
#    initial<-unlist(mapply(function(x,y){if(is.na(x)[1])rep(x,y) else x }, initial.psi, a.ok, SIMPLIFY = TRUE))
    initial<-unlist(mapply(function(x,y){if(is.na(x)[1])rep(x,y) else x }, initial.psi[nomiFINALI], a.ok[a.ok!=0], SIMPLIFY = TRUE))
    ris.psi[,1]<-initial

#    a<-tapply(id.psi.group, id.psi.group, length) #ho sovrascritto "a" di sopra, ma non dovrebbe servire..
#    initial<-unlist(mapply(function(x,y){if(is.na(x)[1])rep(x,y) else x }, initial.psi, a))
#    id <- match(nomiVxb, names(coef(objF)))
#    Cov <-  try(vcov(objF), silent=TRUE)
#    if(class(Cov)!="try-error") {
#        vv <- if (length(id) == 1) Cov[id, id] else diag(Cov[id, id])
#        #if(length(initial)!=length(psi)) initial<-rep(NA,length(psi))
#        psi <- cbind(initial, psi, sqrt(vv))
#        rownames(psi) <- colnames(Cov)[id]
#        colnames(psi) <- c("Initial", "Est.", "St.Err")
#        } else {
#        psi <- cbind(initial, psi)
#        rownames(psi) <- nomiVxb
#        colnames(psi) <- c("Initial", "Est.")
#        }
    objF$Z <- Z
    objF$rangeZ <- rangeZ
    objF$psi.history <- psi.values
    objF$psi <- ris.psi
    objF$it <- (it - 1)
    objF$epsilon <- obj$epsilon
    objF$call <- match.call()
    #objF$nameUV <- list(U = nomiU, V = rownames(psi), Z = name.Z)
    objF$nameUV <- list(U = drop(nomiU), V = rownames(ris.psi), Z = nomiFINALI) #Z = name.Z

    objF$id.group <- if(length(name.Z)<=1) -rowSums(as.matrix(V))
    objF$id.psi.group <- id.psi.group
    objF$id.warn <- id.warn

#    browser()
#    objF$orig.call<-orig.call
#    if (model)  objF$model <- mf #objF$mframe <- data.frame(as.list(KK))

    if(n.boot>0) objF$seed<-employed.Random.seed
    class(objF) <- c("segmented", class(obj0))
    list.obj[[length(list.obj) + 1]] <- objF
    class(list.obj) <- "segmented"
    if (last) list.obj <- list.obj[[length(list.obj)]]
    warning("'segmented.Arima' is at a preliminary stage. Estimates are OK, but the '*.segmented' methods are not expected to work",
      call.=FALSE)
    return(list.obj)
    } #end function
