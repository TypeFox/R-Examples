#`segmented.default` <-
#o1<-segmented(out.lm, seg.Z=~x +z,psi=list(x=c(30,60),z=.3), control=seg.control(display=FALSE, n.boot=20, seed=1515))
#o2<-ss(out.lm, seg.Z=~x +z,psi=list(x=c(30,60),z=.3), control=seg.control(display=FALSE, n.boot=20, seed=1515))

#o2<-ss(out.lm, seg.Z=~x +z,psi=list(x=c(30,60),z=.3), control=seg.control(display=FALSE, n.boot=0))
#o2<-ss(o, seg.Z=~age, psi=41, control=seg.control(display=FALSE, n.boot=0))

segmented.default<-
function(obj, seg.Z, psi, control = seg.control(), model = TRUE, ...) {
#Richiede control$f.obj that should be a string like "sum(x$residuals^2)" or "x$dev"
#-----------------
dpmax<-function(x,y,pow=1){
#deriv pmax
        if(pow==1) ifelse(x>y, -1, 0)
         else -pow*pmax(x-y,0)^(pow-1)
         }
#-----------
#    control$fn.obj<-"sum(x$residuals^2)"
#    control$fn.obj<-"x$dev"
#    control$fn.obj<-"-x$loglik[2]"
#    control$fn.obj<-"x$rho"
#    if(is.null(control$fn.obj)) stop("'segmented.default' needs 'fn.obj' specified in seg.control") else fn.obj<-control$fn.obj
    if(is.null(control$fn.obj)) fn.obj<-"-as.numeric(logLik(x))" else fn.obj<-control$fn.obj

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
    orig.call<-Call<-mf<-obj$call
    orig.call$formula<- mf$formula<-formula(obj) #per consentire lm(y~.)
    m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    if(class(mf$formula)=="name" && !"~"%in%paste(mf$formula)) mf$formula<-eval(mf$formula)
    mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
    mfExt<- mf
    if(!is.null(obj$call$offset) || !is.null(obj$call$weights) || !is.null(obj$call$subset) || !is.null(obj$call$id)){ 
      mfExt$formula <- 
          update.formula(mf$formula, 
          paste(".~.+", paste(
          c(all.vars(obj$call$offset), 
            all.vars(obj$call$weights),
            all.vars(obj$call$subset),
            all.vars(obj$call$id)), 
            collapse = "+")
            ))
          }
    mf <-  eval(mf, parent.frame())
    n<-nrow(mf)
    #questo serve per inserire in mfExt le eventuali variabili contenute nella formula con offset(..)
    nomiOff<-setdiff(all.vars(formula(obj)), names(mf))
    if(length(nomiOff)>=1) mfExt$formula<-update.formula(mfExt$formula,paste(".~.+", paste( nomiOff, collapse="+"), sep=""))
    nomiTUTTI<-all.vars(mfExt$formula) #comprende anche altri nomi (ad es., threshold) "variabili"
    nomiNO<-NULL #dovrebbe contenere
    for(i in nomiTUTTI){
      r<-try(eval(parse(text=i), parent.frame()), silent=TRUE)
      if(class(r)!="try-error" && length(r)==1 && !is.function(r)) nomiNO[[length(nomiNO)+1]]<-i
      }
    if(!is.null(nomiNO)) mfExt$formula<-update.formula(mfExt$formula,paste(".~.-", paste( nomiNO, collapse="-"), sep=""))
    mfExt<-eval(mfExt, parent.frame())
    #apply(mfExt,2,function(x) {if(is.Surv(x)) x[,1:ncol(x)] else x}) 
    if(inherits(obj, "coxph")){
      is.Surv<-NA
      rm(is.Surv)
      for(i in 1:ncol(mfExt)){
            if(is.Surv(mfExt[,i])) aa<-mfExt[,i][,1:ncol(mfExt[,i])]
               }
      mfExt<-cbind(aa,mfExt)
    }      
 
    id.seg<-match(all.vars(seg.Z), names(mfExt))
    name.Z<-names(mfExt)[id.seg]
    Z<-mfExt[,id.seg,drop=FALSE]
#    name.Z <- names(Z)
    
    n.psi<-length(unlist(psi))
    #################
    if(ncol(Z)==1 && length(psi)==1 && n.psi==1 && !any(is.na(psi))) { if(psi==Inf) psi<-median(Z[,1])} #devi selezionare la colonna perche' Z e' un dataframe!
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

    Z<-matrix(unlist(mapply(function(x,y)rep(x,y),Z,a,SIMPLIFY = TRUE)),nrow=n)
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

#CONTROLLARE se e' necessario aggiungere a mfExt le nuove variabili U e V. Non dovrebbe perche'
#   update(, evaluate=FALSE) non "le vuole" e poi comunque vengono calcolate in seg.def.fit()
#   Servono solo se it.max=0..
#----------------------------------------------------------
    U <- pmax((Z - PSI), 0)^pow[1]#U <- pmax((Z - PSI), 0)
    #V <- dpmax(Z,PSI,pow=pow[2])# ifelse((Z > PSI), -1, 0)
    V<-ifelse((Z > PSI), -1, 0)
    for(i in 1:k) {
        mfExt[nomiU[i]] <- U[,i]
        mfExt[nomiV[i]] <- V[,i]
        }

    Fo <- update.formula(formula(obj), as.formula(paste(".~.+", paste(nnomi, collapse = "+"))))
    Fo.noV <- update.formula(formula(obj), as.formula(paste(".~.+", paste(nomiU, collapse = "+"))))
    
    call.ok <- update(obj, formula = Fo,  evaluate=FALSE, data = mfExt) #objF <- update(obj0, formula = Fo, data = KK)
    call.noV <- update(obj, formula = Fo.noV,  evaluate=FALSE, data = mfExt) #objF <- update(obj0, formula = Fo, data = KK)

    if (it.max == 0) {
      if(!is.null(call.noV[["subset"]])) call.noV[["subset"]]<-NULL
      obj1 <- eval(call.noV, envir=mfExt)
      return(obj1)
    }

    #obj1 <- eval(call.ok, envir=mfExt)
    initial <- psi
    obj0 <- obj
    
    #browser()
    
    dev0<- eval(parse(text=fn.obj), list(x=obj))
    if(length(dev0)<=0) stop("error in the objective to be minimized, see 'fn.obj' in ?seg.control") #per fit su cui logLik() does not work..
    if(length(dev0)>1) stop("the objective to be minimized is not scalar, see 'fn.obj' in ?seg.control") #per fit su cui logLik() does not work..
    if(is.na(dev0)) dev0<-10
    
    list.obj <- list(obj)
    nomiOK<-nomiU

    opz<-list(toll=toll,h=h,stop.if.error=stop.if.error,dev0=dev0,visual=visual,it.max=it.max,
        nomiOK=nomiOK, id.psi.group=id.psi.group, gap=gap, visualBoot=visualBoot, pow=pow)

    opz$call.ok<-call.ok
    opz$call.noV<-call.noV
    opz$formula.orig<-formula(obj)
    opz$nomiU<-nomiU
    opz$nomiV<-nomiV
    opz$fn.obj <- fn.obj    

    if(n.boot<=0){
      obj<-seg.def.fit(obj, Z, PSI, mfExt, opz)
      } else {
      obj<-seg.def.fit.boot(obj, Z, PSI, mfExt, opz, n.boot=n.boot, size.boot=size.boot, random=random) #jt, nonParam
      }
#browser()

    if(!is.list(obj)){
        warning("No breakpoint estimated", call. = FALSE)
        return(obj0)
        }
    if(!is.null(obj$obj$df.residual)){
      if(obj$obj$df.residual==0) warning("no residual degrees of freedom (other warnings expected)", call.=FALSE)
      }
    id.psi.group<-obj$id.psi.group
    nomiOK<-obj$nomiOK #sarebbe nomiU
    nomiVxb<-paste("psi",sapply(strsplit(nomiOK,"U"), function(x){x[2]}), sep="")
    #nomiFINALI<-unique(sapply(strsplit(nomiOK, split="[.]"), function(x)x[2])) #nomi delle variabili con breakpoint stimati!
    nomiFINALI<-unique(sub("U[1-9]*[0-9].", "", nomiOK))
    #se e' stata usata una proc automatica "nomiFINALI" sara differente da "name.Z"
    nomiSenzaPSI<-setdiff(name.Z,nomiFINALI)
    if(length(nomiSenzaPSI)>=1) warning("no breakpoints found for: ", paste(nomiSenzaPSI," "), call. = FALSE)
    
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
#        if( (any(diff(sumV)>=2) #se ci sono due breakpoints uguali
#            || any(table(sumV)<=1)) && stop.if.error) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close each other")
        if(any(table(sumV)<=1) && stop.if.error) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close each other")
        }
    rangeZ<-obj$rangeZ
    mfExt<-obj$mfExt
    names(mfExt)[match(obj$nomiV, names(mfExt))]<-nomiVxb
    obj<-obj$obj
    k<-length(psi)
#    beta.c<-if(k == 1) coef(obj)["U"] else coef(obj)[paste("U", 1:ncol(U), sep = "")]
#browser()
    beta.c<- coef(obj)[nomiOK] #nomiOK e' stato estratto da obj e contiene tutti i nomi delle variabili U inserite nel modello
    psi.values[[length(psi.values) + 1]] <- psi
    id.warn <- FALSE
    if (n.boot<=0 && it > it.max) { #it >= (it.max+1)
        warning("max number of iterations attained", call. = FALSE)
        id.warn <- TRUE
    }
    Vxb <- V %*% diag(beta.c, ncol = length(beta.c))

    #se usi una procedura automatica devi cambiare ripetizioni, nomiU e nomiV, e quindi:
#    length.psi<-tapply(as.numeric(as.character(names(psi))), as.numeric(as.character(names(psi))), length)
#    forma.nomiU<-function(xx,yy)paste("U",1:xx, ".", yy, sep="")
#    forma.nomiVxb<-function(xx,yy)paste("psi",1:xx, ".", yy, sep="")
#    nomiU   <- unlist(mapply(forma.nomiU, length.psi, nomiFINALI)) #invece di un ciclo #paste("U",1:length.psi[i], ".", name.Z[i])
#    nomiVxb <- unlist(mapply(forma.nomiVxb, length.psi, nomiFINALI))
#    for(i in 1:ncol(U)) {
#        mfExt[nomiU[i]]<-mf[nomiU[i]]<-U[,i]
#        mfExt[nomiVxb[i]]<-mf[nomiVxb[i]]<-Vxb[,i]
#        }

#    ToDeletenomiU<-nomiU[!id.psi.ok] #salva i nomi delle U per i psi ammissibili
#    ToDeletenomiV<-nomiV[!id.psi.ok] #salva i nomi delle V per i psi ammissibili
#    if(length(ToDeletenomiU)>0 || length(ToDeletenomiV)>0) for(nn in c(ToDeletenomiU, ToDeletenomiV)) mfExt[[nn]]<-NULL

#E' inutile lavorare sul mf, utilizzo quello restituito da seg.def.fit.r
#    for(i in 1:k) {
#        mfExt[nomiOK[i]]<-mf[nomiOK[i]]<-U[,i]
#        mfExt[nomiVxb[i]]<-mf[nomiVxb[i]]<-Vxb[,i]
#        }

    nnomi <- c(nomiOK, nomiVxb)
    Fo <- update.formula(formula(obj0), as.formula(paste(".~.+", paste(nnomi, collapse = "+"))))
    objF <- update(obj0, formula = Fo,  evaluate=FALSE, data = mfExt)
    if(!is.null(objF[["subset"]])) objF[["subset"]]<-NULL
    objF<- eval(objF, envir=mfExt)
    #Puo' capitare che psi sia ai margini e ci sono 1 o 2 osservazioni in qualche intervallo. Oppure ce ne
    #sono di piu' ma hanno gli stessi valori di x
    #objF$coef puo' avere mancanti.. names(which(is.na(coef(objF))))
    objF$offset<- obj0$offset
    isNAcoef<-any(is.na(objF$coefficients))
    if(isNAcoef){
      if(stop.if.error) {stop("at least one coef is NA: breakpoint(s) at the boundary? (possibly with many x-values replicated)", 
        call. = FALSE)} else {
        warning("some estimate is NA: too many breakpoints? 'var(hat.psi)' cannot be computed \n ..returning a 'lm' model", call. = FALSE)
        Fo <- update.formula(formula(obj0), as.formula(paste(".~.+", paste(nomiU, collapse = "+"))))
        objF <- update(obj0, formula = Fo,  evaluate=TRUE, data = mfExt)
        names(psi)<-nomiVxb
        objF$psi<-psi
        return(objF)      
        }
    }
# CONTROLLARE!!!!
#    objF$offset<- obj0$offset
    #sostituire i valori: objF include le U e V, obj solo le U
    if(!gap){
     #names.coef <- names(objF$coefficients)
     #names(obj$coefficients)[match(nomiV, names(coef(obj)))]<-nomiVxb
     #objF$coefficients[names.coef]<-obj$coefficients[names.coef]
     names.coef <- names(obj$coefficients)
     objF$coefficients[names.coef]<-obj$coefficients[names.coef]
     objF$coefficients[nomiVxb]<-rep(0, k)
     if(!is.null(objF$geese$beta)) objF$geese$beta <- objF$coefficients#obj$geese$beta
     if(!is.null(objF$geese$gamma)) objF$geese$gamma <- obj$geese$gamma
     if(!is.null(objF$geese$alpha)) objF$geese$alpha <- obj$geese$alpha
     if(!is.null(objF$fitted.values)) objF$fitted.values<-obj$fitted.values
     if(!is.null(objF$residuals)) objF$residuals<-obj$residuals
     if(!is.null(objF$linear.predictors)) objF$linear.predictors<-obj$linear.predictors
     if(!is.null(objF$deviance)) objF$deviance<-obj$deviance
     if(!is.null(objF$weights)) objF$weights<-obj$weights
     if(!is.null(objF$aic)) objF$aic<-obj$aic + 2*k
     if(!is.null(objF$loglik)) objF$loglik<-obj$loglik #per coxph
     if(!is.null(objF$rho)) objF$rho<-obj$rho #per rq
     if(!is.null(objF$dual)) objF$dual<-obj$dual #per rq
     }
    Cov <-  try(vcov(objF), silent=TRUE)
    idd <- match(nomiVxb, names(coef(objF)))
    if(class(Cov)!="try-error") {
        vv <- if (length(idd) == 1) Cov[idd, idd] else diag(Cov[idd, idd])
        #if(length(initial)!=length(psi)) initial<-rep(NA,length(psi))
#        psi <- cbind(initial, psi, sqrt(vv))
#        rownames(psi) <- colnames(Cov)[idd]
#        colnames(psi) <- c("Initial", "Est.", "St.Err")
        } else {
#        psi <- cbind(initial, psi)
#        rownames(psi) <- nomiVxb
#        colnames(psi) <- c("Initial", "Est.")
      vv<-NA
        }

#browser()
    a<-tapply(id.psi.group, id.psi.group, length) #ho sovrascritto "a" di sopra, ma non dovrebbe servire..
    ris.psi<-matrix(,length(psi),3)
    colnames(ris.psi) <- c("Initial", "Est.", "St.Err")
    rownames(ris.psi) <- nomiVxb
    ris.psi[,2]<-psi
    ris.psi[,3]<-sqrt(vv)
#NB "a" deve essere un vettore che si appatta con "initial.psi" per ottnetere "initial" sotto... Se una variabile alla fine risulta
# senza breakpoint questo non avviene e ci sono problemi nella formazione di "initial". Allora costruisco a.ok
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
    objF$rangeZ <- rangeZ
    objF$psi.history <- psi.values
    objF$psi <- ris.psi
    objF$it <- (it - 1)
    objF$epsilon <- obj$epsilon
    objF$call <- match.call()
    objF$nameUV <- list(U = drop(nomiU), V = rownames(ris.psi), Z = nomiFINALI) #Z = name.Z
    objF$id.group <- if(length(name.Z)<=1) -rowSums(as.matrix(V))
    objF$id.psi.group <- id.psi.group
    objF$id.warn <- id.warn
    objF$orig.call<-orig.call
    if (model)  objF$model <- mf #objF$mframe <- data.frame(as.list(KK))
    if(n.boot>0) objF$seed<-employed.Random.seed
#    class(objF) <- c("segmented", class(obj0))
    list.obj[[length(list.obj) + 1]] <- objF
    class(list.obj) <- "segmented"
    if (last) list.obj <- list.obj[[length(list.obj)]]
    return(list.obj)
    } #end function
