`segmented.glm` <-
function(obj, seg.Z, psi, control = seg.control(), model = TRUE, ...) {
    n.Seg<-1
    if(missing(psi)){if(length(all.vars(seg.Z))>1) stop("provide psi") else psi<-Inf}
    if(length(all.vars(seg.Z))>1 & !is.list(psi)) stop("`psi' should be a list with more than one covariate in `seg.Z'")
    if(is.list(psi)){
      if(length(all.vars(seg.Z))!=length(psi)) stop("A wrong number of terms in `seg.Z' or `psi'")
      if(any(is.na(match(all.vars(seg.Z),names(psi), nomatch = NA)))) stop("Variables in `seg.Z' and `psi' do not match")
      n.Seg <- length(psi)
      }
    if(length(all.vars(seg.Z))!=n.Seg) stop("A wrong number of terms in `seg.Z' or `psi'")
    maxit.glm <- control$maxit.glm
    it.max <- old.it.max<- control$it.max
    toll <- control$toll
    visual <- control$visual
    stop.if.error<-control$stop.if.error
    n.boot<-control$n.boot
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
        if(visual) {visual<-FALSE; visualBoot<-TRUE}#warning("`display' set to FALSE with bootstrap restart", call.=FALSE)}
        if(!stop.if.error) stop("Bootstrap restart only with a fixed number of breakpoints")
     }
    last <- control$last
    K<-control$K
    h<-min(abs(control$h),1)
    if(h<1) it.max<-it.max+round(it.max/2)
#    if(!stop.if.error) objInitial<-obj        
    #-------------------------------
#    #una migliore soluzione.........
#    objframe <- update(obj, model = TRUE, x = TRUE, y = TRUE)
#    y <- objframe$y
#    a <- model.matrix(seg.Z, data = eval(obj$call$data))
#    a <- subset(a, select = colnames(a)[-1])
    orig.call<-Call<-mf<-obj$call
    orig.call$formula<-mf$formula<-formula(obj) #per consentire lm(y~.)
    m <- match(c("formula", "data", "subset", "weights", "na.action","etastart","mustart","offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    #non so a che serva la seguente linea..
    if(class(mf$formula)=="name" && !"~"%in%paste(mf$formula)) mf$formula<-eval(mf$formula)
    #orig.call$formula<-update.formula(orig.call$formula, paste("~.-",all.vars(seg.Z))) #utile per plotting
    
#    nomeRispo<-strsplit(paste(formula(obj))[2],"/")[[1]] #eventuali doppi nomi separati da "/" (tipo "y/n" per GLM binom)
    #la linea sotto aggiunge nel mf anche la variabile offs..
#    if(length(all.vars(formula(obj)))>1){
#      id.rispo<-1
#      if(length(nomeRispo)>=2) id.rispo<-1:2      
#      #questo serve quando formula(obj) ha solo l'intercept
#      agg<-if(length(all.vars(formula(obj))[-id.rispo])==0) "" else "+"
#      mf$formula<-update.formula(mf$formula,paste(paste(seg.Z,collapse=".+"),agg,paste(all.vars(formula(obj))[-id.rispo],collapse="+")))
#    } else {
    #  mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
#    }
    mfExt<- mf
    mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))

#    mfExt$formula<- update.formula(mfExt$formula,paste(paste(seg.Z,collapse=".+"),"+",paste(all.vars(formula(obj)),collapse="+")))   
#    mfExt$formula<- if(!is.null(obj$call$data)) 
#        update.formula(mf$formula,paste(".~",paste(all.vars(obj$call), collapse="+"),"-",obj$call$data,sep=""))
#            else update.formula(mf$formula,paste(".~",paste(all.vars(obj$call), collapse="+"),sep=""))
#-----------
   if(!is.null(obj$call$offset) || !is.null(obj$call$weights) || !is.null(obj$call$subset)){ 
      mfExt$formula <- 
          update.formula(mf$formula, 
          paste(".~.+", paste(
          c(all.vars(obj$call$offset), 
            all.vars(obj$call$weights),
            all.vars(obj$call$subset)), 
            collapse = "+")
            ))
          }
    mf <-  eval(mf, parent.frame())
    n<-nrow(mf)
    #La linea sotto serve per inserire in mfExt le eventuali variabili contenute nella formula con offset(..)
    #   o anche variabili che rientrano in espressioni (ad es., y/n o I(y*n))
    nomiOff<-setdiff(all.vars(formula(obj)), names(mf))
    if(length(nomiOff)>=1) mfExt$formula<-update.formula(mfExt$formula,paste(".~.+", paste( nomiOff, collapse="+"), sep=""))
    
    #ago 2014 c'e' la questione di variabili aggiuntive...
    nomiTUTTI<-all.vars(mfExt$formula) #comprende anche altri nomi (ad es., threshold) "variabili"
    nomiNO<-NULL #dovrebbe contenere
    for(i in nomiTUTTI){
      r<-try(eval(parse(text=i), parent.frame()), silent=TRUE)
      if(class(r)!="try-error" && length(r)==1 && !is.function(r)) nomiNO[[length(nomiNO)+1]]<-i
      }
    #nomiNO dovrebbe contenere i nomi delle "altre variabili" (come th in subset=x<th) 
    if(!is.null(nomiNO)) mfExt$formula<-update.formula(mfExt$formula,paste(".~.-", paste( nomiNO, collapse="-"), sep=""))

    mfExt<-eval(mfExt, parent.frame())

    #id.offs<-pmatch("offset",names(mf)) #questa identifica il nome offset(..). ELiminarlo dal dataframe? non conviene
    #       altrimenti nel model.frame non risulta l'offset

    #mantieni in mfExt solo le variabili che NON ci sono in mf (cosi la funzione occupa meno spazio..)
    #mfExt<-mfExt[,setdiff(names(mfExt), names(mf)),drop=FALSE]

    weights <- as.vector(model.weights(mf))
    offs <- as.vector(model.offset(mf))
    
#    if(!is.null(Call$weights)){ #"(weights)"%in%names(mf)
#      mfExt[all.vars(Call$weights, functions=FALSE)]<- eval(parse(text=all.vars(Call$weights)))
#      #names(mfExt)[which(names(mf)=="(weights)")]<- all.vars(Call$weights, functions=FALSE)#prima di 0.2.9-4 era as.character(Call$weights)
#      # mf["(weights)"]<-weights
#      }
#
#    if(!is.null(Call$offset)){ 
#     mfExt[all.vars(Call$offset, functions=FALSE)]<- eval(parse(text=all.vars(Call$offset)))
##     mfExt["(offset)"]<-eval(parse(text=all.vars(Call$offset)))
##      names(mf)[which(names(mf)=="(offset)")]<-all.vars(Call$offset, functions=FALSE) 
#     }

    mt <- attr(mf, "terms")
    interc<-attr(mt,"intercept")
    y <- model.response(mf, "any")
    XREG <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    
    #il cambio in mf da "offset(_nomevar_)" al "_nomevar_" deve avvenire dopo "model.matrix(mt, mf, contrasts)" 
#    if(!is.null(offs)){
#      #id.offs<-pmatch("offset",names(mf)) #questa identifica il nome offset(..). ELiminarlo dal dataframe? non conviene altrimenti nel model.frame non risulta l'offset
#      id.offs<- which(grepl("(offset)", names(mf))) #per consentire anche offset come argomento di glm()
#      names(mf)[id.offs]<- all.vars(formula(paste("~", names(mf)[id.offs])), functions=FALSE)
#      }
    
    namesXREG0<-colnames(XREG)
    #nameLeftSlopeZero<-setdiff(all.vars(seg.Z), all.vars(formula(obj)))
    nameLeftSlopeZero<-setdiff(all.vars(seg.Z), names(coef(obj))) #in questo modo riconosce che sin(x*pi) NON e' x, ad esempio.
    namesXREG0<-setdiff(namesXREG0, nameLeftSlopeZero)
    
    #dalla 0.3.0-1 eliminati i seguenti (tanto il modello viene stimato su mfExt)
    #nomeRispo<-strsplit(paste(formula(obj))[2],"/")[[1]] #eventuali doppi nomi separati da "/" (tipo "y/n" per GLM binom)
    #nomeRispo<-strsplit(paste(formula(obj))[2],"/")[[1]] #portato sopra
    #if(length(nomeRispo)>=2) mf[nomeRispo[1]]<-weights*y
    
    id.duplic<-match(all.vars(formula(obj)),all.vars(seg.Z),nomatch=0)>0
    if(any(id.duplic)) {
        #new.mf<-mf[,id.duplic,drop=FALSE]
        new.mf<-mf[,all.vars(formula(obj))[id.duplic],drop=FALSE]
        new.XREGseg<-data.matrix(new.mf)
        XREG<-cbind(XREG,new.XREGseg)
        }
    n.psi<- length(unlist(psi))
    id.n.Seg<-(ncol(XREG)-n.Seg+1):ncol(XREG)
    XREGseg<-XREG[,id.n.Seg,drop=FALSE]
    #XREG<-XREG[,-id.n.Seg,drop=FALSE]
    #XREG<-model.matrix(obj0) non va bene perche' non elimina gli eventuali mancanti in seg.Z..
    #Due soluzioni
    #XREG<-XREG[,colnames(model.matrix(obj)),drop=FALSE]
    #XREG<-XREG[,match(c("(Intercept)",all.vars(formula(obj))[-1]),colnames(XREG),nomatch =0),drop=FALSE]
    XREG <- XREG[, match(c("(Intercept)", namesXREG0),colnames(XREG), nomatch = 0), drop = FALSE]
    XREG<-XREG[,unique(colnames(XREG)), drop=FALSE]
    
    #################
    if(ncol(XREGseg)==1 && length(psi)==1 && n.psi==1 && !any(is.na(psi))) { if(psi==Inf) psi<-median(XREGseg)}
    #################

    n <- nrow(XREG)
    #Z <- list(); for (i in colnames(XREGseg)) Z[[length(Z) + 1]] <- XREGseg[, i]
    Z<-lapply(apply(XREGseg,2,list),unlist) #prende anche i nomi!
    name.Z <- names(Z) <- colnames(XREGseg)
    if(length(Z)==1 && is.vector(psi) && (is.numeric(psi)||is.na(psi))){
        psi <- list(as.numeric(psi))
        names(psi)<-name.Z
        }
    if (!is.list(Z) || !is.list(psi) || is.null(names(Z)) || is.null(names(psi))) stop("Z and psi have to be *named* list")
    id.nomiZpsi <- match(names(Z), names(psi))
    if ((length(Z)!=length(psi)) || any(is.na(id.nomiZpsi))) 
        stop("Length or names of Z and psi do not match")
    #dd <- match(names(Z), names(psi))
    nome <- names(psi)[id.nomiZpsi]
    psi <- psi[nome]
    initial.psi<-psi
    for(i in 1:length(psi)) {
        if(any(is.na(psi[[i]]))) psi[[i]]<-quantile(Z[[i]], prob= seq(0,1,l=K+2)[-c(1,K+2)], names=FALSE)
        }
    
    a <- sapply(psi, length)#b <- rep(1:length(a), times = a)
    id.psi.group <- rep(1:length(a), times = a) #identificativo di apparteneza alla variabile
    #Znew <- list()
    #for (i in 1:length(psi)) Znew[[length(Znew) + 1]] <- rep(Z[i], a[i])
    #Z <- matrix(unlist(Znew), nrow = n)
    Z<-matrix(unlist(mapply(function(x,y)rep(x,y),Z,a,SIMPLIFY = TRUE)),nrow=n)
    psi <- unlist(psi)
    psi<-unlist(tapply(psi,id.psi.group,sort))
    k <- ncol(Z)
    PSI <- matrix(rep(psi, rep(n, k)), ncol = k)
    colnames(Z) <- nomiZ <- rep(nome, times = a)
    ripetizioni <- as.numeric(unlist(sapply(table(nomiZ)[order(unique(nomiZ))], function(xxx) {1:xxx})))
    nomiU <- paste("U", ripetizioni, sep = "")
    nomiU <- paste(nomiU, nomiZ, sep = ".")
    nomiV <- paste("V", ripetizioni, sep = "")
    nomiV <- paste(nomiV, nomiZ, sep = ".")
    #forse non serve crearsi l'ambiente KK, usa mf..
    #obj <- update(obj, formula = Fo, data = mf)
    #if (model.frame) obj$model <- mf
  #controlla che model.frame() funzioni sull'oggetto restituito    
#    KK <- new.env()
#    for (i in 1:ncol(objframe$model)) assign(names(objframe$model[i]), objframe$model[[i]], envir = KK)
    if (it.max == 0) {
        #mf<-cbind(mf, mfExt)
        U <- pmax((Z - PSI), 0)
        colnames(U) <- paste(ripetizioni, nomiZ, sep = ".")
        nomiU <- paste("U", colnames(U), sep = "")
        #for (i in 1:ncol(U)) assign(nomiU[i], U[, i], envir = KK)
        #e' necessario il for? puoi usare colnames(U)<-nomiU;mf[nomiU]<-U
        for(i in 1:ncol(U)) mfExt[nomiU[i]]<-mf[nomiU[i]]<-U[,i]
        Fo <- update.formula(formula(obj), as.formula(paste(".~.+", paste(nomiU, collapse = "+"))))
        #obj <- update(obj, formula = Fo, data = KK)
        
        obj <- update(obj, formula = Fo, data = mfExt, evaluate=FALSE)
        if(!is.null(obj[["subset"]])) obj[["subset"]]<-NULL
        obj<-eval(obj, envir=mfExt)
        if (model) obj$model <-mf  #obj$model <- data.frame(as.list(KK))
        names(psi)<-paste(paste("psi", ripetizioni, sep = ""), nomiZ, sep=".")
        obj$psi <- psi
        return(obj)
    }
    #XREG <- model.matrix(obj) creata sopra         
    #o <- model.offset(objframe)
    #w <- model.weights(objframe)
    if (is.null(weights)) weights <- rep(1, n)
    if (is.null(offs)) offs <- rep(0, n)
    fam <- family(obj)
    initial <- psi
    obj0 <- obj
    dev0<-obj$dev
    list.obj <- list(obj)
#    psi.values <- NULL
    nomiOK<-nomiU
    opz<-list(toll=toll,h=h,stop.if.error=stop.if.error,dev0=dev0,visual=visual,it.max=it.max,nomiOK=nomiOK,
        fam=fam, eta0=obj$linear.predictors, maxit.glm=maxit.glm, id.psi.group=id.psi.group, gap=gap,
        pow=pow, visualBoot=visualBoot)   

    if(n.boot<=0){
      obj<-seg.glm.fit(y,XREG,Z,PSI,weights,offs,opz)
    } else {
      obj<-seg.glm.fit.boot(y, XREG, Z, PSI, weights, offs, opz, n.boot=n.boot, size.boot=size.boot, random=random) #jt, nonParam
      }
    if(!is.list(obj)){
        warning("No breakpoint estimated", call. = FALSE)
        return(obj0)
        }
    id.psi.group<-obj$id.psi.group
    nomiOK<-obj$nomiOK
    #nomiFINALI<-unique(sapply(strsplit(nomiOK, split="[.]"), function(x)x[2])) #nomi delle variabili con breakpoint stimati!
    nomiFINALI<-unique(sub("U[1-9]*[0-9].", "", nomiOK)) #nomi originali delle variabili con breakpoint stimati!
    #se e' stata usata una proc automatica "nomiFINALI" sara' differente da "name.Z"
    nomiSenzaPSI<-setdiff(name.Z,nomiFINALI)
    if(length(nomiSenzaPSI)>=1) warning("no breakpoints found for: ", paste(nomiSenzaPSI," "), call. = FALSE)
    it<-obj$it
    psi<-obj$psi
    k<-length(psi)
    psi.values<-if(n.boot<=0) obj$psi.values else obj$boot.restart
    U<-obj$U
    V<-obj$V
    #if(any(table(rowSums(V))<=1)) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close")
    for(jj in colnames(V)) {
        VV<-V[, which(colnames(V)==jj),drop=FALSE]
        sumV<-abs(rowSums(VV))
#        if( (any(diff(sumV)>=2) #se ci sono due breakpoints equivalenti
#            || any(table(sumV)<=1))) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close each other")
       if(any(table(sumV)<=1) && stop.if.error) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close each other")

        }
    rangeZ<-obj$rangeZ 
    obj<-obj$obj
    beta.c<-if(k == 1) coef(obj)["U"] else coef(obj)[paste("U", 1:ncol(U), sep = "")]
    psi.values[[length(psi.values) + 1]] <- psi
    id.warn <- FALSE
    if (n.boot<=0 && it > it.max) { #it >= (it.max+1)
        warning("max number of iterations attained", call. = FALSE)
        id.warn <- TRUE
    }
    Vxb <- V %*% diag(beta.c, ncol = length(beta.c))
    
    #se usi una procedura automatica devi cambiare ripetizioni, nomiU e nomiV, e quindi:
    length.psi<-tapply(as.numeric(as.character(names(psi))), as.numeric(as.character(names(psi))), length)
    forma.nomiU<-function(xx,yy)paste("U",1:xx, ".", yy, sep="")
    forma.nomiVxb<-function(xx,yy)paste("psi",1:xx, ".", yy, sep="")
    nomiU   <- unlist(mapply(forma.nomiU, length.psi, nomiFINALI)) #invece di un ciclo #paste("U",1:length.psi[i], ".", name.Z[i])
    nomiVxb <- unlist(mapply(forma.nomiVxb, length.psi, nomiFINALI))

    #mf<-cbind(mf, mfExt)
    for(i in 1:ncol(U)) {
        mfExt[nomiU[i]]<-mf[nomiU[i]]<-U[,i]
        mfExt[nomiVxb[i]]<-mf[nomiVxb[i]]<-Vxb[,i]
        }
#    for (i in 1:ncol(U)) {
#        assign(nomiU[i], U[, i], envir = KK)
#        assign(nomiVxb[i], Vxb[, i], envir = KK)
#    }
    nnomi <- c(nomiU, nomiVxb)
    Fo <- update.formula(formula(obj0), as.formula(paste(".~.+", 
        paste(nnomi, collapse = "+"))))
    #la seguente linea si potrebbe rimuovere perche' in mfExt c'e' gia' tutto..
    if(is.matrix(y)&& (fam$family=="binomial" || fam$family=="quasibinomial")){
              mfExt<-cbind(mfExt[[1]], mfExt[,-1])
    }
    objF <- update(obj0, formula = Fo, data = mfExt, evaluate=FALSE)
    if(!is.null(objF[["subset"]])) objF[["subset"]]<-NULL
    objF<-eval(objF, envir=mfExt)
    #C'e' un problema..controlla obj (ha due "(Intercepts)" - bhu.. al 27/03/14 non mi sembra!
    #Puo' capitare che psi sia ai margini e ci sono 1 o 2 osservazioni in qualche intervallo. Oppure ce ne 
    #   sono di piu' ma hanno gli stessi valori di x
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

#aggiornare qui i weights???? (piuttosto che sotto)
#------>>>
#------>>>
#------>>>
    if(!gap){
        names.coef<-names(objF$coefficients)
        if(k==1) {names(obj$coefficients)[match(c("U","V"), names(coef(obj)))]<- nnomi
          } else {
      names(obj$coefficients)[match(c(paste("U",1:k, sep=""), paste("V",1:k, sep="")), names(coef(obj)))]<- nnomi  
          }
        objF$coefficients[names.coef]<-obj$coefficients[names.coef]
#        objF$coefficients<- if(sum("(Intercept)"==names(obj$coef))==2) obj$coefficients[-2] else obj$coefficients
        objF$fitted.values<-obj$fitted.values
        objF$linear.predictors<-obj$linear.predictors
        objF$residuals<-obj$residuals
        objF$deviance<-obj$deviance
        objF$aic<-obj$aic + 2*k
        objF$weights<-obj$weights
        }
    Cov <- vcov(objF)
    id <- match(nomiVxb, names(coef(objF)))
    #cat(id,"\n")
    #return(objF)
    vv <- if (length(id) == 1) Cov[id, id] else diag(Cov[id, id])
    #if(length(initial)!=length(psi)) initial<-rep(NA,length(psi))
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
    class(objF) <- c("segmented", class(obj0))
    list.obj[[length(list.obj) + 1]] <- objF
    class(list.obj) <- "segmented"
    if (last) 
        list.obj <- list.obj[[length(list.obj)]]
    return(list.obj)
}
