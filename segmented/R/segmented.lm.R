`segmented.lm` <-
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
        if(visual) {visual<-FALSE; visualBoot<-TRUE}# warning("`display' set to FALSE with bootstrap restart", call.=FALSE)}
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
    orig.call$formula<- mf$formula<-formula(obj) #per consentire lm(y~.)
    m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    if(class(mf$formula)=="name" && !"~"%in%paste(mf$formula)) mf$formula<-eval(mf$formula)
    #orig.call$formula<-update.formula(orig.call$formula, paste("~.-",all.vars(seg.Z)))
#
#genn 2013. dalla versione 0.2.9-4 ho tolto if(length(.. Tra l'altro non capisco perche' lo avevo fatto
#    if(length(all.vars(formula(obj)))>1){
#    mf$formula<-update.formula(mf$formula,paste(paste(seg.Z,collapse=".+"),"+",paste(all.vars(formula(obj))[-1],collapse="+")))
#    } else {
#    mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
#    }
#nov 2013 dalla versione 0.3-0.0 (che dovrebbe essere successiva alla 0.2-9.5) viene creato anche il modelframe esteso che comprende
# termini "originali", prima che fossero trasformati (Ad es., x prima che ns(x) costruisca le basi). Questo permette di avere termini
# ns(), poly(), bs() nel modello di partenza
    mfExt<- mf
    mf$formula<-update.formula(mf$formula,paste(seg.Z,collapse=".+"))
    #mfExt$formula<- update.formula(mfExt$formula,paste(paste(seg.Z,collapse=".+"),"+",paste(all.vars(formula(obj)),collapse="+")))   
#    mfExt$formula<- if(!is.null(obj$call$data)) 
#        update.formula(mf$formula,paste(".~",paste(all.vars(obj$call), collapse="+"),"-",obj$call$data,sep=""))
#            else update.formula(mf$formula,paste(".~",paste(all.vars(obj$call), collapse="+"),sep=""))
#-----------
#   browser()
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
    #questo serve per inserire in mfExt le eventuali variabili contenute nella formula con offset(..)
    nomiOff<-setdiff(all.vars(formula(obj)), names(mf))
    if(length(nomiOff)>=1) mfExt$formula<-update.formula(mfExt$formula,paste(".~.+", paste( nomiOff, collapse="+"), sep=""))
#----------------------------------------------------
#    browser()

#ago 2014 c'e' la questione di variabili aggiuntive...
nomiTUTTI<-all.vars(mfExt$formula) #comprende anche altri nomi (ad es., threshold) "variabili"
nomiNO<-NULL 
for(i in nomiTUTTI){
    r<-try(eval(parse(text=i), parent.frame()), silent=TRUE)
    if(class(r)!="try-error" && length(r)==1 && !is.function(r)) nomiNO[[length(nomiNO)+1]]<-i
    }
#nomiNO dovrebbe contenere i nomi delle "altre variabili" (come th in subset=x<th) 
if(!is.null(nomiNO)) mfExt$formula<-update.formula(mfExt$formula,paste(".~.-", paste( nomiNO, collapse="-"), sep=""))
#funziona ma se c'e' qualche variabile con il nome di una funzione (ad es., "filter")
# anche questa variabile va a finire in nomiNO e quindi poi viene eliminata
# problema risolto con "r<-try(evalq(.."? (prima era "r<-try(eval(..")
#----------------------------------------------------   
    mfExt<-eval(mfExt, parent.frame())
    
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
#      names(mf)[which(names(mf)=="(offset)")]<-all.vars(Call$offset, functions=FALSE) 
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
        if(any(is.na(psi[[i]]))) psi[[i]]<-if(control$quant) {quantile(Z[[i]], prob= seq(0,1,l=K+2)[-c(1,K+2)], names=FALSE)} else {(min(Z[[i]])+ diff(range(Z[[i]]))*(1:K)/(K+1))}
        }

    a <- sapply(psi, length)

    #per evitare che durante il processo iterativo i psi non siano ordinati
    id.psi.group <- rep(1:length(a), times = a) #identificativo di apparteneza alla variabile
    #

    #Znew <- list()
    #for (i in 1:length(psi)) Znew[[length(Znew) + 1]] <- rep(Z[i], a[i])
    #Z <- matrix(unlist(Znew), nrow = n)
    Z<-matrix(unlist(mapply(function(x,y)rep(x,y),Z,a,SIMPLIFY = TRUE)),nrow=n)
    psi <- unlist(psi)
    #se psi e' numerico, la seguente linea restituisce i valori ordinati all'interno della variabile..
    psi<-unlist(tapply(psi,id.psi.group,sort))
    k <- ncol(Z)
    PSI <- matrix(rep(psi, rep(n, k)), ncol = k)
    #controllo se psi e' ammissibile..
    c1 <- apply((Z <= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo <)
    c2 <- apply((Z >= PSI), 2, all) #dovrebbero essere tutti FALSE (prima era solo >)
    if(sum(c1 + c2) != 0 || is.na(sum(c1 + c2)) ) stop("starting psi out of the admissible range")
    
    #questo dovrebbe eliminare i psi non-ammissib. ma non sono sicuro cosa succede se ci sono piu' variabili
#    if(sum(c1 + c2) != 0){
#     id.val.psi<-!((c1+c1)>0) #individua i psi ammissibili (i.e. interni)
#     psi<-psi[id.val.psi]
#     Z<-Z[,id.val.psi]
#     PSI<-PSI[,id.val.psi]
#    } 
#    if(is.na(sum(c1 + c2))) stop("psi out of the range")

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
        obj <- update(obj, formula = Fo, evaluate=FALSE, data=mfExt) #data = mf, 
        if(!is.null(obj[["subset"]])) obj[["subset"]]<-NULL
        obj<-eval(obj, envir=mfExt)
        if (model) obj$model <-mf  #obj$model <- data.frame(as.list(KK))

        psi <- cbind(psi, psi, 0)
        rownames(psi) <- paste(paste("psi", ripetizioni, sep = ""), nomiZ, sep=".")
        colnames(psi) <- c("Initial", "Est.", "St.Err")


        #names(psi)<-paste(paste("psi", ripetizioni, sep = ""), nomiZ, sep=".")
        obj$psi <- psi
        return(obj)
    }
    #XREG <- model.matrix(obj) creata sopra
    #o <- model.offset(objframe)
    #w <- model.weights(objframe)
    if (is.null(weights)) weights <- rep(1, n)
    if (is.null(offs)) offs <- rep(0, n)
    initial <- psi
    obj0 <- obj
    dev0<-sum(obj$residuals^2)
    list.obj <- list(obj)
#    psi.values <- NULL
    nomiOK<-nomiU
    opz<-list(toll=toll,h=h,stop.if.error=stop.if.error,dev0=dev0,visual=visual,it.max=it.max,
        nomiOK=nomiOK,id.psi.group=id.psi.group,gap=gap,visualBoot=visualBoot,pow=pow)
    if(n.boot<=0){
    obj<-seg.lm.fit(y,XREG,Z,PSI,weights,offs,opz)
    } else {
    obj<-seg.lm.fit.boot(y, XREG, Z, PSI, weights, offs, opz, n.boot=n.boot, size.boot=size.boot, random=random) #jt, nonParam
      }
    if(!is.list(obj)){
        warning("No breakpoint estimated", call. = FALSE)
        return(obj0)
        }
    if(obj$obj$df.residual==0) warning("no residual degrees of freedom (other warnings expected)", call.=FALSE)
    id.psi.group<-obj$id.psi.group
    nomiOK<-obj$nomiOK
    #nomiFINALI<-unique(sapply(strsplit(nomiOK, split="[.]"), function(x)x[2])) #nomi delle variabili con breakpoint stimati!
    #nomiFINALI<-sub("U[1-9].", "", nomiOK) #nomi originali delle variabili con breakpoint stimati!
    nomiFINALI<- unique(sub("U[1-9]*[0-9].", "", nomiOK))
    #se e' stata usata una proc automatica "nomiFINALI" sara' differente da "name.Z"
    nomiSenzaPSI<-setdiff(name.Z,nomiFINALI)
    if(length(nomiSenzaPSI)>=1) warning("no breakpoints found for: ", paste(nomiSenzaPSI," "), call. = FALSE)
    it<-obj$it
    psi<-obj$psi
    psi.values<-if(n.boot<=0) obj$psi.values else obj$boot.restart
    U<-obj$U
    V<-obj$V
#browser()
    
    #if(any(table(rowSums(V))<=1)) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close")
    for(jj in colnames(V)) {
        VV<-V[, which(colnames(V)==jj), drop=FALSE]
        sumV<-abs(rowSums(VV))
#        if( (any(diff(sumV)>=2) #se ci sono due breakpoints uguali
#            || any(table(sumV)<=1)) && stop.if.error) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close each other")
# Tolto perche' se la variabile segmented non e' ordinata non ha senso..
#magari potresti fare un abs(diff(psi))<=.0001? ma clusterizzato..
        if(any(table(sumV)<=1) && stop.if.error) stop("only 1 datum in an interval: breakpoint(s) at the boundary or too close each other")

        }
    rangeZ<-obj$rangeZ
    obj<-obj$obj
    k<-length(psi)
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
    
    #nomiU   <- unlist(mapply(forma.nomiU, length.psi, name.Z)) #invece di un ciclo #paste("U",1:length.psi[i], ".", name.Z[i])
    #nomiVxb <- unlist(mapply(forma.nomiVxb, length.psi, name.Z))
    nomiU   <- unlist(mapply(forma.nomiU, length.psi, nomiFINALI)) #invece di un ciclo #paste("U",1:length.psi[i], ".", name.Z[i])
    nomiVxb <- unlist(mapply(forma.nomiVxb, length.psi, nomiFINALI))

    
#    nomiVxb <-paste("psi",ripetizioni, ".",nomiZ ,sep="")
#Dalla 0.2.9-5 eliminati i seguenti. La linea sopra sembra sufficiente
#    colnames(U) <- colnames(Vxb) <-sapply(strsplit(nomiOK,"U"),function(x)x[2])
#    #colnames(U) <- paste(ripetizioni, nomiZ, sep = ".")
#    #colnames(Vxb) <- paste(ripetizioni, nomiZ, sep = ".")
#    nomiU <- paste("U", colnames(U), sep = "")
#    nomiVxb <- paste("psi", colnames(Vxb), sep = "")
#    for (i in 1:ncol(U)) {
#        assign(nomiU[i], U[, i], envir = KK)
#        assign(nomiVxb[i], Vxb[, i], envir = KK)
#    }

    #mf<-cbind(mf, mfExt) #questo creava ripetizioni..
    for(i in 1:ncol(U)) {
        mfExt[nomiU[i]]<-mf[nomiU[i]]<-U[,i]
        mfExt[nomiVxb[i]]<-mf[nomiVxb[i]]<-Vxb[,i]
        }
    nnomi <- c(nomiU, nomiVxb)
#    browser()
    Fo <- update.formula(formula(obj0), as.formula(paste(".~.+", paste(nnomi, collapse = "+"))))
    #objF <- update(obj0, formula = Fo, data = KK)
    objF <- update(obj0, formula = Fo,  evaluate=FALSE, data = mfExt)
    #eliminiamo subset, perche' se e' del tipo subset=x>min(x) allora continuerebbe a togliere 1 osservazione 
    if(!is.null(objF[["subset"]])) objF[["subset"]]<-NULL
    objF<-eval(objF, envir=mfExt)
    #Puo' capitare che psi sia ai margini o molto vicini (e ci sono solo 1 o 2 osservazioni in qualche intervallo. Oppure ce ne 
    #sono di piu' ma hanno gli stessi valori di x. In questo caso objF$coef puo' avere mancanti.. names(which(is.na(coef(objF))))
    
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
    
    if(!gap){
      names.coef<-names(objF$coefficients)
      #questi codici funzionano e si basano sull'assunzioni che le U e le V siano ordinate..
      if(k==1){
        names(obj$coefficients)[match(c("U","V"), names(coef(obj)))]<- nnomi
          } else {
        names(obj$coefficients)[match(c(paste("U",1:k, sep=""), paste("V",1:k, sep="")), names(coef(obj)))]<- nnomi  
          }
      objF$coefficients[names.coef]<-obj$coefficients[names.coef] #sostituisce gli 0 
      #objF$coefficients<-obj$coefficients
      #names(objF$coefficients)<-names.coef
      objF$fitted.values<-obj$fitted.values
      objF$residuals<-obj$residuals
    }
    Cov <- vcov(objF) 
    id <- match(nomiVxb, names(coef(objF)))
    vv <- if (length(id) == 1) Cov[id, id] else diag(Cov[id, id])
    #if(length(initial)!=length(psi)) initial<-rep(NA,length(psi))
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
    #initial<-unlist(mapply(function(x,y){if(is.na(x)[1])rep(x,y) else x }, initial.psi, a.ok, SIMPLIFY = TRUE))
    initial<-unlist(mapply(function(x,y){if(is.na(x)[1])rep(x,y) else x }, initial.psi[nomiFINALI], a.ok[a.ok!=0], SIMPLIFY = TRUE))
    ris.psi[,1]<-initial
    #psi <- cbind(initial, psi, sqrt(vv))
    #rownames(psi) <- colnames(Cov)[id]
    
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
