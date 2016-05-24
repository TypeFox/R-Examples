seg.glm.fit<-function(y,XREG,Z,PSI,w,offs,opz,return.all.sol=FALSE){
#-------------------------
dpmax<-function(x,y,pow=1){
#deriv pmax
        if(pow==1) ifelse(x>y, -1, 0)
         else -pow*pmax(x-y,0)^(pow-1)
         }
#--------------------
    c1 <- apply((Z <= PSI), 2, all) #prima era solo <
    c2 <- apply((Z >= PSI), 2, all) #prima era solo >
    if(sum(c1 + c2) != 0 || is.na(sum(c1 + c2))) stop("psi out of the range")
    pow<-opz$pow
    eta0<-opz$eta0
    fam<-opz$fam
    maxit.glm<-opz$maxit.glm
    #--------------
    nomiOK<-opz$nomiOK
    toll<-opz$toll
    h<-opz$h
    stop.if.error<-opz$stop.if.error
    dev.new<-opz$dev0
    visual<-opz$visual
    id.psi.group<-opz$id.psi.group
    it.max<-old.it.max<-opz$it.max
    id.psi.group<-opz$id.psi.group
    gap<-opz$gap
    rangeZ <- apply(Z, 2, range)
    #k<-ncol(Z)
    psi<-PSI[1,]
    names(psi)<-id.psi.group
    it <- 1
    epsilon <- 10
    dev.values<-psi.values <- NULL
    id.psi.ok<-rep(TRUE, length(psi))
    sel.col.XREG<-unique(sapply(colnames(XREG), function(x)match(x,colnames(XREG))))
    if(is.numeric(sel.col.XREG))XREG<-XREG[,sel.col.XREG,drop=FALSE] #elimina le ripetizioni, ad es. le due intercette..
    while (abs(epsilon) > toll) {
        k<-ncol(Z)
        U <- pmax((Z - PSI), 0)^pow[1]
        V <- dpmax(Z,PSI,pow=pow[2])# ifelse((Z > PSI), -1, 0)
        X <- cbind(XREG, U, V)
        rownames(X) <- NULL
        if (ncol(V) == 1)
            colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c("U", "V")
        else colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c(paste("U",
            1:ncol(U), sep = ""), paste("V", 1:ncol(V), sep = ""))
        #obj <- lm.wfit(x = X, y = y, w = w, offset = o)
        #controlla******************
        obj <- suppressWarnings(glm.fit(x = X, y = y, offset = offs,
            weights = w, family = fam, control = glm.control(maxit = maxit.glm),
            etastart = eta0))
        eta0 <- obj$linear.predictors
        dev.old<-dev.new
        dev.new <- dev.new1<- obj$dev
        if(return.all.sol) dev.new1 <- glm.fit(x=cbind(XREG, U),y=y, family=fam, weights=w, offset=offs, etastart=eta0)$dev
        dev.values[[length(dev.values) + 1]] <- dev.new1
        if (visual) {
            flush.console()
            if (it == 1)
                cat(0, " ", formatC(dev.old, 3, format = "f"),
                  "", "(No breakpoint(s))", "\n")
            spp <- if (it < 10) "" else NULL
            cat(it, spp, "", formatC(dev.new1, 3, format = "f"), "",length(psi),"\n")
            #cat(paste("iter = ", it, spp," dev = ",formatC(dev.new,digits=3,format="f"), " n.psi = ",formatC(length(psi),digits=0,format="f"), sep=""), "\n")

        }
        epsilon <- (dev.new - dev.old)/(dev.old+.001)
#        epsilon <- (dev.new1 - dev.old)/dev.old #se vuoi usare la *vera* (e non la working che tiene conto dei gap) deviance
        obj$epsilon <- epsilon
        it <- it + 1
        obj$it <- it
        #class(obj) <- c("segmented", class(obj))
        #list.obj[[length(list.obj) + ifelse(last == TRUE, 0, 1)]] <- obj
        if (k == 1) {
            beta.c <- coef(obj)["U"]
            gamma.c <- coef(obj)["V"]
        }
        else {
            beta.c <- coef(obj)[paste("U", 1:ncol(U), sep = "")]
            gamma.c <- coef(obj)[paste("V", 1:ncol(V), sep = "")]
        }
        if (it > it.max) break
        psi.values[[length(psi.values) + 1]] <- psi.old <- psi
        #if(it>=old.it.max && h<1) H<-h
        psi <- psi.old + h*gamma.c/beta.c
        PSI <- matrix(rep(psi, rep(nrow(Z), ncol(Z))), ncol = ncol(Z))
        #check if psi is admissible..
        a <- apply((Z <= PSI), 2, all) #prima era solo <
        b <- apply((Z >= PSI), 2, all) #prima era solo >
        if(stop.if.error) {
        isErr <- sum(a + b) != 0 || is.na(sum(a + b))
            if(isErr) {
                if(return.all.sol) return(list(dev.values, psi.values)) else stop("(Some) estimated psi out of its range")
                }
            } else {
            id.psi.ok<-!is.na((a+b)<=0)&(a+b)<=0
            Z <- Z[,id.psi.ok,drop=FALSE]
            psi <- psi[id.psi.ok]
            PSI <- PSI[,id.psi.ok,drop=FALSE]
            nomiOK<-nomiOK[id.psi.ok] #salva i nomi delle U per i psi ammissibili
            id.psi.group<-id.psi.group[id.psi.ok]
            names(psi)<-id.psi.group
            if(ncol(PSI)<=0) return(0)
            } #end else
        #obj$psi <- psi
    } #end while
    #queste due righe aggiunte nella versione 0.2.9-3 (adesso i breakpoints sono sempre ordinati)
    psi<-unlist(tapply(psi, id.psi.group, sort))
    names(psi)<-id.psi.group
    PSI <- matrix(rep(psi, rep(nrow(Z), length(psi))), ncol = length(psi))
    #aggiunto da qua..
    U <- pmax((Z - PSI), 0)
    V <- ifelse((Z > PSI), -1, 0)
    X <- cbind(XREG, U, V)
    rownames(X) <- NULL
    if (ncol(V) == 1) colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c("U", "V")
        else colnames(X)[(ncol(XREG) + 1):ncol(X)] <- c(paste("U", 1:ncol(U), sep = ""), paste("V", 1:ncol(V), sep = ""))
    #tolto il suppressWarnings(
    obj <- glm.fit(x = X, y = y, offset = offs,
            weights = w, family = fam, control = glm.control(maxit = maxit.glm),
            etastart = eta0)
    obj$epsilon <- epsilon
    obj$it <- it
    obj.new <- glm.fit(x = cbind(XREG, U), y = y, offset = offs,
            weights = w, family = fam, control = glm.control(maxit = maxit.glm),
            etastart = eta0)
    SS.new<- obj.new$dev
    if(!gap){
          names.coef<-names(obj$coefficients)
          obj$coefficients<-c(obj.new$coefficients, rep(0,ncol(V)))
          names(obj$coefficients)<-names.coef
          obj$residuals<-obj.new$residuals
          obj$fitted.values<-obj.new$fitted.values
          obj$linear.predictors<-obj.new$linear.predictors
          obj$deviance<-obj.new$deviance
          obj$weights<-obj.new$weights
          obj$aic<-obj.new$aic #+ 2*ncol(V) #ho fatto la modifica in segmented.glm(): "objF$aic<-obj$aic + 2*k"
          }
    #fino a qua..
    obj<-list(obj=obj,it=it,psi=psi,psi.values=psi.values,U=U,V=V,rangeZ=rangeZ,
        epsilon=epsilon,nomiOK=nomiOK, dev.no.gap=SS.new, id.psi.group=id.psi.group)
    return(obj)
    }
