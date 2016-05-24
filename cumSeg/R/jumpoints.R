jumpoints <-
function(y, x, k=min(30,round(length(y)/10)), output="2",
    psi=NULL, round=TRUE, control = fit.control(), selection=sel.control(), ...)  {
#jump-point models
#y: the response; x the explanatory (if missing index integers are assumed)
#psi: the starting values for the breakpoints. If NULL k quantiles are used.
#k: the number of breakpoints to be estimated. This argument is ignored when psi is specified
#output: "1" restituisce tutti i breakpoint individuati. fit.control() può essere usato per modificare
#       it.max. Aumentandolo può ridurre il n. di "putative" (candidate) psi e migliorare la performance
#       della selezione dei psi attraverso un qualche criterio.
#output: "2" applica il criterio in selection per selezionare i jumpoints "significativi"
#output: "3" ri-applica l'algoritmo segmented assumendo come punti di partenza quelli selezionati; in genere
#       se ci sono psi spuri questi possono essere ulteriormente eliminati
#------------
#--- seg.lm.fit0
#funzioni interne
seg.lm.fit0<- function(y, Z, PSI, control, round=FALSE, ...){
#Questa è una versione semplificata di seg.lm.fit() per evitare calcoli inutili
#y: la risposta
#Z: matrice di variabile segmented
#control: lista che controlla il processo di stima
#round: approssimare le soluzioni del punto di svolta?
#----------------
      it.max <- old.it.max <- control$it.max
      toll <- control$toll
      visual <- control$visual
      last <- control$last
      stop.if.error<-control$stop.if.error
      h <- min(abs(control$h), 1)
      if (h < 1)
        it.max <- it.max + round(it.max/2)
      it <- 1
      epsilon<-10
      k <- ncol(PSI)
      psi.values <- NULL
      H <- 1
      psi<-PSI[1,]
      #NB Poichè Z contiene ripetizioni della stessa variabile è sufficiente prendere Z[,1]
      #if (intercept)  XREG <- cbind(1,Z[,1],Xlinear) else XREG <- cbind(Z[,1],Xlinear)
      #se AR
#      n<-length(y)
#      XREG<-cbind(c(y[1],y[1:(n-1)]),Z[,1])
      #
      XREG<-cbind(Z[,1])
      #obj sotto serve solo per la stampare la dev
      #obj<-lm.fit(y=y, x=XREG)
      obj<-list(residuals=rep(10,3))
      while (abs(epsilon) > toll) {
        U <- pmax((Z - PSI), 0)
        V <- ifelse((Z > PSI), -1, 0)
        X<-cbind(XREG, U, V)
        dev.old <- sum(obj$residuals^2)
        rownames(X) <- NULL
        if (ncol(V) == 1) {
            #colnames(X)[ (ncol(XREG) + 1):ncol(X)] <- c("U", "V")
            colnames(X)[ (ncol(XREG) ):ncol(X)] <- c("firstSlope","U", "V")
              } else {
            colnames(X)<-rev(c(paste("V", ncol(V):1, sep = ""),
                paste("U",ncol(U):1, sep = ""),rep("firstSlope",ncol(X)-ncol(U)-ncol(V))))
            }
        obj <- lm.fit(x = X, y = y) #drop(solve(crossprod(X),crossprod(X,y)))
        dev.new <- sum(obj$residuals^2)
        if (visual) {
          if (it == 1) cat(0, " ", formatC(dev.old, 3, format = "f"),"", "(No breakpoint(s))", "\n")
            spp <- if (it < 10) "" else NULL
            cat(it, spp, "", formatC(dev.new, 3, format = "f"), "---",ncol(V),"breakpoints","\n")
        }
        epsilon <- (dev.new - dev.old)/dev.old
        obj$epsilon <- epsilon
        it <- it + 1
        obj$it <- it
        class(obj) <- c("segmented", class(obj))
        if (k == 1) {
            beta.c <- coef(obj)["U"]
            gamma.c <- coef(obj)["V"]
        } else {
          #se ci sono contrasti i beta.c quali sono?
          beta.c <- coef(obj)[paste("U", 1:ncol(U), sep = "")]
          gamma.c <- coef(obj)[paste("V", 1:ncol(V), sep = "")]
        }
        if (it > it.max) break
        psi.values[[length(psi.values) + 1]] <- psi.old <- psi
        if (it >= old.it.max && h < 1) H <- h
        psi <- round(psi.old + H * gamma.c/beta.c,0)
        PSI <- matrix(rep(psi, rep(nrow(Z), ncol(Z))), ncol = ncol(Z))
        #check if psi is admissible..
        a <- apply((Z <= PSI), 2, all)
        b <- apply((Z >= PSI), 2, all)
        if(stop.if.error) {
          if(sum(a + b) != 0 || is.na(sum(a + b))) stop("(Some) estimated psi out of its range")
          } else {
          id.psi.ok<-!is.na((a+b)<=0)&(a+b)<=0
          Z <- Z[,id.psi.ok,drop=FALSE]
          psi <- psi[id.psi.ok]
          PSI <- PSI[,id.psi.ok,drop=FALSE]
          #id.psi.ok<-!a|b #indici di psi validi
#          Z <- Z[,!is.na(id.psi.ok),drop=FALSE]
#          psi <- psi[!is.na(id.psi.ok)]
#          PSI <- PSI[,!is.na(id.psi.ok),drop=FALSE]
          }
    if(ncol(PSI)<=0) {
         warning("No breakpoint estimated", call. = FALSE)
         obj<-lm.fit(x = XREG, y = y)
         obj$fitted.values<-rep(obj$coef,length(y))
         obj$est.means<-obj$coef
         return(obj)
         }
    } # end while
   if(round) {
      psi<-round(psi,0)
      PSI <- matrix(rep(psi, rep(nrow(Z), ncol(Z))), ncol = ncol(Z))
      V <- ifelse((Z > PSI), -1, 0)
      #V serve per i fitted...si può evitare di crearla?
      }
   #obj$psi <- if(round) round(sort(psi),2) else sort(psi)
   obj$psi <- sort(psi)
   obj$beta.c <- beta.c[order(psi)]
   obj$gamma.c <- gamma.c[order(psi)]
   obj$epsilon <- epsilon
   obj$V<- V[,order(psi)]
#un'ultima verifica..
    obj$psi<-obj$psi[!is.na(obj$beta.c)]
    obj$V<-as.matrix(as.matrix(obj$V)[,!is.na(obj$beta.c)])
    obj$beta.c<-obj$beta.c[!is.na(obj$beta.c)]
   return(obj)
  }
#
#-- pen.MDL
pen.MDL<-function(id,n){
#restituisce un vettore (di dim=length(id)) che rappresenta la penalità 2*sum\log n_j per
#ogni partizione (active set)
#length(id) è il num (max) di breakpoints ed n è il vettore delle numerosità della
#partizione.
    do.m<-function(id,n.col){
        blockdiag <- function(...) {
            args <- list(...)
            nc <- sapply(args,ncol)
            cumnc <- cumsum(nc)
            ##  nr <- sapply(args,nrow)
            ## NR <- sum(nr)
            NC <- sum(nc)
            rowfun <- function(m,zbefore,zafter) {
              cbind(matrix(0,ncol=zbefore,nrow=nrow(m)),m,
                    matrix(0,ncol=zafter,nrow=nrow(m)))
            }
            ret <- rowfun(args[[1]],0,NC-ncol(args[[1]]))
            for (i in 2:length(args)) {
              ret <- rbind(ret,rowfun(args[[i]],cumnc[i-1],NC-cumnc[i]))
            }
            ret
            } #end blockgiag
        id<-sort(id) #sort(unlist(id))
        if(length(id)==1) {
            m<-t(rep(1,id))} else {
        m<-do.call(blockdiag,lapply(c(id[1],diff(id)),function(xx)t(rep(1,xx))))
          }
        m<-blockdiag(m,t(rep(1,n.col-ncol(m))))
        m
        } #end do.m
      #inizio codici veri
      if(length(n)!=(length(id)+1)) stop("Errore in 'id' o 'n'")
      A<-matrix(rev(id),length(id),length(id),byrow=FALSE)
      A[col(A)>row(A)]<-NA
      r<-rev(apply(A,2,function(x)(x[!is.na(x)])))
      lista.m<-lapply(r,do.m,n.col=length(n))
      #sapply(lista.m,function(xx)drop(xx%*%n))
      ris<-sapply(lista.m,function(xx)2*sum(log(drop(xx%*%n))))
      ris<-c(2*sum(log(n)),ris)
      ris
      }
#-----------
#------------
     n <- length(y)
#     if(any(is.na(y))) {
#       id<-complete.cases(y,x)
#        x<-x[id]
#        y<-y[id]
#        }
     if(missing(x)) {
            x<-1:n
            Y<-cumsum(y)
          } else {
            if(length(x)!=n) stop("Lengths of x and y differ")
            y<-y[order(x)]
            x<-sort(x)
            diffx <- c(x[1],diff(x))
#            DD <- matrix(diffx, ncol=n, nrow=n, byrow=TRUE)
#            DD[col(DD)>row(DD)]<-0   # matrice di trasformazione di y in Y
#            Y <- drop(DD%*%y)
            Y<-cumsum(y*diffx)
          }
     if(is.null(psi)) psi<-quantile(x, prob= seq(0,1,l=k+2)[-c(1,k+2)], names=FALSE)
      k<-length(psi)
      Z <- matrix(rep(x, k), nrow = n)
      PSI <- matrix(rep(psi, rep(n, k)), ncol = k)

    it.max <- old.it.max <- control$it.max
    if (it.max == 0)  U <- pmax((Z - PSI), 0)

    Xlin<- NULL
    obj<-seg.lm.fit0(y=Y, Z=Z, PSI=PSI, control=control, round=round, ...)
    obj$y<-y

    if(!is.null(obj$psi)){
        obj$fitted.values<-drop(abs(obj$V)%*%obj$beta.c)
        if("firstSlope"%in%names(coef(obj))) obj$fitted.values<- obj$fitted.values + obj$coef["firstSlope"]
        obj$id.group<- -rowSums(obj$V)
    }
    obj$est.means<-cumsum(c(obj$coef["firstSlope"], obj$beta.c))
    obj$n.psi<-length(obj$psi)
##---------------------------------------
##--- primo output
    if(output=="1") {
        class(obj)<-"aCGHsegmented"
        return(obj)
        }
    #se vuoi selezionare anche i psi significativi...
    psi0<-obj$psi
    est.means0<-obj$est.means
    display1<-selection$display
    edf.psi<-selection$edf.psi
    type<-selection$type
    ##nome<-deparse(substitute(type))
    #plot.it<-selection$plot.it
    Cn<-eval(parse(text=selection$Cn))
    S<-selection$S
    tipoAlg<-selection$alg
    require(lars)
    if(is.null(obj$psi)) stop("No estimated breakpoint in obj")
    #-------------------
    #pesi2<-1/abs(diff(c(obj$psi,n))*obj$beta.c)
    #pesi3<-1/sqrt(diff(c(obj$psi,n))^2+obj$beta.c^2)
    #s1<-scale(diff(c(obj$psi,n)))
    #b1<-scale(obj$beta.c)
    #pesi3<-1/sqrt(b1^2+(s1/2)^2)
    #olars<-lars(t(pesi3*t(abs(obj$V))), y=y, type="lasso", normalize=FALSE, intercept=TRUE, trace=display1)
    #
    ##se bic mdl non considerare penaliz.:
    #if(type!="rss") tipoAlg<-"stepwise"
    olars<-lars(abs(obj$V), y=y, type=tipoAlg, normalize=FALSE, intercept=TRUE, trace=display1)
    #type="stepwise" ,"lasso"
    id.var.entry<-(1:ncol(obj$V))[order(olars$entry)]
    edf<- if(edf.psi) (olars$df-1)*2+1 else olars$df
    RSS<-olars$RSS

    #######ottieni RSS unpenalized
#    if(RSS.unp){
#    RSS<-vector(length=length(RSS))
#    RSS[1]<-olars$RSS[1]
#    for(i in 1:(length(RSS)-1)){
#        X.ok<- cbind(1,obj$V[,id.var.entry[1:i]])
#        RSS[i+1]<-sum(lm.fit(y=y,x=X.ok)$residuals^2)
#        }
#    }
    #######################################################
    max.rss<-function(RSS){
          var.mod<-(RSS/n)
          ll<- -(log(2*pi*var.mod)+1)*n/2
          new.r<-((ll[length(ll)]-ll[-1])/(ll[length(ll)]-ll[2]))*(length(ll)-1)+1
          diff2<-diff(new.r,diff=2)>S
          if(!any(diff2)) return(0)
          maxll=max(which(diff2))+1
          return(maxll)
          }
    min.r<-switch(type,
        bic = which.min(log(RSS/n)+log(n)*edf*Cn/n),
        mdl = which.min(n*log(RSS/(n-edf))+Cn*pen.MDL(id.var.entry,as.numeric(table(-rowSums(obj$V))))),
        rss = max.rss(RSS)
        )
    crit<-switch(type,
        bic = (log(RSS/n)+log(n)*edf*Cn/n),
        mdl = (n*log(RSS/(n-edf))+Cn*pen.MDL(id.var.entry,as.numeric(table(-rowSums(obj$V))))),
        rss = (RSS)
        )
#        rss = {
#            diff.rss<-diff(RSS,diff=2)>.75
#            if(any(diff.rss)) max(which(diff.rss)) else 0 #olars$RSS
#              }
        #new.r<-((RSS[length(RSS)]-RSS[-1])/(RSS[length(RSS)]-rss[2]))*(length(RSS)-1)+1
        #max(which(c(10,diff(new.r,diff=2))>.75))
    id<-sort(c(0,id.var.entry)[1:min.r])
    if(length(id)<=1) {
        o<-list(y=y,est.means=mean(y),id.var.entry=id.var.entry,n.psi=0,criterion=crit)
        class(o)<-"aCGHsegmented"
        return(o)
        }
    id<-id[-1] #escludi l'intercetta
    psi1<- obj$psi[id]
    if(output=="2"){
        V.ok<- cbind(1,abs(obj$V[,id]))
        if(tipoAlg!="stepwise"){
            hat.b<-drop(solve(crossprod(V.ok),crossprod(V.ok,y))) #better than hat.b<-lm.fit(x=V.ok,y=y)$coef
            } else {
              hat.b<-olars$beta[min.r,id]
              hat.b<-c(olars$mu-sum(colMeans(as.matrix(V.ok[,-1]))*hat.b),hat.b)
            }
        fittedvalues<-drop(V.ok%*%hat.b)
        #gruppi<-rowSums(V.ok[,-1])
        #fittedvalues<-rep(cumsum(hat.b),tapply(gruppi,gruppi,length))
        ris<-list(id.var.entry=id.var.entry,psi=psi1,n.psi=length(id), psi.order.entry=psi0[id.var.entry],
                  psi0=psi0,est.means0=est.means0,est.means=cumsum(hat.b),criterion=crit,
                  fitted.values=fittedvalues)
        ris$y<-y
        class(ris)<-"aCGHsegmented"
        return(ris)
        }
    k<-length(psi1)
    Z <- matrix(rep(x, k), nrow = n)
    PSI <- matrix(rep(psi1, rep(n, k)), ncol = k)
    obj<-seg.lm.fit0(y=Y, Z=Z, PSI=PSI, control=fit.control(toll = 1e-04, it.max = 10, stop.if.error=FALSE),
        round=round, ...)
#    obj$y<-y
    if(!is.null(obj$est.means)){
        obj$n.psi<-0
        obj$psi0<-psi0
        obj$psi1<-psi1
        obj$criterion<-crit
        return(obj)
      }
    fitted.v<-drop(abs(obj$V)%*%obj$beta.c)
    if("firstSlope"%in%names(coef(obj))) fitted.v<- fitted.v + obj$coef["firstSlope"]
    obj$id.group<- -rowSums(obj$V)
    est.means<-cumsum(c(obj$coef["firstSlope"], obj$beta.c)) #unique(fitted.v)
    ris<-list(id.var.entry=id.var.entry,psi=obj$psi,n.psi=length(obj$psi), psi.order.entry=psi0[id.var.entry],
                  psi0=psi0,psi1=psi1,est.means0=est.means0, est.means=est.means,criterion=crit) #type=r)

    ris$fitted.values<-fitted.v
    ris$y<-y
    class(ris)<-"aCGHsegmented"
    return(ris)
    }

