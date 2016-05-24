`tempeff` <-
function(formula, data, subset, na.action, fcontrol=fit.control(), etastart=NULL, drop.L, ...){
# source("d:/lavori/jss/modtempeff/r/tempeff.r")
#--required functions: bspline() lagged()
bspline<-function(x, ndx, xlr=NULL, deg=3, deriv=0){
#x: vettore di dati
#xlr: il vettore di c(xl,xr)
#ndx: n.intervalli in cui dividere il range
#deg: il grado della spline
#Restituisce ndx+deg basis functions per ndx-1 inner nodi
#require(splines)
    if(is.null(xlr)){
    xl<-min(x)-.01*diff(range(x))
    xr<-max(x)+.01*diff(range(x))
       } else {
         if(length(xlr)!=2) stop("quando fornito, xlr deve avere due componenti")
        xl<-xlr[1]
        xr<-xlr[2]
     }
    dx<-(xr-xl)/ndx
    knots<-seq(xl-deg*dx,xr+deg*dx,by=dx)
    B<-splineDesign(knots,x,ord=deg+1,derivs=rep(deriv,length(x)))
    #B<-spline.des(knots,x,bdeg+1,0*x) #$design
    B #the B-spline base matrix
}#end_fn

lagged<-function(x,lag=1){#by T.Lumley
        if (lag==0) return(x)
            n<-length(x)
            c(rep(NA,lag),x[-( (n-lag+1):n)])
      }
#-----------------------------------
    if (missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m<-match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame") #restituisce model.frame(formula = y ~ x, data = d, drop.unused.levels = TRUE)
    mf <- eval(mf, parent.frame()) #restituisce una dataframe
    y <- model.response(mf, "any")
    mt <- attr(mf, "terms")
    X<-model.matrix(mt, mf, contrasts) 
#    if(!"(Intercept)"%in%colnames(X)) stop("The model intercept is required")
    tf <- terms(formula, specials = c("csdl","seas","dl"))
    id.csdl<-attr(tf,"specials")$csdl #non dipende se c'è intercetta o meno..
    if(length(id.csdl)>1) stop("Only one csdl() term is allowed!")
    id.seas<-attr(tf,"specials")$seas
    id.dl<-attr(tf,"specials")$dl
    
    #NB: (id.csdl) è la posizione della variabile csdl(x) nel model frame (include y ma non dipende da l'intercetta)
    toll<-fcontrol$toll
    visual<-fcontrol$visual
    it.max<-fcontrol$it.max
    GLM<-fcontrol$GLM
    maxit.glm<-fcontrol$maxit.inner
    #    
    ndx.seas<-0
    if(!is.null(id.csdl)){ #se c'è csdl() nella formula..
        if("z"%in%names(match.call(expand.dots = FALSE)$...)){
            new.mf <- match.call(expand.dots = TRUE)
            new.m<-match(c("formula", "data","z","subset", "na.action"), names(new.mf), 0)
            etic<-nomeTemp<-as.character(new.mf[["z"]])
                } else {
            etic<-nomeTemp<-all.vars(formula)[id.csdl]
            }
        only.seas<-FALSE
        testo.csdl<-names(mf)[id.csdl] #colnames(X)[id.csdl] #questa è la stringa del tipo csda(..)
        #testo.csdl<-grep("^csdl\\(",colnames(X),value=TRUE)
        nomiOrig.mf<-names(mf)
        names(mf)<-all.vars(formula)        
        z<-with(mf,eval(parse(text=testo.csdl)))
        names(mf)<-nomiOrig.mf #ripristina i nomi..
#        #estrai tutti gli argomenti per la modellazione della temp.
#        psi<-attr(z,"psi"); L<-attr(z,"L") ecc..
        #
        # se gli argomenti di csdl(), i.e. psi,L,ndx,ridge,DL, diff.varying sono messi come 
        #argomenti di tempeff(), questi sovrascrivono i corrispondenti di csdl(). (utile per update)
        new.mf <- match.call(expand.dots = FALSE)
        new.args<-new.mf$...
        psi<-L<-ndx<-DL<-ridge<-diff.varying<-heat.power<-NULL
        for(j in c("psi","L","ndx","DL","ridge","diff.varying","heat.power")){
          if(j %in% names(new.args)) assign(j,new.args[[j]]) else assign(j,attr(z,j))
        }
        #
        pcontrol<-list(DL=DL,diff.varying=diff.varying, ridge.formulas=ridge)
        #
        add.args.temp<-any(c("psi","L","ridge","ndx","DL","diff.varying","heat.power")%in%names(new.args))
        if(add.args.temp){
          fo.no.csdl<-update.formula(formula,as.formula(paste(".~.-",testo.csdl)))
          #agg.testo<-paste("csdl(",nomeTemp,",psi=",deparse(substitute(psi)), paste(",L=c(",L[1],",",L[2],")",sep=""),sep="")
          #
          agg.testo<-paste("csdl(",nomeTemp,",psi=",deparse(substitute(psi)), ",L=",deparse(substitute(L)),sep="")
          #
          if(length(grep("ridge",testo.csdl))>0 || "ridge"%in%names(new.args)){
            testo.ridge<-if(is.null(ridge)) "NULL" else
               paste("list(cold='", ridge[["cold"]],"'," ,"heat='",ridge[["heat"]] ,"')",sep="")
            agg.testo<-paste(agg.testo,",ridge=",testo.ridge,sep="")
            }
          if(length(grep("ndx",testo.csdl))>0 || "ndx"%in%names(new.args)){
            #agg.testo<-paste(agg.testo,",ndx=c(",ndx[1],",",ndx[2],")",sep="")
            agg.testo<-paste(agg.testo,",ndx=", deparse(substitute(ndx)),sep="")
            }
          if(length(grep("DL",testo.csdl))>0 || "DL"%in%names(new.args)){
            agg.testo<-paste(agg.testo,",DL=",DL,sep="")
            }
          if(length(grep("diff.varying",testo.csdl))>0 || "diff.varying"%in%names(new.args)){
            agg.testo<-paste(agg.testo,",diff.varying=",diff.varying,sep="")
            }
        nuovo.testo.csdl<-paste(agg.testo,")",sep="")
        new.formula<-update.formula(fo.no.csdl, as.formula(paste(".~.+",nuovo.testo.csdl)))
          }
        #
        if("z"%in%names(match.call(expand.dots = FALSE)$...)){
            new.mf <- match.call(expand.dots = TRUE)
            new.m<-match(c("formula", "data","z","subset", "na.action"), names(new.mf), 0)
            new.mf<-new.mf[c(1,new.m)]
            new.mf[[1]]<-as.name("model.frame")
            new.mf <- eval(new.mf, parent.frame())
            z<-new.mf[,length(all.vars(formula))+1]
            #new.mf<-paste("model.frame(~1,", "data=", mf[c(m[2])],",z=",match.call(expand.dots = FALSE)$...[["z"]],")")
            #new.mf <- eval(parse(text=new.mf), parent.frame())
            } else {
        z<-X[,testo.csdl]
        #oppure z<-X[,id.csdl]; oppure: z<-as.numeric(z) #(è la temperatura)
        #oppure: assign(nomeTemp, X[,id.csdl]);z<-eval(parse(text=nomeTemp))
        #
        }
        X<-X[,-match(testo.csdl,colnames(X))] #matrice del modello senza temp
        #
        ndx<-eval(ndx)
        L<-eval(L)
        yy<-y[-(1:max(L))]
        XX<-X[-(1:max(L)),]
        #
        psi<-eval(psi)
        k<-length(psi)
        if(any(psi<=min(z)) || any(psi>=max(z))) stop("Invalid starting values for psi")
        if(length(psi)>2) stop("One or two breakpoints allowed")
        if(length(psi)==1) psi<-c(psi,psi)
        psi.new<-psi
        n<-length(yy)
        PSI1 <- matrix(rep(rep(psi.new[1],(L[1]+1)), rep(n, L[1]+1)), ncol = L[1]+1)
        PSI2 <- matrix(rep(rep(psi.new[2],(L[2]+1)), rep(n, L[2]+1)), ncol = L[2]+1)
        Xlag<-as.matrix(sapply(0:max(L),function(i) lagged(z,i)))[-(1:max(L)),]
        U1<-pmin(Xlag[,1:(L[1]+1)]-PSI1,0)
        U2<-pmax(Xlag[,1:(L[2]+1)]-PSI2,0)^heat.power
        #
        A1<-bspline(0:L[1],ndx=ndx[1]) #B-spline basis for cold
        A2<-bspline(0:L[2],ndx=ndx[2]) #B-spline basis for heat
        colnames(A1)<-colnames(A2)<-NULL
        Xf<-U1%*%A1
        Xc<-U2%*%A2
        colnames(Xf)<-paste("Xf",1:ncol(Xf),sep="")
        colnames(Xc)<-paste("Xc",1:ncol(Xc),sep="")
    } else { #se *non* c'è csdl()
        only.seas<-TRUE
#        if(is.null(id.seas)) stop("without csdl(), at least provide seas()!")
        if(is.null(id.seas)) {#se non c'è smoothing..
        warning("Using glm.fit() for simple Poisson GLM fitting",call.=FALSE)
#        if(!missing(drop.L)) data<- data[-(1:drop.L),]
#        o<-glm(formula,data=data,family=poisson)
        #oppure
        if(!missing(drop.L)) {
            X<-as.matrix(X[-(1:drop.L),])
            y<-y[-(1:drop.L)]
            }
        o<-glm.fit(x=X, y=y, family=poisson())
        o$call<-match.call()
        o$formula<-o$call$formula
        o$df.residual<-length(o$residuals)-length(o$coefficients)
        class(o)<-"modTempEff"
        return(o)        
        }
    testo.seas<-names(mf)[id.seas] #colnames(X)[id.seas] #questa è la stringa del tipo "seas(..)"
    #Non serve valutare `seas()' come csdl(), mi interessa solo ndx.seas
    ndx.seas<-as.numeric(unlist(strsplit(unlist(strsplit(testo.seas,","))[2],"\\)")))
    X<-as.matrix(X[,-match(testo.seas,colnames(X))])
    if(!missing(drop.L)) {
            X<-as.matrix(X[-(1:drop.L),])
            y<-y[-(1:drop.L)]
            }
    o<-tempeff.fit(y,X,etastart=etastart,only.seas=only.seas,ndx.seas=ndx.seas) #tolto ...
    o$call<-match.call()
    o$formula<-o$call$formula
    o$df.residual<-length(o$residuals)-sum(o$edf)
    class(o)<-"modTempEff"
    return(o)
    } #end se *non* c'è csdl()
    if(!is.null(id.seas)){ #se c'è la seas
        names(mf)<-nomiOrig.mf #ripristina i nomi..
        testo.seas<-names(mf)[id.seas]#questa è la stringa del tipo "seas(..)"
        ndx.seas<-as.numeric(unlist(strsplit(unlist(strsplit(testo.seas,","))[2],"\\)")))        
        XX<-as.matrix(XX[,-match(testo.seas,colnames(XX))])
        }
    if(it.max==0){
        #se it.max=0 restituisce un modello con assegnati breakpoints==psi
        o<-tempeff.fit(yy,XX,Af=A1,Ac=A2,Xf,Xc,V=NULL,etastart=etastart,penalty=pcontrol,
          only.seas=only.seas,ndx.seas=ndx.seas) #tolto ...
        o$call<-match.call()
        o$formula<-if(!add.args.temp) formula else new.formula
        o$psi<-unique(psi)
        class(o)<-"modTempEff"
        return(o)
        }
    dev00<-glm.fit(x=XX, y=yy, family=poisson())$dev
    XREG<-cbind(XX,Xf,Xc)
    obj<-glm.fit(XREG, yy, family=poisson())
    id1<-grep("Xf",names(coef(obj)))
    id2<-grep("Xc",names(coef(obj)))
    a1<-coef(obj)[id1]
    a2<-coef(obj)[id2]
    beta1<- A1%*%a1
    beta2<- A2%*%a2
    initial <- psi
    it <- 1
    epsilon <- 10
    while(abs(epsilon)>toll) { #start_while
            eta0<-obj$linear.predictors
            U1<-pmin(Xlag[,1:(L[1]+1)]-PSI1,0)
            U2<-pmax(Xlag[,1:(L[2]+1)]-PSI2,0)^heat.power
            Xf<-U1%*%A1
            Xc<-U2%*%A2
            colnames(Xf)<-paste("Xf",1:ncol(Xf),sep="")
            colnames(Xc)<-paste("Xc",1:ncol(Xc),sep="")
            V1<-ifelse(Xlag[,1:(L[1]+1)] <PSI1, -1, 0)
            V2<-ifelse(Xlag[,1:(L[2]+1)] >PSI2, -1, 0)
            if(heat.power!=1) V2<-2*pmax(Xlag[,1:(L[2]+1)]-PSI2,0)*V2
            V1<-t(t(V1)*as.vector(beta1))
            V2<-t(t(V2)*as.vector(beta2))
            V1<-rowSums(V1)
            V2<-rowSums(V2)
            V<- if(k==1) V1+V2 else cbind(V1,V2)
            XREG<-cbind(XX,Xf,Xc,V)
            dev.old<-obj$dev
            #----stima il modello..
            if(GLM){
              obj<-suppressWarnings(glm.fit(XREG, yy, family=poisson(),
                  etastart=eta0,control=glm.control(maxit=maxit.glm)))
              a1<-obj$coef[grep("Xf",names(obj$coef))]
              a2<-obj$coef[grep("Xc",names(obj$coef))]
              beta1<- A1%*%a1
              beta2<- A2%*%a2
              d<-obj$coef[grep("V",names(obj$coef))]
              } else {
                obj<-tempeff.fit(yy,XX,Af=A1,Ac=A2,Xf,Xc,V=V,gam.fit.it=maxit.glm,etastart=eta0,ndx.seas=ndx.seas,
                  penalty=pcontrol) #tolto ...
                beta1<-obj$betaCold
                beta2<-obj$betaHeat
                d<- obj$delta
                }
            #----
            #controlla bene nel caso di 1-2 psi
            psi.old <- psi
            psi <- psi.old + d
#            if(psi<=min(z) || psi>=max(z)) stop("psi fuori dal range")
            #PSI <- matrix(rep(psi.new, rep(nrow(x), ncol(x))), ncol = ncol(x))
            PSI1 <- matrix(rep(rep(psi[1],(L[1]+1)), rep(n, L[1]+1)), ncol = L[1]+1)
            PSI2 <- matrix(rep(rep(psi[2],(L[2]+1)), rep(n, L[2]+1)), ncol = L[2]+1)
            dev.new <- obj$dev #if(GLM) obj$dev else obj$dev
        if (visual) {
            flush.console()
            if (it == 1)
                cat(0, " ", formatC(dev00, 3, format = "f"),
                  "", "---- without 'csdl' variable", "\n")
            spp <- if (it < 10) "" else NULL
            cat(it, spp, "", formatC(c(dev.new,unique(psi)), 3, format = "f"), "\n")
            }
        #epsilon <- (dev.new - dev.old)/dev.old #
        epsilon <- abs(dev.new-dev.old)/(abs(dev.new)+0.1)
        it <- it + 1
        if (it > it.max)
            break
        } #end_while
    if (it > it.max) warning("max number of iterations attained", call. = FALSE)
    obj<-tempeff.fit(yy,XX,Af=A1,Ac=A2,Xf,Xc,V=V,etastart=eta0,ndx.seas=ndx.seas,penalty=pcontrol) #tolto ...
    var.psi<- obj$Ve[obj$id.d,obj$id.d]
    var.psi.bayes<-  obj$Vp[obj$id.d,obj$id.d]
    if(k==2) {
        var.psi<-diag(var.psi)
        var.psi.bayes<-diag(var.psi.bayes)
          }
    psi<-cbind(unique(psi),sqrt(var.psi),sqrt(var.psi.bayes))
    colnames(psi)<-c("Est","SE.freq","SE.bayes")
    obj$psi<-psi
    obj$formula<-formula
    if(add.args.temp) obj$formula<-new.formula
    obj$df.residual<-n-sum(obj$edf)
    obj$call<-match.call()
    #if(add.args.temp) obj$call[["formula"]]<-new.formula
    class(obj)<-"modTempEff"
    return(obj)
    } #end_function..
