`tempeff.fit` <-
function(y, X, Af=NULL, Ac=NULL, Xf=NULL, Xc=NULL, V=NULL, ndx.seas=0, only.seas=FALSE, 
      penalty=list(DL=FALSE,diff.varying=FALSE,ridge.formulas=NULL), 
      gam.fit.it=NULL, etastart=NULL, spstart=NULL, fit.method="magic"){
#--DUBBI: spstart (valori iniziali per sp) e gam.fit.it (numero max di iterazioni) pare che vengano ignorati da
#       gam.fit..???
#y: la risposta (vettore)
#X: la matrice di esplicative (escluso temperatura)
#Af, Ac: le B-spline basi
#Xf, Xc: le pseudo esplicative per il freddo e caldo (ad es., Xf=Xlagged.f%*%Af)
#V: eventuale matrice V per il punto di svolta
#seasonP: se diverso da 0 la stagionalità è modellata con P-spline con quel dato ndx
#fit.metodo=un carattere "magic" o "mgcv"..Da un errore se è "mgcv" ?
#only.seas: se TRUE restituisce un modello con le esplicative e soltanto la stagionalità
#   modellata come P-splines. In tal caso seasonP viene imposto TRUE.
#penalty: una lista con componenti:
#   1)"DL" specifica se la difference penalty deve riguardare i coeff della spline (DL=FALSE)
#   2)"diff.varying" specifica se considerare la penalità sulla curvatura variabile. Valori ammissibili
#   "no", "spline", "DL"
#   3)"ridge" applicare ridge penalties ai DL coeff?
#ridge.formulas: una list di due caratteri *nominati* per specificare la formula in termini di xlag della varying
#   ridge penalty per freddo e caldo. Ignorato se penalty$ridge=FALSE
#   bspline(xlag,ndx=3,deg=1)[,-1] #restituisce "ndx+deg" nodi per "ndx-1" inner knots
##i mancanti già eliminati..
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

#    require(mgcv)
    if(!is.matrix(X)) stop("X has to be a matrix")
    if(only.seas && ndx.seas==0) stop("if only.seas is TRUE it should be 'ndx.seas>0' ")
    n<-length(y)
#    test.dim<-identical(n,nrow(X))+identical(n,nrow(Xf))+identical(n,nrow(Xc))
#    if(!isTRUE(identical(test.dim,3))) stop("y, X, Xf, Xc have to be of the same dimension")
    S<-list()
    off.seas<- NULL
    nome.seas<- NULL
    if(ndx.seas>0){
        if(isTRUE(all.equal(X[,1],rep(1,nrow(X)),check.attributes = FALSE))) X<-X[,-1]
        B.seas<-bspline(1:n, ndx.seas, deg=3)
        off.seas<-ncol(X)+1
        X<-cbind(X, B.seas)
        S<-list(crossprod(diff(diag(ncol(B.seas)), diff = 2)))
        nome.seas<-"lambda.seas"
        }
    if(only.seas){
      M<-list(y=y,X=X,S=S,off=off.seas,sp=rep(-1,length(S)),
        am=FALSE,intercept=TRUE,fit.method=fit.method,w=rep(1,n),offset=rep(0,n),
        sig2=1,conv.tol=1e-07,max.half=15)
      o<-gam.fit(M,family=poisson())
      id.seas<-seq(off.seas,length.out=ncol(B.seas))
      o$fit.seas<-drop(B.seas%*%o$coef[id.seas])
      o$edf.seas<-o$edf[id.seas]
      o$rank.seas <- ncol(B.seas)
      return(o)
        }
    if(!is.null(V)){
        V<-as.matrix(V)
        n.d<-ncol(V)
        X<-cbind(X,V)
        }
    XREG<-cbind(X,Xf,Xc)
    Df<-diff(diag(ncol(Af)), diff = 2) #difference matrix for freddo
    Dc<-diff(diag(ncol(Ac)), diff = 2) #difference matrix for caldo
    #nomi<-c("l0.Cold","l0.Heat")
    nomi<-c(nome.seas,"lambda.Cold","lambda.Heat")
    if(penalty$DL) {
        Df<-diff(diag(nrow(Af)), diff = 2)%*%Af
        Dc<-diff(diag(nrow(Ac)), diff = 2)%*%Ac
        }
    #S<-list(crossprod(Df),crossprod(Dc))
    S[[length(S)+1]]<-crossprod(Df)
    S[[length(S)+1]]<-crossprod(Dc)
    off.f<-ncol(X)+1
    off.c<-ncol(X)+ncol(Xf)+1
    #tutti.off<-c(off.f,off.c)
    tutti.off<-c(off.seas,off.f,off.c)
    M<-list(y=y,X=XREG,S=S,off=tutti.off,sp=rep(-1,length(S)),
        am=FALSE,intercept=TRUE,fit.method=fit.method,w=rep(1,n),offset=rep(0,n),
        sig2=1,conv.tol=1e-07,max.half=15)

    #-----SE SEVI METTERE LA DIFF PENALTY VARIABILE------------------------------------------
    if(penalty$diff.varying==TRUE){ #if varying diff penalty
    #una Bspline "povera" (lineare con 2 nodi interni)
        xx.f<-1:nrow(Df)
        if(penalty[[2]]!="spline"){
                S<-list(crossprod(diff(diag(length(xx.f)+2), diff = 2)%*%Af),
                    crossprod(diag((xx.f))%*%diff(diag(length(xx.f)+2),diff = 2)%*%Af))
                    } else {
                S<-list()
                basexx.f<-bspline(xx.f,ndx=3,deg=1)
                for(i in 1:ncol(basexx.f)) S[[length(S)+1]]<-crossprod(diag(sqrt(basexx.f[,i]))%*%Df)
                    }
        nomi.freddo<-paste("lambda",1:length(S),".Cold",sep="")
        #adesso il caldo:
        xx.c<-1:nrow(Dc)
        if(penalty[[2]]!="spline"){
                    S[[length(S)+1]]<-crossprod(diff(diag(length(xx.c)+2), diff = 2)%*%Ac)
                    S[[length(S)+1]]<-crossprod(diag((xx.c))%*%diff(diag(length(xx.c)+2),diff = 2)%*%Ac)
                    } else {
                    basexx.c<-bspline(xx.c,ndx=3,deg=1)
                    for(i in 1:ncol(basexx.c)) S[[length(S)+1]]<-crossprod(diag(sqrt(basexx.c[,i]))%*%Dc)
                    }
        #basexx.c<-bspline(xx.c,ndx=3,deg=1)
        #for(i in 1:ncol(basexx.c)) S[[length(S)+1]]<-crossprod(diag(sqrt(basexx.c[,i]))%*%Dc)
        nomi.caldo<-paste("lambda",1:(length(S)-length(nomi.freddo)),".Heat",sep="")
        M$S<-S
        M$off<-if(penalty[[2]]!="spline")
            {c(rep(off.f,2),rep(off.c,2))} else {
                c(rep(off.f,ncol(basexx.f)), rep(off.c,ncol(basexx.c)))}
        M$sp<-rep(-1,length(S))
        nomi<-c(nomi.freddo,nomi.caldo)
        } #end diff-varying
    #-----FINE DIFF PENALTY VARIABILE------

    if(suppressWarnings(!is.null(penalty$ridge.formulas[1]))){ #if ridge varying penalty..
        xx.rid.f<-1:nrow(Af)
        xx.rid.c<-1:nrow(Ac)

        f.freddo<-function(l) eval(parse(text=penalty$ridge.formulas$cold))
        f.caldo<-function(l)  eval(parse(text=penalty$ridge.formulas$heat))
        xx.rid.f.valori<-f.freddo(xx.rid.f)
        xx.rid.c.valori<-f.caldo(xx.rid.c)
        if(is.matrix(xx.rid.f.valori)){
            idfreddo<-seq(length(S)+1,length.out=ncol(xx.rid.f.valori))
            S[idfreddo]<-lapply(data.frame(xx.rid.f.valori),function(yy){crossprod(diag(sqrt(yy))%*%Af)})
            off.f<-rep(off.f,length(idfreddo))
            nomiFr.ridge<-paste("omega.Cold",1:length(idfreddo))
            } else {
              S[[length(S)+1]]<-crossprod(diag(sqrt(xx.rid.f.valori))%*%Af)
              nomiFr.ridge<-"omega.Cold"
              }
        if(is.matrix(xx.rid.c.valori)){
            idcaldo<-seq(length(S)+1,length.out=ncol(xx.rid.c.valori))
            S[idcaldo]<-lapply(data.frame(xx.rid.c.valori),function(yy){crossprod(diag(sqrt(yy))%*%Ac)})
            off.c<-rep(off.c,length(idcaldo))
            nomiCa.ridge<-paste("omega.Heat",1:length(idcaldo))
            } else {
              S[[length(S)+1]]<-crossprod(diag(sqrt(xx.rid.c.valori))%*%Ac)
              nomiCa.ridge<-"omega.Heat"
              }
        M$S<-S
        M$off<-c(M$off,off.f,off.c)
        M$sp<-if(is.null(spstart)) rep(-1,length(S)) else spstart
        nomi<-c(nomi,nomiFr.ridge,nomiCa.ridge)
        }
    o<-if(!is.null(gam.fit.it)) {
            gam.fit(M,family=poisson(),fixedSteps=gam.fit.it,etastart=etastart)
            } else {
            gam.fit(M,family=poisson(),etastart=etastart)
            }
    o$sp.mio<-matrix(o$sp,ncol=1,dimnames=list(nomi,"spar"))
    o$id.f<-id.f<-seq(ncol(X)+1,length=ncol(Af)) #id coeff cold
    o$id.c<-id.c<-seq(ncol(X)+ncol(Xf)+1,length=ncol(Ac))  #id coeff heat
    if(!is.null(V)){
        o$id.d<-id.d<-rev(seq(ncol(X),by=-1,length=n.d)) #id coeff gamma
        o$delta<-o$coefficients[id.d]
        o$Tdelta<-c(coef=o$delta,tvalue=o$delta/sqrt(diag(as.matrix(o$Ve[id.d,id.d]))))
        }
    o$betaCold<-Af%*%o$coef[id.f]
    o$betaHeat<-Ac%*%o$coef[id.c]
    #o$Af<-Af
    #o$Ac<-Ac
    Var<- o$Ve #o Vp is bayesian (larger)
    o$SE.c<-sqrt(diag(Af%*%Var[id.f,id.f]%*%t(Af)))
    o$SE.h<-sqrt(diag(Ac%*%Var[id.c,id.c]%*%t(Ac)))
    o$ToTcold<-c(-sum(o$betaCold),
      sqrt(drop(rep(1,nrow(Af))%*%Af%*%Var[id.f,id.f]%*%t(Af)%*%rep(1,nrow(Af)))))
    o$ToTheat<-c(sum(o$betaHeat),
      sqrt(drop(rep(1,nrow(Ac))%*%Ac%*%Var[id.c,id.c]%*%t(Ac)%*%rep(1,nrow(Ac)))))
    #calcolo var bayesiana
    Var<- o$Vp
    o$SE.c.bayes<-sqrt(diag(Af%*%Var[id.f,id.f]%*%t(Af)))
    o$SE.h.bayes<-sqrt(diag(Ac%*%Var[id.c,id.c]%*%t(Ac)))
    o$ToTcold.bayes<-c(-sum(o$betaCold),
      sqrt(drop(rep(1,nrow(Af))%*%Af%*%Var[id.f,id.f]%*%t(Af)%*%rep(1,nrow(Af)))))
    o$ToTheat.bayes<-c(sum(o$betaHeat),
      sqrt(drop(rep(1,nrow(Ac))%*%Ac%*%Var[id.c,id.c]%*%t(Ac)%*%rep(1,nrow(Ac)))))
    o$edf.cold<-o$edf[id.f]; o$rank.cold<-ncol(Af)
    o$edf.heat<-o$edf[id.c]; o$rank.heat<-ncol(Ac)
    if(ndx.seas>0){
      id.seas<-seq(off.seas,length.out=ncol(B.seas))
      o$fit.seas<-drop(B.seas%*%o$coef[id.seas])
      o$edf.seas<-o$edf[id.seas]
      o$rank.seas <- ncol(B.seas)
      }
    o$call<-match.call()
    o}

