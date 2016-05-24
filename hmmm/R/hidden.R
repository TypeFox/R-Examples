


hidden.emfit<-function(y,model.obs,model.lat=NULL,nlat=1,
      noineq = TRUE, maxit = 10, maxiter=100, norm.diff.conv = 1e-05, 
    norm.score.conv = 1e-05, y.eps = 0,  mup = 1, step = 1,printflag=0,old.tran.p=NULL,
bb=NULL,q.par=1)
{
if(is.null(model.obs$Formula)){
hidden.emfit.S(y=y,model.obs=model.obs,model.lat=model.lat,nlat=nlat,
      noineq = noineq, maxit = maxit, maxiter=maxiter, norm.diff.conv = norm.diff.conv, 
    norm.score.conv = norm.score.conv, y.eps = y.eps,  mup = mup, step = step,
printflag=printflag,old.tran.p=old.tran.p,
bb=bb,q.par=q.par)}
else{hidden.emfit.X(y=y,model.obs=model.obs,model.lat=model.lat,nlat=nlat,
      noineq = noineq, maxit = maxit, maxiter=maxiter, norm.diff.conv = norm.diff.conv, 
    norm.score.conv = norm.score.conv, y.eps = y.eps,  mup = mup, step = step,
printflag=printflag,old.tran.p=old.tran.p,
bb=bb,q.par=q.par)}
}


##################################################################################################
hidden.emfit.X<-function(y,model.obs,model.lat=NULL,nlat=1,
      noineq = TRUE, maxit = 10, maxiter=100, norm.diff.conv = 1e-05, 
    norm.score.conv = 1e-05, y.eps = 0,  mup = 1, step = 1,printflag=0,old.tran.p=NULL,
bb=NULL,q.par=q.par)
{Nobs<-dim(y)[1]

#-------------------------------------------------------------------
#CORPO FUNZIONE
#INIZIALIZZAZIONE
#modalità osservate più latenti
if(!class(model.lat)=="hmmmmod"){
######---------------------------######
#0ld nstrata<-model.obs$modello$livelli[1]

#newNEWS
nstrata<-model.obs$lev.strata
nlat<-length(nstrata)
marglat<-paste(rep("l",2*nlat),collapse="-")
#old marglat<-c("l-l")
#marglat<-marg.list(marglat)
marglat0<-paste(rep(c("marg","l"),each=nlat),collapse="-")
namelat<-paste("lat",(1:nlat),sep="_")


namelat2<-paste("laglat",(1:nlat),sep="_")
namelat<-c(namelat,namelat2)

model.lat<-hmmm.model(marg=marg.list(c(marglat0,marglat)),lev=c(nstrata,nstrata),names=namelat,
#X=diag(1,nstrata*nstrata-1)
)
#togliere X=diag(1,nstrata*nstrata-1)
#print(model.lat)
######################################
}

lev<-model.obs$modello$livelli
levobs<-lev
#p-1 sono osservate attuali
p<-length(lev)

#modalità osservate
################---------------------###############
#old levobs<-lev[-1]

#numero stati latenti
#old strata<-model.obs$modello$livelli[1]
strata<-prod(model.obs$lev.strata)

###########################################################
#inizializzo prob transizione con distrib equiprobabilità
if(is.null(old.tran.p)){ 
#old.tran.p<-matrix(runif(strata^2,0,1),strata,strata)
old.tran.p<-matrix(1,strata,strata)
old.tran.p<-prop.table(old.tran.p,2)
}
A<-diag(1,strata)-old.tran.p
A<-rbind(A,matrix(1,1,strata))
old.iniz.p<-solve(t(A)%*%A)%*%t(A)
old.iniz.p<-old.iniz.p[,strata+1]


#inizializzo working data set  osservazioni complete: 
#latenti attuale e osservate
#r<-apply(y,1,function(x) matrix(rep(x,strata),ncol=length(x)))
#r<-matrix(r,ncol=dim(y)[2],byrow=TRUE)
r<-matrix(t(y),dim(y)[1]*strata,dim(y)[2],byrow=T)
##############-----------------------#################
oggi<-gl(strata,dim(y)[1],dim(r)[1])
#ieri<-gl(strata,strata,dim(r)[1])
r<-cbind(r,oggi)
r<-as.data.frame(r)

#inizializzo working data per latente attuale e ritardata
#la attuale varia più velocemente
oggi<-gl(strata,1,strata^2*dim(y)[1])
ieri<-gl(strata,strata,strata^2*dim(y)[1])
rrr<-cbind(oggi,ieri)
rrr<-as.data.frame(rrr)
#################################################

#inizializzo prob a posteriori degli stati
#ESTEP
#bb<-c(table(as.data.frame(y)))
if(is.null(bb)){
bb<-matrix(runif(prod(levobs)*strata,0,1),strata,prod(levobs))
bb<-prop.table(bb,1)

}
#bb<-c(table(as.data.frame(y)))
#bb<-matrix(bb,strata,length(bb),byrow=TRUE)
#bb<- prop.table(bb,1)

old.Ptobs<-bb
ii<-cumprod(c(1,levobs[-1]))
#ii<-t(apply(iiy,1,cumprod))
#/(max(y[,1]-1,1))
#ii<-cbind(ii,rep(1,dim(ii)[1])
#ii<-(cumprod(rev(ii))


yii<-(y-1)%*%ii+1
bb<-bb[,yii]


#bb<-matrix(1/dim(y)[1],nrow=2,ncol=dim(y)[1],byrow=TRUE)

ffs<-BLHK.filter(bb,t(old.tran.p),old.iniz.p)


#probabilità congiunte stati tempo t e tempo t-1
#smoothed
smooth.p<-array(
apply(X=ffs$pairs,FUN=function(x) t(x),
MARGIN=c(3))
,dim(ffs$pairs))


#inizilizzazioni per  condizioni convergenza
#######################################----------###################
#oldpar<-rep(999,prod(lev)-1+strata^2-1-2*(strata-1))
oldpar<-rep(999,(prod(lev)-1)*strata+strata^2-strata)
Lold<-999999

par.conv<-tran.conv<-L.conv<-99999999
#EM iterations----------------------------------------------
iter<-1
while ( (  ((par.conv> norm.diff.conv)&(tran.conv> norm.diff.conv)  ) ||(L.conv > norm.score.conv))
     &(iter< maxiter))
   {


#MSTEP

#tabella frequenze   congiunte delle osservate e latente attuale  per M step
r$count<-c(t(ffs$margin))

oldm<-xtabs(count~. ,data =r)
#edit(as.data.frame(oldm))

#tabella frequenze congiunte latente attuale e ritardata per M step
rrr$count<-c(smooth.p)
oldmmm<-xtabs(count~. ,data =rrr)





#stima modello hmmm definito da hmm.hmm.model usando come vettore frequenze congiunte oldm 
#sink("o")
fit.obs<-hmmm.mlfit(c(oldm), model.obs, maxit =maxit, norm.diff.conv =  norm.diff.conv, 
    norm.score.conv = norm.score.conv, y.eps =y.eps, 
    mup = mup, step = step) 
fit.lat<-hmmm.mlfit(c(oldmmm), model.lat, maxit =maxit, norm.diff.conv =  norm.diff.conv, 
    norm.score.conv = norm.score.conv, y.eps =y.eps, 
    mup = mup, step = step) 
#sink()
if(is.null(model.lat$Formula)){
p1<-fit.lat$L[-c(strata:(2*strata-2))]}
else{p1<-fit.lat$L}
p2<-fit.obs$L
newpar<-c(p1,p2)
#forma tabellare frequenze teoriche congiunte obs + latente attuale
######CRITICO°°°°°°°°°°°°
lev2<-c(prod(levobs),strata)
newm<-array(fit.obs$m,lev2)
newmmm<-array(fit.lat$m,c(strata,strata))
#da forma tabellare frequenze teoriche congiunte
#ricavo probabilità di transizione
new.tran.p<-prop.table(newmmm,2)
#calcolo distr.  invariante della latente
A<-diag(1,strata)-new.tran.p
A<-rbind(A,matrix(1,1,strata))
new.iniz.p<-solve(t(A)%*%A)%*%t(A)
new.iniz.p<-new.iniz.p[,strata+1]


#ESTEP
#in forma tabellare frequenze teoriche congiunte
#marginalizzo rispetto latente ritardata
#newm<-apply(newm,c(1:p,(p+2):length(lev)),sum)
#calcolo probabilità osservate dato stato latente attuale
########CRITICAL##################
b<-prop.table(newm,2)

#per ogni istante temporale calcolo:
#probabilità osservazioni tempo t 
#dato ogni stato latente
#b<-matrix(b,ncol=strata)
b<-t(b)
Ptobs<-b
ii<-cumprod(c(1,levobs[-1]))
#ii<-t(apply(iiy,1,cumprod))
#/(max(y[,1]-1,1))
#ii<-cbind(ii,rep(1,dim(ii)[1])
#ii<-(cumprod(rev(ii))


yii<-(y-1)%*%ii+1

bb<-b[,yii]
 
#FILTERING + SMOOTHING
# bb[,1]<-1
ffs<-BLHK.filter(bb,t(new.tran.p),new.iniz.p)


#probabilità congiunte stati tempo t e tempo t-1
#smoothed
smooth.p<-array(
apply(X=ffs$pairs,FUN=function(x) t(x),
MARGIN=c(3))
,dim(ffs$pairs))
#verosimiglianza osservazioni
Lnew<-sum(log(ffs$fy))
#CRITERI CONVERGENZA
par.conv<-quantile(abs(newpar-oldpar),q.par)

#par.conv<-mean(abs(newpar-oldpar))
L.conv<-abs(Lnew-Lold)
tran.conv<-quantile(abs(c(new.tran.p-old.tran.p,Ptobs-old.Ptobs)),1)
oldpar<-newpar
old.tran.p<-new.tran.p
Lold<-Lnew
old.Ptobs<-Ptobs
if(printflag > 0){
if(iter%%printflag==0){ cat("loglik= ",Lnew," Lconv= ",L.conv, 
" prob.conv= ",tran.conv, " niter= ",iter)
cat("\n")
}
}

iter<-iter+1
}  
#END EM iterations----------------------------------------------------------------------------------------------------
#restituisco working models, prob transizione, distr iniziale, 
#prob stati smmothed e filtrate, probabilità osservazioni tempo t dato passato , max verosim marginale, condizioni convergenza
#hmlist<-list(model=list(lat=fit.lat,obs=fit.obs),initial=new.iniz.p,
#Ptr<-new.tran.p,filter=ff,smooth=ffs, 
#Gsq=Lnew,conv=c(iter,par.conv,tran.conv,L.conv))
#dl<-hmmm.model.summary(model.lat,fit.lat,printflag=FALSE)
#do<-hmmm.model.summary(model.obs,fit.obs,printflag=FALSE)
#do<-do[-c(1:strata-1),]
#dl<-dl[-c(strata:(2*strata-2)),]
#dof<-(dim(model.obs$matrici$X)[1]+dim(model.lat$matrici$X)[1]-
#dim(model.obs$matrici$X)[2]-dim(model.lat$matrici$X)[2])
####################################dof<-fit.obs$df+fit.lat$df
fit.obs$Gsq<-Lnew
fit.obs$Xsq<-NaN
fit.obs$Wsq<-NaN
################################################fit.obs$df<-dof
fit.lat$Gsq<-NaN
fit.lat$Xsq<-NaN
fit.lat$Wsq<-NaN
#fit.lat$df<-


hmlist<-list(model.lat=fit.lat,model.obs=fit.obs,vecpar=list(lat=p1,obs=p2),initial=new.iniz.p,
Ptr=new.tran.p,Ptobs=Ptobs  ,smooth=ffs,conv=list(niter=iter,par.conv=par.conv,transprob.conv=tran.conv,Lconv=L.conv))
class(hmlist)<-"hidden"
hmlist
}


hidden.emfit.S<-function(y,model.obs,model.lat=NULL,nlat=1,
      noineq = TRUE, maxit = 10, maxiter=100, norm.diff.conv = 1e-05, 
    norm.score.conv = 1e-05, y.eps = 0,  mup = 1, step = 1,printflag=0,old.tran.p=NULL,
bb=NULL,q.par=q.par)
{Nobs<-dim(y)[1]

#-------------------------------------------------------------------
#CORPO FUNZIONE
#INIZIALIZZAZIONE
#modaltà osservate più latenti
if(!class(model.lat)=="hmmmmod"){
######---------------------------######
#0ld nstrata<-model.obs$modello$livelli[1]

#new
nstrata<-model.obs$modello$livelli[1:nlat]
marglat<-paste(rep("l",2*nlat),collapse="-")
#old marglat<-c("l-l")
#marglat<-marg.list(marglat)
marglat0<-paste(rep(c("marg","l"),each=nlat),collapse="-")
namelat<-paste("lat",(1:nlat),sep="_")


namelat2<-paste("laglat",(1:nlat),sep="_")
namelat<-c(namelat,namelat2)

model.lat<-hmmm.model(marg=marg.list(c(marglat0,marglat)),lev=c(nstrata,nstrata),names=namelat,
#X=diag(1,nstrata*nstrata-1)
)
#togliere X=diag(1,nstrata*nstrata-1)

######################################
}

#------------------------------------------------------------------------------
else{nlat<-length(model.lat$modello$livelli)/(1+is.null(model.lat$Formula))}
#------------------------------------------------------------------------------

lev<-model.obs$modello$livelli
#p-1 sono osservate attuali
p<-length(lev)

#modalità osservate
################---------------------###############
#old levobs<-lev[-1]
levobs<-lev[-(1:nlat)]
#numero stati latenti
#old strata<-model.obs$modello$livelli[1]
strata<-prod(model.obs$modello$livelli[1:nlat])

###########################################################
#inizializzo prob transizione con distrib equiprobabilità
if(is.null(old.tran.p)){ 
#old.tran.p<-matrix(runif(strata^2,0,1),strata,strata)
old.tran.p<-matrix(1,strata,strata)
old.tran.p<-prop.table(old.tran.p,2)
}
A<-diag(1,strata)-old.tran.p
A<-rbind(A,matrix(1,1,strata))
old.iniz.p<-solve(t(A)%*%A)%*%t(A)
old.iniz.p<-old.iniz.p[,strata+1]


#inizializzo working data set  osservazioni complete: 
#latenti attuale e osservate
r<-apply(y,1,function(x) matrix(rep(x,strata),ncol=length(x)))
r<-matrix(r,ncol=dim(y)[2],byrow=TRUE)
##############-----------------------#################
oggi<-gl(strata,1,dim(r)[1])
ieri<-gl(strata,strata,dim(r)[1])
r<-cbind(oggi,r)
r<-as.data.frame(r)

#inizializzo working data per letente attuale e ritardata
#la attuale varia più velocemente
oggi<-gl(strata,1,strata^2*dim(y)[1])
ieri<-gl(strata,strata,strata^2*dim(y)[1])
rrr<-cbind(oggi,ieri)
rrr<-as.data.frame(rrr)
#################################################

#inizializzo prob a posteriori degli stati
#ESTEP
#bb<-c(table(as.data.frame(y)))
if(is.null(bb)){
bb<-matrix(runif(prod(lev),0,1),strata,prod(levobs))
bb<-prop.table(bb,1)
}

old.Ptobs<-bb
#bb<-c(table(as.data.frame(y)))
#bb<-matrix(bb,strata,length(bb),byrow=TRUE)
#bb<- prop.table(bb,1)

ii<-cumprod(c(1,levobs[-1]))
#ii<-t(apply(iiy,1,cumprod))
#/(max(y[,1]-1,1))
#ii<-cbind(ii,rep(1,dim(ii)[1])
#ii<-(cumprod(rev(ii))


yii<-(y-1)%*%ii+1
bb<-bb[,yii]


#bb<-matrix(1/dim(y)[1],nrow=2,ncol=dim(y)[1],byrow=TRUE)

ffs<-BLHK.filter(bb,t(old.tran.p),old.iniz.p)


#probabilità congiunte stati tempo t e tempo t-1
#smoothed
smooth.p<-array(
apply(X=ffs$pairs,FUN=function(x) t(x),
MARGIN=c(3))
,dim(ffs$pairs))


#inizilizzazioni per  condizioni convergenza
#######################################----------###################
oldpar<-rep(999,prod(lev)-1+strata^2-1-2*(strata-1))
Lold<-999999

par.conv<-tran.conv<-L.conv<-99999999
#EM iterations----------------------------------------------
iter<-1
while ( (  ((par.conv> norm.diff.conv)&(tran.conv> norm.diff.conv)  ) ||(L.conv > norm.score.conv))
     &(iter< maxiter))
   {


#MSTEP

#tabella frequenze   congiunte delle osservate e latente attuale  per M step
r$count<-c(ffs$margin)

oldm<-xtabs(count~. ,data =r)
#edit(as.data.frame(oldm))

#tabella frequenze congiunte latente attuale e ritardata per M step
rrr$count<-c(smooth.p)
oldmmm<-xtabs(count~. ,data =rrr)





#stima modello hmmm definito da hmm.hmm.model usando come vettore frequenze congiunte oldm 
#sink("o")
fit.obs<-hmmm.mlfit(c(oldm), model.obs, maxit =maxit, norm.diff.conv =  norm.diff.conv, 
    norm.score.conv = norm.score.conv, y.eps =y.eps, 
    mup = mup, step = step) 
fit.lat<-hmmm.mlfit(c(oldmmm), model.lat, maxit =maxit, norm.diff.conv =  norm.diff.conv, 
    norm.score.conv = norm.score.conv, y.eps =y.eps, 
    mup = mup, step = step) 
#sink()
if(is.null(model.lat$Formula)){
p1<-fit.lat$L[-c(strata:(2*strata-2))]}
else{p1<-fit.lat$L}
#p1<-fit.lat$L[-c(strata:(2*strata-2))]
p2<-fit.obs$L[-c(1:strata-1)]
newpar<-c(p1,p2)
#forma tabellare frequenze teoriche congiunte obs + latente attuale
######CRITICO°°°°°°°°°°°°
lev2<-c(strata,lev[-(1:nlat)])
newm<-array(fit.obs$m,lev2)
newmmm<-array(fit.lat$m,c(strata,strata))
#da forma tabellare frequenze teoriche congiunte
#ricavo probabilità di transizione
new.tran.p<-prop.table(newmmm,2)
#calcolo distr.  invariante della latente
A<-diag(1,strata)-new.tran.p
A<-rbind(A,matrix(1,1,strata))
new.iniz.p<-solve(t(A)%*%A)%*%t(A)
new.iniz.p<-new.iniz.p[,strata+1]


#ESTEP
#in forma tabellare frequenze teoriche congiunte
#marginalizzo rispetto latente ritardata
#newm<-apply(newm,c(1:p,(p+2):length(lev)),sum)
#calcolo probabilità osservate dato stato latente attuale
########CRITICAL##################
b<-prop.table(newm,1)

#per ogni istante temporale calcolo:
#probabilità osservazioni tempo t 
#dato ogni stato latente
b<-matrix(b,nrow=strata)
Ptobs<-b
ii<-cumprod(c(1,levobs[-1]))
#ii<-t(apply(iiy,1,cumprod))
#/(max(y[,1]-1,1))
#ii<-cbind(ii,rep(1,dim(ii)[1])
#ii<-(cumprod(rev(ii))


yii<-(y-1)%*%ii+1

bb<-b[,yii]
 
#FILTERING + SMOOTHING
# bb[,1]<-1
ffs<-BLHK.filter(bb,t(new.tran.p),new.iniz.p)


#probabilità congiunte stati tempo t e tempo t-1
#smoothed
smooth.p<-array(
apply(X=ffs$pairs,FUN=function(x) t(x),
MARGIN=c(3))
,dim(ffs$pairs))
#verosimiglianza osservazioni
Lnew<-sum(log(ffs$fy))
#CRITERI CONVERGENZA
#par.conv<-mean(abs(newpar-oldpar))

par.conv<-quantile(abs(newpar-oldpar),q.par)
L.conv<-abs(Lnew-Lold)
tran.conv<-quantile(abs(c(new.tran.p-old.tran.p,Ptobs-old.Ptobs)),1)

#CRITERI CONVERGENZA
#par.conv<-max(abs(newpar-oldpar))
#L.conv<-quantile(abs(c(Lnew-Lold)),0.75)
#tran.conv<-max(abs(c(new.tran.p-old.tran.p)))
oldpar<-newpar
old.tran.p<-new.tran.p
old.Ptobs<-Ptobs
Lold<-Lnew
if(printflag > 0){
if(iter%%printflag==0){ cat("loglik= ",Lnew," Lconv= ",L.conv,
" prob.conv= ",tran.conv, " niter= ",iter)
cat("\n")
}
}

iter<-iter+1
}  
#END EM iterations----------------------------------------------------------------------------------------------------
#restituisco working models, prob transizione, distr iniziale, 
#prob stati smmothed e filtrate, probabilità osservazioni tempo t dato passato , max verosim marginale, condizioni convergenza
#hmlist<-list(model=list(lat=fit.lat,obs=fit.obs),initial=new.iniz.p,
#Ptr<-new.tran.p,filter=ff,smooth=ffs, 
#Gsq=Lnew,conv=c(iter,par.conv,tran.conv,L.conv))
#dl<-hmmm.model.summary(model.lat,fit.lat,printflag=FALSE)
#do<-hmmm.model.summary(model.obs,fit.obs,printflag=FALSE)
#do<-do[-c(1:strata-1),]
#dl<-dl[-c(strata:(2*strata-2)),]
#dof<-(dim(model.obs$matrici$X)[1]+dim(model.lat$matrici$X)[1]-
#dim(model.obs$matrici$X)[2]-dim(model.lat$matrici$X)[2])
#############################################dof<-fit.obs$df+fit.lat$df
fit.obs$Gsq<-Lnew
fit.obs$Xsq<-NaN
fit.obs$Wsq<-NaN
######################################################fit.obs$df<-dof
fit.lat$Gsq<-NaN
fit.lat$Xsq<-NaN
fit.lat$Wsq<-NaN
#fit.lat$df<-


hmlist<-list(model.lat=fit.lat,model.obs=fit.obs,vecpar=list(lat=p1,obs=p2),initial=new.iniz.p,
Ptr=new.tran.p,Ptobs=Ptobs,smooth=ffs,conv=list(niter=iter,par.conv=par.conv,transprob.conv=tran.conv,Lconv=L.conv))
class(hmlist)<-"hidden"
hmlist
}

#Baum-Lindgren-Hamilton-Kim filter and smoother
#Hans-Martin Krolzig (1997) : Markov Switching Vector autoregression-Springe#r

BLHK.filter <- function(pobs,pr,p0)
{
 #forward
  Tt <- ncol(pobs)
  m <- nrow(pobs)
  filterx <- matrix(p0,nrow=m,ncol=Tt)
  fy <- rep(1,Tt)
  filterx[,1] <- pobs[,1]*filterx[,1]
  fy[1] <- sum(filterx[,1])
  filterx[,1] <- filterx[,1]/fy[1]
  for (ti in (2:Tt)) {
    filterx[,ti] <- pobs[,ti]*(filterx[,ti-1]%*%pr)
    fy[ti] <- sum(filterx[,ti])
    filterx[,ti] <- filterx[,ti]/fy[ti]
  }
 
#backward



  r <- matrix(1,nrow=m,ncol=Tt)
  sm2 <- array(1,dim=c(m,m,Tt))
  for (ti in (Tt:2)) {
    mat <- t(matrix(r[,ti]*pobs[,ti],m,m))
    mat <- pr*mat/fy[ti]
    r[,ti-1] <- apply(mat,1,sum)
    sm2[,,ti-1] <- t( matrix(filterx[,ti-1],m,m)*mat)
  }
  sm2[,,1]<-NA
  list(fx=filterx, fy=fy,margin=r*filterx,pairs=sm2)
}











print.hidden<-function(x,printflag=FALSE,...){fitted<-x
#print(fitted$model.lat,printflag=printflag,aname="latent model",printhidden=TRUE)

#print(fitted$model.obs,printflag=printflag,aname="observation model",printhidden=1)
hmmm.model.summary(fitted$model.obs$model,fitted$model.obs,aname="observation model",printflag=printflag,printhidden=1)
if(printflag==TRUE){
cat("  effects relating  to  latent variables only, if printed ,","\n", "must not be considered" , "\n")}
hmmm.model.summary(fitted$model.lat$model,fitted$model.lat,aname="latent model",printflag=printflag,printhidden=2)
if(printflag==TRUE){
cat("  effects relating  to the lagged latent variables only, if printed","\n", "must not be considered" , "\n")}
}
summary.hidden<-function(object,...){
fitted<-object


#model.summary(fitted$model.obs,cell.stats=FALSE,model.info=FALSE)

if(fitted$model.obs$model$modello$strata >1){
cat("\n LATENT STATES  EFFECTS...")
  cat("\n")
  print(fitted$model.obs$beta)
}
#print(fitted$model.lat,printflag=TRUE,aname="latent model")
cat("\n","Transition probabilities","\n")
print(fitted$Ptr)
#print(fitted$model.obs,printflag=TRUE,aname="observation model")
 cat("\n","Probabilities of observations (by rows) given the latent states (by columns)","\n")
print(t(fitted$Ptobs))

}

hmm.hmm.anova<-function(modelloA,modelloB){
Gsq<-2*abs(modelloA$model.obs$Gsq-modelloB$model.obs$Gsq)
a<-modelloA$model.obs
alat<-modelloA$model.lat
if(class(a)=="hmmmfit"){
dfA=a$df+dim(a$Zlist$DMAT)[1]-dim(a$Zlist$DMAT)[2]-a$model$modello$strata+
alat$df+dim(alat$Zlist$DMAT)[1]-dim(alat$Zlist$DMAT)[2]-alat$model$modello$strata}
pA<-signif(1-pchisq(a$Gsq,dfA),5)

b<-modelloB$model.obs
blat<-modelloB$model.lat
if(class(b)=="hmmmfit"){
dfB=b$df+dim(b$Zlist$DMAT)[1]-dim(b$Zlist$DMAT)[2]-b$model$modello$strata+
blat$df+dim(blat$Zlist$DMAT)[1]-dim(blat$Zlist$DMAT)[2]-blat$model$modello$strata}
pB<-signif(1-pchisq(b$Gsq,dfB),5)


dof<-abs(dfA-dfB)
P<-round((1-pchisq(Gsq,dof)),5)
#print(matrix(c(round(a$Gsq,5),round(b$Gsq,5),round(Gsq,5),dfA,dfB,dof,"","",P),3,3,
#dimnames = list( c("model A", "model B","LR test"), c("statistics value", "dof","pvalue")  )    ),quote=FALSE)



anova.table<-matrix(c(round(a$Gsq,5),round(b$Gsq,5),round(Gsq,5),dfA,dfB,dof,NA,NA,P),3,3,
dimnames = list( c("model A", "model B","LR test"), c("statistics value", "dof","pvalue")  )    )


}
anova.hidden<-function(object,objectlarge,...){t<-
hmm.hmm.anova(object,objectlarge)
t}

akaike<-function(...,LRTEST=FALSE,ORDERED=FALSE,NAMES=NULL){
MOD<-list(...)
n<-length(MOD)
VAK<-matrix(0,1,8)
if(LRTEST){
if(class(MOD[[1]])=="hidden")
{
a<-MOD[[1]]$model.obs
b<-MOD[[1]]$model.lat
adf1<-a$df+dim(a$Zlist$DMAT)[1]-dim(a$Zlist$DMAT)[2]-a$model$modello$strata
bdf1<-b$df+dim(b$Zlist$DMAT)[1]-dim(b$Zlist$DMAT)[2]-b$model$modello$strata
df1<-adf1+bdf1
Gsq1<-2*MOD[[1]]$model.obs$Gsq
}

else{
df1=MOD[[1]]$df+dim(MOD[[1]]$Zlist$DMAT)[1]-dim(MOD[[1]]$Zlist$DMAT)[2]-MOD[[1]]$model$modello$strata
y<-MOD[[1]]$y 
Gsq1 <--MOD[[1]]$Gsq+ 2*sum(y[y>0]*log(y[y>0]))
}
}
for(i in 1:n) {
GSQ<-df<-P<-0
if(class(MOD[[i]])=="hidden"){
if(is.null(MOD[[i]])){ VAK<-rbind(VAK,c(i,rep(NA,7))) }
if(!is.null(MOD[[i]]))
{
a<-MOD[[i]]$model.obs
b<-MOD[[i]]$model.lat
adf<-a$df+dim(a$Zlist$DMAT)[1]-dim(a$Zlist$DMAT)[2]-a$model$modello$strata
bdf<-b$df+dim(b$Zlist$DMAT)[1]-dim(b$Zlist$DMAT)[2]-b$model$modello$strata
npar<-length(MOD[[i]]$vecpar$obs)+length(MOD[[i]]$vecpar$lat)-adf-bdf
AK<-2*npar-2*a$Gsq

#################
if(LRTEST){
GSQ<- -2*a$Gsq+Gsq1
df<-adf+bdf-df1
P<-1-pchisq(GSQ,df)
}
#########################
#AK<-c(i,a$Gsq,npar,AK)
VAK<-rbind(VAK,c(i,a$Gsq,adf+bdf,npar,GSQ,df,P,AK))
}
}
else{
dfi=MOD[[i]]$df+dim(MOD[[i]]$Zlist$DMAT)[1]-dim(MOD[[i]]$Zlist$DMAT)[2]-MOD[[i]]$model$modello$strata
npar<-length(MOD[[i]]$y)-MOD[[i]]$model$modello$strata-dfi


y<-MOD[[i]]$y   
Gsqi <--MOD[[i]]$Gsq/2+ sum(y[y>0]*log(y[y>0]))
AK<--2*Gsqi+2*npar
###########################################################
if(LRTEST){
df<-dfi-df1
GSQ<--2*Gsqi+Gsq1
P<-1-pchisq(GSQ,df)
}
########################################################
VAK<-rbind(VAK,c(i,Gsqi,MOD[[i]]$df,npar,GSQ,df,P,AK))

}
}
VAK<-VAK[-1,]

Delta<-VAK[,8]-min(VAK[,8])
VAK<-cbind(VAK,matrix(Delta,length(Delta),1))
VAK[,-1]
if(is.null(NAMES)){
rownames(VAK)<-paste("model",1:n,sep="")}
else{rownames(VAK)<-NAMES}
if(LRTEST){
colnames(VAK)<-c("#model","loglik","dfmodel","npar","LRTEST","dftest","PVALUE","AIC","DeltaAIC")
}
else{

VAK<-VAK[,c(1,2,3,4,8,9)]
colnames(VAK)<-c("#model","loglik","dfmodel","npar","AIC","DeltaAIC")
}
if(ORDERED){r<-order(VAK[,dim(VAK)[2]])
VAK[r,]}
else{VAK}


}

