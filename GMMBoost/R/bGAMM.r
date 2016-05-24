bs.design<-function(x, diff.ord, spline.degree, knots.no)
{

    ## generate a B-Spline-Matrix with equidistant knots (code by Thomas Kneib):
    n<-length(x)
    xl<-min(x)
    xr<-max(x)
    xmin<-xl-(xr-xl)/100
    xmax<-xr+(xr-xl)/100
    dx<-(xmax-xmin)/(knots.no-1)
    knots<-seq(xmin-spline.degree*dx,xmax+spline.degree*dx,by=dx)
    B<-splines::spline.des(knots,x,spline.degree+1)$design
    knots<-knots[1:(knots.no+spline.degree-1)]

    ## generate Penalization-Matrix
    D<-diag(length(knots))
    d<-min(diff.ord,spline.degree)
    if(d<diff.ord) warning(paste("order of differences > degree of
    splines:\n new order of differences=",d,"\n"))
    if(d>0) {for(i in 1:d) D<-diff(D)}

    ## reparametrization: B_unpen = unpenalized part, B_pen=penalized part
    B.unpen.fact<-rep(1,length(knots))
    if(diff.ord>1) {for(i in 2:diff.ord)
                        B.unpen.fact<-cbind(B.unpen.fact,knots^(i-1)) }

    B.unpen<-(B%*%B.unpen.fact)
    B.pen  <-B%*%t(D)%*%solve(D%*%t(D))

    return(list(B=B, X=B.unpen, Z=B.pen))
}


#############################################################################################################################################
###################################################### main boosting function #################################################################
#fix=formula(y~1);add=formula(~u1+u2+u3+u4+u5+u6+u7+u8+u9+u10+u11+u12+u13+u14+u15);rnd=list(Person=~1);data=Hirst;lambda=10;family=poisson(link=log);control=list(diff.ord=4,spline.degree=6,nbasis=15,add.fix="u7")
#fix=formula(points~ball.possession+as.factor(yellow.red.card));add=formula(points~transfer.spendings+transfer.receits+unfair.score+sold.out+ave.attend);rnd=list(team=~1+tackles);data=soccer;lambda=1001;family=poisson(link=log);
#control = list(nbasis=15,spline.degree=2,diff.ord=1,add.fix="ave.attend",overdispersion=TRUE,start=c(5,rep(0,52)),sel.method="bic",method="REML")


est.bGAMM<-function(fix,add,rnd,data,lambda,family=gaussian(link = "identity"),control=list())
{
if(attr(terms(fix), "intercept")==0)
{
variables<-attr(terms(fix),"term.labels")
fix<- paste(rownames((attr(terms(fix),"factors")))[1],"~ +1",sep="")
for (ir in 1:length(variables))
fix <- paste(fix, variables[ir], sep="+")
fix <-formula(fix)
}

y <- model.response(model.frame(fix, data))

X <- model.matrix(fix, data)
X<-as.matrix(X)

lin.old<-dim(X)[2]

if(attr(terms(add), "intercept")==0)
{
variables<-attr(terms(add),"term.labels")
add<- paste(rownames((attr(terms(add),"factors")))[1],"~ +1",sep="")
for (ir in 1:length(variables))
add <- paste(add, variables[ir], sep="+")
add <-formula(add)
}

U <- model.matrix(add, data)
U<-U[,-1]

m.old<-dim(U)[2]


very.old.names<-attr(terms(fix),"term.labels")
factor.names<-character()
for (i in 1:length(very.old.names))
{
if (substr(very.old.names[i],1,9)=="as.factor")
factor.names<-c(factor.names,very.old.names[i])
}

factor.list<-list()
if(length(factor.names)>0)
{
spl<-strsplit(factor.names,"\\(")

categ.names<-character()
for(uz in 1:length(spl))
categ.names<-c(categ.names,spl[[uz]][2])

spl2<-strsplit(categ.names,"\\)")
categ.names2<-character()
for(uz in 1:length(spl2))
categ.names2<-c(categ.names2,spl2[[uz]])
factor.names<-categ.names2

for (i in 1:length(factor.names))
factor.list[[i]]<-levels(as.factor(data[,factor.names[i]]))
}    



old.names<-attr(X,"dimnames")[[2]]

U.names<-attr(U,"dimnames")[[2]]

rndformula <- as.character(rnd)

trmsrnd <- terms(rnd[[1]])
newrndfrml <- "~ -1"
newrndfrml <- paste(newrndfrml,  if(attr(trmsrnd, "intercept")) names(rnd)[1] else "", sep=" + ")

if(length(attr(trmsrnd, "variables"))>1)
{
newrndfrml <- paste(newrndfrml,
         paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                     paste(lbl, names(rnd)[1], sep=":")
                 }), collapse=" + "), sep="+") }


W_start <- model.matrix(formula(newrndfrml), data)

rnlabels<-terms(formula(newrndfrml))
random.labels<-attr(rnlabels,"term.labels")
k<-table(data[,colnames(data)==(names(rnd)[1])])
n<-length(k)
s<-dim(W_start)[2]/n

if(s>1)
{
W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
for (i in 2:n)
W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
}else{
W<-W_start
}

control<-do.call(bGAMMControl, control)

nbasis<-control$nbasis
diff.ord<-control$diff.ord
spline.degree<-control$spline.degree
knots.no<-nbasis-(spline.degree-2)


if(length(attr(trmsrnd,"factors"))!=0 && !is.element(colnames(attr(trmsrnd,"factors")),very.old.names))
stop("Use only fixed effects as random slopes!")

if(!is.null(control$add.fix) && !is.element(control$add.fix,U.names))
stop("Forced additive terms must be also part of the additive formula!")

if(!(diff.ord<spline.degree))
stop("Order of differences must be lower than degree of B-spline polynomials!")


print(paste("Iteration ", 1,sep=""))


nbasis.old<-nbasis


U.fix<-NULL
dim.U.fix<-0
if(!is.null(control$add.fix))
{
add.out<-!is.element(U.names,control$add.fix)
U.fix<-U[,!add.out]
U.fix<-as.matrix(U.fix)
colnames(U.fix)<-U.names[!add.out]
U.names<-U.names[add.out]
dim.U.fix<-dim(U.fix)[2]
U<-U[,add.out]
}

m<-dim(U)[2]

Basis<-list()
Phi<-list()

for (r in 1:m)
{
Basis[[r]]<-bs.design(U[,r],diff.ord=diff.ord,spline.degree=spline.degree,knots.no=knots.no)
Phi[[r]]<-cbind(Basis[[r]]$X[,-1],Basis[[r]]$Z)
colnames(Phi[[r]])<-paste(colnames(U)[r],rep(1:dim(Phi[[r]])[2],each=1), sep=".")
}

if(!is.null(U.fix))
{
Basis.fix<-list()
Phi.fix<-list()
for (rr in 1:dim.U.fix)
{
Basis.fix[[rr]]<-bs.design(U.fix[,rr],diff.ord=diff.ord,spline.degree=spline.degree,knots.no=knots.no)
Phi.fix[[rr]]<-cbind(Basis.fix[[rr]]$X[,-1],Basis.fix[[rr]]$Z)
colnames(Phi.fix[[rr]])<-paste(colnames(U.fix)[rr],rep(1:dim(Phi.fix[[rr]])[2],each=1), sep=".")
}
for (rr in 1:dim.U.fix)
X<-cbind(X,Phi.fix[[rr]])
}


lin<-dim(X)[2]

if(length(control$start)==0)
control$start<-c(rep(0,(lin+m*nbasis+n*s)))

if(length(control$q_start)==0)
{
control$q_start<-rep(0.1,s)
if(s>1)
control$q_start<-diag(control$q_start)
}

beta_null<-control$start[1:lin.old]
if(!is.null(U.fix))
beta_null<-c(beta_null,rep(0,dim.U.fix*nbasis))

ranef_null<-control$start[(lin.old+1):(lin.old+n*s)]
q_start<-control$q_start

k1<-rep(0,lin.old)
if(diff.ord>1)
{
k2<-c(rep(0,diff.ord-1),rep(1,ncol(Basis[[r]]$Z)))
}else{
k2<-rep(1,ncol(Basis[[r]]$Z))
}
k22<-rep(k2,dim.U.fix+1)
k3<-rep(0,s*n)
K<-lambda*diag(c(k1,k22,k3))

N<-length(y)

Z_fastalles<-X

lin0<-sum(beta_null!=0)
if(s==1)
{
Q_start<-diag(q_start^2,s)
p_start<-c(rep(0,lin0),rep((q_start^2)^(-1),n*s))
P_start<-diag(p_start)
}else{
Q_start<-q_start
P_start<-matrix(0,lin0+n*s,lin0+n*s)
for(jf in 1:n)
P_start[(lin0+(jf-1)*s+1):(lin0+jf*s),(lin0+(jf-1)*s+1):(lin0+jf*s)]<-chol2inv(chol(q_start))
}

if(lin>1)
{
Eta_start<-X%*%beta_null+W%*%ranef_null
}else{
Eta_start<-rep(beta_null,N)+W%*%ranef_null
}

D_start<-as.vector(family$mu.eta(Eta_start))
Mu_start<-as.vector(family$linkinv(Eta_start))
Sigma_start<-as.vector(family$variance(Mu_start))
W0_inv<-D_start*1/Sigma_start*D_start

Z_start<-cbind(X[,beta_null!=0],W)
M0<-Z_start%*%chol2inv(chol((t(Z_start)%*%(Z_start*W0_inv))+P_start))%*%t(Z_start*W0_inv)

Mu<-as.vector(family$linkinv(Eta_start))

for (r in 1:m)
Z_fastalles<-cbind(Z_fastalles,Phi[[r]])

Z_alles<-cbind(Z_fastalles,W)

#########################  some definitions ################
Z_r<-list()
for (r in 1:m)
Z_r[[r]]<-cbind(X,Phi[[r]])

Z_r_gross<-list()
for (r in 1:m)
Z_r_gross[[r]]<-cbind(Z_r[[r]],W)

Sigma<-as.vector(family$variance(Mu))

if(control$overdispersion==T)
{
phi_start<-(sum((y-Mu)^2/Mu))/(N-sum(diag(M0)))
Sigma<-Sigma*phi_start
}
########################################################## some definitions ################################################

Delta<-matrix(0,control$steps,(lin+m*nbasis+s*n))

Standard<-array(0,c(control$steps,m,N))
Standard_ma<-list()

for (ijk in 1:m)
Standard_ma[[ijk]]<-as.numeric()

Delta[1,1:lin]<-beta_null
Delta[1,(lin+m*nbasis+1):(lin+m*nbasis+s*n)]<-t(ranef_null)

Delta_start<-Delta[1,]

IC<-matrix(0,nrow=control$steps,ncol=m)
IC_stop<-rep(0,control$steps)
komp<-rep(0,control$steps)

E<-diag(N)

M1<-list()
diag_H1<-rep(0,m)
FinalHat<-list()

phi<-rep(0,control$steps)

Q<-array(0,c(s,s,control$steps+1))
Q[,,1]<-Q_start
######### derive D and score_r and F_r  ###########
D<-as.vector(family$mu.eta(Eta_start))

Delta_r<-matrix(0,nrow=lin+nbasis+s*n,ncol=m)
Eta_r<-matrix(0,N,m)
Mu_r<-matrix(0,N,m)

grp<-rep(1,m*nbasis)
for (r in 2:m)
{
grp[((r-1)*nbasis+1):(r*nbasis)]<-r
}

if(s==1)
{
P_alles<-c(rep(0,lin+m*nbasis),rep((Q_start^(-1)),n*s))
P_alles<-diag(P_alles)
}else{
P_alles<-matrix(0,lin+m*nbasis+n*s,lin+m*nbasis+n*s)
for(jf in 1:n)
P_alles[(lin+m*nbasis+(jf-1)*s+1):(lin+m*nbasis+jf*s),(lin+m*nbasis+(jf-1)*s+1):(lin+m*nbasis+jf*s)]<-chol2inv(chol(Q_start))
}

score_alles<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
F_alles<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P_alles

l=1
for (r in 1:m)
{
Inv_r<-chol2inv(chol(F_alles[c(1:lin,(lin+(r-1)*nbasis+1):(lin+r*nbasis),(lin+m*nbasis+1):(lin+m*nbasis+n*s)),c(1:lin,(lin+(r-1)*nbasis+1):(lin+r*nbasis),(lin+m*nbasis+1):(lin+m*nbasis+n*s))]+K))
Delta_r[,r]<-Inv_r%*%score_alles[c(1:lin,(lin+(r-1)*nbasis+1):(lin+r*nbasis),(lin+m*nbasis+1):(lin+m*nbasis+n*s))]

############# derive predictor Eta_r ########
Eta_r[,r]<-Eta_start+Z_r_gross[[r]]%*%Delta_r[,r]

############# new  Mu_r ##############
Mu_r[,r]<-as.vector(family$linkinv(Eta_r[,r]))

################### derive Hat-Matrix  ##############
M1[[r]]<-(Z_r_gross[[r]]*sqrt(Sigma*D*1/Sigma*D))%*%Inv_r%*%t(Z_r_gross[[r]]*sqrt(D*1/Sigma*D*1/Sigma))
diag_H1[r]<-sum(diag(M0+M1[[r]]%*%(E-M0)))

######## likelihood of information criterion   ########
    if (family$family=="poisson")
    IC[1,r]<-sum(y*log(Mu_r[,r])-Mu_r[,r])#-0.5*pena

    if (family$family=="binomial")
    IC[1,r]<-sum(y*log(Mu_r[,r])+(1-y)*log(1-Mu_r[,r]))#-0.5*pena

    if (family$family=="gaussian")
    IC[1,r]<-sum(y*Mu_r[,r]-0.5*(Mu_r[,r]^2))#-0.5*pena

   ############ AIC ##############
    if (control$sel.method=="aic")
    IC[1,r]<-(-IC[1,r])+diag_H1[r]

   ############ BIC ###############
    if (control$sel.method=="bic")
    IC[1,r]<-(-IC[1,r])+0.5*diag_H1[r]*log(n)
}

   mi_ic<-match(min(IC[1,]),IC[1,])

#################### update ###########################
g<-c(rep(T,lin),grp==mi_ic,rep(T,s*n))

Delta_r[(lin+1):(lin+nbasis),mi_ic]<-control$nue*Delta_r[(lin+1):(lin+nbasis),mi_ic]
Delta[1,g]<-Delta[1,g]+t(Delta_r[,mi_ic])

Nue_ma<-c(rep(1,lin),rep(control$nue,nbasis),rep(1,s*n))
Uu<-(E-(cbind(Z_r[[mi_ic]],W)*sqrt(Sigma*D*1/Sigma*D))%*%(chol2inv(chol(F_alles[c(1:lin,(lin+(mi_ic-1)*nbasis+1):(lin+mi_ic*nbasis),(lin+m*nbasis+1):(lin+m*nbasis+n*s)),c(1:lin,(lin+(mi_ic-1)*nbasis+1):(lin+mi_ic*nbasis),(lin+m*nbasis+1):(lin+m*nbasis+n*s))]+K))*Nue_ma)%*%t(cbind(Z_r[[mi_ic]],W)*sqrt(D*1/Sigma*D*1/Sigma)))%*%(E-M0)
FinalHat[[1]]<-E-Uu

defre<-sum(diag(FinalHat[[1]]))

komp[1]<-mi_ic

X1<-cbind(X,Phi[[mi_ic]])

if(s==1)
{
PP<-rep((Q_start^(-1)),n*s)
PP<-diag(PP)
}else{
PP<-matrix(0,n*s,n*s)
for(jf in 1:n)
PP[((jf-1)*s+1):(jf*s),((jf-1)*s+1):(jf*s)]<-chol2inv(chol(Q_start))
}

k1a<-rep(0,lin.old)
if(diff.ord>1)
{
k2a<-c(rep(0,diff.ord-1),rep(1,ncol(Basis[[r]]$Z)))
}else{
k2a<-rep(1,ncol(Basis[[r]]$Z))
}
k22a<-rep(k2a,dim.U.fix+1)
KK<-lambda*diag(c(k1a,k22a))

R1<-chol2inv(chol((t(X1)%*%(X1*D*1/Sigma*D)+KK)-(t(X1)%*%(W*D*1/Sigma*D)%*%chol2inv(chol(t(W)%*%(W*D*1/Sigma*D)+PP))%*%t(t(X1)%*%(W*D*1/Sigma*D)))))
R2<-(t(X1*D*1/Sigma)-(t(X1)%*%(W*D*1/Sigma*D)%*%chol2inv(chol(t(W)%*%(W*D*1/Sigma*D)+PP)))%*%t(W*D*1/Sigma))

Q_mi<-control$nue*R1[(lin+1):(lin+nbasis),]%*%R2%*%(E-M0)

Eta<-Z_alles%*%Delta[1,]

Mu<-as.vector(family$linkinv(Eta))

Sigma<-as.vector(family$variance(Mu))

Standard_ma[[mi_ic]]<-Phi[[mi_ic]]%*%Q_mi
Standard[1,mi_ic,]<-sqrt(diag(Standard_ma[[mi_ic]]%*%t(Standard_ma[[mi_ic]]*Sigma)))

D<-as.vector(family$mu.eta(Eta))

if (control$method=="EM")
{
aktuell_vec<-!is.element(Delta[1,],0)
X_aktuell<-Z_alles[,aktuell_vec]

m_aktuell<- sum(is.element(1:m,komp))

if(s==1)
{
P_gross<-c(rep(0,lin+m_aktuell*nbasis),rep((Q_start^(-1)),n*s))
P_gross<-diag(P_gross)
}else{
P_gross<-matrix(0,lin+m_aktuell*nbasis+n*s,lin+m_aktuell*nbasis+n*s)
for(jf in 1:n)
P_gross[(lin+m_aktuell*nbasis+(jf-1)*s+1):(lin+m_aktuell*nbasis+jf*s),(lin+m_aktuell*nbasis+(jf-1)*s+1):(lin+m_aktuell*nbasis+jf*s)]<-chol2inv(chol(Q_start))
}

F_gross<-t(X_aktuell)%*%(X_aktuell*D*1/Sigma*D)+P_gross

k1<-rep(0,lin.old)
if(diff.ord>1)
{
k2<-c(rep(0,diff.ord-1),rep(1,ncol(Basis[[r]]$Z)))
}else{
k2<-rep(1,ncol(Basis[[r]]$Z))
}
k22<-rep(k2,dim.U.fix+m_aktuell)
k3<-rep(0,n*s)
K_gross<-lambda*diag(c(k1,k22,k3))

############################# update Q ################
InvFisher<-chol2inv(chol(F_gross+K_gross))

Q1<-InvFisher[(lin+m_aktuell+1):(lin+m_aktuell+s),(lin+m_aktuell+1):(lin+m_aktuell+s)]+Delta[1,(lin+m+1):(lin+m+s)]%*%t(Delta[1,(lin+m+1):(lin+m+s)])
for (i in 2:n)
Q1<-Q1+InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[1,(lin+m+(i-1)*s+1):(lin+m+i*s)]%*%t(Delta[1,(lin+m+(i-1)*s+1):(lin+m+i*s)])
Q1<-1/n*Q1
}else{

Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[1,1:(lin+m*nbasis)]

aktuell_vec<-!is.element(Delta[1,1:(lin+m*nbasis)],0)
X_aktuell<-Z_fastalles[,aktuell_vec]

if(s==1)
{
optim.obj<-nlminb(sqrt(Q_start),likelihood_bobyqa,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, model="no.model",  lower = 1e-14, upper = 20)
Q1<-as.matrix(optim.obj$par)^2
}else{
q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
upp<-rep(up1,length(q_start_vec))
low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
optim.obj<-try(bobyqa(q_start_vec,likelihood,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, model="no.model", lower=low,upper=upp))
Q1<-matrix(0,s,s)
Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
Q1<-Q1+t(Q1)
diag(Q1)<-(optim.obj$par[1:s])

#### Check for positive definitness ########
      for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
      Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
}}
    if (family$family=="poisson")
    IC_stop[1]<-sum(y*log(Mu)-Mu)

    if (family$family=="binomial")
    IC_stop[1]<-sum(y*log(Mu)+(1-y)*log(1-Mu))

    if (family$family=="gaussian")
    IC_stop[1]<-sum(y*Mu-0.5*(Mu^2))

   ############ AIC ##############
    if (control$sel.method=="aic")
    IC_stop[1]<-(-IC_stop[1])+defre

   ############ BIC ###############
    if (control$sel.method=="bic")
    IC_stop[1]<-(-IC_stop[1])+0.5*defre*log(n)

if(control$overdispersion==T)
{
phi[1]<-(sum((y-Mu)^2/Mu))/(N-defre)
Sigma<-Sigma*phi[1]
}

Q[,,2]<-Q1
###############################################################################################################################################
################################################################### Boost ###################################################################
if(control$steps!=1)
{
for (l in 2:control$steps)
{
print(paste("Iteration ", l,sep=""))
if(s==1)
{
P_alles<-c(rep(0,lin+m*nbasis),rep((Q1^(-1)),n*s))
P_alles<-diag(P_alles)
}else{
P_alles<-matrix(0,lin+m*nbasis+n*s,lin+m*nbasis+n*s)
for(jf in 1:n)
P_alles[(lin+m*nbasis+(jf-1)*s+1):(lin+m*nbasis+jf*s),(lin+m*nbasis+(jf-1)*s+1):(lin+m*nbasis+jf*s)]<-chol2inv(chol(Q1))
}

score_alles<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
F_alles<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P_alles

  for (r in 1:m)
  {
  if(r==1)
  Inv_r<-chol2inv(chol(F_alles[c(1:lin,(lin+(r-1)*nbasis+1):(lin+r*nbasis),(lin+m*nbasis+1):(lin+m*nbasis+n*s)),c(1:lin,(lin+(r-1)*nbasis+1):(lin+r*nbasis),(lin+m*nbasis+1):(lin+m*nbasis+n*s))]+K))
  Delta_r[,r]<-Inv_r%*%score_alles[c(1:lin,(lin+(r-1)*nbasis+1):(lin+r*nbasis),(lin+m*nbasis+1):(lin+m*nbasis+n*s))]

  Eta_r[,r]<-Eta_start+Z_r_gross[[r]]%*%Delta_r[,r]
  Mu_r[,r]<-as.vector(family$linkinv(Eta_r[,r]))

  M1[[r]]<-(Z_r_gross[[r]]*sqrt(Sigma*D*1/Sigma*D))%*%Inv_r%*%t(Z_r_gross[[r]]*sqrt(D*1/Sigma*D*1/Sigma))
  diag_H1[r]<-sum(diag(E-((E-M1[[r]])%*%Uu)))

   ######## likelihood of information criterion    ########
      if (family$family=="poisson")
      IC[l,r]<-sum(y*log(Mu_r[,r])-Mu_r[,r])#-0.5*pena

      if (family$family=="binomial")
      IC[l,r]<-sum(y*log(Mu_r[,r])+(1-y)*log(1-Mu_r[,r]))#-0.5*pena

      if (family$family=="gaussian")
      IC[l,r]<-sum(y*Mu_r[,r]-0.5*(Mu_r[,r]^2))#-0.5*pena

   ############ AIC ##############
      if (control$sel.method=="aic")
      IC[l,r]<-(-IC[l,r])+diag_H1[r]

   ############ BIC ###############
      if (control$sel.method=="bic")
      IC[l,r]<-(-IC[l,r])+0.5*diag_H1[r]*log(n)
   }

   mi_ic<-match(min(IC[l,]),IC[l,])

   g<-c(rep(T,lin),grp==mi_ic,rep(T,s*n))

   Delta[l,]<-Delta[l-1,]
   Delta_r[(lin+1):(lin+nbasis),mi_ic]<-control$nue*Delta_r[(lin+1):(lin+nbasis),mi_ic]
   Delta[l,g]<-Delta[l,g]+t(Delta_r[,mi_ic])

komp[l]<-mi_ic

Uu<-(E-(cbind(Z_r[[mi_ic]],W)*sqrt(Sigma*D*1/Sigma*D))%*%(chol2inv(chol(F_alles[c(1:lin,(lin+(mi_ic-1)*nbasis+1):(lin+mi_ic*nbasis),(lin+m*nbasis+1):(lin+m*nbasis+n*s)),c(1:lin,(lin+(mi_ic-1)*nbasis+1):(lin+mi_ic*nbasis),(lin+m*nbasis+1):(lin+m*nbasis+n*s))]+K))*Nue_ma)%*%t(cbind(Z_r[[mi_ic]],W)*sqrt(D*1/Sigma*D*1/Sigma)))%*%(E-M0)
FinalHat[[l]]<-E-Uu
defre<-sum(diag(FinalHat[[l]]))

X1<-cbind(X,Phi[[mi_ic]])

if(s==1)
{
PP<-rep((Q1^(-1)),n*s)
PP<-diag(PP)
}else{
PP<-matrix(0,n*s,n*s)
for(jf in 1:n)
PP[((jf-1)*s+1):(jf*s),((jf-1)*s+1):(jf*s)]<-chol2inv(chol(Q1))
}

R1<-chol2inv(chol((t(X1)%*%(X1*D*1/Sigma*D)+KK)-(t(X1)%*%(W*D*1/Sigma*D)%*%chol2inv(chol(t(W)%*%(W*D*1/Sigma*D)+PP))%*%t(t(X1)%*%(W*D*1/Sigma*D)))))
R2<-(t(X1*D*1/Sigma)-(t(X1)%*%(W*D*1/Sigma*D)%*%chol2inv(chol(t(W)%*%(W*D*1/Sigma*D)+PP)))%*%t(W*D*1/Sigma))

Q_mi<-control$nue*R1[(lin+1):(lin+nbasis),]%*%R2%*%(E-FinalHat[[l-1]])

    if(is.matrix(Standard_ma[[mi_ic]]))
    {
    Standard_ma[[mi_ic]]<-Standard_ma[[mi_ic]]+Phi[[mi_ic]]%*%Q_mi
    }else{
    Standard_ma[[mi_ic]]<-Phi[[mi_ic]]%*%Q_mi
    }

Eta<-Z_alles%*%Delta[l,]

Mu<-as.vector(family$linkinv(Eta))

Sigma<-as.vector(family$variance(Mu))

    for(jh in as.numeric(levels(as.factor(komp[1:l])))[as.numeric(levels(as.factor(komp[1:l])))!=0])
    Standard[l,jh,]<-sqrt(diag(Standard_ma[[jh]]%*%t(Standard_ma[[jh]]*Sigma)))

D<-as.vector(family$mu.eta(Eta))

if (control$method=="EM")
{
aktuell_vec<-!is.element(Delta[l,],0)
X_aktuell<-Z_alles[,aktuell_vec]

m_aktuell<- sum(is.element(1:m,komp))

if(s==1)
{
P_gross<-c(rep(0,lin+m_aktuell*nbasis),rep((Q1^(-1)),n*s))
P_gross<-diag(P_gross)
}else{
P_gross<-matrix(0,lin+m_aktuell*nbasis+n*s,lin+m_aktuell*nbasis+n*s)
for(jf in 1:n)
P_gross[(lin+m_aktuell*nbasis+(jf-1)*s+1):(lin+m_aktuell*nbasis+jf*s),(lin+m_aktuell*nbasis+(jf-1)*s+1):(lin+m_aktuell*nbasis+jf*s)]<-chol2inv(chol(Q1))
}

F_gross<-t(X_aktuell)%*%(X_aktuell*D*1/Sigma*D)+P_gross

k22<-rep(k2,dim.U.fix+m_aktuell)
k3<-rep(0,n*s)
K_gross<-lambda*diag(c(k1,k22,k3))

############################# update Q ################
InvFisher<-chol2inv(chol(F_gross+K_gross))

Q1<-InvFisher[(lin+m_aktuell+1):(lin+m_aktuell+s),(lin+m_aktuell+1):(lin+m_aktuell+s)]+Delta[l,(lin+m+1):(lin+m+s)]%*%t(Delta[l,(lin+m+1):(lin+m+s)])

for (i in 2:n)
Q1<-Q1+InvFisher[(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s),(lin+m_aktuell+(i-1)*s+1):(lin+m_aktuell+i*s)]+Delta[l,(lin+m+(i-1)*s+1):(lin+m+i*s)]%*%t(Delta[l,(lin+m+(i-1)*s+1):(lin+m+i*s)])
Q1<-1/n*Q1
}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[l,1:(lin+m*nbasis)]

aktuell_vec<-!is.element(Delta[l,1:(lin+m*nbasis)],0)
X_aktuell<-Z_fastalles[,aktuell_vec]

if(s==1)
{
optim.obj<-try(nlminb(sqrt(Q1),likelihood_bobyqa,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, model="no.model",  lower = 1e-12, upper = 20))
Q1<-as.matrix(optim.obj$par)^2
}else{
Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
optim.obj<-try(bobyqa(Q1_vec,likelihood,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, model="no.model", lower=low,upper=upp))

Q1<-matrix(0,s,s)
Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
Q1<-Q1+t(Q1)
diag(Q1)<-(optim.obj$par[1:s])

#### Check for positiv definitness ########
for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
       Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
}}
    if (family$family=="poisson")
    IC_stop[l]<-sum(y*log(Mu)-Mu)

    if (family$family=="binomial")
    IC_stop[l]<-sum(y*log(Mu)+(1-y)*log(1-Mu))

    if (family$family=="gaussian")
    IC_stop[l]<-sum(y*Mu-0.5*(Mu^2))

   ############ AIC ##############
    if (control$sel.method=="aic")
    IC_stop[l]<-(-IC_stop[l])+defre

   ############ BIC ###############
    if (control$sel.method=="bic")
    IC_stop[l]<-(-IC_stop[l])+0.5*defre*log(n)

if(control$overdispersion==T)
{
phi[l]<-(sum((y-Mu)^2/Mu))/(N-defre)
Sigma<-Sigma*phi[l]
}
Q[,,l+1]<-Q1

kritval<-sqrt(sum((Delta[l-1,]-Delta[l,])^2))/sqrt(sum(Delta[l,]^2))
if(kritval<1e-8)
break

kritval2<-abs(IC_stop[l-1]-IC_stop[l])/abs(IC_stop[l])
if(kritval2<1e-8)
break
}}
if (l<control$steps)
{
for (wj in (l+1):control$steps)
{
Delta[wj,]<-Delta[l,]
Q[,,wj+1]<-Q1
IC[wj,]<-IC[l,]
IC_stop[wj]<-IC_stop[l]
}}

##########
opt<-match(min(IC_stop),IC_stop)

if (control$OPT==TRUE)
{
Delta_neu<-Delta[opt,]
}else{
Delta_neu<-Delta[l,]
}

Eta_opt<-Z_alles%*%Delta_neu
Mu_opt<-as.vector(family$linkinv(Eta_opt))
Sigma_opt<-as.vector(family$variance(Mu_opt))
D_opt<-as.vector(family$mu.eta(Eta_opt))

aktuell_vec<-!is.element(Delta_neu,0)
X_aktuell<-Z_alles[,aktuell_vec]
m_aktuell<- sum(is.element(1:m,komp[1:opt]))

k1b<-rep(0,lin.old)
if(diff.ord>1)
{
k2b<-c(rep(0,diff.ord-1),rep(1,ncol(Basis[[r]]$Z)))
}else{
k2b<-rep(1,ncol(Basis[[r]]$Z))
}
k22b<-rep(k2b,dim.U.fix+m_aktuell)
k3b<-rep(0,n*s)
K_opt<-lambda*diag(c(k1b,k22b,k3b))

if(s==1)
{
P_opt<-c(rep(0,lin+m_aktuell*nbasis),rep((Q[,,opt+1]^(-1)),n*s))
P_opt<-diag(P_opt)
}else{
P_opt<-matrix(0,lin+m_aktuell*nbasis+n*s,lin+m_aktuell*nbasis+n*s)
for(jf in 1:n)
P_opt[(lin+m_aktuell*nbasis+(jf-1)*s+1):(lin+m_aktuell*nbasis+jf*s),(lin+m_aktuell*nbasis+(jf-1)*s+1):(lin+m_aktuell*nbasis+jf*s)]<-chol2inv(chol(Q[,,opt+1]))
}

F_opt<-t(X_aktuell)%*%(X_aktuell*D_opt*1/Sigma_opt*D_opt)+P_opt

Inv_F_opt<-chol2inv(chol(F_opt+K_opt))
Lin_errors<-sqrt(diag(Inv_F_opt))

Standard_errors<-list()
Standard_errors$smooth<-Standard[opt,,]
rownames(Standard_errors$smooth)<-U.names
if(!is.null(U.fix))
{
Standard_errors$smooth.fix<-matrix(0,dim.U.fix,N)
for (uz in 1:dim.U.fix)
Standard_errors$smooth.fix[uz,]<-sqrt(diag(Phi.fix[[uz]]%*%Inv_F_opt[(lin.old+1):lin,(lin.old+1):lin]%*%t(Phi.fix[[uz]])))

Standard_errors$smooth<-rbind(Standard_errors$smooth.fix,Standard_errors$smooth)
rownames(Standard_errors$smooth)[1:dim.U.fix]<-colnames(U.fix)
}
Standard_errors$fixef<-Lin_errors[1:lin.old]
Standard_errors$ranef<-Lin_errors[(lin+m_aktuell*nbasis+1):(lin+m_aktuell*nbasis+n*s)]

y_dach<-as.vector(family$linkinv(Eta_opt))

if(s==1)
{
Qfinal<-sqrt(Q[,,opt+1])
}else{
Qfinal<-Q[,,opt+1]
}

if(!is.matrix(Qfinal))
Qfinal<-as.matrix(Qfinal)
colnames(Qfinal)<-random.labels
rownames(Qfinal)<-random.labels

names(Delta_neu)[1:lin.old]<-colnames(X)[1:lin.old]

if(!is.null(U.fix))
{
U.fix.names<-character()
for (rr in 1:dim.U.fix)
U.fix.names<-c(U.fix.names,colnames(Phi.fix[[rr]]))
names(Delta_neu)[(lin.old+1):lin]<-U.fix.names
}

U.long.names<-character()
for (r in 1:m)
U.long.names<-c(U.long.names,colnames(Phi[[r]]))
names(Delta_neu)[(lin+1):(lin+m*nbasis)]<-U.long.names

names(Delta_neu)[(lin+m*nbasis+1):(lin+m*nbasis+n*s)]<-colnames(W)

comp<-character()
for(oi in 1:length(komp))
comp<-c(comp,U.names[komp[oi]])

ret.obj=list()
ret.obj$IC<-IC
ret.obj$IC_sel<-IC_stop
ret.obj$opt<-opt
ret.obj$Deltamatrix<-Delta
ret.obj$ranef<-Delta_neu[(lin+m*nbasis+1):(lin+m*nbasis+n*s)]
ret.obj$coefficients<-Delta_neu[1:(lin.old)]
ret.obj$spline.weights<-Delta_neu[(lin.old+1):(lin+m*nbasis)]
ret.obj$fixerror<-Standard_errors$fixef
ret.obj$ranerror<-Standard_errors$ranef
ret.obj$smootherror<-Standard_errors$smooth
ret.obj$lambda<-lambda
ret.obj$components<-comp
ret.obj$Q_long<-Q
ret.obj$Q<-Qfinal
ret.obj$y_hat<-Mu_opt
ret.obj$HatMatrix<-FinalHat[[opt]]
ret.obj$phi<-phi[opt]
ret.obj$family<-family
ret.obj$fix<-fix
ret.obj$add<-add
ret.obj$newrndfrml<-newrndfrml
ret.obj$subject<-names(rnd)[1]
ret.obj$data<-data
ret.obj$U<-cbind(U.fix,U)
ret.obj$U.fix.names<-colnames(U.fix)
ret.obj$nbasis<-nbasis
ret.obj$spline.degree<-spline.degree
ret.obj$diff.ord<-diff.ord
ret.obj$factor.names<-factor.names
ret.obj$factor.list<-factor.list
return(ret.obj)
}



bGAMM <- function(fix=formula, add=formula, rnd=formula, data, lambda, family=NULL, control=list()) UseMethod("bGAMM")

bGAMM.formula <- function(fix,add,rnd,...)
{
est <- est.bGAMM(fix,add,rnd,...)
est$fitted.values <- est$y_hat
est$StdDev <- est$Q
est$call <- match.call()
class(est) <- "bGAMM"
est
}


print.bGAMM <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nFixed Effects:\n")
cat("\nCoefficients:\n")
print(x$coefficients)
cat("\nForced Smooth Effects:\n")
print(x$U.fix.names)
cat("\nSelected Smooth Effects:\n")
print(levels(as.factor(x$comp)))
cat("\nRandom Effects:\n")
cat("\nStdDev:\n")
print(x$StdDev)
}




summary.bGAMM <- function(object, ...)
{
se <- object$fixerror
zval <- coefficients(object) / se
TAB <- cbind(Estimate = coefficients(object),
StdErr = se,
z.value = zval,
p.value = 2*pnorm(-abs(zval)))
res <- list(call=object$call,
coefficients=TAB,StdDev=object$StdDev, comp=object$comp,U.fix.names=object$U.fix.names)
class(res) <- "summary.bGAMM"
res
}


print.summary.bGAMM <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\n")
cat("\nFixed Effects:\n")
cat("\nCoefficients:\n")
printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
cat("\nForced Smooth Effects:\n")
print(x$U.fix.names)
cat("\nSelected Smooth Effects:\n")
print(levels(as.factor(x$comp)))
cat("\nRandom Effects:\n")
cat("\nStdDev:\n")
print(x$StdDev)
}


predict.bGAMM <- function(object,newdata=NULL,...)
{
if(is.null(newdata))
{
y<-fitted(object)
}else{
nbasis<-object$nbasis
diff.ord<-object$diff.ord
spline.degree<-object$spline.degree
knots.no<-nbasis-(spline.degree-2)

if(length(object$factor.names>0))
{
for (i in 1:length(object$factor.names))
{
newdata[,object$factor.names[i]]<-factor(newdata[,object$factor.names[i]],levels=object$factor.list[[i]])
}}

X <- model.matrix(object$fix, newdata)

U <- model.matrix(object$add, newdata)
U<-U[,-1]
if(length(object$U.fix.names)>0)
{
reord<-is.element(colnames(U),object$U.fix.names)
U<-cbind(U[,reord],U[,!reord])
}

m<-dim(U)[2]

Basis<-list()
Phi<-numeric()

for (r in 1:m)
{
Basis[[r]]<-bs.design(U[,r],diff.ord=diff.ord,spline.degree=spline.degree,knots.no=knots.no)
Phi<-cbind(Phi,Basis[[r]]$X[,-1],Basis[[r]]$Z)
}

subj.new<-levels(as.factor(newdata[,object$subject]))
subj.old<-levels(as.factor(object$data[,object$subject]))
subj.test<-is.element(subj.new,subj.old)
subj.ok<-subj.new[subj.test]

krit.random<-!all(!is.element(subj.new,subj.old))
if(krit.random)
{
W_start <- model.matrix(formula(object$newrndfrml), newdata)
}else{
W_start <- NULL
}

rnlabels<-terms(formula(object$newrndfrml))
random.labels<-attr(rnlabels,"term.labels")
s<-length(random.labels)
k<-table(newdata[,colnames(newdata)==(object$subject)])
n<-length(k)

if(s>1)
{for (i in 2:s)
subj.test<-cbind(subj.test,subj.test)
subj.test<-as.vector(t(subj.test))
}

if(krit.random)
{
if(s>1)
{
W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
for (i in 2:n)
W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
}else{
W<-W_start
}
y<- as.vector(object$family$linkinv(cbind(X,Phi)%*%c(object$coef,object$spline.weights)))
rand.ok<-is.element(newdata[,object$subject],subj.ok)
W.neu<-W[,subj.test]
y[rand.ok]<- object$family$linkinv(cbind(X,Phi,W.neu)[rand.ok,]%*%c(object$coef,object$spline.weights,object$ranef[match(colnames(W.neu),names(object$ranef))]))
#y<- as.vector(object$family$linkinv(cbind(X,Phi,W)%*%c(object$coef,object$spline.weights,object$ranef)))
}else{
W<-NULL
y<- as.vector(object$family$linkinv(cbind(X,Phi)%*%c(object$coef,object$spline.weights)))
}
}
y
}




plot.bGAMM <- function(x,which=NULL,...)
{
U<-x$U
m<-dim(U)[2]

if(is.null(which))
which<-1:m

p<-length(which)
if(p>9)
stop("Too many smooth functions! Please specify at maximum nine.")

a<-ceiling(sqrt(p))
b<-floor(sqrt(p))
if(b==0)
b<-1

nbasis<-x$nbasis
diff.ord<-x$diff.ord
spline.degree<-x$spline.degree
knots.no<-nbasis-(spline.degree-2)


spline.ma<-list()
Design<-list()
smooth.ma<-matrix(0,m,dim(U)[1])
smootherror<-x$smootherror

for(i in which)
{
spline.ma[[i]]<-bs.design(sort(U[,i]), diff.ord=diff.ord, spline.degree=spline.degree, knots.no=knots.no)
smootherror[i,]<-x$smootherror[i,order(U[,i])]
Design[[i]]<-cbind(spline.ma[[i]]$X[,-1],spline.ma[[i]]$Z)
smooth.ma[i,]<-Design[[i]]%*%x$spline.weights[((i-1)*nbasis+1):(i*nbasis)]
}

smoothlow<-smooth.ma-2*smootherror
smoothupp<-smooth.ma+2*smootherror

par(mfrow=c(a,b))

for(i in which)
{
plot(sort(U[,i]), smooth.ma[i,], type = "l", lwd=2, xlab=paste(colnames(U)[i]),ylab="",main=" ",cex.lab=2,cex.axis=2,ylim=c(min(smoothlow[i,])-0.1*(max(smoothlow[i,])-min(smoothlow[i,])),max(smoothupp[i,])+0.1*(max(smoothlow[i,])-min(smoothlow[i,]))))
lines(sort(U[,i]),smoothlow[i,], type = "l",lty=2, lwd=2)
lines(sort(U[,i]),smoothupp[i,], type = "l",lty=2, lwd=2)
rug(jitter(U[,i]))
}
}

