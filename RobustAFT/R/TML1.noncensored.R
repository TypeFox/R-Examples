"TML1.noncensored" <-
function(y,errors= c("Gaussian", "logWeibull"), cu=NULL, initial=c("S","input"),otp=c("adaptive","fixed"),
               cov=c("no","parametric","nonparametric"),input=NULL, control = list(), ...)
#               iv=1,ctrg=NULL)
{
	control <- do.call("TML1.noncensored.control", control)
	iv <- control$iv
	gam <- control$gam
        maxit <- control$maxit
        tol <- control$tol
        n <- length(y)
 if (initial=="S")     {old <- comval(); dfcomn2(ipsi=4, xk=1.5477)
                        Beta0 <- integrate(Chiphi, -10, 10)$value
                        dfcomn2(ipsi=old$ipsi, xk=old$xk)
                        ctrS <- c(TML1.noncensored.control.S(...), k0=1.5477, k1=4.6873, Beta0=Beta0)
                       }
if(errors == "Gaussian"){
  if (is.null(cu)) cu <- 2.5
  cl <- -cu
 # Step 1: Initial high bdp estimate
 if (initial=="S")     {z   <- MM.E.gauss(y,onlyS=TRUE,control=ctrS)
                        th0 <- z$lambda; v0 <- z$sigma; nit0 <- z$nit.ref}
 if (initial=="input") {z   <- input
                       th0 <- z$theta; v0 <- z$sigma; nit0 <- 0}
 nares <- list(th0=th0,v0=v0,th1=NA,v1=NA,tl=NA,tu=NA,alpha=NA,beta=NA)
 # Step 2: rejection rule
 yo  <- sort(y); re  <- yo-th0; rs <- re/v0
 tp  <- adaptn(rs,cl,cu,option=otp); if (is.na(tp$tu)) return(nares)
 wi  <- tPsin(rs,tp$tl,tp$tu); yr <- yo[wi!=0 | rs==0]; tp$tn <- length(yr)
 # Step 3: ML-estimate on retained observations

 z   <- MLn1(yr,iv,tp)
 res <- list(th0=th0,v0=v0,th1=z$th1,v1=z$v1,tl=tp$tl,tu=tp$tu,
        alpha=tp$alpha,tn=tp$tn,beta=tp$beta,wi=(wi!=0)*1)
 if (cov!="no") {l <- cl; u <- cu; if (otp=="adaptive") {l <- tp$tl; u <- tp$tu} 
 # if (cov=="halfparametric") K <- Cov2.n1(y,u,z$th1,z$v1,opt="integrals")
   if (cov=="nonparametric")  K <- Cov2.n1(y,u,z$th1,z$v1,opt="averages")
   if (cov=="parametric"   )  K <- CovE.n1(y,u,z$th1,z$v1)
 res  <- c(res,list(CV0=K$CV0,CV1=K$CV1))}
}

if(errors == "logWeibull"){
  if (is.null(cu)) cu <- 1.855356
        cl <- Izero(cu)

# Step 1: Initial high bdp estimate
if (initial=="S")     {z  <- MM.E.gauss(y,control=ctrS,onlyS=TRUE); b0 <- -0.1352
                       v0 <- z$sigma; th0 <- z$lambda-b0*v0; nit0 <- z$nit.ref}
if (initial=="input") {z  <- input;   v0  <- z$v; th0 <- z$tau; nit0 <- 0}
# Step 2: rejection rule
nares <- list(th0=th0,v0=v0,nit0=NA,th1=NA,v1=NA,nit1=NA,tl=NA,tu=NA,alpha=NA,beta=NA)
yo <- sort(y); re  <- yo-th0; rs <- re/v0
tp <- adaptw(rs,cl,cu,otp); if (is.na(tp$tu)) return(nares)
wi <- tPsiw(rs,tp$tl,tp$tu); yr <- yo[wi!=0 | rs==0]; tp$tn <- length(yr)
# Step 3: ML-estimate on retained observations
z  <- MLw1(yr,th0,v0,iv,n,tp,gam,maxit,tol)
yinf<-exp(z$tl*z$v0+z$th0)
ysup<-exp(z$tu*z$v0+z$th0)
res <- list(th0=th0,v0=v0,nit0=nit0,th1=z$th1,v1=z$v1,nit1=z$nit,tl=tp$tl,tu=tp$tu,
       alpha=tp$alpha,tn=tp$tn,beta=tp$beta,yinf=yinf,ysup=ysup,wi=(wi!=0)*1)
if (cov!="no") {l <- cl; u <- cu; if (otp=="adaptive") {l <- tp$tl; u <- tp$tu} 
# if (cov=="halfparametric") K <- Cov2.w1(y,l,u,z$th1,z$v1,opt="integrals")
  if (cov=="nonparametric")  K <- Cov2.w1(y,l,u,z$th1,z$v1,opt="averages")
  if (cov=="parametric"   )  K <- CovE.w1(y,l,u,z$th1,z$v1)
  res  <- c(res,list(CV0=K$CV0,CV1=K$CV1))}
}
res}

