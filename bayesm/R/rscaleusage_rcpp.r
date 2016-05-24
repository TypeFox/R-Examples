rscaleUsage=
function(Data,Prior,Mcmc) 
{
#
# purpose: run scale-usage mcmc
#    draws y,Sigma,mu,tau,sigma,Lambda,e
#                                R. McCulloch 12/28/04
#    added classes 3/07
# 
# arguments:
#    Data:
#     all components are required:
#       k:  integer giving the scale of the responses, each observation is an integer from 1,2,...k
#       x:  data, num rows=number of respondents, num columns = number of questions
#    Prior:
#     all components are optional
#       nu,V: Sigma ~ IW(nu,V)
#       mubar,Am: mu ~N(mubar,Am^{-1})
#       gsigma: grid for sigma
#       gl11,gl22,gl12: grids for ij element of Lambda
#       Lambdanu,LambdaV: Lambda ~ IW(Lambdanu,LambdaV)
#       ge: grid for e
#    Mcmc:
#     all components are optional (but you would typically want to specify R= number of draws)
#       R: number of mcmc iterations
#       keep: frequency with which draw is kept
#       ndghk: number of draws for ghk
#       nprint - print estimated time remaining on every nprint'th draw
#       e,y,mu,Sigma,sigma,tau,Lambda: initial values for the state
#       doe, ...doLambda: indicates whether draw should be made
# output:
#    List with draws of each of Sigma,mu,tau,sigma,Lambda,e
#    eg. result$Sigma is the draws of Sigma
#    Each component is a matrix expept e, which is a vector
#    for the matrices Sigma and Lambda each row transpose of the Vec
#    eg. result$Lambda has rows (Lambda11,Lambda21,Lambda12,Lambda22)

#
# define functions needed
#
# -----------------------------------------------------------------------------------
myin = function(i,ind) {i %in% ind}

ispd = function(mat,d=nrow(mat)) {
if(!is.matrix(mat)) {
res = FALSE
} else if(!((nrow(mat)==d) & (ncol(mat)==d))) {
res = FALSE
} else {
diff = (t(mat)+mat)/2 - mat
perdiff = sum(diff^2)/sum(mat^2)
res = ((det(mat)>0) & (perdiff < 1e-10))
}
res
}
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# print out components of inputs ----------------------------------------------
cat('\nIn function rscaleUsage\n\n')
if(!missing(Data)) {
cat('   Data has components: ')
cat(paste(names(Data),collapse=' ')[1],'\n')
}
if(!missing(Prior)) {
cat('   Prior has components: ')
cat(paste(names(Prior),collapse=' ')[1],'\n')
}
if(!missing(Mcmc)) {
cat('   Mcmc has components: ')
cat(paste(names(Mcmc),collapse=' ')[1],'\n')
}
cat('\n')
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# process Data argument --------------------------
if(missing(Data)) {pandterm("Requires Data argument - list of k=question scale and x = data")}
if(is.null(Data$k)) {
   pandterm("k not specified")
} else {
   k = as.integer(Data$k)
   if(!((k>0) & (k<50))) {pandterm("Data$k must be integer between 1 and 50")}
}
if(is.null(Data$x)) {
   pandterm('x (the data), not specified')
} else {
   if(!is.matrix(Data$x)) {pandterm('Data$x must be a matrix')}
   x = matrix(as.integer(Data$x),nrow=nrow(Data$x))
   checkx = sum(sapply(as.vector(x),myin,1:k))
   if(!(checkx == nrow(x)*ncol(x))) {pandterm('each element of Data$x must be in 1,2...k')}
   p = ncol(x)
   n = nrow(x)
   if((p<2) | (n<1)) {pandterm(paste('invalid dimensions for x: nrow,ncol: ',n,p))}
}
# ++++++++++++++++++++++++++++++++++++++++++++++++

# process Mcmc argument ---------------------

#run mcmc
R = 1000
keep = BayesmConstant.keep
ndghk= 100
nprint = BayesmConstant.nprint
if(!missing(Mcmc)) {
if(!is.null(Mcmc$R))              { R = as.integer(Mcmc$R) }
if(!is.null(Mcmc$keep))           { keep = as.integer(Mcmc$keep) }
if(!is.null(Mcmc$ndghk))          { ndghk = as.integer(Mcmc$ndghk) }
if(!is.null(Mcmc$nprint))         { nprint = as.integer(Mcmc$nprint) }
}
if(R<1) { pandterm('R must be positive')}
if(keep<1) { pandterm('keep must be positive') }
if(ndghk<1) { pandterm('ndghk must be positive') }
if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}

#state
y = matrix(as.double(x),nrow=nrow(x))
mu = apply(y,2,mean)
Sigma = var(y)
tau = rep(0,n)
sigma = rep(1,n)
#Lambda = matrix(c(3.7,-.22,-.22,.32),ncol=2)
#Lambda = matrix(c((k/4)^2,(k/4)*.5*(-.2),0,.25),nrow=2); Lambda[1,2]=Lambda[2,1]
Lambda = matrix(c(4,0,0,.5),ncol=2)
e=0
if(!missing(Mcmc)) {
if(!is.null(Mcmc$y))         { y = Mcmc$y }
if(!is.null(Mcmc$mu))        { mu = Mcmc$mu }
if(!is.null(Mcmc$Sigma))     { Sigma = Mcmc$Sigma }
if(!is.null(Mcmc$tau))       { tau = Mcmc$tau }
if(!is.null(Mcmc$sigma))     { sigma = Mcmc$sigma }
if(!is.null(Mcmc$Lambda))    { Lambda = Mcmc$Lambda }
if(!is.null(Mcmc$e))         { e = Mcmc$e }
}
if(!ispd(Sigma,p)) { pandterm(paste('Sigma must be positive definite with dimension ',p)) }
if(!ispd(Lambda,2)) { pandterm(paste('Lambda must be positive definite with dimension ',2)) }
if(!is.vector(mu)) { pandterm('mu must be a vector') }
if(length(mu) != p) { pandterm(paste('mu must have length ',p)) }
if(length(tau) != n) { pandterm(paste('tau must have length ',n)) }
if(!is.vector(sigma)) { pandterm('sigma must be a vector') }
if(length(sigma) != n) { pandterm(paste('sigma must have length ',n)) }
if(!is.matrix(y)) { pandterm('y must be a matrix') }
if(nrow(y) != n) { pandterm(paste('y must have',n,'rows')) }
if(ncol(y) != p) { pandterm(paste('y must have',p,'columns')) }

#do draws
domu=TRUE
doSigma=TRUE
dosigma=TRUE
dotau=TRUE
doLambda=TRUE
doe=TRUE
if(!missing(Mcmc)) {
if(!is.null(Mcmc$domu))        { domu = Mcmc$domu }
if(!is.null(Mcmc$doSigma))     { doSigma = Mcmc$doSigma }
if(!is.null(Mcmc$dotau))       { dotau = Mcmc$dotau }
if(!is.null(Mcmc$dosigma))     { dosigma = Mcmc$dosigma }
if(!is.null(Mcmc$doLambda))    { doLambda = Mcmc$doLambda }
if(!is.null(Mcmc$doe))         { doe = Mcmc$doe }
}


#++++++++++++++++++++++++++++++++++++++

#process Prior argument ----------------------------------
nu = p+BayesmConstant.nuInc
V= nu*diag(p)
mubar = matrix(rep(k/2,p),ncol=1)
Am = BayesmConstant.A*diag(p)
gs = 200
gsigma = 6*(1:gs)/gs
gl11 = .1 + 5.9*(1:gs)/gs
gl22 = .1 + 2.0*(1:gs)/gs
#gl12 = -.8 + 1.6*(1:gs)/gs
gl12 = -2.0 + 4*(1:gs)/gs
nuL=20
VL = (nuL-3)*Lambda
ge = -.1+.2*(0:gs)/gs

if(!missing(Prior)) {
if(!is.null(Prior$nu))       { nu = Prior$nu; V = nu*diag(p) }
if(!is.null(Prior$V))        { V = Prior$V }
if(!is.null(Prior$mubar))    { mubar = matrix(Prior$mubar,ncol=1) }
if(!is.null(Prior$Am))       { Am = Prior$Am }
if(!is.null(Prior$gsigma))   { gsigma = Prior$gsigma }
if(!is.null(Prior$gl11))     { gl11 = Prior$gl11 }
if(!is.null(Prior$gl22))     { gl22 = Prior$gl22 }
if(!is.null(Prior$gl12))     { gl12 = Prior$gl12 }
if(!is.null(Prior$Lambdanu)) { nuL = Prior$Lambdanu; VL = (nuL-3)*Lambda }
if(!is.null(Prior$LambdaV))  { VL = Prior$LambdaV }
if(!is.null(Prior$ge))       { ge = Prior$ge }
}
if(!ispd(V,p)) { pandterm(paste('V must be positive definite with dimension ',p)) }
if(!ispd(Am,p)) { pandterm(paste('Am must be positive definite with dimension ',p)) }
if(!ispd(VL,2)) { pandterm(paste('VL must be positive definite with dimension ',2)) }
if(nrow(mubar) != p) { pandterm(paste('mubar must have length',p)) }
#++++++++++++++++++++++++++++++++++++++++

#print out run info -------------------------
#
# note in the documentation and in BSM, m is used instead of p
#    for print-out purposes I'm using m   P. Rossi 12/06
cat('   n,m,k: ', n,p,k,'\n')
cat('   R,keep,ndghk,nprint: ', R,keep,ndghk,nprint,'\n')
cat('\n')
cat('   Data:\n')
cat('      x[1,1],x[n,1],x[1,m],x[n,m]: ',x[1,1],x[n,1],x[1,p],x[n,p],'\n\n')
cat('   Prior:\n')
cat('      ','nu: ',nu,'\n')
cat('      ','V[1,1]/nu,V[m,m]/nu: ',V[1,1]/nu,V[p,p]/nu,'\n')
cat('      ','mubar[1],mubar[m]: ',mubar[1],mubar[p],'\n')
cat('      ','Am[1,1],Am[m,m]: ',Am[1,1],Am[p,p],'\n')
cat('      ','Lambdanu: ',nuL,'\n')
cat('      ','LambdaV11,22/(Lambdanu-3): ',VL[1,1]/(nuL-3),VL[2,2]/(nuL-3),'\n')
cat('      ','sigma grid, 1,',length(gsigma),': ',gsigma[1],', ',gsigma[length(gsigma)],'\n')
cat('      ','Lambda11 grid, 1,',length(gl11),': ',gl11[1],', ',gl11[length(gl11)],'\n')
cat('      ','Lambda12 grid, 1,',length(gl12),': ',gl12[1],', ',gl12[length(gl12)],'\n')
cat('      ','Lambda22 grid, 1,',length(gl22),': ',gl22[1],', ',gl22[length(gl22)],'\n')
cat('      ','e grid, 1,',length(ge),': ',ge[1],', ',ge[length(ge)],'\n')
cat('      ','draw e: ',doe,'\n')
cat('      ','draw Lambda: ',doLambda,'\n')
#++++++++++++++++++++++++++++++++++++++++++++

###################################################################
# Wayne Taylor
# 3/14/2015
###################################################################
out = rscaleUsage_rcpp_loop(k,x,p,n,
                            R,keep,ndghk,nprint,
                            y,mu,Sigma,tau,sigma,Lambda,e,
                            domu,doSigma,dosigma,dotau,doLambda,doe,
                            nu,V,mubar,Am, 
                            gsigma,gl11,gl22,gl12,
                            nuL,VL,ge)

R = out$ndpost
###################################################################
attributes(out$drmu)$class=c("bayesm.mat","mcmc")
attributes(out$drmu)$mcpar=c(1,R,keep)
attributes(out$drtau)$class=c("bayesm.mat","mcmc")
attributes(out$drtau)$mcpar=c(1,R,keep)
attributes(out$drsigma)$class=c("bayesm.mat","mcmc")
attributes(out$drsigma)$mcpar=c(1,R,keep)
attributes(out$drLambda)$class=c("bayesm.mat","mcmc")
attributes(out$drLambda)$mcpar=c(1,R,keep)
attributes(out$dre)$class=c("bayesm.mat","mcmc")
attributes(out$dre)$mcpar=c(1,R,keep)
attributes(out$drSigma)$class=c("bayesm.var","bayesm.mat","mcmc")
attributes(out$drSigma)$mcpar=c(1,R,keep)
return(list(Sigmadraw=out$drSigma,mudraw=out$drmu,taudraw = out$drtau,
            sigmadraw=out$drsigma,Lambdadraw=out$drLambda,edraw=out$dre))
}




