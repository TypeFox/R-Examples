# release 1.0.0 : first release
# release 1.0.1 : simplification
# release 1.0.2 : remove use.CV
# release 1.0.3 : correction of the data
# release 1.0.4 : correct export function names

#' @import polynom
#' @import orthopolynom
#' @import parallel
#' @importFrom graphics matplot par
#' @importFrom stats spline toeplitz

require(graphics)
require(stats)
require(parallel)
require(polynom)
require(orthopolynom)

int.trapeze = function(t,y) {
  sum(diff(t)*(y[-1]+y[-length(y)]))/2
}

#' function LaplaceConvolution
#'
#' computes the Laplace convolution of two functions f and g observed at discrete times t. Use trapezoidal formula and spline approximation of f.
#'
#' @param t, numeric vector, the observation times
#' @param g, numeric vector, the observed values of the known Laplace convolution kernel at the observation times
#' @param f, numeric vector, the coefficients the values of the function f to convole with g
#'
#' @return return the Laplace convolution of f and g using Trapezoidal formula and spline approximation for F
#' @author Y. Rozenholc and M. Pensky
#' @export LaplaceConvolution
#'
#' @examples
#'  \dontrun{
#'
#'  library(LaplaceDeconv)
#'  t = seq(0,10,l=100)
#'  g = exp(-5*t)
#'  f = t^2*exp(-t)
#'  # compute the Laplace convolution from functions computed at times t : f and g
#'  fg = LaplaceConvolution(t,g,f)
#'  matplot(t,cbind(f,g,fg),lty=1,type='l')
#'  legend('topright',lty=1,legend=c('f','g','fxg'),col=1:3)
#'  }

LaplaceConvolution=function(t,g,f) {
  ### WORKS FOR IRREGULAR DESIGN - USE SPLINE APPROXIMATION
  # add t=0 if needed
  if (t[1]>0) {t=c(0,t); f=c(0,f); g=c(0,g)}

  dt = diff(t)
  # Laplace convolution of two functions computed at times t
  LapConv=function(j){
    if (j==1) return(0)
    # compute approximation of F at time t[j]-t[i] for i<=j
    # z = approx(t,F,xout=t[j]-t[1:j])$y
    z = spline(t,f,xout=t[j]-t[1:j])$y
    y = g[1:j]*z
    sum(dt[1:(j-1)]*(y[-j]+y[-1]))/2
  }
  # return Laplace convolution approximated by Trapeze
  sapply(1:length(t),LapConv)
}

#' function LaguerreLaplaceConvolution
#'
#' computes the Laplace convolution of two functions f and g observed at discrete times t. Use trapezoidal formula and an expansion of f in the Laguerre function basis.
#'
#' @param t, numeric vector, the observation times
#' @param g, numeric vector, the observed values of the known Laplace convolution kernel at the observation times
#' @param f.coef, numeric vector, the coefficients in the Laguerre function basis of the function f to convole with g
#' @param a, numeric, the scale of the  Laguerre functions basis
#'
#' @return return the Laplace convolution of f and g using Trapezoidal formula and expansion of f in the Laguerre function basis
#' @author Y. Rozenholc and M. Pensky
#' @export LaguerreLaplaceConvolution
#'
#' @examples
#'  \dontrun{
#'
#'  library(LaplaceDeconv)
#'  a = 1/2
#'  t = seq(0,10,l=100)
#'  g = exp(-5*t)
#'  f.coef = c(1,0.25,0.1)
#'  # compute the Laplace convolution from g, kernel computed at times t, and the function described by
#'  # its decomposition in Laguerre function basis with scale a
#'  fg = LaguerreLaplaceConvolution(t,g,f.coef,a)
#'  matplot(t,cbind(MakeLaguerreMatrix(a,3)(t)%*%f.coef,g,fg),lty=1,type='l',ylab='')
#'  legend('topright',lty=1,legend=c('f','g','fxg'),col=1:3)
#'  }

LaguerreLaplaceConvolution=function(t,g,f.coef,a) {
  ### WORKS FOR IRREGULAR DESIGN - NO APPROXIMATION

  n = length(t)
  f.coef = as.matrix(f.coef)
  M = nrow(f.coef)
  dt = diff(t)

  # Laplace convolution of g and f computed at times t
  LapConv=function(j){
    if (j==1) return(rep(0,l=ncol(f.coef)))
    # compute f at time t[j]-t[i] for i<=j from its Laguerre coefficients
    f.rev = MakeLaguerreMatrix(a,M)(t[j]-t[1:j])%*%f.coef
    y = g[1:j]*f.rev
    colSums(dt[1:(j-1)]*(y[-j,,drop=F]+y[-1,,drop=F]))/2
  }
  # return Laplace convolution approximated by Trapeze
  matrix(sapply(1:n,LapConv),ncol=ncol(f.coef),byrow=TRUE)
}



#' function MakeLaguerreMatrix
#'
#' return a constructor of the Laguerre function basis.
#'
#' @param a, numeric, the scale of the  Laguerre functions basis
#' @param M, numeric, the maximal degree
#'
#' @return a constructor of the M-first Laguerre functions with scale a
#' @author Y. Rozenholc and M. Pensky
#' @export MakeLaguerreMatrix
#'
#' @examples
#' \dontrun{
#'
#'  library(LaplaceDeconv)
#'  # build the constructor
#'  bFctLagMat = MakeLaguerreMatrix(M=10)
#'  # compute the Laguerre function basis at times seq(0,15,l=100)
#'  bFctLagMat(seq(0,15,l=100))
#'  }

MakeLaguerreMatrix = function(a=1/2,M){
  # constructor of matrix of M first Laguerre functions (degree from 0 to M-1)
  LaguerreMatx = NULL
  LM = laguerre.polynomials(M-1,normalized=TRUE)
  function(x) {
    LM.x = polynomial.values(LM,2*a*x)
    E = exp(-a*x)*sqrt(2*a)
    LaguerreMatx <<- matrix(nrow=length(x),ncol=M)
    for (i in 1:M) LaguerreMatx[,i] <<- LM.x[[i]]*E
    LaguerreMatx
  }
}

#' function BuildLaguerreSystem
#'
#' build the matrix associated to the Laguerre function basis at given times, its QR decomposition, the coefficients of the kernel function of the Laplace convolution and the Toeplitz matrix associated to the Laplace convolution with the kernel g.
#'
#' @param a, numeric, the scale parameter of the Laguerre functions basis
#' @param g, numeric vector, the observed values of the known Laplace convolution kernel at the observation times
#' @param times, numeric vector, the observation times
#' @param Mmax, numeric, the maximal degree
#' @param tol, tolerance in qr decomposition (default 1e-7)
#' @author Y. Rozenholc and M. Pensky
#'
#' @return a list containing:
#' \itemize{
#' \item \code{gcoef}, numeric vector, the coefficients of g in the Laguerre functions basis
#' \item \code{GM}, numeric matrix, the Toeplitz matrix of the Laplace convolution with g
#' \item \code{Phi}, numeric matrix, the Laguerre function basis
#' \item \code{Phi.qr}, the QR decomposition of Phi
#' \item \code{Mmax}, numeric, the maximal degree
#' \item \code{a}, numeric, the scale of the Laguerre function basis
#' \item \code{g}, numeric vector, the kernel
#' \item \code{times}, numeric vector, the observation times
#' }
#' @export BuildLaguerreSystem
#'
#' @examples
#'  \dontrun{
#'
#'  library(LaplaceDeconv)
#'  t = seq(0,10,l=50)
#'  g = exp(-5*t)
#'  # compute the elements needed for the Laplace deconvolution with kernel observations g at times t
#'  bFctLagMat = BuildLaguerreSystem(a=1/2,g,t,M=10)
#'  }

BuildLaguerreSystem = function(a,g,times,Mmax,tol=1e-7){

  # don't use smaller tolerance than 1e-7 : does not work

  Phi = MakeLaguerreMatrix(a,Mmax)(times)

  Phi.qr = qr(Phi,tol=tol)
  if (Phi.qr$rank<Mmax)  # correct for Phi-rank degeneracy
    return(BuildLaguerreSystem(a,g,times,Mmax=Phi.qr$rank))

  # Find expansion of observations in Laguerre functions
  gcoef = backsolve(qr.R(Phi.qr),t(qr.Q(Phi.qr))%*%g,upper.tri=T)

  # Construction of the lower triangular Toeplitz matrix
  dg = diff(c(0,gcoef))
  GM = toeplitz(dg/sqrt(2*a))
  GM[upper.tri(GM)]=0					##### G_m = GM[1:m,1:m] is lower triangular

  # in order to be able to compute Qm we need GM %*% t(GM) to be of full rank
  rankGG = qr(GM%*%t(GM),tol=tol)$rank
  if (rankGG<Mmax) # correct for G-rank degeneracy
    return(BuildLaguerreSystem(a,g,times,Mmax=rankGG))

  list(gcoef=gcoef,GM=GM,Phi=Phi,Phi.qr=Phi.qr,Mmax=Mmax,a=a,g=g,times=times)
}

Normalizing = function(Y,g,times,sigma,pow=2){
  g.abs = abs(g)^pow
  g.Lpow = sum(diff(times)*(g.abs[-1]+g.abs[-length(g)])/2)^(1/pow)
  list(g=g/g.Lpow,Y=Y/g.Lpow,times=times,sigma=sigma/g.Lpow,g.Lpow=g.Lpow)
}

CheckConditionnementQR = function(a,g,times,Mmax,aleph){

  n = length(times)
  Tmax = times[n]

  LS = BuildLaguerreSystem(a,g,times,Mmax)

  Mmax = LS$Mmax
  while (Mmax>0) {
    PhiPhi.inv = chol2inv(qr.R(LS$Phi.qr)[1:Mmax,1:Mmax,drop=F])
    QM = PhiPhi.inv%*%chol2inv(t(LS$GM[1:Mmax,1:Mmax,drop=F]))*Tmax/n ### formula (2.11)
    if (sum(diag(QM))<=aleph) break
    Mmax = Mmax-1
  }

  if (Mmax>0) {
    LS = BuildLaguerreSystem(a,g,times,Mmax)
    LS$PhiPhi.inv = PhiPhi.inv
    Mmax = LS$Mmax
    LS$GM.inv = forwardsolve(LS$GM,diag(rep(1,Mmax)))
    LS$aleph = aleph
  }

  return(list(a=a,Mmax=Mmax,LS=LS))
}


#' function LagLaplDeconv
#'
#' Main function of this package : computes the Laplace deconvolution with noisy discrete non-equally spaced observations on a finite time interval (see reference).
#'
#' @param Y numeric vector, the observed noisy observations of the Laplace convolution
#' @param g numeric vector, the known kernel of the Laplace convolution
#' @param times numeric vector, the observation times (default \code{1:length(Y)})
#' @param sigma numeric, the noise level
#' @param cpen numeric, the penalization constant (default 2)
#' @param atab numeric vector, an array of value for a (default -1, see Details)
#' @param Mmax integer, the maximum degree (default 25)
#' @param ncores.max numeric, number of used cores (default \code{Inf})
#' @param verbose boolean to control information output (default \code{FALSE})
#' @param withplot boolean to control plot output (default \code{FALSE})
#' @details \code{atab} defines the values of the scale parameter used to build the Laguerre functions basis. If \code{atab} is of length 1 and negative then \code{atab<-seq(0.4,3,by=0.05/abs(atab))/sqrt(Tmax/10)} where \code{Tmax=max(times)}.
#' @author Y. Rozenholc and M. Pensky
#' @references \emph{Laplace deconvolution on the basis of time domain data and its application to Dynamic Contrast Enhanced imaging} by F. Comte, C-A. Cuenod, M. Pensky, Y. Rozenholc (ArXiv http://arxiv.org/abs/1405.7107)
#'
#' @return a list containing:
#' \itemize{
#' \item \code{f.hat}, numeric vector, the estimate at observation times
#' \item \code{q.hat}, numeric vector, the reconstructed convolution of the kernel g and the estimate f.hat, computed at observation times
#' \item \code{f.coef}, numeric vector, the coefficients of f.hat in the selected Laguerre function basis
#' \item \code{info}, only for internal use --- not documented
#' \item \code{a.hat}, numeric, the selected value of the parameter a
#' \item \code{M0}, numeric, the dimension of the selected model
#' \item \code{sigma}, numeric, the noise level
#' \item \code{cpen}, numeric, the penalization constant used in the penalty
#' }
#' @export LagLaplDeconv
#'
#' @aliases LaguerreLaplaceDeconvolution, LagLaplaceDeconvolution, LaplaceDeconvolution, LagLaplaceDeconv, LaplaceDeconv, LaplDeconv, LaguerrePenalizedQR
#' @examples
#'  \dontrun{
#'  #### AN ARTICIAL EXAMPLE ####
#'
#'  library(LaplaceDeconv)
#'  par(mfrow=c(1,1))
#'  set.seed(29102015)
#'
#'  sigma=0.02
#'  a = 1
#'  t = seq(0,5,l=100)
#'  g = 20*t^2*exp(-5*t)
#'  f.coef = c(0.4,0.02,0.01)
#'
#'  # compute the Laplace convolution from g, kernel computed at times t, and the function
#'  # described by its decomposition in Laguerre function basis with scale a :
#'  fg = LaguerreLaplaceConvolution(t,g,f.coef,a)
#'
#'  # the noisy observations :
#'  Y = fg+sigma*rnorm(length(fg))
#'
#'  # estimation of f from the observation and the kernel :
#'  L = LagLaplDeconv(Y,g,t,sigma)
#'  matplot(t,cbind(g,MakeLaguerreMatrix(a,3)(t)%*%f.coef,fg,L$q.hat,L$f.hat,Y),lty=1,
#'    type=c('b',rep('l',4),'p'),ylab='',pch='x')
#'
#'  # display results of estimation
#'  legend('topright',lty=c(rep(1,5),0),pch=c('x',rep('',4),'x'),
#'    legend=c(
#'      'g: partially observed kernel',
#'      'f: unknown',
#'      'q=fxg: unknown convolution',
#'      expression(hat(q)*': plug-in convolution'),
#'      expression(hat(f)*': estimation of f'),
#'      'Y: observations'),
#'    col=1:6)
#'  }
#'
#'  \dontrun{
#'  #### A REAL EXAMPLE USING DCE-MRI DATA FROM A TUMOR ####
#'
#'  library(LaplaceDeconv)
#'  par(mfrow=c(1,2))
#'
#'  # load data from patient before the treatment
#'  data(EX_DCEMRI_t0)
#'
#'  # display AIF and tumoral enhancements
#'  matplot(ex_dcemri$times,
#'    cbind(ex_dcemri$AIF,ex_dcemri$TUM_1,ex_dcemri$TUM_2,ex_dcemri$TUM_3),
#'    ylab='',lty=1,type=c('b',rep('p',3)),pch='+',main='Observations')
#'  legend('topright',pch='+',legend=c('AIF','TUM_1','TUM_2','TUM_3'),col=1:4)
#'
#'  # estimation of the contrast agent survival functions
#'  L1 = LagLaplDeconv(ex_dcemri$TUM_1,ex_dcemri$AIF,ex_dcemri$times,ex_dcemri$sigma)
#'  L2 = LagLaplDeconv(ex_dcemri$TUM_2,ex_dcemri$AIF,ex_dcemri$times,ex_dcemri$sigma)
#'  L3 = LagLaplDeconv(ex_dcemri$TUM_3,ex_dcemri$AIF,ex_dcemri$times,ex_dcemri$sigma)
#'
#'  matlines(ex_dcemri$times,cbind(L1$q.hat,L2$q.hat,L3$q.hat),type='l',lty=1,col=2:4)
#'
#'  # display results of estimation
#'  matplot(ex_dcemri$times,cbind(L1$f.hat,L2$f.hat,L3$f.hat),type='l',lty=1,col=2:4,
#'    ylab='survival',main='Contrast agent survival fcts')
#'  legend('topright',lty=1,col=2:4,
#'    legend=c(
#'      paste0('TUM_1 - a.hat=',round(L1$a.hat,digits=2)),
#'      paste0('TUM_2 - a.hat=',round(L2$a.hat,digits=2)),
#'      paste0('TUM_3 - a.hat=',round(L3$a.hat,digits=2))
#'      )
#'    )
#'  }

LagLaplDeconv = function(Y,g,times=1:length(Y),sigma,cpen=2,atab=-1,Mmax=25,
                         ncores.max=Inf,verbose=FALSE,withplot=FALSE){

  # Implement a choice of parameter a based on minimizing the residuals of the penalized estimate when a varies

  Y = as.matrix(Y)

  n = dim(Y)[1]
  Ntir = dim(Y)[2]

  if (times[1]>0) {
    print('Translate times to have times[1]=0')
    times = times-times[1]
  }

  use.multicore = (ncores.max>1)
  ncores.max = min(ncores.max,detectCores())

  L = Normalizing(Y,g,times,sigma)
  Y = L$Y; g = L$g; times = L$times; sigma = L$sigma

  Tmax = times[n]

  if ((length(atab)==1)&&(atab<0)) atab = seq(0.4,3,by=0.05/(-atab))/sqrt(Tmax/10) # new proposal

  M0tab = NA*atab
  CC.list = list()

  ################# INTERNAL FUNCTIONS ############################################
  DoLaguerrePenalization = function(ia,return_ctr=TRUE,index=1:Ntir,M0=NULL) {

    # Do minimum contrast estimation and return residuals in the observation space

    if (is.na(M0tab[ia])) {
      CC = CheckConditionnementQR(atab[ia],g,times,Mmax,aleph=n^0.8)
      M0tab[ia] <<- CC$Mmax
      CC.list[[ia]] <<- CC$LS
    }

    LaguerreSystem = CC.list[[ia]]
    M0 = M0tab[ia]

    if (M0==0) { # not good conditionning : exit
      if (return_ctr) return(rep(Inf,length(index)))
      else return(list(f.hat=NA*Y[index],q.hat=NA*Y[index],f.coef=NA,
                       info=c(M0=0,Mmax=Mmax,m.hat=NA,R2.hat=Inf,aleph=NA)))
    }

    # Extract informations from LaguerreSystem
    g.coef = LaguerreSystem$gcoef
    G.mat = LaguerreSystem$GM
    G.inv = LaguerreSystem$GM.inv
    Phi = LaguerreSystem$Phi
    Phi.qr = LaguerreSystem$Phi.qr
    PhiPhi.inv = LaguerreSystem$PhiPhi.inv

    # Find expansion of observations in Laguerre functions
    q.coef = backsolve(qr.R(Phi.qr),t(qr.Q(Phi.qr))%*%Y[,index,drop=F],upper.tri=TRUE) # qcoef = qr.coef(Phi.qr, Y)

    # Solve the system G.mat x fcoef = qcoef
    f.coef = backsolve(G.mat,q.coef,upper.tri=FALSE)

    Ctr = -apply(f.coef^2,2,cumsum)

    nu2_rho2 = function(m) {
      # formula (2.15)
      AAm = n/Tmax*G.inv[1:m,1:m,drop=F]%*%PhiPhi.inv[1:m,1:m,drop=F]%*%t(G.inv[1:m,1:m,drop=F])
      nu2 = sum(diag(AAm))
      rho2 = eigen(AAm,symmetric=TRUE,only.values=TRUE)$values[1] ##### Spectral radius of AAm
      return(c(nu2=nu2,rho2=rho2))
    }

    # Compute penalty
    Nu2Rho2 = sapply(1:M0,nu2_rho2)
    qen = Nu2Rho2['rho2',]*(log(1+Nu2Rho2['rho2',]/Nu2Rho2['rho2',1]) + cpen*log(1:M0))
    pen = 8*sigma^2*Tmax/n*(Nu2Rho2['nu2',] + qen)

    # penalized contrast
    CtrPen = Ctr+pen

    if (M0==1) CtrPen = matrix(CtrPen,nrow=1) # force to be a line matrix

    # selected dimension = degree+1
    m.hat = apply(CtrPen,2,which.min)

    # build estimates and reconstructions
    f.hat = q.hat = Y[,index,drop=F]
    tapply(1:length(index),as.factor(m.hat),function(ind) {
      m = m.hat[ind[1]];
      tmp.coef = f.coef[1:m,ind,drop=F]
      f.hat[,ind]<<-Phi[,1:m,drop=FALSE]%*%tmp.coef
      q.hat[,ind]<<-LaguerreLaplaceConvolution(times,g,tmp.coef,atab[ia])
      NULL
    })

    # cross validation criterion
    R2.hat = colMeans((Y[,index]-q.hat)^2)^0.5

    # WAS used for JRSS-v2 but does not change the results and is not in the theory
    # R2.hat = R2.hat + 2*sigma^2*m.hat*(1+log(m.hat)^3)/n # AKAIKE LIKE REGULARIZATION

    if (return_ctr)	return(R2.hat)
    else return(list(f.hat=f.hat,q.hat=q.hat,f.coef=f.coef,
                     info=c(M0=M0,Mmax=Mmax,m.hat=m.hat,R2.hat=R2.hat,aleph=LaguerreSystem$aleph)))

  }

  DoLaguerrePenalizationParallel = function(ind.a,index=1:Ntir,add.plot=FALSE) {
    if (!use.multicore) R2.hat = sapply(ind.a, DoLaguerrePenalization,index=index)
    else R2.hat = simplify2array(mclapply(ind.a, DoLaguerrePenalization,index=index,
                                          mc.cores=min(length(ind.a),ncores.max),mc.preschedule=FALSE))

    if (length(index)==1) R2.hat = matrix(R2.hat,nrow=1)

    if (verbose)
      matplot(atab,t(R2.hat),type='l',log='y',ylab='CV(a)',col=index,xlab='a',add=add.plot,lty=3-2*add.plot)

    ind  = apply(R2.hat,1,which.min)

    ind
  }
  ################# END INTERNAL FUNCTIONS ############################################

  if (withplot) par(mfrow=c(2+verbose,1))

  if (length(atab)>1) {
    # use the grid to optimize
    ind.a.opt = DoLaguerrePenalizationParallel(1:length(atab))
  } else ind.a.opt = 1

  f.hat = q.hat = Y
  m.hat = M0 = vector(length=Ntir)
  for (i in 1:Ntir) {
    D = DoLaguerrePenalization(ind.a.opt[i],return_ctr=FALSE,index=i)
    f.hat[,i] = D$f.hat
    q.hat[,i] = D$q.hat
    m.hat[i] = D$info['m.hat']
    M0[i] = D$info['M0']
    # if (verbose) print(D$Nu2Rho2)
  }

  q.hat = L$g.Lpow*q.hat; Y = L$g.Lpow*Y; g = L$g.Lpow*g; sigma = L$g.pow*sigma;

  if (withplot) {
    matplot(times,cbind(Y,q.hat),type=c('p','l'),pch=c('+',''),lty=1,xlab='times',ylab='q',col=c('red','blue'))
    matplot(times,f.hat,type='l',xlab='times',ylab='f')
  }

  if (Ntir==1)
    return(list(f.hat=f.hat,q.hat=q.hat,f.coef=D$f.coef,info=D$info,a.hat=atab[ind.a.opt],M0=M0,sigma=sigma,cpen=cpen))
  else
    return(list(f.hat=f.hat,q.hat=q.hat,a.hat=atab[ind.a.opt],m.hat=m.hat,M0=M0,Mmax=Mmax,sigma=sigma,cpen=cpen))
}


#' Data from two Dynamical Contrast Enhanced Computer Tomography (DCE-CT) image sequences realized one patient of the cohort REMISCAN before treatment starts (t0) and after 15 days of treatment (t1). These are the real data used to illustrate the paper in References.
#'
#' Each dataset contains a variable \code{ex_dcemri} which is a list containing
#' #' \itemize{
#'   \item \code{AIF}, numeric vector, the averaged enhancement curve recorded in the aorta, the convolution kernel
#'   \item \code{TUM_1}, \code{TUM_2}, \code{TUM_3}, numeric vectors, three averaged enhancement curves recorded in homogenous regions in the tumor, the observations
#'   \item \code{times}, numeric vector, the observation times
#'   \item \code{sigma}, numeric,  the noise level
#'   \item \code{delta}, numeric, equals 0, a time shift assumed to be less than 1 time unit
#' }
#
#' @docType data
#' @usage data(EX_DCEMRI_t0)
#' @usage data(EX_DCEMRI_t1)
#' @name EX_DCEMRI_t0.rda, EX_DCEMRI_t1.rda
#' @author C-A. Cuenod and Y. Rozenholc
#' @aliases EX_DCEMRI_t0, EX_DCEMRI_t1, AIF, TUM_1, TUM_2, TUM_3, times, sigma, delta, ex_dcemri
#' @references "Laplace deconvolution on the basis of time domain data and its application to Dynamic Contrast Enhanced imaging" by F. Comte, C-A. Cuenod, M. Pensky, Y. Rozenholc (ArXiv http://arxiv.org/abs/1405.7107)
#' @references REMISCAN - Project number IDRCB 2007-A00518-45/P060407/STIC 2006; Research Ethics Board (REB) approved- cohort funding by INCa (1M Euros) and promoted by the AP-HP (Assistance Publique Hôpitaux de Paris). Inclusion target: 100 patients. Start in 2007. Closed since July 2015.
#' @keywords data
NULL




