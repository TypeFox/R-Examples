linmod = function(xfdobj, yfdobj, betaList, wtvec=NULL)  {
#  LINMOD  Fits an unrestricted or full functional linear model of the form
#       y(t) = \alpha(t) + \int x(s) \beta(s,t) ds + \epsilon(t),
#  where
#       \beta(s,t) = \phi'(s) B \psi(t)
#  
#  Arguments:
#  XFD       a functional data object for the independent variable 
#  YFD       a functional data object for the   dependent variable 
#  BETACELL  a List array of length 2 containing a functional
#               parameter object for \alpha and a bifdPar object for
#               \beta as a function of s and t, respectively.
#  WTVEC    a vector of weights
#  Returns:  a list object linmodList with fields
#  BETA0ESTFD    a functional parameter object for \alpha
#  BETA1ESTBIFD  a bivariate functional parameter object for \beta
#  YHATFDOBJ     a functional data object for the approximation to y

#  Last modified 15 February 2014

#  check xfdobj and yfdobj

if (!is.fd(xfdobj)) {
    stop("XFD is not a functional data object.")
}

if (!is.fd(yfdobj)) {
    stop("YFD is not a functional data object.")
}

ybasis  = yfdobj$basis
ynbasis = ybasis$nbasis
ranget  = ybasis$rangeval

xbasis  = xfdobj$basis
ranges  = xbasis$rangeval

nfine = max(c(201,10*ynbasis+1))
tfine = seq(ranget[1],ranget[2],len=nfine)

#  get dimensions of data

coefy   = yfdobj$coef
coefx   = xfdobj$coef
coefdx  = dim(coefx)
coefdy  = dim(coefy)
ncurves = coefdx[2]
if (coefdy[2] != ncurves) {
    stop ("Numbers of observations in first two arguments do not match.")
}

#  set up or check weight vector

if (!is.null(wtvec)) wtvec = wtcheck(ncurves, wtvec)

#  get basis parameter objects

if (!inherits(betaList, "list")) stop("betaList is not a list object.")

if (length(betaList) != 2) stop("betaList not of length 2.")

alphafdPar  = betaList[[1]]
betabifdPar = betaList[[2]]

if (!inherits(alphafdPar, "fdPar")) {
    stop("BETACELL[[1]] is not a fdPar object.")
}
if (!inherits(betabifdPar, "bifdPar")) {
    stop("BETACELL[[2]] is not a bifdPar object.")
}

#  get Lfd objects

alphaLfd = alphafdPar$Lfd
betasLfd = betabifdPar$Lfds
betatLfd = betabifdPar$Lfdt

#  get smoothing parameters

alphalambda = alphafdPar$lambda
betaslambda = betabifdPar$lambdas
betatlambda = betabifdPar$lambdat

#  get basis objects

alphafd    = alphafdPar$fd

alphabasis = alphafd$basis
alpharange = alphabasis$rangeval
if (alpharange[1] != ranget[1] || alpharange[2] != ranget[2]) {
    stop("Range of ALPHAFD coefficient and YFD not compatible.")
}

betabifd = betabifdPar$bifd

betasbasis = betabifd$sbasis
betasrange = betasbasis$rangeval
if (betasrange[1] != ranges[1] || betasrange[2] != ranges[2]) {
    stop("Range of BETASFD coefficient and XFD not compatible.")
}

betatbasis = betabifd$tbasis
betatrange = betatbasis$rangeval
if (betatrange[1] != ranget[1] || betatrange[2] != ranget[2]) {
    stop("Range of BETATFD coefficient and YFD not compatible.")
}

#  get numbers of basis functions

alphanbasis = alphabasis$nbasis
betasnbasis = betasbasis$nbasis
betatnbasis = betatbasis$nbasis

#  get inner products of basis functions and data functions

Finprod = inprod(ybasis, alphabasis)
Ginprod = inprod(ybasis, betatbasis)
Hinprod = inprod(xbasis, betasbasis)

ycoef = yfdobj$coef
xcoef = xfdobj$coef
Fmat = t(ycoef) %*% Finprod
Gmat = t(ycoef) %*% Ginprod
Hmat = t(xcoef) %*% Hinprod

if (is.null(wtvec)) {
    HHCP = t(Hmat) %*% Hmat
    HGCP = t(Hmat) %*% Gmat
    H1CP = as.matrix(apply(Hmat,2,sum))
    F1CP = as.matrix(apply(Fmat,2,sum))
} else {
    HHCP = t(Hmat) %*% (outer(wtvec,rep(betasnbasis))*Hmat)
    HGCP = t(Hmat) %*% (outer(wtvec,rep(betatnbasis))*Gmat)
    H1CP = t(Hmat) %*% wtvec
    F1CP = t(Fmat) %*% wtvec
}

#  get inner products of basis functions

alphattmat = inprod(alphabasis, alphabasis)
betalttmat = inprod(betatbasis, alphabasis)
betassmat  = inprod(betasbasis, betasbasis)
betattmat  = inprod(betatbasis, betatbasis)

#  get penalty matrices

if (alphalambda > 0) {
    alphapenmat = eval.penalty(alphabasis, alphaLfd)
} else {
    alphapenmat = NULL
}
if (betaslambda > 0) {
    betaspenmat = eval.penalty(betasbasis, betasLfd)
} else {
    betaspenmat = NULL
}
if (betatlambda > 0) {
    betatpenmat = eval.penalty(betatbasis, betatLfd)
} else {
    betatpenmat = NULL
}

#  set up coefficient matrix and right side for stationary equations

betan = betasnbasis*betatnbasis
ncoef = alphanbasis + betan
Cmat  = matrix(0,ncoef,ncoef)
Dmat  = matrix(0,ncoef,1)

#  rows for alpha

ind1 = 1:alphanbasis
ind2 = ind1
Cmat[ind1,ind2] = ncurves*alphattmat
if (alphalambda > 0) {
    Cmat[ind1,ind2] = Cmat[ind1,ind2] + alphalambda*alphapenmat
}
ind2 = alphanbasis + (1:betan)
Cmat[ind1,ind2] = t(kronecker(H1CP,betalttmat))

Dmat[ind1] = F1CP

#  rows for beta

ind1 = alphanbasis + (1:betan)
ind2 = 1:alphanbasis
Cmat[ind1,ind2] = t(Cmat[ind2,ind1])
ind2 = ind1
Cmat[ind1,ind2] = kronecker(HHCP,betattmat)
if (betaslambda > 0) {
    Cmat[ind1,ind2] = Cmat[ind1,ind2] + 
                      betaslambda*kronecker(betaspenmat,betattmat)
}
if (betatlambda > 0) {
    Cmat[ind1,ind2] = Cmat[ind1,ind2] + 
                      betatlambda*kronecker(betassmat,betatpenmat)
}

Dmat[ind1] = matrix(t(HGCP),betan,1)

#  solve the equations

coefvec = symsolve(Cmat, Dmat)

#  set up the coefficient function estimates

#  functional structure for the alpha function

ind1 = 1:alphanbasis
alphacoef = coefvec[ind1]

alphafdnames = yfdobj$fdnames
alphafdnames[[3]] = "Intercept"
alphafd = fd(alphacoef, alphabasis, alphafdnames)

#  bi-functional structure for the beta function

ind1 = alphanbasis + (1:betan)
betacoef    = t(matrix(coefvec[ind1],betatnbasis,betasnbasis))
betafdnames = xfdobj$fdnames
betafdnames[[3]] = "Reg. Coefficient"
betafd = bifd(t(betacoef), betasbasis, betatbasis, betafdnames)

#  functional data structure for the yhat functions

xbetacoef = t(Hmat %*% betacoef)   
xbetafd   = fd(xbetacoef, betatbasis)
yhatmat   = eval.fd(tfine, alphafd) %*% matrix(1,1,ncurves) + 
            eval.fd(tfine, xbetafd)
yhatfd    = smooth.basis(tfine, yhatmat, ybasis)$fd 

linmodList = list(beta0estfd=alphafd, beta1estbifd=betafd, yhatfdobj=yhatfd)

return(linmodList)

}


