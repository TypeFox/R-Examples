################################################################################
# Functional Response Regression for (fd or fdata class objects)
# The function is a wrapped of linmod function proposed by Ramsay and Silverman (2005) 
# to model the relationship between the functional response and the functional covariate
# by basis representation of both. 
################################################################################
fregre.basis.fr<- function(x,y,basis.s=NULL,basis.t=NULL,
lambda.s=0,lambda.t=0,Lfdobj.s=vec2Lfd(c(0,0),range.s),
Lfdobj.t=vec2Lfd(c(0,0),range.t),weights=NULL,...){
call<-match.call()
isfdx<-is.fd(x)
x.orig<-x
isfdy<-is.fd(y)
y.orig<-y
if (isfdx) {
  xfdobj<-x
  basis.x<- x$basis
  nbasis.x<- basis.x$nbasis
  range.x<- basis.x$rangeval
  if (is.null(basis.s))  basis.s<-basis.x
}
else {
  if (!is.fdata(x))  stop("x is not a functional data object of class fd or fdata.")
    range.x<-x$rangeval
    if (is.null(basis.s))  {
        np<-ncol(x)
        #nbasis.b = min(51,floor(np/10))
        nbasis.b =min(7,floor(np/2))
        basis.s<-basis.x<-create.bspline.basis(rangeval=range.x,nbasis=nbasis.b)
   }
   else   basis.x<-basis.s
   xfdobj<-Data2fd(argvals =x$arg, y = t(x$data), basisobj = basis.s,...)
}
if (isfdy) {
  yfdobj<-y
  basis.y<-y$basis
  nbasis.y<-basis.y$nbasis
  range.y<-basis.y$rangeval
  np.y<-nbasis.y
  nfine = min(c(201, 10 * np.y+ 1))
  tty <- seq(range.y[1],range.y[2],len=nfine)
#tfine = seq(ranget[1], ranget[2], len = nfine)

  if (is.null(basis.t))  basis.t<-basis.y
}
else {
  if (!is.fdata(y))  stop("y is not a functional data object of class fd or fdata.")
  range.y<-y$rangeval
  np.y<-ncol(y)
  nbasis.y =min(7,floor(np.y/2))
#  tty<-y$argvals
  if (is.null(basis.t))  {
        #nbasis.y = min(51,floor(np/10))
        basis.y<-basis.t<-create.bspline.basis(rangeval=range.y,nbasis=nbasis.y)
   }
   else   basis.y<-basis.t
  nfine = min(c(201, 10 * np.y+ 1))
  tty <- seq(range.y[1],range.y[2],len=nfine)
   yfdobj<-Data2fd(argvals =y$argvals, y = t(y$data), basisobj = basis.t,...)
}
#xLfdobj<-vec2Lfd(c(0,0,0), c(range.x[1],range.x[2]))
coefy   = yfdobj$coef
coefx   = xfdobj$coef
coefdx  = dim(coefx)
coefdy  = dim(coefy)
n=ncurves = coefdx[2]
################################################################################
#alphabasis<-basis.t#create.constant.basis(range.y)
range.t = basis.t$rangeval
range.s = basis.s$rangeval
###fdPar(smallbasis3, xLfdobj, xlambda)$Lfd
alphafdPar = fdPar(basis.s,Lfdobj.s, lambda.s)
alphalambda = alphafdPar$lambda
alphabasis<-alphafdPar$fd$basis
alphanbasis = alphabasis$nbasis
Finprod = inprod(basis.y, alphabasis)
alphattmat = diff(range.y) #cambiar
alphattmat = inprod(alphabasis, alphabasis)
################################################################################
if (range.s[1] != range.x[1] || range.s[2] != range.x[2]) {
    stop("Range of beta.s and x not compatible.")
}
nbasis.s = basis.s$nbasis
Hinprod = inprod(basis.x, basis.s)
xcoef = xfdobj$coef
basis.ss  = inprod(basis.s, basis.s)
Q<-S<-R<- NULL
range.t = basis.t$rangeval
if (range.t[1] != range.y[1] || range.t[2] != range.y[2]) {
    stop("Range of BETATFD coefficient and YFD not compatible.")
}
nbasis.t = basis.t$nbasis
Ginprod = inprod(basis.y, basis.t)
ycoef = yfdobj$coef
basis.tt  = inprod(basis.t, basis.t)
basis.talpha= inprod(basis.t, alphabasis)
Fmat = t(ycoef) %*% Finprod
Gmat = t(ycoef) %*% Ginprod
Hmat = t(xcoef) %*% Hinprod
if (is.null(weights)) {
    HHCP = t(Hmat) %*% Hmat
    HGCP = t(Hmat) %*% Gmat
    H1CP = as.matrix(colSums(Hmat))
    F1CP = as.matrix(colSums(Fmat))
} else {
    HHCP = t(Hmat) %*% (outer(weights,rep(nbasis.s))*Hmat)
    HGCP = t(Hmat) %*% (outer(weights,rep(nbasis.t))*Gmat)
    H1CP = t(Hmat) %*% weights
    F1CP = t(Fmat) %*% weights
}
#############################
alphattmat = inprod(alphabasis, alphabasis)
betalttmat = inprod(basis.t, alphabasis)
betassmat = inprod(basis.s, basis.s)
betattmat = inprod(basis.t, basis.t)

Beta1fd  = bifd(matrix(0,nbasis.s,nbasis.t), basis.s, basis.t)#
Beta1Par = bifdPar(Beta1fd, Lfdobj.s, Lfdobj.t, lambda.s, lambda.t)
betaslambda<-Beta1Par$lambdas
betatlambda<-Beta1Par$lambdat

if (alphalambda > 0) {   Q = eval.penalty(alphabasis, alphafdPar$Lfd) }
if (lambda.s > 0) {   R = eval.penalty(basis.s,Beta1Par$Lfds)}
if (lambda.t > 0) {  S = eval.penalty(basis.t,Beta1Par$Lfdt)}

betan = nbasis.s*nbasis.t
ncoef = alphanbasis+betan
Cmat = matrix(0, ncoef, ncoef)
Dmat = matrix(0, ncoef, 1)
ind1 = 1:alphanbasis
ind2 = ind1
Cmat[ind1, ind2] = ncurves * alphattmat
if (alphalambda > 0) {    Cmat[ind1, ind2] = Cmat[ind1, ind2] + alphalambda * Q  }
ind2 = alphanbasis + (1:betan)
Cmat[ind1,ind2] = t(kronecker(H1CP,basis.talpha))
Dmat[ind1] = F1CP
ind1 = alphanbasis + (1:betan)
ind2 = 1:alphanbasis
Cmat[ind1, ind2] = t(Cmat[ind2, ind1])
ind2 = ind1
Cmat[ind1, ind2] = kronecker(HHCP, betattmat)

if (betaslambda> 0) {
  Cmat[ind1,ind2] = Cmat[ind1,ind2] + lambda.s*kronecker(R,basis.tt)
}
if (betatlambda > 0) {
    Cmat[ind1,ind2] = Cmat[ind1,ind2] + lambda.t*kronecker(basis.ss,S)
}
Dmat[ind1] = matrix(t(HGCP), betan, 1)
coefvec = symsolve(Cmat, Dmat)
ind1 = 1:alphanbasis
alpha.est = coefvec[ind1]
alphafdnames = yfdobj$fdnames
alphafdnames[[3]] = "Intercept"
alphafd = fd(alpha.est, alphabasis, alphafdnames)
ind1 = alphanbasis + (1:betan)
beta.est = matrix(coefvec[ind1],nbasis.t,nbasis.s)
betafdnames = xfdobj$fdnames
betafdnames[[3]] = "Reg. Coefficient"
betafd = bifd(t(beta.est), basis.s, basis.t, betafdnames)
beta.xest = beta.est %*% t(Hmat)
xbetafd = fd(beta.xest, basis.t)
if (isfdy) {
   yhatmat = eval.fd(tty, alphafd) %*% matrix(1, 1, ncurves)+   eval.fd(tty,  xbetafd)
   yhatfd = smooth.basis(tty, yhatmat, basis.y)$fd
   fitted.values  <-yhatfd # fd(coef=yhatmat, basisobj=basis.y, fdnames=yfdobj$fdnames)
 }
else {
   yhatmat = eval.fd(y$argvals, alphafd) %*% matrix(1, 1, ncurves)+   eval.fd(y$argvals,  xbetafd)
   fitted.values<-fdata(t(yhatmat),y$argvals,y$rangeval,y$names)
 }
residuals<-y-fitted.values
out = list(call = call, alpha.est = alphafd, coefficients = beta.est, 
             beta.estbifd = betafd, fitted.values = fitted.values, 
             residuals = residuals, lambda.s = lambda.s, lambda.t = lambda.t, 
             Lfdobj.s = Lfdobj.s, Lfdobj.t = Lfdobj.t, weights = weights, 
             x = x, y = y, H = Hinprod, basis.s = basis.s, basis.t = basis.t, 
             argvals.y = tty
            , betaspenmat= R, betatpenmat= S,alphapenmat=Q)
class(out)<-"fregre.fr"
return(out)
}
################################################################################