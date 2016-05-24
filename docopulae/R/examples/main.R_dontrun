\donttest{
library(copula)


## prepare for SparseGrid integration
ncube = function(dimension) {
    SparseGrid::createIntegrationGrid('GQU', dimension, 3)
}
ncube = nint_integrateNCube_SparseGrid(ncube)
unlockBinding('nint_integrateNCube', environment(nint_integrate))
assign('nint_integrateNCube', ncube, envir=environment(nint_integrate))


## general settings
numDeriv = FALSE

## copula
copula = claytonCopula()
alphas = list(alpha=iTau(copula, 0.5))

## thetas
thetas = c(names(alphas), c('beta1', 'beta2', 'beta3',
                            'beta4', 'beta5', 'beta6'))

## eta
eta = function(theta)
    c(c(theta$beta1, theta$beta2, theta$beta3) %*% (theta$x ** c(0, 1, 2)),
      c(theta$beta4, theta$beta5, theta$beta6) %*% (theta$x ** c(1, 3, 4)))


if (numDeriv) {
    ## margins
    margins = function(y, theta, ...) {
        e = eta(theta)
        cbind(dnorm(y, e), pnorm(y, e))
    }

    ## f
    f = buildf(margins, copula, names=names(alphas))

    ## 2nd derivatives
    d2logf = numDeriv2Logf(f)

} else {
    ## margins
    eta1 = quote(beta1   + beta2*x    + beta3*x**2)
    eta2 = quote(beta4*x + beta5*x**3 + beta6*x**4)

    margins = list(list(pdf=substitute(dnorm(y1, e, 1), list(e=eta1)),
                        cdf=substitute(pnorm(y1, e, 1), list(e=eta1))),
                   list(pdf=substitute(dnorm(y2, e, 1), list(e=eta2)),
                        cdf=substitute(pnorm(y2, e, 1), list(e=eta2))))

    ## mappings
    yMap = list(y1=1, y2=2)

    alphaMap = as.list(names(alphas))
    names(alphaMap) = names(alphas)
    thetaMap = c(alphaMap,
                 list(beta1='beta1', beta2='beta2', beta3='beta3',
                      beta4='beta4', beta5='beta5', beta6='beta6',
                      x='x'))

    ## f
    ff = buildf(margins, copula)
    f = expr2f(ff, yMap=yMap, thetaMap=thetaMap)

    ## 2nd derivatives
    cat('building derivatives ...')
    tt = system.time(d2logf <- Deriv2Logf(ff, thetas,
                                          yMap=yMap, thetaMap=thetaMap))
    cat('\n')
    print(tt)
}

f
str(d2logf)


## param
model = function(theta) {
    integrand = function(y, theta, i, j)
        -d2logf(y, theta, i, j) * f(y, theta)

    yspace = nint_space(nint_intvDim(-Inf, Inf),
                        nint_intvDim(-Inf, Inf))

    fisherIf = function(x) {
        theta$x = x

        ## probability integral transform
        e = eta(theta)

        tt = nint_transform(integrand, yspace, 1:2, list(
            g=function(x) pnorm(x, e),
            gij=function(y) {
                t1 = qnorm(y, e)
                cbind(t1, 1/dnorm(t1, e))
            }
        ))

        fisherI(tt$f, theta, thetas, tt$space)
    }

    return(param(fisherIf, 1))
}

theta = c(alphas, list(beta1=1, beta2=1, beta3=1,
                       beta4=1, beta5=1, beta6=1,
                       x=0))
m = model(theta)

## update.param
system.time(m <- update(m, matrix(seq(0, 1, length.out=101), ncol=1)))

## find D-optimal design
D = Dsensitivity(defaults=list(x=m$x, desx=m$x, mod=m))

system.time(d <- FedorovWynn(D, 7.0007, maxIter=1e3))
d$tag$FedorovWynn$tolBreak

getM(m, d)
dev.new(); plot(d, sensTol=7, main='d')

rd = reduce(d, 0.05)

try(getM(m, rd))
m2 = update(m, rd)
getM(m2, rd)

dev.new(); plot(rd, main='rd')
dev.new(); plot(rd, sensx=d$x, sens=D(x=d$x, desx=rd$x, desw=rd$w, mod=m2),
                main='rd + sensitivity')

## find Ds-optimal design
dsNames = c(names(alphas), 'beta1', 'beta2', 'beta3')
Ds = Dsensitivity(dsNames, defaults=list(x=m$x, desx=m$x, mod=m))

system.time(ds <- FedorovWynn(Ds, 4.0004, maxIter=1e3))
ds$tag$FedorovWynn$tolBreak

dev.new(); plot(ds, sensTol=4, main='ds')

## create custom design
n = 4
d2 = design(x=matrix(seq(0, 1, length.out=n), ncol=1), w=rep(1/n, n))

m = update(m, d2)
dev.new(); plot(d2, sensx=d$x, sens=D(x=d$x, desx=d2$x, desw=d2$w, mod=m),
                sensTol=7, main='d2 + sensitivity')

## compare designs
Defficiency(ds, d, m)
Defficiency(d, ds, m, dsNames=dsNames) # Ds-efficiency
Defficiency(d2, d, m)
Defficiency(d2, ds, m) # D-efficiency

## end with nice plot
dev.new(); plot(rd, sensx=d$x, sens=D(x=d$x, desx=rd$x, desw=rd$w, mod=m2),
                main='rd + sensitivity')
}
