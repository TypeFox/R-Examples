data(Chocolate)

## Helper function:
responseName <- function(di) {
  ev <- environment(di)
  if ("lm" %in% class(ev$expr)) {
    deparse(formula(environment(di)$expr)[[2]])
  } else {
    "rt"
  }
}

## Color ramp to visualize desirabilities:
color <- c(rgb(1, 0, 0), grey(seq(0.01, 1, length.out=200)))

##########################################################################################
## Demo starts here:
summary(Chocolate)

oask <- devAskNewPage(TRUE)
## Use linear models as predictor:
m.E     <- lm(E ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)
m.d90   <- lm(d90 ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)
m.Fe    <- lm(Fe ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)
m.etaCa <- lm(etaCa ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)
m.tauCa <- lm(tauCa ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)

## Construct desirabilities
d.rt    <- derringerSuich(c(-Inf, 30, 45, 1, 1))
d.E     <- derringerSuich(c(-Inf, 3, 4, 1, 1))
d.d90   <- derringerSuich(c(20, 21, 23, 1, 1))
d.Fe    <- derringerSuich(c(-Inf, 20, 30, 1, 1))
d.etaCa <- derringerSuich(c(1, 1.5, 2, 1, 1))
d.tauCa <- derringerSuich(c(5, 8, 10, 1, 1))

d.Fe2 <- d.Fe
class(d.Fe2) <- c("dsA1", class(d.Fe2))

## Plot of desirabilities
opar <- par(mfrow=c(3, 2))
plot(d.rt, main="Runtime")
plot(d.E, main="Energy")
plot(d.d90, main="Particle size")
plot(d.Fe, main="Iron content")
plot(d.etaCa, main="Plastic viscosity")
plot(d.tauCa, main="yield value")
par(opar)

## Goal:
## Minimize runtime while keeping the rest of the parameters 'in spec'.
fnrt <- function(x) {
  if (is.data.frame(x)) {
    return (x$rt)
  } else if (is.matrix(x)){
    return (x[1,])
  } else {
    return (x[1])
  }
}

cd.rt     <- compositeDF(fnrt, d.rt)
cd.E      <- compositeDF(m.E, d.E)
cd.d90    <- compositeDF(m.d90, d.d90)
cd.Fe     <- compositeDF(m.Fe, d.Fe)
cd.etaCa  <- compositeDF(m.etaCa, d.etaCa)
cd.tauCa  <- compositeDF(m.tauCa, d.tauCa)

## Combine desirabilities using geometric DI:
di <- geometricDI(cd.rt, cd.E, cd.d90, cd.Fe, cd.etaCa, cd.tauCa,
                  weights=c(10, .1, 5, 1, 1, 1))

## Maximize desirability index:
of <- optim(c(40, 70), di, control=list(fnscale=-1))

## Visualize DI:
rt <- seq(20, 45, length.out=201)
as <- seq(5, 90, length.out=201)
z <- outer(rt, as, function(x,y) di(data.frame(rt=x, as=y)))
image(rt, as, z, col=color)
abline(v=of$par[1], h=of$par[2], col="red", lty=2)

## Individual desirabilities:
opar <- par(mfrow=c(3, 2))
for (dfn in c(cd.rt, cd.E, cd.d90, cd.Fe, cd.etaCa, cd.tauCa)) {
  z <- outer(rt, as, function(x, y) dfn(data.frame(rt=x, as=y)))
  image(rt, as, z, col=color, main=responseName(dfn))
  abline(v=of$par[1], h=of$par[2], col="red", lty=2)
}
par(opar)

##########################################################################################
## Realistic desirabilities:
rcd.E      <- compositeDF(m.E, realisticDF(d.E))
rcd.d90    <- compositeDF(m.d90, realisticDF(d.d90))
rcd.Fe     <- compositeDF(m.Fe, realisticDF(d.Fe))
rcd.etaCa  <- compositeDF(m.etaCa, realisticDF(d.etaCa))
rcd.tauCa  <- compositeDF(m.tauCa, realisticDF(d.tauCa))

rdi <- geometricDI(cd.rt, rcd.E, rcd.d90, rcd.Fe, rcd.etaCa, rcd.tauCa,
                   weights=c(10, .1, 5, 1, 1, 1))

## Maximize desirability index:
og <- optim(c(40, 70), rdi, control=list(fnscale=-1))

## Deterministic/Idealized desirability of solution is 0!
di(og$par)

rz <- outer(rt, as, function(x,y) rdi(data.frame(rt=x, as=y)))

## Response surface:
image(rt, as, rz, col=color)
abline(v=og$par[1], h=og$par[2], col="red", lty=2)
abline(v=of$par[1], h=of$par[2], col="blue", lty=2)

## Individual desirabilities:
opar <- par(mfrow=c(3, 2))
for (dfn in c(cd.rt, rcd.E, rcd.d90, rcd.Fe, rcd.etaCa, rcd.tauCa)) {
  z <- outer(rt, as, function(x, y) dfn(data.frame(rt=x, as=y)))
  image(rt, as, z, col=color, main=responseName(dfn))
  abline(v=of$par[1], h=of$par[2], col="red", lty=2)
}
par(opar)
devAskNewPage(oask)
