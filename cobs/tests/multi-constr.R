#### Examples which use the new feature of more than one 'constraint'.

suppressMessages(library(cobs))

options(digits = 6)
if(!dev.interactive(orNone=TRUE)) pdf("multi-constr.pdf")

source(system.file("util.R", package = "cobs"))

op <- options(warn = 2) ## << all warnings to errors!

set.seed(101)
x <- seq(-2,2, length = 100)
y <- x^2 + 0.5*rnorm(100)
## Constraints -- choosing ones that are true for  f(x) = x^2
PW  <- rbind(
             c(0, -3,9), # f(-3) = 9
             c(0,  3,9), # f(3 ) = 9
             c(2,  0,0)) # f'(0) = 0

mod <- cobs (x,y,constraint = "convex", pointwise = PW)
mod
stopifnot(all.equal(predict(mod, c(-3, 3))[,"fit"], c(9,9), tol = 1e-12))

## derivative 0 at 0 -- we miss a 'deriv = 1' argument [-> see ../TODO]
eps <- 1e-6
stopifnot(abs(diff(predict(mod, c(-eps, eps))[,"fit"])/(2*eps)) < .001 * eps)
options(op)# allow warnings


set.seed(908)
x <- seq(-1,2, len = 50)
f.true <- pnorm(2*x)
y <- f.true + rnorm(50)/10
plot(x,y); lines(x, f.true, col="gray", lwd=2, lty=3)

## constraint on derivative at right end:
(con <- rbind(c(2 , max(x), 0))) # f'(x_n) == 0

## Using 'trace = 3' --> 'trace = 2' inside drqssbc2()

## Regression splines (lambda = 0)
c2   <- cobs(x,y, trace = 3)
c2i  <- cobs(x,y, constraint = c("increase"), trace = 3)
c2c  <- cobs(x,y, constraint = c("concave"), trace = 3)

c2IC <- cobs(x,y, constraint = c("inc", "concave"), trace = 3)
## here, it's the same as just "i":
all.equal(fitted(c2i), fitted(c2IC))

c1   <- cobs(x,y, degree = 1, trace = 3)
c1i  <- cobs(x,y, degree = 1, constraint = c("increase"), trace = 3)
c1c  <- cobs(x,y, degree = 1, constraint = c("concave"), trace = 3)

plot(c1)
lines(predict(c1i), col="forest green")
all.equal(fitted(c1), fitted(c1i), tol = 1e-9)# but not 1e-10

## now gives warning (not error):
c1IC <- cobs(x,y, degree = 1, constraint = c("inc", "concave"), trace = 3)

cp2   <- cobs(x,y,                          pointwise = con, trace = 3)
cp2i  <- cobs(x,y, constraint = "increase", pointwise = con)# warn: check 'ifl'
## when plotting it, we see that it gave a trivial constant!!
cp2c  <- cobs(x,y, constraint = "concave",  pointwise = con, trace = 3)

## now gives warning (not error):
cp2IC <- cobs(x,y, constraint = c("inc", "concave"), pointwise = con, trace = 3)

cp1   <- cobs(x,y, degree = 1,                            pointwise = con, trace = 3)
cp1i  <- cobs(x,y, degree = 1, constraint = "increase",   pointwise = con, trace = 3)
cp1c  <- cobs(x,y, degree = 1, constraint = "concave",    pointwise = con, trace = 3)

cp1IC <- cobs(x,y, degree = 1, constraint = c("inc", "concave"), pointwise = con, trace = 3)


plot(x,y, main = "cobs(*, degree= 1, constraint = *, pointwise= *)")
matlines(x,cbind(fitted(c1),
                 fitted(c1i),
                 fitted(c1c),
                 fitted(cp1),
                 fitted(cp1i),
                 fitted(cp1c)),
        col = 1:6, lty=1)
legend("bottomright", inset = .02, col = 1:6, lty=1,
       legend = c("none", "increase","concave",
       "pt", "pt + incr.", "pt + conc."))

if(dev.interactive()) x11() # cheap way to get plot in new window, when testing

plot(x,y, main = "cobs(*, degree= 2, constraint = *, pointwise= *)")
matlines(x,cbind(fitted(c2),
                 fitted(c2i),
                 fitted(c2c),
                 fitted(cp2),
                 fitted(cp2i),
                 fitted(cp2c)),
        col = 1:6, lty=1)
legend("bottomright", inset = .02, col = 1:6, lty=1,
       legend = c("none", "increase","concave",
       "pt", "pt + incr.", "pt + conc."))

##--> "increase + pointwise" gives constant which seems plain wrong  <<<< BUG ???
