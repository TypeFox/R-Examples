## ----FNver,echo=FALSE,results='asis',cache=FALSE-------------------------
FNver <- packageDescription("FuzzyNumbers")$Version
if (as.integer(substring(FNver, nchar(FNver))) == "0") {
   cat("\\DevelopmentVersiontrue\n")
} else {
   cat("\\DevelopmentVersionfalse\n")
}

## ----echo=FALSE,results='hide',warning=FALSE,message=FALSE,cache=FALSE----
options(digits=7)
options(width=73)
require('knitr')
# require('tikzDevice')
#
# options(tikzDefaultEngine = 'pdftex')
#
# options(tikzLatexPackages = c( # dolaczanie uzywanych pakietow TeX-a
#    '\\usepackage{amsmath,amssymb,amsfonts}', # pakiety AMS
#    '\\usepackage{tikz}',
# #   '\\usepackage[MeX,T1,plmath]{polski}', # obsluga m.in. polskich ogonkow
#    '\\usepackage[utf8]{inputenc}',
#    '\\usepackage[T1]{fontenc}',
#    '\\usetikzlibrary{calc}',
#    '\\usepackage[english]{babel}',
#    '\\selectlanguage{english}',
#    '\\usepackage{standalone}'
# ))
#
# options(tikzMetricsDictionary='~/R/tikzMetrics')
#
# options(tikzDocumentDeclaration = '\\documentclass[11pt]{standalone}\n')
#
# options(tikzMetricPackages = c(
#    '\\usepackage[utf8]{inputenc}',
#    '\\usepackage[T1]{fontenc}',
#    '\\usepackage{amsmath,amssymb,amsfonts}',
#    '\\usetikzlibrary{calc}',
#    '\\usepackage[english]{babel}',
#    '\\selectlanguage{english}'
# ))



# opts_knit$set(progress = TRUE, verbose = TRUE)

opts_chunk$set(
   keep.source=TRUE,
   out.width='4.5in',
   fig.width=6,
   fig.height=6/sqrt(3),
#    fig.path='figures-knitr/',
#    cache.path='cache-knitr/',
   cache=TRUE,
   tidy=FALSE,
#    dev='cairo_pdf',
#    dev.args=list(pointsize=11),
#    dev='tikz',
#    external=TRUE,
   fig.align='center',
   size='small'
)

# knit_theme$set(knit_theme$get('solarized-light'))

## ----eval=FALSE--------------------------------------------------------
#  install.packages('FuzzyNumbers')

## ----eval=FALSE--------------------------------------------------------
#  install.packages('devtools')
#  library('devtools')
#  install_github('FuzzyNumbers', 'Rexamine')

## ----libraryFuzzyNumbers,results='hide',warning=FALSE,message=FALSE----
library('FuzzyNumbers') # Load the package

## ----helpFuzzyNumbers,eval=FALSE---------------------------------------
#  library(help='FuzzyNumbers')

## ----A1ex--------------------------------------------------------------
A1 <- FuzzyNumber(1, 2, 4, 7,
    left=function(x) x,
   right=function(x) 1-x
)

## ----A1ex2,dependson='A1ex'--------------------------------------------
class(A1)

## ----A1ex3,dependson='A1ex'--------------------------------------------
A1

## ----A1ex4,dependson='A1ex',fig.keep='none'----------------------------
plot(A1)

## ----A1ex5,dependson='A1ex',echo=FALSE---------------------------------
par(mar=c(4,4,2,1))
plot(A1, xlab=expression(x), ylab=expression(alpha), main='')

## ----ExFconvert--------------------------------------------------------
f <- splinefun(c(-4,-3.5,-3,-2.2,-2), c(0,0.4,0.7,0.9,1), method='monoH.FC')
g <- splinefun(c(-1,0,10), c(1,0.5,0), method='monoH.FC')

## ----ExFconvert2,dependson='ExFconvert'--------------------------------
convertSide(f, -4, -2)(c(0,1))
convertSide(g, -1, 10)(c(0,1))
convertSide(g, 10, -1)(c(0,1)) # interesting!

## ----ExFconvert3,dependson='ExFconvert2',fig.keep='none'---------------
B <- FuzzyNumber(10,20,20,30,
    left=convertSide(f, -4, -2),
   right=convertSide(g, -1, 10)
)
plot(B, xlab=expression(x), ylab=expression(alpha))

## ----ExFconvert4,dependson='ExFconvert3',echo=FALSE--------------------
par(mar=c(4,4,2,1))
B <- FuzzyNumber(10,20,20,30,
    left=convertSide(f, -4, -2),
   right=convertSide(g, -1, 10)
)
plot(B, xlab=expression(x), ylab=expression(alpha))

## ----alphacutEx,fig.keep='none'----------------------------------------
A1 <- FuzzyNumber(1, 2, 4, 7,
    left=function(x) x,
   right=function(x) 1-x
)
A2 <- FuzzyNumber(1, 3, 4, 7,
   lower=function(alpha) pbeta(alpha, 5, 9), # CDF of a beta distr.
   upper=function(alpha) pexp(1/alpha-1) # transformed CDF of an exp. distr.
)
plot(A1, col='blue')
plot(A2, col='red', lty=2, add=TRUE)
legend('topright', c(expression(mu[A1]), expression(mu[A2])),
   col=c('blue', 'red'), lty=c(1,2))

## ----alphacutEx2,dependson='alphacutEx',echo=FALSE---------------------
par(mar=c(4,4,2,1))
plot(A1, col='blue', xlab=expression(x), ylab=expression(alpha))
plot(A2, col='red', lty=2, add=TRUE)
legend('topright', expression(mu[A[1]], mu[A[2]]),
   col=c('blue', 'red'), lty=c(1,2))

## ----A3def-------------------------------------------------------------
A3 <- FuzzyNumber(1, 2, 4, 5)
A3

## ----A3deffig,dependson='A3def',fig.keep='none'------------------------
plot(A3)

## ----A3deffig2,dependson='A3def',echo=FALSE----------------------------
par(mar=c(4,4,2,1))
plot(A3, xlab=expression(x), ylab=expression(alpha), shadowcol='gray')

## ----alphacutEx3,dependson='alphacutEx'--------------------------------
alphacut(A2, 0.5) # A2 has alpha-cut generators defined
alphacut(A1, 0.5) # A1 hasn't got them

## ----alphacutEx4,dependson='alphacutEx'--------------------------------
evaluate(A1, 6.5) # A1 has side generators defined
evaluate(A2, 6.5) # A2 hasn't got them

## ----alphacutEx5,dependson='alphacutEx'--------------------------------
A2['lower']
A2['upper']
A2['left']
A2['right']

## ----approxfuns--------------------------------------------------------
l <- function(x) pbeta(x, 1, 2)
r <- function(x) 1-pbeta(x, 1, 0.1)
A4 <- FuzzyNumber(-2, 0, 0, 2,
   left  = l,
   right = r,
   lower = approxInvert(l),
   upper = approxInvert(r)
)

x <- seq(0,1,length.out=1e5)
max(abs(qbeta(x, 1, 2) - A4['lower'](x)))     # sup-error estimator
max(abs(qbeta(1-x, 1, 0.1) - A4['upper'](x))) # sup-error estimator

## ----TrapEx1a----------------------------------------------------------
T1 <- TrapezoidalFuzzyNumber(1, 1.5, 4, 7)

## ----TrapEx1f,dependson='TrapEx1e',results='asis',echo=FALSE-----------
cat(as.character(T1, toLaTeX=TRUE, varnameLaTeX='T_1'))

## ----TrapEx1b,dependson='TrapEx1a'-------------------------------------
class(T1)

## ----TrapEx1c,dependson='TrapEx1b',fig.keep='none'---------------------
plot(T1)

## ----TrapEx1d,dependson='TrapEx1c',echo=FALSE--------------------------
par(mar=c(4,4,2,1))
plot(T1, xlab=expression(x), ylab=expression(alpha))

## ----TrapEx1e,dependson='TrapEx1d'-------------------------------------
T1['lower']
T1['upper']
T1['left']
T1['right']

## ----TrapEx2-----------------------------------------------------------
TrapezoidalFuzzyNumber(1,2,2,3)   # triangular FN
TriangularFuzzyNumber(1,2,3)      # the same
TrapezoidalFuzzyNumber(2,2,3,3)   # `crisp' interval
as.TrapezoidalFuzzyNumber(c(2,3)) # the same
TrapezoidalFuzzyNumber(5,5,5,5)   # `crisp' real
as.TrapezoidalFuzzyNumber(5)      # the same

## ----PLFNEx1a,fig.keep='none'------------------------------------------
P1 <- PiecewiseLinearFuzzyNumber(1, 2, 3, 4,
   knot.n=1, knot.alpha=0.25, knot.left=1.5, knot.right=3.25)
class(P1)
P1
P2 <- PiecewiseLinearFuzzyNumber(1, 2, 3, 4,
   knot.n=2, knot.alpha=c(0.25,0.6),
   knot.left=c(1.5,1.8), knot.right=c(3.25, 3.5))
P2
plot(P1, type='b', from=0, to=5, xlim=c(0.5,4.5))
plot(P2, type='b', col=2, lty=2, pch=2, add=TRUE, from=0, to=5)

## ----PLFNEx1b,dependson='PLFNEx1a',echo=FALSE--------------------------
par(mar=c(4,4,2,1))
plot(P1, type='b', xlab=expression(x), ylab=expression(alpha), from=0, to=5, xlim=c(0.5,4.5))
plot(P2, type='b', col=2, lty=2, pch=2, add=TRUE, from=0, to=5)

## ----PLFNEx1c,dependson='PLFNEx1b',------------------------------------
P1['knots']
P1['allknots'] # including a1,a2,a3,a4

## ----PLFNEx1d,dependson='PLFNEx1c',results='asis',echo=FALSE-----------
cat(as.character(P1, toLaTeX=TRUE, varnameLaTeX='P_1'))

## ----PLFNEx3-----------------------------------------------------------
PiecewiseLinearFuzzyNumber(knot.left=c(0,0.5,0.7,1),
                           knot.right=c(2,2.2,2.4,3))['allknots']

## ----PLFNEx2a,fig.keep='none'------------------------------------------
alpha <- c(0.3, 0.5, 0.7)
P3 <- as.PiecewiseLinearFuzzyNumber(
   TrapezoidalFuzzyNumber(1,2.5,4,7),
         knot.n=3, knot.alpha=alpha
)
P3
plot(P3, type='b', from=-1, to=9, xlim=c(0,8))
abline(h=alpha, col='gray', lty=2)
abline(v=P3['knot.left'], col='gray', lty=3)
abline(v=P3['knot.right'], col='gray', lty=3)
text(7.5, alpha, sprintf('a=%g', alpha), pos=3)

## ----PLFNEx2b,dependson='PLFNEx2a',echo=FALSE--------------------------
par(mar=c(4,4,2,1))
plot(P3, type='b', xlab=expression(x), ylab=expression(alpha), from=-1, to=9, xlim=c(0,8))
abline(h=alpha, col='gray', lty=2)
abline(v=P3['knot.left'], col='gray', lty=3)
abline(v=P3['knot.right'], col='gray', lty=3)
text(7.5, alpha, parse(text=sprintf('alpha*"="*%g', alpha)), pos=3)

## ----PLFNEx2c,dependson='PLFNEx2b',------------------------------------
(as.FuzzyNumber(P3))

## ----PowerEx1a,fig.keep='none'-----------------------------------------
X <- PowerFuzzyNumber(-3, -1, 1, 3, p.left=2, p.right=0.1)
class(X)
X
plot(X)

## ----PowerEx1b,dependson='PowerEx1a',echo=FALSE------------------------
par(mar=c(4,4,2,1))
plot(X, type='l', xlab=expression(x), ylab=expression(alpha))

## ----PowerEx1c,dependson='PowerEx1b',results='asis',echo=FALSE---------
cat(as.character(X, toLaTeX=TRUE, varnameLaTeX='X'))

## ----depicting1a,fig.keep='none'---------------------------------------
A <- FuzzyNumber(-5, 3, 6, 20,
    left=function(x) pbeta(x,0.4,3),
   right=function(x) 1-x^(1/4),
   lower=function(alpha) qbeta(alpha,0.4,3),
   upper=function(alpha) (1-alpha)^4
)
plot(A)

## ----depicting1b,dependson='depicting1a',echo=FALSE--------------------
par(mar=c(4,4,2,1))
plot(A, xlab=expression(x), ylab=expression(alpha))

## ----depicting1c,dependson='depicting1b',fig.keep='none'---------------
plot(A, n=3, type='b')
plot(A, n=6, add=TRUE,  lty=2, col=2, type='b', pch=2)
plot(A, n=101, add=TRUE, lty=4, col=4) # default n

## ----depicting1d,dependson='depicting1c',echo=FALSE--------------------
par(mar=c(4,4,2,1))
plot(A, n=3, type='b', xlab=expression(x), ylab=expression(alpha))
plot(A, n=6, add=TRUE, lty=2, col=2, type='b', pch=2)
plot(A, n=101, add=TRUE, lty=4, col=4) # default n

## ----depicting1e,dependson='depicting1d',fig.keep='none'---------------
plot(A, n=3, at.alpha=numeric(0), type='b') # use alpha-cuts
plot(A, n=3, type='b', col=2, lty=2, pch=2, add=TRUE) # use sides

## ----depicting1f,dependson='depicting1e',echo=FALSE--------------------
par(mar=c(4,4,2,1))
plot(A, n=3, at.alpha=numeric(0), type='b', xlab=expression(x), ylab=expression(alpha)) # use alpha-cuts
plot(A, n=3, type='b', col=2, lty=2, pch=2, add=TRUE) # use side generators

## ----depicting1g,dependson='depicting1f',fig.keep='none'---------------
plot(A, draw.alphacuts=TRUE)

## ----depicting1h,dependson='depicting1g',echo=FALSE--------------------
par(mar=c(4,4,2,1))
plot(A, draw.alphacuts=TRUE, ylab=expression(x), xlab=expression(alpha))

## ----A1ex6,dependson='A1ex',eval=FALSE---------------------------------
#  pdf('figure1.pdf', width=8, height=5) # create file
#  plot(A)
#  dev.off() # close graphical device and save the file

## ----depicting1i,dependson='depicting1h',eval=FALSE--------------------
#  cat(as.character(A, toLaTeX=TRUE, varnameLaTeX='A'))

## ----depicting1i2,dependson='depicting1i',results='asis',echo=FALSE----
cat(as.character(A, toLaTeX=TRUE, varnameLaTeX='A'))

## ----depicting1i3,dependson='depicting1i',results='asis',echo=FALSE----
cat(as.character(A, toLaTeX=TRUE, varnameLaTeX='A'))

## ----depicting2a,fig.keep='none'---------------------------------------
X <- PiecewiseLinearFuzzyNumber(0, 1, 2, 5, knot.n=1,
   knot.alpha=0.6, knot.left=0.3, knot.right=4)

plot.default(NA, xlab=expression(x), ylab=expression(mu[S](x)),
   xlim=c(-0.3,5.3), ylim=c(0,1)) # empty window

xpos <- c(X['a1'], X['knot.left'],  X['a2'],
          X['a3'], X['knot.right'], X['a4'])
xlab <- expression(s[1], s[2], s[3], s[4], s[5], s[6])
abline(v=xpos, col='gray', lty=3)
text(xpos, 1.05, xlab, pos=3, xpd=TRUE)

abline(h=c(0, X['knot.alpha'], 1), col='gray', lty=2)
text(5.1, X['knot.alpha'], expression(alpha[0]), pos=4, xpd=TRUE)

plot(X, add=TRUE, type='l', from=-1, to=6)
plot(X, add=TRUE, type='p', from=-1, to=6)

## ----depicting2b,dependson='depicting2a',echo=FALSE--------------------
par(mar=c(4,4,2,1))
X <- PiecewiseLinearFuzzyNumber(0, 1, 2, 5, knot.n=1,
   knot.alpha=0.6, knot.left=0.3, knot.right=4)

plot.default(NA, xlab=expression(x), ylab=expression(mu[S](x)),
   xlim=c(-0.3,5.3), ylim=c(0,1)) # empty window

xpos <- c(X['a1'], X['knot.left'],  X['a2'],
          X['a3'], X['knot.right'], X['a4'])
xlab <- expression(s[1], s[2], s[3], s[4], s[5], s[6])
abline(v=xpos, col='gray', lty=3)
text(xpos, 1.05, xlab, pos=3, xpd=TRUE)

abline(h=c(0, X['knot.alpha'], 1), col='gray', lty=2)
text(5.1, X['knot.alpha'], expression(alpha[0]), pos=4, xpd=TRUE)

plot(X, add=TRUE, type='l', from=-1, to=6)
plot(X, add=TRUE, type='p', from=-1, to=6)

## ----comput1a,fig.keep='none'------------------------------------------
A <- FuzzyNumber(-5, 3, 6, 20,
    left=function(x) pbeta(x,0.4,3),
   right=function(x) 1-x^(1/4),
   lower=function(alpha) qbeta(alpha,0.4,3),
   upper=function(alpha) (1-alpha)^4
)

## ----comput1b,dependson='comput1a'-------------------------------------
supp(A)

## ----comput1c,dependson='comput1b'-------------------------------------
core(A)

## ----comput1d,dependson='comput1c'-------------------------------------
alphacut(A, 0) # same as supp(A) (if alpha-cut generators are defined)
alphacut(A, 1) # same as core(A)
(a <- alphacut(A, c(0, 0.5, 1)))
a[1, ]
a[2, 2]
a[, "L"]

## ----comput1e,dependson='comput1d'-------------------------------------
evaluate(A, 1)
evaluate(A, c(-3,0,3))
evaluate(A, seq(-1, 2, by=0.5))

## ----comput1f,dependson='comput1e'-------------------------------------
expectedInterval(A)

## ----comput1g,dependson='comput1f'-------------------------------------
expectedValue(A)

## ----comput1h,dependson='comput1g'-------------------------------------
weightedExpectedValue(A, 0.5) # equivalent to expectedValue(A)
weightedExpectedValue(A, 0.25)

## ----comput1i,dependson='comput1h'-------------------------------------
value(A)

## ----comput1j,dependson='comput1i'-------------------------------------
width(A)

## ----comput1k,dependson='comput1j'-------------------------------------
ambiguity(A)

## ----comput1l,dependson='comput1k'-------------------------------------
diff(supp(A))

## ----addition,fig.keep='none'------------------------------------------
A <- TrapezoidalFuzzyNumber(0, 1, 1, 2)
B <- TrapezoidalFuzzyNumber(1, 2, 2, 3)
plot(A, xlim=c(0,6))
plot(B, add=TRUE, col=2, lty=2)
plot(A+B, add=TRUE, col=4, lty=4)

## ----additionFig,dependson='addition',echo=FALSE,results='hide'--------
par(mar=c(4,4,2,1))
plot(A, xlim=c(0,6), xlab=expression(x), ylab=expression(alpha))
plot(B, add=TRUE, col=2, lty=2)
plot(A+B, add=TRUE, col=4, lty=4)
legend('topright', expression(A, B, A+B), lty=c(1,2,4), col=c(1,2,4))

## ----ops2--------------------------------------------------------------
A <- piecewiseLinearApproximation(PowerFuzzyNumber(1,2,3,4,p.left=2,p.right=0.5),
   method="Naive", knot.n=20)
B <- piecewiseLinearApproximation(PowerFuzzyNumber(2,3,4,5,p.left=0.1,p.right=3),
   method="Naive", knot.n=40)
A+A # the same as 2*A
A+B # note the number of knots has increased

## ----fapply,fig.keep='none'--------------------------------------------
A <- as.PiecewiseLinearFuzzyNumber(TrapezoidalFuzzyNumber(0,1,2,3), knot.n=100)
plot(fapply(A, function(x) sqrt(log(x+1))))

## ----fapplyFig,dependson='fapply',echo=FALSE,results='hide'------------
par(mar=c(4,4,2,1))
plot(fapply(A, function(x) log(x+1)^0.5), xlab=expression(x), ylab=expression(alpha))
legend('topleft', expression(sqrt(log(A+1))), lty=1)

## ----exponential,fig.keep='none'---------------------------------------
A <- as.PiecewiseLinearFuzzyNumber(TrapezoidalFuzzyNumber(-2,-1,-1,2), knot.n=10)
plot(A, xlim=c(-8,8))
plot(A^2, add=TRUE, col=2, lty=2)
plot(A^3, add=TRUE, col=4, lty=4)

## ----exponentialFig,dependson='exponential',echo=FALSE,results='hide'----
par(mar=c(4,4,2,1))
plot(A, xlim=c(-8,8), xlab=expression(x), ylab=expression(alpha))
plot(A^2, add=TRUE, col=2, lty=2)
plot(A^3, add=TRUE, col=4, lty=4)
legend('topright', expression(A, A^2, A^3), lty=c(1,2,4), col=c(1,2,4))

## ----TrapDist----------------------------------------------------------
T1 <- TrapezoidalFuzzyNumber(-5, 3, 6, 20)
T2 <- TrapezoidalFuzzyNumber(-4, 4, 7, 21)
distance(T1, T2, type='Euclidean') # L2 distance /default/
distance(T1, T2, type='EuclideanSquared') # Squared L2 distance

## ----ApproxExA---------------------------------------------------------
A <- FuzzyNumber(-5, 3, 6, 20,
   left=function(x) pbeta(x,0.4,3),
   right=function(x) 1-x^(1/4),
   lower=function(alpha) qbeta(alpha,0.4,3),
   upper=function(alpha) (1-alpha)^4
)

## ----ApproxExA_naive,dependson='ApproxExA',fig.keep='none'-------------
(T1 <- trapezoidalApproximation(A, method='Naive'))
distance(A, T1)

## ----ApproxExA_naive2,dependson='ApproxExA_naive',echo=FALSE-----------
par(mar=c(4,4,2,1))
plot(A, xlab=expression(x), ylab=expression(alpha))
plot(T1, col='red', lty=2, add=TRUE)
legend('topright', legend=expression(A, 'Naive approx.'),
   col=c('black', 'red'), lty=c(1,2))

## ----ApproxExA_L2n,dependson='ApproxExA',fig.keep='none'---------------
(T2 <- trapezoidalApproximation(A, method='NearestEuclidean'))
distance(A, T2)

## ----ApproxExA_naiveL2n2,dependson='ApproxExA_L2n',echo=FALSE----------
par(mar=c(4,4,2,1))
plot(A, xlab=expression(x), ylab=expression(alpha))
plot(T2, col='red', lty=2, add=TRUE)
legend('topright', legend=expression(A, L[2]-"nearest approx."),
   col=c('black', 'red'), lty=c(1,2))

## ----ApproxExA_ExpInt,dependson='ApproxExA',fig.keep='none'------------
ambiguity(A)
width(A)/3
(T3 <- trapezoidalApproximation(A, method='ExpectedIntervalPreserving'))
distance(A, T3)
expectedInterval(A)
expectedInterval(T3)

## ----ApproxExA_ExpInt2,dependson='ApproxExA_ExpInt',eval=FALSE,echo=FALSE----
#  par(mar=c(4,4,2,1))
#  plot(A, xlab=expression(x), ylab=expression(alpha))
#  plot(T3, col='red', lty=2, add=TRUE)

## ----ApproxExA_ExpInt3,fig.keep='none'---------------------------------
(B  <- FuzzyNumber(1, 2, 3, 45,
   lower=function(x) sqrt(x),
   upper=function(x) 1-sqrt(x)))
(TB1 <- trapezoidalApproximation(B, 'NearestEuclidean'))
(TB2 <- trapezoidalApproximation(B, 'ExpectedIntervalPreserving'))
distance(B, TB1)
distance(B, TB2)

## ----ApproxExA_ExpInt3Fig,dependson='ApproxExA_ExpInt3',echo=FALSE-----
par(mar=c(4,4,2,1))
plot(B, xlab=expression(x), ylab=expression(alpha), log='x', xlim=c(0.9,46))
plot(TB1, col='red',  lty=2, add=TRUE)
plot(TB2, col='blue', lty=3, add=TRUE)
legend('topright', expression(mu[B], mu[TB[1]], mu[TB[2]]),
   col=c('black', 'red', 'blue'), lty=c(1,2,3))

## ----ApproxExA_RestrSuppCore,dependson='ApproxExA',fig.keep='none'-----
(T4 <- trapezoidalApproximation(A, method='SupportCoreRestricted'))
distance(A, T4)

## ----ApproxExA_RestrSuppCore2,dependson='ApproxExA_RestrSuppCore',echo=FALSE----
par(mar=c(4,4,2,1))
plot(A, xlab=expression(x), ylab=expression(alpha))
plot(T4, col='red', lty=2, add=TRUE)
legend('topright', legend=expression(A, 'Supp&core restr.'),
   col=c('black', 'red'), lty=c(1,2))

## ----ApproxExAPLFN-----------------------------------------------------
A <- FuzzyNumber(-5, 3, 6, 20,
   left=function(x) pbeta(x,0.4,3),
   right=function(x) 1-x^(1/4),
   lower=function(alpha) qbeta(alpha,0.4,3),
   upper=function(alpha) (1-alpha)^4
)

## ----ApproxPLFNNaive,dependson='ApproxExAPLFN',fig.keep='none'---------
P1 <- piecewiseLinearApproximation(A, method='Naive',
         knot.n=1, knot.alpha=0.5)
P1['allknots']
print(distance(A, P1), 8)

## ----ApproxPLFNNaive2,dependson='ApproxPLFNNaive',echo=FALSE-----------
par(mar=c(4,4,2,1))
plot(A, xlab=expression(x), ylab=expression(alpha))
plot(P1, col='red', lty=2, add=TRUE)
legend('topright', legend=expression(A, 'Naive approx.'),
   col=c('black', 'red'), lty=c(1,2))

## ----ApproxPLFNNearest,dependson='ApproxExAPLFN',fig.keep='none'-------
P2 <- piecewiseLinearApproximation(A,
   method='NearestEuclidean', knot.n=3, knot.alpha=c(0.25,0.5,0.75))
print(P2['allknots'], 6)
print(distance(A, P2), 12)

## ----ApproxPLFNNearest2,dependson='ApproxPLFNNearest',echo=FALSE-------
par(mar=c(4,4,2,1))
plot(A, xlab=expression(x), ylab=expression(alpha))
plot(P2, col='red', lty=2, add=TRUE)
legend('topright', legend=expression(A, L[2]-"nearest approx."),
   col=c('black', 'red'), lty=c(1,2))

## ----ApproxPLFNNearestConv1,dependson='ApproxPLFNNearest2',fig.keep='none'----
n <- 1:27
d <- matrix(NA, ncol=4, nrow=length(n))
# d[,1] - Naive approximator's error for given knot.n
# d[,2] - Best L2 approximator's error
# d[,3] - theoretical upper bound
for (i in seq_along(n))
{
   P1 <- piecewiseLinearApproximation(A, method='Naive',
            knot.n=n[i]) # equidistant knots
   P2 <- piecewiseLinearApproximation(A, method='NearestEuclidean',
            knot.n=n[i]) # equidistant knots
   d[i,1] <- distance(A, P1)
   d[i,2] <- distance(A, P2)

   acut <- alphacut(A, seq(0, 1, length.out=n[i]+2))
   # d[i,3] <- sqrt(sum((c(diff(acut[,1]), diff(acut[,2]))^2)/(n[i]+1))) # beter ubound
   d[i,3] <- sqrt(2)*max(abs(c(diff(acut[,1]), diff(acut[,2]))))
}
matplot(n, d, type='l', log="y", lty=c(1,2,4), col=c(1,2,4))

## ----ApproxPLFNNearestConv2,dependson='ApproxPLFNNearestConv1',echo=FALSE----
par(mar=c(4,4,2,1))
matplot(n, d, type='l', xlab=expression(n), ylab=expression(d(A, '*')), log="y", lty=c(1,2,4), col=c(1,2,4))
legend(x=19, y=3.2, expression(d(A, {N^n}(A)), d(A, {Pi^n}(A)), 'upper bound'),
   lty=c(1,2,4), col=c(1,2,4))

## ----ApproxPLFNNearest3,dependson='ApproxPLFNNearest2',fig.keep='none'----
a <- seq(1e-9, 1-1e-9, length.out=100) # many alphas from (0,1)
d <- numeric(length(a)) # distances /to be calculated/
for (i in seq_along(a))
{
   P1 <- piecewiseLinearApproximation(A, method='NearestEuclidean',
            knot.n=1, knot.alpha=a[i])
   d[i] <- distance(A, P1)
}
plot(a, d, type='l', xlab=expression(alpha), ylab=expression(D[A](alpha)))

## ----ApproxPLFNNearest4,dependson='ApproxPLFNNearest3',echo=FALSE------
par(mar=c(4,4,2,1))
plot(a, d, type='l', xlab=expression(alpha), ylab=expression(D[A](alpha)))

## ----ApproxPLFNNearest5,dependson='ApproxPLFNNearest4',warning=FALSE----
for (i in 1:5) # 5 iterations
{
   a0 <- runif(1,0,1) # random starting point
   optim(a0,
      function(a)
      {
         P1 <- piecewiseLinearApproximation(A, method='NearestEuclidean',
                                             knot.n=1, knot.alpha=a)

         distance(A, P1)
      }, method='L-BFGS-B', lower=1e-9, upper=1-1e-9) -> res
   cat(sprintf('%.9f %6g *%.9f* %.9f\n', a0, res$counts[1], res$par, res$value))
}

## ----suppcoreplfn1a----------------------------------------------------
A <- FuzzyNumber(0, 3, 4, 5,
   lower=function(x) qbeta(x, 2, 1),
   upper=function(x) 1-x^3
)

## ----suppcoreplfn1a2,dependson='suppcoreplfn1a'------------------------
knot.alpha <- 0.2
P1 <- piecewiseLinearApproximation(A, knot.alpha=knot.alpha)
P2 <- piecewiseLinearApproximation(A, method="SupportCorePreserving",
   knot.alpha=knot.alpha)
distance(A, P1)
print(alphacut(P1, c(0, knot.alpha, 1)))
distance(A, P2)
print(alphacut(P2, c(0, knot.alpha, 1)))

## ----suppcoreplfn1a3,dependson='suppcoreplfn1a2',echo=FALSE------------
par(mar=c(4,4,2,1))
plot(A, xlab=expression(x), ylab=expression(alpha))
plot(P1, col=2, add=TRUE, lty=2)
plot(P2, col=4, add=TRUE, lty=4, lwd=2)
abline(h=knot.alpha, col="grey", lty=3)
legend("topleft", legend=expression(A, P[1], P[2]), lty=c(1,2,4), col=c(1,2,4), lwd=c(1,1,2))

## ----suppcoreplfn1a4,dependson='suppcoreplfn1a'------------------------
D <- function(a) distance(A,
   piecewiseLinearApproximation(A, method="SupportCorePreserving", knot.alpha=a))
optimize(D, lower=0, upper=1)

## ----suppcoreplfn1a5,dependson='suppcoreplfn1a4',echo=FALSE------------
par(mar=c(4,4,2,1))
a <- seq(0, 1, length.out=1001)
Da <- sapply(a, D)
plot(a[-c(1,length(a))], Da[-c(1,length(a))], xlab=expression(alpha), ylab=expression(D(alpha)),
   type='l', xlim=c(0,1), ylim=range(Da))
points(c(0, 1), Da[c(2, length(a)-1)])
points(c(0, 1), Da[c(1,length(a))], pch=16)
opt <- optimize(D, interval=c(0,1))
abline(v=opt$minimum, col="grey", lty=3)
abline(h=opt$objective, col="grey", lty=3)

## ----compareXY1--------------------------------------------------------
x = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(0.2, 1.0, 2.8))
y = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(0, 1.8, 2.2))

## ----compareXY1plot,dependson='compareXY1'-----------------------------
plot(x, col=2)
plot(y, col=4, add=TRUE)

## ----compareXY1c1,dependson='compareXY1'-------------------------------
possibilityExceedance(x,y)
necessityExceedance(x,y)

## ----compareXY1c2,dependson='compareXY1'-------------------------------
possibilityStrictExceedance(x,y)
necessityStrictExceedance(x,y)

## ----compareXY1c3,dependson='compareXY1'-------------------------------
possibilityUndervaluation(x,y)
necessityUndervaluation(x,y)
possibilityStrictUndervaluation(x,y)
necessityStrictUndervaluation(x,y)

## ----compareXY2--------------------------------------------------------
x = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(1.7, 2.7, 2.8), knot.n = 9)
y = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(0, 1.8, 2.2), knot.n = 9)
min = min(x@a1,y@a1)
max = max(x@a4,y@a4)
plot(x, col=2, xlim = c(min,max))
plot(y, col=4, add=TRUE)
possibilityExceedance(x,y)
necessityExceedance(x,y)
possibilityStrictExceedance(x,y)
necessityStrictExceedance(x,y)
possibilityUndervaluation(x,y)
necessityUndervaluation(x,y)
possibilityStrictUndervaluation(x,y)
necessityStrictUndervaluation(x,y)

## ----minmaxXY,fig.keep='none'------------------------------------------
x = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-4.8, -3 , -1.5), knot.n = 9)
y = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-5.5, -2.5, -1.1), knot.n = 9)
min = min(x@a1,y@a1)
max = max(x@a4,y@a4)
plot(x, col=1, xlim = c(min,max))
plot(y, col=2, add=TRUE)
maxFN = maximum(x,y)
minFN = minimum(x,y)
plot(minFN, col=4)
plot(maxFN, col=6, add=TRUE)

## ----minmaxXYFig,dependson='minmaxXY',echo=FALSE,results='hide'--------
par(mar=c(4,4,2,1))
plot(x, col=1, xlim = c(min,max))
plot(y, col=2, add=TRUE)
legend('topright', expression(x, y), lty=c(1,1), col=c(1,2))
plot(minFN, col=4)
plot(maxFN, col=6, add=TRUE)
legend('topright', expression(min, max), lty=c(1,1), col=c(4,6))

