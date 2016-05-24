##
##  q u a d r a t u r e . R  Test suite
##


quad <- pracma::quad
quadl <- pracma::quadl
quadgk <- pracma::quadgk
quadgr <- pracma::quadgr
quadinf <- pracma::quadinf
quad2d <- pracma::quad2d
dblquad <- pracma::dblquad

simpson2d <- pracma::simpson2d
simpadpt <- pracma::simpadpt
gauss_kronrod <- pracma::gauss_kronrod
clenshaw_curtis <- pracma::clenshaw_curtis
romberg <- pracma::romberg

gaussLegendre <- pracma::gaussLegendre
gaussHermite <- pracma::gaussHermite
gaussLaguerre <- pracma::gaussLaguerre

##  F i n i t e  I n t e r v a l s

f1 <- function(x) exp(x)*sin(x)      # [0, pi]  12.0703463163896 = 1/2*(1+e^pi)
f2 <- pracma::runge                  # [-1, 1]   0.549360306778006
f3 <- function(x) 1/(x^3 - 2*x - 5)  # [0, 2]   -0.460501533846733
f4 <- function(x) abs(sin(10*x))     # [0, pi]   2.0

# quad (Adaptive Simpson)
all.equal(quad(f1, 0, pi, tol=1e-12),   12.0703463163896,
          tolerance = 1e-12)
all.equal(quad(f2, -1, 1, tol=1e-12),    0.549360306778006,
          tolerance = 1e-12)
all.equal(quad(f3, 0, 2, tol=1e-12),    -0.460501533846733,
          tolerance = 1e-12)
all.equal(quad(f4, 0, pi, tol=1e-12),    2.0,
          tolerance = 1e-12)

# quadl (Adaptive Lobatto)
all.equal(quadl(f1, 0, pi, tol=1e-9),   12.0703463163896,
          tolerance = 1e-12)
all.equal(quadl(f2, -1, 1, tol=1e-9),    0.549360306778006,
          tolerance = 1e-12)
all.equal(quadl(f3, 0, 2, tol=1e-9),    -0.460501533846733,
          tolerance = 1e-12)
all.equal(quadl(f4, 0, pi, tol=1e-12),   2.0,
          tolerance = 1e-12)

# quadgr (Gauss-Richardson)
all.equal(quadgr(f1, 0, pi, tol=1e-12)$value, 12.0703463163896,
          tolerance = 1e-13)
all.equal(quadgr(f2, -1, 1, tol=1e-12)$value,  0.549360306778006,
          tolerance = 1e-15)
all.equal(quadgr(f3, 0, 2, tol=1e-12)$value,  -0.460501533846733,
          tolerance = 1e-15)
all.equal(quadgr(f4, 0, pi, tol=1e-12)$value,  2.0,
          tolerance = 1e-15)

# quadgk (Adaptive Gauss-Kronrod)
all.equal(quadgk(f1, 0, pi), 12.0703463163896,
          tolerance = 1e-13)
all.equal(quadgk(f2, -1, 1),  0.549360306778006,
          tolerance = 1e-15)
all.equal(quadgk(f3, 0, 2),  -0.460501533846733,
          tolerance = 1e-13)
all.equal(quadgk(f4, 0, pi, tol = 1e-12),  2.0,
          tolerance = 1e-12)

# Adaptive Simpson (simpadpt)
all.equal(simpadpt(f1, 0, pi, tol=1e-12),     12.0703463163896,
          tolerance = 1e-13)
all.equal(simpadpt(f2, -1, 1, tol=1e-12),      0.549360306778006,
          tolerance = 1e-12)
all.equal(simpadpt(f3, 0, 2, tol=1e-12),      -0.460501533846733,
          tolerance = 1e-13)
all.equal(simpadpt(f4, 0, pi, tol=1e-12),      2.0,
          tolerance = 1e-14)

# Gauss-Kronrod
all.equal(gauss_kronrod(f1, 0, pi)$value, 12.0703463163896,
          tolerance = 1e-13)
all.equal(gauss_kronrod(f2, -1, 1)$value,  0.549360306778006,       # BAD
          tolerance = 1e-2)
all.equal(gauss_kronrod(f3, 0, 2)$value,  -0.460501533846733,       # Bad
          tolerance = 1e-5)
all.equal(gauss_kronrod(f4, 0, pi)$value,  2.0,                     # BAD
          tolerance = 1e-0)

# Clenshaw-Curtis
all.equal(clenshaw_curtis(f1, 0, pi, n = 128), 12.0703463163896,
          tolerance = 1e-12)
all.equal(clenshaw_curtis(f2, -1, 1, n = 128),  0.549360306778006,
          tolerance = 1e-12)
all.equal(clenshaw_curtis(f3, 0, 2, n = 128),  -0.460501533846733,
          tolerance = 1e-12)
all.equal(clenshaw_curtis(f4, 0, pi, n = 1024),  2.0,               # Bad
          tolerance = 2e-5)

# romberg
all.equal(romberg(f1, 0, pi, tol=1e-12)$value, 12.0703463163896,
          tolerance = 1e-12)
all.equal(romberg(f2, -1, 1, tol=1e-12)$value,  0.549360306778006,
          tolerance = 1e-12)
all.equal(romberg(f3, 0, 2, tol=1e-12)$value,  -0.460501533846733,  # BAD
          tolerance = 1e-3)
all.equal(romberg(f4, 0, pi, tol=1e-12)$value,  2.0,
          tolerance = 1e-12)

f5 <- function(x) log(x)*sin(x)/x       # pi/2 * gamma , cannot be computed !
f6 <- function(x) sin(x)^2 * exp(-x)    # [0, Inf] , 0.4
f7 <- function(x) sin(x)^2 * exp(-x^2)  # [-Inf, Inf] , (e-1)*sqrt(pi)/(4*e)
x7 <- (exp(1)-1) * sqrt(pi) / (2*exp(1))

# quadinf
all.equal(quadinf(f6, 0, Inf), 0.4, tolerance = 1e-15)
all.equal(quadinf(f7, -Inf, Inf), x7, tolerance = 1e-15)

all.equal(quadgr(f6, 0, Inf)$value, 0.4, tolerance = 1e-11)
all.equal(quadgr(f7, -Inf, Inf)$value, x7, tolerance = 1e-9)

gL <- gaussLaguerre(64)
all.equal(sum(gL$w * sin(gL$x)^2), 0.4, tolerance = 1e-15)
gH <- gaussHermite(64)
all.equal(sum(gH$w * sin(gH$x)^2), x7, tolerance = 1e-14)

f8 <- function(x, y) y * sin(x)  # [0, pi/2]x[0, 1] , 1/2
f9 <- function(x, y) ifelse(x^2 + y^2 <= 1, 1-x^2-y^2, 0)

# quad2d
all.equal(quad2d(f8, 0, pi/2, 0, 1), 0.5,         tolerance = 1e-15)
all.equal(quad2d(f9, -1, 1, 0, 1, n = 128), pi/4, tolerance = 1e-6)

# dblquad
all.equal(dblquad(f8, 0, pi/2, 0, 1), 0.5, tolerance = 1e-15)
#all.equal(dblquad(f9, -1, 1, 0, 1), pi/4,  tolerance = 1e-6)
    # disabled because of problems with Fedora and Solaris

# simpson2d
all.equal(simpson2d(f8, 0, pi/2, 0, 1), 0.5, tolerance = 1e-9)
all.equal(simpson2d(f9, -1, 1, 0, 1), pi/4,  tolerance = 1e-5)

# Integrals with singularities at boundaries:
f11 <- function(t) log(1-t) / t  # [1, 0]  pi^2/6 , dilogarithm
f12 <- function(t) log(-log(t))  # [0, 1]  gamma = 0.57721 56649 01532 ...
f13 <- function(t) 1 / sqrt(t)   # [0, 1]  2.0

all.equal(quad(f11, 1, 0, tol = 1e-12),  1.64493406684823,
          tolerance = 1e-10)
all.equal(quad(f12, 0, 1, tol = 1e-12), -0.577215664901533,
          tolerance = 5e-10)
all.equal(quad(f13, 0, 1, tol = 1e-12),  2.0,
          tolerance = 1e-4)                                         # Bad
                                       
all.equal(quadl(f11, 1, 0, tol = 1e-12),  1.64493406684823,
          tolerance = 1e-12)           
all.equal(quadl(f12, 0, 1, tol = 1e-12), -0.577215664901533,
          tolerance = 5e-12)           
all.equal(quadl(f13, 0, 1, tol = 1e-12),  2.0,
          tolerance = 1e-7)                                         # Bad

all.equal(quadgr(f11, 1, 0, tol = 1e-12)$value,  1.64493406684823,
          tolerance = 1e-12)
all.equal(quadgr(f12, 0, 1, tol = 1e-12)$value, -0.577215664901533,
          tolerance = 5e-12)
all.equal(quadgr(f13, 0, 1, tol = 1e-12)$value,  2.0,
          tolerance = 1e-12)

all.equal(simpadpt(f11, 1, 0, tol = 1e-12),  1.64493406684823,
          tolerance = 1e-11)           
all.equal(simpadpt(f12, 0, 1, tol = 1e-12), -0.577215664901533,
          tolerance = 5e-11)           
all.equal(simpadpt(f13, 0, 1, tol = 1e-10),  2.0,
          tolerance = 1e-7)                                         # Bad

##  E o F
