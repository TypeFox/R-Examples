# 24-12-2013 MRC-Epid JHZ

# The following is from experiments with GCTA on heritability estimates

VG <- 0.017974; VVG <- 0.003988^2
VGE <- 0.002451; VVGE <- 0.005247^2
Ve <- 0.198894; VVe <- 0.005764^2

cVGVGE <- -7.93348e-06
cVGVe <- -5.54006e-06
cVGEVe <- -1.95297e-05

Vp <- VG + VGE + Ve
s <- VVG + VVGE + VVe
s12 <- 2*cVGVGE
s13 <- 2*cVGVe
s23 <- 2*cVGEVe
VVp <- s + s12 + s13 + s23
cVpVG <- VVG + cVGVGE + cVGVe
cVpVGE <- cVGVGE + VVGE + cVGEVe

library(gap)

Vp
sqrt(VVp)
VG/Vp
sqrt(VR(VG, VVG, Vp, VVp, cVpVG))
VGE/Vp
sqrt(VR(VGE,VVGE,Vp,VVp, cVpVGE))

K <- 0.05
x <- qnorm(1-K)
z <- dnorm(x)
1/sqrt(2*pi)*exp(-x^2/2)
P <- 0.496404
fK <- (K*(1-K)/z)^2
fP <- P*(1-P)
f <- fK/fP
ho <- 0.274553
vo <- 0.067531
f*ho
f*vo
hl <- 0.232958
vl <- 0.057300^2
r1 <- hl/ho
r2 <- vl/vo
r1==r2
z2 <- K^2*(1-K)^2/(f*fP)
x2 <- -log(2*pi*z2)
sqrt(x2)

V <- c(0.017974, 0.002451, 0.198894)
VCOV <- matrix(0,3,3)
diag(VCOV) <- c(0.003988, 0.005247, 0.005764)^2
VCOV[2,1] <- -7.93348e-06
VCOV[3,1] <- -5.54006e-06
VCOV[3,2] <- -1.95297e-05
z <- h2GE(V,VCOV)

h2 <- 0.274553
se <- 0.067531
P <- 0.496404
hl <- 0.232958
vl <- 0.057300^2

z <- h2l(P=P,h2=h2,se=se)
