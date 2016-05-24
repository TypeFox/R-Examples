 times = seq(0,20,0.5)
 FHN.xnames = c('V','R')
 FHN.fn = make.fhn()
 FHN.knots  = seq(0,20,2)
 FHN.norder  = 3
 FHN.nbasis = length(FHN.knots) + FHN.norder - 2
 FHN.range  = c(0,20)
 FHN.basis  = create.bspline.basis(FHN.range, FHN.nbasis,
                                   FHN.norder, FHN.knots)
 FHN.fdnames = list(NULL, NULL, FHN.xnames)
 FHN.data = matrix(c(
   -0.9508,    0.1941,
   -0.1886,    1.0975,
    1.4686,   -0.0005,
    1.9981,    1.1561,
    2.0249,    0.7667,
    2.0081,    0.0546,
    1.4857,   -0.2347,
    1.4184,   -1.6662,
    2.3930,   -0.3252,
    0.0610,   -1.8704,
    2.0342,   -2.0501,
    0.5523,   -0.9018,
   -0.7991,   -1.2990,
   -2.7791,   -0.2625,
   -2.1212,    0.2284,
   -1.8151,    0.6424,
   -1.2904,    0.1545,
   -2.1196,    0.9703,
   -0.7323,    0.9908,
   -0.8960,    1.4679,
    1.7302,    1.2176,
    2.3323,    1.0775,
    2.0664,    0.2541,
    2.3299,   -0.0338,
    2.1683,    0.2920,
    1.2035,   -1.4959,
    1.5007,   -0.8775,
    0.7062,   -1.4168,
    0.2399,   -1.0865,
    0.7983,   -0.6172,
   -1.4062,   -0.3967,
   -1.9705,   -0.9559,
   -1.3631,   -0.3275,
   -1.3702,    0.2982,
   -1.3165,    0.3683,
   -0.6459,    0.9082,
   -0.4809,    1.1491,
   -0.1056,    0.5921,
    1.4264,    0.8459,
    1.6757,   -0.4546,
    2.4930,    0.8772),41,2,byrow=TRUE)
 colnames(FHN.data) = c("V","R")
 FHN.xfd   = smooth.basis(times, FHN.data, FHN.basis,
                          fdnames=FHN.fdnames)$fd
 FHN.coefs = FHN.xfd$coefs
 FHN.pars0        = c(0.2, 0.2, 3.0)
 names(FHN.pars0) = c("a", "b", "c")
 FHN.lambda       = 1e4*c(1,1)
 FHN.x0           = c(-1,1)
 profile.obj = LS.setup(FHN.pars0,FHN.coefs,FHN.fn,FHN.basis,FHN.lambda,
                        times=times,names=c("V","R"))
 lik   = profile.obj$lik     #  lik object
 proc  = profile.obj$proc    #  proc object
 
f      = SplineCoefsErr(coefs,times,FHN.data,lik,proc,FHN.pars0)

dfdc   = SplineCoefsDC(coefs,times,FHN.data,lik,proc,FHN.pars0)

d2fdc2 = SplineCoefsDC2(coefs,times,FHN.data,lik,proc,FHN.pars0)

d2fdcdp = SplineCoefsDCDP(coefs,times,FHN.data,lik,proc,FHN.pars0)

#  Profile.LS argument correspondences

data      = FHN.data
fn        = FHN.fn
pars      = FHN.pars0
coefs     = FHN.coefs
coefs.0   = FHN.coefs
basisvals = FHN.basis
lambda    = FHN.lambda
fdobj       = NULL
more        = NULL
weights     = NULL 
quadrature  = NULL
active      = NULL
in.method   = NULL
out.method  = NULL
options.in  = NULL
options.out = NULL
eps         = 1e-6
poslik      = 0 
posproc     = 0 
discrete    = 0 
sgn         = 1


profile.obj = Profile.LS(FHN.fn,FHN.data,times,FHN.pars0,FHN.coefs,basisvals=FHN.basis,FHN.lambda)
