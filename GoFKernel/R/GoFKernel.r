# GoFKernel version 2.1-0
# January 2016 - Jose M. Pavia
# Tested in R version 3.3.2
# The  fan.test function depends on KernSmooth package
#
# FUNCTIONS:
#      inverse
#      random.function
#      area.between
#      support.facto
#      density.reflected
#      dgeometric.test
#      fan.test
#
#      inverse: computes the inverse function of any (cumulative distribution) function.
#      random.function: generates draws of any random variable given its distribution function.
#      area.between: calculates the area between a theoretical density and a kernel estimate.
#      support.facto: computes the numerical limits of a density with infinite support.
#      density.reflected: computes a kernel estimate using reflection in the borders.
#      dgeometric.test: test to  compare geometrically a theoretical distribution with an empirical kernel density one.
#      fan.test: Implementation of the Fan's test.
#
# DATA:
#      risk76.1929
#      A vector containing the time exposed to the risk of death with 76 years during
#      2006 for the 2006 registered Spanish immigrants born in 1929.
#      Under the null hypothesis of uniform distribution of date of birth and date of
#      migration, this time exposed to risk is distributed as a f(x)=2-2x 0<x<1.


inverse <- function (f, lower=-Inf, upper=Inf) {
# DESCRIPTION: Function to calculate the inverse function of a function
#
# INPUTS:
#       f: function for which we want to obtain its inverse
#       lower: lower limit of f domain (support of the random variable)
#       upper: upper limit of f domain (support of the random variable)
#
# OUTPUT:
#       A function, the inverse of f
    if (!is.numeric(lower) || !is.numeric(upper) || lower >=upper)
       stop("lower < upper is not fulfilled")

    if (!is.finite(f(lower)) & is.finite(lower)) lower <- lower+.Machine$double.eps
    if (!is.finite(f(upper)) & is.finite(upper)) upper <- upper-.Machine$double.eps

    if (is.infinite(lower) & is.infinite(upper)) {
       function(y) optim(0,(function (x) abs(f(x) - y)), method="L-BFGS-B")$par
   } else if (is.infinite(lower) & is.finite(upper)) {
       function(y) optim(upper,(function (x) abs(f(x) - y)), lower = lower, upper = upper, method="L-BFGS-B")$par
   } else if (is.finite(lower) & is.infinite(upper)) {
       function(y) optim(lower,(function (x) abs(f(x) - y)), lower = lower, upper = upper, method="L-BFGS-B")$par
   } else {
      if (f(lower)<f(upper)) {
          function(y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)$root
      } else {
          function(y) optim((lower+upper)/2,(function (x) (f(x) - y)^2), lower = lower, upper = upper, method="L-BFGS-B")$par
      }
   }
}
#................................................


random.function <- function (n=1, f, lower=-Inf, upper=Inf, kind="density") {
# DESCRIPTION: Function to generate random draws of a continuos random variable
#              given its density or cumulative distribution function
#
# INPUTS:
#        n: number of draws
#        f: density function (default) of the random variable
#        lower: lower limit of the support of the random variable
#        upper: ipper limit of the support of the random variable
#        kind:  character string with the kind of function used to identify the distribution.
#               either "density" or "cumulative" as alternative
#
# OUTPUT:
#        A vector of random draws of f
#
   if (!is.finite(f(lower))) lower <- lower+.Machine$double.eps
   if (!is.finite(f(upper))) upper <- upper-.Machine$double.eps

   if (!is.numeric(lower) || !is.numeric(upper) || lower >=upper)
       stop("lower < upper is not fulfilled")
   cumulative <- function (f, lower) {
       function(z) integrate(f,lower,upper=z,rel.tol=1.e-10,subdivisions=1000)$value
   }
   if (kind=="density"){
      f.distr <- cumulative(f=f,lower=lower)
      limits <- support.facto(f,lower=lower,upper=upper)
      lower <- max(c(lower,limits[1]))
      upper <- min(c(upper,limits[2]))
   } else {
      f.distr <- f
   }
#   if (f.distr(lower)!=0)
#       stop("f.distr(lower) must be 0")
#   if (f.distr(upper)!=1)
#       stop("f.distr(upper) must be 1")
   cum.inverse <- inverse(f.distr, lower=lower, upper=upper)
   sample.unif <- runif(n)
   sample.function <- sapply(sample.unif,cum.inverse)
   return(sample.function)
}
#................................................


support.facto <- function(f,lower=-Inf,upper=Inf){
# DESCRIPCION: For a density function of theoretical infinite support
#              computes its de facto numerical limits
#
1# INPUT:
#       f: a density function. It requires that its first two ordinary moments exist.
#          Otherwise, the limits are not modified.
#       lower: theoretical lower limit of the support of the random variable
#       upper: theoretical lower limit of the support of the random variable
#
# OUTPUT:
#       A two components vector with the de facto lower and upper limits

  f.mean<-try(integrate(function(x) x*f(x), lower=lower, upper=upper, subdivisions=10000)$value,silent=T)
  f.x2<-try(integrate(function(x) x^2*f(x), lower=lower, upper=upper, subdivisions=10000)$value, silent=T)
  if (is.numeric(f.x2) & f.x2!=0){
     f.sd<-sqrt(f.x2-f.mean^2)
     k.f<-sqrt(1/10^-5) # Chebyshev
     k.min<-sqrt(1/10^-21) # Chebyshev

     g<-function(y){ if (f(y)<10^(-20) & y<f.mean) { f.mean-y } else { 10^50 } }
     lower.facto<-optim(f.mean-k.f*f.sd,g, lower=f.mean-k.min*f.sd, upper=f.mean, method="Brent")$par

     g<-function(y){ if (f(y)<10^(-20) & y>f.mean) { y-f.mean } else { 10^50 } }
     upper.facto<-optim(f.mean+k.f*f.sd, g, lower=f.mean, upper=f.mean+k.min*f.sd, method="Brent")$par
   } else {
     lower.facto <- lower
     upper.facto <- upper
   }
 output<-c(lower.facto,upper.facto)
 return(output)
}
#................................................


area.between <- function(f, kernel.density, lower=-Inf, upper=Inf){
# DESCRIPTION: function to calculate the area between a theoretical density and
#              a kernel empirical estimate in a given interval
#
# INPUTS:
#       f: theorical density funtion
#       kernel.density: empirical kernel estimate, object of the class density
#       lower: lower limit of the support of f
#       upper: upper limit of the support of f
#
# OUTPUT:
#       A number with the numerical value of the area between both functions

  bw<-kernel.density$x[2]-kernel.density$x[1] # bandwidth
  selected<-which(kernel.density$x>=lower & kernel.density$x<=upper)
  area.bt<-sum(abs(kernel.density$y[selected]-f(kernel.density$x[selected]))*bw)
  return(area.bt)
}
#................................................


density.reflected <- function(x, lower=-Inf, upper=Inf, weights=NULL, ...){
# DESCRIPTION: function to compute from a random sample of data in an interval
#              a kernel estimate using reflection in the borders
#
# INPUTS:
#       x: a numeric vector of data values
#       lower: lower limit of the interval to which x belongs to
#       upper: upper limit of the interval to which x belongs to
#       ...: further density arguments.
#
# OUTPUT:
#       Object of the class density with borders correction
#
# OBSERVATIONS:
#       To estimate the density the default option of fucntion density is used
#
mantener<- !is.na(x)
x <- x[mantener]
if (upper < max(x)) warning("There are values in the sample higher than the upper limit")
if (lower > min(x)) warning("There are values in the sample smaller than the lower limit")
# Degenerate sample: all observations equal
if (sd(x)==0){
    dx <- density(c(x,x[1]+.Machine$double.eps,x[1]-.Machine$double.eps))
} else {
   # weights
   if(is.null(weights)){
      pesos <- rep(1/length(x),length(x))
   } else {
      pesos <- weights[mantener]
   }
   argumentos <- list(...)
          # edge distance of observations to reflect
   if("bw" %in% names(argumentos)){
       if(is.numeric(argumentos$bw)) broad<-4*argumentos$bw
   } else {
       pesos <- pesos/sum(pesos)
       broad <- 4*density(x,weights=pesos,...)$bw
   }
   # density
   if (is.infinite(lower) & is.infinite(upper)) {
         dx <- density(x, weights=pesos, ...)
   } else if (is.infinite(lower) & is.finite(upper)) {
     reflected<-which(x >= (upper-broad)) # Observations to reflect
     x.reflect<-c(x,2*upper-x[reflected])
     p.reflect<-c(pesos,pesos[reflected])
     p.reflect <- p.reflect/sum(p.reflect)
     dx <- density(x.reflect, weights=p.reflect, ...)
     # Estimation is restricted to the interval
     dx$y<-(dx$y[dx$x>=lower & dx$x<=upper])
     dx$x<-(dx$x[dx$x>=lower & dx$x<=upper])
     # Adjustment of density estimates
     bw<-dx$x[2]-dx$x[1]
     area.under <- sum(dx$y)*bw
     dx$y<-dx$y/area.under
   } else if (is.finite(lower) & is.infinite(upper)) {
     reflected<-which(x <= (lower+broad)) # Observations to reflect
     x.reflect<-c(x,-x[reflected]+2*lower)
     p.reflect<-c(pesos,pesos[reflected])
     p.reflect <- p.reflect/sum(p.reflect)
     dx <- density(x.reflect, weights=p.reflect, ...)
     # Estimation is restricted to the interval
     dx$y<-dx$y[dx$x>=lower & dx$x<=upper]
     dx$x<-dx$x[dx$x>=lower & dx$x<=upper]
     # Adjustment of density estimates
     bw<-dx$x[2]-dx$x[1]
     area.under <- sum(dx$y)*bw
     dx$y<-dx$y/area.under
   } else {
     reflected.inf<-which(x <= (lower+broad)) # Observations to reflect
     reflected.sup<-which(x >= (upper-broad))
     x.reflect<-c(x,-x[reflected.inf]+2*lower)
     p.reflect<-c(pesos,pesos[reflected.inf])
     x.reflect<-c(x.reflect,2*upper-x[reflected.sup])
     p.reflect<-c(p.reflect,pesos[reflected.sup])
     p.reflect <- p.reflect/sum(p.reflect)
     dx <- density(x.reflect, weights=p.reflect, ...)
     # Estimation is restricted to the interval
     dx$y<-dx$y[dx$x>=lower & dx$x<=upper]
     dx$x<-dx$x[dx$x>=lower & dx$x<=upper]
     # Adjustment of density estimates
     bw<-dx$x[2]-dx$x[1]
     area.under <- sum(dx$y)*bw
     dx$y<-dx$y/area.under
   }
}
return(dx)
}
#................................................


dgeometric.test <- function(x, fun.den, par=NULL, lower=-Inf, upper=Inf, n.sim=101, bw=NULL){
# DESCRIPTION:
#             Test to compare a theoretical distribution with an empirical kernel estimate
#             The estatistic measures the distance (in area) between the empirical
#             kernel estimate and a theoretical density function
#
# INPUTS:
#        x: a numeric vector of data values
#        fun.den: an actual density distribution function, such as dnorm. Only continuous
#                densities are valid.
#        par: list of parameters of the density function under the null hypothesis, default NULL
#        lower: lower end point of the support of the variable defined by fun.den, default -Inf
#        upper: upper end point of the support of the variable defined by fun.den, default Inf
#        n.sim: number of iterations performed to calculate the p.value, default 100
#        bw: a numeric value giving the bandwidth to be used in the test, default NULL. If it
#            is not provided, the bandwidth is not constant and is the one provided by the
#            function density under hypothesis of a Gaussian kernel.
#
# OUTPUT:
#       The output is an object of the class htest exactly like for the Kolmogorov-Smirnov
#       test, ks.test. A list containing as most important elements the value of the
#       statistic and its p.value

   data.name <- deparse(substitute(x)) # Name of the variable

   x <- c(na.omit(x))                  # Handling missing values

   if (length(x)<2)
        stop("Not enough observations (at least after removing missing values)")

   if (min(x)<lower || max(x)>upper){
      output <- list(statistic = Inf, p.value = 0, method = "Geometric test", iterations = 0, data.name = data.name)
      warning("There is at least a sample observation out of the theoretical support")
      return(output)
   }

   f<-fun.den
   if (length(par)!=0) f<-function(m) do.call(fun.den,c(list(m),par))

   if (is.null(bw)) bw <- density(x)$bw
   emp.density <- density.reflected(x, lower, upper, bw=bw) # Empirical density
   statistic <- area.between(f, emp.density, lower, upper)

   # p.value estimate by simulation
   n.x<-length(x)
   name.function <- deparse(substitute(fun.den))
   if (nchar(name.function)==1){
       name.function <- NULL
   } else {
       name.function <- paste("r",substr(name.function,2,nchar(name.function)),sep="")
   }
   if (is.null(name.function) || !exists(name.function, mode="function")){
      simulations <- replicate(n.sim, random.function(n=n.x, f=f, lower=lower, upper=upper))
   } else {
      if (is.null(par)){
          simulations <- replicate(n.sim,eval(do.call(name.function,list(n.x))))
      } else {
          simulations <- replicate(n.sim,eval(do.call(name.function,c(n.x,par))))
      }
   }
   density.sim<-apply(simulations, 2, density.reflected, lower=lower, upper=upper, bw=bw)
   statistic.sim <- c(do.call("cbind",lapply(density.sim, area.between, f=f, lower=lower, upper=upper)))
   p.value<-sum(statistic.sim>statistic)/n.sim
   names( statistic ) <- "Tn"
   output <- list(statistic = statistic, p.value = p.value, method = "Geometric test", iterations=n.sim, data.name = data.name)
   class(output) <- "htest"
   return(output)
}
#................................................


fan.test <- function(x, fun.den, par=NULL, lower=-Inf, upper=Inf, kernel="normal", bw=NULL){
# DESCRIPTION:
#             Given a sample of a continuous univariate random variable and a density
#             function fun.den (with suppot in the interval (lower, upper)), fan.test computes
#             the test statistic and the corresponding p-value using a nonparametric
#             kernel approximation to the test with null hypothesis that the sample follows the fun.den
#
# Reference: Fan, Y (1994) "Testing the goodness-of-fit of a parametric density function
#                           by kernel method", Econometric Theory, 10, 316-356.
#
# OBSERVATIONS:
#               -To proper run the function requires the package KernSmooth, to estimate the
#                bandwidth.
#               -The version of the test programed is the one available in Li and Racine (2007, pp.380-1,
#                equation 12.28).
#                Li, O. and Racine, J.F. (2007) Nonparametric Econometrics, Princeton
#                                         University Press, New Jersey.
#               -Tested using R version 3.1.0 and KernSmooth 2.23-8
#               -In its default option, the dpik function available in the library KernSmooth is used.
#                The method for selecting the bandwidth of a kernel density estimate in dpik was proposed by
#                Sheather and Jones (1991) and is described in Section 3.6 of Wand and Jones (1995).
#               -Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
#                              selection method for kernel density estimation.
#                              Journal of the Royal Statistical Society, B, 53.
#               -Wand, M. P. and Jones, M. C. (1995). Kernel Smoothing. Chapman and Hall, London.
#
# INPUTS:
#        x: a numeric vector of data values.
#        fun.den: an actual density distribution function, such as dnorm. Only continuous
#                densities are valid.
#        par: a list of additional parameters of the distribution specified, default NULL.
#      	 lower: lower end point of the support of the variable characterized by fun.den,
#               default -Inf.
#      	 upper: upper end point of the support of the variable characterized by fun.den,
#               default Inf.
#      	 kernel: character string with the kernel to be used, either "normal" (a N(0,1) density),
#               "box" (a uniform in -1 to 1) or "epanech" (a Epanechnikov quadratic kernel), default
#               "normal".
#        bw: a numeric value giving the bandwidth to be used in the test. If it is not provided,
#            the bandwidth is estimated from the data using the dpik function available in the
#            library KernSmooth.
#
# OUTPUT:
#       The output is an object of the class htest exactly like for the Kolmogorov-Smirnov
#       test, ks.test. A list containing as most important elements the value of the
#       statistic and its p.value

  data.name <- deparse(substitute(x))  # Name of the variable

  x <- c(na.omit(x))                   # Handling missing values
  if (length(x)<2)
      stop("Not enough observations (at least after removing missing values)")

  if (!is.numeric(lower) || !is.numeric(upper) || lower >=upper)
      stop("lower < upper is not fulfilled")

  if (min(x)<lower || max(x)>upper){
      output <- list(statistic = -Inf, p.value = 0, method = "Fan's test", kernel=kernel, data.name = data.name)
      warning("There is at least a sample value out of the theoretical support")
      return(output)
  }

  if (is.null(bw)){
    #  require(KernSmooth)
      h.min <- dpik(x,kernel=kernel)       # bandwidth
  } else {
      h.min <- bw
  }

  f<-fun.den
  if (length(par)!=0) f<-function(m) do.call(fun.den,c(list(m),par))

  n.x<-length(x)                       # sample size

  if (kernel=="normal") nucleo <- dnorm
  if (kernel=="box") nucleo <- function(x) ifelse(x>=-1 & x<=1, 0.5, 0)
  if (kernel=="epanech") nucleo <- function(x) ifelse(x>=-1 & x<=1, (3/4)*(1-x^2), 0)

  if (!is.finite(f(lower))) lower <- lower+h.min
  if (!is.finite(f(upper))) upper <- upper-h.min

  # Computing Ig,n
  # First component: Sampling component
  samp1.no.diag <- matrix(nrow=1,ncol=n.x)
  for(j in 1:n.x) samp1.no.diag[1,j] <- sum(nucleo((x-x[j])/h.min))-nucleo((x[j]-x[j])/h.min)
  primera.comp <- sum(samp1.no.diag,na.rm=TRUE)/(h.min*n.x*(n.x-1))

  # Second component: integrate of the square of the kernel convolution of the density null function
  convolution <- function(y,v) (1/h.min)*nucleo((y-v)/h.min)*f(v)
  int.conv <- function(y) sapply(y, function(v) integrate(convolution, lower=lower, upper=upper, v=v, rel.tol=1.e-10)$value)
  fun.cuad <- function(y){int.conv(y)^2}
  segunda.comp <- integrate(fun.cuad,lower=lower,upper=upper,rel.tol=1.e-10,subdivisions=1000)$value

  # Third component: sum of the convolution of the density in the sampled values
  tercera.comp <- (2/n.x)*sum(int.conv(x))

  Ign <- primera.comp + segunda.comp - tercera.comp

  # sigma.cuadrado estimate:
  samp1<-matrix(nrow=1,ncol=n.x)
  for(j in 1:n.x){samp1[1,j] <- sum(nucleo((x-x[j])/h.min)^2)}
  sigma.cuadrado<-2*(1/(h.min*n.x^2))*sum(samp1,na.rm=TRUE)
  
  est<-n.x*sqrt(h.min)*Ign/sqrt(sigma.cuadrado)  # statistic

  # p.value
  p.value<-pnorm(est, lower.tail = FALSE)
  names( est ) <- "Ig"
  # output list
  output <- list(statistic = est, p.value = p.value, method = "Fan's test", kernel=kernel, data.name = data.name)
  class(output) <- "htest"
  return(output)
}
#................................................

# Vector containing the time exposed to risk of death with 76 years during 2006
# for 2006 Spanish immigrants born in 1929
risk76.1929 <- c(0.93975974, 0.8850083, 0.86310107, 0.85215347, 0.79739828, 0.79465792,
                 0.79465448, 0.79192045, 0.78096862, 0.77548848, 0.77001306, 0.76727587,
                 0.76453888, 0.74811108, 0.74811229, 0.73716469, 0.73716057, 0.73716467,
                 0.73442562, 0.73168194, 0.72894432, 0.72620919, 0.72621445, 0.72347073,
                 0.70704561, 0.70704968, 0.70430765, 0.69882752, 0.69609542, 0.69335875,
                 0.68788458, 0.68240598, 0.6769247, 0.67145643, 0.66871539, 0.66324277,
                 0.66324243, 0.65776635, 0.65775908, 0.64954808, 0.6468077, 0.64407303,
                 0.64407156, 0.64134081, 0.64133594, 0.63586141, 0.63312249, 0.63312176,
                 0.63038174, 0.62490714, 0.62217206, 0.62217165, 0.61395989, 0.61395921,
                 0.61122046, 0.61121966, 0.60848053, 0.60848309, 0.60574155, 0.60300606,
                 0.60026798, 0.60026544, 0.59753473, 0.58384606, 0.58384383, 0.58384362,
                 0.58110029, 0.58110727, 0.57562873, 0.57563101, 0.57289158, 0.57288841,
                 0.56467757, 0.56193403, 0.56194171, 0.55920316, 0.55646635, 0.55098816,
                 0.55098505, 0.54551477, 0.5400368, 0.53729881, 0.53182034, 0.52908044,
                 0.52634385, 0.51813598, 0.50444628, 0.50444871, 0.50444261, 0.50444679,
                 0.5017054, 0.50170573, 0.49348819, 0.48528201, 0.48253992, 0.4825382,
                 0.47979858, 0.47706891, 0.47158611, 0.47158945, 0.46885327, 0.46611169,
                 0.46611208, 0.45790069, 0.45789701, 0.44968548, 0.44695285, 0.4469521,
                 0.43873317, 0.43325804, 0.42778055, 0.41957103, 0.41683335, 0.4113554,
                 0.40588583, 0.40588571, 0.40587923, 0.40588465, 0.40314811, 0.40313943,
                 0.40314676, 0.40313879, 0.40313978, 0.40040783, 0.39767138, 0.39492854,
                 0.39493291, 0.39219244, 0.38671393, 0.3867187, 0.38398134, 0.3812414,
                 0.37850661, 0.37576959, 0.37576654, 0.37303017, 0.3730276, 0.37028784,
                 0.37029093, 0.37029104, 0.36755638, 0.367554, 0.3648163, 0.36481084,
                 0.36481439, 0.36207379, 0.35933724, 0.35660073, 0.35386096, 0.35112205,
                 0.34838872, 0.34564714, 0.34291229, 0.34290685, 0.33743961, 0.32922277,
                 0.32374295, 0.32101118, 0.32100545, 0.32100735, 0.32101102, 0.31827017,
                 0.31279774, 0.31279666, 0.31005757, 0.307317, 0.30731882, 0.30184612,
                 0.29636687, 0.29362987, 0.29362652, 0.29089113, 0.2908921, 0.29089085,
                 0.29089252, 0.28267419, 0.28268022, 0.28267379, 0.27446314, 0.27172218,
                 0.26625457, 0.26625436, 0.26077997, 0.26077668, 0.25804207, 0.25803847,
                 0.25529969, 0.25530385, 0.25529909, 0.25256134, 0.25255975, 0.25255828,
                 0.24981942, 0.24982092, 0.24982081, 0.24708856, 0.24160736, 0.24161044,
                 0.24161497, 0.23613479, 0.23613898, 0.23339618, 0.23340153, 0.23066202,
                 0.22244297, 0.21970527, 0.21970559, 0.21970978, 0.21423243, 0.21423196,
                 0.2114982, 0.21149037, 0.21149137, 0.21148972, 0.20875686, 0.20876074,
                 0.20876066, 0.20601671, 0.20601383, 0.20327656, 0.20054677, 0.19506789,
                 0.19506601, 0.19233313, 0.19232527, 0.19232573, 0.1868525, 0.18685787,
                 0.18684835, 0.18684978, 0.18138131, 0.18137531, 0.18137683, 0.17864324,
                 0.17864087, 0.17864223, 0.16768794, 0.16494937, 0.15673341, 0.15400111,
                 0.15399411, 0.15400046, 0.15125699, 0.15126385, 0.14852516, 0.14578728,
                 0.14579006, 0.1430459, 0.14304958, 0.14304549, 0.14031451, 0.14031127,
                 0.14030554, 0.13757326, 0.13483732, 0.13483353, 0.13209225, 0.12935738,
                 0.12662357, 0.12661942, 0.12388228, 0.12388532, 0.12114817, 0.12114219,
                 0.11840873, 0.11841149, 0.11840378, 0.11566763, 0.11567133, 0.11566692,
                 0.10197974, 0.10198221, 0.09650407, 0.09650663, 0.09650818, 0.09650691,
                 0.09377124, 0.09376507, 0.09377005, 0.08828613, 0.08829514, 0.0882919,
                 0.08555037, 0.08554954, 0.08281718, 0.08008109, 0.08007756, 0.0800729,
                 0.07733819, 0.0773428, 0.07734233, 0.07460097, 0.07460242, 0.0746005,
                 0.07186557, 0.07186081, 0.07186548, 0.06912449, 0.06912977, 0.06912403,
                 0.06639034, 0.06364774, 0.06090969, 0.06091666, 0.0609082, 0.05817205,
                 0.055433, 0.05270056, 0.0526969, 0.04995732, 0.04995942, 0.04996077,
                 0.04996539, 0.04995996, 0.04722328, 0.04448264, 0.04448406, 0.04448119,
                 0.04448971, 0.04174433, 0.03627417, 0.03627589, 0.03627122, 0.03627587,
                 0.03627574, 0.03353526, 0.03353788, 0.03079664, 0.03079083, 0.03079667,
                 0.03079845, 0.02805408, 0.02531573, 0.0253152, 0.02258661, 0.02258142,
                 0.02258499, 0.02258552, 0.02258428, 0.01984428, 0.01711121, 0.01710496,
                 0.01710158, 0.0143732, 0.01163319, 0.01163078, 0.01163578, 0.00888996,
                 0.00889585, 0.00888818, 0.00615601, 0.00615705, 0.00615132, 0.00615461,
                 0.00341586, 0.0034219, 0.00341765, 0.00342142, 0.00341658, 0.0006752,
                 0.00067748, 0.00067793)
#................................................