.coda_diags.stats = list(
    studentbridge=c("E","var","autocov","loglik_mean","loglik_extremum","extremum","ratio_extremum","ratio_loglik_extremum"),
    brownianbridge=c("E","var","autocov","loglik_mean","loglik_extremum","extremum","ratio_extremum","ratio_loglik_extremum"),
    loglikbridge=c("E","var","autocov","extremum","ratio_extremum"))
.coda_diags.rhos = seq(from=0.0,to=.9,by=.1)
.coda_diags.Ns = c(10,30,100,300,1000,3000)#,10000)

.onLoad <- function(libname, pkgname) {
    load(system.file("coda_diags.ecdf.Rdata",package="codadiags"), envir=parent.env(environment()))
    #data("coda_diags.ecdf", package="codadiags", envir=parent.env(environment()))
}
#load(system.file("coda_diags.ecdf.Rdata",package="codadiags"))

#' Build the null CDF (cumulative density function) for a given statistic, for arbitrary length and autocorrelation sequence.
#' @param stat statistic used
#' @param N length of the target sequence to be tested
#' @param rho autocorrelation (1st coeff) of the target sequence to test
null.param.cdf <- function(stat, N,rho) {
    .coda_diags.rhos.ix1 = which.min(abs(.coda_diags.rhos-rho))
    rho1 = .coda_diags.rhos[.coda_diags.rhos.ix1]
    if (rho1 >= max(.coda_diags.rhos) || rho1 <= min(.coda_diags.rhos) || rho == rho1) {
        rho2 = rho1
    } else {
        if (rho>rho1) {
            .coda_diags.rhos.ix2 = .coda_diags.rhos.ix1 + 1
        } else {
            .coda_diags.rhos.ix2 = .coda_diags.rhos.ix1 - 1
        }
        rho2 = .coda_diags.rhos[.coda_diags.rhos.ix2]
    }
    
    .coda_diags.Ns.ix1 = which.min(abs(.coda_diags.Ns-N))
    N1 = .coda_diags.Ns[.coda_diags.Ns.ix1]
    if (N1 >= max(.coda_diags.Ns) || N1 <= min(.coda_diags.Ns) || N1 == N) {
        N2 = N1
    } else {
        if (N>N1) {
            .coda_diags.Ns.ix2 = .coda_diags.Ns.ix1 + 1
        } else {
            .coda_diags.Ns.ix2 = .coda_diags.Ns.ix1 - 1
        }
        N2 = .coda_diags.Ns[.coda_diags.Ns.ix2]
    }    
    
    ecdf.N1.rho1 = get(paste(".coda_diags.",(stat),"_N",N1,"_rho",rho1,sep=""))
    if (N1 != N2) {
        ecdf.N2.rho1 = get(paste(".coda_diags.",(stat),"_N",N2,"_rho",rho1,sep=""))
        if (rho1 != rho2) {
            ecdf.N1.rho2 = get(paste(".coda_diags.",(stat),"_N",N1,"_rho",rho2,sep=""))
            ecdf.N2.rho2 = get(paste(".coda_diags.",(stat),"_N",N2,"_rho",rho2,sep=""))
            
            null.ecdf = function(v){
                ((N-N2)*((rho-rho2)*ecdf.N1.rho1(v) + (rho-rho1)*ecdf.N1.rho2(v)) / (2*rho-rho1-rho2) +(N-N1)*((rho-rho2)*ecdf.N2.rho1(v) + (rho-rho1)*ecdf.N2.rho2(v)) / (2*rho-rho1-rho2)) / (2*N-N1-N)
            }
        } else {
            null.ecdf = function(v){
                ((N-N2)*ecdf.N1.rho1(v) + (N-N1)*ecdf.N2.rho1(v)) / (2*N-N1-N2)
            }
        }
    } else {
        if (rho1 != rho2) {
            ecdf.N1.rho2 = get(paste(".coda_diags.",(stat),"_N",N1,"_rho",rho2,sep=""))
            
            null.ecdf = function(v){
                ((rho-rho2)*ecdf.N1.rho1(v) + (rho-rho1)*ecdf.N1.rho2(v)) / (2*rho-rho1-rho2)
            }
        } else {
            null.ecdf = ecdf.N1.rho1
        }
    }
    
    return(null.ecdf)
    
    #rho = .coda_diags.rhos[which.min(abs(.coda_diags.rhos-rho))]
    #N = .coda_diags.Ns[which.min(abs(.coda_diags.Ns-N))]
    #return(get(paste((stat),"_N",N,"_rho",rho,sep="")))
}

#' Asymptotic CDF for a given statistic
#' @param stat statistic
#' @param forceUseECDF if true, 
#' - if stat is loglik_mean.brownianbridge, use the Anderson-Darling CDF
#' #NO - if stat is var.brownianbridge, use the Cramer von Mises CDF
#' - if stat is extremum.brownianbridge, use the Kolmogorov-Smirnov CDF
#' #NO - if stat is ratio_loglik_extremum.brownianbridge, use the chi-square (3 freedom degrees) CDF
#' - if stat is ratio_extremum.brownianbridge, use the Bay CDF
#' else, use tabulated empirical CDF built on white noise process (length 10000)
null.lim.cdf <- function(stat, forceUseECDF=FALSE) {
    if (!isTRUE(forceUseECDF)) {
        if (stat=="loglik_mean.brownianbridge") {return(function(x)1-ad.cdf(-x))}
        # if (stat=="var.brownianbridge") {return(cvm.cdf)} # no: \int_0^1 B^2 != var(B)
        if (stat=="extremum.brownianbridge") {return(ks.cdf)}
        # if (stat=="ratio_loglik_extremum.brownianbridge") {return(maxinv.chi23.cdf)}# no: chi23 is correspondig to loglik_maximum.bronwiabridge stat. , not ratio_...
        if (stat=="ratio_extremum.brownianbridge") {return(maxinv.bay.cdf)}
    }
    return(get(paste(".coda_diags.",(stat),"_asymptotic",sep="")))    
}

#' Anderson-Darling cumulative density function, copy from ADGofTest package.
#' @author Carlos J. Gil Bellosta
#' @param x value to test
#' @param n sample size for Anderson-Darling statistic
#' @references G. and J. Marsaglia, "Evaluating the Anderson-Darling Distribution", Journal of Statistical Software, 2004
#' @examples
#' require(codadiags)
#' plot(null.lim.cdf("loglik_mean.brownianbridge",forceUseECDF=TRUE),col='blue',xlim=c(-4,0))
#' plot(Vectorize(function(x)1-ad.cdf(-x)),col='green',add=TRUE,xlim=c(-4,0))
ad.cdf <- function( x, n=1000 ) {
    if( x < 2 )
        x <- exp(-1.2337141/x)/sqrt(x)*(2.00012+(.247105- (.0649821-(.0347962-(.011672-.00168691*x)*x)*x)*x)*x)
    else
        x <- exp(-exp(1.0776-(2.30695-(.43424-(.082433-(.008056 -.0003146*x)*x)*x)*x)*x))
    
    
    if( x > 0.8 )
        return( x + (-130.2137+(745.2337-(1705.091-(1950.646-(1116.360-255.7844*x)*x)*x)*x)*x)/n )
    
    z <- 0.01265 + 0.1757 / n
    
    if( x < z ){ 
        v <- x / z
        v <- sqrt(v)*(1.-v)*(49*v-102)
        return ( x + v * (.0037/(n*n)+.00078/n+.00006)/n )
    }
    
    v <- (x-z) / (0.8-z)
    v <- -0.00022633+(6.54034-(14.6538-(14.458-(8.259-1.91864*v)*v)*v)*v)*v
    x + v * (.04213+.01365/n)/n
}

#' Cramer von Mises cumulative density function, import from coda package.
#' @param x value to test
#' @seealso coda::pcramer
#' @references Csorgo S. and Faraway, JJ. The exact and asymptotic distributions of the Cramer-von Mises statistic. J. Roy. Stat. Soc. (B), 58, 221-234 (1996).
#' @examples
#' require(codadiags)
#' plot(null.lim.cdf("var.brownianbridge",forceUseECDF=TRUE),col='blue')
#' plot(Vectorize(cvm.cdf),col='green',add=TRUE)
cvm.cdf <- function ( x ) {
   coda::pcramer(x)
}

#' Kolmogorov-Smirnov cumulative density function, copy from stats::ks.test.
#' @param x value to test
#' @param n sample size for Kolmogorov-Smirnov statistic
#' @useDynLib codadiags
#' @seealso package stats, ks.c
#' @references George Marsaglia, Wai Wan Tsang and Jingbo Wang (2003), Evaluating Kolmogorov's distribution. Journal of Statistical Software, 8/18. http://www.jstatsoft.org/v08/i18/.
#' @examples 
#' require(codadiags)
#' plot(null.lim.cdf("extremum.brownianbridge",forceUseECDF=TRUE),col='blue',xlim=c(0.01,4))
#' plot(Vectorize(ks.cdf),col='green',add=TRUE,xlim=c(0.01,4))
ks.cdf <- function( x, n=100 ) {
    .Call("pKolmogorov2x", x/sqrt(n), n)
}

#' Bay cumulative density function, corresponding to -B(t+)/B(t-), where B(t+) (resp. B(t-)) is the maximum (resp.minimum) of B(t)/(t*(1-t)).
#' @param x value to test
#' @author X. Bay
bay.cdf <- function( x ) {
    if(x==0) return(0)
    y=pi/(1+x)
    return(1-x*y^2/pi*(1/y-cos(y)/sin(y)))
}

#' CDF of max(x,1/x) (=cdf(x)-cdf(1)+cdf(1)-cdf(1/x)) where x is 'Bay' distributed
#' @param x value to test
#' @examples
#' require(codadiags)
#' plot(null.lim.cdf("ratio_extremum.brownianbridge",forceUseECDF=TRUE),col='blue',xlim=c(0,10))
#' plot(Vectorize(maxinv.bay.cdf),col='green',add=TRUE,xlim=c(0,10))
maxinv.bay.cdf <- function( x ) {
    if (x<1) return(0)
    return(bay.cdf(x)-bay.cdf(1/x))
}
