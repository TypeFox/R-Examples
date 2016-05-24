#############################################################################
##                                                                         ##
##   Runuran                                                               ##
##                                                                         ##
##   (c) 2007 Josef Leydold and Wolfgang Hoermann                          ##
##   Department for Statistics and Mathematics, WU Wien                    ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Wrapper for standard distributions                                    ##
##                                                                         ##
#############################################################################
##   Remark: Please sort function calls alphabetically!                    ##
#############################################################################


#############################################################################
## Continuous univariate Distributions                                      #
#############################################################################

## -- Beta distribution - (replacement for rbeta) ---------------------------
urbeta <- function (n,shape1,shape2,lb=0,ub=1) {
        unr<-new("unuran", paste("beta(",shape1,",",shape2,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)	
}

udbeta <- function (shape1,shape2,lb=0,ub=1) {
  if (missing (shape1) || missing (shape2))
    stop ("argument 'shape1' or 'shape2' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "beta", c(shape1,shape2), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Burr distribution -----------------------------------------------------
urburr <- function (n,a,b,lb=0,ub=Inf) {
## works in theory for a >= 1 and b >= 2 
## numerical problems for a*b > 175 or so
        unr <- new("unuran", paste("distr=cont;pdf='",a*(b-1),"*x^(",a-1,")/(1+x^",a,")^",b,"'; domain=(",lb,",",ub,")"),"TDR")
        unuran.sample(unr,n)
}

## TODO

## -- Cauchy distribution - (replacement for rcauchy) -----------------------
urcauchy <- function (n,location=0,scale=1,lb=-Inf,ub=Inf) {
  unr<-new("unuran",paste("cauchy(",location,",",scale,"); domain=(",lb,",",ub,")"),"HINV")
  unuran.sample(unr,n)	
}

udcauchy <- function (location=0,scale=1,lb=-Inf,ub=Inf) {
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "cauchy", c(location,scale), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Chi distribution ------------------------------------------------------
urchi <- function (n,df,lb=0,ub=Inf) {
        unr <- new("unuran", paste("chi(",df,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udchi <- function (df,lb=0,ub=Inf) {
  if (missing (df))
    stop ("argument 'df' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "chi", c(df), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Chi^2 distribution - (replacement for rchisq) -------------------------
urchisq <- function (n,df,lb=0,ub=Inf) {
        unr <- new("unuran", paste("chisquare(",df,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udchisq <- function (df,lb=0,ub=Inf) {
  if (missing (df))
    stop ("argument 'df' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "chisquare", c(df), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Exponential distribution - (replacement for rexp) ---------------------
urexp <- function (n,rate=1,lb=0,ub=Inf) {
        unr <- new("unuran", paste("exponential(",1./rate,"); domain=(",lb,",",ub,")"), "CSTD")
        unuran.sample(unr,n)
}

udexp <- function (rate=1,lb=0,ub=Inf) {
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "exponential", c(1./rate), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- F distribution  - (replacement for rf) --------------------------------
urf <- function (n,df1,df2,lb=0,ub=Inf) {
        unr <- new("unuran", paste("F(",df1,",",df2,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udf <- function (df1,df2,lb=0,ub=Inf) {
  if (missing (df1) || missing (df2) )
    stop ("argument 'df1' or 'df2' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "F", c(df1,df2), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Frechet (Extreme Value type II) distribution --------------------------
urextremeII <- function (n,shape,location=0,scale=1,lb=location,ub=Inf) {
        unr <- new("unuran", paste("extremeII(",shape,",",location,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udfrechet <- function (shape,location=0,scale=1,lb=location,ub=Inf) {
  if (missing (shape))
    stop ("argument 'shape' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "extremeII", c(shape,location,scale), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Gamma distribution  - (replacement for rgamma) ------------------------
urgamma <- function (n,shape,scale=1,lb=0,ub=Inf) {
        unr<-new("unuran", paste("gamma(",shape,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udgamma <- function (shape,scale=1,lb=0,ub=Inf) {
  if (missing (shape))
    stop ("argument 'shape' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "gamma", c(shape,scale), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Generalized hyperbolic distribution -----------------------------------
udghyp <- function (lambda,alpha,beta,delta,mu, lb=-Inf,ub=Inf) {
  if (missing (lambda) || missing (alpha) || missing (beta)
      || missing (delta) || missing (mu))
    stop ("argument 'lambda', 'alpha', 'beta', 'delta', or 'mu' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "ghyp", c(lambda,alpha,beta,delta,mu),
                      c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Generalized inverse Gaussian ------------------------------------------
urgig <- function (n,lambda,omega,lb=1.e-12,ub=Inf) { 
        ## works for lambda>=1 and omega>0 and for lambda>0 and omega>=0.5
        unr<-new("unuran",
                 paste("cont; pdf='x^(",lambda-1,")*exp(-(",omega/2,")*(x+1/x))'; domain=(",lb,",",ub,")"),"TDR")
        unuran.sample(unr,n)
}

udgig <- function (theta,psi,chi, lb=0,ub=Inf) {
  if (missing (theta) || missing (psi) || missing (chi) )
    stop ("argument 'theta', 'psi' or 'chi' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "gig2", c(theta,psi,chi), c(lb,ub), PACKAGE="Runuran")
  distr
}

udgiga <- function (theta,omega,eta=1, lb=0,ub=Inf) {
  if (missing (theta) || missing (omega) )
    stop ("argument 'theta' or 'omega' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "gig", c(theta,omega,eta), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Gumbel (Extreme Value type I) distribution ----------------------------
urextremeI <- function (n,location=0,scale=1,lb=-Inf,ub=Inf) {
        unr <- new("unuran", paste("extremeI(",location,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udgumbel <- function (location=0,scale=1,lb=-Inf,ub=Inf) {
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "extremeI", c(location,scale), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Hyperbolic distribution -----------------------------------------------
urhyperbolic <- function (n,shape,scale=1,lb=-Inf,ub=Inf) {
        unr <- new("unuran",
                   paste("cont; pdf='exp(-",shape,"*sqrt(1.+x*x/",(scale*scale),"))'; domain=(",lb,",",ub,")"),
                   "TDR")
	 unuran.sample(unr,n)
}

udhyperbolic <- function (alpha,beta,delta,mu, lb=-Inf,ub=Inf) {
  if (missing (alpha) || missing (beta) || missing (delta) || missing (mu))
    stop ("argument 'alpha', 'beta', 'delta',or 'mu' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "hyperbolic", c(alpha,beta,delta,mu), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Inverse Gaussian ------------------------------------------------------
udig <- function (mu,lambda, lb=0,ub=Inf) {
  if (missing (mu) || missing (lambda) )
    stop ("argument 'mu' or 'lambda' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "ig", c(mu,lambda), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Laplace (double exponential) distribution -----------------------------
urlaplace <- function (n,location=0,scale=1,lb=-Inf,ub=Inf) {
        unr <- new("unuran", paste("laplace(",location,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udlaplace <- function (location=0,scale=1,lb=-Inf,ub=Inf) {
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "laplace", c(location,scale), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Lognormal distribution  - (replacement for rlnorm) --------------------
urlnorm <- function (n,meanlog=0,sdlog=1,lb=0,ub=Inf) {
        exp(urnorm(n,meanlog,sdlog,log(lb),log(ub)))
}

udlnorm <- function (meanlog=0,sdlog=1, lb=0,ub=Inf) {
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "lognormal", c(meanlog,sdlog), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Logistic distribution - (replacement for rlogistic) -------------------
urlogis <- function (n,location=0,scale=1,lb=-Inf,ub=Inf) {
        unr <- new("unuran", paste("logistic(",location,",",scale,"); domain=(",lb,",",ub,")"), "CSTD")
        unuran.sample(unr,n)
}

udlogis <- function (location=0,scale=1,lb=-Inf,ub=Inf) {
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "logistic", c(location,scale), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Lomax distribution (Pareto distribution of second kind) ---------------
urlomax <- function (n,shape,scale=1,lb=0,ub=Inf) {
        unr <- new("unuran", paste("lomax(",shape,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udlomax <- function (shape,scale=1,lb=0,ub=Inf) {
  if (missing (shape))
    stop ("argument 'shape' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "lomax", c(shape,scale), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Meixner distribution --------------------------------------------------
udmeixner <- function (alpha,beta,delta,mu, lb=-Inf,ub=Inf) {
  if (missing (alpha) || missing (beta) || missing (delta) || missing (mu))
    stop ("argument 'alpha', 'beta', 'delta', or 'mu' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "meixner", c(alpha,beta,delta,mu),
                      c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Normal (Gaussian) distribution - (replacement for rnorm) --------------
urnorm <- function (n,mean=0,sd=1,lb=-Inf,ub=Inf) {
        unr<-new("unuran",paste("normal(",mean,",",sd,"); domain=(",lb,",",ub,")"),"HINV")
        unuran.sample(unr,n)
}

udnorm <- function (mean=0,sd=1,lb=-Inf,ub=Inf) {
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "normal", c(mean,sd), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Pareto distribution ---------------------------------------------------
urpareto <- function (n,k,a,lb=k,ub=Inf) {
        unr <- new("unuran", paste("pareto(",k,",",a,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udpareto <- function (k,a,lb=k,ub=Inf) {
  if (missing (k) || missing (a))
    stop ("argument 'k' or 'a' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "pareto", c(k,a), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Planck distribution ---------------------------------------------------
urplanck <- function (n,a,lb=1.e-12,ub=Inf) { 
        ## works for a>=1 
        unr <- new("unuran", paste("cont; pdf='x^",a,"/(exp(x)-1)'; domain=(",lb,",",ub,")"), "TDR")
        unuran.sample(unr,n)
}

#udplanck <- function (a,lb=1.e-12,ub=Inf) { 
#  distr <- new ("unuran.cont",empty=TRUE)
#  distr@distr <-.Call("Runuran_std_cont", distr, "planck", c(a), c(lb,ub), PACKAGE="Runuran")
#  distr
#}
## TODO

## -- Powerexponential (Subbotin) distribution ------------------------------
urpowerexp <- function (n,shape,lb=-Inf,ub=Inf) {
        unr <- new("unuran", paste("powerexponential(",shape,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udpowerexp <- function (shape,lb=-Inf,ub=Inf) {
  if (missing (shape))
    stop ("argument 'shape' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "powerexponential", c(shape), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Rayleigh distribution -------------------------------------------------
urrayleigh <- function (n,scale=1,lb=0,ub=Inf) {
        unr <- new("unuran", paste("rayleigh(",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udrayleigh <- function (scale=1,lb=0,ub=Inf) {
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "rayleigh", c(scale), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Slash distribution ----------------------------------------------------
udslash <- function (lb=-Inf,ub=Inf) {
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "slash", numeric(0), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Student's t distribution - (replacement for rt) -----------------------
urt <- function (n,df,lb=-Inf,ub=Inf) { 
        unr <- new("unuran", paste("student(",df,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udt <- function (df,lb=-Inf,ub=Inf) { 
  if (missing (df))
    stop ("argument 'df' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "student", c(df), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Triangular distribution with lower border a, mode m and upper border b 
urtriang <- function (n,a,m,b,lb=a,ub=b) {
        if (a>=b || m<=a || m>=b)
                stop("Invalid arguments for a,m,b")
        l <- c(a*a, -2*a, 1) / ((a-b)*(a-m))
        r <- c(a*b-a*m+b*m, -2*b, 1) / ((a-b)*(b-m))
        cdfstring <- paste("'(x<=",m,")*(",l[1],"+(",l[2],")*x+(",l[3],")*x*x)+(x>",
                           m,")*(",r[1],"+(",r[2],")*x+(",r[3],")*x*x)';", sep="")
        domainstring <- paste("domain=(",max(lb,a),",",min(ub,b),")", sep="")
        unr <- new("unuran", paste("cont; cdf=",cdfstring,domainstring), "HINV")
        unuran.sample(unr,n)
}

#udtriang <- function (df,lb=-Inf,ub=Inf) { 
#  distr <- new ("unuran.cont",empty=TRUE)
#  distr@distr <-.Call("Runuran_std_cont", distr, "student", c(df), c(lb,ub), PACKAGE="Runuran")
#  distr
#}
## TODO

## -- Variance gamma distribution -------------------------------------------
udvg <- function (lambda, alpha, beta, mu, lb=-Inf, ub=Inf) {
  if (missing (lambda) || missing (alpha) || missing (beta) || missing (mu))
    stop ("argument 'lambda', 'alpha', 'beta', or 'mu' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "vg", c(lambda,alpha,beta,mu),
                      c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Weibull distribution - (replacement for rweibull) ---------------------
urweibull <- function (n,shape,scale=1,lb=0,ub=Inf) {
        unr <- new("unuran", paste("weibull(",shape,",",scale,"); domain=(",lb,",",ub,")"), "HINV")
        unuran.sample(unr,n)
}

udweibull <- function (shape,scale=1,lb=0,ub=Inf) {
  if (missing (shape))
    stop ("argument 'shape' missing")
  distr <- new ("unuran.cont",empty=TRUE)
  distr@distr <-.Call("Runuran_std_cont", distr, "weibull", c(shape,scale), c(lb,ub), PACKAGE="Runuran")
  distr
}


#############################################################################
## Discrete univariate Distributions                                        #
#############################################################################


## -- Binomial distribution - (replacement for rbinom) ----------------------
urbinom <- function (n,size,prob,lb=0,ub=size) { 
        unr <- new("unuran", paste("binomial(",size,",",prob,"); domain=(",lb,",",ub,")"), "DGT")
        unuran.sample(unr,n)
}

udbinom <- function (size,prob,lb=0,ub=size) {
  if (missing (size) || missing (prob))
    stop ("argument 'size' or 'prob' missing")
  distr <- new ("unuran.discr",empty=TRUE)
  distr@distr <-.Call("Runuran_std_discr", distr, "binomial", c(size,prob), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Geometric distribution - (replacement for rgeom) ----------------------
urgeom <- function (n,prob,lb=0,ub=Inf) {
        if (prob > 0.02) {
                ub  <- min(ub,2000);
                unr <- new("unuran", paste("geometric(",prob,"); domain=(",lb,",",ub,")"), "DGT");
	}
        else {
                unr <- new("unuran", paste("geometric(",prob,"); domain=(",lb,",",ub,")"), "DARI");
        }
        unuran.sample(unr,n)
}
 
udgeom <- function (prob,lb=0,ub=Inf) {
  if (missing (prob))
    stop ("argument 'prob' missing")
  distr <- new ("unuran.discr",empty=TRUE)
  distr@distr <-.Call("Runuran_std_discr", distr, "geometric", c(prob), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Hypergeometric distribution - (replacement for rhyper) ----------------
urhyper <- function (nn,m,n,k,lb=max(0,k-n),ub=min(k,m)) {
        unr <- new("unuran", paste("hypergeometric(",m+n,",",m,",",k,"); domain=(",lb,",",ub,")"), "DGT")
        unuran.sample(unr,nn)
}

udhyper <- function (m,n,k,lb=max(0,k-n),ub=min(k,m)) {
  if ( missing (m) || missing (n) || missing(k))
    stop ("argument 'm', 'n' or 'k' missing")
  distr <- new ("unuran.discr",empty=TRUE)
  distr@distr <-.Call("Runuran_std_discr", distr, "hypergeometric", c(m+n,m,k), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Logarithmic distribution ----------------------------------------------
urlogarithmic <- function (n,shape,lb=1,ub=Inf) {
        if(shape<0.98) {
                ub  <- min(ub,2000);
                unr <- new("unuran", paste("logarithmic(",shape,"); domain=(",lb,",",ub,")"), "DGT")
        }
        else {
                unr <- new("unuran", paste("logarithmic(",shape,"); domain=(",lb,",",ub,")"), "DARI")
        }
        unuran.sample(unr,n)
}

udlogarithmic <- function (shape,lb=1,ub=Inf) {
  if (missing (shape))
    stop ("argument 'shape' missing")
  distr <- new ("unuran.discr",empty=TRUE)
  distr@distr <-.Call("Runuran_std_discr", distr, "logarithmic", c(shape), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Negative binomial distribution - (replacement for rnbinom) ------------
urnbinom <- function (n,size,prob,lb=0,ub=Inf) {
        if (pnbinom(1000,size,prob,lower.tail=F) < 1.e-10){
                ub  <- min(ub,1000);
                unr <- new("unuran", paste("negativebinomial(",prob,",",size,"); domain=(",lb,",",ub,")"), "DGT")
        }
        else {
                unr <- new("unuran", paste("negativebinomial(",prob,",",size,"); domain=(",lb,",",ub,")"), "DARI")
        }
        unuran.sample(unr,n)
}

udnbinom <- function (size,prob,lb=0,ub=Inf) {
  if (missing (size) || missing (prob))
    stop ("argument 'size' or 'prob' missing")
  distr <- new ("unuran.discr",empty=TRUE)
  distr@distr <-.Call("Runuran_std_discr", distr, "negativebinomial", c(prob,size), c(lb,ub), PACKAGE="Runuran")
  distr
}

## -- Poisson distribution - (replacement for rpois) ------------------------
urpois <- function (n,lambda,lb=0,ub=Inf) {
        if (ppois(1000,lambda,lower.tail=F) < 1.e-10) {
                ub <- min(ub,1000);
                unr <- new("unuran", paste("poisson(",lambda,"); domain=(",lb,",",ub,")"), "DGT")
        }
        else {
                unr <- new("unuran", paste("poisson(",lambda,"); domain=(",lb,",",ub,")"), "DARI")
        }
        unuran.sample(unr,n)
}

udpois <- function (lambda,lb=0,ub=Inf) {
  if (missing (lambda))
    stop ("argument 'lambda' missing")
  distr <- new ("unuran.discr",empty=TRUE)
  distr@distr <-.Call("Runuran_std_discr", distr, "poisson", c(lambda), c(lb,ub), PACKAGE="Runuran")
  distr
}


#############################################################################
## Continuous multivariate Distributions:                                   #
##    Yet net implemented.                                                  #
#############################################################################
