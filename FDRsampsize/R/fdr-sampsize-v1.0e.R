#> BEGIN DESCRIPTION
#> Package: FDRsampsize
#> Type: Package
#> Title: Compute Sample Size that Meets Requirements for Average Power and FDR
#> Version: 1.0
#> Author: Stan Pounds <stanley.pounds@stjude.org>
#> Maintainer: Stan Pounds <stanley.pounds@stjude.org>
#> Depends: R (>= 2.15.1)
#> Imports: stats
#> Date: 2016-01-06
#> Description: Defines a collection of functions to compute average power and sample size for studies that use the false discovery rate as the final measure of statistical significance.
#> License: GPL-2
#> END DESCRIPTION


#############################################################
#> BEGIN NAMESPACE
#> export(afdr,alpha.fdr,alpha.power,average.power,fdr.power,fdr.sampsize,pi.star,power.cox,power.hart,power.li,power.onesampt,power.oneway,power.ranksum,power.signtest,power.tcorr,power.twosampt)
#> importFrom("stats","pf","pnorm","qf","qnorm")
#> S3method(print,fdr.sampsize)
#> END NAMESPACE

#############################################################
#> BEGIN FDRsampsize-package
#> \alias{FDRsampsize-package}
#> \alias{FDRsampsize}
#> \docType{package}
#> \title{An R package to Perform Power and Sample Size Calculations for Microarray Studies}
#> \description{ A general approach to performing power and sample size calculations for microarray
#> studies has been developed in the literature.  However, the software associated with those articles
#> implements the approach only for studies that will perform the t-test or one-way ANOVA to compare 
#> gene expression across two or more groups.  Here, we describe a set of R routines that implement 
#> the general method for power and sample size calculations for a wider variety of statistical tests.  
#> These routines accept the name of a function that computes the power for the statistical test of 
#> interest and thus have the flexibility to perform calculations for virtually any statistical test 
#> with a known power formula.}
#> \details{
#> \tabular{ll}{
#> Package: \tab FDRsampsize\cr
#> Type: \tab Package\cr
#> Version: \tab 1.0\cr
#> Date: \tab 2016-01-06\cr
#> License: \tab GPL(>=2) \cr }
#> }
#> \author{Stan Pounds <stanley.pounds@stjude.org> }
#> \references{A Onar-Thomas, S Pounds. FDRsampsize: An R package to Perform Generalized Power and Sample Size Calculations for Planning Studies that use the False Discovery Rate to Measure Significance, Manuscript 2016.
#>
#> Pounds, Stan, and Cheng Cheng. "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271.
#>
#>  Jung, Sin-Ho. "Sample size for FDR-control in microarray data analysis." Bioinformatics 21.14 (2005): 3097-3104.}
#> END FDRsampsize-package


######################################
#> BEGIN fdr.sampsize
#> \alias{print.fdr.sampsize}
#> \title{Determine sample size required to achieve a desired average power while controlling the FDR at a specified level.}
#> \description{Determines the sample size needed to achieve a desired average power while controlling the FDR at a specified level.}
                                        #> \arguments{
fdr.sampsize=function(fdr,             #> \item{fdr}{Desired FDR (scalar numeric)}
                       ave.pow,         #> \item{ave.pow}{Desired average power (scalar numeric)}
                       pow.func,        #> \item{pow.func}{Character string name of function to compute power; must accept n, alpha, and eff.size as its first three arguments.  Other optional arguments may also be provided.}
                       eff.size,        #> \item{eff.size}{Numeric vector of effect sizes; interally, this will be provided as the third argument of pow.func}
                       null.effect,     #> \item{null.effect}{Scalar value of the effect size under the null hypothesis.  This may be 0 or 1 for tests that respectively use differences or ratios for comparisons.}
                       max.n=500,       #> \item{max.n}{Maximum n to consider}
                       min.n=5,         #> \item{min.n}{Minimum n to consider}
                       tol=0.00001,     #> \item{tol}{Tolerance for bisection calculations}
                       eps=0.00001,     #> \item{eps}{Epsilon for numerical differentiation}
                       lam=0.95,        #> \item{lam}{Lambda for computing anticipated pi0 estimate, see Storey 2002.}
                       ...)             #> \item{...}{additional arguments for pow.func} }

{
	# Compute alpha needed to give desired power at maximum sample size
  maxpow=alpha.power(ave.pow,
	                    max.n,
	                    pow.func,
	                    eff.size,
	                    null.effect,
	                    tol,...)
	
  # Determine the anticipated FDR with that sample size and alpha
  minfdr=afdr(max.n,
	             maxpow$alpha,
	             pow.func,
	             eff.size,
	             lam,eps,...)
  
	pi0=mean(eff.size==null.effect) # Get the true value of pi0  

	
	if (minfdr>fdr) # Indicates that max.n cannot achieve specified FDR level
	{
		pi.hat=pi.star(max.n,pow.func,eff.size,lam=lam,eps=eps,...)
		act.fdr=pi0*maxpow$alpha/(pi0*maxpow$alpha+(1-pi0)*maxpow$ave.pow)
		res=list(OK=F,
		          n=max.n,
		          alpha=maxpow$alpha,
		          fdr.hat=minfdr,
		          act.fdr=act.fdr,
		          ave.pow=maxpow$ave.pow,
		          act.pi=pi0,
		          pi.hat=pi.hat,
		          eff.size=eff.size)
		class(res)="fdr.sampsize"
		return(res)	 # Return what is achievable with maximum sample size	          
	}
	
	# Determine FDR and power for minimum sample size
	minpow=alpha.power(ave.pow,min.n,pow.func,eff.size,null.effect,tol,...)
	maxfdr=afdr(min.n,minpow$alpha,pow.func,eff.size,lam,eps,...)

	if (maxfdr<fdr)
	{
		pi.hat=pi.star(min.n,pow.func,eff.size,lam=lam,eps=eps,...)
		act.fdr=pi0*minpow$alpha/(pi0*minpow$alpha+(1-pi0)*minpow$ave.pow)
		res=list(OK=T,
		          n=min.n,
		          alpha=minpow$alpha,
		          fdr.hat=maxfdr,
		          act.fdr=act.fdr,
		          ave.pow=minpow$ave.pow,
		          act.pi=pi0,
		          pi.hat=pi.hat,
		          eff.size=eff.size)
		class(res)="fdr.sampsize"
		return(res)			
	}
	res0=n.fdr(ave.pow,fdr,pow.func,eff.size,null.effect,
	            lam=lam,eps=eps,n0=min.n,n1=max.n,tol=tol,...)
	res=list(OK=T,
	          n=res0$n,
	          alpha=res0$alpha,
	          fdr.hat=res0$fdr.hat,
	          act.fdr=res0$act.fdr,
	          ave.pow=res0$ave.pow,
	          act.pi=res0$act.pi,
	          pi.hat=res0$pi.hat,
	          eff.size=eff.size)
	class(res)="fdr.sampsize"
	return(res)
}
#> \details{This function checks the technical conditions regarding whether the desired FDR can be
#>          achieved by min.n or max.n before calling n.fdr.  Thus, for most applications,
#>          fdr.sampsize should be used instead of n.fdr.}
#> \value{\code{fdr.sampsiz}e returns an object of class 'FDRsampsize', which is a list with the following components:
#>            \item{OK}{indicates if the requirement is met}
#>            \item{n}{the computed sample size}
#>            \item{alpha}{the p-value threshold that gives the desired FDR}
#>            \item{fdr.hat}{anticipated value of the FDR estimate given n and effect size }
#>            \item{act.fdr}{actual expected FDR given n and effect size}
#>            \item{ave.pow}{average power}
#>            \item{act.pi}{actual value of pi0, the proportion of tests with a true null hypothesis.}
#>            \item{pi.hat}{expected value of the pi0 estimate}
#>            \item{eff.size}{input effect size vector}
#>           }

#> \references{A Onar-Thomas, S Pounds. "FDRsampsize: An R package to Perform Generalized Power and Sample Size Calculations for Planning Studies that use the False Discovery Rate to Measure Significance", Manuscript 2015.
#> 
#> Pounds, Stan, and Cheng Cheng. "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271.
#> 
#>             Jung, Sin-Ho. "Sample size for FDR-control in microarray data analysis." Bioinformatics 21.14 (2005): 3097-3104.}

#> \examples{
#  # Load the package
#> power.twosampt             # show the power.cox function
#> res=fdr.sampsize(fdr=0.1,
#>                  ave.pow=0.8,
#>                  pow.func=power.twosampt,
#>                  eff.size=rep(c(1,0),c(10,990)),
#>                  null.effect=0)
#> res
#> }
#> END fdr.sampsize

# BEGIN print.FDRsampsize
# \title{Print the result object returned by fdr.sampsize.}
# \description{Print an object of class 'FDRsampsize', which is returned by fdr.sampsize.}
                                           # \arguments{
print.fdr.sampsize=function(x,...) # \item{x}{an FDRsampsize object, which is the result of \link{fdr.sampsize}}
                                           # \item{...}{additional arguments for print}
                                           # }
{
  object=x
  res=c(n=round(object$n,2),
        alpha=object$alpha,
        ave.pow=object$ave.pow,
        fdr.hat=object$fdr.hat,
        act.fdr=object$act.fdr)
  print(res,...)
}
# \value{prints a vector with the sample size estimate, 
#         the p-value threshold alpha, the average power,
#         the anticipated FDR estimate, and the 
#         anticipated actual FDR}

# END print.FDRsampsize

######################################
#> BEGIN fdr.power
#> \title{Compute the average power at a specific FDR control level}
#> \description{Compute the average power at a specific level of FDR control for a given effect size and sample size}
                                       #> \arguments{
fdr.power=function(fdr,                #> \item{fdr}{Desired FDR, scalar}
                      n,               #> \item{n}{sample size}
               pow.func,               #> \item{pow.func}{name of R function to compute statistical power}
               eff.size,               #> \item{eff.size}{effect size vector; will be provided as the third argument of pow.func}
            null.effect,               #> \item{null.effect}{value of effect size that corresponds to null hypothesis}
               lam=0.95,               #> \item{lam}{name of R function to compute statistical power}
             eps=0.0001,               #> \item{eps}{epsilon for numerical differentiation}
             tol=0.0001,               #> \item{tol}{tolerance for bisection solution to alpha}
                    ...)               #> \item{...}{additional agruments for the functions} }

{
	a0=alpha.fdr(fdr,n,pow.func,eff.size,null.effect,lam,eps,tol,...)
	if (!a0$OK) return(0)
	res=average.power(n,a0$alpha,pow.func,eff.size,null.effect=null.effect,...)
	return(res)
}
#> \value{average power (scalar) of the tests with a false null hypothesis}
#> \examples{
#> fdr.power(fdr=0.10,n=50,pow.func=power.twosampt,
#>           eff.size=rep(0:1,c(900,100)),null.effect=0)
#>}
#> \references{A Onar-Thomas, S Pounds "FDRsampsize: An R package to Perform Generalized Power and Sample Size Calculations for Planning Studies that use the False Discovery Rate to Measure Significance", Manuscript 2016.
#>
#> Gadbury GL, et al. (2004) Power and sample size estimation in high dimensional biology.  Statistical Methods in Medical Research 13(4):325-38.
#>
#>             Pounds S and Cheng C (2005) Sample size determination for the false discovery rate.  Bioinformatics 21(23): 4263-71.}
#> END fdr.power


######################################
#> BEGIN alpha.fdr
#> \title{Find the fixed p-value threshold that controls the FDR at a specified level}
#> \description{Find the p-value threshold that satisfies an FDR requirement (if such a threshold exists)}
                                       #> \arguments{
alpha.fdr=function(fdr,               #> \item{fdr}{Desired FDR, scalar}
                      n,               #> \item{n}{sample size}
               pow.func,               #> \item{pow.func}{an R function to compute statistical power}
               eff.size,               #> \item{eff.size}{effect size vector}
            null.effect,               #> \item{null.effect}{value of effect size that corresponds to the null hypothesis}
               lam=0.95,               #> \item{lam}{the lambda parameter in computing the pi0 (proportion of tests with a true null) estimate of Storey (2002)}
             eps=0.0001,               #> \item{eps}{epsilon for numerical differentiation}
             tol=0.0001,               #> \item{tol}{tolerance for bisection solution to alpha}
                    ...)               #> \item{...}{additional agruments for the functions} }

{
	pi.hat=pi.star(n,pow.func,eff.size,lam,eps,...)
	
	pow.temp=pow.func(n,eps,eff.size,...)
	pr.rej=mean(pow.temp)
	pdf0=exp(log(pr.rej)-log(eps))
	minpossfdr=pi.hat/pdf0
	if (minpossfdr>fdr) return(list(alpha=0,fdr=minpossfdr,OK=F))
	
	a0=0
	a1=1
	fdr0=minpossfdr
	fdr1=mean(eff.size==null.effect)
	
	k=ceiling(-log(tol)/log(2))
	for (i in 1:k)
	{
		a.mid=(a0+a1)/2
		fdr.mid=afdr(n,a.mid,pow.func,eff.size,lam=lam,eps=eps,...)
		if (fdr.mid>fdr)
		{
			a1=a.mid
			fdr1=fdr.mid
		}
		else
		{
			a0=a.mid
			fdr0=fdr.mid
		}
	}
	return(list(alpha=a0,fdr=fdr0,OK=T))
}
#> \value{a list with the following components:
#>            \item{fdr}{the FDR at that alpha}
#>            \item{alpha}{the determined alpha}
#>            \item{OK}{indicates if the requirement is met}
#>           }
#> \examples{
#> alpha.fdr(fdr=0.1,n=50,pow.func=power.twosampt,
#>           eff.size=rep(0:1,c(900,100)),null.effect=0)
#> }
#> \references{A Onar-Thomas, S Pounds "FDRsampsize: An R package to Perform Generalized Power and Sample Size Calculations for Planning Studies that use the False Discovery Rate to Measure Significance", Manuscript 2015.
#>
#>             Pounds, Stan, and Cheng Cheng. "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271.
#> 
#>             Jung, Sin-Ho. "Sample size for FDR-control in microarray data analysis." Bioinformatics 21.14 (2005): 3097-3104.}
#> END alpha.fdr

######################################
#> BEGIN n.fdr
#> \title{Find the sample size that meets desired FDR and power criteria}
#> \description{Find smallest sample size that meets requirements for average power and FDR}
                                 #> \arguments{
n.fdr=function(ave.pow,         #> \item{ave.pow}{required average power (scalar)}
                fdr,             #> \item{fdr}{required FDR (scalar)}
                pow.func,        #> \item{pow.func}{name of R function that computes statistical power}
                eff.size,        #> \item{eff.size}{effect size vector}
                null.effect,     #> \item{null.effect}{Value of effect size that indicates null}
                lam=0.95,        #> \item{lam}{p-value at which to evaluate ensemble PDF}
                eps=0.0001,      #> \item{eps}{epsilon for numerical differentiation}
                n0=5,            #> \item{n0}{smallest sample size to be considered for bisection}
                n1=500,          #> \item{n1}{maximum sample size to be considered for bisection}
                tol=0.000001,    #> \item{tol}{tolerance for solving for alpha in iterations}
                ...)             #> \item{...}{additional agruments for the functions} }
{
	a0=alpha.power(ave.pow,n0,pow.func,eff.size,null.effect,tol,...)
	a1=alpha.power(ave.pow,n1,pow.func,eff.size,null.effect,tol,...)
	afdr0=afdr(n0,a0$alpha,pow.func,eff.size,lam=lam,eps=eps,...)
	afdr1=afdr(n1,a0$alpha,pow.func,eff.size,lam=lam,eps=eps,...)
	#print(afdr0)
	k=ceiling((log(n1-n0)-log(tol))/log(2))
	i=0
	cl.eq=(ceiling(n1)==ceiling(n0))
	while((i<k)&(!cl.eq))
	{
		i=i+1
		n.mid=(n1+n0)/2
		a.mid=alpha.power(ave.pow,n.mid,pow.func,eff.size,null.effect,tol,...)
		afdr.mid=afdr(n.mid,a.mid$alpha,pow.func,eff.size,lam=lam,eps=eps,...)
		if (afdr.mid>fdr)
		{
			n0=n.mid
			afdr0=afdr.mid
		}
		else
		{
			n1=n.mid
			afdr1=afdr.mid
		}
		cl.eq=(ceiling(n1)==ceiling(n0))
	}
	n=ceiling(n1)
	a=alpha.power(ave.pow,n,pow.func,eff.size,null.effect,tol,...)
	pow.res=average.power(n,a$alpha,pow.func,eff.size,null.effect=null.effect,...)
	fdr.hat=afdr(n,a$alpha,pow.func,eff.size,lam=lam,eps=eps,...)
	pi.hat=pi.star(n,pow.func,eff.size,lam=lam,eps=eps,...)
	act.pi0=mean(eff.size==null.effect)
	act.fdr=act.pi0*a$alpha/(act.pi0*a$alpha+(1-act.pi0)*pow.res)
	res=list(n=n,
	            alpha=a$alpha,
	            fdr.hat=fdr.hat,
	            ave.pow=pow.res,
	            act.fdr=act.fdr,
	            pi.hat=pi.hat,
	            act.pi=act.pi0)
	class(res)="FDRsampsize"
	return(res)
}
#> \details{This performs the sample size calculation without checking whether the
#>          minimum or maximum sample size satisfy the desired requirements.  The fdr.sampsize
#>          function checks these conditions and then calls n.fdr.  Thus, many users will
#>          may prefer to use the sampsize.fdr function instead of n.fdr.}
#> \value{a list with the following components:
#>            \item{n}{a sample size estimate}
#>            \item{alpha}{the p-value cut-off}
#>            \item{fdr.hat}{an approximation to the expected value of the FDR estimate given n }
#>            \item{ave.pow}{the average power}
#>            \item{fdr.act}{the actual FDR given n}
#>            \item{pi.hat}{expected value of the pi.hat estimator given n}
#>            \item{act.pi}{actual pi0}
#>           }
#> \references{A Onar-Thomas, S Pounds. "FDRsampsize: An R package to Perform Generalized Power and Sample Size Calculations for Planning Studies that use the False Discovery Rate to Measure Significance", Manuscript 2015.
#> 
#>             Pounds, Stan, and Cheng Cheng. "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271.
#> 
#>             Jung, Sin-Ho. "Sample size for FDR-control in microarray data analysis." Bioinformatics 21.14 (2005): 3097-3104.}
#> END n.fdr



######################################
#> BEGIN average.power
#> \title{Compute average power for a given sample size}
#> \description{Compute average power for given sample size, effect size, and p-value threshold}
                                          #> \arguments{
average.power=function(n,                #> \item{n}{sample size}
                     alpha,               #> \item{alpha}{p-value cut off (scalar)}
                     pow.func,            #> \item{pow.func}{an R function to compute statistical power}
                     eff.size,            #> \item{eff.size}{effect size vector}
                     null.effect,         #> \item{null.effect}{value of effect size that corresponds to null hypothesis}
                     ...)                 #> \item{...}{additional agruments for the functions} }
{
  pi=mean(eff.size==null.effect)
	pr.rej=mean(pow.func(n,alpha,eff.size,...))
	res=exp(log(pr.rej-pi*alpha)-log(1-pi))
	return(res)
}
#> \value{average power (scalar)}
#> \examples{
#> average.power(n=50,alpha=0.01,pow.func=power.twosampt,
#>               eff.size=rep(0:1,c(900,100)),null.effect=0)
#> }
#> \references{
#> Pounds, Stan, and Cheng Cheng. "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271.
#> Gadbury GL, et al. (2004) Power and sample size estimation in high dimensional biology.  Statistical Methods in Medical Research 13(4):325-38.
#> Jung, Sin-Ho. "Sample size for FDR-control in microarray data analysis." Bioinformatics 21.14 (2005): 3097-3104.}
#> END average.power


######################################
#> BEGIN afdr
#> \title{Compute the anticipated FDR}
#> \description{Compute the anticipated FDR for given sample size, p-value threshold, and effect size.}
                              #> \arguments{
afdr=function(n,             #> \item{n}{sample size (scalar)}
              alpha,         #> \item{alpha}{p-value cut-off (scalar)}
              pow.func,      #> \item{pow.func}{an R function that computes statistical power}
              eff.size,      #> \item{eff.size}{effect size vector}
              lam=0.95,      #> \item{lam}{p-value at which to evaluate ensemble PDF (for pi.star)}
              eps=0.0001,    #> \item{eps}{epsilon for numerical differentiation}
              ...)          #> \item{...}{additional agruments for the functions} }
{
  pi.hat=pi.star(n,pow.func,eff.size,lam,eps,...)
	pow.temp=pow.func(n,alpha,eff.size,...)
	pr.rej=mean(pow.temp)
	return(pi.hat*alpha/pr.rej)
}

#> \value{the aFDR}
#> \details{The aFDR is defined by Pounds and Cheng (2005) as the anticipated false discovery rate
#>          incurred by performing all tests with p-value threshold alpha given the same size
#>          effect size and power function.}
#> \examples{
#> afdr(n=50,alpha=0.01,pow.func=power.twosampt,eff.size=rep(c(1,0),c(100,900)))
#> }
#> \references{Pounds, Stan, and Cheng Cheng. "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271.
#> 
#>             Jung, Sin-Ho. "Sample size for FDR-control in microarray data analysis." Bioinformatics 21.14 (2005): 3097-3104.}
#> END afdr



######################################
#> BEGIN alpha.power
#> \title{Find the p-value threshold that gives a specified average power}
#> \description{Find p-value cut-off that yields desired average power given n and effect size}
                                             #> \arguments{
alpha.power=function(ave.pow,               #> \item{ave.pow}{desired average power (scalar)}
                            n,               #> \item{n}{sample size}
                     pow.func,               #> \item{pow.func}{name of R function to compute statistical power}
                     eff.size,               #> \item{eff.size}{effect size vector}
                  null.effect,               #> \item{null.effect}{value of effect size that corresponds to null hypothesis}
                 tol=0.000001,               #> \item{tol}{tolerance for bisection solution to alpha}
                          ...)               #> \item{...}{additional agruments for the functions} }
{
  pow0=0
	pow1=1
	a0=0
	a1=1
	k=ceiling(-log(tol)/log(2))
	for (i in 1:k)   # Bisection loop
	{
		a.mid=(a0+a1)/2
		pow.mid=average.power(n,a.mid,pow.func,eff.size,null.effect=null.effect,...)
		#print(pow.mid)
		if (pow.mid>ave.pow)
		{
			a1=a.mid
			pow1=pow.mid
		}
		else
		{
			a0=a.mid
			pow0=pow.mid
		}
	}
	return(list(alpha=a1,ave.pow=pow1))
}
#> \value{a list with the following components:
#>            \item{alpha}{desired value of alpha}
#>            \item{ave.pow}{average power at that alpha}
#>           }
#> \examples{
#> alpha.power(ave.pow=0.8,n=50,pow.func=power.twosampt,
#>             eff.size=rep(0:1,c(900,100)),null.effect=0)
#> }
#> \references{A Onar-Thomas, S Pounds. "FDRsampsize: An R package to Perform Generalized Power and Sample Size Calculations for Planning Studies that use the False Discovery Rate to Measure Significance", Manuscript 2015.
#> Pounds, Stan, and Cheng Cheng. "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271.
#> 
#>             Jung, Sin-Ho. "Sample size for FDR-control in microarray data analysis." Bioinformatics 21.14 (2005): 3097-3104.}
#> END alpha.power


######################################
#> BEGIN pi.star
#> \title{Compute the anticipated null proportion estimate}
#> \description{Compute an approximation of the expected value of the null proportion estimate given the sample size and effect size.}
                                     #> \arguments{
pi.star=function(n,                 #> \item{n}{sample size}
           pow.func,                 #> \item{pow.func}{an R function to compute statistical power}
           eff.size,                 #> \item{eff.size}{effect size vector}
           lam=0.95,                 #> \item{lam}{p-value at which to numerically evaluate p-value pdf (scalar)}
         eps=0.0001,                 #> \item{eps}{epsilon for numerical differentiation}
                ...)                 #> \item{...}{additional agruments for the functions} }
{
  pow.lam=pow.func(n,lam,eff.size,...)
	pow.lam.eps=pow.func(n,lam-eps,eff.size,...)
  res=exp(log(mean(pow.lam)-
	             mean(pow.lam.eps))-
	         log(eps))
	return(res)
}
#> \value{scalar value for approximated E(pi.hat)}
#> \references{#> Pounds, Stan, and Cheng Cheng. "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271.}
#> END pi.star



######################################
#> BEGIN power.cox
#> \title{Compute the power of a single-predictor Cox regression model}
#> \description{Use the formula of Hseih and Lavori (2000) to compute the power of a single-predictor Cox model.}
                          #> \arguments{
power.cox=function(n,    #> \item{n}{number of events (scalar)}
                alpha,    #> \item{alpha}{p-value threshold (scalar)}
                logHR,    #> \item{logHR}{log hazard ratio (vector)}
                    v)    #> \item{v}{variance of predictor variable (vector)}
                          #> }
{
   pnorm(qnorm(alpha/2),sqrt(n*v)*logHR)+
   1-pnorm(qnorm(1-alpha/2),sqrt(n*v)*logHR)
}
#> \value{vector of power estimates for two-sided test}
#> \references{
#> Hsieh, FY and Lavori, Philip W (2000) Sample-size calculations for the Cox proportional hazards regression model with nonbinary covariates.  Controlled Clinical Trials 21(6):552-560.
#> }
#> \examples{
#  # Load the package
#> power.cox             # show the power.cox function
#> res=fdr.sampsize(fdr=0.1,
#>                  ave.pow=0.8,
#>                  pow.func=power.cox,
#>                  eff.size=rep(c(log(2),0),c(100,900)),
#>                  null.effect=0,
#>                  v=1)
#> res



 
#>}
#> END power.cox

######################################
#> BEGIN power.onesampt
#> \title{Compute power of the one-sample t-test}
#> \description{Estimate power of the one-sample t-test;Uses classical power formula for one-sample t-test}
                               #> \arguments{
power.onesampt=function(n,    #> \item{n}{number of events (scalar)}
                     alpha,    #> \item{alpha}{p-value threshold (scalar)}
                     delta,    #> \item{delta}{difference of actual mean from null mean (vector)}
                   sigma=1)    #> \item{sigma}{standard deviation (vector or scalar, default=1)}
                               #> }

{
   1-pf(qf(1-alpha,1,n-1),1,n-1,ncp=n*(delta/sigma)^2)
}
#> \value{vector of power estimates for two-sided test}
#> \examples{
#  # Load the package
#> power.onesampt        # show the power function
#> res=fdr.sampsize(fdr=0.1,
#>                  ave.pow=0.8,
#>                  pow.func=power.onesampt,
#>                  eff.size=rep(c(2,0),c(100,900)),
#>                  null.effect=0,
#>                  sigma=1)
#> res



 
#>}
#> END power.onesampt


######################################
#> BEGIN power.twosampt
#> \title{Compute power of the two-samples t-test}
#> \description{Estimate power of the two-samples t-test;Uses classical power formula for two-sample t-test;Assumes equal variance and sample size }
                               #> \arguments{
power.twosampt=function(n,    #> \item{n}{per-group sample size (scalar)}
                     alpha,    #> \item{alpha}{p-value threshold (scalar)}
                     delta,    #> \item{delta}{difference between population means (vector)}
                   sigma=1)    #> \item{sigma}{standard deviation (vector or scalar)}
                               #> }
{
  1-pf(qf(1-alpha,1,2*(n-1)),1,2*(n-1),ncp=n*(delta/sigma)^2/2)
}
#> \value{vector of power estimates for two-sided test}
#> \details{For many applications, the null.effect is zero difference of means.}
#> \examples{
#  # Load the package
#> power.twosampt        # show the power function
#> res=fdr.sampsize(fdr=0.1,
#>                  ave.pow=0.8,
#>                  pow.func=power.twosampt,
#>                  eff.size=rep(c(2,0),c(100,900)),
#>                  null.effect=0,
#>                  sigma=1)
#> res



 
#>}
#> END power.twosampt



######################################
#> BEGIN power.oneway
#> \title{Compute power of one-way ANOVA}
#> \description{Compute power of one-way ANOVA;Uses classical power formula for ANOVA;Assumes equal variance and sample size }
                             #> \arguments{
power.oneway=function(n,    #> \item{n}{per-group sample size (scalar)}
                   alpha,    #> \item{alpha}{p-value threshold (scalar)}
                   theta,    #> \item{theta}{sum of ((group mean - overall mean)/stdev)^2 across all groups for each hypothesis test (vector)}
                   k=2)      #> \item{k}{the number of groups to be compared, default k=2}
                             #> }
{
   1-pf(qf(1-alpha,k-1,k*(n-1)),k-1,k*(n-1),ncp=n*theta)
}
#> \value{vector of power estimates for test of equal means}
#> \details{For many applications, the null effect is zero for the parameter theta described above.}
#> \examples{
#  # Load the package
#> power.oneway        # show the power function
#> res=fdr.sampsize(fdr=0.1,
#>                  ave.pow=0.8,
#>                  pow.func=power.oneway,
#>                  eff.size=rep(c(2,0),c(100,900)),
#>                  null.effect=0,
#>                  k=3)
#> res



 
#>}
#> END power.oneway



######################################
#> BEGIN power.tcorr
#> \title{Compute Power of the t-test for non-zero correlation}
#> \description{Estimate power of t-test for non-zero correlation;Uses classical power formula for t-test}
                             #> \arguments{
power.tcorr=function(n,     #> \item{n}{sample size (scalar)}
                   alpha,    #> \item{alpha}{p-value threshold (scalar)}
                    rho)     #> \item{rho}{population correlation coefficient (vector)}
                             #> }
{
   1-pf(qf(1-alpha,1,n-2),1,n-2,ncp=(n-2)*rho^2/(1-rho^2))
}
#> \value{vector of power estimates for two-sided tests}
#> \details{For many applications, the null.effect is rho=0.}
#> \examples{
#  # Load the package
#> power.tcorr        # show the power function
#> res=fdr.sampsize(fdr=0.1,
#>                  ave.pow=0.8,
#>                  pow.func=power.tcorr,
#>                  eff.size=rep(c(0.3,0),c(100,900)),
#>                  null.effect=0)
#> res



 
#>}
#> END power.tcorr


######################################
#> BEGIN power.ranksum
#> \title{Compute power of the rank-sum test}
#> \description{Compute power of rank-sum test;Uses formula of Noether (JASA 1987)}
                               #> \arguments{
power.ranksum=function(n,     #> \item{n}{sample size (scalar)}
                   alpha,      #> \item{alpha}{p-value threshold (scalar)}
                    p)         #> \item{p}{Pr (Y>X), as in Noether (JASA 1987)}
                               #> }
{  mu0=0.5*n*n
   mu1=p*n*n
   delta=(mu1-mu0)
   sig2=n*n*(2*n+1)/12
   sig=sqrt(sig2)
   z.reject=-abs(qnorm(alpha/2,0,sig))
   pow=pnorm(z.reject,delta,sig)+pnorm(-z.reject,delta,sig,lower.tail=F)
   return(pow)
}

#> \value{vector of power estimates for two-sided tests}
#> \details{In most applications, the null effect size will be designated by p = 0.5, which indicates that
#           the probability that a randomly selected subject from group A will have a
#           greater value of the variable than a randomly selected subject from group B.  
#>          Thus, in the example below, the argument null.effect=0.5 is specified in the call to fdr.sampsize.}

#> \examples{
#  # Load the package
#> power.ranksum        # show the power function
#> res=fdr.sampsize(fdr=0.1,
#>                  ave.pow=0.8,
#>                  pow.func=power.ranksum,
#>                  eff.size=rep(c(0.8,0.5),c(100,900)),
#>                  null.effect=0.5)
#> res



 
#>}
#> \references{Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests.  
#>             Journal of the American Statistical Association, 82:645-647.}
#> END power.ranksum

######################################
#> BEGIN power.signtest
#> \title{Compute power of the sign test}
#> \description{Use the Noether (1987) formula to compute the power of the sign test}
                               #> \arguments{
power.signtest=function(n,    #> \item{n}{sample size (scalar)}
                   alpha,      #> \item{alpha}{p-value threshold (scalar)}
                    p)         #> \item{p}{Pr (Y>X), as in Noether (JASA 1987)}
                               #> }
{
    mu0=0.5*n
    mu1=p*n
    sig0=sqrt(n*0.25)
    sig1=sqrt(n*p*(1-p))
    pow=pnorm(qnorm(alpha/2,mu0,sig0),mu1,sig1)+
      pnorm(qnorm(1-alpha/2,mu0,sig0),mu1,sig1,lower.tail=F)
    return(pow)
}
#> \value{vector of power estimates for two-sided tests}
#> \details{In most applications, the null effect size will be designated by p = 0.5 instead of p = 0.
#>          Thus, in the call to fdr.sampsize, we specify null.effect=0.5 in the example below.}

#> \references{Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests.  
#>             Journal of the American Statistical Association, 82:645-647.}
#> \examples{
#  # Load the package
#> power.signtest        # show the power function
#> res=fdr.sampsize(fdr=0.1,
#>                  ave.pow=0.8,
#>                  pow.func=power.signtest,
#>                  eff.size=rep(c(0.8,0.5),c(100,900)),
#>                  null.effect=0.5)
#> res



 
#>}
#> END power.signtest

#> BEGIN power.hart
#> \title{Compute Power for RNA-seq Experiments Assuming Negative Binomial Distribution.}
#> \description{Use the formula of Hart et al (2013) to compute power for comparing RNA-seq expression across two groups assuming a negative binomial distribution.}
                             #> \arguments{
power.hart=function(n,       #> \item{n}{per-group sample size (scalar)}
                    alpha,   #> \item{alpha}{p-value threshold (scalar)}
                    log.fc,  #> \item{log.fc}{log fold-change (vector), usual null hypothesis is log.fc=0}
                    mu,      #> \item{mu}{read depth per gene (vector, same length as log.fc)}
                    sig)     #> \item{sig}{coefficient of variation (CV) per gene (vector, same length as log.fc)}
                             #> }
{
  z.alpha=qnorm(alpha/2)
  res=pnorm(z.alpha,sqrt(n*log.fc^2/(2*(1/mu+sig^2))))+
      pnorm(-z.alpha,sqrt(n*log.fc^2/(2*(1/mu+sig^2))),lower.tail=F)
  return(res)
}
#> \value{vector of power estimates for the set of two-sided tests}
#> \references{SN Hart, TM Therneau, Y Zhang, GA Poland, and J-P Kocher (2013).  
#>             Calculating Sample Size Estimates for RNA Sequencing Data.  
#>             Journal of Computational Biology 20: 970-978.}
#> \details{This function is based on equation (1) of Hart et al (2013).  It
#>          assumes a negative binomial model for RNA-seq read counts and 
#>          equal sample size per group.}
#> \examples{
#  # Load the package
#> power.hart       # show the power function
#> n.hart=2*(qnorm(0.975)+qnorm(0.9))^2*(1/20+0.6^2)/(log(2)^2) # Equation 6 of Hart et al
#> power.hart(n.hart,0.05,log(2),20,0.6)                        # Recapitulate 90% power  
#> res=fdr.sampsize(fdr=0.1,
#>                  ave.pow=0.8,
#>                  pow.func=power.hart,
#>                  eff.size=rep(c(log(2),0),c(100,900)),
#>                  null.effect=0,mu=5,sig=1)
#> res



 
#>}

#> END power.hart

#> BEGIN power.li
#> \title{Compute Power for RNA-Seq Experiments Assuming Poisson Distribution}
#> \description{Use the formula of Li et al (2013) to compute power for comparing RNA-seq expression across two groups assuming the Poisson distribution.}
#>                \arguments{
power.li=function(n,         #> \item{n}{per-group sample size}
                  alpha,     #> \item{alpha}{p-value threshold}  
                  rho,       #> \item{rho}{fold-change, usual null hypothesis is that rho=1}
                  mu0,       #> \item{mu0}{average count in control group}
                  w=1,       #> \item{w}{ratio of total number of }  
                  type="w")  #> \item{type}{type of test: "w" for Wald, "s" for score, "lw" for log-transformed Wald, "ls" for log-transformed score.}
#>                           }
{
  z.alpha=qnorm(alpha/2)
  if (!is.element(type,c("w","lw","s","ls")))
    stop("type must be 'w', 's', 'ls', or 'lw'.")
  if (any(w<=0)) stop("w must be >0.")
  if (any(rho<=0)) stop("rho must be >0.")
  if (type=="w")  
  {
    mu.z=sqrt(n*mu0*(rho-1)^2/(1+rho/w))
    res=pnorm(z.alpha,mu.z)+pnorm(-z.alpha,mu.z,lower.tail=F)
    return(res)
  }
  if (type=="lw")
  {
    mu.z=sqrt(n*mu0*(log(rho)^2)/(1+1/(rho*w)))
    res=pnorm(z.alpha,mu.z)+pnorm(-z.alpha,mu.z,lower.tail=F)
    return(res)
  }
  if (type=="s")
  {
    z.cut=z.alpha*sqrt((1+w*rho)/(w+rho))
    mu.z=sqrt(n*mu0*(rho-1)^2/(1+rho/w))
    res=pnorm(z.cut,mu.z)+pnorm(-z.cut,mu.z,lower.tail=F)
    return(res)
  }
  if (type=="ls")
  {
    z.cut=z.alpha*sqrt((2+w+1/w)/((1+w*rho)*(1+1/(w*rho))))
    mu.z=sqrt(n*mu0*(log(rho)^2)/(1+1/(rho*w)))
    res=pnorm(z.cut,mu.z)+pnorm(-z.cut,mu.z,lower.tail=F)
    return(res)
  }

}
#> \value{vector of power estimates for two-sided tests}
#> \references{C-I Li, P-F Su, Y Guo, and Y Shyr (2013).  
#>             Sample size calculation for differential expression analysis of
#>             RNA-seq data under Poisson distribution.
#>             Int J Comput Biol Drug Des 6(4).  doi:10.1504/IJCBDD.2013.056830}
#> \details{This function computes the  power for each of a series of two-sided
#>          tests defined by the input parameters.  The power is based on the
#>          sample size formulas in equations 10-13 of Li et al (2013).
#>          Also, note that the null.effect is set to 1 in the examples because
#>          the usual null hypothesis is that the fold-change = 1.}
#> \examples{
#  # Load the package
#> power.li      # show the power function
#> power.li(88,0.05,1.25,5,0.5,"w")  # recapitulate 80% power in Table 1 of Li et al (2013)
#> res=fdr.sampsize(fdr=0.1,
#>                  ave.pow=0.8,
#>                  pow.func=power.li,
#>                  eff.size=rep(c(1.5,1),c(100,900)),
#>                  null.effect=1,
#>                  mu0=5,w=1,type="w")
#> res



 
#>}

#> END power.li

