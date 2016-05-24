### fdrtool.R  (2013-09-15)
###
###    Estimate (Local) False Discovery Rates For Diverse Test Statistics
###    
###
### Copyright 2007-13 Korbinian Strimmer
###
### This file is part of the `fdrtool' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA



 
fdrtool = function(x, 
  statistic=c("normal", "correlation", "pvalue"),
  #statistic=c("normal", "correlation", "pvalue", "studentt"),
  plot=TRUE, color.figure=TRUE, verbose=TRUE, 
  cutoff.method=c("fndr", "pct0", "locfdr"),
  pct0=0.75)
{
  statistic = match.arg(statistic)
  cutoff.method = match.arg(cutoff.method)
  
  if ( is.vector(x) == FALSE )
  	stop("input test statistics must be given as a vector!")
 
  if ( length(x) < 200 ) warning("There may be too few input test statistics for reliable FDR calculations!")
  if (statistic=="pvalue")
  {
    if (max(x) > 1 | min(x) < 0) 
      stop("input p-values must all be in the range 0 to 1!")
  }

  
#### step 1 ####

  if(verbose) cat("Step 1... determine cutoff point\n")

  # determine cutoff point for censoring

  if (cutoff.method=="pct0")
  {
    # use specified quantile

    if(statistic=="pvalue") x0 = quantile(x, probs=1-pct0)
    else x0 = quantile(abs(x), probs=pct0)
  }
  else if ( cutoff.method=="locfdr" & (statistic=="normal" | statistic=="correlation") )
  {
    # use same procedure as in locfdr R package (due to Brit Katzen-Turnbull)

    if(statistic=="normal") z = x
    if(statistic=="correlation") z = atanh(x)

    iqr = as.double(diff(quantile(z, probs=c(.25, .75))))
    sdhat = iqr/(2*qnorm(.75))
    N = length(z)
    # b = 3.55-.44*log(N, 10)                               # locfdr 1.1-3
    b = ifelse(N > 500000, 1,  4.3 * exp(-0.26*log(N,10)) ) # locfdr 1.1-6
    z0 = b*sdhat
    
    if(statistic=="normal") x0 = z0
    if(statistic=="correlation") x0 = tanh(z0)
  }
  else
  {
    if(cutoff.method=="locfdr")
    warning("cutoff.method=\"locfdr\" only available for normal and correlation statistic.")

    # control false nondiscovery rate

    x0 = fndr.cutoff(x, statistic)
  }



#### step 2 ####

  if(verbose) cat("Step 2... estimate parameters of null distribution and eta0\n")

  cf.out <- censored.fit(x=x, cutoff=x0, statistic=statistic)
# cf.out looks as follows for p-values
#     cutoff N0      eta0    eta0.var
#[1,]   0.96 64 0.3730473 0.002141996
# the other test statistics have two more colums containing the scale parameter

  if (statistic=="pvalue")
    scale.param = NULL
  else
    scale.param <- cf.out[1,5] # variance parameter

  eta0 = cf.out[1,3]

#### step 2 ####

  if(verbose) cat("Step 3... compute p-values and estimate empirical PDF/CDF\n")

  nm = get.nullmodel(statistic)
  pval = nm$get.pval(x, scale.param)

  # determine cumulative empirical distribution function (pvalues)
  ee <- ecdf.pval(pval, eta0=eta0)

  g.pval <- grenander(ee)

  #cat("DEBUG: Grenander eta0=", g.pval$f.knots[length(g.pval$f.knots)], "\n")
  #cat("DEBUG: estimated eta0=", eta0 , "\n\n")

  # mixture density and CDF  
  f.pval = approxfun( g.pval$x.knots,  g.pval$f.knots, method="constant", rule=2)
  f0.pval = function(x) return( ifelse(x > 1 | x < 0, 0, rep(1, length(x))) )  

  F.pval = approxfun( g.pval$x.knots,  g.pval$F.knots, method="linear", 
           yleft=0, yright=g.pval$F.knots[length(g.pval$F.knots)])
  F0.pval = function(x) return( ifelse(x > 1, 1, ifelse(x < 0, 0, x )) )

  #fdr.pval = function(p) pmin( eta0   / f.pval(p), 1) # eta0*f0/ f
  fdr.pval = function(p)
  {
    p[ p == .Machine$double.eps ] = 0
    pmin( eta0   / f.pval(p), 1) # eta0*f0/ f
  }

  Fdr.pval = function(p) pmin( eta0*p / F.pval(p), 1) # eta0*F0/ F
  

#### step 4 ####

  if(verbose) cat("Step 4... compute q-values and local fdr\n")

  qval <- Fdr.pval(pval) 
  lfdr <- fdr.pval(pval)


 
#### return results ####


  result = list(pval=pval, qval=qval, lfdr=lfdr, 
             statistic=statistic, param=cf.out)

  if (plot)
  {
    if(verbose) cat("Step 5... prepare for plotting\n")

    ##############
    # zeta > 0 in the following
     

    
 
    if(statistic=="pvalue")
    {
      f0 <- function(zeta) return( nm$f0(zeta, scale.param) ) 
      F0 <- function(zeta) return( nm$F0(zeta, scale.param) )
      get.pval <- function(zeta) return( nm$get.pval(1-zeta, scale.param) )
      x0 = 1-x0
    } 
    else
    {
      f0 <- function(zeta) return( 2*nm$f0(zeta, scale.param)  ) 
      F0 <- function(zeta) return( 2*nm$F0(zeta, scale.param)-1  )
      get.pval <- function(zeta) return( nm$get.pval(zeta, scale.param) )
    }

    fdr = function(zeta)  fdr.pval(get.pval(zeta)) 
    Fdr = function(zeta)  Fdr.pval(get.pval(zeta)) 
     
    F   = function(zeta)  1-eta0*get.pval(zeta)/Fdr(zeta) 
    FA  = function(zeta)  (F(zeta)-eta0*F0(zeta))/(1-eta0)		

    f   = function(zeta)  eta0*(f0(zeta))/fdr(zeta) 
    fA  = function(zeta)  (f(zeta)-eta0*f0(zeta))/(1-eta0)		

   
    ##############

    ax = abs(x) 
    if (statistic=="pvalue") ax = 1-ax  # reverse p-val plot 

    xxx = seq(0, max(ax), length.out=500)
    
    ll = pvt.plotlabels(statistic, scale.param, eta0)

    par(mfrow=c(3,1))

    if (color.figure)
      cols = c(2,4) # colors for f0,F0 and fA,FA
    else
      cols = c(1,1)

    hist(ax, freq=FALSE, bre=50,
      main=ll$main, xlab=ll$xlab, cex.main=1.8)
    lines(xxx, eta0*f0(xxx), col=cols[1], lwd=2, lty=3 )
    lines(xxx, (1-eta0)*fA(xxx), col=cols[2], lwd=2 )
    if (statistic=="pvalue") 
      pos1 = "topleft" else pos1="topright"
    legend(pos1, 
      c("Mixture", "Null Component", "Alternative Component"), 
      lwd=c(1, 2, 2), col=c(1,cols), lty=c(1,3,1), bty="n", cex=1.5)
 
    plot(xxx, F(xxx), lwd=1, type="l", ylim=c(0,1),
      main="Density (first row) and Distribution Function (second row)",
      xlab=ll$xlab, ylab="CDF", cex.main=1.5)
    lines(xxx, eta0*F0(xxx), col=cols[1], lwd=2, lty=3)
    lines(xxx, (1-eta0)*FA(xxx), col=cols[2], lwd=2)
    
    # DEBUG show cutoff in green line
    #lines(c(x0,x0),c(0,1), col=3)

    plot(xxx, Fdr(xxx), type="l", lwd=2, ylim=c(0,1),
      main="(Local) False Discovery Rate", ylab="Fdr and fdr",
      xlab=ll$xlab, lty=3, cex.main=1.5)
    lines(xxx, fdr(xxx), lwd=2)
    if (eta0 > 0.98) 
      pos2 = "bottomleft" else pos2="topright"
    legend(pos2, 
      c("fdr (density-based)", "Fdr (tail area-based)"), 
      lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)

    # DEBUG show cutoff in green line
    #lines(c(x0,x0),c(0,1), col=3)

    par(mfrow=c(1,1))

    rm(ax)
  }

  if(verbose) cat("\n")

  return(result)
}


#####

## create labels for plots
pvt.plotlabels <- function(statistic, scale.param, eta0)
{
   if (statistic=="pvalue")
   {
     main = paste("Type of Statistic: p-Value (eta0 = ", round(eta0, 4), ")", sep="")
     xlab ="1-pval"
   }

   if (statistic=="studentt")
   {
     df = scale.param 
     main = paste("Type of Statistic: t-Score (df = ", round(df,3),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(t)"
   }

   if (statistic=="normal")
   {
     sd = scale.param 
     main = paste("Type of Statistic: z-Score  (sd = ", round(sd,3),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(z)"
   }

   if (statistic=="correlation")
   {
     kappa =scale.param      
     main = paste("Type of Statistic: Correlation (kappa = ", round(kappa,1),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(r)"
   }

   return(list(main=main, xlab=xlab))
}

