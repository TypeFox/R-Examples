#******************************************************************************* 
#
# Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
# Copyright (C) 2005, University of California
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
#
#*******************************************************************************


## tgp.choose.as:
##
## pick which type of "errors" to be returned, either for
## plotting purposes or for adaprive sampling purposes

"tgp.choose.as" <-
function(out, as)
{
  ## choose AS stats to plot

  ## default quantile diffs (as=NULL), or predictive variance (as="s2")
  if(is.null(as) || as == "s2" || as == "ks2") {

    ## use fulle data, XX & X
    X <- out$XX

    ## choose quantile diffs or s2
    if(is.null(as)) { 
      criteria <- c(out$Zp.q, out$ZZ.q)
      name <- "quantile diff (error)"
      if(!is.null(out$Zp.q)) X <- rbind(out$X, X)
    } else if(as == "ks2") {
      criteria = c(out$Zp.ks2, out$ZZ.ks2)
      name <- "kriging var"
      if(!is.null(out$Zp.ks2)) X <- rbind(out$X, X)
    } else {
      if(is.matrix(out$Zp.s2)) criteria <- c(diag(out$Zp.s2), diag(out$ZZ.s2))
      else criteria <- c(out$Zp.s2, out$ZZ.s2)
      name <- "pred var"
      if(!is.null(out$Zp.s2)) X <- rbind(out$X, X)
    }
    
  } else {

    ## only use predictive data
    X <- out$XX

    ## default choice is ALM stats (quantile diffs)
    criteria <- out$ZZ.q
    name <- "ALM stats"

    ## choose ALC or EGO stats
    if(as == "alc") {
      if(is.null(out$Ds2x)) cat("NOTICE: out$Ds2x is NULL, using ALM\n")
      else { criteria <- out$Ds2x; name <- "ALC stats" }
    } else if(as == "improv") {
      if(is.null(out$improv)) cat("NOTICE: out$improv is NULL, using ALM\n")
      else {
        criteria <- out$improv[,1];
        name <- paste("Improv stats (g=", out$g[1], ")", sep="")
      }
    } else if(as != "alm")
      warning(paste("as criteria \"", as, "\" not recognized; defaulting to \"alm\"",
                    sep=""))
  }

  ## there might be nothing to plot
  if(is.null(criteria)) stop("no predictive data, so nothing to plot")
  
  ## return
  return(list(X=X, criteria=criteria, name=name))
}


## tgp.choose.center:
##
## pick which type of center (mean, median, kriging mean, etc)
## to be returned, mostly for plotting purposes

"tgp.choose.center" <-
function(out, center)
{
  X <- out$XX

  ## check center description
  if(center != "mean" && center != "med"  && center != "km") {
    warning(paste("center = \"", center, "\" invalid, defaulting to \"mean\"\n",
                  sep=""))
    center <- "mean"
  }
  
  ## choose center as median or mean
  if(center == "med") {
    name <- "median";
    Z <- c(out$Zp.med, out$ZZ.med)
    if(!is.null(out$Zp.med)) X <- rbind(out$X, X)
  } else if(center == "km") {
    name <- "kriging mean";
    Z <- c(out$Zp.km, out$ZZ.km)
    if(!is.null(out$Zp.km)) X <- rbind(out$X, X)
  } else {
    name <- "mean";
    Z <- c(out$Zp.mean, out$ZZ.mean)
    if(!is.null(out$Zp.mean)) X <- rbind(out$X, X)
  }

  ## there might be nothing to plot
  if(is.null(Z)) stop("no predictive data, so nothing to plot")
  
  ## return
  return(list(X=X, Z=Z, name=name))
}
