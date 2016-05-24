### nullmodel.R  (2009-11-19)
###
###     Details on the FDR Null Model
###
### Copyright 2007-2009 Korbinian Strimmer 
###
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



# this function specifies all that is needed to define a null model
get.nullmodel = function(
  statistic=c("normal", "correlation", "pvalue", "studentt")
)
{
  statistic <- match.arg(statistic)

  ###### normal z-scores ######

  if (statistic=="normal")
  {
    # untruncated density
    f0 = function(x, param, log=FALSE)
    {
      return( dnorm(x, sd=param, log=log) )
    }

    # corresponding distribution function 
    F0 = function(x, param)
    {
      return( pnorm(x, sd=param) )
    }

    # interquartile range
    iqr = function(param)
    {
       return( qnorm(.75, sd=param)-qnorm(.25, sd=param) )
    }

    # parameter support
    get.support = function() return( c(1e-9,Inf) )
  }

  ###### p-values ######

  if (statistic=="pvalue")
  {
    # untruncated density
    f0 = function(x, param, log=FALSE)
    {
      if(log==TRUE) 
        return( rep(0, length(x)) )
      else 
        return( rep(1, length(x)) )
    }

    # corresponding distribution function 
    F0 = function(x, param)
    {
      return( x )
    }

    # interquartile range
    iqr = function(param)
    {
       return( 0.5 )
    }

    # parameter support
    get.support = NULL
  }

  ###### correlation coefficients ######

  if (statistic=="correlation")
  {
    # untruncated density
    f0 = function(x, param, log=FALSE)
    {
       return( dcor0(x, kappa=param, log=log) ) 
    }

    # corresponding distribution function 
    F0 = function(x, param)
    {
      return( pcor0(x, kappa=param) )
    }

    # interquartile range
    iqr = function(param)
    {
       return( qcor0(.75, kappa=param)-qcor0(.25, kappa=param) )
    }

    # parameter support
    get.support = function() return( c(3,1e9) )
  }

  ###### t scores ######

  if (statistic=="studentt")
  {
    # untruncated density
    f0 = function(x, param, log=FALSE)
    {
       return( dt(x, df=param, log=log) )
    }

    # corresponding distribution function 
    F0 = function(x, param)
    {
      return( pt(x, df=param) )
    }

    # interquartile range
    iqr = function(param)
    {
       return( qt(.75, df=param)-qt(.25, df=param) )
    }

    # parameter support
    get.support = function() return( c(1,1000) )
  }

  ###### corresponding p-values ######

  if (statistic=="pvalue")
  {
    # one-sided p-values
    get.pval = function(x, param)
    {
      #return( F0(x) )
      return( x )  
    }
  }
  else
  {
    # two-sided p-values
    get.pval = function(x, param)
    {
      ax = abs(x)

      # return pval=0 only if abs(x)=Inf
      return( ifelse(ax==Inf, 0, 
                pmax(.Machine$double.eps, 
                  2-2*F0(ax, param))) )   
    }
  }


  return(list(
         f0 = f0, 
         F0=F0,
         iqr=iqr,
         get.pval = get.pval,
         get.support = get.support
  ))
}
