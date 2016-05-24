### monoreg.R  (2014-04-15)
###
###     Monotone Regression
###
### Copyright 2006-2014 Korbinian Strimmer 
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


# monotone regression

monoreg = function(x, y=NULL, w=rep(1, length(x)), 
   type=c("isotonic", "antitonic"))
{
    # get x-y coordinates
    xy = xy.coords(x,y)
    x = xy$x
    y = xy$y
    
    
    # remove duplicated x values
    xvals = unique(x)
    lx = length(xvals)
    
    if (lx == length(x)) # x values are all unique
    {
       wts   = w
       yvals = y 
    }
    else # merge duplicated x values
    {   
       wts = rep(NA, lx)
       yvals = rep(NA, lx)
       for (i in 1:lx)
       {
         idx = which(xvals[i]==x)
	 
	 if (length(idx) > 1)
	 {
	    warning("Duplicated x value (x=", xvals[i], ") detected!\nThe corresponding weights and y values will be merged.")
	    wi = w[idx]
            wts[i] = sum(wi) # add all weights
            yvals[i] = sum(y[idx]*wi)/wts[i]  # weighted mean
	 }
	 else
	 {
	    wts[i] = w[idx]
	    yvals[i] = xy$y[idx]
	 }
	 
       }
    } 
    
    # sort
    
    ord = order(xvals)
    xvals = xvals[ord]
    yvals = yvals[ord]
    wts = wts[ord]
       
    # perform monotonic regression      
          
    type = match.arg(type)
        
    if (type=="isotonic")
    {
        yf = pvt.isoMean(yvals, wts)
    }
    
    if (type=="antitonic")
    {
	yf = -pvt.isoMean(-yvals, wts)
    }
	
    # output results
    
    out = list(x=xvals, 
               y=yvals, 
	       w=wts,
               yf=yf,
	       type=type,
	       call = match.call())
    class(out) = c("monoreg")
    
    return(out)
}

fitted.monoreg = function(object, ...) object$yf
residuals.monoreg = function(object, ...) object$y - fitted(object)


plot.monoreg = function(x, main, main2,
     plot.type = c("single", "row.wise", "col.wise"), ...)
{
  plot.type = match.arg(plot.type)
  if (plot.type=="row.wise") par(mfrow=c(2,1))
  if (plot.type=="col.wise") par(mfrow=c(1,2))
  
  
  if (missing(main))
  {
     if (x$type=="isotonic")
        main = paste("Isotonic Regression:", deparse(x$call))
     else
        main = paste("Antitonic Regression:", deparse(x$call))
  }

  plot(x$x, x$y, main=main, xlab="x", ylab="y", ...)
  lines(x$x, x$yf, col=4, lty=3)
  points(x$x, x$yf, col=4, pch=20)
  
  if (plot.type != "single")
  { 
  
    if (missing(main2))
    {
      if(x$type == "isotonic")
        main2="Cumulative Sum Diagram and Greatest Convex Minorant"
      else
        main2="Cumulative Sum Diagram and Least Concave Majorant"  
    }
  
    G = c(0, cumsum(x$y*x$w))
    M = c(0, cumsum(x$yf*x$w))
    W = c(0, cumsum(x$w))
  
    plot(W, G, main=main2, 
       xlab="cumsum(w)", ylab="cumsum(y*w)", type="l")
    points(W, G)
    lines(W, M, col=4, lty=3)
    points(W, M, col=4, pch=20)
  
    par(mfrow=c(1,1))
  }
}


######## internal function ###################

pvt.isoMean = function(y, w)
{
  # Input:	y: measured values in a regression setting
  #		w: weights
  # Output: 	vector containing estimated (isotonic) values

  n = length(y)

  if(n == 1)
  {
    return(y)
  }
  else
  {
    ghat = .C("C_isomean",
            as.double(y),
            as.double(w),
            as.integer(n),
            ghat=double(n), PACKAGE="fdrtool")$ghat

    return(ghat)
  }
}
