### wt.scale.R  (2008-10-14)
###
###    Weighted Expectations and Variances
###    
###
### Copyright 2006-2008 Korbinian Strimmer
###
### This file is part of the `corpcor' library for R and related languages.
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


# in all the following functions, 
# w is a vector of weights with sum(w)=1


# mean and variance


# this exists already in R 
#weighted.mean = function(xvec, w)
#{
#  return( sum(w*xvec) )
#}

wt.var = function(xvec, w) 
{
  w = pvt.check.w(w, length(xvec))

  # bias correction factor
  h1 = 1/(1-sum(w*w))   # for w=1/n this equals the usual h1=n/(n-1)
       
  xc = xvec-weighted.mean(xvec, w)
  s2 = h1*weighted.mean(xc*xc, w)

  return( s2 ) 
}


wt.moments = function(x, w)
{
  x = as.matrix(x)
  w = pvt.check.w(w, nrow(x))
  # bias correction factor
  h1 = 1/(1-sum(w*w))   # for w=1/n this equals the usual h1=n/(n-1)
 
     
  # m = apply(x, 2, weighted.mean, w=w)
  m = colSums(w*x)  # same as above, but much faster
  
  # v = apply(x, 2, wt.var, w=w)
  v = h1*(colSums(w*x^2)-m^2) # same as above, but much faster
 
  
  # set small values of variance exactly to zero
  v[v < .Machine$double.eps] = 0
  
  return( list(mean=m, var=v) )
}


# scale using the weights
wt.scale = function(x, w, center=TRUE, scale=TRUE)
{
  x = as.matrix(x)
  w = pvt.check.w(w, nrow(x))
  
  # compute column means and variances
  wm = wt.moments(x, w)

  if (center==TRUE)
  {
      x = sweep(x, 2, wm$mean, "-")	
      attr(x, "scaled:center") = wm$mean
  } 
  
  if (scale==TRUE)
  {
      sc = sqrt(wm$var)
         
      x = sweep(x, 2, sc, "/")
      attr(x, "scaled:scale") = sc
      
      zeros = (sc == 0.0)
      x[,zeros] = 0
      
      if (any(zeros))
      {
        warning(paste(sum(zeros), "instances of variables with zero scale detected!"),
	 call. = FALSE)
      }  
  } 

  return(x)
}

