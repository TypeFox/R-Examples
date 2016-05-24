### dyn.scale.R  (2008-11-14)
###
###    Dynamical Scale, Moments, and Weights
###
### Copyright 2005-2008 Rainer Opgen-Rhein and Korbinian Strimmer
###
###
###
### This file is part of the `GeneTS' library for R and related languages.
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




# mean and variance of longitudinal data matrix
dyn.moments = function(x)
{
  w = dyn.weights(x)
  
  return( wt.moments(x, w) )
}

# scale longitudinal data matrix
dyn.scale = function(x, center=TRUE, scale=TRUE)
{
  w = dyn.weights(x)
  x = wt.scale(x, w, center, scale)

  return( x )
}

# compute weights from 
dyn.weights = function(x)
{
  if ( is.longitudinal(x) )
  {
    tvec = get.time.repeats(x)$time
    rvec = get.time.repeats(x)$repeats
        
    tw = time2weights(tvec)
    w = rep(tw, rvec) / rep(rvec, rvec)
  }
  else if (is.matrix(x))
  {
    n = dim(x)[1]
    w = rep(1/n, n)
  }
  else
  {
     stop("x has to be a matrix, or a longitudinal object")
  }
  
  return(w)
}




#  turn time vector into weight vector
# 
# t:   a vector with the time points
# w:   the corresponding weight vector (sums up to 1)
#
# weights are assume a linear spline
# example: for 6 equidistant points the weights
#          are  c(0.1, 0.2, 0.2, 0.2, 0.2, 0.1) 

time2weights = function(t)
{
  # number of time points
  u = length(t) 

  if (u == 1) return(1)

  dt = diff(t)
  w = (c(0, dt) + c(dt, 0))/(2*sum(dt))

  return(w)
}



