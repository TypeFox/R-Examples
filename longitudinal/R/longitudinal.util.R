

### longitudinal.util.R  (2005-05-07)
###
###    Utility functions for longitudinal data
###
### Copyright 2005 Korbinian Strimmer
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



###### PRIVATE ########


# indixes: where do time points end and begin (given repeats)
get.time.idx = function(r)
{
   cr = cumsum(c(1,r))
   start = cr[1:length(r)]
   end = cr[1:length(r)+1]-1
   
   idx = matrix(c(start,end),nrow=length(r),ncol=2)
   
   return( idx  )
}

###### PUBLIC ########

get.time.repeats = function(x)
{
  if (!is.longitudinal(x))
    stop("argument is not a longitudinal object")

  
  return( list(time=attr(x, "time"), repeats=attr(x, "repeats")) )
}



## utility functions

condense.longitudinal = function(x, s, func=median)
{
  if (!is.longitudinal(x))
    stop("argument is not a longitudinal object")


   if(missing(s))
     s = 1:dim(x)[2]

   tr = get.time.repeats(x)
   tp = length(tr$time)
   t.idx = get.time.idx(tr$repeats)

   mat = matrix(NA, nrow = tp, ncol = length(s))
   colnames(mat) = colnames(x)[s]

   for (i in 1:tp)
   {
     tt = (t.idx[i,1]):(t.idx[i,2])
     xtt = x[tt, s, drop=FALSE]  # keep results as matrix
     mat[i,] = apply(xtt, 2, func)
   }

   return(mat)
}


combine.longitudinal = function(x1, x2)
{
  # some basic checks
  if (!is.longitudinal(x1))
    stop("argument 1 is not a longitudinal object")
    
  if (!is.longitudinal(x2))
    stop("argument 2 is not a longitudinal object")

   if (dim(x1)[2] != dim(x2)[2])
     stop("different number of variables in both longitudinal objects")

   tr1 = get.time.repeats(x1)
   tr2 = get.time.repeats(x2)
   tidx1 = get.time.idx(tr1$repeats)
   tidx2 = get.time.idx(tr2$repeats)
   
   # combine the time series objects
   time3 = sort(union(tr1$time, tr2$time))   
   repeats3 = rep(0, length(time3))
     
   x3 = NULL
   for (i in 1:length(time3))
   {
      overlap1 = (tr1$time == time3[i])
      if (any(overlap1)) 
      {
        k1 = which(overlap1)
	t1 = seq(tidx1[k1,1], tidx1[k1,2])
	x3 = rbind(x3, x1[t1,])
        repeats3[i] = repeats3[i] + length(t1)
      } 
    
      overlap2 = (tr2$time == time3[i])
      if (any(overlap2)) 
      {
        k2 = which(overlap2)
	t2 = seq(tidx2[k2,1], tidx2[k2,2])
	x3 = rbind(x3, x2[t2,])
        repeats3[i] = repeats3[i] + length(t2)
      } 
   }
   
   return( as.longitudinal(x3, repeats=repeats3, time=time3) )
}



is.equally.spaced = function(x)
{
  if (!is.longitudinal(x))
    stop("argument is not a longitudinal object")


  t = get.time.repeats(x)$time
  
  if (length(t)==1) return(NA)
  
  dt = diff(t)
  if (all(dt[1] == dt))
    return(TRUE)
  else
    return(FALSE)
}


is.regularly.sampled = function(x)
{
  if (!is.longitudinal(x))
    stop("argument is not a longitudinal object")

  r = get.time.repeats(x)$repeats
  
  if (all(r[1] == r))
    return(TRUE)
  else
    return(FALSE)
}


has.repeated.measurements = function(x)
{
  if (!is.longitudinal(x))
    stop("argument is not a longitudinal object")
  
  if (length(get.time.repeats(x)$time) == dim(x)[1])
     return(FALSE)
  else
     return(TRUE)
}
