### dcor0.R  (2007-01-09)
###
###    Distribution of the Correlation Coefficient (rho=0) 
###    and Related Functions
###    
###
### Copyright 2003-07 Korbinian Strimmer
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



# density function
dcor0 <- function (x, kappa, log = FALSE)
{
  t <- r2t(x, kappa)
  df <- kappa-1
  vol <- sqrt(df)*(1-x^2)^(-3/2)
  if (log) 
    d <-  dt(t, df=df, log=log) + log(vol)
  else
    d <-  dt(t, df=df, log=log) * vol
 
  return(d) 
}

# distribution function
pcor0 <- function(q, kappa, lower.tail=TRUE, log.p=FALSE)
{
  t <- r2t(q, kappa)
  df <- kappa-1
  p <- pt(t, df=df, lower.tail = lower.tail, log.p = log.p)

  return(p)
}

# quantile function
qcor0 <- function(p, kappa, lower.tail=TRUE, log.p=FALSE)
{
  df <- kappa-1
  r <- t2r(qt(p, df=df, lower.tail = lower.tail, log.p = log.p), df)
  
  return(r)
}

# random number generator
rcor0 <- function(n, kappa)
{
  df <- kappa-1
  r <- t2r(rt(n, df),df)
  
  return(r)
}


### conversion from r to t statistic (and vice versa)

r2t <- function(r, kappa)
{
  t = r*sqrt((kappa-1)/(1-r*r)) 

  return(t) # df = kappa-1
} 

t2r <- function(t, df)
{
  r = t/sqrt(t*t+df)  

  return(r) # kappa = df+1
} 

