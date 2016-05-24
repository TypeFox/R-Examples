### halfnorm.R  (2007-01-01)
###
###     Half-Normal Distribution
###
### Copyright 2006-2007 Korbinian Strimmer 
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



dhalfnorm <- function(x, theta=sqrt(pi/2), log=FALSE)
{
  sd.norm = theta2sd(theta)
  
  if (log)
    d = ifelse(x<0, -Inf, log(2)+dnorm(x, mean=0, sd=sd.norm, log=TRUE))
  else
    d = ifelse(x<0, 0,  2*dnorm(x, mean=0, sd=sd.norm)) 

  return(d)
}

phalfnorm <- function(q, theta=sqrt(pi/2), lower.tail=TRUE, log.p=FALSE)
{
  sd.norm = theta2sd(theta)
  
  p = ifelse(q < 0, 0, 2*pnorm(q, mean=0, sd=sd.norm)-1)
  if (lower.tail == FALSE) p = 1-p
  if(log.p) p = log(p)
  
  return( p )
}

qhalfnorm <- function(p, theta=sqrt(pi/2), lower.tail=TRUE, log.p=FALSE)
{
  sd.norm = theta2sd(theta)
  
  if (log.p) p = exp(p)
  if (lower.tail == FALSE) p = 1-p
  q = ifelse(p < 0, NaN, qnorm((p+1)/2, mean=0, sd=sd.norm))
  
  return(q)
}

rhalfnorm <- function(n, theta=sqrt(pi/2))
{
  sd.norm = theta2sd(theta)
  return( abs(rnorm(n, mean=0, sd=sd.norm)) )
}


# conversion between standard deviation of normal distribution
# and theta parameter of corresponding half-normal distribution

sd2theta <- function(sd)
{
  return(sqrt(pi/2)/sd)
}

theta2sd <- function(theta)
{
  return(sqrt(pi/2)/theta)
}

