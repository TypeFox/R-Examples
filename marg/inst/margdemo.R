## file marg/inst/margdemo.R, v 1.0.1 (2013-05-14)
##
##  Copyright (C) 2000-2013 Alessandra R. Brazzale 
##
##  This file is part of the "marg" package for R.  This program is 
##  free software; you can redistribute it and/or modify it under the 
##  terms of the GNU General Public License as published by the Free 
##  Software Foundation; either version 2 of the License, or (at your 
##  option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
##  MA 02111-1307 USA or look up the web page 
##  http://www.gnu.org/copyleft/gpl.html.
##
##  Please send any comments, suggestions or errors found to:
##  Alessandra R. Brazzale, Department of Statistics, University of
##  Padova, Via C. Battisti 241/243, 35121 Padova (PD), Italy.
##  Email: alessandra.brazzale@unipd.it 
##  Web: http://www.stat.unipd.it/~brazzale


## ===================================================================
##
##                         User-defined "rsm" families
##                                             
## ===================================================================


## This file gives two examples of how to define a new "rsm" family,
## that is, a family of distributions which is not included in the
## "rsm.distributions" object defined in the "marg" package of the
## "hoa" package bundle.  The two distributions we consider are
## obtained by contaminating the standard Gumbel distributions at the 
## 10% rate with:
##
## 1) a centered Gumbel distribution with scale parameter 10, i.e.
##     0.9*EV(0,1) + 0.1*EV(0,10);
## 2) a standard Gumbel distribution re-centered at the value log(4),
##    i.e. 0.9*logExp(1) + 0.1*logExp(4).

## NOTE: This functions will be used in the demonstration file
##       associated with the "sampling" package of the "hoa" bundle.



## Example 1) 0.9*EV(0,1) + 0.1*EV(0,10)
## =====================================

## You need:


## 1) A family generator function
##    ---------------------------

contam.EV <- function()
{
  make.family.rsm("contam.EV")
}


## 2) A `.distributions' object (of the same name than the generator function!)
##    -----------------

contam.EV.distributions <- structure(
                      .Data = list(
	g0 = function(y,...)
             { 
               dens <- 0.9*exp(-exp(y)+y) + 
                         0.1*exp(-exp(y/10)+y/10)/10 
               -log(dens)
             },
        g1 = function(y,...) 
             { 
               dens <- 0.9*exp(-exp(y)+y) + 
                         0.1*exp(-exp(y/10)+y/10)/10 
               dens.1 <- 0.9*exp(-exp(y)+y)*(1-exp(y)) + 
                           0.1*exp(-exp(y/10)+y/10)*
                             (1/10-exp(y/10)/10)/10
               -dens.1/dens
             }, 
        g2 = function(y,...)  
             {
               dens <- 0.9*exp(-exp(y)+y) + 
                         0.1*exp(-exp(y/10)+y/10)/10 
               dens.1 <- 0.9*exp(-exp(y)+y)*(1-exp(y)) + 
                           0.1*exp(-exp(y/10)+y/10)*
                             (1/10-exp(y/10)/10)/10
               dens.2 <- 0.9*(exp(-exp(y)+y)*(1-exp(y))^2-exp(-exp(y)+
                           2*y)) + 0.1*(exp(-exp(y/10)+y/10)*
                                     (1/10-exp(y/10)/10)^2 -
                         exp(-exp(y/10)+2*y/10)/10^2)/10
               -(dens.2*dens-(dens.1)^2)/dens^2
             } ),          
                      .Dim = c(3,1),
                      .Dimnames = list(c("g0","g1","g2"), 
                                       c("contam.EV")))


## 3) The corresponding d-, p-, q- and r- distribution functions 
##    (if available)
##    ----------------------------------------------------------

dcontam.EV <- function(x)
{
     0.9*dweibull(exp(x), shape=1, scale=1)*exp(x) +
	0.1*dweibull(exp(x), shape=1/10, scale=1)*exp(x)
}


rcontam.EV <- function(n)                          
{
    if( is.na(n) )
        return(NA)
    val1 <- rweibull(n, shape=1, scale=1)
    val2 <- rweibull(n, shape=1/10, scale=1)
    alpha <- runif(n)
    log(ifelse(alpha < 0.1, val2, val1))
}


pcontam.EV <- function(q)
{
     0.9*pweibull(exp(q), shape=1, scale=1) +
	0.1*pweibull(exp(q), shape=1/10, scale=1)
}


## Example 2) 0.9*logExp(1) + 0.1*logExp(4)
## ========================================

## You need:


## 1) A family generator function
##    ---------------------------

contam.EV <- function()
{
  make.family.rsm("contam.EV")
}


## 2) A `.distributions' object (of the same name than the generator function!)
##    -----------------

contam.EV.distributions <- structure(
                      .Data = list(
	g0 = function(y,...)
             { 
               dens <- 0.9*exp(-exp(y)+y) + 0.1*exp(-exp(y)*4+y)*4 
               -log(dens)
             },
        g1 = function(y,...) 
             { 
               dens <- 0.9*exp(-exp(y)+y) + 0.1*exp(-exp(y)*4+y)*4 
               dens.1 <- 0.9*exp(-exp(y)+y)*(1-exp(y)) + 
                         0.1*exp(-exp(y)*4+y)*(1-exp(y)*4)*4
               -dens.1/dens
             }, 
        g2 = function(y,...)  
             {
               dens <- 0.9*exp(-exp(y)+y) + 0.1*exp(-exp(y)*4+y)*4 
               dens.1 <- 0.9*exp(-exp(y)+y)*(1-exp(y)) + 
                         0.1*exp(-exp(y)*4+y)*(1-exp(y)*4)*4
               dens.2 <- 0.9*(exp(-exp(y)+y)*(1-exp(y))^2-
                           exp(-exp(y)+2*y)) +
                         0.1*(exp(-exp(y)*4+y)*(1-exp(y)*4)^2 -
                             exp(-exp(y)*4+2*y)*4)*4
               ret <- -(dens.2*dens-(dens.1)^2)/dens^2
	       ret[is.na(ret)] <- 0
	       ret
             } ),          
                      .Dim = c(3,1),
                      .Dimnames = list(c("g0","g1","g2"), 
                                       c("contam.EV")))

	
## 3) The corresponding d-, p-, q- and r- distribution functions 
##    (if available)
##    ----------------------------------------------------------


dcontam.EV <- function(x)
{
     0.9*dweibull(exp(x), shape=1, scale=1)*exp(x) +
	0.1*dexp(exp(x), rate=4)*exp(x)
}

rcontam.EV <- function(n)                          
{
    if( is.na(n) )
        return(NA)
    val1 <- rweibull(n, shape=1, scale=1)
    val2 <- rexp(n, rate=4)
    alpha <- runif(n)
    log(ifelse(alpha < 0.1, val2, val1))
}

pcontam.EV <- function(q)
{
     0.9*pweibull(exp(q), shape=1, scale=1) +
	0.1*dexp(exp(q), rate=4)
}

