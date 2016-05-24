
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:             DESCRIPTION:
#  adapt                 Integrates over a two dimensional unit square
################################################################################


adapt <-
  function(ndim=NULL, lower, upper, functn, ...)
{
  ans <- cubature::adaptIntegrate(
    f=functn, lowerLimit=lower, upperLimit=upper, ...)
  
  # Return Value
  ans
}


# -----------------------------------------------------------------------------
# replaced by the cubature adaptIntegrate package
# you will find the code in the deprecated folder


# adapt <- function (ndim, lower, upper, minpts = 100, maxpts = NULL, 
#                   functn, eps = 0.01, ...)


# Title: adapt -- multidimensional numerical integration
# Package: adapt
# Version: 1.0-4
# Author: FORTRAN by Alan Genz, 
#     S by Mike Meyer, R by Thomas Lumley and Martin Maechler
# Description: Adaptive Quadrature in up to 20 dimensions
# Depends:
# License: Unclear (Fortran) -- code in Statlib's ./S/adapt
# Maintainer: Thomas Lumley <tlumley@u.washington.edu>
# Packaged: Fri Apr 20 11:38:07 2007; thomas


# [from Statlib's original  http://lib.stat.cmu.edu/S/adapt ]
# This code contains an S function and supporting C and Fortran code for
# adaptive quadrature.  The underlyling fortran code is purported to
# work in from 2 to 20 dimensions.  The code is set up to dynamically
# load from a central library area.  If you can not do dynamic loading,
# you may need to build a staticly loaded version.  The adapt S function
# calls load.if.needed to do the dynamic loading.  You will have to
# change the functions used here (probably to call library.dynam).
# S code written by Michael Meyer (mikem@andrew.cmu.edu).
# October, 1989.


# 2002-03-14  Martin Maechler  <maechler@stat.math.ethz.ch>
# * DESCRIPTION (Version): 1.0-3 --> CRAN
# * R/adapt.R (adapt): use defaults for minpts, maxpts, eps;
#   more logical maxpts default (for ndim >= 7) using rulcls
# * man/adapt.Rd: extended example
# 2002-03-13  Martin Maechler  <maechler@stat.math.ethz.ch>
# * DESCRIPTION (Version): 1.0-2
# * man/adapt.Rd: indentation, using \code{.}, etc;
#   example also tries p=5 dimensions
# * R/adapt.R: clean up (spaces)
# 2002-01-09  Martin Maechler  <maechler@stat.math.ethz.ch>
# * R/adapt.R: do not use .Alias anymore
# 2001-06-29  Thomas Lumley <tlumley@u.washington.edu>
# * move (improved!) integrate() into base, using .Call() etc.


# Message-ID: <4AD7A74B.3020108@math.wsu.edu>
# Date: Thu, 15 Oct 2009 15:50:51 -0700
# From: Alan Genz <genz@math.wsu.edu>
# User-Agent: Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.21) 
#     Gecko/20090402 SeaMonkey/1.1.16
# MIME-Version: 1.0
# To: Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
# CC: Alan C Genz <alangenz@wsu.edu>
# Subject: Re: adapt
# References: <4AD3032B.4090801@itp.phys.ethz.ch>
# In-Reply-To: <4AD3032B.4090801@itp.phys.ethz.ch>
# Content-Type: text/plain; charset=ISO-8859-1; format=flowed
# Content-Transfer-Encoding: 7bit
# Status:  O
# Dear Prof. Wuertz,
# Thank you for your message and your interest in my adaptive integration
# Fortran code. I would be pleased if you included my code in your open
# source R package under the Gnu GPL2 license. You have my permission to 
# do this.
# Sincerely,
# Alan Genz


################################################################################

