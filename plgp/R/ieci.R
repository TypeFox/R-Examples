#******************************************************************************* 
#
# Particle Learning of Gaussian Processes
# Copyright (C) 2010, University of Cambridge
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
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## calc.eis:
##
## wrapper used to calculate the ECIs (expected conditional improvements)
## in C, and then calculate the IECI by a (weighted) mean

calc.eis <- function(tmat, fmin, w=NULL)
  {
    n <- nrow(tmat)

    return(.C("calc_eis_R",
              tmat = as.double(t(tmat)),
              n = as.integer(n),
              fmin = as.double(fmin),
              bw = as.integer(length(w)),
              w = as.double(w),
              eis = double(n),
              PACKAGE="plgp")$eis)
  }


## calc.vars:
##
## wrapper used to calculate the predictive variances

calc.vars <- function(tmat, w=NULL)
  {
    if(is.null(w)) return(tmat$df*tmat$s2/(tmat$df-2))
    else return(w*tmat$df*tmat$s2/(tmat$df-2))
  }


## calc.iecis:
##
## wrapper used to calculate the ECIs (expected conditional improvements)
## in C, and then calculate the IECI by a (weighted) mean

calc.iecis <- function(ktKik, k, Xcand, X, Ki, Xref, d, g, s2p, phi,
                       badj, tm, tdf, fmin, w, verb)
  {
    m <- length(ktKik)
    I <- nrow(Xcand)

    return(.C("calc_iecis_R",
              ktKik = as.double(ktKik),
              m = as.integer(m),
              k = as.double(k),
              n = as.integer(nrow(k)),
              Xcand = as.double(t(Xcand)),
              I = as.integer(I),
              col = as.integer(ncol(Xcand)),
              X = as.double(t(X)),
              Ki = as.double(Ki),
              Xref = as.double(t(Xref)),
              d = as.double(d),
              dlen = as.integer(length(d)),
              g = as.double(g),
              s2p = as.double(s2p),
              phi = as.double(phi),
              badj = as.double(badj),
              tm = as.double(tm),
              tdf = as.integer(tdf),
              fmin = as.double(fmin),
              w = as.double(w),
              w.null = as.integer(is.null(w)),
              verb = as.integer(verb),
              iecis = double(I),
              PACKAGE="plgp")$iecis)
  }


## calc.alcs:
##
## wrapper used to calculate the ALCs in C

calc.alcs <- function(k, Xcand, X, Ki, Xref, d, g, s2p, phi,
                      badj, tdf, w, verb)
  {
    m <- nrow(Xref)
    I <- nrow(Xcand)

    return(.C("calc_alcs_R",
              m = as.integer(m),
              k = as.double(k),
              n = as.integer(nrow(k)),
              Xcand = as.double(t(Xcand)),
              I = as.integer(I),
              col = as.integer(ncol(Xcand)),
              X = as.double(t(X)),
              Ki = as.double(Ki),
              Xref = as.double(t(Xref)),
              d = as.double(d),
              dlen = as.integer(length(d)),
              g = as.double(g),
              s2p = as.double(s2p),
              phi = as.double(phi),
              badj = as.double(badj),
              tdf = as.integer(tdf),
              w = as.double(w),
              w.null = as.integer(is.null(w)),
              verb = as.integer(verb),
              alcs = double(I),
              PACKAGE="plgp")$alcs)
  }


## calc.ieci:
##
## wrapper used to calculate the ECIs (expected conditional improvements)
## in C, and then calculate the IECI by a (weighted) mean

calc.ieci <- function(ktKik, k, x, X, Ki, Xref, d, g, s2p, phi,
                     badj, tm, tdf, fmin, w)
  {
    m <- length(ktKik)

    return(.C("calc_ieci_R",
              ktKik = as.double(ktKik),
              m = as.integer(m),
              k = as.double(k),
              n = as.integer(nrow(k)),
              x = as.double(x),
              col = as.integer(length(x)),
              X = as.double(t(X)),
              Ki = as.double(Ki),
              Xref = as.double(t(Xref)),
              d = as.double(d),
              dlen = as.integer(length(d)),
              g = as.double(g),
              s2p = as.double(s2p),
              phi = as.double(phi),
              badj = as.double(badj),
              tm = as.double(tm),
              tdf = as.integer(tdf),
              fmin = as.double(fmin),
              w = as.double(w),
              ieci = double(1),
              PACKAGE="plgp")$ieci)
  }


## calc.ecis:
##
## wrapper used to calculate the ECIs (expected conditional improvements)
## in C

calc.ecis <- function(ktKik, k, x, X, Ki, Xref, d, g, s2p, phi,
                     badj, tm, tdf, fmin)
  {
    m <- length(ktKik)

    ## current version writes over ktKiK 
    return(.C("calc_ecis_R",
              ktKik = as.double(ktKik),
              m = as.integer(m),
              k = as.double(k),
              n = as.integer(nrow(k)),
              x = as.double(x),
              col = as.integer(length(x)),
              X = as.double(t(X)),
              Ki = as.double(Ki),
              Xref = as.double(t(Xref)),
              d = as.double(d),
              dlen = as.integer(length(d)),
              g = as.double(g),
              s2p = as.double(s2p),
              phi = as.double(phi),
              badj = as.double(badj),
              tm = as.double(tm),
              tdf = as.integer(tdf),
              fmin = as.double(fmin),
              PACKAGE="plgp")$ktKik)
  }


## calc.ktKik.x:
##
## no longer used; testing interface for "_R" C function

calc.ktKik.x <- function(ktKik, k, g, mui, kxy)
  {
    m <- length(ktKik)

    ## current version writes over ktKiK 
    r <- .C("calc_ktKikx_R",
            ktKik = as.double(ktKik),
            m = as.integer(m),
            k = as.double(k),
            n = as.integer(nrow(k)),
            g = as.double(g),
            mui = as.double(mui),
            kxy = as.double(kxy),
            PACKAGE="plgp")
   
    return(r$ktKik)
  }


## calc2.ktKik.x:
##
## no longer used; testing interface for "_R" C function

calc2.ktKik.x <- function(ktKik, k, x, X, Ki, Xref, d, g)
  {
    m <- length(ktKik)

    ## current version writes over ktKiK 
    r <- .C("calc2_ktKikx_R",
            ktKik = as.double(ktKik),
            m = as.integer(m),
            k = as.double(k),
            n = as.integer(nrow(k)),
            x = as.double(x),
            col = as.integer(length(x)),
            X = as.double(t(X)),
            Ki = as.double(Ki),
            Xref = as.double(t(Xref)),
            d = as.double(d),
            dlen = as.integer(length(d)),
            g = as.double(g),
            PACKAGE="plgp")
    
    return(r$ktKik)
  }

