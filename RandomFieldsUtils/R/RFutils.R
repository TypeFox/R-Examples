## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



.struve <- function(x, nu, sign, expon.scaled) {
  storage.mode(x) <- "double"
  storage.mode(nu) <- "double"
  storage.mode(expon.scaled) <- "logical"
  storage.mode(sign) <- "double"
#  res <- double(max(length(x), length(nu)))
 .Call("struve", x, nu, sign, expon.scaled)
}


struveH <- function(x, nu)  .struve(x, nu, -1, FALSE)
struveL <- function(x, nu, expon.scaled=FALSE)  .struve(x, nu, 1, expon.scaled)
I0L0 <- function(x) {
  storage.mode(x) <- "double"
#  res <- double(length(x))
  .Call("I0ML0", x)
}


solvePosDef <- function(a, b=NULL, logdeterminant=FALSE) {
  if (logdeterminant) {
    logdet <- double(1)
    res <- .Call("solvePosDef", a, b, logdet)
    return(list(inv=res, logdet=logdet))
  } else {
    .Call("solvePosDef", a, b, double(0))
  }
}


Print <- function(..., digits=6, empty.lines=2) { #
  ## ?"..1"
#  print(..1)
#  print(substitute(..1))
#   print(missing(..100))
   
  max.elements <- 99
  l <- list(...)
  n <- as.character(match.call())[-1]
  cat(paste(rep("\n", empty.lines), collapse="")) #
  for (i in 1:length(l)) {
    cat(n[i]) #
    if (!is.list(l[[i]]) && is.vector(l[[i]])) {
      L <- length(l[[i]])
      if (L==0) cat(" = <zero>")#
      else {
        cat(" [", L, "] = ", sep="")
        cat(if (is.numeric(l[[i]]))
            round(l[[i]][1:min(L , max.elements)], digits=digits)#
            else l[[i]][1:min(L , max.elements)]) #
        if (max.elements < L) cat(" ...")
      }
    } else {
       if (is.list(l[[i]])) {
        cat(" =  ") #
        str(l[[i]], digits.d=digits) #
      } else {
        cat(" =")
        if (length(l[[i]]) <= 100 && FALSE) {
          print(if (is.numeric(l[[i]])) round(l[[i]], digits=digits)#
                else l[[i]])
        } else {
          if (length(l[[i]]) > 1 && !is.vector(l[[i]]) && !is.matrix(l[[i]])
              && !is.array(l[[i]])) cat("\n")
          str(l[[i]]) #
        }
      }
    }
    cat("\n")
  }
}

