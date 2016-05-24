# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"calc.npar" <- 
  function(object)
{

  ##attach(object)

  ## Thetas
  if(!any(is.null(object$thetas)||is.null(object$thetas))) nth <- length(object$thetas)

  ## SE of the thetas
  #if(!any(is.null(object$sethetas)||is.na(object$sethetas))) {
  if(!any(is.null(object$sethetas))){
    nseth <- length(object$sethetas[!is.na(object$sethetas)])
  }
  else {
    nseth <- 0
  }

  if(!any(is.null(object$omega)||is.na(object$omega))) {
    nom <- 0
    for(i in 1:length(object$omega)) {
      ## This only gives back the number of non-zero elements.
      ## Handling of initial zeros is done in plotsum2
      sel <- object$omega[[i]] != 0
      if(!any(sel==TRUE)){
        nom <- nom + 1
      }
      ##nom <- length(omega[[i]]) + nom
      nom <- length(object$omega[[i]][sel]) + nom
    }
  }
  ## SE of the omegas
  #if(!any(is.null(object$seomegas)||is.na(object$seomegas))) {
  if(!any(is.null(object$seomegas))) {
    nseom <- 0
    for(i in 1:length(object$seomegas)) {
      ## This only gives back the number of non-zero elements.
      ## Handling of initial zeros is done in plotsum2
      sel <- object$seomegas[[i]] != 0
      sel2 <- !is.na(object$seomegas[[i]])
      sel3 <- sel & sel2
      #if(!any(sel==TRUE)){
      #  nom <- nom + 1
      #}
      ##nom <- length(omega[[i]]) + nom
      nseom <- length(object$seomegas[[i]][sel3]) + nseom
    }
  }
  else {
    nseom <- 0
  }
  ## Sigmas
  nsi <- 0
  if(!any(is.null(object$sigma)||is.na(object$sigma))) {
    for(i in 1:length(object$sigma)) {
      sel <- object$sigma[[i]] != 0
      nsi <- length(object$sigma[[i]][sel]) + nsi
    }
  }
  ## SE of the sigmas
  #if(!any(is.null(object$sesigmas)||is.na(object$sesigmas))) {
  if(!any(is.null(object$sesigmas))) {
    nsesi <- 0
    for(i in 1:length(object$sesigmas)) {
      ## This only gives back the number of non-zero elements.
      ## Handling of initial zeros is done in plotsum2
      sel <- object$sesigmas[[i]] != 0
      sel2 <- !is.na(object$sesigmas[[i]])
      sel3 <- sel & sel2
      #if(!any(sel==TRUE)){
      #  nom <- nom + 1
      #}
      ##nom <- length(omega[[i]]) + nom
      nsesi <- length(object$sesigmas[[i]][sel3]) + nsesi
    }
  }
  else {
    nsesi <- 0
  }
  npar <- nth + nom + nsi
  if(length(nseth) > 0 || length(nseom) > 0 || length(nsesi) > 0) {
    ret.list <- list(npar = npar, nth = nth, nseth = nseth, nom = 
                     nom, nseom = nseom, nsi = nsi, nsesi = nsesi)
  }
  else {
    ret.list <- list(npar = npar, nth = nth, nom = nom, nsi = nsi)
  }
  #invisible(ret.list)
  return(ret.list)
}
