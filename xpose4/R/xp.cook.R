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

"xp.cook"<-
  function(gam.object){
    
    ##assign(pos = 1, "data", gam.object$data)
    fit.s <- summary.glm(gam.object)
    fit.infl <- lm.influence(gam.object)
    R <- gam.object$R
    I <- t(R) %*% R
    Iinv <- fit.s$cov.unscaled
    ass <- gam.object$assign
    names(ass) <- names(gam.object$coefficients)
    D <- matrix(0, length(gam.object$residuals), length(ass))
    dimnames(D) <- list(names(gam.object$residuals), names(gam.object$coefficients))
    Dcoef <- scale(fit.infl$coefficients, center = gam.object$coefficients, scale= F)
    for(subname in names(ass)) {
      sub <- ass[[subname]]
      Dcoefi <- Dcoef[, sub, drop = F] %*% t(R[, sub, drop = F])
      denom <- I[sub, sub, drop = F] %*% Iinv[sub, sub, drop = F]
      denom <- sum(diag(denom)) * fit.s$dispersion
      D[, subname] <- apply(Dcoefi^2, 1, sum)/denom
    }
    ##remove("data",fr=0)
    D
  }
