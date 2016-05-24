# Copyright (C) 2003 Roderick D. Ball
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
# Or see www.gnu.org/copyleft/gpl.htlm


smspline <- function(formula, data){
# Copyright (C) 2003 Roderick D. Ball
# generate the Z matrix
#browser()
if(is.vector(formula)){
  x <- formula
}else{
  if(missing(data)){
    mf <- model.frame(formula)
  }else{
    mf <- model.frame(formula,data)
  }
  if(ncol(mf) != 1)stop("formula can have only one variable")
  x <- mf[,1]
}

x.u <- sort(unique(x))
Zx <- smspline.v(x.u)$Zs
Zx[match(x,x.u),]
}

smspline.v <-
function(time)
{
# Copyright (C) 2003 Roderick D. Ball
### generate the X and Z matrices for a mixed model
### formulation of a smoothing spline as per Verbyla'a notes 5.3.2
### Verbyla, A.P. Mixed models for prectitioners. University of Adelaide and
### The South Australian Research and Development Institute, 1999, 115pp.
### limitation: knot points identical to data
### smoothing penalty \lambda_s \int g''(t) dt
### lambda_s = sigma^2/sigma^2_s
### y = X_s beta_s + Z_s u_s + e
### X_s = [1 | t]
### Z_s = Q (t(Q) %*%Q)^-1
### u_s ~ N(0, G_s)
### y[i] = g(t[i]) + e[i]
### let u_s = L v_s
### so cov(u) = G_s = L L'
	t1 <- sort(unique(time))
	p <- length(t1)
	h <- diff(t1)
	h1 <- h[1:(p - 2)]
	h2 <- h[2:(p - 1)]
	Q <- matrix(0, nrow = p, ncol = p - 2)
	Q[cbind(1:(p - 2), 1:(p - 2))] <- 1/h1
	Q[cbind(1 + 1:(p - 2), 1:(p - 2))] <- -1/h1 - 1/h2
	Q[cbind(2 + 1:(p - 2), 1:(p - 2))] <- 1/h2
	Gs <- matrix(0, nrow = p - 2, ncol = p - 2)
	Gs[cbind(1:(p - 2), 1:(p - 2))] <- 1/3 * (h1 + h2)
	Gs[cbind(1 + 1:(p - 3), 1:(p - 3))] <- 1/6 * h2[1:(p - 3)]
	Gs[cbind(1:(p - 3), 1 + 1:(p - 3))] <- 1/6 * h2[1:(p - 3)]
	Gs
	# Zus = Q %*% (t(Q)%*%Q)^-1
	Zus <- t(solve(t(Q) %*% Q, t(Q)))
	#R _ choleski(Gs)
	R <- chol(Gs, pivot = FALSE)
#Splus
#	if(attr(R, "rank") < nrow(R))
#		stop("singular G matrix")
        tol <- max(1e-12,1e-8*mean(diag(R)))
        if(sum(abs(diag(R))) < tol)
		stop("singular G matrix")
	Zvs <- Zus %*% t(R)
	list(Xs = cbind(rep(1, p), t1), Zs = Zvs, Q = Q, Gs = Gs, R = R)
}

approx.Z <- function(Z,oldtimes,newtimes){
# Copyright (C) 2003 Roderick D. Ball
# linear interpolation of Z matrix by column.
oldt.u <- sort(unique(oldtimes))
if(length(oldt.u) !=length(oldtimes) ||
   any(oldt.u != oldtimes)){
   Z <- Z[match(oldt.u,oldtimes),]
   oldtimes <- oldt.u
}
apply(Z,2,function(u,oldt,newt){
   approx(oldt,u,xout=newt)$y},
   oldt=oldtimes, newt=newtimes)
}

