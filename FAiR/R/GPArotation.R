#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich
#
#     These functions are all slightly modified from those in GPArotation/R/GPArotation.R,
#     which is Copyright 2005-2006 Coen Baarnards and Robert Jennrich and licensed under
#     the GPL V2+ as of version 2008.05-1 of the GPArotation package.
#
#     FAiR is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     FAiR is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with FAiR.  If not, see <http://www.gnu.org/licenses/>.

vgQ.geomin <- function(L, delta = .01) {
	k <- ncol(L)
	p <- nrow(L)
	L2 <- L^2 + delta
	pro <- exp(rowSums(log(L2))/k)
	f <- sum(pro)
	return(f)
}

vgQ.quartimin <- function(L) {
	X <- L^2 %*% (!diag(TRUE, ncol(L)))
	f <- sum(L^2 * X)/4
	return(f)
}

vgQ.target <- function(L, Target) {
	f <- sum((L - Target)^2)
	return(f)
}

vgQ.pst <- function(L, Target) {
	f <- sum( (L - Target)^2, na.rm = TRUE)
	return(f)
}

vgQ.oblimax <- function(L) {
	f <- -(log(sum(L^4)) - 2 * log(sum(L^2)))
	return(f)
}

vgQ.simplimax <- function(L, k = nrow(L)) {
	Imat <- sign(L^2 <= sort(L^2)[k])
	f <- sum(Imat*L^2)
	return(f)
}

vgQ.bentler <- function(L) {
	L2 <- L^2
	M <- crossprod(L2)
	D <- diag(diag(M))
	f <- -(log(det(M))-log(det(D)))/4
	return(f)
}

vgQ.cf <- function(L, kappa = 0) {
	k <- ncol(L)
	p <- nrow(L)
	N <- matrix(1,k,k) - diag(k)
	M <- matrix(1,p,p) - diag(p)
	L2 <- L^2
	f1 <- (1 - kappa) * sum(diag(crossprod(L2, L2 %*% N)))/4
	f2 <- kappa * sum(diag(crossprod(L2, M %*% L2)))/4
	f <- f1 + f2
	return(f)
}

vgQ.infomax <- function(L) {
	k <- ncol(L)
	p <- nrow(L)
	S <- L^2
	s <- sum(S)
	s1 <- rowSums(S)
	s2 <- colSums(S)
	E <- S / s
	e1 <- s1 / s
	e2 <- s2 / s
	Q0 <- sum(-E * log(E))
	Q1 <- sum(-e1 * log(e1))
	Q2 <- sum(-e2 * log(e2))
	f <- log(k) + Q0 - Q1 - Q2
	return(f)
}

vgQ.mccammon <- function(L) {
	k <- ncol(L)
	p <- nrow(L)
	S <- L^2
	M <- matrix(1,p,p)
	s2 <- colSums(S)
	P <- S / matrix(rep(s2,p), ncol = k, byrow=TRUE)
	Q1 <- -sum(P * log(P))
	H <- -(log(P) + 1)
	R <- M %*% S
	G1 <- H / R - M %*% (S * H / R^2)
	s <- sum(S)
	p2 <- s2 / s
	Q2 <- -sum(p2 * log(p2))
	f <- log(Q1) - log(Q2)
	return(f)
}

vgQ.oblimin <- function(L, gam = 0) {
	X <- L^2 %*% (!diag(TRUE, ncol(L))) 
	if(0 != gam) {
		p <- nrow(L)
		X <- (diag(1, p) - matrix(gam / p, p, p)) %*% X
	}
	f <- sum(L^2 * X) / 4
	return(f)
}

NormalizingWeight <- function (A, normalize = FALSE) {
	if(is.function(normalize)) normalize <- normalize(A)
	else if(is.character(normalize)) {
		normalize <- match.arg(tolower(normalize), c("kaiser", "cureton-mulaik"))
		if(normalize == "kaiser") normalize <- sqrt(rowSums(A^2))
		else if(normalize == "cureton-mulaik") {
			reduced <- tcrossprod(A)
			pc <- svd(reduced, 1, 0)$u
			h  <- sqrt(diag(reduced))
			a  <- abs(pc / h)
			k  <- ncol(A)^(-1/2)
			w  <- ( acos(k) -  acos(a) ) * pi / (2 * 
				acos(k) - pi * (a < k))
			w  <- cos(w)^2
			normalize <- h / w
		}
	}
	if(nrow(A) != length(normalize)) stop("normalize length wrong")
	return(normalize)
}
