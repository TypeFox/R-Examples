## This file is part of the 'agop' library.
##
## Copyright 2013 Marek Gagolewski, Anna Cena
##
## 'agop' is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## 'agop' is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with 'agop'. If not, see <http://www.gnu.org/licenses/>.



#' @title
#' WAM and OWA Operators
#' 
#' @description
#' Computes the Weghted Arithmetic Mean or the
#' Ordered Weighted Averaging aggregation operator.
#' 
#' @details
#' The OWA operator is given by
#' \deqn{
#' \mathsf{OWA}_\mathtt{w}(\mathtt{x})=\sum_{i=1}^{n} w_{i}x_{\{i\}}
#' }{
#' OWA_w(x) = sum_i(w_i * x_{i})
#' }
#' where \eqn{x_{\{i\}}}{x_{i}} denotes the \eqn{i}-th greatest
#' value in \code{x}.
#' 
#' The WAM operator is given by
#' \deqn{
#' \mathsf{WAM}_\mathtt{w}(\mathtt{x})=\sum_{i=1}^{n} w_{i}x_{i}
#' }{
#' WAM_w(x) = sum_i(w_i * x_i)
#' }
#' 
#' If the elements of \code{w} does not sum up to \eqn{1}, then
#' they are normalized and a warning is generated.
#' 
#' Both functions return the ordinary arithmetic mean by default.
#' Special cases of OWA include the trimmed mean (cf. \code{\link{mean}})
#' and winsorized mean.
#' 
#' There is a strong connection between the OWA operators
#' and the Choquet integrals.
#' 
#' @param x numeric vector to be aggregated
#' @param w numeric vector of the same length as \code{x}, with elements in \eqn{[0,1]},
#' and such that \eqn{\sum_i w_i=1}{sum(x)=1}; weights
#' @return single numeric value
#' 
#' @rdname owa
#' @name owa
#' @export
#' @family aggregation_operators
#' @references
#' Yager R.R., On ordered weighted averaging aggregation operators in multicriteria decision making, IEEE Transactions on Systems, Man, and Cybernetics 18(1), 1988, pp. 183-190.\cr
owa <- function(x, w=rep(1/length(x), length(x))) {
   .Call("owa", x, w, PACKAGE="agop")
}


#' @rdname owa
#' @export
wam <- function(x, w=rep(1/length(x), length(x))) {
   .Call("wam", x, w, PACKAGE="agop")
}


#' @title
#' WMax, WMin, OWMax, and OWMin Operators
#' 
#' @description
#' Computes the (Ordered) Weighted Maximum/Minimum.
#' 
#' @details
#' The OWMax operator is given by
#' \deqn{
#' \mathsf{OWMax}_\mathtt{w}(\mathtt{x})=\bigvee_{i=1}^{n} w_{i}\wedge x_{\{i\}}
#' }{
#' OWMax_w(x) = max_i{ min{w_i, x_{i}} }
#' }
#' where \eqn{x_{\{i\}}}{x_{i}} denotes the \eqn{i}-th greatest
#' value in \code{x}.
#' 
#' The OWMin operator is given by
#' \deqn{
#' \mathsf{OWMin}_\mathtt{w}(\mathtt{x})=\bigwedge_{i=1}^{n} w_{i}\vee x_{\{i\}}
#' }{
#' OWMin_w(x) = min_i{ max{w_i, x_{i}} }
#' }
#' 
#' The WMax operator is given by
#' \deqn{
#' \mathsf{WMax}_\mathtt{w}(\mathtt{x})=\bigvee_{i=1}^{n} w_{i}\wedge x_{i}
#' }{
#' WMax_w(x) = max_i{ min{w_i, x_i} }
#' }
#' 
#' The WMin operator is given by
#' \deqn{
#' \mathsf{WMin}_\mathtt{w}(\mathtt{x})=\bigwedge_{i=1}^{n} w_{i}\vee x_{i}
#' }{
#' WMin_w(x) = min_i{ max{w_i, x_i} }
#' }
#' 
#' \code{OWMax} and \code{WMax} return the greatest value in \code{x}
#' by default, and \code{OWMin} and \code{WMin} - the smallest value in \code{x}.
#' 
#' Note that e.g. in the case of OWMax operator
#' the aggregation w.r.t. \code{w} gives the same result as
#' that of w.r.t. \code{sort(w)}.
#' Moreover, classically, it is assumed that if we agregate
#' vectors with elements in \eqn{[a,b]}, then
#' the largest weight should be equal to \eqn{b}.
#' 
#' There is a strong connection between the OWMax/OWMin operators
#' and the Sugeno integrals. Additionally, it may be shown
#' that the OWMax and OWMin classes are equivalent.
#' 
#' Moreover, \code{\link{index_h}} for integer data
#' is a particular OWMax operator.
#' 
#' @param x numeric vector to be aggregated
#' @param w numeric vector of the same length as \code{x}; weights
#' @return single numeric value
#' 
#' @rdname owmax
#' @name owmax
#' @export
#' @family aggregation_operators
#' @references
#' Dubois D., Prade H., Testemale C., Weighted fuzzy pattern matching, Fuzzy Sets and Systems 28, 1988, pp. 313-331.\cr
owmax <- function(x, w=rep(Inf, length(x))) {
   .Call("owmax", x, w, PACKAGE="agop")
}


#' @rdname owmax
#' @export
owmin <- function(x, w=rep(-Inf, length(x))) {
   .Call("owmin", x, w, PACKAGE="agop")
}


#' @rdname owmax
#' @export
wmax <- function(x, w=rep(Inf, length(x))) {
   .Call("wmax", x, w, PACKAGE="agop")
}


#' @rdname owmax
#' @export
wmin <- function(x, w=rep(-Inf, length(x))) {
   .Call("wmin", x, w, PACKAGE="agop")
}



# #' Computes the S-statistic2 w.r.t. to the identity function for data transformed by the inverse of a control function.
# #'
# #' Given a sequence of \eqn{n} non-negative numbers \eqn{x=(x_1,\dots,x_n)},
# #' where \eqn{x_i\ge x_j} for \eqn{i\le j},
# #' and a nondecreasing function \eqn{\gamma: R\to[0,1]}{\gamma: R->[0,1]},
# #' the \dfn{S-statistic2} (Gagolewski, Grzegorzewski, 2010) for \eqn{x} is defined as
# #' \deqn{V_n(x)=\max_{i=1,\dots,n}\{\min\{\gamma(x_i), i/n \}\}}{V_n(x)=max{ min{i/n, \gamma(x_i)} } for i=1,\dots,n}
# #'
# #' If \code{disable.check} is set to \code{FALSE}, then
# #' eventual \code{NA} values are removed from the input vector.
# #'
# #' If a non-increasingly sorted vector is given as input (set \code{sorted.dec} to \code{TRUE})
# #' the result is computed in linear time.
# #'
# #' @references Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties,
# #' In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical Methods in Data Analysis, (SMPS 2010), Springer-Verlag, 2010, 281-288.
# #'
# #' @title S-statistic2
# #' @param x a vector of real numbers.
# #' @param kappaInv a nondecreasing function ranging on [0,1], \eqn{\gamma} (see Details), the inverse of a so-called control function.
# #' @param sorted.dec logical; \code{TRUE} if the vector has already been sorted non-increasingly; defaults \code{FALSE}.
# #' @param disable.check logical; \code{TRUE} to disable some validity checks on the input vector; defaults \code{FALSE}.
# #' @return The function returns a single number or NA if improper input has been given.
# #' @seealso \code{\link{index.h}}, \code{\link{index.g}}, \code{\link{index.rp}}, \code{\link{index.lp}}, \code{\link{Sstat}}, \code{\link{psstat}}, \code{\link{dsstat}}
# #' @examples
# #' x <- rpareto2(25, 1.05, 1);
# #' kappaInv <- function(x) { pmax(0,pmin(1,x/25)); }
# #' Sstat2(x, kappaInv, FALSE, TRUE);
# #' @export
# Sstat2 <- function(x, kappaInv, sorted.dec=FALSE, disable.check=FALSE)
# {
# 	if (!disable.check)
# 	{
# 		if (length(x) == 0) return(0);
# 		if (mode(x) != "numeric") return(NA);
# 		if (any(x < 0)) return(NA);
# 		x <- x[!is.na(x)];
# 	}
# 
# 	if (!sorted.dec)
# 		x <- sort(x, decreasing=TRUE);
# 
# 	# internal method needs O(log n) time, however kappaInv(x) is O(n)
# 	.C("Sstat2", as.double(kappaInv(x)), as.integer(length(x)), out=double(1), DUP=FALSE, PACKAGE="agop")$out;
# }


# #' Computes the S-statistic w.r.t. to a control function.
# #'
# #' Given a sequence of \eqn{n} non-negative numbers \eqn{x=(x_1,\dots,x_n)},
# #' where \eqn{x_i\ge x_j} for \eqn{i\le j},
# #' and an increasing function \eqn{\kappa: [0,1]\to[a,b]}{\kappa: [0,1]->[a,b]} for some \eqn{a,b},
# #' the \dfn{S-statistic} (Gagolewski, Grzegorzewski, 2010) w.r.t. \eqn{\kappa} for \eqn{x} is defined as
# #' \deqn{V_n(x)=\max_{i=1,\dots,n}\{\min\{x_i, \kappa(i/n) \}\}}{V_n(x)=max{ min{\kappa(i/n), x_i} } for i=1,\dots,n}
# #'
# #' If \code{disable.check} is set to \code{FALSE}, then
# #' eventual \code{NA} values are removed from the input vector.
# #'
# #' If a non-increasingly sorted vector is given as input (set \code{sorted.dec} to \code{TRUE})
# #' the result is computed in linear time.
# #'
# #' @references Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties,
# #' In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical Methods in Data Analysis, (SMPS 2010), Springer-Verlag, 2010, 281-288.
# #'
# #' @title S-statistic
# #' @param x a vector of real numbers.
# #' @param kappa an increasing function, \eqn{\kappa} (see Details), a so-called control function.
# #' @param sorted.dec logical; \code{TRUE} if the vector has already been sorted non-increasingly; defaults \code{FALSE}.
# #' @param disable.check logical; \code{TRUE} to disable some validity checks on the input vector; defaults \code{FALSE}.
# #' @return The function returns a single number or NA if improper input has been given.
# #' @examples
# #' x <- rpareto2(25, 1.05, 1);
# #' kappa <- function(x) { pmax(0,pmin(1,x))*25; }
# #' Sstat(x, kappa, FALSE, TRUE);
# #' @seealso \code{\link{index.h}}, \code{\link{index.g}}, \code{\link{index.rp}}, \code{\link{index.lp}}, \code{\link{Sstat2}}, \code{\link{psstat}}, \code{\link{dsstat}}
# #' @export
# Sstat <- function(x, kappa, sorted.dec=FALSE, disable.check=FALSE)
# {
# 	if (!disable.check)
# 	{
# 		if (length(x) == 0) return(0);
# 		if (mode(x) != "numeric") return(NA);
# 		if (any(x < 0)) return(NA);
# 		x <- x[!is.na(x)];
# 	}
# 
# 	n <- length(x);
# 
# 	if (!sorted.dec)
# 		x <- sort(x, decreasing=TRUE);
# 
# 	return (max(pmin(x,kappa((1:n)/n))));
# }


