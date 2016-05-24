## This file is part of the 'agop' library.
##
## Copyright 2013 Marek Gagolewski, Anna Cena
##
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



#' @title Aggregation Operators Package for R
#'
#' @description
#' ``The process of combining several numerical values into a single
#' representative one is called aggregation, and the numerical function
#' performing this process is called aggregation function.
#' This simple definition demonstrates the size of the field of application of aggregation:
#'    applied mathematics (e.g. probability, statistics, decision theory), computer science
#'    (e.g. artificial  intelligence, operation research), as well as many applied fields
#'    (economics and finance, pattern recognition and image processing, data fusion,
#'     multicriteria decision making, automated reasoning etc.). Although history of aggregation is probably
#'     as old as mathematics (think of the arithmetic mean), its existence
#'     has reminded underground till only recent (...).''
#'     (Grabisch et al, 2009, p. xiii)
#'
#' @details
#' \pkg{agop} is an open source (LGPL 3) package for R,
#' to which anyone can contribute.
#' It started as a fork of the \pkg{CITAN}
#' package (Gagolewski, 2011).
#'
#'
#' For more infrmation refer to the Package Vignette.
#' Its most recent version is available at
#' \url{http://github.com/Rexamine/agop/raw/master/inst/doc/agop-Tutorial.pdf}.
#'
#' @author
#' Marek Gagolewski \email{gagolews@@rexamine.com} [aut,cre],\cr
#' Anna Cena \email{cena@@rexamine.com} [ctb]
#'
#' \bold{Keywords}: aggregation, bibliometrics, scientometrics, scientific impact,
#' webometrics, preorders, means, OWA, OWMax, OWMin, Hirsch's h-index,
#' Egghe's g-index.
#'
#' \bold{Acknowledgments}:
#'  The development of the package in March-June 2013 was partially supported
#'  by the European Union from resources of the European Social Fund, Project PO KL
#'  ``Information technologies: Research and their interdisciplinary
#'  applications'', agreement UDA-POKL.04.01.01-00-051/10-00.
#'
#' @useDynLib agop
#' @name agop-package
#' @docType package
#' @import igraph
#' @import Matrix
#' @references
#' Beliakov G., Pradera A., Calvo T., Aggregation Functions: A Guide for Practitioners, Springer-Verlag, 2007.\cr
#' Cena A., Gagolewski M., OM3: ordered maxitive, minitive, and modular aggregation operators
#'  - Part I: Axiomatic analysis under arity-dependence, In: Bustince H. et al (Eds.),
#'  Aggregation Functions in Theory and in Practise (AISC 228), Springer-Verlag, Heidelberg, 2013, pp. 93-103. \cr
#' Cena A., Gagolewski M., OM3: ordered maxitive, minitive, and modular aggregation operators
#' - Part II: A simulation study, In: Bustince H. et al (Eds.),
#' Aggregation Functions in Theory and in Practise (AISC 228), Springer-Verlag, Heidelberg, 2013, pp. 105-115.\cr
#' Dubois D., Prade H., Testemale C., Weighted fuzzy pattern matching,
#' Fuzzy Sets and Systems 28, 1988, pp. 313-331.\cr
#' Gagolewski M., On the Relationship Between Symmetric Maxitive, Minitive,
#' and Modular Aggregation Operators, Information Sciences 221, 2013, pp. 170-180. \cr
#' Gagolewski M., Grzegorzewski P., Possibilistic Analysis of Arity-Monotonic
#' Aggregation Operators and Its Relation to Bibliometric Impact Assessment of Individuals,
#' International Journal of Approximate Reasoning 52(9), 2011, pp. 1312-1324. \cr
#' Gagolewski M., Mesiar R., Aggregating Different Paper Quality Measures
#' with a Generalized h-index, Journal of Informetrics 6(4), 2012, pp. 566-579. \cr
#' Gagolewski M., Bibliometric Impact Assessment with R and the CITAN Package,
#' Journal of Informetrics 5(4), 2011, pp. 678-692. \cr
#' Gagolewski M., Grzegorzewski P., A Geometric Approach to the Construction
#' of Scientific Impact Indices, Scientometrics 81(3), 2009, pp. 617-634. \cr
#' Gagolewski M., Statistical Hypothesis Test for the Difference between
#' Hirsch Indices of Two Pareto-Distributed Random Samples,
#' In: Kruse R. et al (Eds.), Synergies of Soft Computing and Statistics
#' for Intelligent Data Analysis (AISC 190), Springer-Verlag, Heidelberg, 2013, pp. 359-367.\cr
#' Gagolewski M., On the Relation Between Effort-Dominating and Symmetric
#'  Minitive Aggregation Operators, In: Greco S. et al (Eds.),
#'  Advances in Computational Intelligence, Part III (CCIS 299),
#'  Springer-Verlag, Heidelberg, 2012, pp. 276-285.\cr
#' Gagolewski M., Grzegorzewski P., Axiomatic Characterizations
#' of (quasi-) L-statistics and S-statistics and the Producer Assessment
#'  Problem, In: Galichet S., Montero J., Mauris G. (Eds.), Proc. EUSFLAT/LFA 2011,
#'   Atlantic Press, 2011, pp. 53-58.\cr
#' Gagolewski M., Grzegorzewski P., S-Statistics and Their Basic Properties,
#' In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical
#' Methods in Data Analysis (AISC 77), Springer-Verlag,
#' Heidelberg, 2010, pp. 281-288.\cr
#' Gagolewski M., Grzegorzewski P., Arity-Monotonic Extended Aggregation Operators,
#'  In: Hullermeier E., Kruse R., Hoffmann F. (Eds.),
#'  Information Processing and Management of Uncertainty in Knowledge-Based
#'  Systems (CCIS 80), Springer-Verlag, Heidelberg, 2010, pp. 693-702.\cr
#' Grabisch M., Marichal J.-L.,  Mesiar R., Pap E., Aggregation functions,
#' Cambridge University Press, 2009.\cr
#' Hirsch J.E., An index to quantify individual's scientific research output,
#'  Proceedings of the National Academy of Sciences 102(46), 2005, pp. 16569-16572.\cr
#' Shilkret, N., Maxitive measure and integration, Indag. Math. 33, 1971, pp. 109-116.\cr
#' Yager R.R., On ordered weighted averaging aggregation operators
#' in multicriteria decision making, IEEE Transactions on Systems,
#' Man, and Cybernetics 18(1), 1988, pp. 183-190.\cr
invisible(NULL)
