#' @title Package for calculation of classic and Fuzzy AHP
#'
#' @description
#' \pkg{FuzzyAHP} is an open source (LGPL 3) package for R.
#' The package is only suitable for AHP that uses categorical rating of criteria
#' for alternatives instead of pairwise comparison of alternatives according to
#' each criteria. This adaptation of AHP is common in situations when the number
#' of alternatives is hight and the pairwise comparison is thus inadequate or
#' impossible to construct. The weights for criteria are, however, still determined
#' from the pairwise comparison matrix. This approach towards AHP is common in
#' Geosciences as well as other fields.
#'
#' The determination of criteria weights is done according to process described
#' by Krejčí, Pavlačka, and Talašová (2016), which yelds significantly narrower
#' fuzzy numbers than previously used approaches.
#'
#' @details
#' Please see vignettes for more details about the package and examples of use.
#'
#' Complete list of classes and methods call \code{help(package="FuzzyAHP")}.
#' \cr\cr
#'
#' @references
#' Krejčí, Jana, Ondřej Pavlačka, and Jana Talašová. 2016. “A fuzzy extension of Analytic Hierarchy Process based on the constrained fuzzy arithmetic.” Fuzzy Optimization and Decision Making. doi:10.1007/s10700-016-9241-0.
#'
#' @name FuzzyAHP-package
#' @docType package
#' @import methods
#' @author Jan Caha \email{cahik@@atlas.cz}, with contributions from Aneta Drážná
#'
#' @encoding UTF-8
invisible(NULL)
