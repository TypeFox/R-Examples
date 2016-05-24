#' Datasets from Hosmer, Lemeshow and Sturdivant, "Applied Logistic
#' Regression" (3rd ed.)
#'
#' This package is a unofficial companion to the textbook "Applied
#' Logistic Regression" by D.W. Hosmer, S. Lemeshow and
#' R.X. Sturdivant (3rd ed.).
#'
#' It includes all the datasets used in the book, both for easy reproducibility
#' (assuming a didactic point of view) and algorithms benchmarking purposes.
#'
#' Some analysis proposed in the text are reproduced in the examples,
#' in order to provide data testing and code demos at the same time.
#'
#' The vignette includes all the examples (with graphics too); therefore is
#' organized per-dataset.
#'
#' Datasets and variables have lower-case name with respect to the
#' original sources. Categorical data were packaged as factor.
#'
#' Regarding data coding, help pages lists the internal/factor representation
#' of the data (eg 1: No, 2: Yes), not the original one (eg 0: No, 1:
#' Yes). This is intended to allow easier/safer recoding based on
#' as.integer, especially for multinomial variables. 
#' 
#' @name aplore3
#' @docType package
#' @source Hosmer, D.W., Lemeshow, S. and Sturdivant, R.X. (2013) Applied
#'         Logistic Regression, 3rd ed., New York: Wiley
NULL
