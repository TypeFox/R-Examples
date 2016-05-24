#' Miscellaneous utils for Italian fiscal codes
#'
#' These package provide tools to deal with fiscal codes: a function
#' to check if the fiscal code is valid (\code{\link{wrong_fc}}) and
#' another to guess the code given personal data (\code{\link{guess_fc}}).
#'
#' Fiscal codes are far from perfect but ubiquitous personal id codes in
#' Italy, especially useful for merge purposes from a data analyst
#' standpoint.
#' 
#' \bold{WARNING}: provided routines aim to be reasonably good. Unfortunately
#' can't be perfect since, as wikipedia puts it, "On the internet, there
#' are several web sites and software packages that allow one to calculate
#' one's tax code. However, the only official tax code is the one provided
#' by the tax office, which handles cases of identical tax codes (which is
#' a pretty frequent case for people not born in Italy, as in this case the
#' 4-characters town code in the codice fiscale is replace by a 3-digit
#' country code), an arbitrary changing of a tax code, as well as cases
#' where a code is incorrect, but still valid (because provided by the tax
#' office)."
#'
#' For more info see
#' \url{http://en.wikipedia.org/wiki/Italian_fiscal_code_card}.
#'
#' @name ifctools
#' @docType package
#' @useDynLib ifctools
#' @importFrom stats var
NULL
