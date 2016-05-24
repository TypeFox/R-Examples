#' cymruservices is an R package that provides interfaces to various
#' \href{http://www.team-cymru.org/services.html}{Team Cymru Services} including
#' The Bogon Refrerence, The IP to ASN Mapping Project and The Malware Hash Registry
#'
#' @name cymruservices
#' @note A direct connection to TCP Port 43 (WHOIS) is required for most of these
#'       API functions to work properly.
#' @docType package
#' @author Bob Rudis (@@hrbrmstr)
#' @import utils
#' @importFrom purrr safely
#' @importFrom stringr str_match_all
NULL
