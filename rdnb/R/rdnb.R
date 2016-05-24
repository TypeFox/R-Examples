#' @title R Interface to the Deutsche Nationalbibliothek (German National Library) API
#' @description A wrapper for the Deutsche Nationalbibliothek (German National Library) API, available at \url{http://www.dnb.de}. The German National Library is the German central archival library, collecting, archiving, bibliographically classifying all German and German-language publications from 1913, foreign publications about Germany, translations of German works, and the works of German-speaking emigrants published abroad between 1933 and 1945.
#' A personal access token is required for usage.
#' @name rdnb
#' @docType package
#' @details All bibliographic data of the German National Library are provided free of charge and can be freely re-used under "Creative Commons Zero" (CC0 1.0) terms. The metadata and online interfaces are provided with no guarantee of their being continuous, punctual, error-free or complete, or of their not infringing the rights of third parties (e.g. personal rights and copyright).
#' A personal access token is required for usage and can be requested by sending an e-mail to the Interface Service (\email{schnittstellen-service@@dnb.de}). The e-mail must include the required catalogue "Catalogue of German National Library (DNB) / Katalog der Deutschen Nationalbibliothek (DNB)" and the access option "via access token / ueber Zugangscode".
#' If you do not want to enter your token for each R session, put the following in your .Renviron or .Rprofile file: \code{DNB_TOKEN=PUTYOURTOKENHERE}.
#' @references About the DNB: \url{http://www.dnb.de/EN/Wir/wir_node.html}; about the interface and access requirements: \url{http://www.dnb.de/EN/Service/DigitaleDienste/SRU/sru_node.html}; the DNB web search: \url{http://dnb.dnb.de}
#' @importFrom httr GET
#' @importFrom httr content
#' @importFrom xml2 as_list
#' @importFrom xml2 read_xml
#' @importFrom stats setNames
#' @import brew
#' @import grDevices
#' @import methods
#' @import utils
#' @aliases rdnb rdnb-package
NULL
