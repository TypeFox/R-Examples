momIntegrated <- function(densFn, order, param = NULL, about = 0,
                          absolute = FALSE){
     if (missing(densFn) || !(is.function(densFn) || is.character(densFn))){
         stop("'densFn' must be supplied as a function or name")
     }
     ## Set default integration limits
     low <- -Inf
     high <- Inf
     if (is.character(densFn)) {
         distname <- tolower(densFn)
         if (is.null(densFn)){
             stop("unsupported distribution")
         }
         if (distname == "ghyp" || distname == "generalized hyperbolic") {
             if (is.null(param)) {
                 param <- c(1,1,0,1,0)
             }
             if (absolute == FALSE) {
                 ddist <- function(x, order, param, about) {
                     (x - about)^order * dghyp(x, Theta = param)
                 }
             } else {
                 ddist <- function(x, order, param, about) {
                     abs(x - about)^order * dghyp(x, Theta = param)
                 }
             }
         } else if (distname == "hyperb" || distname == "hyperbolic") {
             if (is.null(param)) {
                 param <- c(0,1,1,0)
             }
             if (absolute == FALSE) {
                 ddist <- function(x, order, param, about) {
                     (x - about)^order * dhyperb(x, Theta = param)
                 }
             } else {
                 ddist <- function(x, order, param, about) {
                     abs(x - about)^order * dhyperb(x, Theta = param)
                 }
             }
         } else  if (distname == "gig" ||
                     distname == "generalized inverse gaussian") {
             if (is.null(param)) {
                 param <- c(1,1,1)
             }
             low <- 0
             if (absolute == FALSE) {
                 ddist <- function(x, order, param, about) {
                     (x - about)^order * dgig(x, Theta = param)
                 }
             } else {
                 ddist <- function(x, order, param, about) {
                     abs(x - about)^order * dgig(x, Theta = param)
                 }
             }
         } else  if (distname == "gamma") {
             if (order <= -param[1]){
                 stop("Order must be greater than shape parameter for gamma")
             }
             if (is.null(param)) {
                 param <- c(1,1)
             }
             low <- 0
             if (absolute == FALSE) {
                 ddist <- function(x, order, param, about) {
                     (x - about)^order *
                         dgamma(x, shape = param[1], rate = param[2])
                 }
             } else {
                 ddist <- function(x, order, param, about) {
                     abs(x - about)^order *
                         dgamma(x, shape = param[1], rate = param[2])
                 }
             }
         } else  if (distname == "invgamma" ||
                     distname == "inverse gamma"){
             if (is.null(param)) {
                 param <- c(-1,1)
             }
             if (param[1] <= order){
                 stop("Order must be less than shape parameter for inverse gamma")
             }
             low <- 0
             dinvgamma <- function(x, shape, rate = 1, scale = 1/rate){
                 dens <- ifelse(x <= 0, 0,
                                (scale/x)^shape*exp(-scale/x)/(x*gamma(shape)))
                 return(dens)
             }

             if (absolute == FALSE) {
                 ddist <- function(x, order, param, about) {
                     (x - about)^order *
                         dinvgamma(x, shape = param[1], rate = param[2])
                 }
             } else {
                 ddist <- function(x, order, param, about) {
                     abs(x - about)^order *
                         dinvgamma(x, shape = param[1], rate = param[2])
                 }
             }


         } else if (distname == "vg" || distname == "variance gamma") {
             if (!exists("dvg", mode = "function")){
                 stop("package Variance Gamma needs to be loaded")
             }
             if (is.null(param)) {
                 param <- c(0,1,0,1)
             }
             if (absolute == FALSE) {
                 ddist <- function(x, order, param, about) {
                     (x - about)^order * dvg(x, param = param)
                 }
             } else {
                 ddist <- function(x, order, param, about) {
                     abs(x - about)^order * dvg(x, param = param)
                 }
             }
         }
     } else {
         if (is.function(densFn)) {
             stop("general density functions code not yet implemented")
         }
     }
     mom <- integrate(ddist, low, high, param = param,
                      order = order, about = about,
                      subdivisions = 1000,
                      rel.tol = .Machine$double.eps^0.5)[[1]]

     ## Return Value:
     return(mom)
 }
