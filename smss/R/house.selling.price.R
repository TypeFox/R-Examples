##' @name house.selling.price
##' @title House Selling Price Data
##'
##' @description Selling price of homes in Gainesville, Florida, fall 2006, from Alachua County public records. 
##' Excerpt in Table 9.4 of 4th edition.
##' 
##' @docType data
##' @format \Sexpr[stage=build,results=rd]{data(house.selling.price); smss:::describe_df(house.selling.price)}
##' \describe{
##' \item{\code{case}}{observation id}
##" \item{\code{Taxes}}{annual tax bill (dollars)}
##" \item{\code{Beds}}{number of bedrooms}
##' \item{\code{Baths}}{number of bathrooms}
##' \item{\code{New}}{whether new (1 = yes, 0 = no)}
##' \item{\code{Price}}{selling price (dollars)}
##' \item{\code{Size}}{size of home (square feet)}
##' }
##'
##' @examples
##' data(house.selling.price)
##' summary(house.selling.price)
##' @source \url{http://www.stat.ufl.edu/~aa/social/data.html}
NULL
