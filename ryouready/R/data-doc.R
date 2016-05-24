# package data documentation

#### Kählers NGO Datensatz ####

#' NGO Dataset
#'
#' Data set used by Kähler (2008).
#'
#' @name d.ngo
#' @docType data
#' @references  
#' Kähler, W.-M. (2008). \emph{Statistische Datenanalyse: Verfahren verstehen
#' und mit SPSS gekonnt einsetzen}. Wiesbaden: Vieweg.
#' @keywords data
#'
NULL


#### Daten für Eta Beispiele ####

#' Sample data set for eta function examples
#'
#' Data set for eta examples.
#'
#' @name d.eta
#' @docType data
#' @keywords data
#'
NULL


#### Daten für t-Test und Wilcoxon-Test ####

#' Student self assessment data
#'
#' The participipoants were asked to assess in h
#'
#' @name d.superiority
#' @docType data
#' @keywords data
#'
NULL


#### data for examples in eta coefficient ####

# for demonstration of eta coefficient
# x   interval data
# x1  group factor
# y   dependent variable
#
# set.seed(5)
# x <- sample(10:100, 50, rep=T) / 10 
# y <- -(x-5.5)^2 + rnorm(length(x), sd=3) 
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# y <- round(range01(y) *9 + 1, 1)
# x1 <- rep(2, length(x))
# x1[x < 4 | x > 7] <- 1
# d <- data.frame(x, x1, y)
# write.table(d, "data/eta.txt", row.names=FALSE)
# write.csv2(d, "data/eta.csv", row.names=FALSE)
# plot(x,y)
# plot(x1,y)

