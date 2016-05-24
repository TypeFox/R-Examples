# ciao<-read.table(file="Z:/aon.txt",header=FALSE,sep=" ")
# names(ciao)<-c("class","loss","sumInsured")
# aon<-transform(ciao, sumInsured=sumInsured*1000)
# devtools::use_data(aon)

#' AON Re Belgium data.
#'
#' A dataset containing losses and sum insured for building property damage of AON Re Belgium portfolio.
#' Claims are split by building category.
#'
#' @format A data frame with 1823 rows and 3 variables:
#' \describe{
#'   \item{class}{building category}
#'   \item{loss}{loss size}
#'   \item{sumInsured}{sum insured}
#'   ...
#' }
#' @source \url{http://lstat.kuleuven.be/Wiley/}
"aon"

# ciao<-read.table(file="Z:/lossdata.txt",header=FALSE,sep="\t")
# names(ciao)<-c("loss","alae","limit","censored")
# loss<-ciao
# devtools::use_data(loss)

#' Loss-ALAE data of Freez and Valdez
#'
#' A dataset containing losses (claim amount), alae and policy limit.
#' Also censoring information is reported (actual claim amount to exceed policy limit).
#'
#' @format A data frame with 1500 rows and 4 variables:
#' \describe{
#'   \item{loss}{actual claim amount for the claim}
#'   \item{alae}{allocated loss adjustment expense for the claim}
#'   \item{limit}{policy limit}
#'   \item{censored}{0 not censored, 1 censored}
#'   ...
#' }
#' @source \url{http://lstat.kuleuven.be/Wiley/}
"loss"
