
#' @title
#' oin_group_element
#' 
#' @description
#' Calculates the organization group from the organization number for one oin.
#' 
#' @param one_oin Character elemen with oin.
#' 
#' @references 
#' http://www.skatteverket.se/download/18.70ac421612e2a997f85800040284/1302507382017/70909.pdf
#' 
#' @keywords internal
#' 
#' @return
#' Character categegory of organisational group.

oin_group_element <- function(one_oin){
  switch(substr(one_oin, 1, 1),
         "1" = "D\u00F6dsbo",
         "2" = "Stat, landsting, kommuner, f\u00F6rsamlingar",
         "3" = "Utl\u00E4ndska f\u00F6retag som bedriver n\u00E4ringsverksamhet eller \u00E4ger fastigheter i Sverige",
         "5" = "Aktiebolag",
         "6" = "Enkelt bolag",
         "7" = "Ekonomiska f\u00F6reningar",
         "8" = "Ideella f\u00F6reningar och stiftelser",
         "9" = "Handelsbolag, kommanditbolag och enkla bolag",
         as.character(NA))
}
