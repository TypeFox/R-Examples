#' GLSZM Features
#'
#' @param glszm A matrix of class "glszm" produced by \code{glszm}.
#' @references \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0102107#s5} 
#' @name glszm_features
NULL
#> NULL

#' @describeIn glszm_features Small Area Emphasis
#' 
#
glszm_SAE <- function(glszm){
  dat <- t(glszm) / as.numeric(colnames(glszm))^2
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum(dat) / sum(glszm)
}

#' @describeIn glszm_features Large Area Emphasis
#' 
#
glszm_LAE <- function(glszm){
  sum(t(glszm) * as.numeric(colnames(glszm))^2) / sum(glszm)
}

#' @describeIn glszm_features Intensity Variability
#' 
#
glszm_IV <- function(glszm){
  sum(apply(glszm, 1, sum)^2)/sum(glszm)
}

#' @describeIn glszm_features Size Zone Variability
#' 
#
glszm_SZV <- function(glszm){
  sum(apply(glszm, 2, sum)^2)/sum(glszm)
}

#' @describeIn glszm_features Zone percentage
#' 
#
glszm_ZP <- function(glszm){
  n_voxels <- sum(sapply(1:ncol(glszm), function(i) sum(glszm[,i])* as.numeric(colnames(glszm)[i])))
  sum(apply(glszm/n_voxels, 2, sum))
}

#' @describeIn glszm_features Low intensity emphasis
#' 
#
glszm_LIE <- function(glszm){
  dat <- glszm / as.numeric(rownames(glszm))^2
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum(dat) / sum(glszm)
}

#' @describeIn glszm_features High intensity emphasis
#' 
#
glszm_HIE <- function(glszm){
  dat <- glszm * as.numeric(rownames(glszm))^2
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum(dat) / sum(glszm)
}

#' @describeIn glszm_features Low intensity small area emphasis
#' 
#
glszm_LISAE <- function(glszm){
  dat <- t( t(glszm) / as.numeric(colnames(glszm))^2 ) / (as.numeric(rownames(glszm))^2)
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum( dat ) / sum(glszm)
}

#' @describeIn glszm_features High intensity small area emphasis
#' 
#
glszm_HISAE <- function(glszm){
  dat <- t(glszm * as.numeric(rownames(glszm))^2) / as.numeric(colnames(glszm))^2
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum(dat) /sum(glszm)
}

#' @describeIn glszm_features Low intensity large area emphasis
#' 
#
glszm_LILAE <- function(glszm){
  dat <- t(as.numeric(colnames(glszm))^2 * t(glszm)) / as.numeric(rownames(glszm))^2
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum(dat) /sum(glszm)
}

#' @describeIn glszm_features High intensity Large area emphasis
#' 
#
glszm_HILAE <- function(glszm){
  dat <- t(as.numeric(colnames(glszm))^2 * t(glszm)) * (as.numeric(rownames(glszm))^2)
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum(dat /sum(glszm))
}