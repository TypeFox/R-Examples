#' GLRLM Features
#'
#' @param glrlm A matrix of class "glrlm" produced by \code{glrlm}.
#' @references \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0102107#s5} 
#' @name glrlm_features
NULL
#> NULL



#' @describeIn glrlm_features Grey level non-uniformity
#' 
glrlm_GLN <- function(glrlm){
  sum(apply(glrlm, 1, sum)^2)/sum(glrlm)
}

#' @describeIn glrlm_features High Gray level run emphasis
#' 
glrlm_HGLRE <- function(glrlm){
  sum(as.numeric(rownames(glrlm))^2 * glrlm)/sum(glrlm)
}


#' @describeIn glrlm_features Long Run Emphasis
#' 
##
glrlm_LRE <- function(glrlm){
  sum(as.numeric(colnames(glrlm))^2 * t(glrlm))/sum(glrlm)
}

#' @describeIn glrlm_features Long run high gray level emphasis
#' 
##
glrlm_LRHGLE <- function(glrlm){
  sum(t(as.numeric(colnames(glrlm))^2 * t(glrlm)) * as.numeric(rownames(glrlm))^2) /sum(glrlm)
}

#' @describeIn glrlm_features Long Run Low Gray Level Emphasis
#' 
##
glrlm_LRLGLE <- function(glrlm){
  dat <- t(as.numeric(colnames(glrlm))^2 * t(glrlm)) / as.numeric(rownames(glrlm))^2
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum(dat) /sum(glrlm)
}

#' @describeIn glrlm_features Low gray level run emphasis
#' 
## 
glrlm_LGLRE <- function(glrlm){
  dat <- glrlm / as.numeric(rownames(glrlm))^2
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum(dat) / sum(glrlm)
}

#' @describeIn glrlm_features un length non-uniformity
#' 
##R
glrlm_RLN <- function(glrlm){
  sum(apply(glrlm, 2, sum)^2)/sum(glrlm)
}

#' @describeIn glrlm_features Run Percentage
#' 
## 
glrlm_RP <- function(glrlm){
  n_voxels <- sum(sapply(1:ncol(glrlm), function(i) sum(glrlm[,i])* as.numeric(colnames(glrlm)[i])))
  sum(apply(glrlm/n_voxels, 2, sum))
}

#' @describeIn glrlm_features Short run emphasis
#' 
## 
glrlm_SRE <- function(glrlm){
  dat <- t(glrlm) / as.numeric(colnames(glrlm))^2
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum(dat) / sum(glrlm)
}


#' @describeIn glrlm_features rt run high gray level emphasis
#' 
# Sho
glrlm_SRHGLE <- function(glrlm){
  dat <- t(glrlm * as.numeric(rownames(glrlm))^2) / as.numeric(colnames(glrlm))^2
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum(dat) /sum(glrlm)
}

#' @describeIn glrlm_features Short run low grey emphasis
#' 
# 
glrlm_SRLGLE <- function(glrlm){
  dat <- t( t(glrlm) / as.numeric(colnames(glrlm))^2 ) / (as.numeric(rownames(glrlm))^2)
  dat[is.infinite(dat) | is.na(dat)] <- 0
  sum( dat ) / sum(glrlm)
}





