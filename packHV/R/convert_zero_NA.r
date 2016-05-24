#' Convert 0s in NA
#'
#' Converts 0s in \code{NA}
#'
#' @param data the data frame in which we can find \code{vars}
#' @param vars a character vector of covariates for which to transform 0s in \code{NA}
#' @return The modified data frame
#' @author Hugo Varet
#' @examples
#' my.data=data.frame(x=rbinom(20,1,0.5),y=rbinom(20,1,0.5),z=rbinom(20,1,0.5))
#' my.data=convert_zero_NA(my.data,c("y","z"))

convert_zero_NA=function(data,vars){
  for (var in vars){
    data[,var]=ifelse(data[,var]==0,NA,data[,var])
  }
  data
}

#my.data=data.frame(x=rbinom(20,1,0.5),y=rbinom(20,1,0.5),z=rbinom(20,1,0.5))
#my.data=convert_zero_NA(my.data,c("y","z"))

