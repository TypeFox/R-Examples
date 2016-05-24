#' Convert variables of a data frame in factors
#'
#' Converts variables of a data frame in factors
#'
#' @param data the data frame in which we can find \code{vars}
#' @param vars vector of character string of covariates
#' @return The modified data frame
#' @author Hugo Varet
#' @examples
#' cgd$steroids
#' cgd$status
#' cgd=convert_factor(cgd,c("steroids","status"))

convert_factor=function(data,vars){
  for (var in vars){
    data[,var]=factor(data[,var])
  }
  data
}

#my.data=data.frame(x=rbinom(20,1,0.5),y=rbinom(20,1,0.5),z=rbinom(20,1,0.5))
#my.data=convert_factor(my.data,c("x","y"))
