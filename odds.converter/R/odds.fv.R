#' Calculate the fair ( vigorish free) odds for a vector of vigged odds. 
#'
#' @param ... Odds of the format defined in input, the odds can be submitted as a data.frame or as a vector
#' @param input Default to US Odds. Valid selections are  "us" for US odds,
#' "prob" for Probability ,"dec" for Decimal odds,"malay" for Malaysian odds,
#' "hk" for Hong Kong odds, "indo" for Indonesian odds
#' @param output Default to US Odds. Valid selections are  "us" for US odds,
#' "prob" for Probability ,"dec" for Decimal odds,"malay" for Malaysian odds,
#' "hk" for Hong Kong odds, "indo" for Indonesian odds 
#'
#' @return A vector or a matrix of the fair values of the the given vigged values in the format chosen in output
#'
#' @examples
#' odds.fv(-120,110)
#' odds.fv(2.05,3.12,2.90,input="dec",output="prob")
#' df <- data.frame(Home = c(1.5,1.8,1.9),
#'                 Away = c(2.9,2.2,2.05))
#' odds.fv(df,input = "dec",output = "prob")
#' odds.fv(df$Home,df$Away,input="dec",output="prob")
odds.fv <- function (...,input="us",output="us") {
  
  if ( !is.element(input,c("us","prob","dec","malay","hk","indo")) ||
       !is.element(output,c("us","prob","dec","malay","hk","indo"))){
    stop("Valid Oddsformats are us,prob,dec,malay,hk,indo")}
  
  
  odds <- data.frame(...)
  
  if(ncol(odds)==1){
    odds <- t(odds)
  }
  
  if(input=="dec"){
    vig.prob.odds <-  mapply(odds.dec2prob, odds )
  }
  if (input=="indo"){
    vig.prob.odds <-  mapply(odds.indo2prob, odds )
  }
  if (input=="malay"){
    vig.prob.odds <-  mapply(odds.malay2prob, odds )
  }
  if (input=="hk"){
    vig.prob.odds <-  mapply(odds.hk2prob, odds )
  }
  if (input=="us"){
    vig.prob.odds <-  mapply(odds.us2prob, odds )
  }
  if (input=="prob"){
    vig.prob.odds <-  odds 
  }
  if(inherits(vig.prob.odds,"matrix")){
    vig <- rowSums(vig.prob.odds)
    vig[vig < 1] <- NA
  }else{
    vig <- sum(vig.prob.odds) 
    if(is.na(vig) | vig < 1 ){
      return(NA)
    }
  }
  
  fair.prob.odds <- vig.prob.odds/vig
  
  if (output=="prob"){
    return(fair.prob.odds)
  }
  
  if (output=="dec"){
    fair.prob.dec <- odds.prob2dec(fair.prob.odds)
    return(fair.prob.dec)
  }
  if (output=="indo"){
    fair.prob.indo <- odds.prob2indo(fair.prob.odds)
    return(fair.prob.indo)
  }
  if (output=="hk"){
    fair.prob.hk <- odds.prob2hk(fair.prob.odds)
    return(fair.prob.hk)
  }
  if (output=="prob"){
    return(fair.prob.odds)
  }
  if (output=="us"){
    fair.prob.us <- odds.prob2us(fair.prob.odds)
    return(fair.prob.us)
  }
  if (output=="malay"){
    fair.prob.malay <- odds.prob2malay(fair.prob.odds)
    return(fair.prob.malay)
  }
}


