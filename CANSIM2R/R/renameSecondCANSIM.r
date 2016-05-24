renameSecondCANSIM <- function(df, nColumnsInFirstCansim){

  for(i in 1:ncol(df)){
    if(i > 2) names(df)[i] <- paste0("V", nColumnsInFirstCansim+i-2-2)
  }

  return(df)
}
