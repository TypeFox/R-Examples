labelCANSIM <- function(df){
  for(i in 1:ncol(df)){
    label(df[[i]]) <- names(df[i])
    if(i > 2) names(df)[i] <- paste0("V", i-2)

    label(df[[1]]) <- "period"
    label(df[[2]]) <- "id"
  }

  return(df)
}
