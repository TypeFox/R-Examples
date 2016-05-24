createStatCanVariables <- function(df){
  VectorPosition <- match("Vector",names(df))

  #Only create new variable if there is more than one column from StatCan
  if(VectorPosition > 4) df$StatCanVariable <- apply(df[,c(3:(VectorPosition-1))], 1, function(x) paste(x, collapse = "; "))
  else df$StatCanVariable <- df[,3]

  return(df)
}
