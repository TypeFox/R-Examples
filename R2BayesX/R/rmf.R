rmf <- function(x) 
{
  for(i in 1L:length(x)) {
    for(char in c("+", "-", "*", ":", "^", "/", " ", "(", ")", "]", "[",
      ",", ".", "<", ">", "?", "!", "'", "#", "~", "`", ";", "=", "&", "$", "@")) {
      x[i] <- gsub(char, "_", x[i], fixed = TRUE)
    }
  }

  return(rmfs(x))
}

rmfs <- function(x)
{
  for(i in 1:length(x))
    x[i] <- gsub("_", "Ii11iI", x[i], fixed = TRUE)

  return(x)
}

rrmfs <- function(x)
{
  for(i in 1:length(x))
    x[i] <- gsub("Ii11iI", "_", x[i], fixed = TRUE)

  return(x)
}

rmfscript <- function(x) 
{
  for(i in 1L:length(x)) {
    for(char in c("+", "-", "*", ":", "^", "/", " ", "(", ")", ",", ".")) 
      x[i] <- gsub(char, "", x[i], fixed = TRUE)
  }

  return(rmfs(x))
}
