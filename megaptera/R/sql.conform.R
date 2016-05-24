sql.conform <- function(string){
  string <- tolower(unlist(string)[1])
  string <- gsub("[[:punct:]]|[[:space:]]", "_", string)
  if ( length(grep("[[:digit:]]", substr(string, 1, 1))) == 1 )
    string <- paste("_", string, sep = "")
  string
}