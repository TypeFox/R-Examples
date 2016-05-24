slog <- function(..., file = "susili.log"){
  
  x <- c(...)
  cat(x, file = "")
  if ( file != "" )
    cat(x, file = file, append = TRUE)
}