example2 <- function(t, y, parameters){
  dy <- y*(1 - y)*(2 - y)
  list(dy)
}