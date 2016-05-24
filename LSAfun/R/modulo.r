modulo <- function(x,y){
  
  mod <- x%%y
  
  if(mod == 0){mod1 <- y}
  if(mod != 0){mod1 <- mod}
  
  mod1
}