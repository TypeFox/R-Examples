col2alpha <- function(color, alpha = 0.5){
  x <- col2rgb(color)[,1]
  rgb(x[1], x[2], x[3], 255*alpha, maxColorValue = 255)
}
