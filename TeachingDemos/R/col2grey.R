col2grey <- function(cols){
    rgb <- col2rgb(cols)
    gry <- rbind( c(0.3, 0.59, 0.11) ) %*% rgb
    rgb(gry,gry,gry, maxColorValue=255)
}

col2gray <- function(cols){
    col2grey(cols)
}
