setalpha <- function(Col, alpha) {

  if (is.null(Col))
    return(Col)

  ii <- which (! is.na(Col))
  if (length(ii) > 0) {   
    pcol <- alpha.col (Col, alpha)
    Col[ii] <- pcol[ii] 
  }    
  Col
}


## =============================================================================
## Suitable colors
## =============================================================================

jet.col <- function (n = 100, alpha = 1) {

 # red-green-blue colors on scale of 0 to 1
  red    <- c(0,   0,   0,   255, 255, 128)
  green  <- c(0,   0,   255, 255, 0,   0  )
  blue   <- c(143, 255, 255, 0,   0,   0  )

  x.from <- c(0.0, seq(0.125, 1, by = 0.25), 1)  # scale from 0-1
  x.to   <- seq (0, 1, length.out = n)

  expand <- function(col)
    approx(x = x.from, y = col, xout = x.to)$y

  return (rgb(expand(red), expand(green), expand(blue), 
               maxColorValue = 255, alpha = alpha*255))
}

jet2.col <- function (n = 100, alpha = 1) {

  red    <- c(0,   0,   255, 255, 210)
  green  <- c(78,  255, 255, 0,   0  )
  blue   <- c(255, 255, 0,   0,   0  )

  x.from <- seq (0, 1, length.out = 5)  
  x.to   <- seq (0, 1, length.out = n)

  expand <- function(col)
    approx(x = x.from, y = col, xout = x.to)$y

  return (rgb(expand(red), expand(green), expand(blue), 
               maxColorValue = 255, alpha = alpha*255))
}

alpha.col <- function (col = "grey", alpha = 0.5) {

  RGBini <- col2rgb(col)

  return( rgb(t(RGBini), 
          maxColorValue = 255, alpha = alpha*255))
}

ramp.col <- function (col = c("grey", "black"), 
                      n = 100, alpha = 1) {

  RGBini <- col2rgb(col)

  x.from <- seq(0, 1, length.out = length(col)) 
  x.to   <- seq(0, 1, length.out = n)

  expand <- function(col)
    approx(x = x.from, y = col, xout = x.to)$y

  return (rgb(expand(RGBini["red",]), 
                 expand(RGBini["green",]), 
                 expand(RGBini["blue",]), 
                 maxColorValue = 255, alpha = alpha*255))
}

gg.col <-  function (n = 100, alpha = 1) {
 
  ramp.col(col = c("#0072B2", "#56B4E9", "#009E73", 
        "#CC79A7", "#D55E00", "#000000"), n = n, alpha = alpha)
} 

gg2.col <-  function (n = 100, alpha = 1) {
 
  ramp.col(col = c("#0072B2", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", 
               "#D55E00", "#000000"), n = n, alpha = alpha)
} 

MeanColors <- function(col) {           
  rgb(t(rowMeans(col2rgb(col))), maxColorValue = 255)
}

