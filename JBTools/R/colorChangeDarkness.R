colorChangeDarkness <- function(
  ##title<< Change the dark/brightness of a color
  col ##<< color to change. Can be a R color name or a hexadecimal color specification
  , factor ##<< decimal factor to darken/brighten the color with. Values < 1 lead to darker
           ##   colors,  values > 1 to brighter colors.
  )
  ##description<< Function to yield darker or brighter variants of the same color
  ##              by mixing this color with white or black
  ##seealso<<
  ##\code{\link{mixcolor}}
{
  if (inherits(col, 'character')) 
    col = col2hex(col)
  
  ##TODO make loop easier
  n.cols.out <- max(c(length(col), length(factor)))
  col        <- rep(col, length.out = n.cols.out)
  factor     <- rep(factor, length.out = n.cols.out)
  cols.out   <- character(n.cols.out)
  for (i in 1:n.cols.out) {
    if (factor[i] < 1) {
      col.mix <- hex2RGB(col2hex('black'))
      alpha   <- 1 - (1- factor[i])
    }  else {
      col.mix <- hex2RGB(col2hex('white'))
      alpha   <- 1 - (factor[i] - 1)
    }
    cols.out[i] <- hex(mixcolor(alpha, col.mix, hex2RGB(col[i])))
  }
  ##value<< hexadecimal code for the new color.
  return(cols.out)
}


