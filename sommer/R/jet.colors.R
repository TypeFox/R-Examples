jet.colors <- function(n, alpha=1) {
  if(n > 0) {
    if(length(alpha) != 1 & length(alpha) != n) {
      print('Warning: using only first alpha value')
      alpha <- alpha[1]
    }
    if(length(alpha) == 1) {
      alpha <- rep(alpha, n)
    }
    ## TODO Include alpha values
    return(colorRampPalette(c('#000066', 'blue', 'cyan', 'yellow',
                              'red', '#660000'))(n))
  } else {
    ## Return an empty character string if they requested nothing.
    character()
  }
}

bathy.colors <- function(n, alpha=1)
  return(rgb(seq(0.9,0,len=n), seq(0.9,0,len=n), 1, alpha))
