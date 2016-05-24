################################################################################
# Color schemes
################################################################################
# apply color scheme to a numeric vector
# decreasing: relationship gets stronger with decreasing values?
# TRUE for, typically, e-values, gaps, mism, ... 190
# FALSE for bit scores, per_id, aln_length
apply_color_scheme <- function(x, direction=NULL, color_scheme="grey",
                               decreasing=FALSE, rng=NULL, transparency=0.5){
  # check arguments
  # if length is 0, return length 0
  if (!is.null(direction) && length(direction) == 0) return (character(0))
  # if x is null and direction is not, get x to 1s (mainly for blue/red)
  if (is.null(x) && !is.null(direction)) {
    x <- rep(1, length(direction))
  }
  if (!is.numeric(x)) stop("Color column is not numeric")
  if (is.null(rng)) rng <- range(x)
  if (!(is.logical(transparency) && !(transparency))
      && !is.numeric(transparency))
    stop ("transparency should be FALSE or numeric")
  col <- rep(grey(0.5), length(x))
  # red blue
  if (any(color_scheme %in% c("red_blue", "blue_red"))){
    if (is.null(direction)) direction <- rep(1, length(x))
    blues <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
               "#4292C6", "#2171B5", "#08519C", "#08306B")
    reds  <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A",
               "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
    # case: only one value:
    if (diff(rng) == 0){
      level <- rep(9, length(x))
    } else { # case: several values
      level <- round(((x-rng[1])/diff(rng))*8+1)
    }
    if (decreasing) level <- -level + 10
    col[direction==1] <- reds[level[direction==1]] 
    col[direction==-1] <- blues[level[direction==-1]]
  }
  # grey: between 0.25 and 0.75
  else if (any(color_scheme %in% c("grey", "gray", "grays", "greys")))
    # case: only one value:
    if (diff(rng) == 0){
      col <- rep(grey(0.5), length(x))
    } else {
      level <- 0.75-((x-rng[1])/diff(rng))*0.5
      if (decreasing) level <- -level+1
      col <- grey(level)
    }
  else {
    stop("Color scheme name invalid, choose between red_blue or grey")
  }
  if (transparency && transparency != 1){
    # Convert ratio into hexadecimal
    if (transparency > 1 || transparency < 0)
      stop("transparency should be between 0 and 1")
    tpc <- floor(transparency*256)
    tpc <- sprintf("%X", tpc)
    if (nchar(tpc) == 1) tpc <- paste("0", tpc, sep="")
    col <- paste(col, tpc, sep="")
  }
  col
}
