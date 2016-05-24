dopos <- function(pos, limits, width, height, side.legend, shift)
{
  if(side.legend > 1L)
    limits <- rev(limits)
  shift <- c(diff(limits[[1]]), diff(limits[[2]])) * shift
  if(pos == "bottomleft") {
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + shift[1]
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + shift[2]
  }
  if(pos == "topleft") {
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + shift[1]
    ylim <- c(max(limits[[2L]], na.rm = TRUE) - height, max(limits[[2L]], na.rm = TRUE)) - shift[2]
  }
  if(pos == "topright") {
    xlim <- c(max(limits[[1L]], na.rm = TRUE) - width, max(limits[[1L]], na.rm = TRUE)) - shift[1]
    ylim <- c(max(limits[[2L]], na.rm = TRUE) - height, max(limits[[2L]], na.rm = TRUE)) - shift[2]
  }
  if(pos == "bottomright") {
    xlim <- c(max(limits[[1L]], na.rm = TRUE) - width, max(limits[[1L]], na.rm = TRUE)) - shift[1]
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + shift[2]
  }
  if(pos == "bottom") {
    m <- mean(limits[[1]] - min(limits[[1]], na.rm = TRUE), na.rm = TRUE) - 0.5 * width
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + m
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + shift[2]
  }
  if(pos == "top") {
    m <- mean(limits[[1]] - min(limits[[1]], na.rm = TRUE), na.rm = TRUE) - 0.5 * width
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + m
    ylim <- c(max(limits[[2L]], na.rm = TRUE) - height, max(limits[[2L]], na.rm = TRUE)) - shift[2]
  }
  if(pos == "left") {
    m <- mean(limits[[2]] - min(limits[[2]], na.rm = TRUE), na.rm = TRUE) - 0.5 * height
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + shift[1]
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + m
  }
  if(pos == "right") {
    m <- mean(limits[[2]] - min(limits[[2]], na.rm = TRUE), na.rm = TRUE) - 0.5 * height
    xlim <- c(max(limits[[1L]], na.rm = TRUE) - width, max(limits[[1L]], na.rm = TRUE)) - shift[1]
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + m
  }
  if(pos == "center") {
    mx <- mean(limits[[1]] - min(limits[[1]], na.rm = TRUE), na.rm = TRUE) - 0.5 * width
    my <- mean(limits[[2]] - min(limits[[2]], na.rm = TRUE), na.rm = TRUE) - 0.5 * height
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + mx
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + my
  }

  return(list(xlim = xlim, ylim = ylim))
}

