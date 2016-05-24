###'
###' A text progress bar with label
###'
###' This is the base txtProgressBar but with a little modification to
###' implement the label parameter for style=3. For full info see txtProgressBar
###' 
###' @param min min value for bar
###' @param max max value for bar
###' @param initial initial value for bar
###' @param char the character (or character string) to form the progress bar.
###' @param width progress bar width
###' @param label text to put at the end of the bar
###' @param title ignored
###' @param style bar style
###' @export


txtProgressBar2 <-
  function (min = 0, max = 1, initial = 0, char = "=", width = NA, 
            title="", label="", style = 1) 
{
  force(label);force(title)
  if (!style %in% 1L:3L) 
    style <- 1
  .val <- initial
  .killed <- FALSE
  .nb <- 0L
  .pc <- 0L
  nw <- nchar(char, "w")
  if (is.na(width)) {
    width <- getOption("width")
    if (style == 3L) 
      width <- width - 10L
    width <- trunc(width/nw)
  }
  up1 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    if (.nb < nb) {
      cat(paste(rep.int(char, nb - .nb), collapse = ""))
      flush.console()
    }
    else if (.nb > nb) {
      cat(paste(c("\r", rep.int(" ", .nb * nw)), collapse = ""))
      cat(paste(c("\r", rep.int(char, nb)), collapse = ""))
      flush.console()
    }
    .nb <<- nb
  }
  up2 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    if (.nb <= nb) {
      cat(paste("\r", rep.int(char, nb), collapse = ""))
      flush.console()
    }
    else {
      cat(paste(c("\r", rep.int(" ", .nb * nw)), collapse = ""))
      cat(paste(c("\r", rep.int(char, nb)), collapse = ""))
      flush.console()
    }
    .nb <<- nb
  }
  up3 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    pc <- round(100 * (value - min)/(max - min))
    if (nb == .nb && pc == .pc) 
      return()
    cat(paste(c("\r  |", rep.int(" ", nw * width + 6)), collapse = ""))
    cat(paste(c("\r  |", rep.int(char, nb), rep.int(" ", 
                                                    nw * (width - nb)), sprintf("| %3d%%", pc)), collapse = ""))
    if(nchar(label))cat(" : ",label,sep="")
    flush.console()
    .nb <<- nb
    .pc <<- pc
  }
  lab <- function(label){label<<-label}
  getVal <- function() .val
  kill <- function() if (!.killed) {
    cat("\n")
    flush.console()
    .killed <<- TRUE
  }
  up <- switch(style, up1, up2, up3)
  if (initial > min) 
    up(initial)
  structure(list(getVal = getVal, up = up, kill = kill, label=lab), class = "txtProgressBar")
}

###' set the progress bar
###'
###' update a text progress bar. See help(txtProgressBar) for more info.
###'
###' @param pb text progress bar object
###' @param value new value
###' @param title ignored
###' @param label text for end of progress bar
###' @export

setTxtProgressBar2 <- 
function (pb, value, title = NULL, label = NULL) 
{
    if (!inherits(pb, "txtProgressBar")) 
        stop("'pb' is not from class \"txtProgressBar\"")
    oldval <- pb$getVal()
    if(!is.null(label)){
      pb$label(label)
    }
    pb$up(value)
    invisible(oldval)
}
