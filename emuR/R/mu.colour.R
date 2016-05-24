##' Function for specifying color, linetype, and line-widths in EMU plotting
##' functions.
##' 
##' The function specifies color, linetype and linewidths in EMU plotting
##' functions as is used mostly in calls from within plot.trackdata,
##' plot.spectral, eplot, and dplot
##' 
##' Parameters are also supplied for use with the function 'legend'
##' 
##' @param labs A vector of character labels
##' @param col A code passed to the 'col' argument in plotting functions. There
##' are four possibilities. Either logical, a character vector, or a numeric
##' vector.  In the first case, if TRUE, then a different numeric code is given
##' for each unique label type. For example, if labs is c("a", "b", "a", "c"),
##' then the output is c(1, 2, 1, 3). If F, then for this example, the output
##' is c(1, 1, 1, 1). In the second case, the character vector can be either a
##' single element specifying a character, or there can be as many elements as
##' there are unique colors. Thus if col = "red", then for the example c("a",
##' "b", "a", "c"), the output is c("red", "red", "red", "red"). Alternatively,
##' since there are three unique labels for this example, then the user could
##' specify col = c("green", "red", "blue") and the output is c("green", "red",
##' "green", "blue") if labs is c("a", "b", "a", "c").  In the third case,
##' 'col'. can be either a single element numeric vector, or its length must be
##' equal to the number of unique types in labs.  For example, if col=3 and if
##' labs = c("a", "b", "a", "c"), then the output is c(3, 3, 3, 3).
##' Alternatively, if col = c(2,3,1), then the output is c(2, 3, 2, 1) for the
##' same example. Finally, col can be specified as a character or numeric
##' vector that is the same length as labs, allowing the user to choose the
##' color in which each line should be drawn.  The default is col = TRUE.
##' @param linetype A code specifying linetypes, i.e. the values passed to lty
##' in plotting functions.There are 2 possibilities.  Either logical, a
##' character vector, or a numeric vector.  In the first case, if TRUE, then a
##' different numeric code is given for each unique label type. For example, if
##' labs is c("a", "b", "a", "c"), then the output is c(1, 2, 1, 3). If F, then
##' for this example, the output is c(1, 1, 1, 1).  In the second case,
##' 'linetype' can be either a single element numeric vector, or its length
##' must be equal to the number of unique types in labs.  For example, if
##' linetype=3 and if labs = c("a", "b", "a", "c"), then the output is c(3, 3,
##' 3, 3). Alternatively, if linetype = c(2,3,1), then the output is c(2, 3, 2,
##' 1) for the same example. Finally, linetype can be specified as a numeric
##' vector that is the same length as labs, allowing the user to choose the
##' linetype in which each line should be drawn.  The default is linetype=F
##' @param lwd A code passed to the lwd argument in plotting functions.  'lwd'
##' can be either a single element numeric vector, or its length must be equal
##' to the number of unique types in labs.  For example, if lwd=3 and if labs =
##' c("a", "b", "a", "c"), then the output is c(3, 3, 3, 3). Alternatively, if
##' lwd = c(2,3,1), then the output is c(2, 3, 2, 1) for the same example. The
##' default is NULL in which case all lines are drawn with lwd=1
##' @param pch A code passed to the pch argument in plotting functions.
##' Functions in the same way as lwd above
##' @return If it is a LISTRUE, use \item{colour}{A code for the color'}
##' \item{linetype}{A code for the linetype} \item{lwd}{A code for the line
##' width} \item{legend}{A list consisting of \$legend\$lab, \$legend\$lty and
##' \$legend\$lwd that specify the parameters for the 'legend' function.
##' 
##' ...  }
##' @author Steve Cassidy, modified by Jonathan Harrington
##' @seealso \code{\link{plot.trackdata}} \code{\link{dplot}}
##' \code{\link{eplot}} \code{\link{plot.spectral}}
##' @keywords utilities
##' @examples
##' 
##' # examples will be given using the above functions
##' # b/w but with different linetypes
##' eplot(vowlax.fdat.5[,1:2], vowlax.l, col=FALSE, lty=TRUE)
##' 
##' # user-defined colors
##' eplot(vowlax.fdat.5[,1:2], vowlax.l, col=c("green", "blue", "red", "orange"))
##' 
##' # spectral plot, user-defined colors, the last one is dotted
##' # and with a line-thickness of 2
##' plot(vowlax.dft.5[1:20,], vowlax.l[1:20], 
##' col=c("green", "blue", "red", "orange"), 
##' fun=mean, lty=c(1, 1, 1, 2), lwd=c(1, 1, 1, 2))
##' 
##' # similar but using dplot()
##' dplot(vowlax.fdat[1:20,2], vowlax.l, 
##' col=c("green", "blue", "red", "orange"), 
##' lwd=c(1, 1, 1, 2), lty=c(1, 1, 1, 2))
##' 
##' # the default except plot everything with a dotted line and plotting symbol 4
##' dplot(vowlax.fdat[,2], vowlax.l, average=TRUE, lty=2, pch=4, type="b", xlim=c(40, 60))
##' 
##' # the default except plot everything with a dotted line and
##' # with double line thickness
##' eplot(vowlax.fdat.5[,1:2], vowlax.l, lty=2, lwd=2)
##' 
##' @export mu.colour
`mu.colour` <- function (labs, col = TRUE, linetype = FALSE, 
                         lwd = NULL, pch=NULL) 
{
  result <- NULL
  if (is.logical(col)) {
    if (col) 
      result$colour <- label_num(labs)
    else result$colour <- rep(1, length(labs))
  }
  else if (length(col) == length(labs)) 
    result$colour <- col
  else if (length(col) == length(unique(labs))) {
    k <- 1
    result$colour <- labs
    for (j in unique(labs)) {
      temp <- labs == j
      result$colour[temp] = col[k]
      k <- k + 1
    }
  }
  else if (length(col) == 1) 
    result$colour <- rep(col, length(labs))
  if (is.logical(linetype)) {
    if (linetype) 
      result$linetype <- label_num(labs)
    else result$linetype <- rep(1, length(labs))
  }
  else if (length(linetype) == length(labs)) 
    result$linetype <- linetype
  else if (length(linetype) == length(unique(labs))) {
    k <- 1
    result$linetype <- labs
    for (j in unique(labs)) {
      temp <- labs == j
      result$linetype[temp] = linetype[k]
      k <- k + 1
    }
  }
  else if (length(linetype) == 1) 
    result$linetype <- rep(linetype, length(labs))
  if (is.null(lwd)) 
    result$lwd <- rep(1, length(labs))
  else if (length(lwd) == length(labs)) 
    result$lwd <- lwd
  else if (length(lwd) == length(unique(labs))) {
    k <- 1
    result$lwd <- labs
    for (j in unique(labs)) {
      temp <- labs == j
      result$lwd[temp] = lwd[k]
      k <- k + 1
    }
  }
  else if (length(lwd) == 1) 
    result$lwd <- rep(lwd, length(labs))
  
  if (is.null(pch)) 
    result$pch <- rep(1, length(labs))
  else if (length(pch) == length(labs)) 
    result$pch <- pch
  else if (length(pch) == length(unique(labs))) {
    k <- 1
    result$pch <- labs
    for (j in unique(labs)) {
      temp <- labs == j
      result$pch[temp] = pch[k]
      k <- k + 1
    }
  }
  else if (length(pch) == 1) 
    result$pch <- rep(pch, length(labs))
  
  
  p1 <- paste(labs, result$colour, result$linetype, result$lwd, result$pch)
  p1.temp <- duplicated(p1)
  result$legend$lab <- labs[!p1.temp]
  result$legend$col <- result$colour[!p1.temp]
  result$legend$lty <- result$linetype[!p1.temp]
  result$legend$lwd <- result$lwd[!p1.temp]
  result$legend$pch <- result$pch[!p1.temp]
  result
}





## return the colour for a given label via the colour object








##' get a EMU color
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export mu.colour.get
mu.colour.get <- function(col.lty, label) {
  
  colour <- col.lty$legend$col[match(label, col.lty$legend$lab)]
  return( colour )
  
}









##' mu linetype get
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export mu.linetype.get
mu.linetype.get <- function(col.lty, label) {
  
  lty <- col.lty$legend$lty[match(label, col.lty$legend$lab)]
  return( lty )
  
}
