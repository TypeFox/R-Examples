#' Derivative Titration Curves
#' 
#' Given a two-column data frame, with volumes of titrant in the first
#' column and pH, pAnalyte, pTitrant, pMetal, or solution potential in 
#' the second column, this function calculates and plots the associated
#' first-derivative and second derivative titration curves.
#' 
#' @param df A data frame with two columns; the first column must
#' contain the volumes of titrant and the second column must contain
#' values for the associated pH, pAnalyte, pTitrant, pMetal, or 
#' solution potential. A typical object to pass to this function is 
#' that created by the other functions in this package; however, the 
#' data frame can be one prepared separately, provided that it matches 
#' the structure defined above.
#' 
#' @param plot Logical; if TRUE, plots the titration curve.
#'  
#' @param \dots Additional arguments to pass to \code{plot()} function.
#' 
#' @return A list that consists of two data frames, one for the first
#' derivative titration curve and one for the second derivative 
#' titration curve.
#' 
#' @author David T. Harvey, DePauw University. \email{harvey@@depauw.edu}
#' 
#' @export
#' 
#' @importFrom graphics plot lines par
#' 
#' @examples
#' ### Derivative weak acid/strong base titration curves
#' ab = wa_sb()
#' ex16 = derivative(ab)
#' str(ex16)


derivative = function(df, plot = TRUE, ...) {
  n = length(df[ , 1])
  x = df[ , 1]
  y = df[ , 2]
  x1 = seq(1,n-1)
  x2 = seq(1,n-2)
  y1 = seq(1,n-1)
  y2 = seq(1,n-2)
  for (i in 1:(n-1)) {
    x1[i] = (x[i] + x[i+1])/2
    y1[i] = (y[i+1] - y[i])/(x[i+1] - x[i])
  }
  for (j in 1:(n-2)) {
    x2[j] = (x1[j] + x1[j+1])/2
    y2[j] = (y1[j+1] - y1[j])/(x1[j+1] - x1[j])
  }
  if (plot == TRUE) {
  opt = par(mfrow = c(1,2))
  plot(x1, y1, type = "l", col = "blue", lwd = 2, 
       xlab = "volume of titrant (mL)", ylab = "first derivative",
       xaxs = "i", ...)
  plot(x2, y2, type = "l", col = "blue", lwd = 2, 
       xlab = "volume of titrant (mL)", ylab = "second derivative", 
       xaxs = "i",  ...)
  par(opt)
  }
  df.f = data.frame(x1,y1)
  df.s = data.frame(x2,y2)
  obj.ret = list("first_deriv" =df.f,"second_deriv"=df.s)
  invisible(obj.ret)
}

