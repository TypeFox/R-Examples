plotAdditionalAxis <- function(
##title<< Add a second plot axis with transformed label values
    side = 1      ##<< integer: which axis to use as a basis for the second one.
    ,trans.fun    ##<< function: the transfer function to use between the two axis values.
    ,label = c()  ##<< character: labels of the axis.
    ,...          ##<< further arguments to pass to the axis call.
    )
  ##description<<
  ## This function adds a second axis additional labels to a plot. It uses the
  ## axis values of the opposite side and mathematically transforms these into the
  ## values added. This can be used for example to indicate the period and frequency
  ## of a periodic signal. 
  ##seealso<<
  ## \code{\link{axis}}
{
  ##TODO move to package JBTools
  xaxis.values <- axTicks(side = (2 - side %% 2))
  side.new     <- (side + 1) %% 4 + 1
  axis(side = side.new, at = xaxis.values, labels = round(trans.fun(xaxis.values), digits = 2))
  mtext(side = side.new, text = label, line = 3)
  ##value<<
  ## Nothing is returned.
}
