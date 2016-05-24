
#' Create a Schumaker spline
#' @export
#' @param x A vector of x coordinates
#' @param y A corresponding vector of y coordinates
#' @param ff (Optional) A corresponding vector of gradiants at the data points. If not supplied this is estimated.
#' @param Vectorised This is a boolean parameter. Set to TRUE if you want to be able to input vectors to the created spline. If you will only input single values set this to FALSE as it is a bit faster.
#' @param Extrapolation This determines how the spline function responds when an input is recieved outside the domain of x. The options are "Curve" which outputs the result of the point on the quadratic curve at the nearest interval, "Constant" which outputs the y value at the end of the x domain and "Linear" which extends the spline using the gradiant at the edge of x.
#'
#' @return A list with 3 spline functions and a table with spline intervals and coefficients. The first spline is the schumaker spline, the second spline is the first derivative of the schumaker spline, the third spline is the second derivative of the schumaker spline. Each function takes an x value (or vector if Vectorised = TRUE) and outputs the interpolated y value (or relevant derivative).
#' @references Judd (1998). Numerical Methods in Economics. MIT Press
#' @examples
#' x = seq(1,6)
#' y = log(x)

#' SSS = schumaker::Schumaker(x,y, Vectorised = TRUE)

#' xarray = seq(1,6,0.01)
#' Result = SSS$Spline(xarray)
#' Result2 = SSS$DerivativeSpline(xarray)
#' Result3 = SSS$SecondDerivativeSpline(xarray)

#' plot(xarray, Result, ylim=c(-0.5,2))
#' lines(xarray, Result2, col = 2)
#' lines(xarray, Result3, col = 3)

Schumaker <- function(x,y, ff = "Not-Supplied", Vectorised = TRUE, Extrapolation = c("Curve", "Constant", "Linear")){
  Extrapolation = Extrapolation[1]
  if (!(Extrapolation %in% c("Constant", "Linear", "Curve"))){stop("The extrapolation parameter defines what the function returns when evaluated
                                                                      outside the domain of the interpolation data. \n Choose 'Constant' for constant
                                                                      extrapolation. This returns the value at the nearest edge of the domain. \n 'Linear'
                                                                      extends out a line from the edge of the domain with a slope of the gradiant at
                                                                      that point. \n 'Curve' extrapolation uses the parabolic arc for the last interval.")}

# Schumaker shape-preserving quadratic interpolation spline.

n = length(x)

if (ff == "Not-Supplied"){
  # Judd (1998), page 233, second last equation
  L = sqrt( (x[2:n]-x[1:(n-1)])^2 + (y[2:n]-y[1:(n-1)])^2)
  # Judd (1998), page 233, last equation
  d = (y[2:n]-y[1:(n-1)])/(x[2:n]-x[1:(n-1)])
  # Judd (1998), page 234, Eqn 6.11.6
  Conditionsi = (d[1:(n-2)]*d[2:(n-1)] > 0)
  MiddleSiwithoutApplyingCondition = (L[1:(n-2)]*d[1:(n-2)]+L[2:(n-1)] * d[2:(n-1)]) / (L[1:(n-2)]+L[2:(n-1)])
  sb = Conditionsi * MiddleSiwithoutApplyingCondition
  # Judd (1998), page 234, Second Equation line plus 6.11.6 gives this array of slopes.
  ff = c(((-sb[1]+3*d[1])/2),  sb,  ((3*d[n-1]-sb[n-2])/2))
}

NumberOfIntervalsWithKnots = 2*(n-1)

Intervals = 1:(n-1)
IntervalTab = data.frame(IntervalNum = sort(rep(Intervals,2)),
                                    SubIntervalNum = rep(c(1,2),n-1),
                                    StartOfInterval = numeric(NumberOfIntervalsWithKnots),
                                    EndOfInterval = numeric(NumberOfIntervalsWithKnots)
                                    )

Evals = do.call(rbind, lapply(Intervals, function(IntervalNum) SchumakerIndInterval(c(y[IntervalNum], y[IntervalNum+1]), c(ff[IntervalNum], ff[IntervalNum+1]), c(x[IntervalNum], x[IntervalNum+1]))))

IntervalTab = cbind(IntervalTab, Evals)
rm(Evals)


IntervalTab[IntervalTab$SubIntervalNum == 1, "StartOfInterval"] = x[1:n-1]
IntervalTab[IntervalTab$SubIntervalNum == 2, "EndOfInterval"] = x[2:n]
IntervalTab[IntervalTab$SubIntervalNum == 2, "StartOfInterval"] = IntervalTab$tsi[IntervalTab$SubIntervalNum == 2]
IntervalTab[IntervalTab$SubIntervalNum == 1, "EndOfInterval"] = IntervalTab$tsi[IntervalTab$SubIntervalNum == 1]
# This gets rid of empty intervals. The 1e-10 is there in case of numerical imprecision.
IntervalTab <- IntervalTab[which(IntervalTab$EndOfInterval + 1e-10 > IntervalTab$StartOfInterval),]

if ((Extrapolation %in% c("Constant", "Linear"))){
  Botx = min(IntervalTab$StartOfInterval)
  Boty   = y[1]

  BotInt = findInterval(Botx, IntervalTab$StartOfInterval)
  BotB = IntervalTab[BotInt ,c("B")]
  BotC   = Boty - BotB

    if (Extrapolation == "Constant"){

      BotB = 0
      BotC   = Boty
    }

  BotRow = data.frame(IntervalNum = 0, SubIntervalNum = 0, StartOfInterval = Botx-1, EndOfInterval = Botx, tsi = 0, C = 0, B =  BotB, A = BotC)


  Topx = max(IntervalTab$EndOfInterval)
  Topy   = y[n]

  TopInt = findInterval(Topx, IntervalTab$StartOfInterval)

  TopB = IntervalTab[TopInt ,c("B")]
  TopC   = Topy

  if (Extrapolation == "Constant"){
    TopB = 0
    TopC   = Topy
  }

  TopRow = data.frame(IntervalNum = 0, SubIntervalNum = 0,StartOfInterval = Topx, EndOfInterval = Topx + 1, tsi = 0, C = 0, B =  TopB, A = TopC)

  IntervalTab = rbind(BotRow, IntervalTab, TopRow)
}

# It is important use individual vectors and matrices rather than datatables or data.frames for speed
IntStarts = c(IntervalTab$StartOfInterval, Inf)
SpCoefs = data.matrix(IntervalTab[,c("C", "B", "A")])

# This is the end spline which looks up the correct interval and evaluates with the correpsonding coefficients
Spline0 = ppmak(IntStarts, SpCoefs, Vectorised )
Spline1 = ppmakDeriv(IntStarts, SpCoefs, Vectorised )
Spline2 = ppmak2Deriv(IntStarts, SpCoefs, Vectorised )
# This just boosts the speed of evaluation by 5ish percent. Not essential.
CompiledSpline0 = compiler::cmpfun(Spline0)
CompiledSpline1 = compiler::cmpfun(Spline1)
CompiledSpline2 = compiler::cmpfun(Spline2)
return(list(Spline = CompiledSpline0, DerivativeSpline = CompiledSpline1, SecondDerivativeSpline = CompiledSpline2 , IntervalTab = IntervalTab))
}
