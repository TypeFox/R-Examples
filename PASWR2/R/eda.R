#' @title Exploratory Data Analysis
#' 
#' @description Function that produces a histogram, density plot, boxplot, and Q-Q plot
#' 
#' @details The function \code{eda()} will not return console window information on data sets containing more than 5000 observations. It will, however, still produce graphical output for data sets containing more than 5000 observations.
#' 
#' @param x is a numeric vector where \code{NA}s and \code{Inf}s are allowed but will be removed.
#' @param trim is a fraction (between 0 and 0.5, inclusive) of values to be trimmed from each end of the ordered data such that if \code{trim = 0.5}, the result is the median.
#' @param dec is a number specifying the number of decimals
#' 
#' @return Function returns various measures of center and location. The values returned for the quartiles are based on the default \pkg{R} definitions for quartiles. For more information on the definition of the quartiles, type \code{?quantile} and read about the algorithm used by \code{type = 7}.
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu> 
#' 
#' @export
#' 
#' @examples
#' eda(x = rnorm(100))
#' # Produces four graphs for the 100 randomly
#' # generated standard normal variates.
#' 
#' @keywords hplot
#######################################################################
eda <- function(x, trim = 0.05, dec = 3)
{
  # require(e1071)
  # rgb(0, 128/255, 1, names="Adkblue") #Alan's dark blue
  # rgb(169/255, 226/255, 1, names="Altblue") #Alan's light blue
  opar <- par(no.readonly = TRUE)
  Altblue <- "#A9E2FF"
  Adkblue <- "#0080FF"
  Ared <- "#C51111"
  varname <- deparse(substitute(x))
  N <- length(x)
  UM <- round(sum(is.na(x)), 0)
  x <- x[is.finite(x)]
  n <- length(x)
  Q1 <- round(quantile(x)[2], digits = dec)
  Q3 <- round(quantile(x)[4], digits = dec)
  IQR <- round(Q3 - Q1, digits = dec)
  Min <- round(min(x), digits = dec)
  Max <- round(max(x), digits = dec)
  Stdev <- round(sd(x), digits = dec)
  Mean <- round(mean(x), digits = dec)
  Median <- round(median(x), digits = dec)
  TrMean <- round(mean(x, trim = trim), digits = dec)
  Var <- round(var(x), digits = dec)
  SE <- round(Stdev/sqrt(n), digits = dec)
  Range <- round(Max - Min, digits = dec)
  par(oma = c(0, 1, 4, 1))
  par(mfrow = c(2, 2))
  par(mar = c(1, 1, 2, 1))
  hist(x, freq = FALSE, col = Adkblue, xlab = "", ylab = "", axes = FALSE,
       main = paste("Histogram of", varname), breaks = "Scott" )
  box()
  iqd <- summary(x)[5] - summary(x)[2]
  plot(density(x), xlab = "", ylab = "", 
       axes = FALSE, type = "n", main = paste("Density of", varname))
  lines(density(x), col=Ared)
  box()
  LH <- fivenum(x)[2]
  UH <- fivenum(x)[4]
  HS <- UH - LH
  l.out <- x[x < (LH - 1.5 * HS)]
  r.out <- x[x > (UH + 1.5 * HS)]
  outliers <- c(l.out, r.out)
  rest <- x[x > (LH - 1.5 * HS) & x < (UH + 1.5 * HS)]
  Minrest <- min(rest)
  Maxrest <- max(rest)
  plot(x, x, main = paste("Boxplot of", varname), xlab =
         "", ylab = "", axes = FALSE, type = "n", xlim = c(min(x), max(
           x)), ylim = c(0, 1))
  box()
  polygon(c(LH, LH, UH, UH), c(0.3, 0.7, 0.7, 0.3), density = -1, col= Altblue)
  points(outliers, c(rep(0.5, length(outliers))), col = Ared)
  lines(c(min(rest), LH), c(0.5, 0.5), lty = 1)
  lines(c(UH, max(rest)), c(0.5, 0.5), lty = 1)
  lines(c(min(rest), min(rest)), c(0.4, 0.6))
  lines(c(max(rest), max(rest)), c(0.4, 0.6))
  lines(c(LH, LH), c(0.3, 0.7))
  lines(c(UH, UH), c(0.3, 0.7))
  lines(c(Median, Median), c(0.3, 0.7))
  lines(c(LH, UH), c(0.3, 0.3))
  lines(c(LH, UH), c(0.7, 0.7))
  points(Mean, 0.5, pch = 16, col = "black")
  qqnorm(x, col = "black", main = paste("Q-Q Plot of", varname), xlab = "",
         ylab = "", axes = FALSE)
  qqline(x, col = Ared)
  box()
  mtext("EXPLORATORY  DATA  ANALYSIS", side = 3, outer = TRUE, cex = 1.5,
        col = Adkblue, line = 1)
  on.exit(par(opar))
  SW <- shapiro.test(x)
  K <- round(kurtosis(x), digits = dec)
  S <- round(skewness(x), digits = dec)
  SWpval <- round(SW$p.value, digits = dec)
  TOT <- c(n, UM, Min, Q1, Mean, Median, TrMean, Q3, Max, Stdev, Var,
           SE, IQR, Range, K, S, SWpval)
  names(TOT) <- c("Size (n)", "Missing", "Minimum", " 1st Qu", "   Mean",
                  " Median", "TrMean", " 3rd Qu", "   Max", " Stdev",
                  "   Var", "SE Mean", " I.Q.R.", "  Range", "Kurtosis",
                  "Skewness", "SW p-val")
  return(TOT)
}
