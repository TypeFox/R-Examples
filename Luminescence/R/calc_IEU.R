#' Apply the internal-external-uncertainty (IEU) model after Thomsen et al.
#' (2007) to a given De distribution
#'
#' Function to calculate the IEU De for a De data set.
#'
#' This function uses the equations of Thomsen et al. (2007).  The parameters a
#' and b are estimated from dose-recovery experiments.
#'
#' @param data \code{\linkS4class{RLum.Results}} or \link{data.frame}
#' (\bold{required}): for \code{data.frame}: two columns with De
#' \code{(data[,1])} and De error \code{(values[,2])}
#'
#' @param a \code{\link{numeric}}: slope
#'
#' @param b \code{\link{numeric}}: intercept
#'
#' @param interval \code{\link{numeric}}: fixed interval (e.g. 5 Gy) used for
#' iteration of Dbar, from the mean to Lowest.De used to create Graph.IEU
#' [Dbar.Fixed vs Z]
#'
#' @param decimal.point \code{\link{numeric}} (with default): number of decimal
#' points for rounding calculations (e.g. 2)
#'
#' @param plot \code{\link{logical}} (with default): plot output
#'
#' @param \dots further arguments (\code{trace, verbose}).
#'
#' @return Returns a plot (optional) and terminal output. In addition an
#' \code{\linkS4class{RLum.Results}} object is returned containing the
#' following element:
#'
#' \item{summary}{\link{data.frame} summary of all relevant model results.}
#' \item{data}{\link{data.frame} original input data} \item{args}{\link{list}
#' used arguments} \item{call}{\link{call} the function call}
#' \item{tables}{\link{list} a list of data frames containing all calculation
#' tables}
#'
#' The output should be accessed using the function
#' \code{\link{get_RLum}}.
#'
#' @section Function version: 0.1.0
#'
#' @author Rachel Smedley, Geography & Earth Sciences, Aberystwyth University
#' (United Kingdom) \cr Based on an excel spreadsheet and accompanying macro
#' written by Kristina Thomsen.
#'
#' @seealso \code{\link{plot}}, \code{\link{calc_CommonDose}},
#' \code{\link{calc_CentralDose}}, \code{\link{calc_FiniteMixture}},
#' \code{\link{calc_FuchsLang2001}}, \code{\link{calc_MinDose}}
#'
#' @references Smedley, R.K., 2015. A new R function for the Internal External Uncertainty (IEU) model.
#' Ancient TL 33, 16-21.
#'
#' Thomsen, K.J., Murray, A.S., Boetter-Jensen, L. & Kinahan, J.,
#' 2007. Determination of burial dose in incompletely bleached fluvial samples
#' using single grains of quartz. Radiation Measurements 42, 370-379.
#'
#' @examples
#'
#' ## load data
#' data(ExampleData.DeValues, envir = environment())
#'
#' ## apply the IEU model
#' ieu <- calc_IEU(ExampleData.DeValues$CA1, a = 0.2, b = 1.9, interval = 1)
#'
#' @export
calc_IEU <- function(
  data,
  a,
  b,
  interval,
  decimal.point = 2,
  plot = TRUE,
  ...
) {

  ##==========================================================================##
  ## CONSISTENCY CHECK OF INPUT DATA
  ##==========================================================================##
  if(missing(data)==FALSE){
    if(is(data, "data.frame") == FALSE & is(data,"RLum.Results") == FALSE){
      stop("[calc_IEU] Error: 'data' object has to be of type
           'data.frame' or 'RLum.Results'!")
    }else{
      if(is(data, "RLum.Results") == TRUE){
        data <- get_RLum(data, signature(object = "De.values"))
      }
    }
  }

  ##==========================================================================##
  ## ... ARGUMENTS
  ##==========================================================================##
  extraArgs <- list(...)
  ## console output
  if ("verbose" %in% names(extraArgs)) {
    verbose <- extraArgs$verbose
  } else {
    verbose <- TRUE
  }
  # trace calculations
  if ("trace" %in% names(extraArgs)) {
    trace <- extraArgs$trace
  } else {
    trace <- FALSE
  }
  # TODO: main, xlab, ylab, xlim, ylim, pch, col


  ##============================================================================##
  ## CALCULATIONS
  ##============================================================================##
  empty <- NULL
  Table.Fixed.Iteration <- data.frame(matrix(nrow = 0, ncol = 9))
  colnames(data) <- c("De", "De.Error")
  data <- data[order(data$De), ]
  Mean <- mean(data$De)
  Dbar <- round(Mean, decimal.point)
  Lowest.De <- round(data$De[1], decimal.point)

  # (a) Calculate IEU at fixed intervals of Dbar starting from the Mean and
  # subtracting the interval until Dbar is < Lowest.De; this creates a plot
  N <- nrow(data)
  Rank.number <- t(c(1:N))
  De.Total.Error <- sqrt((data$De.Error^2) + (((a * Dbar) + b)^2))
  Table.Calculations <- data.frame(Rank.number = c(Rank.number),
                                   De = c(data$De),
                                   De.Total.Error = c(De.Total.Error))
  Z.top <- cumsum(Table.Calculations$De/(Table.Calculations$De.Total.Error^2))
  Z.bottom <- cumsum(1/(Table.Calculations$De.Total.Error^2))
  Z <- Z.top/Z.bottom
  Table.Calculations["Z"] <- Z

  temp <- NULL
  for (j in 1:N) {
    for (i in j) {
      Z <- Table.Calculations$Z[j]
      x <- ((Table.Calculations$De[1:i] - Z)^2)/((Table.Calculations$De.Total.Error[1:i])^2)
      y <- (sum(x))
      temp <- rbind(temp, data.frame(y))
    }
  }

  EXT.top <- temp
  EXT.bottom <- (Table.Calculations$Rank.number - 1) * Z.bottom
  EXT <- EXT.top/EXT.bottom
  INT <- 1/Z.bottom
  R <- sqrt(INT/EXT)
  R.Error <- (2 * (Table.Calculations$Rank.number - 1))^(-0.5)

  Table.IEU <- data.frame(Table.Calculations$Rank.number, Table.Calculations$De,
                          Table.Calculations$De.Total.Error, Table.Calculations$Z,
                          EXT.top, EXT, INT, R, R.Error)

  colnames(Table.IEU) <- c("Rank.number", "De", "De.Error", "Z", "EXT.top",
                           "EXT", "INT", "R", "R.Uncertainty")

  Unity <- Table.IEU[R >= 1, ]
  Max <- max(Unity$Rank.number, na.rm = TRUE)
  Above.Z <- Table.IEU[Max, 4]
  Above.Error <- Table.IEU[Max, 6]
  Below.Z <- Table.IEU[Max + 1, 4]
  Below.Error <- Table.IEU[Max + 1, 6]
  Above.R <- Table.IEU[Max, 8]
  Below.R <- Table.IEU[Max + 1, 8]
  Slope <- (Above.R - Below.R)/(Above.Z - Below.Z)
  Intercept <- Above.R - (Slope * Above.Z)
  IEU.De <- round(((1 - Intercept)/Slope), decimal.point)
  IEU.Error <- max(sqrt(Above.Error), sqrt(Below.Error))
  IEU.Error <- round(IEU.Error, decimal.point)
  n <- Max + 1

  Dbar.Fixed <- Dbar - interval
  Dbar.Mean <- c(1, Dbar, Dbar.Fixed, IEU.De, IEU.Error, n, Below.R, a, b)

  repeat {
    if (Dbar.Fixed < Lowest.De) {
      break
    } else {
      Dbar <- Dbar.Fixed
    }
    De.Total.Error <- sqrt((data$De.Error^2) + (((a * Dbar) + b)^2))
    Table.Calculations <- data.frame(Rank.number = c(Rank.number),
                                     De = c(data$De),
                                     De.Total.Error = c(De.Total.Error))
    Z.top <- cumsum(Table.Calculations$De/(Table.Calculations$De.Total.Error^2))
    Z.bottom <- cumsum(1/(Table.Calculations$De.Total.Error^2))
    Z <- Z.top/Z.bottom
    Table.Calculations["Z"] <- Z

    temp <- NULL
    for (j in 1:N) {
      for (i in j) {
        Z <- Table.Calculations$Z[j]
        x <- ((Table.Calculations$De[1:i] - Z)^2)/((Table.Calculations$De.Total.Error[1:i])^2)
        y <- (sum(x))
        temp <- rbind(temp, data.frame(y))
      }
    }

    EXT.top <- temp
    EXT.bottom <- (Table.Calculations$Rank.number - 1) * Z.bottom
    EXT <- EXT.top/EXT.bottom
    INT <- 1/Z.bottom
    R <- sqrt(INT/EXT)
    R.Error <- (2 * (Table.Calculations$Rank.number - 1))^(-0.5)

    Table.IEU <- data.frame(Table.Calculations$Rank.number, Table.Calculations$De,
                            Table.Calculations$De.Total.Error, Table.Calculations$Z,
                            EXT.top, EXT, INT, R, R.Error)

    colnames(Table.IEU) <- c("Rank.number", "De", "De.Error", "Z", "EXT.top",
                             "EXT", "INT", "R", "R.Uncertainty")

    Unity <- Table.IEU[R >= 1, ]
    Max <- max(Unity$Rank.number, na.rm = TRUE)
    Above.Z <- Table.IEU[Max, 4]
    Above.Error <- Table.IEU[Max, 6]
    Below.Z <- Table.IEU[Max + 1, 4]
    Below.Error <- Table.IEU[Max + 1, 6]
    Above.R <- Table.IEU[Max, 8]
    Below.R <- Table.IEU[Max + 1, 8]
    Slope <- (Above.R - Below.R)/(Above.Z - Below.Z)
    Intercept <- Above.R - (Slope * Above.Z)
    Zbar <- round(((1 - Intercept)/Slope), decimal.point)
    Zbar.Error <- max(sqrt(Above.Error), sqrt(Below.Error))
    Zbar.Error <- round(IEU.Error, decimal.point)
    n <- Max + 1
    Dbar.Fixed <- Dbar - interval
    Table.Fixed.Iteration <- rbind(Table.Fixed.Iteration,
                                   cbind(1, Dbar, Dbar.Fixed, Zbar, Zbar.Error,
                                         n, Below.R, a, b))
  }

  Table.Fixed.Iteration <- rbind(Dbar.Mean, Table.Fixed.Iteration)
  colnames(Table.Fixed.Iteration) <- c(FALSE, "Dbar", "Dbar.Fixed", "Zbar",
                                       "Zbar.Error", "n", "Below.R", "a", "b")

  if (plot) {
    plot(Table.Fixed.Iteration$Dbar,
         Table.Fixed.Iteration$Zbar,
         type = "b",
         ylab = "Zbar, weighted mean  (Gy)",
         xlab = "Dbar (Gy)",
         asp = 1/1)

    arrows(Table.Fixed.Iteration$Dbar, Table.Fixed.Iteration$Zbar + Table.Fixed.Iteration$Zbar.Error,
           Table.Fixed.Iteration$Dbar, Table.Fixed.Iteration$Zbar - Table.Fixed.Iteration$Zbar.Error,
           col = 1, angle = 90, length = 0.05, code = 3)

    abline(0, 1, untf = FALSE, lty = 3)
  }

  # (b) Calculate Dbar by iteration from [Dbar = Lowest.De] until [IEU.De = Dbar];
  # this calculates the IEU De
  Dbar <- Lowest.De
  N <- nrow(data)
  Rank.number <- t(c(1:N))
  De.Total.Error <- sqrt((data$De.Error^2) + (((a * Dbar) + b)^2))
  Table.Calculations <- data.frame(Rank.number = c(Rank.number),
                                   De = c(data$De),
                                   De.Total.Error = c(De.Total.Error))
  Z.top <- cumsum(Table.Calculations$De/(Table.Calculations$De.Total.Error^2))
  Z.bottom <- cumsum(1/(Table.Calculations$De.Total.Error^2))
  Z <- Z.top/Z.bottom
  Table.Calculations["Z"] <- Z

  temp <- NULL
  for (j in 1:N) {
    for (i in j) {
      Z <- Table.Calculations$Z[j]
      x <- ((Table.Calculations$De[1:i] - Z)^2)/((Table.Calculations$De.Total.Error[1:i])^2)
      y <- (sum(x))
      temp <- rbind(temp, data.frame(y))
    }
  }

  EXT.top <- temp
  EXT.bottom <- (Table.Calculations$Rank.number - 1) * Z.bottom
  EXT <- EXT.top/EXT.bottom
  INT <- 1/Z.bottom
  R <- sqrt(INT/EXT)
  R.Error <- (2 * (Table.Calculations$Rank.number - 1))^(-0.5)

  Table.IEU <- data.frame(Table.Calculations$Rank.number, Table.Calculations$De,
                          Table.Calculations$De.Total.Error, Table.Calculations$Z,
                          EXT.top, EXT, INT, R, R.Error)

  colnames(Table.IEU) <- c("Rank.number", "De", "De.Error", "Z",
                           "EXT.top", "EXT", "INT", "R", "R.Uncertainty")

  Unity <- Table.IEU[R >= 1, ]
  Max <- max(Unity$Rank.number, na.rm = TRUE)
  Above.Z <- Table.IEU[Max, 4]
  Above.Error <- Table.IEU[Max, 6]
  Below.Z <- Table.IEU[Max + 1, 4]
  Below.Error <- Table.IEU[Max + 1, 6]
  Above.R <- Table.IEU[Max, 8]
  Below.R <- Table.IEU[Max + 1, 8]
  Slope <- (Above.R - Below.R)/(Above.Z - Below.Z)
  Intercept <- Above.R - (Slope * Above.Z)
  IEU.De <- round(((1 - Intercept)/Slope), decimal.point)
  IEU.Error <- max(sqrt(Above.Error), sqrt(Below.Error))
  IEU.Error <- round(IEU.Error, decimal.point)
  n <- Max + 1

  repeat {
    if (IEU.De <= Dbar) {
      break
    } else {
      Dbar <- IEU.De
    }
    De.Total.Error <- sqrt((data$De.Error^2) + (((a * Dbar) + b)^2))
    Table.Calculations <- data.frame(Rank.number = c(Rank.number),
                                     De = c(data$De),
                                     De.Total.Error = c(De.Total.Error))
    Z.top <- cumsum(Table.Calculations$De/(Table.Calculations$De.Total.Error^2))
    Z.bottom <- cumsum(1/(Table.Calculations$De.Total.Error^2))
    Z <- round((Z.top/Z.bottom), decimal.point)
    Table.Calculations["Z"] <- Z

    temp <- NULL
    for (j in 1:N) {
      for (i in j) {
        Z <- Table.Calculations$Z[j]
        x <- ((Table.Calculations$De[1:i] - Z)^2)/((Table.Calculations$De.Total.Error[1:i])^2)
        y <- (sum(x))
        temp <- rbind(temp, data.frame(y))
      }
    }

    EXT.top <- temp
    EXT.bottom <- (Table.Calculations$Rank.number - 1) * Z.bottom
    EXT <- EXT.top/EXT.bottom
    INT <- 1/Z.bottom
    R <- sqrt(INT/EXT)
    R.Error <- (2 * (Table.Calculations$Rank.number - 1))^(-0.5)

    Table.IEU <- data.frame(Table.Calculations$Rank.number, Table.Calculations$De,
                            Table.Calculations$De.Total.Error, Table.Calculations$Z,
                            EXT.top, EXT, INT, R, R.Error)

    colnames(Table.IEU) <- c("Rank.number", "De", "De.Error", "Z", "EXT.top",
                             "EXT", "INT", "R", "R.Error")

    # to reduce the number of plots and increase perfomance
    # intermediate calculations are only plotted when trace = TRUE
    if (plot && trace) {
      ymin <- min(Table.IEU$R[2:nrow(Table.IEU)] - Table.IEU$R.Error[2:nrow(Table.IEU)])
      ymax <- max(Table.IEU$R[2:nrow(Table.IEU)] + Table.IEU$R.Error[2:nrow(Table.IEU)])
      ylim <- c(ifelse(ymin > 0, 0, ymin), ymax)

      plot(Table.IEU$Z, Table.IEU$R,
           type = "b",
           ylab = expression(paste("R = [", alpha["in"], "/", alpha["ex"],"]")),
           xlab = "Z [Gy]",
           ylim = ylim)

      arrows(Table.IEU$Z, Table.IEU$R + Table.IEU$R.Error,
             Table.IEU$Z, Table.IEU$R - Table.IEU$R.Error,
             col = 1, angle = 90,
             length = 0.05, code = 3)

      abline(1, 0, untf = FALSE, lty = 3)
    }

    Unity <- Table.IEU[R >= 1, ]
    Max <- max(Unity$Rank.number, na.rm = TRUE)
    Above.Z <- Table.IEU[Max, 4]
    Above.Error <- Table.IEU[Max, 6]
    Below.Z <- Table.IEU[Max + 1, 4]
    Below.Error <- Table.IEU[Max + 1, 6]
    Above.R <- Table.IEU[Max, 8]
    Below.R <- Table.IEU[Max + 1, 8]
    Slope <- (Above.R - Below.R)/(Above.Z - Below.Z)
    Intercept <- Above.R - (Slope * Above.Z)
    IEU.De <- round(((1 - Intercept)/Slope), decimal.point)
    IEU.Error <- max(sqrt(Above.Error), sqrt(Below.Error))
    IEU.Error <- round(IEU.Error, decimal.point)
    n <- Max + 1

    if (trace) {
      message(sprintf("[Iteration of Dbar] \n Dbar: %.4f \n IEU.De: %.4f \n IEU.Error: %.4f \n n: %i \n R: %.4f \n",
                      Dbar, IEU.De, IEU.Error, n, Below.R))
    }

  }

  # final plot
  if (plot) {
    ymin <- min(Table.IEU$R[2:nrow(Table.IEU)] - Table.IEU$R.Error[2:nrow(Table.IEU)])
    ymax <- max(Table.IEU$R[2:nrow(Table.IEU)] + Table.IEU$R.Error[2:nrow(Table.IEU)])
    ylim <- c(ifelse(ymin > 0, 0, ymin), ymax)

    plot(Table.IEU$Z, Table.IEU$R,
         type = "b",
         ylab = expression(paste("R = [", alpha["in"], "/", alpha["ex"],"]")),
         xlab = "Z [Gy]",
         ylim = ylim)

    arrows(Table.IEU$Z, Table.IEU$R + Table.IEU$R.Error,
           Table.IEU$Z, Table.IEU$R - Table.IEU$R.Error,
           col = 1, angle = 90,
           length = 0.05, code = 3)

    abline(1, 0, untf = FALSE, lty = 3)
  }


  Table.Results <- data.frame(Dbar, IEU.De, IEU.Error, n, a, b)
  colnames(Table.Results) <- c("Dbar", "IEU.De (Gy)", "IEU.Error (Gy)",
                               "Number of De", "a", "b")

  ##==========================================================================##
  ## TERMINAL OUTPUT
  ##==========================================================================##
  if (verbose) {
    message(sprintf(
      "\n [calc_IEU] \n\n Dbar: %.2f \n IEU.De (Gy): %.2f \n IEU.Error (Gy): %.2f Number of De: %.0f \n a: %.4f \n b: %.4f",
      Table.Results[1], Table.Results[2], Table.Results[3],
      Table.Results[4], Table.Results[5], Table.Results[6]))
  }

  ##==========================================================================##
  ## RETURN VALUES
  ##==========================================================================##
  summary <- Table.Results[ ,c(-1, -5, -6)]
  colnames(summary) <- c("de", "de_err", "n")

  call <- sys.call()
  args <- list(a = a, b = b, interval = interval,
               decimal.point = decimal.point, plot = plot)

  newRLumResults.calc_IEU <- set_RLum(
    class = "RLum.Results",
    data = list(summary = summary,
                data = data,
                args = args,
                call = call,
                tables = list(
                  Table.IEUCalculations = Table.IEU,
                  Table.Fixed.Iteration = Table.Fixed.Iteration,
                  Table.IEUResults = Table.Results
                )))

  invisible(newRLumResults.calc_IEU)

}
