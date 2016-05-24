#' Creates a Curve Object containing pairs of Tenors with relevant rates and the interpolation function.
#' Also, methods for populating the object via a .csv file and the generation of the interpolation function via cubic splines are included.
#' @title  Curve Class
#' @param Tenors          The Tenors of the curve
#' @param Rates           The rates on the corresponding tenors
#' @param interp_function (Optional) The interpolation function of the curve. Can be populated via the 'CalcInterpPoints' method
#' @return An object of type Curve
#' @export
#' @author Tasos Grivas <tasos@@openriskcalculator.com>

Curve = setRefClass("Curve",
                  # the timebuckets grouping is only relevant for IRDs
                  fields = list(Tenors = "numeric",
                                Rates  = "numeric",
                                interp_function = "function"
                                ),

                  methods = list(
                    CalcInterpPoints = function(time_points)
                    {
                      interp_function <<- splinefun(Tenors, Rates, method="natural")
                      return(interp_function(time_points))
                    },
                    PopulateViaCSV = function(csvfilename)
                    {
                      raw_data <- read.csv(system.file("extdata", csvfilename, package = "xVA"),header=TRUE,stringsAsFactors = FALSE,strip.white=TRUE)
                      Tenors <<- as.numeric(raw_data[,1])
                      Rates <<- as.numeric(raw_data[,2])/100
                    }
                  )
                  )
