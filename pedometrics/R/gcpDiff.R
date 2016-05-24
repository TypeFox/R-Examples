#' Difference on xyz coordinates between ground control points
#' 
#' This function estimates the difference, absolute difference, and squared
#' difference on x, y and z coordinates of two sets of ground control points
#' (GCP). It also estimates the module (difference vector), its square and
#' azimuth. The result is a data frame ready to be used to define a object of
#' class \code{spsurvey.object}.
#' 
#' @details
#' \subsection{Type of data}{
#' Two types of validation data that can be submitted to function
#' \code{gcpDiff()}: those coming from horizontal (positional) validation
#' exercises (\code{type = "xy"}), and those coming from vertical validation
#' exercises (\code{type = "z"}).
#' 
#' Horizontal (positional) validation exercises compare the position of
#' \code{measured} point data with the position of \code{predicted} point data.
#' Horizontal displacement (error) is measured in both \sQuote{x} and
#' \sQuote{y} coordinates, and is used to calculate the error vector (module)
#' and its azimuth.  Both objects \code{measured} and \code{predicted} used
#' with function \code{gcpDiff()} must be of class
#' \code{SpatialPointsDataFrame}. They must have at least one column named
#' \sQuote{siteID} giving the identification of every case. Matching of case
#' IDs is mandatory. Other columns are discarded.
#' 
#' Vertical validation exercises are interested in comparing the
#' \code{measured} value of a variable at a given location with that
#' \code{predicted} by some model. In this case, error statistics are
#' calculated only for the the vertical displacement (error) in the \sQuote{z}
#' coordinate. Both objects \code{measured} and \code{predicted} used with
#' function \code{gcpDiff()} must be of class \code{SpatialPointsDataFrame}.
#' They also must have a column named \sQuote{siteID} giving the identification
#' of evary case. Again, matching of case IDs is mandatory. However, both
#' objects must have a column named \sQuote{z} which contains the values of the
#' \sQuote{z} coordinate. Other columns are discarded.
#' }
#' \subsection{Data aggregation}{
#' Validation is sometimes performed using cluster or transect sampling. Before
#' estimation of error statistics, the data needs to be aggregated by cluster
#' or transect. The function \code{gcpDiff()} aggregates validation data of
#' \code{type = "z"} calculating the mean value per cluster. Thus, aggregation
#' can only be properly done if the \sQuote{siteID} column of both objects
#' \code{measured} and \code{predicted} provides the identification of
#' clusters.  Setting \code{aggregate = TRUE} will return aggregated estimates
#' of error statistics. If the data has been aggregated beforehand, the
#' parameter \code{aggregate} can be set to \code{FALSE}.
#' }
#' \subsection{Case matching}{
#' There are circumstances in which the number of cases in the object
#' \code{measured} is larger than that in the object \code{predicted}. The
#' function \code{gcpDiff()} compares the number of cases in both objects and
#' automatically drops those cases of object \code{measured} that do not match
#' the cases of object \code{predicted}. However, case matching can only be
#' done if case IDs are exactly the same for both objects. Otherwise, estimated
#' error statistics will have no meaning at all.
#' }
#' 
#' @param measured Object of class \code{SpatialPointsDataFrame}
#' with the reference GCP. A column named \sQuote{siteID} giving case names is
#' mandatory. See \sQuote{Details}, item \sQuote{Type of data}.
#' 
#' @param predicted An object of class
#' \code{SpatialPointsDataFrame} with the point data being
#' validated. A column named \sQuote{siteID} giving case names is mandatory.
#' See \sQuote{Details}, item \sQuote{Type of data}.
#' 
#' @param type Type of data under analysis. Defaults to \code{type = "xy"}.
#' \sQuote{Details}, item \sQuote{Type of data}.
#' 
#' @param aggregate Logical for aggregating the data when it comes from cluster
#' sampling. Used only when \code{type = "z"}. Defaults to \code{aggregate =
#' FALSE}. See \sQuote{Details}, item \sQuote{Data aggregation}.
#' 
#' @param rounding Rounding level of the data in the output data frame.
#' 
#' @return An object of class \code{data.frame} ready to be used to feed the
#' argument \code{data.cont} when creating a \code{spsurvey.analysis} object.
#' @note Data of \code{type = "xy"} cannot be submitted to cluster aggregation
#' in the present version.
#' @author Alessandro Samuel-Rosa <\email{alessandrosamuelrosa@@gmail.com}>
#' @seealso \code{\link[pedometrics]{coordenadas}},
#' \code{\link[pedometrics]{gcpVector}},
#' \code{\link[spsurvey]{spsurvey.analysis}}.
#' @references Kincaid, T. M. and Olsen, A. R. (2013). spsurvey: Spatial Survey
#' Design and Analysis. R package version 2.6. URL:
#' \url{http://www.epa.gov/nheerl/arm/}.
#' @keywords methods
#' @export
#' @examples
#' 
#' \dontrun{
#' ## Create an spsurvey.analysis object
#' my.spsurvey <- 
#'   spsurvey.analysis(design = coordenadas(my.data),
#'                     data.cont = delta(ref.data, my.data),
#'                     popcorrect = TRUE, pcfsize = length(my.data$id),
#'                     support = rep(1, length(my.data$id)),
#'                     wgt = rep(1, length(my.data$id)), vartype = "SRS")
#' }
#' 
gcpDiff <- 
  function(measured, predicted, type = "xy", aggregate = FALSE, rounding = 0) {
    if(type == "xy") {  # difference in the geographic space
      measured <- coordenadas(measured)
      predicted <- coordenadas(predicted)
      # Because the number of GCPs of 'measured' and 'predicted'
      # may not be the same
      measured <- measured[match(predicted[, 1], measured[, 1]), ][, 2:3]
      diff <- predicted[, 2:3] - measured
      names(diff) <- c("dx", "dy")
      azim <- gcpVector(dx = diff$dx, dy = diff$dy)
      abso <- abs(diff)
      colnames(abso) <- c("abs.dx", "abs.dy")
      quad <- diff * diff
      colnames(quad) <- c("sq.dx", "sq.dy")
      diff <- round(diff, rounding)
      abso <- round(abso, rounding)
      quad <- round(quad, rounding)
      azim <- round(azim, rounding)
      siteID <- predicted$siteID
      erro <- data.frame(siteID, diff, abso, quad, azim)
      erro <- erro[order(as.numeric(row.names(erro))), ]
      row.names(erro) <- NULL
      return(erro)
    }
    if(type == "z") {  # difference in the attribute space
      measured <- measured
      predicted <- predicted
      dz <- predicted$z - measured$z
      abs.dz <- abs(dz)
      sq.dz <- dz*dz
      dz <- round(dz, rounding)
      abs.dz <- round(abs.dz, rounding)
      sq.dz <- round(sq.dz, rounding)
      siteID <- predicted$siteID
      if(aggregate == FALSE) {
        erro <- data.frame(siteID, dz, abs.dz, sq.dz)
        row.names(erro) <- NULL
        return(erro)
      }
      if(aggregate == TRUE) {
        erro <- data.frame(dz, abs.dz, sq.dz)
        by <- list(siteID)
        names(by) <- "siteID"
        erro <- aggregate(erro, by, FUN = mean)
        #erro = data.frame(siteID = seq(1:length(dz)), erro)
        row.names(erro) <- NULL
        return(erro)
      }
    }
    print(erro)
  }
# End!
