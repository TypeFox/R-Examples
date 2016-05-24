#' Prepare HMD data for \code{lexis.hmd()}
#' 
#' \code{prepare.hmd()} prepares the raw 'Deaths by Lexis triangles' HMD data for further use by \code{lexis.hmd}.
#' 
#' @param file, the name of the 'Deaths by Lexis triangles' file downloaded from the Human Mortality Database.
#' @details This function reads the raw data into R and transforms data to \code{numeric} and \code{Date}.
#' Furthermore seven columns (\code{upper, x1, x2, x3, y1, y2, y3}) that contain the coordinates of the triangles will be added.
#' The age group \code{110+} will be removed from the data.
#' @author Philipp Ottolinger
#' @importFrom utils read.csv
#' @importFrom stats complete.cases
#' @export prepare.hmd
#' @examples
#' library(LexisPlotR)
#' # Load sample data
#' path <- system.file("extdata", "Deaths_lexis_sample.txt", package = "LexisPlotR")
#' deaths.triangles <- prepare.hmd(path)

prepare.hmd <- function(file) {
  data <- read.csv(file, sep="", skip = 2)
  data$Year <- as.numeric(as.character(data$Year))
  data$Age <- as.numeric(as.character(data$Age))
  data$Cohort <- as.numeric(as.character(data$Cohort))
  data$upper <- ifelse(data$Year - data$Age > data$Cohort, TRUE, FALSE)
  data <- data[complete.cases(data),]
  data$x1 <- ifelse(data$upper == FALSE, data$Year, data$Year)
  data$x2 <- ifelse(data$upper == FALSE, data$Year + 1, data$Year)
  data$x3 <- ifelse(data$upper == FALSE, data$Year + 1, data$Year + 1)
  data$y1 <- ifelse(data$upper == FALSE, data$Age, data$Age)
  data$y2 <- ifelse(data$upper == FALSE, data$Age, data$Age + 1)
  data$y3 <- ifelse(data$upper == FALSE, data$Age + 1, data$Age + 1)
  data$x1 <- as.Date(paste(data$x1, "-01-01", sep = ""), origin = "1970-01-01")
  data$x2 <- as.Date(paste(data$x2, "-01-01", sep = ""), origin = "1970-01-01")
  data$x3 <- as.Date(paste(data$x3, "-01-01", sep = ""), origin = "1970-01-01")
  return(data)
}