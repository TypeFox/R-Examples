##' Create a lag 
##'
##' Lag variables by an arbitrary number of periods even if the data is grouped
##'
##' @param df A dataframe with groups, time periods, and a variable to be lagged
##' @param group The grouping factor in the dataframe
##' @param values The names of the variables to be lagged
##' @param time The variable representing time periods
##' @param periods A scalar for the number of periods to be lagged in the data. Can 
##' be negative to indicate leading variable.
##' @return A dataframe with a newly created variable lagged
##' @export
##' @examples
##' 
##' test_data <- expand.grid(id = sample(letters, 10), 
##'                         time = 1:10)
##' test_data$value1 <- rnorm(100)
##' test_data$value2 <- runif(100)
##' test_data$value3 <- rpois(100, 4)
##' group <- "id"
##' time <- "time"
##' values <- c("value1", "value2")
##' vars <- c(group, time, values)
##' periods <- 2
##' newdat <- lag_data(test_data, group="id", time="time", 
##'                  values=c("value1", "value2"), periods=3)
lag_data <- function(df, group, time, periods, values) {
  if(any(group %in% values) | any(time %in% values)){
    stop("Cannot lag a grouping or a time variable")
  }
  if(!class(df[, time]) %in% c("integer", "numeric")){
    stop("Time must be a numeric or integer value")
  }
  if(!class(periods) %in% c("integer", "numeric")){
    stop("Periods must be a numeric or integer value")
  }
  if(periods - as.integer(periods) != 0){
    periods <- as.integer(periods)
    warning(paste0("parameter periods has been forced to integer of value ", 
                   periods))
  }
  vars <- c(group, time, values)
  tmp.data <- df[, vars]
  tmp.data[, time] <- tmp.data[, time] + periods
  lab <- ifelse(periods < 0, ".lead", ".lag")
  new.vals <- paste0(values, lab, abs(periods))
  names(tmp.data) <- c(group, time, new.vals)
  newdat <- merge(df, tmp.data, all.x=TRUE, by=c(group,time))
  return(newdat)
}

