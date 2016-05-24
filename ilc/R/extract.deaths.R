extract.deaths <-
function(data, ages=data$age, years=data$year,
                           combine.upper = T, fill.method=NULL, series = names(data$rate)[1]){
    if (class(data) != "demogdata" || data$type != "mortality")
        stop("Not mortality data")
    data <- extract.years(data, years)
    data <- extract.ages(data, ages=ages, combine.upper)
    if (!is.null(fill.method)) data <- fill.demogdata(data, method=fill.method)
    # correction for NA rates which are due to 0 exposures:
    na.ind <- !bool(data$pop[[series]])
    data$rate[[series]][na.ind] <- 0
    # Note: NA exposures will simply result in NA deaths.
    data$rate[[series]] <- data$rate[[series]]*data$pop[[series]]
    data$type <- 'deaths'  # to avoid conflicts with demogdata functions
    # though it will stop most demography functions working
    data
}
