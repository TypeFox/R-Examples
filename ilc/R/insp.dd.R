insp.dd <-
function(data, what=c('rate','pop', 'deaths'),
                    ages=data$age, years=data$year, series = names(data$rate)[1]){
    if (class(data) != "demogdata" || data$type != "mortality")
        stop("Not mortality data")
    what <- match.arg(what)
    if (what == 'deaths'){
        data <- extract.deaths(data, ages, years, combine.upper=F, series=series)
        what <- 'rate'
    }
    else {
        data <- extract.years(data, years)
        data <- extract.ages(data, ages=ages, combine.upper=F)
    }
    data[[what]][[series]]
}
