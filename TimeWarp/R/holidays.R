.holidays <- new.env()

holidays <- function(years, type, silent = FALSE)
{
    if (!all(type %in% c("NONE", allHolidays())))
        stop(paste("No '", type[!(type %in% allHolidays())][1],
                   "' holidays exist.", sep = ""))

    extractHolidays <- function(type)
    {
        if (type == "NONE")
            return(emptyDate())

        hDays <- .holidays[[type]]

        if (!silent && any(years < min(hDays$years)))
            warning(paste("Do not have any '", type,
                          "' holiday data for year(s) ",
                          paste(unique(years[years < min(hDays$years)])
                                , collapse = ", "), ".", sep = ""))

        if (!silent && any(years > max(hDays$years)))
            warning(paste("Do not have any '", type,
                          "' holiday data for year(s) ",
                          paste(unique(years[years > max(hDays$years)])
                                , collapse = ", "), ".", sep = ""))

        hDays[hDays$years %in% years, 'days']
    }

    structure(unique(sort(unlist(sapply(type, extractHolidays,
                                        USE.NAMES = FALSE)))),
              class = 'Date')
}

.registerHolidays <- function(type,dates){
    d <- unique(sort(dates))
    .holidays[[type]] <- data.frame(days=d,years=years(d))
}

registerHolidays <- function(type,dates){
    if (type %in% ls(envir=.holidays)) warning(paste('Overwriting',type,'holidays.'))

    d <- NULL
    if (inherits(dates,'Date')) d <- dates
    else if (is.character(dates)) d <- dateParse(dates)

    .registerHolidays(type,d)
}

addToHolidays <- function(type,dates){

    d <- NULL
    if (inherits(dates,'Date')) d <- dates
    else if (is.character(dates)) d <- dateParse(dates)

    if (!(type %in% ls(envir=.holidays)))
        warning(paste('No',type,'holidays exist. Registering.'))
    else
        d <- unique(sort(c(.holidays[[type]]$days,d)))

    .registerHolidays(type,d)
}

unregisterHolidays <- function(type,dates){
    if (!(type %in% ls(envir=.holidays))) {
        warning(paste('No',type,'holidays exist.'))
        return()
    }
    remove(list=type, envir=.holidays)
    invisible(NULL)
}

allHolidays <- function(){
    sort(c(ls(envir=.holidays), "NONE"))
}

isHoliday <- function(dates, type)
{
    if (!inherits(dates, 'Date'))
        if (is.character(dates))
            dates <- dateParse(dates)
        else
            dates <- as.Date(dates)

    if (inherits(type, 'Date')) {
        return(dates %in% type)
    } else {
        if (type == "NONE")
            return(rep(FALSE, length(dates)))
        else
            if (!(type %in% ls(envir=.holidays)))
                stop(paste('no', type, 'holidays exist.'))

        return(dates %in% .holidays[[type]]$days)
    }
}


