
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C) for this R-port:
#   1999 - 2012 Diethelm Wuertz, Zurich, <wuertz@itp.phys.ethz.ch>
#   2009 - 2012 Rmetrics Association, Zurich, www.rmetrics.org


################################################################################
# FUNCTION:             DESCRIPTION:
#  yahooImport           Downloads market data from chart.yahoo.com
#  yahooSeries           Easy to use download from chart.yahoo.com
################################################################################


yahooImport <-
    function (query, file = "tempfile", source = NULL,
    frequency = c("daily", "weekly", "monthly"),
    from = NULL, to = Sys.timeDate(), nDaysBack = 366,
    save = FALSE, sep = ";", try = TRUE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Downloads market data from Yahoo's web site

    # Notes:
    #   Requires: fields() cuts a string in fields
    #   Yahoo Token Description:
    #   s     Selected Ticker-Symbol
    #   a     First Quote starts with Month (mm): 0-11, Jan-Dec
    #   b     First Quote starts with Day (dd)
    #   c     First Quote starts with Year: as CCYY
    #   d     Last Quote ends with Month (mm): 0-11, Jan-Dec
    #   e     Last Quote ends with Day (dd)
    #   f     Last Quote ends with Year (yy): as CCYY
    #   r     Aggregation Level
    #   z     Selected Ticker-Symbol [optional]
    #   IBM SHARES, test 19/20 century change 01-12-1999 -- 31-01-2000:
    #   "s=IBM&a=11&b=1&c=1999&d=0&e=31&f=2000&g=d&x=.csv"

    # Examples:
    #   yahooImport("IBM", nDaysBack = 10)
    #   yahooImport(symbols = c("^DJI", "IBM"))
    #   yahooImport(symbols = c("^DJI", "IBM"), frequency = "weekly")
    #   yahooImport(frequency = "monthly", nDaysBack = 366)

    # FUNCTION:

    # Check:
    stopifnot(length(query) == 1)

    # Match Arguments:
    freq = match.arg(frequency)
    aggregation = substr(freq, 1, 1)

    # Automatic Selection of From / To:
    if (is.null(from) & is.null(to)) {
        to = Sys.timeDate()
        from = as.character(to - nDaysBack*24*3600)
        to = as.character(to) }

    # Extract Atoms - From:
    yearFrom = substring(from, 1, 4)
    monthFrom = as.character(as.integer(substring(from, 6, 7))-1)
    dayFrom = substring(from, 9, 10)

    # Extract Atoms - To:
    yearTo = substring(to, 1, 4)
    monthTo = as.character(as.integer(substring(to, 6, 7))-1)
    dayTo = substring(to, 9, 10)

    # Compose Query:
    Query = paste("s=", query, "&a=", monthFrom, "&b=", dayFrom,
        "&c=", yearFrom, "&d=", monthTo, "&e=", dayTo, "&f=", yearTo,
        "&g=", aggregation, "&x=.csv", sep = "")

    # Source:
    if (is.null(source))
        source = "http://chart.yahoo.com/table.csv?"

    # Download:
    if (try) {
        # First try if the Internet can be accessed:
        z = try(yahooImport(query, file, source, frequency, from, to,
            nDaysBack, save, sep, try = FALSE))
        if (inherits(z, "try-error") || inherits(z, "Error")) {
            return("No Internet Access")
        } else {
            return(z)
        }
    } else {
        # Download File:
        frequency = freq
        url = paste(source, Query, sep = "")
        tmp <- tempfile()
        download.file(url = url, destfile = tmp)

        # Read data and revert:
        X = as.timeSeries(read.table(tmp, header = TRUE, sep = ","))
    }

    # Save to file:
    if (save) {
        write.table(as.data.frame(X), file = file, sep = sep)
    } else {
        unlink(file)
    }

    # Result:
    ans = new("fWEBDATA",
        call = match.call(),
        param = c(
            "Instrument" = query,
            "Frequency " = frequency),
        data = X,
        title = "Data Import from www.yahoo.com",
        description = description() )

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


yahooSeries <-
    function(symbols, from = NULL, to = Sys.timeDate(), nDaysBack = 366, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Downloads easily time series data from Yahoo

    # Arguments:
    #   symbols - a character vector of symbol names
    #   from - from date
    #   to - to date
    #   nDaysBack - number of n-days back
    #   ... - arguments passed to the *Import()

    # Examples:
    #   yahooSeries("IBM", nDaysBack = 10)
    #   yahooSeries("IBM", frequency = "weekly")
    #   yahooSeries(c("^DJI", "IBM"))
    #   yahooSeries(c("^DJI", "IBM"), frequency = "monthly")

    # FUNCTION:

    # Download:
    X = yahooImport(query = symbols[1], ...)@data
    colnames(X) <- paste(symbols[1], colnames(X), sep = ".")
    N = length(symbols)
    if (N > 1) {
        for (i in 2:N) {
            Y = yahooImport(query = symbols[i], ...)@data
            colnames(Y) <- paste(symbols[i], colnames(Y), sep = ".")
            X = merge(X, Y)
        }
    }

    # Time Window:
    if (is.null(from)) from = to - nDaysBack*24*3600
    X = window(X, from, to)

    # Return Value:
    X
}


################################################################################

