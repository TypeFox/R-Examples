
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
#  fredImport            Downloads market data from research.stlouisfed.org
#  fredSeries            Easy to use download from research.stlouisfed.org
################################################################################


fredImport <-
    function(query, file = "tempfile", source = NULL,
    frequency = "daily",
    from = NULL, to = Sys.timeDate(), nDaysBack = NULL,
    save = FALSE, sep = ";", try = TRUE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Downloads Monthly Market Data, Indices and Benchmarks from
    #   St. Louis FED, "research.stlouisfed.org".

    # Value:
    #   An One Column data frame with row names denoting the dates
    #   given in the POSIX format "%Y%m%d". The column lists the
    #   downloaded data records.

    # Examples:
    #   fredImport("DPRIME")

    # Notes:
    #   This function is written for one-column daily data sets.
    #   Some example data sets include:
    #     DEXUSEU   U.S. / Euro Foreign Exchange Rate
    #     DEXSZUS   Switzerland / U.S. Foreign Exchange Rate
    #     DGS1      1-Year Treasury Constant Maturity Rate
    #     DPRIME    Bank Prime Loan Rate

    # FUNCTION:

    # Settings:
    stopifnot(length(query) == 1)

    # Source"
    if (is.null(source))
        source = "http://research.stlouisfed.org/fred2/series/"

    # Check:
    if (frequency != "daily")
        stop("Only daily data records are supported!")

    # Download:
    if (try) {
        # Try for Internet Connection:
        z = try(fredImport(query, file, source, frequency, from, to,
            nDaysBack, save, sep, try = FALSE))
        if (inherits(z, "try-error") || inherits(z, "Error")) {
            return("No Internet Access")
        } else {
            return(z)
        }
    } else {
        # Download File:
        queryFile = paste(query, "/downloaddata/", query, ".txt", sep = "")
        url = paste(source, queryFile, sep = "")
        tmp = tempfile()
        download.file(url = url, destfile = tmp)

        # Scan the file:
        x1 = scan(tmp, what = "", sep = "\n")

        # Extract dates ^19XX and ^20XX:
        x2 = x1[regexpr("^[12][90]", x1) > 0]
        x1 = x2[regexpr(" .$", x2) < 0]

        # Compose Time Series:
        data = matrix(
            as.numeric(substring(x1, 11, 999)), byrow = TRUE, ncol = 1)
        charvec = substring(x1, 1, 10)
        X = timeSeries(data, charvec, units = query)
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
        title = "Data Import from research.stlouisfed.org",
        description = description() )

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


fredSeries <-
    function(symbols, from = NULL, to = Sys.timeDate(), nDaysBack = 366, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Downloads easily time series data from St. Louis FRED

    # Arguments:
    #   symbols - a character vector of symbol names
    #   from - from date
    #   to - to date
    #   nDaysBack - number of n-days back
    #   ... - arguments passed to the *Import()

    # Examples:
    #   fredSeries("DPRIME")[1:10, ]

    # FUNCTION:

    # Download:
    X = fredImport(query = symbols[1], ...)@data
    N = length(symbols)
    if (N > 1) {
        for (i in 2:N) {
            X = merge(X, fredImport(query = symbols[i], ...)@data)
        }
    }

    # Time Window:
    if (is.null(from)) from = to - nDaysBack*24*3600
    X = window(X, from, to)

    # Return Value:
    X
}


################################################################################


