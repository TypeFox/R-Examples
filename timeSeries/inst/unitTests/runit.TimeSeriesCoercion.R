
# Rmetrics is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# Rmetrics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################


test.asTimeSeries =
function()
{
    # as.timeSeries.default - Returns the input
    # as.timeSeries.numeric - Transforms a numeric vector into a 'timeSeries'
    # as.timeSeries.data.frame - Transformas a 'data.frame' into a 'timeSeries'
    # as.timeSeries.matrix - Trasformas a 'matrix' into a 'timeSeries'
    # as.timeSeries.ts - Tranf orms a 'ts' object into a 'timeSeries'
    # as.timeSeries.character - Loads and transformas from a demo file
    # as.timeSeries.zoo - Transforms a 'zoo' object into a 'timeSeries'

    # Create timeSeries Object:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    data = round(rnorm(12), 3)
    charvec = timeCalendar(2006)
    uTS = timeSeries(data, charvec, units = "uTS")
    uTS
    checkTrue(inherits(uTS, "timeSeries"))
    checkTrue(is.timeSeries(uTS))

    # Check Positions:
    positions = timeCalendar()
    class(positions)
    whichFormat(format(positions))
    whichFormat(as.character(positions))

    # Data Input is a Vector - Returns a timeSeries with dummy positions:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rnorm(12)

    # as.numeric - add dummy dates:
    data = as.numeric(x)
    tS = as.timeSeries(data)
    head(tS)

    # as. numeric [as.vector] - add dummy dates:
    data = as.vector(x)
    tS = as.timeSeries(data)
    head(tS)

    # Data Inpiut is a data.frame:
    data(MSFT)
    x.df = as.data.frame(MSFT)
    head(x.df)
    # First Column holds Positions:
    tS = MSFT
    head(tS)
    # Missing Positions - return signal series
    # x.df = msft.dat[, -1]
    # head(x.df)
    # tS = as.timeSeries(x.df)
    # head(tS)

    # Data Input is a Matrix:
    data(MSFT)
    x.mat = as.matrix(MSFT)
    # tS = as.timeSeries(x.mat)
    # head(tS)                              # CHECK

    # Data Input is an Univariate/Muiltivariate timeSeries:
    x = MSFT
    class(x)
    tS = as.timeSeries(x)
    head(tS)

    # Note, data is a demo file ...
    tS = MSFT
    head(tS)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.asTimeSeriesDJ1 =
function()
{
    # Load Data:
    # use instead dummy data set just for testing ...
    Data = matrix(exp(cumsum(rnorm(30*100, sd = 0.1))), ncol = 30)
    Positions = format(timeSequence("2006-01-01", length.out = 100))
    DowJones30 = data.frame(Positions, Data)

    # Taking Dates from First Column:
    DJ = DowJones30[21:30, c(1, 11:15)]
    DJ
    class(DJ)
    as.timeSeries(DJ)

    # Adding Dates through Rownames Assignment:
    DJ = DowJones30[21:30, c(11:15)]
    rownames(DJ)<-DowJones30[21:30, 1]
    DJ
    as.timeSeries(DJ)

    # Missing Dates - Using Dummy Dates:
    DJ = DowJones30[21:30, c(11:15)]
    DJ
    class(DJ)
    as.timeSeries(DJ)

    # With recordIDs:
    if (FALSE) {
        DJ = DowJones30[21:30, c(1,11:15)]
        DJ = cbind(DJ, LETTERS[1:10])
        class(DJ)
        tsDJ = as.timeSeries(DJ)
        tsDJ
        tsDJ@recordIDs
    }

    DJ = DowJones30[21:30, c(11:15)]
    rownames(DJ) = DowJones30[21:30, 1]
    DJ = cbind(DJ, LETTERS[1:10])
    tsDJ = as.timeSeries(DJ)
    tsDJ
    tsDJ@recordIDs

    DJ = DowJones30[21:30, c(11:15)]
    DJ =cbind(DJ, LETTERS[1:10])
    tsDJ = as.timeSeries(DJ)
    tsDJ
    tsDJ@recordIDs

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.fromTimeSeriesUV =
function()
{
    if (FALSE) { 
    # DW has to be fixed ...
    
    # as.vector.timeSeries - Converts a univariate 'timeSeries' to a vector
    # as.matrix.timeSeries - Converts a 'timeSeries' to a 'matrix'
    # as.data.frame.timeSeries - Converts a 'timeSeries' to a 'data.frame'
    # as.ts.timeSeries - Converts a 'timeSeries' to a 'ts'

    # Univariate Case:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    data = round(rnorm(12), 3)
    charvec = timeCalendar(2006)
    uTS = timeSeries(data, charvec, units = "uTS")
    uTS

    # Vector:
    VEC = as.vector(uTS)
    head(VEC)
    class(VEC)
    checkIdentical(class(VEC), "numeric")

    # Numeric:
    # VEC = as.numeric(uTS)                         # Not implemented !
    # head(VEC)
    # class(VEC)
    # checkIdentical(class(VEC), "numeric")

    # Matrix:
    MAT = as.matrix(uTS)
    head(MAT)
    class(MAT)
    checkIdentical(class(MAT), "matrix")

    # Data Frame:
    DF = as.data.frame(uTS)
    head(DF)
    checkIdentical(class(DF), "data.frame")

    # Time Series:
    TS = as.ts(uTS)
    head(TS)
    class(TS)
    checkIdentical(class(TS), "ts")
    }

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.fromTimeSeriesMV =
function()
{
    if (FALSE) { 
    # DW has to be fixed ...
    
    # as.vector.timeSeries - Converts a univariate 'timeSeries' to a vector
    # as.matrix.timeSeries - Converts a 'timeSeries' to a 'matrix'
    # as.data.frame.timeSeries - Converts a 'timeSeries' to a 'data.frame'
    # as.ts.timeSeries - Converts a 'timeSeries' to a 'ts'

    # Multivariate Case:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    data = matrix(round(rnorm(24), 3), ncol = 2)
    charvec = timeCalendar(2006)
    mTS = timeSeries(data, charvec)
    mTS

    # Matrix:
    MAT = as.matrix(mTS)
    head(MAT)
    class(MAT)
    checkIdentical(
        target = class(MAT),
        current = "matrix")
    checkIdentical(
        target = as.vector(MAT[, 1]),
        current = as.numeric(MAT)[1:12])

    # Data Frame:
    DF = as.data.frame(mTS)
    head(DF)
    class(DF)
    checkIdentical(
        target = class(DF),
        current = "data.frame")

    # Time Series:
    TS = as.ts(mTS)
    head(TS)
    class(TS)
    checkIdentical(
        target = class(TS),
        current = c("mts", "ts"))
    }
    
    # Return Value:
    return()
}


################################################################################

