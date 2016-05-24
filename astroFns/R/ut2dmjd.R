ut2dmjd <-
function(yr=2012, mo=1, dy=1, hr=0, mi=0, se=0) {
# Convert from UT day and time to decimal modified Julian date
#
# A. Harris, U. Maryland Astronomy, 4/18/2008, 1/8/2009, 6/12/2012
#
# Uses function ymd2jd to compute Julian date

    if (dy <= 0) stop('*** Day must be positive ***')
    if (mo < 1 | mo > 12) stop('*** Month must be within a year ***')

    ymd2jd(yr, mo, dy) + (hr + (mi + se/60)/60)/24 - 2400000.5

}

