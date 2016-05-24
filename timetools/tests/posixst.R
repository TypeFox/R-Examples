#===========================================================
# The following code test if the POSIXst constructors are ok
# first serie of tests check the ranges of definition
library( timetools )

# all theses must raise an error
try( POSIXst(62, 'second', 'minute') )
try( POSIXst(3602, 'second', 'hour') )
try( POSIXst(86402, 'second', 'day') )
try( POSIXst(604802, 'second', 'week') )
try( POSIXst(2678402, 'second', 'month') )
try( POSIXst(31622402, 'second', 'year') )

try( POSIXst(60, 'minute', 'hour') )
try( POSIXst(1440, 'minute', 'day') )
try( POSIXst(10080, 'minute', 'week') )
try( POSIXst(44640, 'minute', 'month') )
try( POSIXst(527040, 'minute', 'year') )

try( POSIXst(24, 'hour', 'day') )
try( POSIXst(168, 'hour', 'week') )
try( POSIXst(744, 'hour', 'month') )
try( POSIXst(8784, 'hour', 'year') )

try( POSIXst(7, 'day', 'week') )
try( POSIXst(0, 'day', 'month') )
try( POSIXst(32, 'day', 'month') )
try( POSIXst(366, 'day', 'year') )

try( POSIXst(1, 'week', 'month') )
try( POSIXst(1, 'week', 'year') )

try( POSIXst(12, 'month') )

# all theses must succeed
POSIXst(c(0, 61), 'second', 'minute')
POSIXst(c(0, 3601), 'second', 'hour')
POSIXst(c(0, 86401), 'second', 'day')
POSIXst(c(0, 604801), 'second', 'week')
POSIXst(c(0, 2678401), 'second', 'month')
POSIXst(c(0, 31622401), 'second', 'year')

POSIXst(c(0, 59), 'minute', 'hour')
POSIXst(c(0, 1439), 'minute', 'day')
POSIXst(c(0, 10079), 'minute', 'week')
POSIXst(c(0, 44639), 'minute', 'month')
POSIXst(c(0, 527039), 'minute', 'year')

POSIXst(c(0, 23), 'hour', 'day')
POSIXst(c(0, 167), 'hour', 'week')
POSIXst(c(0, 743), 'hour', 'month')
POSIXst(c(0, 8783), 'hour', 'year')

format(POSIXst(c(0, 6), 'day', 'week'), '%s %v of %m')
POSIXst(c(1, 31), 'day', 'month')
POSIXst(c(0, 365), 'day', 'year')

format(POSIXst(c(0, 11), 'month'), '%s %v of %m')

#==============================================
# the following tests check posixst 'extractor'

cet <- as.POSIXct('2012-09-08 02:00', 'CET')
utc <- as.POSIXct('2012-09-08', 'UTC')

# in UTC
POSIXst(utc, 'second', 'minute')
POSIXst(utc, 'second', 'hour')
POSIXst(utc, 'second', 'day')
POSIXst(utc, 'second', 'week')
POSIXst(utc, 'second', 'month')
POSIXst(utc, 'second', 'year')

POSIXst(utc, 'minute', 'hour')
POSIXst(utc, 'minute', 'day')
POSIXst(utc, 'minute', 'week')
POSIXst(utc, 'minute', 'month')
POSIXst(utc, 'minute', 'year')

POSIXst(utc, 'hour', 'day')
POSIXst(utc, 'hour', 'week')
POSIXst(utc, 'hour', 'month')
POSIXst(utc, 'hour', 'year')

format(POSIXst(utc, 'day', 'week'), '%s %v of %m')
POSIXst(utc, 'day', 'month')
POSIXst(utc, 'day', 'year')

format(POSIXst(utc, 'month'), '%s %v of %m')
POSIXst(utc, 'year')

# in CET
POSIXst(cet, 'second', 'minute')
POSIXst(cet, 'second', 'hour')
POSIXst(cet, 'second', 'day')
POSIXst(cet, 'second', 'week')
POSIXst(cet, 'second', 'month')
POSIXst(cet, 'second', 'year')

POSIXst(cet, 'minute', 'hour')
POSIXst(cet, 'minute', 'day')
POSIXst(cet, 'minute', 'week')
POSIXst(cet, 'minute', 'month')
POSIXst(cet, 'minute', 'year')

POSIXst(cet, 'hour', 'day')
POSIXst(cet, 'hour', 'week')
POSIXst(cet, 'hour', 'month')
POSIXst(cet, 'hour', 'year')

format(POSIXst(cet, 'day', 'week'), '%s %v of %m')
POSIXst(cet, 'day', 'month')
POSIXst(cet, 'day', 'year')

format(POSIXst(cet, 'month'), '%s %v of %m')
POSIXst(cet, 'year')

# in UTC from CET. All must be TRUE
all.equal( POSIXst(cet, 'second', 'minute', 'UTC'), POSIXst(utc, 'second', 'minute') )
all.equal( POSIXst(cet, 'second', 'hour', 'UTC'), POSIXst(utc, 'second', 'hour') )
all.equal( POSIXst(cet, 'second', 'day', 'UTC'), POSIXst(utc, 'second', 'day') )
all.equal( POSIXst(cet, 'second', 'week', 'UTC'), POSIXst(utc, 'second', 'week') )
all.equal( POSIXst(cet, 'second', 'month', 'UTC'), POSIXst(utc, 'second', 'month') )
all.equal( POSIXst(cet, 'second', 'year', 'UTC'), POSIXst(utc, 'second', 'year') )

all.equal( POSIXst(cet, 'minute', 'hour', 'UTC'), POSIXst(utc, 'minute', 'hour') )
all.equal( POSIXst(cet, 'minute', 'day', 'UTC'), POSIXst(utc, 'minute', 'day') )
all.equal( POSIXst(cet, 'minute', 'week', 'UTC'), POSIXst(utc, 'minute', 'week') )
all.equal( POSIXst(cet, 'minute', 'month', 'UTC'), POSIXst(utc, 'minute', 'month') )
all.equal( POSIXst(cet, 'minute', 'year', 'UTC'), POSIXst(utc, 'minute', 'year') )

all.equal( POSIXst(cet, 'hour', 'day', 'UTC'), POSIXst(utc, 'hour', 'day') )
all.equal( POSIXst(cet, 'hour', 'week', 'UTC'), POSIXst(utc, 'hour', 'week') )
all.equal( POSIXst(cet, 'hour', 'month', 'UTC'), POSIXst(utc, 'hour', 'month') )
all.equal( POSIXst(cet, 'hour', 'year', 'UTC'), POSIXst(utc, 'hour', 'year') )

all.equal( POSIXst(cet, 'day', 'week', 'UTC'), POSIXst(utc, 'day', 'week') )
all.equal( POSIXst(cet, 'day', 'month', 'UTC'), POSIXst(utc, 'day', 'month') )
all.equal( POSIXst(cet, 'day', 'year', 'UTC'), POSIXst(utc, 'day', 'year') )

all.equal( POSIXst(cet, 'month', tz='UTC'), POSIXst(utc, 'month') )
all.equal( POSIXst(cet, 'year', tz='UTC'), POSIXst(utc, 'year') )

# in CET from UTC. All must be TRUE
all.equal( POSIXst(utc, 'second', 'minute', 'CET'), POSIXst(cet, 'second', 'minute') )
all.equal( POSIXst(utc, 'second', 'hour', 'CET'), POSIXst(cet, 'second', 'hour') )
all.equal( POSIXst(utc, 'second', 'day', 'CET'), POSIXst(cet, 'second', 'day') )
all.equal( POSIXst(utc, 'second', 'week', 'CET'), POSIXst(cet, 'second', 'week') )
all.equal( POSIXst(utc, 'second', 'month', 'CET'), POSIXst(cet, 'second', 'month') )
all.equal( POSIXst(utc, 'second', 'year', 'CET'), POSIXst(cet, 'second', 'year') )

all.equal( POSIXst(utc, 'minute', 'hour', 'CET'), POSIXst(cet, 'minute', 'hour') )
all.equal( POSIXst(utc, 'minute', 'day', 'CET'), POSIXst(cet, 'minute', 'day') )
all.equal( POSIXst(utc, 'minute', 'week', 'CET'), POSIXst(cet, 'minute', 'week') )
all.equal( POSIXst(utc, 'minute', 'month', 'CET'), POSIXst(cet, 'minute', 'month') )
all.equal( POSIXst(utc, 'minute', 'year', 'CET'), POSIXst(cet, 'minute', 'year') )

all.equal( POSIXst(utc, 'hour', 'day', 'CET'), POSIXst(cet, 'hour', 'day') )
all.equal( POSIXst(utc, 'hour', 'week', 'CET'), POSIXst(cet, 'hour', 'week') )
all.equal( POSIXst(utc, 'hour', 'month', 'CET'), POSIXst(cet, 'hour', 'month') )
all.equal( POSIXst(utc, 'hour', 'year', 'CET'), POSIXst(cet, 'hour', 'year') )

all.equal( POSIXst(utc, 'day', 'week', 'CET'), POSIXst(cet, 'day', 'week') )
all.equal( POSIXst(utc, 'day', 'month', 'CET'), POSIXst(cet, 'day', 'month') )
all.equal( POSIXst(utc, 'day', 'year', 'CET'), POSIXst(cet, 'day', 'year') )

all.equal( POSIXst(utc, 'month', tz='CET'), POSIXst(cet, 'month') )
all.equal( POSIXst(utc, 'year', tz='CET'), POSIXst(cet, 'year') )

#===================================
# test for Time*DataFrame extraction

ti <- RegularTimeInstantDataFrame('2012-01-01', '2012-01-01 02:23', 'second')
POSIXst(ti, 'hour', 'day', tz='CET')
POSIXst(ti, 'hour', 'day')

POSIXstti <- RegularTimeIntervalDataFrame('2012-01-01', '2012-02-03', 'day')
POSIXst(ti, 'hour', 'day', cursor=0.5)

ti <- RegularTimeIntervalDataFrame('2012-01-01', '2012-01-03', 'day')
POSIXst(ti, 'minute', 'day')
lapply(POSIXst(ti, 'day', 'week'), format, '%s %v of %m')

#========================
# test for the Ops method

one <- second(utc, 'minute', tz='CET')
two <- second(utc+POSIXctp(4, 'second'), 'minute', tz='CET')
three <- POSIXst(cet, 'second', 'minute')

# ==

c(one, two, three) == c(one, two, three)
three == one
c(three, two) == c(one, two)
c(one, two) == c(two, three)
c(one, two, two) == c(one, two, three)
c(one, two) == c(one, two, three)

# !=

c(one, two, three) != c(one, two, three)
three != one
c(three, two) != c(one, two)
c(one, two) != c(two, three)
c(one, two, two) != c(one, two, three)
c(one, two) != c(one, two, three)

# <

c(one, two, three) < c(one, two, three)
three < one
c(three, two) < c(one, two)
c(one, two) < c(two, three)
c(one, two, two) < c(one, two, three)
c(one, two) < c(one, two, three)

# <=

c(one, two, three) <= c(one, two, three)
three <= one
c(three, two) <= c(one, two)
c(one, two) <= c(two, three)
c(one, two, two) <= c(one, two, three)
c(one, two) <= c(one, two, three)

# >

c(one, two, three) > c(one, two, three)
three > one
c(three, two) > c(one, two)
c(one, two) > c(two, three)
c(one, two, two) > c(one, two, three)
c(one, two) > c(one, two, three)

# >=

c(one, two, three) >= c(one, two, three)
three >= one
c(three, two) >= c(one, two)
c(one, two) >= c(two, three)
c(one, two, two) >= c(one, two, three)
c(one, two) >= c(one, two, three)

