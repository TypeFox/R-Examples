## ----, message=FALSE-----------------------------------------------------
library(bizdays)
cal <- Calendar(holidays=holidaysANBIMA, name='ANBIMA', weekdays=c('saturday', 'sunday'), dib=252)
bizdays.options$set(default.calendar=cal)

## ------------------------------------------------------------------------
from_dates <- c('2013-07-12', '2012-06-13')
to_dates <- seq(as.Date('2014-02-17'), as.Date('2016-07-21'), by='months')
bizdays(from_dates, to_dates)

