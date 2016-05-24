# Create a new date object
date(); Sys.Date(); Sys.time()
date1 <- as.Date(x = c("1/15/2013", "10/29/2014"), format = "%m/%d/%Y")
date2 <- as.Date(x = c("15JAN13", "29OCT14"), format = "%d%b%y")
date3 <- seq(from = as.Date("2013-1-15"), by = "days", length = 3)
c(class(date1), class(date3)); date1; date3

# Create a new date/time object
my.time <- c("5/16/2016 09:45:04", "12/31/2018 18:45:04")
time1 <- as.POSIXct(my.time, format = "%m/%d/%Y %H:%M:%S")
time2 <- as.POSIXlt(my.time, format = "%m/%d/%Y %H:%M:%S")
time3 <- ISOdate(year = c(2014, 2016), month = c(5, 12),
  day = c(16, 30), hour = 9)                              # POSIXct only
time4 <- strptime(my.time, format = "%m/%d/%Y %H:%M:%S")  # POSIXlt only
class(time1); class(time2); class(time3); class(time4)
attributes(time2); names(time2); attr(time2, "names")[1]
time1; time2; time3; time4
unclass(time1) 
unclass(time2)

# Five ways of extracting a date/time component
a1 <- weekdays(date1); a2 <- months(time1); a3 <- quarters(time2)
b1 <- strftime(date1, format = "%Y")
b2 <- strftime(time2, format = "myDate = %m/%d/%Y; myTime = %H:%M:%S")
b3 <- strftime(time3, format = "%H")
c1 <- as.POSIXlt(time1)$year + 1900  # need POSIXlt class
c2 <- time2$year + 1900
d1 <- format(date1, format = "%b")
e1 <- substr(as.character(time1), start = 1, stop = 4)
c(class(a1), class(b1), class(c1), class(d1), class(e1))
a1; b1; c1; d1; e1

# Operation on date/time objects
date1[2] - date2[1]; date1 + 3; date1[1]  < date1[2]
mean(time1); range(date3)
dif <- difftime(date2, date1, units = "days"); dif; str(dif)

# Application: results sensitive to time zone or daylight saving time
x <- c("2005-04-02 19:03:00", "2005-04-03 02:00:00",
  "2005-04-03 14:25:00")
y <- strptime(x, format = "%Y-%m-%d %H:%M:%S")
z <- strptime(x, format = "%Y-%m-%d %H:%M:%S", tz = "EST")
dif.y <- difftime(y[2], y[1], units = "mins")  # Doest not work
dif.z <- difftime(z[2], z[1], units = "mins")  # Yes, work!
y; is.na(y); dif.y
z; is.na(z); dif.z