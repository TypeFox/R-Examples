decyear2date <- 
function(x) {
	yr <- floor(x)
	len <- ifelse(leapYear(yr), 366, 365)
	julday <- floor((x - yr) * len)
	as.Date(julday, origin = as.Date(paste(yr, 1, 1, sep = '-')))
}
