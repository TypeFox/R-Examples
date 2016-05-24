NCEP.aggregate <-
function(wx.data, YEARS=TRUE, MONTHS=TRUE, DAYS=TRUE, HOURS=TRUE, fxn='%>0'){

## 'wx.data' is a three dimentional array as returned by the function gather.wx.data()

## YEARS, MONTHS, DAYS, and HOURS should be marked TRUE if they are to be retained or FALSE if to be aggregated

## 'fxn' can be one of internal R functions such as 'mean', 'median', 'sd', or '%>0' to calculate the percentage of observations greater than zero ##



##################################
## If the fxn to perform is '%>0', create the function to calculate the percent of observations greater than zero ##

if(fxn == '%>0') {
fxn <- function(x){
	
	return(length(which(x > 0))/(length(x) - sum(is.na(x))))
}}


####################################################################
## First determine if any aggregation has been previously applied ##
## Use semi-arbitrary values to replace missing datetime elements ##
## We choose 2008 for the year because it is a leap year and can account
## for instances of February 29.  We choose October because it has 31
## days and can account for the longest months.
orig.dt <- strsplit(dimnames(wx.data)[[3]][1], split="_")[[1]]
if(orig.dt[1] == "XXXX") { substr(dimnames(wx.data)[[3]], start=1, stop=4) <- "2008" 
	YEARS = FALSE }
if(orig.dt[2] == "XX") { substr(dimnames(wx.data)[[3]], start=6, stop=7) <- "10"
	MONTHS = FALSE }
if(orig.dt[3] == "XX") { substr(dimnames(wx.data)[[3]], start=9, stop=10) <- "01"
	DAYS = FALSE }
if(orig.dt[4] == "XX") { substr(dimnames(wx.data)[[3]], start=12, stop=13) <- "01"
	HOURS = FALSE }


## Create a vector of dates from the names used in the data matrix ##
dt <- strptime(dimnames(wx.data)[[3]], "%Y_%m_%d_%H", tz='UTC')

## Determine which unique time increments should be included ##
if(YEARS == TRUE) dt.years <- unique(format(dt, "%Y")) else dt.years <- NA
if(MONTHS == TRUE) dt.months <- unique(format(dt, "%m")) else dt.months <- NA
if(DAYS == TRUE) dt.days <- unique(format(dt, "%d")) else dt.days <- NA
if(HOURS == TRUE) dt.hours <- unique(format(dt, "%H")) else dt.hours <- NA


dt4ag <- expand.grid(list(dt.years, dt.months, dt.days, dt.hours))

dt4ag.f <- c()
for(i in 1:length(dt4ag[,1])){
dt4ag.f[i] <- paste(ifelse(YEARS == TRUE, as.character(dt4ag$Var1[i]), "XXXX"), ifelse(MONTHS == TRUE, as.character(dt4ag$Var2[i]), "XX"),
	ifelse(DAYS == TRUE, as.character(dt4ag$Var3[i]), "XX"), ifelse(HOURS == TRUE, as.character(dt4ag$Var4[i]), "XX"), sep="_")
}

## Sort these dates ##
dt4ag.f.s <- sort(dt4ag.f)

## Determine which group each observation belongs to for aggregating ##
grp <- c()
for(i in 1:length(dimnames(wx.data)[[3]])){
grp[i] <- which(dt4ag.f.s == 	format(dt[i], paste(ifelse(YEARS == TRUE, "%Y", "XXXX"), ifelse(MONTHS == TRUE, "%m", "XX"), ifelse(DAYS == TRUE, "%d", "XX"), ifelse(HOURS == TRUE, "%H", "XX"), sep='_')))
}


## Create an empty array for output ##
out.wx.data <- array(NA, dim=c(length(dimnames(wx.data)[[1]]), length(dimnames(wx.data)[[2]]), length(dt4ag.f.s[unique(grp)])), 
		dimnames=list(dimnames(wx.data)[[1]], dimnames(wx.data)[[2]], dt4ag.f.s[unique(grp)])) 



## Apply the aggrigation and store the output in the appropriate location in the ouput array ##
for(i in 1:length(dimnames(wx.data)[[1]])){
for(j in 1:length(dimnames(wx.data)[[2]])){
	out.wx.data[i,j,] <- as.numeric(aggregate(wx.data[i,j,], by=list(grp), FUN=fxn)$x)
}}


return(out.wx.data)
}

