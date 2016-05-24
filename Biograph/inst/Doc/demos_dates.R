# code to test the date functions
rm(list=ls())
library (Biograph)
data(GLHS)
data(NLOG98)

# convert cmc
z1 <- cmc_as_Date(GLHS$marriage)
z2 <- cmc_as_year (GLHS$marriage)

# Date to year
jj <- Date_as_year (z1)

# year as date
cc.dec <- year_as_Date (z2)

# Date as cmc
Date_as_cmc ("2012-08-30",format="%Y-%m-%d")

date_convert (d="01/01/2011",format.in="%m/%d/%Y",selectday=1,format.out="day-month-year",format.born="%m/%d/%Y")
date_convert (d="03/05/2011",format.in="%m/%d/%Y",selectday=1,format.out="age",born="03/05/1946",format.born="%m/%d/%Y")
date_convert (d="03/03/2011",format.in="%m/%d/%Y",selectday=1,format.out="age",born="03/05/1946",format.born="%m/%d/%Y")
date_convert (d=2010.333,format.in="year",format.out="age",born=1980.00,format.born="year")
date_convert (d=2010.330,format.in="year",format.out="%d-%B-%Y",format.born="year")
year_as_Date(x=2010.330,format.out="%d-%B-%Y")
date_convert (d=2010.330,format.in="year",format.out="age",born=1980.22,format.born="year")
date_convert (d=899,format.in="CMC",format.out="%d-%B-%Y",format.born="CMC")
date_convert(d=555,format.in="CMC",format.out="%d%b%y",format.born="CMC")
date_convert(d=555,format.in="CMC",format.out="%d-%B-%Y",format.born="CMC")
d<- "1946-03-01"
z<- date_convert(d,format.in="%Y-%m-%d",format.out="year",format.born="%Y-%m-%d") 
d <- "1Mar46"
z<- date_convert(d,format.in="%d%b%y",format.out="year",format.born="%d%b%y")  # requires chron ERROR 2046
d <- "1Mar1946"
z<- date_convert(d,format.in="%d%b%Y",format.out="year",format.born="%d%b%y") 
as.Date("1Mar1946",format="%d%b%Y")
d<- "01/March/1946"
z<- date_convert(d,format.in="%d/%B/%Y",format.out="year") 
d<- "01-March-1946"
z<- date_convert(d,format.in="%d-%B-%Y",format.out="year") 
a <- date_convert(d,format.in="%d-%B-%Y",format.out="age",born="15-April-1900") 
d="1March1950"
as.Date(d,format="%d%b%Y")
date_convert(d,format.in="%d%b%Y",format.out="year")
bb<- date_convert (d=GLHS$marriage,format.in="CMC",selectday=1,format.out="day-mon-year",born=Bdata$born)

# Convert entire Biograph object
GLHS.b <- date_b(GLHS,selectday=1,format.out="%d%b%Y",covs=c("marriage","LMentry"))
GLHS.yr <- date_b(GLHS,selectday=1,format.out="year",covs=c("marriage","LMentry"))
GLHS.a <- date_b(GLHS,format.out="age",covs=c("marriage","LMentry"))               



