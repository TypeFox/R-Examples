holiday.fixed <- function( years, month, day )
{
  # function to return a fixed-date holiday that appears always
  # on the same month and day of the given year(s)
  len <- length( years )
  timeCalendar( m = rep( month[1], len), d = rep( day[1], len ), y = years )
}

holidays <- function( years, type = "USFederal", move = FALSE )
{
  # function that finds all the holidays in the given years of
  # the given types (by calling the holiday.xxx functions, where
  # xxx is one of the character strings in the type vector).
  #
  # If move is TRUE, it moves the holidays to the nearest weekday;
  # alternatively, move can be given as a vector the same length
  # as type, which will move the corresponding holidays only.
  #
  # The result will be a single sorted time vector object

  if( length( move ) == 1 )
    move <- rep( move, length( type ))

  # calculate the holidays

  if( length( type ) >= 1 )
    hollist <- lapply( 1:length( type ),
		 function( i, type, years, move )
		 {
		   x <- do.call( paste( "holiday.", type[[i]], sep = "" ),
		                 args = list( years = years ))
		   if( move[[i]] ) x <- holiday.nearest.weekday( x )
		   x
		 }, type, years, move )
  else
    hollist = list()

  # flatten the list and return it
  sort( unique( do.call( "c", hollist )))
}

holiday.weekday.number <- function( years, month, weekday, index )
{
# This is a function for calculating holidays that fall on e.g. the
# 3rd Monday of January.  The month is given as the month argument,
# the day of the week (Monday in this example, which is 1) as the
# weekday argument, and n (3 in this example) as the index argument.
# To get the last one in the month, pass -1 for index.
#
# month: 1-12  weekday 0-6, for Sunday - Saturday; index: -1 or 1 to 5
# years: e.g. 1997:2007
#
# Function will omit dates that don't exist, e.g. if you ask
# for the 5th Friday in March and there isn't one that year

  ret <- .Call( "time_from_month_day_index", as.integer(month),
               as.integer(weekday), as.integer(index), as.integer(years))
  ret <- timeZoneConvert(ret, as(timeDateOptions("time.zone")[[1]], "character"))
  ret@format <- timeDateFormatChoose(ret@columns[[2]], ret@time.zone)
  ret[ !is.na( ret ) ]
}

holiday.nearest.weekday <- function( dates. )
{
# Function that moves input dates falling on Saturday to the day
# before, and Sundays to the day after
  wkdays <- .Call( "time_to_weekday", as(dates., "timeDate"),
                  timeZoneList())
  addto <- -1 * ( wkdays == 6 ) + ( wkdays == 0 )
  dates. + addto
}


holiday.NewYears <- function( years )
{
  # return a date vector of January 1st
  holiday.fixed( years, month = 1, day = 1 )
}

holiday.MLK <- function( years )
{
# 3rd Monday in January, since 1986, according to
# http://www.usis.usemb.se/Holidays/celebrate/intro.htm
# and verified by other sites
   holiday.weekday.number( years, month = 1, weekday = 1, index = 3 )
}

holiday.Presidents <- function( years )
{
# 3rd Monday in February, since 1971, according to
# http://www.usis.usemb.se/Holidays/celebrate/intro.htm
# and verified by other sites
   holiday.weekday.number( years, month = 2, weekday = 1, index = 3 )
}


holiday.Easter <- function( years )
{
  # calculate when easter is
  ret <- .Call( "time_easter", as.integer(years))
  ret <- timeZoneConvert(ret, as(timeDateOptions("time.zone")[[1]], "character"))
  ret@format <- timeDateFormatChoose(ret@columns[[2]], ret@time.zone)
  ret
}

holiday.GoodFriday <- function(years) holiday.Easter(years) - 2

holiday.Memorial <- function( years )
{
# last Monday in May, since 1971, according to
# http://www.usis.usemb.se/Holidays/celebrate/intro.htm
# and verified by other sites
   holiday.weekday.number( years, month = 5, weekday = 1, index = -1 )
}


holiday.Independence <- function( years )
{
  # return a date vector of July 4th
  holiday.fixed( years, month = 7, day = 4 )
}

holiday.Labor <- function( years )
{
# First Monday in September, since 1882, according to
# http://www.usis.usemb.se/Holidays/celebrate/intro.htm
# and verified by other sites
   holiday.weekday.number( years, month = 9, weekday = 1, index = 1 )
}

holiday.Columbus <- function( years )
{
# 2nd Monday in October, since 1971, according to
# http://www.usis.usemb.se/Holidays/celebrate/intro.htm
# and verified by other sites
   holiday.weekday.number( years, month = 10, weekday = 1, index = 2 )
}

holiday.Veterans <- function( years )
{
  # return a date vector of November 11th
  holiday.fixed( years, month = 11, day = 11 )
}

holiday.Remembrance <- holiday.Veterans

holiday.Thanksgiving <- function( years )
{
# 4th Thursday in November, according to
# http://www.usis.usemb.se/Holidays/celebrate/intro.htm
# and verified by other sites
   holiday.weekday.number( years, month = 11, weekday = 4, index = 4 )
}

holiday.Christmas <- function( years )
{
  # return a date vector of Christmas
  holiday.fixed( years, month = 12, day = 25 )
}

holiday.USFederal <- function( years )
{
# New Years, MLK, Presidents, Memorial, Independence, Labor, Columbus,
# Veterans, Thanksgiving, Christmas, according to
# http://www.usis.usemb.se/Holidays/celebrate/intro.htm
# and verified by other sites
  holidays( years, c( "NewYears", "MLK", "Presidents", "Memorial",
	       "Independence", "Labor", "Columbus", "Veterans",
	       "Thanksgiving", "Christmas" ), move = TRUE )
}

holiday.NYSE <- function( years )
{
# Full-day holidays for the New York Stock Exchange since 1885, using
# information from their web site, www.nyse.com, under Data Library.
  years <- years[years >= 1885]
  # Every year they have observed New Years, Independence Day,
  # Thanksgiving and Christmas
  ret <- holidays( years, c("NewYears", "Independence", "Thanksgiving",
			    "Christmas"))
  # Labor day since 1887
  ret <- c( ret, holiday.Labor( years[ years >= 1887 ]))
  # Good Friday except 1898, 1906, 1907
  ret <- c(ret,
		holiday.GoodFriday(years[years != 1898 & years != 1906 &
					 years != 1907]))
  # Columbus Day on October 12, 1909 - 1953
  ret <- c(ret,
		holiday.fixed(years[ years >= 1909 & years <= 1953], 10, 12))
  # Martin Luther King Day since 1998
  ret <- c(ret, holiday.MLK(years[ years >= 1998 ]))
  # Lincoln's Birthday on Feb 12, 1896 - 1953
  ret <- c(ret,
		holiday.fixed(years[ years >= 1896 & years <= 1953], 2, 12))
  # Washington's Birthday on Feb 22 through 1970, and President's Day 1971-
  ret <- c(ret,
		holiday.fixed(years[ years <= 1970], 2, 22),
		holiday.Presidents( years[ years >= 1971]))
  # Decoration/Memorial Day May 30 through 1970, Memorial 1971-
  ret <- c(ret,
		holiday.fixed(years[ years <= 1970], 5, 30),
		holiday.Memorial( years[ years >= 1971]))
  # Veterans Day 1918, 1921, 1934-1953
  ret <- c(ret,
		holiday.Veterans(years[(years == 1918) | (years == 1921) |
					((years >= 1934) & (years <= 1953))]))
  # Election Day, Tuesday after 1st Mon in November, thru 1968 and then
  # presidential year election years only 1972 to 1980
  ret <- c(ret,
		holiday.weekday.number(years[years <= 1968 | years == 1972 |
					     years == 1976 | years == 1980],
				       11, 1, 1) + 1)
  # Holiday list is complete. Now Sunday holidays move to Monday.
  wkdays <- .Call( "time_to_weekday", ret,
		  timeZoneList())
  sun <- wkdays == 0
  ret[sun] <- ret[sun] + 1
  # after July 3, 1959, Saturday holidays move to Friday except if
  # at the end of monthly/yearly accounting period (last biz day of a month)
  sat <- wkdays == 6 & ret >= timeCalendar(d=3, m=7, y=1959) &
     (ret + timeRelative( "-1day +a1mth -1wkd" )) != (ret - 1)
  ret[sat] <- ret[sat] -1
  # take out remaining weekend dates, sort, and return
  wkdays <- .Call( "time_to_weekday", ret, timeZoneList())
  unique(sort(ret[ wkdays != 0 & wkdays != 6]))
}

holiday.Anzac <- function( years )
{
# 25th of April according to
# http://www.adfa.oz.au/~awm/anzacday/traditio.htm
   holiday.fixed( years, month = 4, day = 25 )
}

holiday.Australia <- function( years )
{
# January 26th according to
# http://www.effect.net.au/cuddlyk/myworld/history/ausday/ausday.html
  holiday.fixed( years, month = 1, day = 26 )
}

holiday.May <- function( years )
{
# May 1st, May Day, International labour day
  holiday.fixed( years, month = 5, day = 1 )
}

holiday.VE <- function( years )
{
# May 8th, V-E day, according to
# http://shoga.wwa.com/~android7/holidays.htm
  holiday.fixed( years, month = 5, day = 8 )
}


holiday.Canada <- function( years )
{
# July 1st according to
# http://fas.sfu.ca/canheritage/homepage/canid_hp/theme.htm
  holiday.fixed( years, month = 7, day = 1 )
}

holiday.Bastille <- function( years )
{
# July 14th, French Bastille day
  holiday.fixed( years, month = 7, day = 14 )
}

holiday.AllSaints <- function( years )
{
# November 1st, All Saints Day
  holiday.fixed( years, month = 11, day = 1 )
}

holiday.Thanksgiving.Canada <- function( years )
{
# 2nd Monday in October, according to
# http://www.smiley.cy.net/bdecie/Canada.html
  holiday.weekday.number( years, month = 10, weekday = 1, index = 2 )
}


holiday.Victoria <- function( years )
{
# Monday on or preceeding the 24th, Canadian holiday, according to
# http://pages.citenet.net/users/ctmx1108/webcalend/web-cal-top.html
# Can be calculated as 3 days after the 3rd Friday of the month
   holiday.weekday.number( years, month = 5, weekday = 5, index = 3 ) + 3
}

holiday.StPatricks <- function( years )
{
# March 17th
  holiday.fixed( years, month = 3, day = 17 )
}


