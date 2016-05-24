###################################################################
#
# This function is part of WACSgen V1.0
# Copyright © 2013,2014,2015, D. Allard, BioSP,
# and Ronan Trépos MIA-T, INRA
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. http://www.gnu.org
#
###################################################################



wacs.selectDates  = function(data, from, to)
{
  select = rep(TRUE,dim(data)[1])
  if (!is.null(from) | !is.null(to)){  
    for (i in 1:dim(data)[1]){
      currentdate = wacs.ymd2date(data$year[i],data$month[i],data$day[i])
      select[i] = (currentdate >= from) &  (currentdate <= to)
    }
  }
  return(select)
}
  

wacs.ymd2date = function(year,month,day){
  ddate = as.Date(paste(year,"/",month,"/",day,sep=""))
  return(ddate)
}


wacs.season = function(month ,day, seasons)
{
    for (s in 1:nrow(seasons)) {
        if (month < seasons$month[s]) {
            return (s);
        }
        if ((month == seasons$month[s]) && (day < seasons$day[s])) {
            return (s);
        }
    }
    return (1);
}
  
wacs.year = function(date)
{
    return(as.numeric(format.Date(date, '%Y')));
}

wacs.month = function(date)
{
    return(as.numeric(format.Date(date, '%m')));
}

wacs.day = function(date)
{
    return(as.numeric(format.Date(date, '%d')));
}

wacs.dayOfYear = function(date)
{
    return(as.numeric(strftime(date, format = "%j")));
}

wacs.leapYear = function(year)
{
    if (year %% 400 == 0) {
        return (TRUE);
    }
    if (year %% 100 == 0) {
        return (FALSE);
    }
    if (year %% 4 == 0) {
        return (TRUE);
    }
    return (FALSE);
}
