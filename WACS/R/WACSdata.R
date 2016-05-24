
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



#' Format data for WACS
#' 
#' \code{WACSdata} Builds a data structure compatible with WACS functions
#' 
#' @export
#' 
#' @param data    A dataframe containing series of values for each variable 
#' @param mapping The names of special variables: year, month, day, rain, 
#'                tmin and tmax. Eg. list(RR = "rain", Tmin = "tmin") [optional; default is NULL]
#' @param bounds  A list of lists indicating the bounds for some variables
#'                eg. list(rain=list(min=0, max=7)) [optional; default is NULL]
#'                If not provided is set automatically according to data
#' @param from    Date at which the estimation should begin [optional; default is NULL]
#'                
#' @param to      Date at which the estimation should stop [optional; default is NULL]
#'                
#' @param skip    Vector of column names to skip[optional; default is NULL]
#' 
#' @param Trange Boolean value. When 
#'  \code{Trange=TRUE}, the couple \code{(tmin, trange=tmax-tmin)} is modeled. When
#'  \code{Trange=FALSE}, he couple \code{(tmin, tmax)} is modeled. Default is \code{Trange=FALSE}
#'  
#' @param seasons Vector of string of format 'mm-dd', gives the dates
#'                of change of seasons 
#'                (default: is c("03-01", "06-01", "09-01","12-01"))
#' @return A data frame structure, which will be used to call \link{WACSestim}, the function that estimates the parameters of the statistical model.
#' 
#' @author D. Allard, BioSP, Ronan Trépos MIA-T, INRA
#' 
#' @note
#' 
#' \code{bounds} can be provided as a list, as shown above.  If \code{bounds=NULL}, bounds are computed from the data. Some variables will have minimal values set automatically to 0
#' (\code{trange,V,RG,ETPP}) and maximal values to 100 (\code{ETPP}). Other minimum (resp. maximum) values are computed by adding (resp. subtracting) to 
#' the maximum (resp. minimum value) its difference to the 10th largest (resp. lowest) value.
#' 
#' 
#' \code{from} and \code{to} must be provided with format 'yyyy-mm-dd' (e.g. '2012-01-30').
#' 
#' There can be as many seasons as desired, with unequal length. There can also be one single season, in which case a single date is entered. 
#' 
#' Default is \code{seasons = c("03-01","06-01","09-01","12-01")}.
#' 
#' 
#' @examples
#' \dontrun{
#'   ## Simple example
#'   ThisData = WACSdata(ClimateSeries,from="1995-01-01",to="1999-12-31",
#'                       Trange=F,seasons=c("03-01","06-01","09-01","12-01"))
#'  }
#'  
 

WACSdata=function(data, 
                  mapping=NULL, 
                  bounds=NULL, 
                  from = NULL, 
                  to = NULL, 
                  skip=NULL, 
                  Trange=FALSE, 
                  seasons = c("03-01", "06-01", "09-01","12-01"))
{
    ## skip useless columns
    for (toskip in skip){
        data[[toskip]] <- NULL
    }
    ## check variable names
    mapping = wacs.buildMapping(mapping_user=mapping, names(data),Trange);
    data = wacs.renameData(data, mapping)
    if (! all(data$tmin <= data$tmax)) {
        stop ("[WACSdata] found tmin > tmax");
    }
    if (! all(data$rain >= 0)) {
        stop ("[WACSdata] found rain < 0");
    }
    ## set seasons
    l = length(seasons);
    seasonsRes = data.frame(month=rep(0,l), day=rep(0,l))
    for (s in 1:l) {
        tmp = strsplit(seasons[s],"-")[[1]]
        seasonsRes$month[s] = as.numeric(tmp[1]);
        seasonsRes$day[s] = as.numeric(tmp[2]);
    }
    seasons=seasonsRes;
    data$season = unlist(lapply(1:nrow(data), function(x){
       return(wacs.season(data$month[x], data$day[x], seasons));
    }))
    
    #reorder data
    data = wacs.reorderData(data,mapping,Trange);

    # Select period of time between 'from' and 'to'
    sel  = wacs.selectDates(data, from, to)
    data = data[sel,]
    
    #build default bounds
    bounds = wacs.getBounds(data, bounds, mapping);
    
    # remove 29 feb 
    data = data[-(which(data$month == 2 & data$day == 29)), ]

    # Removes very small amount of rain
    data$rain[data$rain <0.05] = 0
    
    res = list(data=data, mapping=mapping, 
               bounds = bounds, seasons=seasons, Trange=Trange)
    class(res) = "WACSdata"
    return(res)
}

