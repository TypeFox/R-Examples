#' @import RCurl
#' @import httr
#' @import jsonlite
library(RCurl)
library(httr)
library(jsonlite)

#' get current weather data for one location
#'
#' @param api_key Open weather map API key
#' @param cityID city ID
#' @param city name of city
#' @param country name of country
#' @param coordinates (lat,lon) coordinates of the location of your interest
#' @param zip_code zip code
#' @return data frame giving current weather data for one location
#' @export
#' @examples
#' data=get_current_weather(api_key,city="guwahati")

get_current_weather=function(api_key,cityID=NA,city="",country="",coordinates=NA,zip_code=NA)
{
  url="http://api.openweathermap.org/data/2.5/weather?"

  city=gsub(" ","+",city)
  country=gsub(" ","+",country)
  
  if(!is.na(cityID))
  {
    url=paste(url,"id=",cityID,sep="")
  }

  if(city != "")
  {
    url=paste(url,"q=",city,sep="")

    if(country!="")
    {
     
      url=paste(url,",",country,sep="")
    }
  }

  if(!is.na(coordinates[1]))
  {
    url=paste(url,"lat=",coordinates[1],"&lon=",coordinates[2],sep="")
  }

  if(!is.na(zip_code))
  {
    url=paste(url,"zip=",zip_code,",",country,sep="")
  }

  url=paste(url,"&APPID=",api_key,sep="")
  d=getURL(url)
  d=fromJSON(d)
  d
}

#' get current weather data for multiple cities
#'
#' @param api_key Open weather map API key
#' @param bbox bounding box [lat of the top left point, lon of the top left point, lat of the bottom right point, lon of the bottom right point, map zoom]
#' @param coordinates (lat,lon) coordinates of the location of your interest
#' @param count number of cities around the point that should be returned
#' @param cityIDs city IDs
#' @param cluster use server clustering of points. Possible values are [yes, no]
#' @param units metric unit
#' @return data frame giving current weather data for several locations
#' @export
#' @examples
#' data=get_multiple_cities(api_key,cityIDs =c(524901,703448,2643743))

get_multiple_cities=function(api_key,bbox=NA,coordinates=NA,count=NA,cityIDs=NA,cluster="yes",units="metric")
{
  url=""
  if(!is.na(bbox[1]))
  {
    url="http://api.openweathermap.org/data/2.5/box/city?"
    url=paste(url,"bbox=",bbox[1],",",bbox[2],",",bbox[3],",",bbox[4],",",bbox[5],"&cluster=",cluster,sep="")
  }

  if(!is.na(coordinates[1]))
  {
    url="http://api.openweathermap.org/data/2.5/find?"

    url=paste(url,"lat=",coordinates[1],"&lon=",coordinates[2],"&cnt=",count,"&cluster=",cluster,sep="")
  }

  if(!is.na(cityIDs[1]))
  {
    url="http://api.openweathermap.org/data/2.5/group?id="
    len=length(cityIDs)
    len=len-1
    for(i in 1:len)
    {
      url=paste(url,cityIDs[i],",",sep="")
    }
    url=paste(url,cityIDs[len+1],"&units=",units,sep="")
  }

  url=paste(url,"&APPID=",api_key,sep="")
  d=getURL(url)
  d=fromJSON(d)
  d
}

#' get weather forecast data for one location
#'
#' @param api_key Open weather map API key
#' @param cityID city ID
#' @param city name of city
#' @param country name of country
#' @param coordinates (lat,lon) coordinates of the location of your interest
#' @return data frame giving weather forecast data for one location
#' @export
#' @examples
#' data=get_weather_forecast(api_key,city="guwahati")

get_weather_forecast=function(api_key,cityID=NA,city="",country="",coordinates=NA)
{
  url="http://api.openweathermap.org/data/2.5/forecast?q="

  city=gsub(" ","+",city)
  country=gsub(" ","+",country)
  
  if(!is.na(cityID))
  {
    url=paste(url,"id=",cityID,sep="")
  }

  if(city != "")
  {
    url=paste(url,"q=",city,sep="")

    if(country!="")
    {
      url=paste(url,",",country,sep="")
    }
  }

  if(!is.na(coordinates[1]))
  {
    url=paste(url,"lat=",coordinates[1],"&lon=",coordinates[2],sep="")
  }

  url=paste(url,"&APPID=",api_key,sep="")
  d=getURL(url)
  d=fromJSON(d)
  d
}
