if (base::getRversion() >= "2.15.1") {
  utils::globalVariables(c("long", "lat", "group", "admin1.map", "admin1.regions"))
}
#' An Administrative Level 1 map of every country in the world
#' 
#' "Administration Level 1" is the generic term for the largest subnational administrative unit of a country. This
#' unit has different names depending on the country: for example, "state" in the USA and 
#' "prefecture" in Japan. In this data.frame the country name is in the column 
#' "admin" and the admin1 region name is in the column "region". Rather than working with this 
#' object directly, consider using the helper functions listed below.
#'  
#' @references Taken from http://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/. 
#' This is version 3.0.0 of the map and is considered to be beta.
#' The wikipedia page on "Administrative division": http://en.wikipedia.org/wiki/Administrative_division
#' @docType data
#' @name admin1.map
#' @usage data(admin1.map)
#' @note This map is too large to efficiently render by itself with ggplot2. You should subset it by country before
#' attempting to render. Please see the helper functions.
#' @seealso \code{\link{admin1.regions}}, \code{\link{get_admin1_regions}}, \code{\link{admin1_map}} and \code{\link{get_admin1_map}}
NULL

#' Names of all (country, region) pairs on the admin1.map data.frame. Here "region" means "Administrative Level 1 Region".
#' @name admin1.regions
#' @usage data(admin1.regions)
#' @docType data
#' @examples
#' data(admin1.regions)
#' head(admin1.regions)
#' @seealso \code{\link{admin1.map}}, \code{\link{get_admin1_regions}}, \code{\link{admin1_map}} and \code{\link{get_admin1_map}}
NULL

#' Render an Administrative Level 1 map for a specified country
#' 
#' Uses the map ?admin1.map.
#' 
#' @param country.name The name of the country you want to render.  
#' 
#' @examples
#' 
#' \dontrun{
#' admin1_map("japan")
#' 
#' admin1_map("canada")
#' }
#' @importFrom ggplot2 ggplot aes geom_polygon ggtitle
#' @importFrom utils data
#' @export
#' @seealso \code{\link{admin1.map}}, \code{\link{admin1.regions}}, \code{\link{get_admin1_regions}}, and \code{\link{get_admin1_map}}
admin1_map = function(country.name)
{
  data(admin1.map, package="choroplethrAdmin1", envir=environment())
  stopifnot(country.name %in% unique(admin1.map$admin))
  
  country.map = admin1.map[admin1.map$admin==country.name,]
  title = paste0("Administrative Level 1 Map of ", country.name)
  ggplot(country.map, aes(long, lat, group=group)) + 
    geom_polygon() +
    ggtitle(title)
}

#' Get all countries on the admin1 map
#' 
#' Uses ?admin1.regions
#' @export
#' @importFrom utils data
#' @examples 
#' get_admin1_countries()
get_admin1_countries = function()
{
  data(admin1.regions, package="choroplethrAdmin1", envir=environment())
  sort(unique(admin1.regions$country))
}

#' Get all admin1 region names for a given country
#' 
#' @param country.name The name of the country you want the admin1 region names of.
#' @export
#' @importFrom utils data
#' @examples
#' get_admin1_regions("japan")
#' get_admin1_regions("canada")
#' @seealso \code{\link{admin1.map}}, \code{\link{admin1.regions}}, \code{\link{admin1_map}} and \code{\link{get_admin1_map}}
get_admin1_regions = function(country.name)
{
  data(admin1.regions, package="choroplethrAdmin1", envir=environment())
  stopifnot(country.name %in% unique(admin1.regions$country))
  
  admin1.regions[admin1.regions$country==country.name,]  
}

#' Get an admin1 map for a country
#' 
#' Uses ?admin1.map. See ?admin1.regions for how countries are spelled in this map.
#' 
#' @param country.name The name of the country you want the admin1 map for.
#' @export
#' @importFrom utils data
#' @examples
#' \dontrun{
#'  japan.map = get_admin1_map("japan")
#' 
#'  ggplot(japan.map, aes(long, lat, group=group)) + 
#'    geom_polygon() +
#'    ggtitle("An admin1 map of Japan")
#' }
#' @seealso \code{\link{admin1.map}}, \code{\link{admin1.regions}}, \code{\link{get_admin1_regions}} and \code{\link{admin1_map}}
get_admin1_map = function(country.name)
{
  data(admin1.regions, package="choroplethrAdmin1", envir=environment())
  stopifnot(country.name %in% unique(admin1.regions$country))
  
  data(admin1.map, package="choroplethrAdmin1", envir=environment())
  admin1.map[admin1.map$admin==country.name,]  
}