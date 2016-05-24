#' Locate map, geocode data, then plot both.
#' 
#' createMap is a smart function that places data artifact on the map. If 
#' necessary it geocodes the data, locates map that fits all data artifacts, 
#' and plots the map with the data shapes sized and colored using metrics.
#' 
#' Geocoding:
#' If parameter \code{locationName} is missing then no geocoding is possible.
#' In that case parameters \code{lonName} and \code{latName} must contain 
#' names of columns with longitude and latitude information assigned to 
#' each data artifact (data point). 
#' If parameter \code{locationName} is defined then geocoding attempts 
#' to use values from the column with this name. Function \code{geocodeFun}
#' specifies geocoding function (with default \code{\link{geocode}}
#' from \link{ggmap} package). To speed up processing and avoid hitting 
#' global limit on Google Map API use memoised version of this function:
#' \code{memoise(geocode)} (see \code{\link[memoise]{memoise}}).
#' 
#' Map Locating:
#' Function operates in 2 modes: explicit map location mode and implicit mode.
#' In explicit mode value \code{location} locates the map using one
#' of two supported formats. If it is a 2-value vector then it contains a 
#' center of the map. If it is 4-value vector then it contains bounding box
#' coordinates: left/bottom/right/top.
#' In implicit mode, when \code{location} is missing, fuction uses parameters
#' \code{locator} and \code{data} to locate the map. If \code{locator} is
#' equal to \code{'center'} then it centers map by averaging longitude and
#' latitude values of all data artifacts. If \code{locator} is equal to \code{'box'}
#' then it determines min/max values of longitutude and latitude of all data 
#' artifacts and locates the map by corresponding bounding box.
#' Note that both modes support require explicit parameter \code{zoom} if 
#' applicable.
#' 
#' Map Types: variety of map avaiable are from several public sources: google, 
#' OpenStreetMap, Stamen, and CloudMade maps. The options and terms for each are 
#' different. For example, not all sources support both color and black-and-white options, 
#' or map types terrain, satellite, roadmap or hybrid.  
#' Note that in most cases by using Google source you are agreeing to the Google Maps API 
#' Terms of Service at https://developers.google.com/maps/terms.
#' 
#' Shapes: data artifacts are shapes placed over the map. Their size and fill are scaled using
#' values in \code{metrics} columns and their location is determined either by 
#' geocoding values from \code{locationName} column or with longitude and latitude values
#' stored in \code{lonName} and \code{latName} columns. 
#' 
#' Labels: If \code{labelName} is specified then column with such name contains text 
#' labels to place on the map (using the same locations as for the shapes).
#' 
#' @param data data frame with artifacts and their locations and metric(s) to be placed 
#'   on the map. If location name is provided (with \code{locationName}) then it is used to 
#'   gecode artifacts first. If not location then longitude and latitude must be provided.
#'   It is caller's responsibility adjust locations with value of \code{zoom} parameter to
#'   fit artifacts on the map.
#' @param maptype map theme as defined in \code{\link{get_map}}. options available are 'terrain', 
#'   'satellite', 'roadmap', and 'hybrid'
#' @param mapColor color (\code{'color'}) or black-and-white (\code{'bw'})
#' @param source Google Maps ('google'), OpenStreetMap ('osm'), Stamen Maps ('stamen'), or 
#'   CloudMade maps ('cloudmade')
#' @param location location of the map: longitude/latitude pair (in that order), or 
#'   left/bottom/right/top bounding box: 'center' uses 2 value vector for the center of the map, 
#'   while 'box' uses 4 value vector as left/bottom/right/top. If missing then function will 
#'   derive map location using parameter \code{locator} and the \code{data}.
#' @param locator in absence of \code{location} specifies how to use data to determine map 
#'   location: when 'center' then function averages out data point longitude and latitude values 
#'   to get approximate cneter for the map; when 'box' it will use min/max of longitude and 
#'   latitude values to determine bounding box: left/bottom/right/top. 
#'   If parameter \code{locationName} is specified then function will geocode values from this 
#'   column first. If paramter \code{locationName} is missing then it assumes that data is already 
#'   geocoded and stored in the columns with the names \code{lonName} and \code{latName}.
#' @param boxBorderMargin margin size in percent of box sizes to increase box when computed 
#'   from data locations.
#' @param zoom map zoom as defined in \code{\link{get_map}}: an integer from 3 (continent) 
#'   to 21 (building), default value 10 (city). Properly setting \code{zoom} for each map is 
#'   responsibility of a caller. Zoom is optional when using bounding box location specification. 
#' @param locationName vector of the column names with address or name to geocode its location (find 
#'   latitude and longitude) using \code{\link{geocode}} (see package \pkg{ggmap}). 
#'   When \code{locationName} is specified then parameters \code{lonName} and \code{latName} are ignored. 
#'   Multiple column names are used in order of appearance: geocoding tries 1st column's values first, 
#'   then, for the data points that didn't get resolved, it tries the 2d column's values, and so on.
#' @param lonName name of the column with longitude value. This value (in combination with value 
#'   from column \code{latName}) is used to place each data point on the map. This parameter is 
#'   ignored if \code{locationName} is defined. 
#' @param latName name of the column with latitude value. This value (in combination with value 
#'   from column \code{lonName}) is used to place each data point on the map. This parameter is 
#'   ignored if \code{locationName} is defined.
#' @param metricName (deprecated) Use parameter \code{metrics} instead.
#' @param metrics character vector of column names with metric values to scale shapes placed on map. First 
#'   metric corresponds to the size (or area depending on \code{scaleSize}), second to the fill gradient.
#'   See also \code{scaleSize} and \code{shapeStroke}.
#' @param scaleRange a numeric vector of lenght 2 that specifies the minimum and maximum size 
#'   of the plotting symbol after transformation (see parameter \code{range} of \code{\link{scale_size}}).
#' @param labelName name of the column to use for the artifact label text when displaying data. 
#' @param shape type of shape to use.
#' @param shapeColour color of metric artifacts placed on map.
#' @param shapeAlpha transparency of an artifact shape expressed as a fraction between 0 (complete 
#'   transparency) and 1 (complete opacity).
#' @param shapeStroke border width of an artifact shape. Remember, that in \code{ggplot2} \code{size} and \code{stroke}
#'   are additive so a point with \code{size = 5} and \code{stroke = 5} will have a diameter of 10mm. \code{\link{createMap}} 
#'   maps \code{metrics[[1]]} to shape size.
#' @param scaleSize logical if TRUE then scale artifact shapes by size (radius), otherwise scale shape's 
#'   area (artifact shapes scaling always uses \code{metrics[[1]]} values).
#' @param textColour color of artifact labels on map.
#' @param textFamily font family (when available) to use for artfiact labels.
#' @param textFace font style to apply to artifact labels: 'plain' (default), 'bold', 'italic', or 
#'   'bold.italic'.
#' @param textSize font size of artifact labels.
#' @param facet name of a column to divide plot into facets for specificed parameter (defualt is 
#'   NULL - no facets). If facet is single value then facet wrap applied (see \code{\link{facet_wrap}}), 
#'   otherwise facet grid (see \code{\link{facet_grid}} with 1st 2 values of the vector.
#' @param ncol number of facet columns (applies when single facet column supplied only - see 
#'   parameter \code{facet}).
#' @param facetScales Are scales shared across all facets: "fixed" - all are the same, "free_x" - vary 
#'   across rows (x axis), "free_y" - vary across columns (Y axis) (default), "free" - both rows and 
#'   columns (see in \code{facet_wrap} parameter \code{scales}).
#' @param geocodeFun geocode function. Default is \code{\link{geocode}} but due to Google API 
#'   restrictions use memoised version, e.g. \code{memoise(geocode)}, instead (see package \pkg{memoise}).
#' @param getmapFun get map function. Defayult is \code{\link{get_map}} but due to map APIs restrictions 
#'   use memoised version, e.g. \code{memose(get_map)}, instead (see package \pkg{memoise}).
#' @param urlonly return url only.
#' @param api_key an api key for cloudmade maps.
#' @param baseSize base font size.
#' @param baseFamily base font family.
#' @param title plot title.
#' @param legendPosition the position of metric guide ("left", "right", "bottom", "top", or two-element 
#'   numeric vector; "none" is no legend).
#' @param metricGuides list or vector with names of guide objects, or objects themselves, for up to 2 metrics. Typical 
#'   guides are \code{"legend"} or \code{"colorbar"} names and \code{\link[ggplot2]{guide_legend}} or
#'   \code{\link[ggplot2]{guide_colorbar}} objects.
#' @param defaultTheme plot theme to use, default is \code{theme_bw}.
#' @param themeExtra any additional \code{ggplot2} theme attributes to add.
#' @return a ggplot object
#' @export  
#' @examples
#' if(interactive()){
#' # initialize connection to Lahman baseball database in Aster 
#' conn = odbcDriverConnect(connection="driver={Aster ODBC Driver};
#'                          server=<dbhost>;port=2406;database=<dbname>;uid=<user>;pwd=<pw>")
#' 
#' data = computeAggregates(channel = conn, "teams_enh",
#'                 aggregates = c("min(name) name", "min(park) park", "avg(rank) rank", 
#'                                "avg(attendance) attendance"),
#'                 by = c("name || ', ' || park teamname", "lgid", "teamid", "decadeid"))
#'                
#' geocodeFun = memoise::memoise(ggmap::geocode)
#' getMapFun = memoise::memoise(ggmap::get_map)
#' 
#' createMap(data=data[data$decadeid>=2000,], 
#'           source = "stamen", maptype = "watercolor", zoom=4, 
#'           facet=c("lgid", "decadeid"),
#'           locationName=c('teamname','name'), 
#'           metrics=c('rank', 'attendance'), shape = 21,
#'           labelName='name', shapeColour="blue", scaleRange = c(2,12), textColour="black",
#'           title='Game Attendance by Decade and League (yearly, 2000-2012)',
#'           geocodeFun=geocodeFun, getmapFun = getMapFun)
#' }
createMap <- function(data,  
                      maptype = "terrain", 
                      mapColor = c("color", "bw"), 
                      source = c("google", "osm", "stamen", "cloudmade"),
                      location = NULL, locator = 'center', boxBorderMargin = 10,
                      zoom = NULL,
                      locationName = NULL,
                      lonName = "LONGITUDE", latName = "LATITUDE",
                      metricName = NULL, metrics = metricName, labelName = NULL, 
                      scaleRange = c(1,6),
                      shape = 19,
                      shapeColour = "red",
                      shapeAlpha = 0.5,
                      shapeStroke = 0.5,
                      scaleSize = TRUE,
                      textColour = "black", textFamily='mono' , textFace="plain", textSize=4,
                      facet = NULL, ncol = 1, facetScales = "fixed",
                      geocodeFun = memoise::memoise(geocode), getmapFun = get_map,
                      urlonly = FALSE, api_key = NULL,  
                      baseSize = 12, baseFamily = "sans", 
                      title = NULL,
                      legendPosition = "right",
                      metricGuides = c("legend", "colorbar"),
                      defaultTheme = theme_bw(base_size = baseSize),
                      themeExtra = NULL) {
  
  # match argument values
  maptype = match.arg(maptype, c('terrain', 'satellite', 'roadmap', 'hybrid', 'watercolor', 'toner'))
  mapColor = match.arg(mapColor, c('color', 'bw'))
  source = match.arg(source, c("google", "osm", "stamen", "cloudmade"))
  locator = match.arg(locator, c('center', 'box'))
  
  if (!is.null(metricName)) {
    toa_dep("0.4.1", "\"metricName\" argument in createMap is deprecated. Use \"metrics\" for column names with values used to scale shapes placed on map.")
  }
  
  if (length(metrics) > 2) {
    stop("createMap supports 2 or fewer metrics.")
  }
  
  if (!all(metrics %in% names(data))) {
    stop(paste("Some of the metrics", paste0("'", metrics, "'", collapse=", "), "are missing from the data."))
  }

  # geocode locations
  if (missing(location)) {
    if (!missing(locationName)) {
      geocodes = adply(data[ , locationName, drop=FALSE], 1, function(x) {
        x = as.data.frame(x, stringsAsFactors=FALSE)
        for(i in length(x)) {
          lookup_value = x[[1,i]]
          if (is.na(lookup_value) || is.null(lookup_value)) break
          latlon_lookup = geocodeFun(lookup_value, output="latlon")
          if (!any(is.na(latlon_lookup))) break
        }
        if(exists("latlon_lookup"))
          latlon_lookup
        else
          data.frame(lon=NA, lat=NA)
      })
      geocodes = geocodes[ , -which(names(geocodes) %in% locationName)]
    }else {
      geocodes = data[, c(lonName, latName)]
    }
    
    # calculate map location using data latitude and longitude
    if (locator == 'center') {
      # calculate center of the map
      ll = colwise(mean, na.rm = TRUE)(geocodes[, 1:2])
      location =c (ll[[1,1]], ll[[1,2]])
    }else if (locator == 'box') {
      # calculate min and max values
      mins = colwise(min, na.rm = TRUE)(geocodes[, 1:2])
      maxs = colwise(max, na.rm = TRUE)(geocodes[, 1:2])
      margin10percentLon = boxBorderMargin/100. * (maxs[1,1] - mins[1,1])
      margin10percentLat = boxBorderMargin/100. * (maxs[1,2] - mins[1,2])
      location = c(mins[1,1] - margin10percentLon, mins[1,2] - margin10percentLat, 
                      maxs[1,1] + margin10percentLon, maxs[1,2] + margin10percentLat)
    }
  }else {
    if (!typeof(location) %in% c('numeric','integer','double')) {
      stop("Parameter location is not numeric.")
    }
    
    if (!length(location) %in% c(2,4)) {
      stop("Length of parameter location must be 2 or 4.")
    }
  }
  
  # Set zoom if missing
  if (missing(zoom))
    if (source == 'google') zoom = 10
    else
      if (length(location) == 4 ) zoom = 'auto'
      else zoom = 10
  
  # normalize parameters for stamen
  if (source == 'stamen' && !(maptype %in% c('terrain', 'watercolor', 'toner'))) {
    warning(paste('Changed maptype for stamen to terrain instead of unsupported', maptype))
    maptype = 'terrain'
  }
  
  # Load map data
  m = getmapFun(location=location, zoom=zoom, scale=2, maptype=maptype, color=mapColor, source=source,
              urlonly = urlonly, api_key = api_key)
  
  if (urlonly) {
    return(m)
  }
  
  # Add geocodes to data
  if (!missing(locationName)) {
    data[, lonName] = geocodes$lon
    data[, latName] = geocodes$lat
    # remove data that didn't get geocoded successfully
    data = data[stats::complete.cases(data[,c(lonName,latName)]),]
  }
  
  # Create map with data
  p = ggmap(m) +
    labs(title=title)
  
  if (length(metrics) == 0) {
    p = p + geom_point(data=data, aes_string(x=lonName, y=latName, size=1), shape=shape, colour=shapeColour,
                       stroke=shapeStroke, alpha=shapeAlpha)
    
  }else if (length(metrics) == 1) {
    metric1 = metrics[[1]]
    p = p + geom_point(data=data, aes_string(x=lonName, y=latName, size=metric1), shape=shape, colour=shapeColour, 
               stroke=shapeStroke, alpha=shapeAlpha) +
    (if (scaleSize)
       scale_radius(metric1, range=scaleRange, guide=metricGuides[[1]])
    else
       scale_size(metric1, range=scaleRange, guide=metricGuides[[1]]))

  }else {
    metric1 = metrics[[1]]
    metric2 = metrics[[2]]
    p = p + geom_point(data=data, aes_string(x=lonName, y=latName, size=metric1, fill=metric2), shape=shape,
                               stroke=shapeStroke, alpha=shapeAlpha) +
      (if (scaleSize)
         scale_radius(metric1, range=scaleRange, guide=metricGuides[[1]])
       else
         scale_size(metric1, range=scaleRange, guide=metricGuides[[1]])) +
      scale_fill_gradient(metric2, low="grey", high=shapeColour, guide=metricGuides[[2]])
  }
  
  if (!missing(labelName)) {
    p = p +
      geom_text(data=data, aes_string(label=labelName, x=lonName, y=latName), colour=textColour, 
                family=textFamily , fontface=textFace, size=textSize, hjust=0.5, vjust=-0.5)
  }
  
  p = applyFacet(p, facet, facetScales, ncol)
  
  # apply themes
  p = p +
    theme(legend.position = legendPosition, 
          plot.title = element_text(family = baseFamily, face = "bold", size = baseSize * 1.4, vjust = 1),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  return(p)
}