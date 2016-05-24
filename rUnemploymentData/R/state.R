if (base::getRversion() >= "2.15.1") {
  utils::globalVariables(c("df_state_unemployment", "state.regions"))
}

#' Get state unemployment data from the US Bureau of Labor Statistics (BLS) website
#' 
#' This function is included to allow you to verify the integrity of ?df_state_unemployment.
#' @param year A year (integer) between 2000 and 2013
#' @importFrom rvest html html_nodes html_text 
#' @export
get_state_unemployment_df = function(year=2013)
{
  stopifnot(is.numeric(year))
  stopifnot(year >= 2000 && year <= 2013) # at the time of this writing, the only dates available
  
  # BLS uses final 2 years for string - URL is formulaic
  year_string = as.character(year)
  year_string = substr(year_string, 3, 4)
  url = paste0(
          "http://www.bls.gov/lau/lastrk",
          year_string,
          ".htm")
  data = rvest::html(url)
  
  # the data is in a table where the states are in ascending order of value
  # the first element is always the average, which we don't care about.
  
  # get state names
  state_names = character(0)
  data_state_names = html_nodes(data, "#bodytext td:nth-child(2)")
  for (i in 2:length(data_state_names)) # 1 is always "UNITED STATES", which we don't want
  {
    name = tolower(html_text(data_state_names[[i]]))
    state_names = c(state_names, name)
  }
  state_names
  
  # get values
  values = numeric(0)
  data_values = html_nodes(data, "td~ td+ td")
  for (i in 2:length(data_values)) # 1 is always the average, which we don't want
  {
    value = as.numeric(html_text(data_values[[i]]))
    values = c(values, value)
  }
  values
  
  data.frame(region=state_names, value=values)
}

#' Build the data object ?df_state_unemployment
#' 
#' This function is included to allow you to verify the integrity of ?df_state_unemployment.
#' This will scrape the Bureau of Labor Statistics Website to get the data.
#' @export
build_state_df = function()
{
  data(state.regions, package="choroplethrMaps", envir=environment())
  df_state_unemployment = data.frame(region=state.regions$region)
  for (year in 2000:2013)
  {
    df = get_state_unemployment_df(year)
    colnames(df) = c("region", eval(year))
    df_state_unemployment = merge(df_state_unemployment, df, all=TRUE)
  }
  df_state_unemployment
}

#' Render Choropleth Map of US State Unemployment Rates
#' 
#' Data comes from ?df_state_unemployment. The choropleth is rendered with the function
#' ?state_choropleth in the choroplethr package.
#' 
#'  @param year The year of data to use. Must be between 2000 and 2013.
#'  @param buckets The number of equally sized buckets to places the values in. 
#'  A value of 1 will use a continuous scale, and a value in [2, 9] will use that many buckets.
#'  @param zoom An optional vector of states to zoom in on. Elements of this vector 
#'  must exactly match the names of states as they appear in the "region" column of 
#'  ?state.regions in the choroplethrMaps package.
#' @importFrom choroplethr state_choropleth
#' @export
state_unemployment_choropleth = function(year = 2013, buckets = 7, zoom = NULL)
{
  # validate input
  stopifnot(is.numeric(year))
  stopifnot(year >= 2000 && year <= 2013)
  
  # get data into right format for choroplethr
  data(df_state_unemployment, package="rUnemploymentData", envir=environment())
  df = df_state_unemployment[, c("region", eval(year))]
  colnames(df) = c("region", "value")  
  
  # sensible defaults
  title  = paste0("State Unemployment Rates: Year ", year)
  legend = "Unemployment Rate"
  state_choropleth(df, title, legend, buckets, zoom)
}

#' Create an animated choropleth of US State Unemployment Data
#' 
#' Data comes from ?df_state_unemployment. The choropleth is rendered with the function
#' ?state_choropleth in the choroplethr package. Note that this command will write files to 
#' your local file system - see ?choroplethr_animate in the choroplethr package for details.
#' @export
#' @importFrom choroplethr state_choropleth choroplethr_animate
animated_state_unemployment_choropleth = function()
{
  data(df_state_unemployment, package="rUnemploymentData", envir=environment())
  
  frames = list()
  for (i in 2:ncol(df_state_unemployment))
  {
    year = colnames(df_state_unemployment)[i]
    df = df_state_unemployment[,c(1,i)]
    colnames(df) = c("region", "value")
    frames[[i-1]] = state_choropleth(df, 
                                     title=paste0("State Unemployment Map: ", year),
                                     legend = "Unemployment Rate")
  }
  choroplethr_animate(frames)
}