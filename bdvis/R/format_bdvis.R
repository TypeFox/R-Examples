#' Prepare data frame for flagging functions
#'
#' \code{format_bdvis} renames certain fields in the data frame to make sure the 
#' other package functions knows how to use them. This step is highly recommended 
#' for the proper working of the functions.
#'
#' When invoked, there are three ways of indicating the function how to
#' transform the data.frame: using the \code{source} parameter, providing a
#' \code{config} object with field mapping, or passing individual values to the
#' mapping function. This is the order in which the function will parse
#' arguments; \code{source} overrides \code{config}, which overrides other
#' mapping arguments.
#'
#' \code{source} refers to the package that was used to retrieve the data.
#' Currently, three values are supported for this argument: "\code{rgbif}",
#' "\code{rvertnet}" and "\code{rinat}", but many more are on their way.
#'
#' \code{config} asks for a configuration object holding the mapping of the
#' field names. This option is basically a shortcut for those users with
#' custom-formatted data.frames who will use the same mapping many times, to
#' avoid having to type them each time. In practice, this object is a named list
#' with the following four fields: \code{Latitude}, \code{Longitude}, 
#' \code{Date_collected} and \code{Scientific_name}. Each element must have a 
#' string indicating the name of the column in the data.frame holding the values 
#' for that element. If the data.frame doesn't have one or more of these fields, 
#' put \code{NA} in that element; otherwise, the function will throw an error. 
#' See the examples section.
#'
#' If none of the two is provided, the function expects the user to provide the
#' mapping by passing the individual column names associated with the right term.
#' See the examples section.
#'
#' @param indf Required. The data.frame on which to operate.
#' @param source Optional. Indicates the package that was used to retrieve the
#'   data. Currently accepted values are "rvertnet", "rgbif" or "rinat". Either
#'   \code{source}, \code{config} or individual parameters must be present (see
#'   details).
#' @param config Optional. Configuration object indicating mapping of field
#'   names from the data.frame to the DarwinCore standard. Useful when importing
#'   data multiple times from a source not available via the \code{source}
#'   argument. Either \code{source}, \code{config} or individual parameters must
#'   be present (see details).
#' @param quiet Optional. Don't show any logging message at all. Defaults to
#'   FALSE.
#' @param gettaxo optional. Call function \link{gettaxo} to build higher level 
#'   taxanony. Defauls to FALSE.
#' @param ... Optional. If none of the previous is present, the four key
#'   arguments (\code{Latitude}, \code{Longitude},
#'   \code{Date_collected}, \code{Scientific_name}) can be put here. See examples.
#'
#' @return The provided data frame, with field names changed to suite the functioning
#' of further visulization functions.
#'
#'@examples \dontrun{
#' # Using the rgbif package and the source argument
#' if (requireNamespace("rinat", quietly=TRUE)) {
#'  d <- get_inat_obs_project("reptileindia") 
#'  d <- format_bdvis(d, source="rinat")
#'
#'  # Using a configuration object, matches 'rinat' schema
#'  conf <- list(Latitude="latitude",
#'               Longitude="longitude",
#'               Date_collected="Observed.on",
#'               Scientific_name="Scientific.name")
#'  d <- format_bdvis(d, config=conf)
#'
#'  # Passing individual parameters, all optional
#'  d <- format_bdvis(d,
#'                 Latitude="lat",
#'                 Longitude="lng",
#'                 Date_collected="ObservedOn",
#'                 Scientific_name="sciname")
#' }
#'}
#'
#' @export
format_bdvis <- function(indf, source=NULL, config=NULL, quiet=FALSE, gettaxo=F, ...) {
  
  # Parse input object type
  bd_check_df(indf)
  
  # Mapping via 'source'
  if (!(is.null(source))) {
    match.arg(source, sources_list)
    if (!(quiet)) message(c("Mapping according to ", source, " format"))
    new_fields = bd_get_source(source)
    # Mapping via 'config'
  } else if (!(is.null(config))) {
    if (!(quiet)) message("Mapping according to config object")
    new_fields = bd_parse_config(config)
    # Mapping via individual parameters
  } else {
    if (!(quiet)) message("Mapping via individual parameters")
    new_fields = bd_parse_args(list(...))
  }
  # Apply the transformation
  if (!(is.null(new_fields$Latitude)) && new_fields$Latitude != "Latitude") {
    if (new_fields$Latitude %in% names(indf)) {
      if ("Latitude" %in% names(indf)) {
        names(indf)[names(indf)=="Latitude"] <- "Latitude::original"
        if (!(quiet)) message("Changed \"Latitude\" to \"Latitude::original\"")
      }
      names(indf)[names(indf)==new_fields$Latitude] <- "Latitude"
      if (!(quiet)) message(c("Changed \"",new_fields$Latitude,"\" to \"Latitude\""))
    }
  }
  if (!(is.null(new_fields$Longitude)) && new_fields$Longitude != "Longitude") {
    if (new_fields$Longitude %in% names(indf)) {
      if ("Longitude" %in% names(indf)) {
        names(indf)[names(indf)=="Longitude"] <- "Longitude::original"
        if (!(quiet)) message("Changed \"Longitude\" to \"Longitude::original\"")
      }
      names(indf)[names(indf)==new_fields$Longitude] <- "Longitude"
      if (!(quiet)) message(c("Changed \"",new_fields$Longitude,"\" to \"Longitude\""))
    }
  }
  if (!(is.null(new_fields$Date_collected)) && new_fields$Date_collected != "Date_collected") {
    if (new_fields$Date_collected %in% names(indf)) {
      if ("Date_collected" %in% names(indf)) {
        names(indf)[names(indf)=="Date_collected"] <- "Date_collected::original"
        if (!(quiet)) message("Changed \"Date_collected\" to \"Date_collected::original\"")
      }
      names(indf)[names(indf)==new_fields$Date_collected] <- "Date_collected"
      if (!(quiet)) message(c("Changed \"",new_fields$Date_collected,"\" to \"Date_collected\""))
    }
  }
  if (!(is.null(new_fields$Scientific_name)) && new_fields$Scientific_name != "Scientific_name") {
    if (new_fields$Scientific_name %in% names(indf)) {
      if ("Scientific_name" %in% names(indf)) {
        names(indf)[names(indf)=="Scientific_name"] <- "Scientific_name::original"
        if (!(quiet)) message("Changed \"Scientific_name\" to \"Scientific_name::original\"")
      }
      names(indf)[names(indf)==new_fields$Scientific_name] <- "Scientific_name"
      if (!(quiet)) message(c("Changed \"",new_fields$Scientific_name,"\" to \"Scientific_name\""))
    }
  }
  indf$Latitude=as.numeric(indf$Latitude)
  indf$Longitude=as.numeric(indf$Longitude)
  indf=getcellid(indf)
  if(gettaxo){
    indf=gettaxo(indf)
  }
  
  indf
}

sources_list <- c(
  "rgbif",
  "rvertnet",
  "rinat"
)

bd_get_source <- function(source) {
  bd_sources <- list(
    rgbif=list(
      Latitude="decimalLatitude",
      Longitude="decimalLongitude",
      Date_collected="eventDate",
      Scientific_name="name"
    ),
    rvertnet=list(
      Latitude="decimallatitude",
      Longitude="decimallongitude",
      Date_collected="eventdate",
      Scientific_name="scientificname"
    ),
    rinat=list(
      Latitude="Latitude",
      Longitude="Longitude",
      Date_collected="Observed.on",
      Scientific_name="Scientific.name"
    )
  )
  return(bd_sources[[source]])
}

bd_parse_config <- function(config){
  if (!("Latitude" %in% names(config))) stop("\"Latitude\" missing from configuration object")
  if (!("Longitude" %in% names(config))) stop("\"Longitude\" missing from configuration object")
  if (!("Date_collected" %in% names(config))) stop("\"Date_collected\" missing from configuration object")
  if (!("Scientific_name" %in% names(config))) stop("\"Scientific_name\" missing from configuration object")
  return(config)
}

bd_parse_args <- function(args) {
  bd_args = list()
  if ("Latitude" %in% names(args)) {
    bd_args$Latitude=args$Latitude
  } else {
    bd_args$Latitude=NULL
  }
  if ("Longitude" %in% names(args)) {
    bd_args$Longitude=args$Longitude
  } else {
    bd_args$Longitude=NULL
  }
  if ("Date_collected" %in% names(args)) {
    bd_args$Date_collected=args$Date_collected
  } else {
    bd_args$Date_collected=NULL
  }
  if ("Scientific_name" %in% names(args)) {
    bd_args$Scientific_name=args$Scientific_name
  } else {
    bd_args$Scientific_name=NULL
  }
  return(bd_args)
}

bd_check_df <- function(indf) {
  if(is.na(indf) || (is.data.frame(indf) && nrow(indf) == 0)) stop("Input data frame missing or empty")
  if(!(is.data.frame(indf))) stop("Provided argument is not a data.frame")
  return(invisible())
}
