#' @title Return Site Information.
#' @description Return site information from the Neotoma Paleoecological Database.
#'
#' \code{get_site} returns site information from the Neotoma Paleoecological Database
#'    based on parameters defined by the user.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr content GET
#' @param sitename character string representing the full or partial site name, or an object of class \code{dataset}, \code{dataset_list}, \code{download} or \code{download_list}
#' @param altmin Minimum site altitude  (in m).
#' @param altmax Maximum site altitude (in m).
#' @param loc A numeric vector c(lonW, latS, lonE, latN) representing the bounding box within which to search for sites.  The convention here is to use negative values for longitudes west of Grewnwich or longitudes south of the equator.
#' @param gpid A character string or numeric value, must correspond to a valid geopolitical identity in the Neotoma Database.  Use get.tables('GeoPoliticalUnits') for a list of acceptable values, or link here: http://api.neotomadb.org/apdx/geopol.htm
#' @param  ... Optional additional arugments
#'
#' @author Simon J. Goring \email{simon.j.goring@@gmail.com}
#' @return A data frame:
#'
#'  \item{\code{siteid}}{Unique database record identifier for the site.}
#'  \item{\code{sitename}}{Name of the site.}
#'  \item{\code{long}}{Mean longitude, in decimal degrees, for a site (-180 to 180).}
#'  \item{\code{lat}}{Mean latitude, in decimal degrees, for a site (-90 to 90).}
#'  \item{\code{elev}}{Elevation in meters.}
#'  \item{\code{description}}{Free form description of a site, including such information as physiography and vegetation around the site.}
#'  \item{\code{long_acc}}{If the site is described by a bounding box this is the box width.}
#'  \item{\code{lat_acc}}{If the site is described by a bounding box this is the box height.}
#'
#' @examples \dontrun{
#' #  What is the distribution of site elevations in Neotoma?
#' all.sites <- get_site()  #takes a bit of time.
#'
#' plot(density(all.sites$elev, from = 0, na.rm=TRUE),
#' main = 'Altitudinal Distribution of Neotoma Sites', xlab = 'Altitude (m)', log='x')
#'
#' #  Get site information from a dataset:
#' nw.datasets <- get_dataset(loc = c(-140, 50, -110, 65), datasettype='pollen',taxonname='Pinus*')
#' nw.sites <- get_site(nw.datasets)
#'
#' }
#' @references
#' Neotoma Project Website: http://www.neotomadb.org
#' API Reference:  http://api.neotomadb.org/doc/resources/sites
#' @keywords IO connection
#' @export
get_site <- function(sitename, altmin, altmax, loc, gpid, ...) {
  UseMethod('get_site')
  
}

#' @title Return Site Information.
#' @description Return site information from the Neotoma Paleoecological Database.
#' @importFrom jsonlite fromJSON
#' @importFrom httr content GET

#' @param sitename A character string representing the full or partial site name.
#' @param ... Arguments passed from the generic method, not used.
#' 
#' @export
get_site.default <- function(sitename, ...) {

  base.uri <- 'http://api.neotomadb.org/v1/data/sites'

  cl <- as.list(match.call())
  
  cl[[1]] <- NULL
  cl <- lapply(cl, eval, envir = parent.frame())

  #  Pass the parameters to param_check to make sure everything is kosher.
  error_test <- param_check(cl)
  if (error_test[[2]]$flag == 1) {
    stop(paste0(unlist(error_test[[2]]$message), collapse='\n  '))
  } else{
    cl <- error_test[[1]]
  }
  
  neotoma_content <- content(GET(base.uri, query = cl), as = "text")
  
  if (identical(neotoma_content, "")) stop("")
  
  aa <- jsonlite::fromJSON(neotoma_content, simplifyVector = FALSE)
  
  if (aa[[1]] == 0) {
    stop(paste('Server returned an error message:\n', aa[[2]]), call. = FALSE)
  }
  if (aa[[1]] == 1) {
    aa <- aa[[2]]
    
    rep_NULL <- function(x) { 
      if (is.null(x)) {NA}
      else{
        if (class(x) == 'list') {
          lapply(x, rep_NULL)
        } else {
          return(x)
        }
      }
    }
    
    aainsta <- rep_NULL(aa)
    
    if (length(aa) == 0) {
      cat('The API call was successful, but no records were returned.\n')
      return()
    }
    
    cat('The API call was successful, you have returned ',
        length(aa), 'records.\n')
    
  }

  if (class(aa) == 'try-error') {
     output <- aa
  } else {
    
    # replace NULL values:
    aa <- lapply(aa, function(x) ifelse(x == "NULL", NA, x))
    
    output <- data.frame(site.id = sapply(aa, '[[', 'SiteID'),
                         site.name = sapply(aa, '[[', 'SiteName'),
                         long = rowMeans(data.frame(sapply(aa, '[[', 'LongitudeWest'), 
                                                    sapply(aa, '[[', 'LongitudeEast')),
                                         na.rm = TRUE),
                         lat = rowMeans(data.frame(sapply(aa, '[[', 'LatitudeNorth'),
                                                   sapply(aa, '[[', 'LatitudeSouth')),
                                        na.rm = TRUE),
                         elev = sapply(aa, '[[', 'Altitude'),
                         description = sapply(aa, '[[', 'SiteDescription'),
                         long.acc = abs(sapply(aa, '[[', 'LongitudeWest') - 
                                          sapply(aa, '[[', 'LongitudeEast')),
                         lat.acc = abs(sapply(aa, '[[', 'LatitudeNorth') - 
                                         sapply(aa, '[[', 'LatitudeSouth')),
                         stringsAsFactors = FALSE)
    
  }

  class(output) <- c('site', 'data.frame')
  output

}


#' @title Return Site Information from a \code{dataset}
#' @description Return site information from the Neotoma Paleoecological Database.
#'
#' @param sitename An object of class \code{dataset}.
#' @param ... Arguments passed from the generic method, not used.
#' @export
get_site.dataset <- function(sitename, ...) {
  site <- sitename$site.data
  class(site) <- c('site', 'data.frame')
  site
}

#' @title Return Site Information from a \code{dataset_list}
#' @description Return site information from the Neotoma Paleoecological Database.
#'
#' @param sitename An object of class \code{dataset_list}.
#' @param ... Arguments passed from the generic method, not used.
#' @export
get_site.dataset_list <- function(sitename, ...) {
  site <- do.call(rbind.data.frame,lapply(sitename, '[[', 'site.data'))
  class(site) <- c('site', 'data.frame')
  site
}

#' @title Return Site Information from a \code{download}
#' @description Return site information from the Neotoma Paleoecological Database.
#'
#' @param sitename An object of class \code{download}.
#' @param ... Arguments passed from the generic method, not used.
#' @export
get_site.download <- function(sitename, ...) {

  site <- sitename$dataset$site.data
  
  class(site) <- c('site', 'data.frame')
  site
}

#' @title Return Site Information from a \code{download_list}
#' @description Return site information from the Neotoma Paleoecological Database.
#'
#' @param sitename An object of class \code{download_list}.
#' @param ... Arguments passed from the generic method, not used.
#' @export
get_site.download_list <- function(sitename, ...) {
  
  site <- do.call(rbind.data.frame,lapply(lapply(sitename, '[[', 'dataset'), '[[', 'site.data'))
  
  class(site) <- c('site', 'data.frame')
  site
}

#' @title Return Site Information from a \code{geochronologic}
#' @description Return site information from the Neotoma Paleoecological Database.
#'
#' @param sitename An object of class \code{geochronologic}.
#' @param ... Arguments passed from the generic method, not used.
#' @export
get_site.geochronologic <- function(sitename, ...) {
  
  site <- sitename[[1]]$site.data
  
  class(site) <- c('site', 'data.frame')
  site
}

#' @title Return Site Information from a \code{geochronologic_list}
#' @description Return site information from the Neotoma Paleoecological Database.
#'
#' @param sitename An object of class \code{geochronologic_list}.
#' @param ... Arguments passed from the generic method, not used.
#' @export
get_site.geochronologic_list <- function(sitename, ...) {
  
  site <- do.call(rbind.data.frame,lapply(sitename, function(y)y[[1]]$site.data))
  
  class(site) <- c('site', 'data.frame')
  site
}
