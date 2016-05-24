#' Make map to visualize BISON data.
#'
#' @importFrom ggplot2 map_data ggplot aes geom_polygon coord_map scale_fill_gradient2 geom_path 
#' theme_bw labs scale_x_continuous scale_y_continuous guides guide_legend geom_point theme %+%
#' element_blank
#' @importFrom grid grid.newpage viewport unit
#' @importFrom sp point.in.polygon
#' @import mapproj
#' @export
#' @param input Input bison object.
#' @param tomap One of points (occurrences), county (counts by county), or state
#'    (counts by state).
#' @param geom geom_point or geom_jitter, not quoted.
#' @param jitter jitter position, see ggplot2 help.
#' @param customize Pass in more to the plot.
#' @return Map (using ggplot2 package) of points on a map.
#' @examples \dontrun{
#' # Using function bison
#' library("ggplot2")
#' out <- bison(species="Accipiter", type="scientific_name", count=300)
#' bisonmap(input=out)
#' bisonmap(input=out, geom=geom_jitter, jitter=position_jitter(width = 0.3, height = 0.3))
#'
#' # Using function bison_solr
#' out <- bison_solr(scientificName='Ursus americanus', rows=200)
#' bisonmap(out)
#' }

bisonmap <- function(input = NULL, tomap="points", geom = geom_point,
                     jitter = NULL, customize = NULL) {
  
  UseMethod("bisonmap")
}

#' @method bisonmap bison
#' @export
#' @rdname bisonmap
bisonmap.bison <- function(input = NULL, tomap="points", geom = geom_point,
                           jitter = NULL, customize = NULL) {
  
  long = lat = group = total = NULL

  if (!is.bison(input)) {
    stop("Input is not of class bison")
  }

  if (tomap == 'points') {
    bison_map_maker(x = input, geom = geom, jitter = jitter, customize = customize)
  } else
    if (tomap == 'county') {
      bycounty <- input$counties
      bycounty$state <- tolower(bycounty$state)
      bycounty$county_name <- gsub("\\scounty", "", tolower(bycounty$county_name))

      counties <- map_data("county")
      counties_plus <- merge(counties, bycounty, by.x = 'subregion', by.y = 'county_name', all.x = TRUE)
      counties_plus <- counties_plus[order(counties_plus$order),]
      counties_plus$total <- as.numeric(counties_plus$total)

      states <- map_data("state")

      ggplot(counties_plus, aes(long, lat, group = group)) +
        geom_polygon(aes(fill = total)) +
        coord_map(projection = "azequalarea") +
        scale_fill_gradient2("", na.value = "white", low = "white", high = "steelblue") +
        geom_path(data = counties, colour = "grey", size = .3, alpha = .4) +
        geom_path(data = states, colour = "grey", size = .4) +
        theme_bw(base_size = 14) +
        labs(x = "", y = "") +
        bison_blanktheme() +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        theme(legend.position = "top") +
        guides(guide_legend(direction = "horizontal")) +
        customize
    } else if (tomap == 'state') {
        bystate <- input$states
        bystate$record_id <- tolower(bystate$record_id)

        states <- map_data("state")
        states_plus <- merge(states, bystate, by.x = 'region', by.y = 'record_id', all.x = TRUE)
        states_plus <- states_plus[order(states_plus$order),]
        states_plus$total <- as.numeric(states_plus$total)

        ggplot(states_plus, aes(long, lat, group = group)) +
          geom_polygon(aes(fill = total)) +
          coord_map(projection = "azequalarea") +
          scale_fill_gradient2("", na.value = "white", low = "white", high = "steelblue") +
          geom_path(data = states, colour = "grey", size = .4) +
          theme_bw(base_size = 14) +
          labs(x = "", y = "") +
          bison_blanktheme() +
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0)) +
          theme(legend.position = "top") +
          guides(guide_legend(direction = "horizontal")) +
          customize
    } else { 
      stop("tomap must be one of points, county, or state") 
    }
}


#' @method bisonmap bison_solr
#' @export
#' @rdname bisonmap
bisonmap.bison_solr <- function(input = NULL, tomap="points", geom = geom_point,
                                jitter = NULL, customize = NULL) {
  
  if (!is.bison_solr(input)) {
    stop("Input is not of class bison_solr")
  }
  if (!tomap == 'points') {
    stop("tomap must equal 'points'")
  }
  # remove NA's due to /null,null/specimen in pointPath field
  input$points <- input$points[!input$points$pointPath == "/null,null/specimen",]
  bison_map_maker(x = input, geom = geom, jitter = jitter, customize = customize)
}
