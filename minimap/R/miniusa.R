#' Make a tile grid map of The United States of America
#'
#' @param states A vector of US state postal abbreviations. This vector must be
#'  some permutation of \code{usa_abb}.
#' @param state_colors A vector of "colors" in the R sense. For example strings
#'  (\code{"blue"}), hex codes (\code{"#D0C7B9"}), etc. The ith color in this
#'  vector will be the color of square that represents the ith element of
#'  \code{states}.
#' @param border_colors Like \code{state_colors} but specifying the border
#'  of the square.
#' @param state_names Should the postal codes for each state be displayed in the
#'  center of the state? The default value is \code{TRUE}.
#' @param state_name_colors Like \code{state_colors} but specifying the color
#'  of the text displayed in each state.
#' @param state_name_cex The size of the text displayed inside of each state.
#' @param font The font of the text displayed inside of each state. The values
#'  \code{"serif"}, \code{"sans"}, and \code{"mono"} are safest to use. Use
#'  other fonts at your own risk. If \code{NULL} a sans-style font will be used.
#'
#' @export
#' @examples
#' \dontrun{
#'  miniusa(state_abb, 1:51)
#' }
miniusa <- function(states, state_colors, border_colors = rep("white", 51),
                    state_names = TRUE, state_name_colors = rep("white", 51),
                    state_name_cex = 1, font = NULL){

  # Make sure all parameters are specified
  if(any(unlist(lapply(list(
    states, state_colors, border_colors, state_name_colors
  ), length)) != 51)){
    stop("Please make sure parameters contain exactly 51 elements.")
  }

  map <- data.frame(State = c('ME','AK','VT','NH','MA','WA','MT','ND','SD',
                              'MN','WI','MI','NY','CT','RI','OR','ID','WY',
                              'NE','IA','IL','IN','OH','PA','NJ','CA','NV',
                              'UT','CO','KS','MO','KY','WV','DC','MD','DE',
                              'AZ','NM','OK','AR','TN','VA','NC','HI','TX',
                              'LA','MS','AL','GA','SC','FL'),
                    X = c(11,0,9,10,11,1,2,3,4,5,6,7,9,10,11,1,2,3,4,5,6,7,8,
                          9,10,0,1,2,3,4,5,6,7,8,9,10,2,3,4,5,6,7,8,0,3,4,5,
                          6,7,8,7),
                    Y = c(7,6,6,6,6,5,5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4,4,3
                          ,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,1,1,1,1,1,1,1,0),
                    stringsAsFactors = FALSE)

  # Make sure all states are present
  if(!all(map$State %in% states)){
    stop("It appears some state names are repeated or missing.")
  }

  user_map <- data.frame(State = states, scol = state_colors, bcol =
                           border_colors, sncol = state_name_colors,
                         stringsAsFactors = FALSE)

  map <- merge(map, user_map, by = "State")

  text_ <- NULL
  if(state_names){
    text_ <- map$State
  }

  plotbox(map$X, map$Y, map$scol, map$bcol,
          n_xboxes = 12, n_yboxes = 8, text_ = text_,
          text_col = map$sncol, text_cex = state_name_cex, text_font = font)
}
