#' Make a tile grid map of Mexico
#'
#' @param estados A vector of Mexican state postal abbreviations. This vector must be
#'  some permutation of \code{mexico_abb}.
#' @param estados_colors A vector of "colors" in the R sense. For example strings
#'  (\code{"blue"}), hex codes (\code{"#D0C7B9"}), etc. The ith color in this
#'  vector will be the color of square that represents the ith element of
#'  \code{estados}.
#' @param border_colors Like \code{estados_colors} but specifying the border
#'  of the square.
#' @param estados_names Should the postal codes for each state be displayed in the
#'  center of the state? The default value is \code{TRUE}.
#' @param estados_name_colors Like \code{estados_colors} but specifying the color
#'  of the text displayed in each state.
#' @param estados_name_cex The size of the text displayed inside of each state.
#' @param font The font of the text displayed inside of each state. The values
#'  \code{"serif"}, \code{"sans"}, and \code{"mono"} are safest to use. Use
#'  other fonts at your own risk. If \code{NULL} a sans-style font will be used.
#'
#' @export
#' @examples
#' \dontrun{
#'  minimexico(mexico_abb, 1:32)
#' }
minimexico <- function(estados, estados_colors, border_colors = rep("white", 32),
                    estados_names = TRUE, estados_name_colors = rep("white", 32),
                    estados_name_cex = 1, font = NULL){

  # Make sure all parameters are specified
  if(any(unlist(lapply(list(
    estados, estados_colors, border_colors, estados_name_colors
  ), length)) != 32)){
    stop("Please make sure parameters contain exactly 32 elements.")
  }

  map <- data.frame(State = c("BCN", "SON", "CHH", "COA", "NLE",
                              "BCS", "SIN", "DUR", "ZAC", "SLP", "TAM",
                              "NAY", "AGU", "GUA", "QUE", "HID",
                              "JAL", "DIF", "TLA", "YUC",
                              "COL", "MOR", "VER", "TAB", "CAM", "ROO",
                              "MIC", "MEX", "PUE", "GRO", "OAX", "CHP"),
                    X = c(c(0, 2:5), 1:6, 2:6, c(3:5, 10), c(4:6, 9:11), 5:10),
                    Y = c(rep(5, 5), rep(4, 6), rep(3, 5), rep(2, 4), rep(1, 6),
                          rep(0, 6)),
                    stringsAsFactors = FALSE)

  # Make sure all states are present
  if(!all(map$State %in% estados)){
    stop("It appears some state names are repeated or missing.")
  }

  user_map <- data.frame(State = estados, scol = estados_colors, bcol =
                           border_colors, sncol = estados_name_colors,
                         stringsAsFactors = FALSE)

  map <- merge(map, user_map, by = "State")

  text_ <- NULL
  if(estados_names){
    text_ <- map$State
  }

  plotbox(map$X, map$Y, map$scol, map$bcol,
          n_xboxes = 12, n_yboxes = 6, text_ = text_,
          text_col = map$sncol, text_cex = estados_name_cex, text_font = font)
}
