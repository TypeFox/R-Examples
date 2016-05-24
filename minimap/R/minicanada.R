#' Make a tile grid map of Canada
#'
#' @param pt A vector of Canadian province and territory postal abbreviations.
#' This vector must be some permutation of \code{canada_abb}.
#' @param pt_colors A vector of "colors" in the R sense. For example strings
#'  (\code{"blue"}), hex codes (\code{"#D0C7B9"}), etc. The ith color in this
#'  vector will be the color of square that represents the ith element of
#'  \code{pt}.
#' @param border_colors Like \code{pt_colors} but specifying the border
#'  of the square.
#' @param pt_names Should the postal codes for each province or territory be
#' displayed in the center of the province or territory? The default value is
#' \code{TRUE}.
#' @param pt_name_colors Like \code{pt_colors} but specifying the color
#'  of the text displayed in each province or territory.
#' @param pt_name_cex The size of the text displayed inside of each province or
#'  territory.
#' @param font The font of the text displayed inside of each province or
#'  territory. The values \code{"serif"}, \code{"sans"}, and \code{"mono"} are
#'  safest to use. Use other fonts at your own risk. If \code{NULL} a
#'  sans-style font will be used.
#'
#' @export
#' @examples
#' \dontrun{
#'  minicanada(canada_abb, 1:13)
#' }
minicanada <- function(pt, pt_colors, border_colors = rep("white", 13),
                    pt_names = TRUE, pt_name_colors = rep("white", 13),
                    pt_name_cex = 1, font = NULL){
  # Make sure all parameters are specified
  if(any(unlist(lapply(list(
    pt, pt_colors, border_colors, pt_name_colors
  ), length)) != 13)){
    stop("Please make sure parameters contain exactly 13 elements.")
  }

  map <- data.frame(State = c("BC", "AB", "SK", "MB", "ON", "QC", "YT", "NT",
                              "NU", "NB", "NL", "NS", "PE"),
                    X = c(0:5, 1:3, 6, 6, 7, 7),
                    Y = c(rep(2, 6), rep(3, 3), 1, 3, 0, 1),
                    stringsAsFactors = FALSE)

  # Make sure all provinces and territories are present
  if(!all(map$State %in% pt)){
    stop("It appears some province or territory names are repeated or missing.")
  }

  user_map <- data.frame(State = pt, scol = pt_colors, bcol =
                           border_colors, sncol = pt_name_colors,
                         stringsAsFactors = FALSE)

  map <- merge(map, user_map, by = "State")

  text_ <- NULL
  if(pt_names){
    text_ <- map$State
  }

  plotbox(map$X, map$Y, map$scol, map$bcol,
          n_xboxes = 8, n_yboxes = 4, text_ = text_,
          text_col = map$sncol, text_cex = pt_name_cex, text_font = font)

}
