#' Generate random fake full names.
#'
#' @param n number of observations.
#' @return A character vector of \code{n} fake randomly generated full names.
#' @examples
#' r_full_names(10)
#' @export
r_full_names <- function(n) {
  genders <- sample(c("male", "female"), size = n, replace = TRUE)
  genders[genders == "male"] <- male_names[sample(1:nrow(male_names), size = sum(genders == "male"), replace = TRUE), ]
  genders[genders == "female"] <- female_names[sample(1:nrow(female_names), size = sum(genders == "female"), replace = TRUE), ]
  return(paste(genders, surnames[sample(1:nrow(surnames), size = n, replace = TRUE), "surname"], sep = " "))
}
