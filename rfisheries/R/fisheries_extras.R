#' Plots data for an rfisheries result
#'
#'@importFrom ggplot2 ggplot theme_get theme_update ggtitle geom_line labs theme_set element_line element_blank aes
#' @param x A landings dataset belonging to either a species or a country.
#' @param linecolor Default line color is steelblue
#' @param linesize Default line size is 0.9
#' @param title Plot title. Title is generated based on species or country code. Specify one here only if you need something else.
#' @param ... additional arguments
#' @importFrom assertthat not_empty are_equal assert_that
#' @export
#' @examples \dontrun{
#' fish_plot(of_landings(country = 'CAN'))
#' fish_plot(of_landings(species = 'COD'))
#'}
fish_plot <- function(x, linecolor = "steelblue", linesize = 0.9, title = NULL, ...) {

assert_that(class(x) == "data.frame")
# This weird step is just to satisfy the notes in check()
year <- NA
catch <- NA
x[,3] <- toupper(x[, 3])
# Both datasets really should have 3 columns.
# Otherwise something is wrong

assert_that(not_empty(x))
assert_that(are_equal(ncol(x), 3))
# This needs to be a new assertion but I haven't figured out how to write new assertions for assert_that.
stopifnot(class(x) == "data.frame")

# Allows to check which type of landings data we're working with (country or species)
species_dataset <- c("year","catch", "species")
country_dataset <- c("year","catch", "country")

# Update ggplot theme for plots
old <- theme_get()


theme_update(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

if(identical(species_dataset, names(x))) {

if(is.null(title)) {
    english_name <- species_code_data[which(species_code_data$a3_code == unique(x$species)), ]$english_name
    title <- paste0("Landings for ", english_name, " (", unique(x$species), ")")
}
fish_plot <- ggplot(x, aes(year, catch)) +
geom_line(color = linecolor, size = linesize) +
labs(x = "Year", y = "Catch (in tonnes)") +
ggtitle(title)
}

# Make an identical plot but this time if a country is specified instead of species.

if(identical(country_dataset, names(x))) {
    if(is.null(title)) {
     country_name <- country_code_data[which(country_code_data$iso3c == unique(x$country)), ]$country
    title <- paste0("Landings for ", country_name, " (", unique(x$country), ")")
    }
fish_plot <-  ggplot(x, aes(year, catch)) +
geom_line(color = linecolor, size = linesize) +
labs(x = "Year", y = "Catch (in tonnes)") +
ggtitle(title)
}
# Reset theme back to old default.
theme_set(old)
# Return plot object
if(!identical(class(fish_plot), "function")) {
    return(fish_plot)
    }
}
