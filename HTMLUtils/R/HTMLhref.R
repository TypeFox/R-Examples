`HTMLhref` <- function#adds an href item to the current HTML page
###adds an href item to the current HTML page
(href, ##<< HTML reference/URL
 txt, ##<< text to display
 file = get(".HTML.file"),##<< file to write to
 append = TRUE##<< append to file (default TRUE)
 ) {
    cat(paste('<a href=\"', href, '\"> ', txt, ' </a>', sep = ""), append = append,  file = file)
}
