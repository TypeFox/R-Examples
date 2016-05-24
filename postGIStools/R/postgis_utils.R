### Utility functions used by package (not exported)

# Find SRID corresponding to proj4string p4s
#  first search spatial_ref_sys table, then try showEPSG. return 0 if not found
find_srid <- function(conn, p4s) {
    if (is.na(p4s)) return("0")
    srs <- tryCatch(RPostgreSQL::dbGetQuery(conn, "SELECT * FROM spatial_ref_sys"),
                    error = function(e) NULL)
    srid <- NULL
    if (!is.null(srs))
        srid <- srs$srid[srs$proj4text == p4s]
    if (length(srid) == 0) {
        srid <- tryCatch(rgdal::showEPSG(p4s), error = function(e) "0")
        if (grepl("OGRERR", srid)) srid <- "0"
    }
    srid
}

# Convert JSON string representation to hstore
json_to_hstore <- function(str) {
    # deal with NULL fields
    str <- stringr::str_replace_all(str, "\\:\\{\\}", "=>NULL")
    # replace punctuation
    stringr::str_replace_all(str, c("\\{" = "", "\\}" = "", "\\[" = "",
                                    "\\]" = "", "\\:" = "=>", "," = ", "))
}

# Make wrapper functions to quote identifiers and strings from given connection
make_id_quote <- function(conn) {
    function(s) DBI::dbQuoteIdentifier(conn, s)
}
make_str_quote <- function(conn) {
    function(s) DBI::dbQuoteString(conn, s)
}

# Test if argument is a single character string
test_single_str <- function(s) {
    if (!is.character(s) | length(s) != 1) {
        stop(paste(deparse(substitute(s)), "is not a single character string"))
    }
}
