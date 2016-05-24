#' Send SELECT query and parse geometry, hstore columns
#'
#' This function is an extension of \code{\link[DBI]{dbGetQuery}} that is useful
#' in cases where selected columns include a PostgreSQL hstore, which is parsed
#' as a list-column, and/or a PostGIS geometry, in which case the output is a
#' spatial data frame (from the \code{\link[sp]{sp}} package).
#'
#' Column names must be explicitly listed in the query (i.e. no \code{"SELECT *"})
#' and the geometry or hstore column should not be aliased (\code{AS}) or transformed.
#' The function issues a warning if a specified \code{geom_name} or
#' \code{hstore_name} does not appear in \code{statement}.
#'
#' Conversion to spatial data frame objects will fail if there are \code{NULL}
#' values in the geometry column, so these should be filtered out in the provided
#' query statement.
#'
#' @param conn A \code{\link[RPostgreSQL]{PostgreSQLConnection-class}} object,
#'   such as the output of \code{\link[DBI]{dbConnect}}.
#' @param statement Character string for a SQL SELECT query.
#' @param geom_name Name of the geometry column (\code{NA} if none).
#' @param hstore_name Name of the hstore column (\code{NA} if none).
#' @return Either a data frame (if \code{geom_name = NA}) or a
#'   Spatial[Points/Lines/Polygons]DataFrame containing the query result. If a
#'   hstore column is present, it appears as a list-column in the data frame,
#'   i.e. each cell is a named list of key-value pairs.
#'
#' @examples
#' \dontrun{
#' library(RPostgreSQL)
#' con <- dbConnect(PostgreSQL(), dbname = "my_db")
#'
#' # If geom column holds points, returns a SpatialPointsDataFrame
#' cities <- get_postgis_query(con, "SELECT name, geom, datalist FROM city",
#'                             geom_name = "geom", hstore_name = "datalist")
#'
#' # Get the populations (part of datalist hstore) as a vector
#' pop <- cities@data$datalist %->% "population"
#' }
#'
#' @references The code for importing geom fields is based on a blog post by
#'   Lee Hachadoorian: \href{http://www.r-bloggers.com/load-postgis-geometries-in-r-without-rgdal/}{Load PostGIS geometries in R without rgdal}.
#' @seealso The \code{\link{\%->\%}} operator for working with hstore columns;
#'   \code{\link{postgis_insert}} and \code{\link{postgis_update}} for writing
#'   to a PostgreSQL connection.
#' @export
get_postgis_query <- function(conn, statement, geom_name = NA_character_,
                              hstore_name = NA_character_) {
    # Check inputs
    if (!is(conn, "PostgreSQLConnection")) {
        stop("conn is not a valid PostgreSQL connection")
    }
    if (length(statement) != 1 | !grepl("^select", tolower(statement))) {
        stop("statement does not appear to be a SELECT query")
    }
    test_single_str(geom_name)
    test_single_str(hstore_name)

    new_query <- edit_select_query(conn, statement, geom_name, hstore_name)

    # Get query output (suppress warnings that json is unrecognized)
    res <- suppressWarnings(RPostgreSQL::dbGetQuery(conn, new_query))

    process_select_result(res, geom_name, hstore_name)
}


## Non-exported functions called by get_postgis_query

edit_select_query <- function(conn, statement, geom_name, hstore_name) {
    # Make shortcut functions for SQL quoting
    quote_id <- make_id_quote(conn)

    # Convert geom field to WKT text format
    #  and join with spatial_ref_sys for projection information
    new_query <- statement
    if (!is.na(geom_name)) {
        if (!grepl(geom_name, new_query, fixed = TRUE)) {
            warning(paste("geom. field", geom_name,
                          "absent from query statement."))
            geom_name <- NA
        } else {
            geom_sub <- paste0("ST_AsText(", quote_id(geom_name), ") AS ",
                               quote_id(paste0(geom_name, "_wkt")), ", ",
                               "ST_SRID(", quote_id(geom_name), ") AS ",
                               quote_id(paste0(geom_name, "_srid")))
            new_query <- gsub(geom_name, geom_sub, new_query, fixed = TRUE)
            new_query <- paste0("SELECT * FROM (", new_query, ") AS base_qry ",
                                "LEFT JOIN spatial_ref_sys ON ",
                                quote_id(paste0(geom_name, "_srid")), " = srid")
        }
    }

    # Add conversion of hstore field to JSON
    if (!is.na(hstore_name)) {
        if (!grepl(hstore_name, new_query, fixed = TRUE)) {
            warning(paste("hstore field", hstore_name,
                          "absent from query statement."))
            hstore_name <- NA
        } else {
            hstore_sub <- paste0("hstore_to_json(", quote_id(hstore_name),
                                 ") AS ", quote_id(paste0(hstore_name, "_json")))
            new_query <- gsub(hstore_name, hstore_sub, new_query, fixed = TRUE)
        }
    }

    new_query
}


process_select_result <- function(res, geom_name, hstore_name) {
    # Convert hstore (now JSON) column into list of lists
    if (!is.na(hstore_name)) {
        hs <- paste0(hstore_name, "_json")
        res[is.na(res[[hs]]), hs] <- "{}" # Replace SQL NULLs with empty lists
        res[[hstore_name]] <- lapply(as.character(res[[hs]]), jsonlite::fromJSON)
        res[[hs]] <- NULL
    }

    # If it has geom, convert into Spatial*DataFrame
    if (!is.na(geom_name)) {
        if (length(unique(res$proj4text)) > 1)
            stop("Returned geometries do not share the same projection.")
        proj <- unique(res$proj4text)
        geom_wkt <- paste0(geom_name, "_wkt")
        sp_obj <- do.call(rbind, Map(rgeos::readWKT, text = res[[geom_wkt]],
                                     id = 1:nrow(res), USE.NAMES = FALSE))
        proj4string(sp_obj) <- CRS(proj)
        # Spatial columns to be discarded
        sp_cols <- c(geom_wkt, paste0(geom_name, "_srid"), "srid", "auth_name",
                     "auth_srid", "srtext", "proj4text")
        dat <- res[!(names(res) %in% sp_cols)]
        if (is(sp_obj, "SpatialPoints")) {
            res <- SpatialPointsDataFrame(sp_obj, dat)
        } else if (is(sp_obj, "SpatialLines")) {
            res <- SpatialLinesDataFrame(sp_obj, dat)
        } else if (is(sp_obj, "SpatialPolygons")) {
            res <- SpatialPolygonsDataFrame(sp_obj, dat)
        } else {
            stop("geom. field cannot be mapped to Point, Line or Polygon type.")
        }
    }
    res
}
