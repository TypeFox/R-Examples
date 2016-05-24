#' Insert or update records in a PostgreSQL table from a R data frame.
#'
#' These functions produce INSERT (\code{postgis_insert}) or UPDATE (\code{postgis_update})
#' queries to write data from a R data frame to a PostgreSQL table, with options
#' to include a geometry layer and a list-column of key-value pairs (as a PostgreSQL hstore).
#' The queries are passed to the database with \code{\link[DBI]{dbSendQuery}}.
#'
#' All column names used in the query must match between the input data frame and
#' the target database table (except for \code{geom_name} which only applies to
#' the table).
#'
#' \code{postgis_update} creates an \emph{UPDATE ... SET ... FROM ...} query,
#' which effectively joins the original table and input data frame based on matching
#' values in \code{id_cols}, then updates the values in \code{update_cols}.
#' The combination of \code{id_cols} must be unique in \code{df}, but they can
#' be duplicated in the database table, in which case multiple rows are updated
#' from a single row in \code{df}. Neither the geometry nor the hstore column can be
#' used in \code{id_cols}.
#'
#' Note that if \code{hstore_concat = TRUE} (the default), hstore columns are updated
#' by \emph{concatenation}, i.e. new keys are added, values associated with
#' existing keys are updated, no keys are deleted. To overwrite whole hstore
#' "cells", potentially deleting keys absent in \code{df}, set \code{hstore_concat = FALSE}.
#'
#' @param conn A \code{\link[RPostgreSQL]{PostgreSQLConnection-class}} object,
#'   such as the output of \code{\link[DBI]{dbConnect}}.
#' @param df A data frame (if \code{geom_name = NA}) or
#'   Spatial[Points/Lines/Polygons]DataFrame.
#' @param tbl Name of the PostgreSQL table to write to.
#' @param geom_name Name of the geometry column in the database table
#'   (\code{NA} if none).
#' @param hstore_name Name of the hstore column in both \code{df} and the
#'   database table (\code{NA} if none).
#' @param write_cols A character vector, corresponding to the columns in
#'   \code{df} to insert in the database table. If \code{NA}, inserts all columns.
#' @return The result of \code{\link[DBI]{dbSendQuery}}.
#'
#' @examples
#' \dontrun{
#' library(RPostgreSQL)
#' con <- dbConnect(PostgreSQL(), dbname = "my_db")
#'
#' # Returns a SpatialPointsDataFrame
#' cities <- get_postgis_query(con, "SELECT name, geom, datalist FROM city",
#'                             geom_name = "geom", hstore_name = "datalist")
#'
#' # Create a new field in hstore and update DB
#' cities@data$datalist %->% "pop_density" <-
#'    cities@data$datalist %->% "population" / cities@data$datalist %->% "area"
#' postgis_update(con, cities, "city",
#'                id_cols = "name", update_cols = "datalist",
#'                geom_name = "geom", hstore_name = "datalist")
#'
#' # Add rows to DB with postgis_insert
#' # (new_cities is a SpatialPointsDataFrame with same columns as cities)
#' postgis_insert(con, new_cities, "city",
#'                geom_name = "geom", hstore_name = "datalist")
#' }
#'
#' @seealso \code{\link{get_postgis_query}} for the inverse operation
#'   (read from database to R).
#' @rdname postgis_insert_update
#' @export
postgis_insert <- function(conn, df, tbl, write_cols = NA,
                           geom_name = NA_character_,
                           hstore_name = NA_character_) {
    query_text <- prep_write_query(conn, df, tbl, mode = "insert", write_cols,
                                   NA, NA, geom_name, hstore_name, NA)
    RPostgreSQL::dbSendQuery(conn, query_text)
}


#' @param id_cols A character vector, corresponding to the columns in \code{df}
#'   used to match records between \code{df} and the database table.
#' @param update_cols A character vector, corresponding to the columns that
#'   must be updated in the database table based on values in \code{df}.
#' @param hstore_concat If TRUE, hstore columns are updated by concatenation.
#'
#' @rdname postgis_insert_update
#' @export
postgis_update <- function(conn, df, tbl, id_cols, update_cols,
                           geom_name = NA_character_,
                           hstore_name = NA_character_, hstore_concat = TRUE) {
    query_text <- prep_write_query(conn, df, tbl, mode = "update", NA, id_cols,
                             update_cols, geom_name, hstore_name, hstore_concat)
    RPostgreSQL::dbSendQuery(conn, query_text)
}


# Function to build query string for INSERT or UPDATE query (based on "mode")
#  (not exported)
prep_write_query <- function(conn, df, tbl, mode, write_cols, id_cols,
                           update_cols, geom_name, hstore_name, hstore_concat) {
    # Check inputs
    if (!is(conn, "PostgreSQLConnection")) {
        stop("conn is not a valid PostgreSQL connection")
    }
    if (!is(df, "data.frame") & !is(df, "Spatial")) {
        stop("df must be a data.frame or Spatial*DataFrame")
    }
    test_single_str(tbl)
    test_single_str(geom_name)
    test_single_str(hstore_name)
    if (mode == "update") {
        if (length(intersect(id_cols, update_cols)) > 0) {
            stop("the same column cannot appear in id_cols and update_cols")
        }
        if (geom_name %in% id_cols || hstore_name %in% id_cols) {
            stop("geometry and hstore columns cannot be used in id_cols")
        }
    }

    # Make shortcut functions for quoting
    quote_id <- make_id_quote(conn)
    quote_str <- make_str_quote(conn)

    # Subset columns of df based on write_cols (if not NA);
    #  if spatial data frame, convert spatial objects into WKT
    if (!is.na(geom_name)) {
        if (!is(df, "Spatial")) {
            stop("geom_name specified but df is not a spatial object")
        }
        if (all(is.na(write_cols))) {
            write_cols <- colnames(df@data)
        } else if (!all(write_cols %in% colnames(df@data))) {
            stop(paste("columns not found in df:",
                       paste(setdiff(write_cols, colnames(df@data)),
                             collapse = ", ")))
        }
        srid <- find_srid(conn, proj4string(df))
        # Note that writeWKT(..., byid = TRUE) behaves differently if only 1 row
        geom_wkt <- rgeos::writeWKT(df, byid = nrow(df@data) > 1)
        df <- cbind(df@data[, write_cols, drop = FALSE], geom_wkt,
                    stringsAsFactors = FALSE)
        igeom <- ncol(df)
        colnames(df)[igeom] <- geom_name
    } else {
        if (is(df, "Spatial")) {
            stop("geom_name must be specified since df is as spatial object")
        }
        if (all(is.na(write_cols))) {
            write_cols <- colnames(df)
        } else if (!all(write_cols %in% colnames(df))) {
            stop(paste("columns not found in df:",
                       paste(setdiff(write_cols, colnames(df)),
                             collapse = ", ")))
        }
        df <- df[, write_cols, drop = FALSE]
    }

    # Convert any factor columns to characters
    fact_cols <- vapply(df, is.factor, FALSE)
    df[, fact_cols] <- lapply(df[, fact_cols], as.character)

    # If mode == update, check that no id_cols or update_cols are missing
    #   and that id_cols uniquely identify records
    if (mode == "update") {
        if (!all(c(id_cols, update_cols) %in% colnames(df)))
            stop(paste("columns not found in df:",
                       paste(setdiff(c(id_cols, update_cols), colnames(df)),
                             collapse = ", ")))
        if (any(duplicated(df[, id_cols])))
            stop("id_cols do not uniquely identify rows in df")
    }

    # Convert hstore column (if present) to JSON text, then hstore text
    if (!is.na(hstore_name)) {
        if (!all(vapply(df[[hstore_name]], is.list, FALSE))) {
            stop(paste("column", deparse(substitute(hstore_name)),
                       "is not a valid hstore"))
        }
        df[[hstore_name]] <- vapply(df[[hstore_name]], jsonlite::toJSON, "")
        df[[hstore_name]] <- json_to_hstore(df[[hstore_name]])
    }

    # Build multirow INSERT or UPDATE query
    if (mode == "update") {
        tbl_q <- quote_id(tbl)
        tbl_tmp <- quote_id(paste0(tbl, "_tmp"))
        id_q <- quote_id(id_cols)
        update_q <- quote_id(update_cols)
        update_tmp <- paste0(tbl_tmp, ".", update_q)
        if (!is.na(hstore_name)) {
            i_hs <- which(update_cols == hstore_name)
            update_tmp[i_hs] <- paste0(update_tmp[i_hs], "::hstore")
            if (hstore_concat) {
                # hstore updated by concatenation (||), with COALESCE in case of NULLs
                update_tmp[i_hs] <- paste0("COALESCE(", tbl_q, ".", update_q[i_hs],
                                           ", hstore('')) || ", update_tmp[i_hs])
            }
        }
        query_text <- paste("UPDATE", tbl_q, "SET",
                            paste(update_q, "=", update_tmp, collapse = ", "),
                            "FROM (VALUES")
    } else {  # mode == "insert"
        col_ids <- quote_id(colnames(df))
        query_text <- paste("INSERT INTO", quote_id(tbl), "(",
                            paste(col_ids, collapse = ", "), ") VALUES")
    }

    # Convert data frame to character array
    #  quoting character columns and replacing NAs with NULLs
    df_q <- do.call(cbind, lapply(df, function(var) {
        if(is.character(var)) quote_str(var)
        else var
    }))
    na_inds <- do.call(cbind, lapply(df, is.na))
    df_q[na_inds] <- "NULL"

    if (!is.na(geom_name)) {
        df_q[, igeom] <- paste0("ST_GeomFromText(", df_q[, igeom], ", ", srid, ")")
    }

    query_values <- apply(df_q, 1, function(row) {
        paste("(", paste(row, collapse = ", "), ")")
    })
    query_text <- paste(query_text, paste(query_values, collapse = ", "))

    if (mode == "update") {
        query_text <- paste(query_text, ") AS", tbl_tmp, "(",
                            paste(quote_id(colnames(df)), collapse = ", "), ")",
                            "WHERE", paste(paste0(tbl_q, ".", id_q), "=",
                                           paste0(tbl_tmp, ".", id_q),
                                           collapse = " AND "))
    }

    query_text
}
