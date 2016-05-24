library(postGIStools)
context("postgis_insert and update")

# Connect to test database
con <- tryCatch(RPostgreSQL::dbConnect(RPostgreSQL::PostgreSQL(),
                             dbname = "d2u06to89nuqei", user = "mzcwtmyzmgalae",
                             host = "ec2-107-22-246-250.compute-1.amazonaws.com",
                             password = "UTv2BuwJUPuruhDqJthcngyyvO"),
                error = function(e) NULL)

if (!is.null(con)) {
    # Read in info for 10 countries (saved in file)
    country_sp <- readRDS("country_sp.rds")

    # Create new temporary table in DB
    RPostgreSQL::dbSendQuery(con, paste("CREATE TEMP TABLE cty_tmp (name text,",
                                        "iso2 text PRIMARY KEY, capital text,",
                                        "population integer, translations hstore,
                                        geom geometry)"))

    # Convenience function to re-import table from DB
    import_cty_tmp <- function() {
        qry <- get_postgis_query(con, paste("SELECT name, iso2, capital, population,",
                                            "translations, geom FROM cty_tmp"),
                                 geom_name = "geom", hstore_name = "translations")
        qry[order(match(qry$iso2, country_sp$iso2)), ]
    }


    # Insert first two rows and re-import
    postgis_insert(con, country_sp[1:2, ], "cty_tmp",
                   geom_name = "geom", hstore_name = "translations")
    qry <- import_cty_tmp()
}

test_that("postgis_insert correctly inserts full rows", {
    if (is.null(con)) skip("PostgreSQL connection unavailable")
    expect_equal(country_sp@polygons[1:2], qry@polygons)
    expect_equal(country_sp@data[1:2, -5], qry@data[, -5])
    expect_equal(country_sp$translations[1:2] %->% "it",
                 qry$translations %->% "it")
})


if (!is.null(con)) {
    # Insert partial information for a few rows
    postgis_insert(con, country_sp[3, ], "cty_tmp",
                   write_cols = c("name", "iso2"), geom_name = "geom")
    postgis_insert(con, country_sp[4, ], "cty_tmp",
                   write_cols = c("name", "iso2", "translations"),
                   geom_name = "geom", hstore_name = "translations")
    qry <- import_cty_tmp()
}

test_that("postgis_insert correctly inserts partial rows", {
    if (is.null(con)) skip("PostgreSQL connection unavailable")
    expect_equal(country_sp@polygons[1:4], qry@polygons)
    expect_equal(country_sp$translations[c(1, 2, 4)] %->% "es",
                 qry$translations[c(1, 2, 4)] %->% "es")
})


if (!is.null(con)) {
    # Use postgis_update to fill in missing fields in rows 3-4
    #  also add, delete and replace values in row 4 hstore
    country_sp$translations[4] %->% "en" <- "Algeria"
    country_sp$translations[4] %->% "fr" <- "AlgErie"
    country_sp$translations[4] %->% "it" <- NULL
    postgis_update(con, country_sp[3:4, ], "cty_tmp", id_cols = "iso2",
                   update_cols = c("capital", "population", "translations"),
                   geom_name = "geom", hstore_name = "translations")
    qry <- import_cty_tmp()

}

test_that("postgis_update works with basic data types", {
    if (is.null(con)) skip("PostgreSQL connection unavailable")
    expect_equal(country_sp@data[1:4, 1:4], qry@data[, 1:4])
})

test_that("postgis_update correctly inserts new hstore", {
    if (is.null(con)) skip("PostgreSQL connection unavailable")
    expect_equal(country_sp$translations[3] %->% "it",
                 qry$translations[3] %->% "it")
})

test_that("postgis_update correctly updates existing hstore", {
    if (is.null(con)) skip("PostgreSQL connection unavailable")
    expect_equal(length(qry$translations[[4]]), 6)
    expect_equal(country_sp$translations[4] %->% "en",
                 qry$translations[4] %->% "en")
    expect_equal(country_sp$translations[4] %->% "fr",
                 qry$translations[4] %->% "fr")
})


if (!is.null(con)) {
    postgis_update(con, country_sp[4, ], "cty_tmp", id_cols = "iso2",
                   update_cols = "translations", geom_name = "geom",
                   hstore_name = "translations", hstore_concat = FALSE)
    qry <- import_cty_tmp()
}

test_that("postgis_update deletes key if hstore_concat is FALSE", {
    if (is.null(con)) skip("PostgreSQL connection unavailable")
    expect_true(is.na(qry$translations[4] %->% "it"))
})


if (!is.null(con)) {
    # Import data with unknown projection
    RPostgreSQL::dbSendQuery(con, paste("CREATE TEMP TABLE cty_tmp2 (name text,",
                                        "iso2 text PRIMARY KEY, capital text,",
                                        "population integer, translations hstore,
                                        geom geometry)"))
    country_no_proj <- country_sp
    proj4string(country_no_proj) <- NA_character_

    postgis_insert(con, country_no_proj, "cty_tmp2", geom_name = "geom",
                   hstore_name = "translations")
    qry <- get_postgis_query(con, "SELECT iso2, geom FROM cty_tmp2",
                             geom_name = "geom")
    qry <- qry[order(match(qry$iso2, country_no_proj$iso2)), ]
}

test_that("postgis_insert works with unknown projection", {
    if (is.null(con)) skip("PostgreSQL connection unavailable")
    expect_equal(country_no_proj@polygons, qry@polygons)
})


test_that("postgis_insert and update fail on bad inputs", {
    if (is.null(con)) skip("PostgreSQL connection unavailable")
    expect_error(postgis_insert(con, country_sp@data[5:6, ], "cty_tmp",
                                write_cols = "currency"))
    expect_error(postgis_update(con, country_sp@data, "cty_tmp",
                                id_cols = "iso2", update_cols = "test"))
    # No geom_name for spatial and vice versa
    expect_error(postgis_insert(con, country_sp[7:8, ], "cty_tmp"))
    expect_error(postgis_insert(con, country_sp@data[9:10, ], "cty_tmp",
                                geom_name = "geom"))
    # Unallowed id columns
    expect_error(postgis_update(con, country_sp@data, "cty_tmp",
                                id_cols = c("iso2", "name"),
                                update_cols = c("iso2", "capital")))
    expect_error(postgis_update(con, country_sp, "cty_tmp", id_cols = "geom",
                                update_cols = "name", geom_name = "geom"))
})


# Disconnecting deletes TEMP tables
if (!is.null(con)) RPostgreSQL::dbDisconnect(con)
