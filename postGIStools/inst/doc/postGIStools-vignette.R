## ----echo=FALSE----------------------------------------------------------
# Don't run this vignette on CRAN
NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(purl = NOT_CRAN, eval = NOT_CRAN)
knitr::opts_chunk$set(warning=FALSE, message=FALSE, collapse=TRUE, comment = "#>")

## ----get_postgis_query---------------------------------------------------
library(RPostgreSQL)
library(postGIStools)

con <- dbConnect(PostgreSQL(), dbname = "d2u06to89nuqei", user = "mzcwtmyzmgalae",
                 host = "ec2-107-22-246-250.compute-1.amazonaws.com",
                 password = "UTv2BuwJUPuruhDqJthcngyyvO")

countries <- get_postgis_query(con, "SELECT name, iso2, capital, population,
                               translations, geom FROM country 
                               WHERE population > 1000000",
                               geom_name = "geom", hstore_name = "translations")

class(countries)

## ----data_str------------------------------------------------------------
str(countries@data[1:2,])

## ----hstore_select-------------------------------------------------------
head(countries$translations %->% "es")

## ----hstore_subset-------------------------------------------------------
countries$translations[5:7] %->% "fr"

## ----hstore_assign-------------------------------------------------------
countries$translations[2] %->% "nl" <- "Albanië"
countries$translations[3] %->% "fr" <- NULL
countries$translations[2:3]

## ----create_tmp_table, results = "hide"----------------------------------
dbSendQuery(con, paste("CREATE TEMP TABLE cty_tmp (name text,", 
                       "iso2 text PRIMARY KEY, capital text,",
                       "translations hstore, geom geometry)"))

## ----postgis_insert------------------------------------------------------
postgis_insert(con, countries[1:10,], "cty_tmp",
               write_cols = c("name", "iso2", "translations"),
               geom_name = "geom", hstore_name = "translations")

# Reimport to check
cty_tmp <- get_postgis_query(con, paste("SELECT name, iso2, capital,",
                                        "geom, translations FROM cty_tmp"),
                             geom_name = "geom", hstore_name = "translations")
head(cty_tmp@data)

## ----postgis_update------------------------------------------------------
postgis_update(con, countries[1:10,], "cty_tmp", id_cols = "iso2", 
               update_cols = "capital", geom_name = "geom", 
               hstore_name = "translations")

cty_tmp <- get_postgis_query(con, paste("SELECT name, iso2, capital,",
                                        "geom, translations FROM cty_tmp"),
                             geom_name = "geom", hstore_name = "translations")
head(cty_tmp@data)

## ----update_hstore-------------------------------------------------------
countries$translations[2] %->% "nl" <- NULL
countries$translations[3] %->% "fr" <- "Algérie"
 
postgis_update(con, countries[1:10,], "cty_tmp", id_cols = "iso2", 
               update_cols = "translations", geom_name = "geom", 
               hstore_name = "translations")

cty_tmp <- get_postgis_query(con, paste("SELECT name, iso2, capital,",
                                        "geom, translations FROM cty_tmp"),
                             geom_name = "geom", hstore_name = "translations")
cty_tmp@data[cty_tmp$iso2 %in% c("AL", "DZ"), ]

# Key deletion not reflected in database unless hstore_concat = FALSE
postgis_update(con, countries[1:10,], "cty_tmp", id_cols = "iso2", 
               update_cols = "translations", geom_name = "geom", 
               hstore_name = "translations", hstore_concat = FALSE)

cty_tmp <- get_postgis_query(con, paste("SELECT name, iso2, capital,",
                                        "geom, translations FROM cty_tmp"),
                             geom_name = "geom", hstore_name = "translations")
cty_tmp@data[cty_tmp$iso2 %in% c("AL", "DZ"), ]

