## ------------------------------------------------------------------------
con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

## ------------------------------------------------------------------------
table <- "mydata"
sql <- c(sprintf("CREATE TABLE IF NOT EXISTS %s", table),
         "(name string PRIMARY KEY,",
         "value blob)")
DBI::dbGetQuery(con, paste(sql, collapse=" "))

## ------------------------------------------------------------------------
value <- mtcars
name <- "mtcars"

dat <- data.frame(name=name,
                  value=I(list(serialize(value, NULL))),
                  stringsAsFactors=FALSE)
sql <- sprintf("INSERT into %s (name, value) values (:name, :value)", table)
RSQLite::dbGetPreparedQuery(con, sql, bind.data=dat)

## ---- eval=FALSE---------------------------------------------------------
#  DBI::dbGetQuery(con, sql, bind.data=dat)

## ------------------------------------------------------------------------
sql <- sprintf('SELECT value FROM %s WHERE name == "%s"', table, name)
x <- unserialize(DBI::dbGetQuery(con, sql)[[1]][[1]])
identical(x, value)

## ------------------------------------------------------------------------
driver_sqlite <- function(path, tbl_data="storr_data", tbl_keys="storr_keys") {
  .R6_driver_sqlite$new(path, tbl_data, tbl_keys)
}

## ------------------------------------------------------------------------
.R6_driver_sqlite <- R6::R6Class(
  "driver_sqlite",
  public=list(
    ## Public data members
    con=NULL,
    tbl_data=NULL,
    tbl_keys=NULL,

    ## On initialisation we'll create the two tables but only if they
    ## do not exist.  We can enforce the constraint that hash must be
    ## unique within tbl_data and key/`namespace pairs must be unique
    ## within tbl_keys.
    initialize=function(path, tbl_data, tbl_keys) {
      self$con <- DBI::dbConnect(RSQLite::SQLite(), path)
      self$tbl_data <- tbl_data
      self$tbl_keys <- tbl_keys

      sql <- c(sprintf("CREATE TABLE if NOT EXISTS %s", tbl_data),
               "(hash string PRIMARY KEY NOT NULL,",
               "value blob NOT NULL)")
      DBI::dbGetQuery(self$con, paste(sql, collapse=" "))

      sql <- c(sprintf("CREATE TABLE IF NOT EXISTS %s", tbl_keys),
               "(namespace string NOT NULL,",
               "key string NOT NULL,",
               "hash string NOT NULL,",
               "PRIMARY KEY (namespace, key))")
      DBI::dbGetQuery(self$con, paste(sql, collapse=" "))
    },

    ## This is purely for identification later.
    type=function() {
      "DBI/sqlite"
    },

    ## Total destruction of the driver; delete all data stored in both
    ## tables, then delete our database connection to render the
    ## driver useless.
    destroy=function() {
      DBI::dbRemoveTable(self$con, self$tbl_data)
      DBI::dbRemoveTable(self$con, self$tbl_keys)
      self$con <- NULL
    },

    ## Return the hash value given a key/namespace pair
    get_hash=function(key, namespace) {
      sql <- sprintf('SELECT hash FROM "%s" WHERE namespace="%s" AND key="%s"',
                     self$tbl_keys, namespace, key)
      DBI::dbGetQuery(self$con, sql)[[1]]
    },
    ## Set the key/namespace pair to a hash
    set_hash=function(key, namespace, hash) {
      sql <- c(sprintf("INSERT OR REPLACE INTO %s", self$tbl_keys),
               sprintf('(namespace, key, hash) VALUES ("%s", "%s", "%s")',
                       namespace, key, hash))
      DBI::dbGetQuery(self$con, paste(sql, collapse=" "))
    },

    ## Return a (deserialised) R object, given a hash
    get_object=function(hash) {
      sql <- c(sprintf("SELECT value FROM %s", self$tbl_data),
               sprintf('WHERE hash = "%s"', hash))
      value <- DBI::dbGetQuery(self$con, paste(sql, collapse=" "))[[1]]
      unserialize(value[[1]])
    },

    ## Set a (serialised) R object against a hash.  This would be
    ## considerably simpler (but probably slower and less accurate) if we
    ## serialised to string with:
    ##   rawToChar(serialize(value, NULL, TRUE))
    set_object=function(hash, value) {
      dat <- data.frame(hash=hash,
                        value=I(list(serialize(value, NULL))),
                        stringsAsFactors=FALSE)
      sql <- c(sprintf("INSERT OR REPLACE INTO %s", self$tbl_data),
               "(hash, value) VALUES (:hash, :value)")
      RSQLite::dbGetPreparedQuery(self$con, paste(sql, collapse=" "),
                                  bind.data=dat)
    },

    ## Check if a key/namespace pair exists.
    exists_hash=function(key, namespace) {
      sql <- sprintf('SELECT 1 FROM %s WHERE namespace = "%s" AND key = "%s"',
                     self$tbl_keys, namespace, key)
      nrow(DBI::dbGetQuery(self$con, sql)) > 0L
    },
    ## Check if a hash exists
    exists_object=function(hash) {
      sql <- sprintf('SELECT 1 FROM %s WHERE hash = "%s"',
                     self$tbl_data, hash)
      nrow(DBI::dbGetQuery(self$con, sql)) > 0L
    },

    ## Delete a key.  Because of the requirement to return TRUE/FALSE on
    ## successful/unsuccessful key deletion this includes an exists_hash()
    ## step first.
    del_hash=function(key, namespace) {
      if (self$exists_hash(key, namespace)) {
        sql <- sprintf('DELETE FROM %s WHERE namespace = "%s" AND key = "%s"',
                       self$tbl_keys, namespace, key)
        DBI::dbGetQuery(self$con, sql)
        TRUE
      } else {
        FALSE
      }
    },
    ## Delete a hash
    del_object=function(hash) {
      if (self$exists_object(hash)) {
        sql <- sprintf('DELETE FROM %s WHERE hash = "%s"', self$tbl_data, hash)
        DBI::dbGetQuery(self$con, sql)
        TRUE
      } else {
        FALSE
      }
    },

    ## List hashes, namespaces and keys.  Because the SQLite driver seems to
    ## return numeric(0) if the result set is empty, we need as.character here.
    list_hashes=function() {
      sql <- sprintf("SELECT hash FROM %s", self$tbl_data)
      as.character(DBI::dbGetQuery(self$con, sql)[[1]])
    },
    list_namespaces=function() {
      sql <- sprintf("SELECT DISTINCT namespace FROM %s", self$tbl_keys)
      as.character(DBI::dbGetQuery(self$con, sql)[[1]])
    },
    list_keys=function(namespace) {
      sql <- sprintf('SELECT key FROM %s WHERE namespace="%s"',
                     self$tbl_keys, namespace)
      as.character(DBI::dbGetQuery(self$con, sql)[[1]])
    }
  ))

## ------------------------------------------------------------------------
dr <- driver_sqlite(":memory:")

## ------------------------------------------------------------------------
dr$list_hashes()

## ------------------------------------------------------------------------
hash <- digest::digest(mtcars)
dr$exists_object(hash)

## ------------------------------------------------------------------------
dr$set_object(hash, mtcars)

## ------------------------------------------------------------------------
dr$exists_object(hash)

## ------------------------------------------------------------------------
head(dr$get_object(hash))

## ------------------------------------------------------------------------
dr$list_hashes()

## ------------------------------------------------------------------------
dr$del_object(hash)

## ------------------------------------------------------------------------
dr$list_hashes()
dr$exists_object(hash)

## ------------------------------------------------------------------------
key <- "aaa"
namespace <- "ns"
dr$set_hash(key, namespace, hash)

## ------------------------------------------------------------------------
dr$exists_hash(key, namespace)

## ------------------------------------------------------------------------
dr$list_keys(namespace)

## ------------------------------------------------------------------------
dr$get_hash(key, namespace)

dr$del_hash(key, namespace)
dr$exists_hash(key, namespace)
dr$list_keys(namespace)

## ------------------------------------------------------------------------
storr::test_driver(function() driver_sqlite(":memory:"))

## ------------------------------------------------------------------------
storr_sqlite <- function(path,
                         tbl_data="storr_data", tbl_keys="storr_keys",
                         default_namespace="objects") {
  storr::storr(driver_sqlite(path, tbl_data, tbl_keys),
               default_namespace)
}

## ------------------------------------------------------------------------
st_sql <- storr_sqlite(":memory:")

## ------------------------------------------------------------------------
st_sql$list()

## ------------------------------------------------------------------------
st_sql$set("foo", runif(10))

## ------------------------------------------------------------------------
st_sql$get("foo")

## ------------------------------------------------------------------------
st_sql$del("foo")
st_sql$list()

## ------------------------------------------------------------------------
st_sql$list_hashes()
st_sql$get_value(st_sql$list_hashes())

## ------------------------------------------------------------------------
st_sql$gc()
st_sql$list_hashes()

## ------------------------------------------------------------------------
st_rds <- storr::storr_rds(tempfile())
microbenchmark::microbenchmark(
  st_sql$set(key, runif(10), use_cache=FALSE),
  st_rds$set(key, runif(10), use_cache=FALSE))

microbenchmark::microbenchmark(
  st_sql$get(key, use_cache=FALSE),
  st_rds$get(key, use_cache=FALSE))

st_sql$destroy()
st_rds$destroy()

