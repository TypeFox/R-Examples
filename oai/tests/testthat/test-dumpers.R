
# valid_dumper ------------------------------------------------------------

context("Dumper validity checker")

dummy_dumper <- function(res, args, as, x=1, ...) {
  c(list(res=res, args=args, as=as, x=x), list(...))
}


test_that("supplied all obligatory arguments works", {
  skip_on_cran()

  expect_true( valid_dumper( function(res, args, as) {}, NULL))
} )

test_that("missing obligatory argument throws error", {
  skip_on_cran()

  expect_error( valid_dumper(  function(res, args) {}, NULL ) )
  expect_error( valid_dumper(  function(res, as) {}, NULL ) )
  expect_error( valid_dumper(  function(args, as) {}, NULL ) )
} )

test_that("supplying additional arguments unaccepted by dumper", {
  skip_on_cran()

  expect_error( valid_dumper( function(res, args, as) {}, list(a=1) ) )
  expect_true( valid_dumper( function(res, args, as, ...) {}, list(a=1) ) )
} )

test_that("missing necessary extra argument", {
  skip_on_cran()

  expect_error( valid_dumper(
    function(res, args, as, x) {},
    NULL
  ) )
  expect_error( valid_dumper(
    function(res, args, as, x, ...) {},
    list(a=1)
  ) )
} )

test_that("user cannot provide res/args/as arguments", {
  skip_on_cran()

  expect_error( valid_dumper( function(res, args, as) {}, list(res=1) ) )
  expect_error( valid_dumper( function(res, args, as) {}, list(args=1) ) )
  expect_error( valid_dumper( function(res, args, as) {}, list(as=1) ) )
} )


test_that("user did not supply extra argument but dumper has a default", {
  skip_on_cran()

  expect_true( valid_dumper(dummy_dumper, NULL) )
} )




# dump_raw_to_txt ---------------------------------------------------------



context("Testing text file dumper")

test_that("list_identifiers saves raw XML to text files", {
  skip_on_cran()

  fnames <- list_identifiers(from = '2014-06-01T', until = '2014-06-01T', as="raw",
                             dumper=dump_raw_to_txt,
                             dumper_args=list(file_dir=tempdir()))
  expect_true(is.character(fnames))
  expect_true(all(sapply(fnames, file.exists)))

  xmls <- lapply(fnames, xml2::read_xml)
  expect_true( all(sapply(xmls, inherits, "xml_document")) )

  unlink(fnames)
} )


test_that("list_records saves raw XML to text files", {
  skip_on_cran()

  fnames <- list_records(from = '2014-06-01T', until = '2014-06-02T',
                         as = "raw",
                         dumper=dump_raw_to_txt,
                         dumper_args=list(file_dir=tempdir()))
  expect_true(is.character(fnames))
  expect_true(all(sapply(fnames, file.exists)))

  xmls <- lapply(fnames, xml2::read_xml)
  expect_true( all(sapply(xmls, inherits, "xml_document")) )

  unlink(fnames)
} )





# dump_to_rds -------------------------------------------------------------




context("Testing RDS file dumper")

test_that("list_identifiers saves raw XML to RDS files", {
  skip_on_cran()

  fnames <- list_identifiers(from="2014-06-01T", until="2014-06-02T", as="raw",
                             dumper=dump_to_rds,
                             dumper_args=list(file_dir=tempdir()))
  expect_true(is.character(fnames))
  expect_true(all(sapply(fnames, file.exists)))

  # read and parse
  xmls <- lapply(fnames, function(fn) xml2::read_xml(readRDS(fn)))
  expect_true( all(sapply(xmls, inherits, "xml_document")) )

  unlink(fnames)
} )


test_that("list_records saves raw XML to RDS files", {
  skip_on_cran()

  fnames <- list_records(from = '2014-06-01T', until = '2014-06-01T',
                         as = "raw",
                         dumper=dump_to_rds,
                         dumper_args=list(file_dir=tempdir()))
  expect_true(is.character(fnames))
  expect_true(all(sapply(fnames, file.exists)))

  # read and parse
  xmls <- lapply(fnames, function(fn) xml2::read_xml(readRDS(fn)))
  expect_true( all(sapply(xmls, inherits, "xml_document")) )

  unlink(fnames)
} )



# raw_to_db ---------------------------------------------------------------



context("Testing raw_to_db dumper with SQLite")


test_that("list_identifiers dumps results to SQLite", {
  skip_on_cran()

  con <- DBI::dbConnect(RSQLite::SQLite(), dbname=":memory:")
  dumprval <- list_identifiers(from = '2014-06-01T', until = '2014-06-01T', as="raw",
                             dumper=dump_raw_to_db,
                             dumper_args=list(dbcon=con, table_name="foo",
                                              field_name="bar") )
  expect_null(dumprval)
  expect_identical(DBI::dbListTables(con), "foo")
  expect_identical(DBI::dbListFields(con, "foo"), "bar")


  # parse
  xmldf <- DBI::dbReadTable(con, "foo")
  xmls <- lapply( xmldf$bar, xml2::read_xml )
  expect_true( all(sapply(xmls, inherits, "xml_document"))  )

  DBI::dbDisconnect(con)
} )


test_that("list_records dumps results to SQLite", {
  skip_on_cran()

  con <- DBI::dbConnect(RSQLite::SQLite(), dbname=":memory:")
  dumprval <- list_records(from = '2014-06-01T', until = '2014-06-01T',
                           as="raw",
                           dumper=dump_raw_to_db,
                           dumper_args=list(dbcon=con, table_name="foo",
                                            field_name="bar") )
  expect_null(dumprval)
  expect_identical(DBI::dbListTables(con), "foo")
  expect_identical(DBI::dbListFields(con, "foo"), "bar")

  # parse
  xmldf <- DBI::dbReadTable(con, "foo")
  xmls <- lapply( xmldf$bar, xml2::read_xml )
  expect_true( all(sapply(xmls, inherits, "xml_document"))  )

  DBI::dbDisconnect(con)
} )
