
# Test esearch() ---------------------------------------------------------

context("Testing 'esearch()'")

if (getOption('reutils.test.remote')) {
  ## retmode = 'xml'
  a <- esearch(term = "cancer", db = "pubmed", reldate = 60, datetype = "edat",
               retmax = 6, usehistory = TRUE)
  b <- esearch(term = "cancer", db = "pubmed", reldate = 60, datetype = "edat",
               retmax = 6, usehistory = FALSE)
  
  test_that("esearch() returns an 'esearch' object", {
    expect_is(a, "esearch")
    expect_is(b, "esearch")
  })
  
  test_that("'content()' returns a character vector or an XMLInternalDocument", {
    expect_that(content(a, "text"), is_a("character"))
    expect_that(content(b, "text"), is_a("character"))
    expect_that(content(a, "xml"), is_a("XMLInternalDocument"))
    expect_that(content(b, "xml"), is_a("XMLInternalDocument"))
    expect_that(content(a, "json"), throws_error("Cannot return data of retmode.+"))
    expect_that(content(b, "json"), throws_error("Cannot return data of retmode.+"))
    expect_that(content(a, 'parsed'), is_a("entrez_uid"))
    expect_that(content(b, 'parsed'), is_a("entrez_uid"))
  })
  
  test_that("Subsetting an 'esearch' returns an 'esearch' object", {
    expect_that(a[1:2], is_a("entrez_uid"))
    expect_that(b[1:2], is_a("entrez_uid"))
    expect_that(length(b[1:2]), equals(2))
  })
  
  test_that("'querykey', 'webenv', and 'database' return the appropriate results", {
    expect_equal(querykey(a), 1)
    expect_match(webenv(a), "NCID_+")
    expect_equal(database(a), "pubmed")
    
    expect_equal(querykey(b), NA_integer_)
    expect_equal(webenv(b), NA_character_)
    expect_equal(database(b), 'pubmed')
  })
  
  test_that("'rettype', and 'retmode' return the appropriate results", {
    expect_equal(rettype(a), "uilist")
    expect_match(retmode(a), "xml")
    
    expect_equal(rettype(b), "uilist")
    expect_match(retmode(b), "xml")
  })
  
  test_that("'uid' returns a character vector for esearch objdect", {
    expect_equal(uid(a), NA_character_)
    expect_is(uid(b), "character")
    expect_equal(length(uid(b)), 6)
  })
  
  ## retmode = 'json'
  a <- esearch(term = "cancer", db = "pubmed", reldate = 60, datetype = "edat",
               retmax = 6, usehistory = TRUE, retmode = 'json')
  b <- esearch(term = "cancer", db = "pubmed", reldate = 60, datetype = "edat",
               retmax = 6, usehistory = FALSE, retmode = 'json')
  
  test_that("'content()' returns a character vector or a json object", {
    expect_that(content(a, "text"), is_a("character"))
    expect_that(content(b, "text"), is_a("character"))
    expect_that(content(a, "xml"), throws_error("Cannot return data of retmode.+"))
    expect_that(content(b, "xml"), throws_error("Cannot return data of retmode.+"))
    expect_that(content(a, "json"), is_a("json"))
    expect_that(content(b, "json"), is_a("json"))
    expect_that(content(a, 'parsed'), is_a("entrez_uid"))
    expect_that(content(b, 'parsed'), is_a("entrez_uid"))
  })
  
  test_that("'retmode' returns the appropriate results", {
    expect_match(retmode(a), "json")
    expect_match(retmode(b), "json")
  })
  
  ## rettype = 'count'
  a <- esearch(term = "cancer", db = "pubmed", reldate = 60, datetype = "edat",
               rettype = "count", retmode = 'xml')
  b <- esearch(term = "cancer", db = "pubmed", reldate = 60, datetype = "edat",
               rettype = "count", retmode = 'json')
  
  test_that("'content()' returns a numeric vector", {
    expect_that(content(a, "parsed"), is_a("numeric"))
    expect_that(content(b, "parsed"), is_a("numeric"))
    expect_that(content(a, "xml"), is_a("XMLInternalDocument"))
    expect_that(content(b, "json"), is_a("json"))
  })
}




  

  


  

  

