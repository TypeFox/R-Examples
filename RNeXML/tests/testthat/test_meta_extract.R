context("extract_metadata")

nex <- add_basic_meta(
            title = "My test title",
            description = "A description of my test",
            creator = "Carl Boettiger <cboettig@gmail.com>",
            publisher = "unpublished data",
            pubdate = "2012-04-01",
            citation = citation("ape"))

test_that("we can extract metadata using the dedicated functions", {

  get_citation(nex)
  get_license(nex)
  get_metadata(nex)
  summary(nex)

  unlink("example.xml")
})



test_that("we can extract all available metadata at a specified level of the DOM", {
 get_metadata(nex) 
 get_metadata(nex, "trees") 
})


