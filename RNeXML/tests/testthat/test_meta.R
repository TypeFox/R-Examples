context("meta")


data(bird.orders)

test_that("We can add additional metadata", {
  ## The short version using an RNeXML API
  nex <- add_basic_meta(
              title = "My test title",
              description = "A description of my test",
              creator = "Carl Boettiger <cboettig@gmail.com>",
              publisher = "unpublished data",
              pubdate = "2012-04-01")
  write.nexml(nex, file = "meta_example.xml")

  expect_true_or_null(nexml_validate("meta_example.xml"))
  expect_is(nexml_read("meta_example.xml"), "nexml")

  unlink("meta_example.xml") # cleanup

})


test_that("We can add R bibentry type metadata", {
  ## The short version using an RNeXML API

  nex <- add_trees(bird.orders)
  nex <- add_basic_meta(nexml=nex, citation=citation("ape")) 
  write.nexml(nex, file = "meta_example.xml")

  expect_true_or_null(nexml_validate("meta_example.xml"))
  expect_is(nexml_read("meta_example.xml"), "nexml")
  unlink("meta_example.xml") # cleanup

})



test_that("We can add additional metadata", {
  ## The short version using an RNeXML API
  nex <- add_trees(bird.orders)
  nex <- add_basic_meta(nexml = nex, citation=citation("ape")) 

  history <- meta(property = "skos:historyNote",
                  content = "Mapped from the bird.orders data in the ape package using RNeXML",
                  id = "meta5144")
  modified <- meta(property = "prism:modificationDate",
                  content = "2013-10-04")
  website <- meta(href = "http://carlboettiger.info", 
                 rel = "foaf:homepage")

  nex <- add_meta(list(history, modified, website), 
                  nex, 
                  namespaces = c(skos = "http://www.w3.org/2004/02/skos/core#",
                                 prism = "http://prismstandard.org/namespaces/1.2/basic/", # check and remove duplicates
                                 foaf = "http://xmlns.com/foaf/0.1/"))

  nexml_write(nex, file = "meta_example.xml")  

  expect_true_or_null(nexml_validate("meta_example.xml"))
  expect_is(nexml_read("meta_example.xml"), "nexml")
  unlink("meta_example.xml") # cleanup

})


test_that("We can directly add additional metadata at arbitrary level", {
  nex <- add_trees(bird.orders)
  modified <- meta(property = "prism:modificationDate",
                  content = "2013-10-04")

  nex@trees[[1]]@meta <- new("ListOfmeta", list(modified))
  get_metadata(nex, "trees") %>%
    dplyr::filter(property == "prism:modificationDate") %>% 
    dplyr::select(content) -> 
    tmp
  expect_identical(tmp[[1]], modified@content)
})



test_that("We can directly add additional metadata using concatenation notation", {
  nex <- add_trees(bird.orders)
  modified <- meta(property = "prism:modificationDate",
                  content = "2013-10-04")
  website <- meta(href = "http://carlboettiger.info", 
                 rel = "foaf:homepage")

  nex@trees[[1]]@meta <- c(modified)         # we can add just one element 
  nex@trees[[1]]@meta <- c(modified,website) # or more than one element

  get_metadata(nex, "trees") %>%
    dplyr::filter(property == "prism:modificationDate") %>% 
    dplyr::select(content) -> 
    tmp
  
  expect_identical(tmp[[1]], modified@content)
})





test_that("We can add arbitrary metadata", {
  rdfa <- '<xhtml:div typeof="foaf:Person" about="http://carlboettiger.info#me">
             <a rel="foaf:account" href="https://twitter.com/cboettig">twitter</a> 
             <a rel="foaf:account" href="https://github.com/cboettig">github</a>
           </xhtml:div>'
  parsed <- xmlRoot(xmlParse(rdfa))
  arbitrary_rdfa <- meta(property="eml:additionalMetadata", content="additional metadata", children = parsed)

  nex <- add_meta(arbitrary_rdfa, 
                  namespaces = c(foaf = "http://xmlns.com/foaf/0.1/", 
                                 eml = "eml://ecoinformatics.org/eml-2.1.1", 
                                 xhtml = "http://www.w3.org/1999/xhtml"))

  nexml_write(nex, file = "example.xml")


  
  expect_is(nexml_read("example.xml", "nexml"), "nexml")
  unlink("example.xml") # cleanup

})


test_that("we can write numeric types of meta elements and get correct datatype", {
          m <- meta(property="numericTest", content = 3.141)
          expect_is(m@content, "character")
          expect_match(m@datatype, ".*:decimal")
})
