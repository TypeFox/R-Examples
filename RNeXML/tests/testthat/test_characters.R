context("character matrices")


## All tests will use this data file
f <- system.file("examples", "comp_analysis.xml", package="RNeXML")

test_that("we can parse XML to S4 and serialize S4 to XML for the basic character classes", {

  doc <- xmlParse(f)
  root <- xmlRoot(doc)

  char <- as(root[["characters"]][["format"]][["char"]], "char")
  out <- as(char, "XMLInternalElementNode")
  expect_is(char, "char") # not as dumb as it looks, at least we're checking our own method here 
  expect_is(out, "XMLInternalElementNode") # dumb check, but provides a dot to show the code above executed successfully 

  format <- as(root[["characters"]][["format"]], "format")
  out <- as(format, "XMLInternalElementNode")
  expect_is(format, "format") 
  expect_is(out, "XMLInternalElementNode") 

  matrix <- as(root[["characters"]][["matrix"]], "obsmatrix")
  out <- as(matrix, "XMLInternalElementNode")
  expect_is(matrix, "obsmatrix") 
  expect_is(out, "XMLInternalElementNode") 

  characters <- as(root[["characters"]], "characters")
  out <- as(characters, "XMLInternalElementNode")
  expect_is(characters, "characters") 
  expect_is(out, "XMLInternalElementNode") 

})

test_that("we can actually parse NeXML files containing character data", {
  nex <- read.nexml(f)
  expect_is(nex, "nexml")
})


## Now that we tested this, store the result so we can use it in later tests 
nex <- read.nexml(f)

test_that("we can extract character matrix with get_characters", {
  x <- get_characters(nex)
  expect_is(x, "data.frame")
## FIXME add additional and more precise expect_ checks 
})




test_that("we can extract a list of character matrices with get_characters_list", {
  x <- get_characters_list(nex)
  expect_is(x, "list")
  expect_is(x[[1]], "data.frame")
## FIXME add additional and more precise expect_ checks 

})


## 
test_that("add_otu can append only unmatched taxa to an existing otus block", {
  orig <- get_taxa(nex)
  x <- get_characters_list(nex)
  nex@otus[[1]]@otu <- new("ListOfotu", nex@otus[[1]]@otu[1:5]) # chop off some of the otu values

  new_taxa <- rownames(x[[1]]) 
  nex2 <- RNeXML:::add_otu(nex, new_taxa, append=TRUE) # add them back 
## should have same contents as orig... 
  get_taxa(nex2)
  expect_identical(sort(orig$label), sort(get_taxa(nex2)$label))

## Note that otu ids are not unique when we chop them off ...
})





## FIXME add_characters needs a method to add character names of states
## and then we need a test for that method


test_that("we can add characters to a nexml file using a data.frame", {
  f <- system.file("examples", "comp_analysis.xml", package="RNeXML")
  
  nex <- read.nexml(f)
  x <- get_characters(nex)
  nexml <- add_characters(x)

  ##  Can we write it out and read it back? 
  nexml_write(nexml, file = "chartest.xml")
  tmp <- nexml_read("chartest.xml")
  tmp_x <- get_characters(tmp)

  expect_is(tmp_x, "data.frame")
  expect_is(tmp, "nexml")
  unlink("chartest.xml")
})




## based on bug on 2014-03-12 65ae459523c529452adb699c3d5d118c0a207402 
test_that("we can add multiple character matrices to a nexml file", {
          library("geiger")
          data(geospiza)
          data(primates)
          nex <- add_characters(geospiza$dat)
          nex <- add_characters(primates$dat, nex)
          expect_is(nex, "nexml")
})

