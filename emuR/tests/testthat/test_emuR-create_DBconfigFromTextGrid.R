##' testthat tests for convert_TextGridCollection_to_emuDB
##'

context("testing create_DBconfigFromTextGrid function")

path2demoData = file.path(tempdir(), "emuR_demoData")
path2testData = file.path(tempdir(), "emuR_testthat")
path2tg = file.path(path2demoData, "TextGrid_collection/msajc003.TextGrid")

dbName = 'test12'

# tmp project base path
basePath=file.path(tempdir(), dbName)

##############################
test_that("test that correct values are set for msajc003", {
  conf = create_DBconfigFromTextGrid(path2tg, dbName,basePath)
  expect_equal(length(conf$linkDefinitions), 0)
  expect_equal(length(conf$ssffTrackDefinitions), 0)
  expect_equal(length(conf$levelDefinitions), 11)
  expect_equal(conf$mediafileExtension, 'wav')
  
  expect_equal(conf$levelDefinitions[[1]]$name, 'Utterance')
  expect_equal(conf$levelDefinitions[[1]]$type, 'SEGMENT')
  expect_equal(conf$levelDefinitions[[1]]$attributeDefinitions[[1]]$name, 'Utterance')
  expect_equal(conf$levelDefinitions[[1]]$attributeDefinitions[[1]]$type, 'STRING')
  
})

##############################
test_that("test only correct tiers are extracted if tierNames is set", {
  conf = create_DBconfigFromTextGrid(path2tg, dbName, basePath, c("Phonetic", "Tone"))

  expect_equal(conf$levelDefinitions[[1]]$name, 'Phonetic')
  expect_equal(conf$levelDefinitions[[1]]$type, 'SEGMENT')
  expect_equal(conf$levelDefinitions[[1]]$attributeDefinitions[[1]]$name, 'Phonetic')
  expect_equal(conf$levelDefinitions[[1]]$attributeDefinitions[[1]]$type, 'STRING')

  expect_equal(conf$levelDefinitions[[2]]$name, 'Tone')
  expect_equal(conf$levelDefinitions[[2]]$type, 'EVENT')
  expect_equal(conf$levelDefinitions[[2]]$attributeDefinitions[[1]]$name, 'Tone')
  expect_equal(conf$levelDefinitions[[2]]$attributeDefinitions[[1]]$type, 'STRING')
  
})
