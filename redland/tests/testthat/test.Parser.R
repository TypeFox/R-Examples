context("Parser tests")
test_that("redland library loads", {
  library(redland)
})
test_that("Parser constructor", {
  library(redland)
  world <- new("World")
  expect_false(is.null(world))

  # Test creating the Storage system
  storage <- new("Storage", world, "hashes", name="", options="hash-type='memory'")
  expect_false(is.null(storage))
  expect_that(class(storage@librdf_storage), matches("_p_librdf_storage_s"))

  # Test creating the Model
  model <- new("Model", world, storage, options="")
  expect_false(is.null(model))
  expect_that(class(model@librdf_model), matches("_p_librdf_model_s"))

  # Test that model creation fails if world is not provided or is null
  err <- try(model <- new("Model", world=NULL, storage, options=""), silent=TRUE)
  expect_that(class(err), matches("try-error"))

  expect_false(is.null(model))
  expect_that(class(model@librdf_model), matches("_p_librdf_model_s"))

  # Test parsing an RDF document into a Model
  parser <- new("Parser", world)
  expect_false(is.null(parser))
  expect_that(class(parser@librdf_parser), matches("_p_librdf_parser_s"))

  parseFileIntoModel(parser, world, system.file('extdata/example.rdf', package='redland'), model)

  # Test creating a Serializer and serializing the content just parsed into the model
  serializer <- new("Serializer", world)
  expect_false(is.null(serializer))
  expect_that(class(serializer@librdf_serializer), matches("_p_librdf_serializer_s"))

  # Test performing a serialization on an RDF model
  rdf <- serializeToCharacter(serializer, world, model)
  expect_that(rdf, matches("John Smith"))

  err <- try(freeModel(model), silent=TRUE)
  expect_false(class(err) ==  "try-error")

  err <- try(freeStorage(storage), silent=TRUE)
  expect_false(class(err) == "try-error")

  err <- try(freeParser(parser), silent=TRUE)
  expect_false(class(err) == "try-error")

  err <- try(freeSerializer(serializer), silent=TRUE)
  expect_false(class(err) == "try-error")

  err <- try(freeWorld(world), silent=TRUE)
  expect_false(class(err) == "try-error")
})
