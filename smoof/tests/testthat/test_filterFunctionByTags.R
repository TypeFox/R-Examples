context("function filtering works well")

test_that("filterFunctionByTags returns reasonable results", {
  expected_tags = c("multimodal", "scalable")
  # check if function works well
  fun.names = filterFunctionsByTags(expected_tags)
  expect_character(fun.names, any.missing = FALSE, all.missing = FALSE)

  # get tags of functions
  fun.generators = lapply(fun.names, getGeneratorByName)

  for (fun.generator in fun.generators) {
    tags = getTags(fun.generator)
    expect_true(isSubset(expected_tags, tags), info =
      sprintf("Generator '%s' filtered, but tags are missing: %s",
        attr(fun.generator, "name"), collapse(setdiff(tags, expected_tags)), sep = ","))
  }

  # now check independently
  fun.names2 = filterFunctionsByTags(expected_tags, or = TRUE)
  expect_subset(fun.names, fun.names2)
})
