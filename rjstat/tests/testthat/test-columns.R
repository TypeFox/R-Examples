context("Columns")

dataset <- readLines("dataset.json")

test_that("column names are correct", {
    expect_named(fromJSONstat(dataset, naming = "label")[[1]],
                 c("A dimension with array index",
                   "A dimension with object index",
                   "A dimension without index", "value"))
    expect_named(fromJSONstat(dataset, naming = "id")[[1]],
                 c("testdimension1", "testdimension2",
                   "testdimension3", "value"))
})

test_that("columns are correct", {
    expect_equal(fromJSONstat(dataset, naming = "label")[[1]][[1]],
                 c("Category 11", "Category 11", "Category 12", "Category 12"))
    expect_equal(fromJSONstat(dataset, naming = "id")[[1]][[1]],
                 c("testcategory11", "testcategory11", "testcategory12",
                   "testcategory12"))

    expect_equal(fromJSONstat(dataset, naming = "label")[[1]][[2]],
                 c("Category 21", "Category 22", "Category 21", "Category 22"))
    expect_equal(fromJSONstat(dataset, naming = "id")[[1]][[2]],
                 c("testcategory21", "testcategory22", "testcategory21",
                   "testcategory22"))

    expect_equal(fromJSONstat(dataset, naming = "label")[[1]][[3]],
                 c("Category 3", "Category 3", "Category 3", "Category 3"))
    expect_equal(fromJSONstat(dataset, naming = "id")[[1]][[3]],
                 c("testcategory3", "testcategory3", "testcategory3",
                   "testcategory3"))
})

test_that("column names are correct for missing labels", {
    expect_named(fromJSONstat(dataset[-24], naming = "label")[[1]],
                 c("testdimension1", "A dimension with object index",
                   "A dimension without index", "value"))
})

test_that("columns are correct for missing labels", {
    expect_equal(fromJSONstat(dataset[-31], naming = "label")[[1]][[1]],
                 c("testcategory11", "testcategory11", "testcategory12",
                   "testcategory12"))
    expect_equal(fromJSONstat(dataset[-44], naming = "label")[[1]][[2]],
                 c("testcategory21", "testcategory22", "testcategory21",
                   "testcategory22"))
})
