context("Datasets")

dataset <- readLines("dataset.json")

test_that("dataset names are correct", {
    expect_named(fromJSONstat(dataset, naming = "label"),
                 c("A dataset with array value", "A dataset with object value"))
    expect_named(fromJSONstat(dataset, naming = "id"), c("dataset", "dataset2"))
    d <- data.frame(V1 = "a", value = 1)
    expect_named(fromJSONstat(toJSONstat(d)), "dataset")
    expect_named(fromJSONstat(toJSONstat(list(d, d))), c("dataset", "dataset2"))
    expect_named(fromJSONstat(toJSONstat(list(a = d, d))), c("a", "dataset2"))
    expect_named(fromJSONstat(toJSONstat(list(d, b = d))), c("dataset", "b"))
    expect_named(fromJSONstat(toJSONstat(list(a = d, b = d))), c("a", "b"))
    expect_named(fromJSONstat(toJSONstat(list(a = d, a = d))),
                 c("a", "a (dataset2)"))
})

test_that("dataset names are correct for missing labels", {
    expect_named(fromJSONstat(dataset[-3], naming = "label"),
                 c("dataset", "A dataset with object value"))
})
