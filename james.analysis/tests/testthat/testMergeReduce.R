context("Merge and reduce")

# load data
data(james)

# equality test for objects of class james
# (based on comparison of summary output)
expect_equal.james <- function(data1, data2, ...){
  data1.summary <- capture.output(summary(data1))
  data2.summary <- capture.output(summary(data2))
  return(all.equal(data1.summary, data2.summary))
}

test_that("merge and reduce retain class", {
  # duplicate data
  james.dup <- mergeJAMES(james, james)
  # extract coconut problem (from original data)
  james.coco <- reduceJAMES(james, problems = "coco")
  # check classes
  expect_is(james.dup, "james")
  expect_is(james.coco, "james")
})

test_that("merge after disjunct reduces yields original data", {
  # split PT and RD search results
  rd.data <- reduceJAMES(james, searches = "Descent")
  pt.data <- reduceJAMES(james, searches = "Tempering")
  # merge and compare with original data
  merged <- mergeJAMES(rd.data, pt.data)
  expect_equal.james(merged, james)
})

test_that("merge of data with itself doubles search runs", {
  # merge data with itself
  james.dup <- mergeJAMES(james, james)
  # verify
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      orig.num.runs <- getNumSearchRuns(james, p, s)
      dup.num.runs <- getNumSearchRuns(james.dup, p, s)
      expect_equal(dup.num.runs, 2*orig.num.runs)
    }
  }
})

test_that("reduce works with both regexes and lists of strings", {
  # extract coconut data using regex
  james.coco.1 <- reduceJAMES(james, problems = ".*coco.*")
  # extract coconut data using list
  james.coco.2 <- reduceJAMES(james, problems = list("coconut"))
  # compare
  expect_equal.james(james.coco.2, james.coco.1)
  expect_equal(getProblems(james.coco.1), c("coconut"))
  # same for maize sets
  james.maize.1 <- reduceJAMES(james, problems = "maize")
  james.maize.2 <- reduceJAMES(james, problems = list("maize-bulk", "maize-accession"))
  # verify
  expect_equal.james(james.maize.2, james.maize.1)
  expect_equal(getProblems(james.maize.1), c("maize-accession", "maize-bulk"))
  # extract random descent results
  james.rd.1 <- reduceJAMES(james, searches = "descent", ignore.case = TRUE)
  james.rd.2 <- reduceJAMES(james, searches = list("Random Descent"))
  # compare
  expect_equal.james(james.rd.1, james.rd.2)
  expect_equal(getProblems(james.rd.1), c("coconut", "maize-accession", "maize-bulk", "pea-small"))
  for(p in getProblems(james.rd.1)){
    expect_equal(getSearches(james.rd.1, p), "Random Descent")
  }
})

test_that("merge throws an error if second argument is not of the same class", {
  expect_error(mergeJAMES(james, "I am not of class JAMES"), "should also be of class \"james\"")
})

test_that("reduce complains when no data is retained", {
  expect_error(reduceJAMES(james, problems = "foo"), "no data is retained")
  expect_error(reduceJAMES(james, searches = "bar"), "no data is retained")
  # case sensitiveness
  expect_error(reduceJAMES(james, searches = "descent"), "no data is retained")
  reduceJAMES(james, searches = "descent", ignore.case = TRUE)
})








