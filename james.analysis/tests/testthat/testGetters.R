context("Extract data")

# load data
data(james)

# test getProblems

test_that("getProblems returns correct ordered problem names", {
  expect_equal(getProblems(james), c("coconut", "maize-accession", "maize-bulk", "pea-small"))
})

test_that("getProblems respects applied filter", {
  expect_equal(getProblems(james, filter = "maize"), c("maize-accession", "maize-bulk"))
})

test_that("getProblems passes additional arguments to grep", {
  expect_equal(getProblems(james, filter = "Coco"), character(0))
  expect_equal(getProblems(james, filter = "Coco", ignore.case = TRUE), "coconut")
})

# test getSearches

test_that("getSearches returns correct ordered search names", {
  for(p in getProblems(james)){
    expect_equal(getSearches(james, p), c("Parallel Tempering", "Random Descent"))
  }
})

test_that("getSearches respects applied filter", {
  for(p in getProblems(james)){
    expect_equal(getSearches(james, p, filter = "Descent"), "Random Descent")
  }
})

test_that("getSearches passes additional arguments to grep", {
  for(p in getProblems(james)){
    expect_equal(getSearches(james, p, filter = "descent"), character(0))
    expect_equal(getSearches(james, p, filter = "descent", ignore.case = TRUE), "Random Descent")
  }
})

test_that("getSearches allows to omit problem name in case of a single problem", {
  expect_error(getSearches(james), "more than one problem")
  # now for coconut problem only
  james.coco <- reduceJAMES(james, problems = "coco")
  expect_equal(getSearches(james.coco), c("Parallel Tempering", "Random Descent"))
})

test_that("getSearches complains if the given problem is not found", {
  expect_error(getSearches(james, "foo"), "does not contain results for problem")
})

# test getSearchRuns

test_that("getSearchRuns allows to omit problem name in case of a single problem", {
  expect_error(getSearchRuns(james, search = "Random Descent"), "more than one problem")
  # now for coconut problem only
  james.coco <- reduceJAMES(james, problems = "coco")
  getSearchRuns(james.coco, search = "Random Descent")
})

test_that("getSearchRuns allows to omit search name in case of a single search for the given problem", {
  expect_error(getSearchRuns(james, problem = "coconut"), "more than one search for problem")
  # now for random descent only
  james.rd <- reduceJAMES(james, searches = "Descent")
  getSearchRuns(james.rd, problem = "coconut")
  # omit both problem and search name
  james.rd.coco <- reduceJAMES(james.rd, problems = "coco")
  getSearchRuns(james.rd.coco)
})

test_that("getSearchRuns complain when given an unknown problem", {
  expect_error(getSearchRuns(james, problem = "foo"), "does not contain results for problem")
})

test_that("getSearchRuns complain when given an unknown search", {
  expect_error(getSearchRuns(james, problem = "coconut", search = "bar"), "does not contain results for search")
})

test_that("getSearchRuns returns a list of the expected format", {
  runs <- getSearchRuns(james, problem = "coconut", search = "Random Descent")
  # check size
  expect_equal(length(runs), 10)
  # check each run
  for(run in runs){
    # check times
    expect_false(is.null(run$times)) # present
    expect_true(is.vector(run$times)) # vector
    expect_true(is.numeric(run$times)) # numeric
    expect_true(all(run$times[run$times != -1] > 0)) # all positive or -1
    expect_true(is.sorted(run$times)) # increasing
    # check values
    expect_false(is.null(run$values)) # present
    expect_true(is.vector(run$values)) # vector
    expect_true(is.numeric(run$values)) # numeric
    expect_true(is.sorted(run$values) || is.sorted(rev(run$values))) # increasing or decreasing
    # compare number of times and values, should match
    expect_equal(length(run$times), length(run$values))
    # check presence of best solution
    expect_false(is.null(run$best.solution))
  }
})

test_that("getSearchRuns returns the correct values", {
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      file <- sprintf("files/search-runs-%s-%s.rds", gsub(" ", "-", p), gsub(" ", "-", s))
      expect_equal_to_reference(getSearchRuns(james, p, s), file)
    }
  }
})

# test getNumSearchRuns

test_that("number of search runs is correctly reported", {
  expect_equal(getNumSearchRuns(james, "coconut", "Random Descent"),
               length(getSearchRuns(james, "coconut", "Random Descent")))
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      expect_equal(getNumSearchRuns(james, p, s), 10)
    }
  }
})

# test getBestSolutionValues

test_that("getBestSolutionValues allows to omit problem name in case of a single problem", {
  expect_error(getBestSolutionValues(james, search = "Random Descent"), "more than one problem")
  # now for coconut problem only
  james.coco <- reduceJAMES(james, problems = "coco")
  getBestSolutionValues(james.coco, search = "Random Descent")
})

test_that("getBestSolutionValues allows to omit search name in case of a single search for the given problem", {
  expect_error(getBestSolutionValues(james, problem = "coconut"), "more than one search for problem")
  # now for random descent only
  james.rd <- reduceJAMES(james, searches = "Descent")
  getBestSolutionValues(james.rd, problem = "coconut")
  # omit both problem and search name
  james.rd.coco <- reduceJAMES(james.rd, problems = "coco")
  getBestSolutionValues(james.rd.coco)
})

test_that("getBestSolutionValues returns a numeric vector", {
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      best.values <- getBestSolutionValues(james, p, s)
      expect_true(is.vector(best.values))
      expect_true(is.numeric(best.values))
    }
  }
})

test_that("number of best solution values corresponds to number of runs", {
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      expect_equal(length(getBestSolutionValues(james, p, s)), getNumSearchRuns(james, p, s))
    }
  }
})

test_that("getBestSolutionValues returns the correct values", {
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      file <- sprintf("files/best-solution-values-%s-%s.rds", gsub(" ", "-", p), gsub(" ", "-", s))
      expect_equal_to_reference(getBestSolutionValues(james, p, s), file)
    }
  }
})

# test getBestSolutions

test_that("getBestSolutions allows to omit problem name in case of a single problem", {
  expect_error(getBestSolutions(james, search = "Random Descent"), "more than one problem")
  # now for coconut problem only
  james.coco <- reduceJAMES(james, problems = "coco")
  getBestSolutions(james.coco, search = "Random Descent")
})

test_that("getBestSolutions allows to omit search name in case of a single search for the given problem", {
  expect_error(getBestSolutions(james, problem = "coconut"), "more than one search for problem")
  # now for random descent only
  james.rd <- reduceJAMES(james, searches = "Descent")
  getBestSolutions(james.rd, problem = "coconut")
  # omit both problem and search name
  james.rd.coco <- reduceJAMES(james.rd, problems = "coco")
  getBestSolutions(james.rd.coco)
})

test_that("number of best solutions corresponds to number of runs", {
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      expect_equal(length(getBestSolutions(james, p, s)), getNumSearchRuns(james, p, s))
    }
  }
})

test_that("getBestSolutions returns the correct solutions", {
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      file <- sprintf("files/best-solutions-%s-%s.rds", gsub(" ", "-", p), gsub(" ", "-", s))
      expect_equal_to_reference(getBestSolutions(james, p, s), file)
    }
  }
})

# test getConvergenceTimes

test_that("getConvergenceTimes allows to omit problem name in case of a single problem", {
  expect_error(getConvergenceTimes(james, search = "Random Descent"), "more than one problem")
  # now for coconut problem only
  james.coco <- reduceJAMES(james, problems = "coco")
  getConvergenceTimes(james.coco, search = "Random Descent")
})

test_that("getConvergenceTimes allows to omit search name in case of a single search for the given problem", {
  expect_error(getConvergenceTimes(james, problem = "coconut"), "more than one search for problem")
  # now for random descent only
  james.rd <- reduceJAMES(james, searches = "Descent")
  getConvergenceTimes(james.rd, problem = "coconut")
  # omit both problem and search name
  james.rd.coco <- reduceJAMES(james.rd, problems = "coco")
  getConvergenceTimes(james.rd.coco)
})

test_that("number of convergence times corresponds to number of runs", {
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      expect_equal(length(getConvergenceTimes(james, p, s)), getNumSearchRuns(james, p, s))
    }
  }
})

test_that("getConvergenceTimes complains if convergence ratio is not within [0,1]", {
  expect_error(getConvergenceTimes(james, "coconut", "Random Descent", r = -1), "value in [0,1]", fixed = TRUE)
  expect_error(getConvergenceTimes(james, "coconut", "Random Descent", r = 1.001), "value in [0,1]", fixed = TRUE)
  expect_error(getConvergenceTimes(james, "coconut", "Random Descent", r = "abc"), "value in [0,1]", fixed = TRUE)
})

test_that("getConvergenceTimes returns a numeric vector", {
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      best.values <- getConvergenceTimes(james, p, s)
      expect_true(is.vector(best.values))
      expect_true(is.numeric(best.values))
    }
  }
})

test_that("convergence times are greater than or equal to -1", {
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      conv.times <- getConvergenceTimes(james, p, s)
      expect_false(any(conv.times < -1))
    }
  }
})

test_that("getConvergenceTimes returns the correct values", {
  for(p in getProblems(james)){
    for(s in getSearches(james, p)){
      file <- sprintf("files/convergence-times-%s-%s.rds", gsub(" ", "-", p), gsub(" ", "-", s))
      expect_equal_to_reference(getConvergenceTimes(james, p, s), file)
    }
  }
})











