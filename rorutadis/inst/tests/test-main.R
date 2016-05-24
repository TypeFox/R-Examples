context("Main functionality")

test_that("calculateExtremeClassCardinalities", {
  perf <- matrix(c(200,100,400,200,100,500,400,400,0.2,0.3,0.1,0.4,0.6,0.4,0.5,0.6), 8)
  nrClasses = 4
  strictVF = FALSE
  criteria = c("g", "g")
  characteristicPoints = c(3, 2)
  
  problem <- buildProblem(perf, nrClasses, strictVF, criteria, characteristicPoints)
  
  problem <- addAssignmentsLB(problem, c(5, 3))
  problem <- addAssignmentsUB(problem, c(5, 3))
  problem <- addAssignmentsLB(problem, c(6, 2))
  
  
  problem <- addMinimalClassCardinalities(problem, c(1, 2))
  problem <- addMinimalClassCardinalities(problem, c(4, 1))
  problem <- addMaximalClassCardinalities(problem, c(1, 6))
  
  problem <- addAssignmentPairwiseAtLeastComparisons(problem, c(5, 3, 1))
  problem <- addAssignmentPairwiseAtLeastComparisons(problem, c(7, 1, 2))
  problem <- addAssignmentPairwiseAtLeastComparisons(problem, c(3, 1, 1))
  
  problem <- addAssignmentPairwiseAtMostComparisons(problem, c(7, 4, 2))
  
  calculateAssignments(problem, FALSE)
  
  expect_true(all.equal(calculateExtremeClassCardinalities(problem),  matrix(data = c(2,1,1,1,3,2,4,4), ncol=2)))
})