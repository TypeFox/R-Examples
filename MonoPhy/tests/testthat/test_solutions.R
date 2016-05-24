test_that("Regular solution is accurate",{
  data(Ericactree)
  controlsolution <- AssessMonophyly(Ericactree)
  #  test solution object
  expect_equal(length(controlsolution[[1]]), 7)
	expect_is(controlsolution, "list")
  # Test summary table
  expect_equal(dim(controlsolution[[1]]$summary)[1], 6)
  expect_equal(dim(controlsolution[[1]]$summary)[2], 2)
  #  Test summary entries
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[1,1]))[controlsolution[[1]]$summary[1,1]], 32)
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[2,1]))[controlsolution[[1]]$summary[2,1]], 9)
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[3,1]))[controlsolution[[1]]$summary[3,1]], 6)
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[4,1]))[controlsolution[[1]]$summary[4,1]], 17)
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[5,1]))[controlsolution[[1]]$summary[5,1]], 7)
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[6,1]))[controlsolution[[1]]$summary[6,1]], 1)
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[1,2]))[controlsolution[[1]]$summary[1,2]], 77)
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[2,2]))[controlsolution[[1]]$summary[2,2]], 25)
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[3,2]))[controlsolution[[1]]$summary[3,2]], 35)
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[4,2]))[controlsolution[[1]]$summary[4,2]], 17)
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[5,2]))[controlsolution[[1]]$summary[5,2]], 11)
  expect_equal(as.numeric(levels(controlsolution[[1]]$summary[6,2]))[controlsolution[[1]]$summary[6,2]], 2)
  #  Test result table
  expect_equal(dim(controlsolution[[1]]$result)[1], 32)
  expect_equal(dim(controlsolution[[1]]$result)[2], 8)
  #  Test lenghts of taxa and tip lists
  expect_equal(length(controlsolution[[1]]$IntruderTaxa), 6)
  expect_equal(length(controlsolution[[1]]$IntruderTips), 6)
  expect_equal(length(controlsolution[[1]]$OutlierTaxa), 1)
  expect_equal(length(controlsolution[[1]]$OutlierTips), 1)
  #  Test Tipstates
  expect_equal(dim(controlsolution[[1]]$TipStates)[1], 77)
  expect_equal(dim(controlsolution[[1]]$TipStates)[2], 3)
})

test_that("Multifurcating solution1 is accurate",{
  data(Ericactree)
  library(spacodiR)
  multifurcEric1 <- phy.deresolve(Ericactree, time.range = c(50, 25), relative = FALSE)
  solution1 <- AssessMonophyly(multifurcEric1)
  #  test solution object
  expect_equal(length(solution1[[1]]), 7)
	expect_is(solution1, "list")
  # Test summary table
  expect_equal(dim(solution1[[1]]$summary)[1], 6)
  expect_equal(dim(solution1[[1]]$summary)[2], 2)
  #  Test summary entries
  expect_equal(as.numeric(levels(solution1[[1]]$summary[1,1]))[solution1[[1]]$summary[1,1]], 32)
  expect_equal(as.numeric(levels(solution1[[1]]$summary[2,1]))[solution1[[1]]$summary[2,1]], 8)
  expect_equal(as.numeric(levels(solution1[[1]]$summary[3,1]))[solution1[[1]]$summary[3,1]], 7)
  expect_equal(as.numeric(levels(solution1[[1]]$summary[4,1]))[solution1[[1]]$summary[4,1]], 17)
  expect_equal(as.numeric(levels(solution1[[1]]$summary[5,1]))[solution1[[1]]$summary[5,1]], 25)
  expect_equal(as.numeric(levels(solution1[[1]]$summary[6,1]))[solution1[[1]]$summary[6,1]], 0)
  expect_equal(as.numeric(levels(solution1[[1]]$summary[1,2]))[solution1[[1]]$summary[1,2]], 77)
  expect_equal(as.numeric(levels(solution1[[1]]$summary[2,2]))[solution1[[1]]$summary[2,2]], 23)
  expect_equal(as.numeric(levels(solution1[[1]]$summary[3,2]))[solution1[[1]]$summary[3,2]], 37)
  expect_equal(as.numeric(levels(solution1[[1]]$summary[4,2]))[solution1[[1]]$summary[4,2]], 17)
  expect_equal(as.numeric(levels(solution1[[1]]$summary[5,2]))[solution1[[1]]$summary[5,2]], 47)
  expect_equal(as.numeric(levels(solution1[[1]]$summary[6,2]))[solution1[[1]]$summary[6,2]], 0)
  #  Test result table
  expect_equal(dim(solution1[[1]]$result)[1], 32)
  expect_equal(dim(solution1[[1]]$result)[2], 8)
  #  Test lenghts of taxa and tip lists
  expect_equal(length(solution1[[1]]$IntruderTaxa), 7)
  expect_equal(length(solution1[[1]]$IntruderTips), 7)
  expect_equal(length(solution1[[1]]$OutlierTaxa), 0)
  expect_equal(length(solution1[[1]]$OutlierTips), 0)
  #  Test Tipstates
  expect_equal(dim(solution1[[1]]$TipStates)[1], 77)
  expect_equal(dim(solution1[[1]]$TipStates)[2], 3)
})

test_that("Multifurcating solution2 is accurate",{
  data(Ericactree)
  library(spacodiR)
  multifurcEric2 <- phy.deresolve(Ericactree, time.range = c(75, 70), relative = FALSE)
  solution2 <- AssessMonophyly(multifurcEric2)
  #  test solution object
  expect_equal(length(solution2[[1]]), 7)
	expect_is(solution2, "list")
  # Test summary table
  expect_equal(dim(solution2[[1]]$summary)[1], 6)
  expect_equal(dim(solution2[[1]]$summary)[2], 2)
  #  Test summary entries
  expect_equal(as.numeric(levels(solution2[[1]]$summary[1,1]))[solution2[[1]]$summary[1,1]], 32)
  expect_equal(as.numeric(levels(solution2[[1]]$summary[2,1]))[solution2[[1]]$summary[2,1]], 9)
  expect_equal(as.numeric(levels(solution2[[1]]$summary[3,1]))[solution2[[1]]$summary[3,1]], 6)
  expect_equal(as.numeric(levels(solution2[[1]]$summary[4,1]))[solution2[[1]]$summary[4,1]], 17)
  expect_equal(as.numeric(levels(solution2[[1]]$summary[5,1]))[solution2[[1]]$summary[5,1]], 7)
  expect_equal(as.numeric(levels(solution2[[1]]$summary[6,1]))[solution2[[1]]$summary[6,1]], 1)
  expect_equal(as.numeric(levels(solution2[[1]]$summary[1,2]))[solution2[[1]]$summary[1,2]], 77)
  expect_equal(as.numeric(levels(solution2[[1]]$summary[2,2]))[solution2[[1]]$summary[2,2]], 25)
  expect_equal(as.numeric(levels(solution2[[1]]$summary[3,2]))[solution2[[1]]$summary[3,2]], 35)
  expect_equal(as.numeric(levels(solution2[[1]]$summary[4,2]))[solution2[[1]]$summary[4,2]], 17)
  expect_equal(as.numeric(levels(solution2[[1]]$summary[5,2]))[solution2[[1]]$summary[5,2]], 11)
  expect_equal(as.numeric(levels(solution2[[1]]$summary[6,2]))[solution2[[1]]$summary[6,2]], 2)
  #  Test result table
  expect_equal(dim(solution2[[1]]$result)[1], 32)
  expect_equal(dim(solution2[[1]]$result)[2], 8)
  #  Test lenghts of taxa and tip lists
  expect_equal(length(solution2[[1]]$IntruderTaxa), 6)
  expect_equal(length(solution2[[1]]$IntruderTips), 6)
  expect_equal(length(solution2[[1]]$OutlierTaxa), 1)
  expect_equal(length(solution2[[1]]$OutlierTips), 1)
  #  Test Tipstates
  expect_equal(dim(solution2[[1]]$TipStates)[1], 77)
  expect_equal(dim(solution2[[1]]$TipStates)[2], 3)
})

test_that("Multifurcating solution3 is accurate",{
  data(Ericactree)
  library(spacodiR)
  multifurcEric3 <- phy.deresolve(Ericactree, time.range = c(20, 40), relative = FALSE)
  solution3 <- AssessMonophyly(multifurcEric3)
  #  test solution object
  expect_equal(length(solution3[[1]]), 7)
	expect_is(solution3, "list")
  # Test summary table
  expect_equal(dim(solution3[[1]]$summary)[1], 6)
  expect_equal(dim(solution3[[1]]$summary)[2], 2)
  #  Test summary entries
  expect_equal(as.numeric(levels(solution3[[1]]$summary[1,1]))[solution3[[1]]$summary[1,1]], 32)
  expect_equal(as.numeric(levels(solution3[[1]]$summary[2,1]))[solution3[[1]]$summary[2,1]], 7)
  expect_equal(as.numeric(levels(solution3[[1]]$summary[3,1]))[solution3[[1]]$summary[3,1]], 8)
  expect_equal(as.numeric(levels(solution3[[1]]$summary[4,1]))[solution3[[1]]$summary[4,1]], 17)
  expect_equal(as.numeric(levels(solution3[[1]]$summary[5,1]))[solution3[[1]]$summary[5,1]], 15)
  expect_equal(as.numeric(levels(solution3[[1]]$summary[6,1]))[solution3[[1]]$summary[6,1]], 1)
  expect_equal(as.numeric(levels(solution3[[1]]$summary[1,2]))[solution3[[1]]$summary[1,2]], 77)
  expect_equal(as.numeric(levels(solution3[[1]]$summary[2,2]))[solution3[[1]]$summary[2,2]], 20)
  expect_equal(as.numeric(levels(solution3[[1]]$summary[3,2]))[solution3[[1]]$summary[3,2]], 40)
  expect_equal(as.numeric(levels(solution3[[1]]$summary[4,2]))[solution3[[1]]$summary[4,2]], 17)
  expect_equal(as.numeric(levels(solution3[[1]]$summary[5,2]))[solution3[[1]]$summary[5,2]], 35)
  expect_equal(as.numeric(levels(solution3[[1]]$summary[6,2]))[solution3[[1]]$summary[6,2]], 2)
  #  Test result table
  expect_equal(dim(solution3[[1]]$result)[1], 32)
  expect_equal(dim(solution3[[1]]$result)[2], 8)
  #  Test lenghts of taxa and tip lists
  expect_equal(length(solution3[[1]]$IntruderTaxa), 8)
  expect_equal(length(solution3[[1]]$IntruderTips), 8)
  expect_equal(length(solution3[[1]]$OutlierTaxa), 1)
  expect_equal(length(solution3[[1]]$OutlierTips), 1)
  #  Test Tipstates
  expect_equal(dim(solution3[[1]]$TipStates)[1], 77)
  expect_equal(dim(solution3[[1]]$TipStates)[2], 3)
})
