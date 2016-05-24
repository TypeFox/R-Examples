context('binary')

test_that('twoalt', {
  data(WorldEvents)
  scores <- calcscore(answer ~ forecast, fam="beta",
                    param=c(1,1), data=WorldEvents,
                    bounds=c(0,1))
  expect_is(scores, 'numeric')
  expect_false(any(is.na(scores)))

  scores2 <- calcscore(answer ~ forecast, fam="pow",
                       param=2, data=WorldEvents,
                       bounds=c(0,1))
  expect_is(scores2, 'numeric')
  expect_false(any(is.na(scores2)))

  scores.man <- with(WorldEvents, (forecast - answer)^2)
  expect_equal(scores, scores.man)
  expect_equal(scores2, scores.man)

  scores3 <- calcscore(answer ~ forecast, fam="sph",
                       param=2, data=WorldEvents,
                       bounds=c(0,1))
  expect_is(scores3, 'numeric')
  expect_false(any(is.na(scores3)))

  scores.man2 <- 1 - with(WorldEvents, (forecast * answer + (1 - forecast)*(1 - answer))/
                      sqrt(forecast^2 + (1-forecast)^2))
  expect_equal(scores3, scores.man2)

  scores.brier <- brierscore(answer ~ forecast, data=WorldEvents)
  expect_is(scores.brier, 'numeric')
  expect_equal(scores2, scores.brier)

  scores.grp.brier <- brierscore(answer ~ forecast, data=WorldEvents, group="forecast")$mnbrier
  expect_equal(scores.grp.brier, tapply(scores2, WorldEvents$forecast, mean))
  
})
