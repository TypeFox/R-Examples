context("createBoxplot")

test_that("createBoxplot with single box works", {
  # data = computePercentiles(conn, "batting_enh", columns="ba")
  data = structure(list(percentile = c(0L, 5L, 10L, 25L, 50L, 75L, 90L, 
                                       95L, 100L), 
                        value = c(0, 0, 0, 0.147540983606557, 0.230769230769231, 
                                  0.273684210526316, 0.30873786407767, 0.333907056798623, 1), 
                        column = c("ba", "ba", "ba", "ba", "ba", "ba", "ba", "ba", "ba")), 
                   row.names = c(NA, 9L), .Names = c("percentile", "value", "column"), class = "data.frame")
  
  p = createBoxplot(data, fill=NULL, title="BA Boxplot", useIQR=TRUE)
  expect_is(p,"ggplot")
  expect_equal(nrow(p$data), length(unique(data$column)))
  expect_equal(length(p$layers), length(unique(data$column)))
  expect_equal(p$labels$title, "BA Boxplot")
  expect_true(p$facet$shrink)
  expect_equal(length(p$facet), 1)
  expect_true(p$data$upper_bound < 1, info="When useIQR=T upper_bound is less than 1")
  
  p = createBoxplot(data, fill=NULL, title="BA Boxplot", useIQR=FALSE)
  expect_equal(p$data$upper_bound, 1, info="when useIQR=F upper_bound is 1")
  
})

test_that("createBoxplot with facet wrap works", {
  # data = computePercentiles(conn, "batting_enh", columns = c("ba", "slg", "ta"),
  #                           by=c('lgid'), where="yearid >= 1980")
  data = structure(list(lgid = c("NL", "NL", "NL", "NL", "NL", "NL", "NL", 
                                 "NL", "NL", "AL", "AL", "AL", "AL", "AL", "AL", "AL", "AL", "AL", 
                                 "NL", "NL", "NL", "NL", "NL", "NL", "NL", "NL", "NL", "AL", "AL", 
                                 "AL", "AL", "AL", "AL", "AL", "AL", "AL", "AL", "AL", "AL", "AL", 
                                 "AL", "AL", "AL", "AL", "AL", "NL", "NL", "NL", "NL", "NL", "NL", 
                                 "NL", "NL", "NL"), 
                        percentile = c(0L, 5L, 10L, 25L, 50L, 75L,
                                       90L, 95L, 100L, 0L, 5L, 10L, 25L, 50L, 75L, 90L, 95L, 100L, 0L, 
                                       5L, 10L, 25L, 50L, 75L, 90L, 95L, 100L, 0L, 5L, 10L, 25L, 50L, 
                                       75L, 90L, 95L, 100L, 0L, 5L, 10L, 25L, 50L, 75L, 90L, 95L, 100L, 
                                       0L, 5L, 10L, 25L, 50L, 75L, 90L, 95L, 100L), 
                        value = c(0, 0, 0, 0.107142857142857, 0.217391304347826, 0.267818574514039, 0.300518134715026, 
                                  0.333333333333333, 0.8, 0, 0, 0, 0.194444444444444, 0.245901639344262, 
                                  0.277777777777778, 0.305128205128205, 0.333333333333333, 1, 0, 
                                  0, 0, 0.125, 0.3, 0.402597402597403, 0.486862442040185, 0.524714828897338, 
                                  2, 0, 0, 0, 0.257142857142857, 0.36038961038961, 0.434146341463415, 
                                  0.5, 0.536178107606679, 4, 0, 0, 0, 0.419047619047619, 0.59375, 
                                  0.729957805907173, 0.859756097560976, 0.984168865435356, 9, 0, 
                                  0, 0, 0.2, 0.5, 0.682539682539683, 0.836283185840708, 1, 11), 
                        column = c("ba", "ba", "ba", "ba", "ba", "ba", "ba", "ba", 
                                   "ba", "ba", "ba", "ba", "ba", "ba", "ba", "ba", "ba", "ba", 
                                   "slg", "slg", "slg", "slg", "slg", "slg", "slg", "slg", "slg", 
                                   "slg", "slg", "slg", "slg", "slg", "slg", "slg", "slg", "slg", 
                                   "ta", "ta", "ta", "ta", "ta", "ta", "ta", "ta", "ta", "ta", 
                                   "ta", "ta", "ta", "ta", "ta", "ta", "ta", "ta")), 
                   row.names = c(NA, 54L), .Names = c("lgid", "percentile", "value", "column"), 
                   class = "data.frame")
  
  p = createBoxplot(data, x='column', facet=c('lgid'), useIQR=TRUE, coordFlip=TRUE, 
                    title="Batting by Leagues and Decades", legendPosition="none")
  expect_is(p,"ggplot")
  expect_equal(nrow(p$data), length(unique(data$column)) * length(unique(data$lgid)))
  expect_equal(length(p$facet$facets), 1)
  expect_equal(names(p$facet$facets), "lgid")
  expect_false(any(p$data$upper_bound == c(1.0,0.8,4.0,2.0,9.0,11.0)), info="when useIQR=T upper bound is different from max values")
  expect_equal(p$theme$legend.position, "none")
  expect_is(p$coordinates, c("CoordFlip","CoordCartesian","Coord"))
  expect_equal(p$labels$title, "Batting by Leagues and Decades")
  
  p = createBoxplot(data, x='column', facet=c('lgid'), useIQR=FALSE, coordFlip=TRUE, 
                    title="Batting by Leagues and Decades", legendPosition="none")
  expect_equal(p$data$upper_bound, c(1.0,0.8,4.0,2.0,9.0,11.0), info="when useIQR=F upper_bound is 1")
  
})




