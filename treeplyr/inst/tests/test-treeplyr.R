context("treeplyr functions work")
test_that("treedata can handle matrix/dataframe input", {
  require(testthat)
  data(anolis)
  td <- make.treedata(anolis$phy, anolis$dat, name_column=1)
  originaldat <- anolis$dat[,-1]
  rownames(originaldat) <- anolis$dat[,1]
  jacknife <- sample(1:nrow(anolis$dat), 10, replace=FALSE)
  td_filtered <- filter(td, SVL > 3.5, island=="Cuba")
  td_filtered_ <- filter_(td, "SVL > 3.5", "island=='Cuba'")
  td_selected <- select(td, SVL, ecomorph, island)
  td_selected_ <- select_(td, "SVL", "ecomorph", "island")
  td_mutated <- mutate(td, logSVL = log(SVL), discrete_awesome = as.numeric(awesomeness>0))
  td_mutated_ <- mutate_(td, logSVL = "log(SVL)", discrete_awesome = "as.numeric(awesomeness>0)")
  td_reorder <- reorder(td, "postorder")
  td_treeply <- treeply(td, rescale, model="OU", 10)
  
  ##Make sure that names don't get mixed up from original data
  expect_identical(td$dat$SVL[match(anolis$phy$tip.label[jacknife], attributes(td)$tip.label)], originaldat[anolis$phy$tip.label[jacknife], "SVL"])
  expect_identical(td_filtered$dat$SVL, originaldat[td_filtered$phy$tip.label, "SVL"])
  expect_identical(td_selected$dat$SVL[match(anolis$phy$tip.label[jacknife], attributes(td)$tip.label)], originaldat[anolis$phy$tip.label[jacknife], "SVL"])
  expect_identical(td_mutated$dat$SVL[match(anolis$phy$tip.label[jacknife], attributes(td)$tip.label)], originaldat[anolis$phy$tip.label[jacknife], "SVL"])
  expect_identical(td_reorder$dat$SVL[match(anolis$phy$tip.label[jacknife], attributes(td)$tip.label)], originaldat[anolis$phy$tip.label[jacknife], "SVL"])
  expect_identical(td_treeply$dat$SVL[match(anolis$phy$tip.label[jacknife], attributes(td)$tip.label)], originaldat[anolis$phy$tip.label[jacknife], "SVL"])
  
  ##Make sure that lazy and standard eval both produce the same results
  expect_identical(td_selected, td_selected_)
  expect_identical(td_mutated, td_mutated_)
  expect_identical(td_filtered, td_filtered_)
  
  ##Make sure that filter works right
  expect_true(min(filter(td, SVL > 3.5)$dat$SVL) > 3.5)
  expect_identical(td_filtered$phy$tip.label, td_filtered$phy$tip.label[as.numeric(rownames(td_filtered$dat))])
  ##Make sure that treeply applies the function correctly
  expect_identical(td_treeply$phy$edge.length, rescale(td$phy, "OU", 10)$edge.length)  

  
})
