context("showGraph")

test_that("showGraph throws errors", {
  
  expect_error(showGraph(channel=NULL, toaGraph(), test=TRUE),
               "Must provide allTables when test==TRUE.")
  
  expect_error(showGraph(channel=NULL, toaGraph("vs","edges"), 
                         allTables = data.frame(TABLE_NAME=c(""), stringsAsFactors = FALSE), 
                         test=TRUE),
               ".*Both vertices and edges must exist as tables or views.*")
})