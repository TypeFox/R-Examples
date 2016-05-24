context("utils-failures")

pitching_info = dget("_pitchingInfo.dat")

pitching_info2 = duplicateSchema(pitching_info)


test_that("exceptions are properly handled by utility functions", {
  
  expect_equal(isTable(NULL, NULL), logical(0))
  
  expect_equal(isTable(NULL, character(0)), logical(0))
  
  expect_null(viewTableSummary())
  
  expect_error(getTableSummary(NULL, tableName="abcdefg", parallel=TRUE),
               "Please register parallel backend appropriate to your platform to run with parallel=TRUE")
  
  expect_error(getTableSummary(NULL, tableName = "abcdefg"),
               "Connection is not valid RODBC object.")
  
  expect_error(getTableSummary(NULL, tableName="fielding", mock=TRUE),
               "Test sql with 'getTableSummary' only for 'batting' or 'pitching' tables.")
  
  expect_error(getTableSummary(NULL, tableName="pitching", except=pitching_info$COLUMN_NAME, mock=TRUE), 
               "No columns specified found in the table 'pitching'")
  
  expect_error(getTableSummary(NULL, tableName="pitching", percentiles=c(TRUE,TRUE), mock=TRUE),
               ".*Percentiles must be either a numeric vector or a logical flag.*")
  
  expect_error(getTableSummary(NULL, tableName="pitching", percentiles="FALSE", mock=TRUE),
               ".*Percentiles must be either a numeric vector or a logical flag.*")
  
  expect_error(getTableSummary(NULL, tableName="pitching", percentiles=-1, mock=TRUE),
               "Invalid percentile value\\(s\\) passed \\(below 0 or above 100\\): -1")
  
})


test_that("getNullCounts throws errors", {
  
  expect_error(getNullCounts(NULL, output='not-an-option'),
               ".*'arg' should be one of \"long\", \"wide\", \"matrix\".*")
  
  expect_error(getNullCounts(NULL),
               "Table name must be specified.")
  
  expect_error(getNullCounts(NULL, NULL, test=TRUE),
               "Must provide tableInfo when test==TRUE.")
  
  expect_error(getNullCounts(NULL, NULL),
               "first argument is not an open RODBC channel")
  
  expect_error(getNullCounts(NULL, 'pitching', pitching_info2),
               "Table name is not uqique - must provide schema using either parameter 'schema' or 'tableName'.")
  
  expect_error(getNullCounts(NULL, "pitching", include=c('lgid','yearid'), tableInfo = pitching_info2, test=TRUE),
               "Table name is not uqique - must provide schema using either parameter 'schema' or 'tableName'.")
  
})


test_that("getNullCounts sql is correct", {
  
  expect_equal_normalized(getNullCounts(NULL, "pitching", include=c('lgid','yearid'), tableInfo = pitching_info, test=TRUE),
                          "SELECT COUNT(1) - COUNT(lgid) AS lgid, COUNT(1) - COUNT(yearid) AS yearid FROM pitching  ")
  
  expect_equal_normalized(getNullCounts(NULL, "pitching", include=c('lgid','yearid'), tableInfo = pitching_info, 
                                        where='yearid > 2000', test=TRUE),
                          "SELECT COUNT(1) - COUNT(lgid) AS lgid, COUNT(1) - COUNT(yearid) AS yearid FROM pitching 
                             WHERE yearid > 2000 ")
  
  expect_equal_normalized(getNullCounts(NULL, "pitching", include=c('lgid','yearid'), tableInfo = pitching_info, 
                                        percent=TRUE, test=TRUE),
                          "SELECT (COUNT(1) - COUNT(lgid) + 1.0)/(COUNT(1) + 1.0) AS lgid, 
                                  (COUNT(1) - COUNT(yearid) + 1.0)/(COUNT(1) + 1.0) AS yearid FROM pitching  ")
  
  expect_equal_normalized(getNullCounts(NULL, "pitching", include=c('lgid','yearid'), tableInfo = pitching_info, 
                                        percent=TRUE, where='yearid > 2000', test=TRUE),
                          "SELECT (COUNT(1) - COUNT(lgid) + 1.0)/(COUNT(1) + 1.0) AS lgid, 
                                  (COUNT(1) - COUNT(yearid) + 1.0)/(COUNT(1) + 1.0) AS yearid FROM pitching
                             WHERE yearid > 2000  ")

})


all_tables = data.frame(TABLE_CAT=c("beehive","beehive","beehive","beehive"),
                        TABLE_SCHEM=c("public","graph","graph","public"),
                        TABLE_NAME=c("batting","edges","vertices","pitching"),
                        TABLE_TYPE=c("TABLE","TABLE","TABLE","TABLE"),
                        stringsAsFactors = FALSE)

test_that("isTable is correct", {
  
  expect_equal(isTable(NULL, c(edges="graph.edges", vertices="graph.vertices"), 
                       allTables = all_tables),
               c(edges=TRUE, vertices=TRUE))
  
  expect_equal(isTable(NULL, c("graph.edges", "graph.vertices", "pitching", "non.existing", "select * from batting"),
                       allTables = all_tables),
               c(graph.edges=TRUE, graph.vertices=TRUE, pitching=TRUE, non.existing=FALSE, `select * from batting`=NA))
  
  expect_equal(isTable(NULL, c('1'="pitching", '2'="batting", '3'="SELECT source, target, COUNT(*) FROM data GROUP BY 1,2"),
                       allTables = all_tables),
               c('1'=TRUE, '2'=TRUE, '3'=NA))
  
})