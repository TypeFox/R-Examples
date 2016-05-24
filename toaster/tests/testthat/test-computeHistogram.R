context("computeHistogram")

pitching_info = dget("_pitchingInfo.dat")
start_end_test_info = dget("_startEndTestInfo.dat")

test_that("computeHistogram throws errors", {
  
  expect_error(computeHistogram(binMethod = 'not a method'),
               ".*'arg' should be one of \"manual\", \"Sturges\", \"Scott\"")
  
  expect_error(computeHistogram(channel=NULL),
               "Must provide table and column names.")
  
  expect_error(computeHistogram(channel=NULL, tableName="fake"),
               "Must provide table and column names.")
  
   expect_error(computeHistogram(channel=NULL, tableName="fake", columnName="nocolumn"),
               "Connection is not valid RODBC object.")
  
  expect_error(computeHistogram(channel=NULL, tableName="fake", columnName="nocolumn", test=TRUE),
               "Must provide tableInfo when test==TRUE")
  
  expect_error(computeHistogram(channel=NULL, tableName="pitching", columnName="nocolumn", tableInfo=pitching_info, test=TRUE),
               "No columns specified found in the table")
  
  expect_error(computeHistogram(channel=NULL, tableName="texts.containertrailerplanpaths", columnName="missentcount", 
                                tableInfo=start_end_test_info, test=TRUE),
               "Start value should not be greater than or equal to end value. Try to run with useIQR=FALSE or check that data is not constant.")
  
  expect_error(computeHistogram(channel=NULL, tableName="texts.containertrailerplanpaths", columnName="missentcount", 
                                tableInfo=start_end_test_info, useIQR=TRUE, test=TRUE),
               "Start value should not be greater than or equal to end value. Try to run with useIQR=FALSE or check that data is not constant.")
  
  
})


test_that("computeHistogram SQL is correct", {
  
  expect_equal_normalized(computeHistogram(channel=NULL, tableName="pitching", columnName="h",
                                           tableInfo=pitching_info, binMethod="manual", 
                                           test=TRUE),
                          "SELECT * FROM hist_reduce(
                                      ON hist_map(
                                           ON (SELECT   cast(h as numeric) h FROM pitching ) as data_input PARTITION BY ANY  
                                           binsize('6.1')
                                           startvalue('0')
                                           endvalue('183')    
                                           VALUE_COLUMN('h')     
                                      ) 
                                      partition by  1 )"
    )
  
  expect_equal_normalized(computeHistogram(channel=NULL, tableName="pitching", columnName="h",
                                             tableInfo=pitching_info, binMethod="manual", by="lgid",
                                             test=TRUE),
                          "SELECT * FROM hist_reduce(
                                      ON hist_map(
                                           ON (SELECT lgid, cast(h as numeric) h FROM pitching ) as data_input PARTITION BY ANY  
                                           binsize('6.1')
                                           startvalue('0')
                                           endvalue('183')    
                                           VALUE_COLUMN('h')     
                                           GROUP_COLUMNS('lgid')
                                      ) 
                                      partition by  lgid)"
  )
  
  expect_equal_normalized(computeHistogram(channel=NULL, tableName="pitching", columnName="h",
                                             tableInfo=pitching_info, binMethod="manual", by=c("lgid","teamid"),
                                             test=TRUE),
                          "SELECT * FROM hist_reduce(
                                      ON hist_map(
                                           ON (SELECT lgid, teamid, cast(h as numeric) h FROM pitching ) as data_input PARTITION BY ANY  
                                           binsize('6.1')
                                           startvalue('0')
                                           endvalue('183')    
                                           VALUE_COLUMN('h')     
                                           GROUP_COLUMNS('lgid', 'teamid')
                                      ) 
                                      partition by  lgid, teamid)"
  )
  
  expect_equal_normalized(computeHistogram(channel=NULL, tableName="pitching", columnName="h",
                                             tableInfo=pitching_info, binMethod="Sturges",
                                             test=TRUE),
                          "SELECT * FROM hist_reduce(
                                      ON hist_map(
                                           ON (SELECT cast(h as numeric) h FROM pitching ) as data_input PARTITION BY ANY
                                           ON hist_prep( ON (SELECT cast(h as numeric) h FROM pitching ) VALUE_COLUMN('h') ) as data_stat DIMENSION
                                           BIN_SELECT('Sturges')
                                           VALUE_COLUMN('h')
                                      )
                                      partition by 1 )"
  )
  
  expect_equal_normalized(computeHistogram(channel=NULL, tableName="pitching", columnName="h",
                                             tableInfo=pitching_info, binMethod="Sturges", by=c("lgid","teamid"),
                                             test=TRUE),
                          "SELECT * FROM hist_reduce(
                                      ON hist_map(
                                           ON (SELECT lgid, teamid, cast(h as numeric) h FROM pitching ) as data_input PARTITION BY ANY
                                           ON hist_prep( ON (SELECT cast(h as numeric) h FROM pitching ) VALUE_COLUMN('h') ) as data_stat DIMENSION
                                           BIN_SELECT('Sturges')
                                           VALUE_COLUMN('h')
                                           GROUP_COLUMNS('lgid', 'teamid')
                                      )
                                      partition by lgid, teamid)"
  )
  
  expect_equal_normalized(computeHistogram(channel=NULL, tableName="pitching", columnName="h",
                                             tableInfo=pitching_info, binMethod="Scott", by=c("lgid", "teamid"),
                                             where="yearid >= 1980",
                                             test=TRUE),
                          "SELECT * FROM hist_reduce(
                                      ON hist_map(
                                           ON (SELECT lgid, teamid, cast(h as numeric) h FROM pitching WHERE yearid >= 1980 ) as data_input PARTITION BY ANY
                                           ON hist_prep( ON (SELECT cast(h as numeric) h FROM pitching WHERE yearid >= 1980 ) VALUE_COLUMN('h') ) as data_stat DIMENSION
                                           BIN_SELECT('Scott')
                                           VALUE_COLUMN('h')
                                           GROUP_COLUMNS('lgid', 'teamid')
                                      )
                                      partition by lgid, teamid)"
  )
  
  expect_equal_normalized(computeHistogram(channel=NULL, tableName="texts.containertrailerplanpaths", columnName="missentcount", 
                                           tableInfo=start_end_test_info, useIQR=FALSE, test=TRUE),
                          "SELECT * FROM hist_reduce(
                             ON hist_map(
                               ON (SELECT   cast(missentcount as numeric) missentcount FROM texts.containertrailerplanpaths   ) as data_input PARTITION BY ANY  
                               binsize('1.6')
                               startvalue('0')
                               endvalue('48')    
                               VALUE_COLUMN('missentcount')     
                             ) 
                             partition by  1 )")
  
})

test_that("computeHistogram with barplot SQL is correct", {
  
  expect_equal_normalized(computeHistogram(channel=NULL, tableName="pitching", columnName="playerid",
                                           tableInfo=pitching_info, by="lgid", where="decadeid > 1900",
                                           test=TRUE),
                          "SELECT playerid, lgid, COUNT(*) cnt FROM pitching 
                            WHERE decadeid > 1900
                            GROUP BY playerid, lgid "
                          )
})


#test_that("computeHistogram with column frequency SQL is correct", {
#  
#  expect_equal_normalized(computeHistogram(channel=NULL, tableName="pitching", ))
#})