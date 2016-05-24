context("computeAggregates")

test_that("computeAggregates throws errors", {
  
  expect_error(computeAggregates(channel=NULL), 
               "Must have table name.")
  
  expect_error(computeAggregates(channel=NULL, tableName="fake"),
               "Must have one or more columns.")
  
  expect_error(computeAggregates(channel=NULL, tableName="table_name", by="column1",
                       aggregates = vector()),
               "Must have at least one aggregate defined.")
  
  expect_error(computeAggregates(NULL, tableName="table_name", by="column1"),
               "Connection is not valid RODBC object.")
  
})


test_that("computeAggregates SQL is correct", {
  
  expect_equal_normalized(computeAggregates(channel=NULL, "teams_enh", 
                                  by = c("name || ', ' || park teamname", "lgid", "teamid", "decadeid"),
                                  aggregates = c("min(name) name", "min(park) park", "avg(rank) rank", 
                                                 "avg(attendance) attendance"),
                                  test = TRUE),
                          "SELECT name || ', ' || park teamname, lgid, teamid, decadeid, min(name) name, 
                                  min(park) park, avg(rank) rank, avg(attendance) attendance 
                             FROM teams_enh  
                            GROUP BY name || ', ' || park, lgid, teamid, decadeid"                          
                          )
  
  expect_equal_normalized(computeAggregates(channel=NULL, "teams_enh",
                                  by = c("teamid", "decadeid"),
                                  aggregates = c("min(rank) minrank", "max(rank) maxrank"),
                                  where = "lgid = 'AL'",
                                  test = TRUE),
                          "SELECT teamid, decadeid, min(rank) minrank, max(rank) maxrank
                             FROM teams_enh
                            WHERE lgid = 'AL'
                            GROUP BY teamid, decadeid"
                          )
  
  expect_equal_normalized(computeAggregates(channel=NULL, "pitching_enh",
                                  by = c("teamid", "decadeid"), 
                                  aggregates = c("sum(so) so", 
                                                 "sum(so)/(sum(sum(so)) over (partition by decadeid)) percent"),
                                  where = "decadeid >= 1980",
                                  test = TRUE),
                          "SELECT teamid, decadeid, sum(so) so, 
                                  sum(so)/(sum(sum(so)) over (partition by decadeid)) percent
                             FROM pitching_enh
                            WHERE decadeid >= 1980
                            GROUP BY teamid, decadeid"
                          )
})