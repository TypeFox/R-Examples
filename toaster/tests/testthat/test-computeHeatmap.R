context("computeHeatmap")

test_that("computeHeatmap throws errors", {
  
  expect_error(computeHeatmap(channel=NULL),
               "Must have table name.")
  
  expect_error(computeHeatmap(channel=NULL, table='fake', dimension2='dim2'),
               "Must have all 2 heatmap dimensions defined to compute.")
  
  expect_error(computeHeatmap(channel=NULL, table='fake', dimension1='dim1', aggregates='count(*)'),
               "Must have all 2 heatmap dimensions defined to compute.")
  
  expect_error(computeHeatmap(channel=NULL, table='fake', dimension1='dim1', dimension2='dim2',
                              aggregateFun=c("COUNT(*)", "COUNT(*)/(sum(count(*)) over ())"), 
                              aggregateAlias="alias"),
               '(Defunct; last used in version 0.2.4)')

  expect_error(computeHeatmap(channel=NULL, table='fake', dimension1='dim1', dimension2='dim2',
                              aggregates=vector()),
               "Must have at least one aggregate defined.")
  
  expect_error(computeHeatmap(channel=NULL, table='fake', dimension1='dim1', dimension2='dim2'),
               "Connection is not valid RODBC object.")
  
})

test_that("computeHeatmap SQL is correct", {
  
  expect_equal_normalized(computeHeatmap(channel=NULL, tableName="teams_enh", 
                                         dimension1='franchid', dimension2='decadeid', 
                                         aggregates='avg(w) w', where="decadeid >= 1950",
                                         test=TRUE),
                          "SELECT franchid, decadeid, avg(w) w
                             FROM teams_enh
                            WHERE decadeid >= 1950
                           GROUP BY 1, 2"
  )
  
  expect_equal_normalized(computeHeatmap(channel=NULL, tableName="teams_enh", 
                                         dimension1='franchid', dimension2='decadeid', 
                                         aggregates='avg(w) w', 
                                         where="decadeid >= 1950",
                                         test=TRUE),
                          "SELECT franchid, decadeid, avg(w) w
                             FROM teams_enh
                            WHERE decadeid >= 1950
                           GROUP BY 1, 2"
  )
  
  expect_equal_normalized(computeHeatmap(channel=NULL, tableName="teams_enh",
                                         dimension1='franchid', dimension2='decadeid', 
                                         aggregates=c('avg(w-l) wl', 'avg(r) r', 'avg(h) h'), 
                                         where="decadeid >= 1950",
                                         test=TRUE),
                          "SELECT franchid, decadeid, avg(w-l) wl, avg(r) r, avg(h) h
                             FROM teams_enh
                            WHERE decadeid >= 1950
                            GROUP BY 1, 2"
  )
  
  expect_equal_normalized(computeHeatmap(channel=NULL, tableName="teams_enh",
                                         dimension1='franchid', dimension2='decadeid', 
                                         aggregates=c('avg(w-l) wl', 'avg(r) r', 'avg(h) h'), 
                                         where="decadeid >= 1950",
                                         test=TRUE),
                          "SELECT franchid, decadeid, avg(w-l) wl, avg(r) r, avg(h) h
                             FROM teams_enh
                            WHERE decadeid >= 1950
                            GROUP BY 1, 2"
  )
  
  expect_equal_normalized(computeHeatmap(channel=NULL, tableName="teams_enh", 
                                         dimension1='franchid', dimension2='decadeid', 
                                         aggregates='avg(w) w', by='lgid',
                                         test=TRUE),
                          "SELECT lgid, franchid, decadeid, avg(w) w
                             FROM teams_enh
                            GROUP BY 1, 2, 3"
  )
  
  expect_equal_normalized(computeHeatmap(channel=NULL, tableName="teams_enh", 
                                         dimension1='franchid', dimension2='decadeid', 
                                         aggregates='avg(w) w',
                                         by='lgid',
                                         test=TRUE),
                          "SELECT lgid, franchid, decadeid, avg(w) w
                             FROM teams_enh
                            GROUP BY 1, 2, 3"
  )
  
  expect_equal_normalized(computeHeatmap(channel=NULL, tableName="teams_enh", 
                                         dimension1='franchid', dimension2='decadeid', 
                                         aggregates='avg(w) w', where="decadeid >= 1950", 
                                         by='lgid', test=TRUE),
                          "SELECT lgid, franchid, decadeid, avg(w) w
                             FROM teams_enh
                            WHERE decadeid >= 1950
                            GROUP BY 1, 2, 3"
  )
  
  expect_equal_normalized(computeHeatmap(channel=NULL, tableName="teams_enh", 
                                         dimension1='franchid', dimension2='decadeid', 
                                         aggregates='avg(w) w', 
                                         where="decadeid >= 1950", by='lgid',
                                         test=TRUE),
                          "SELECT lgid, franchid, decadeid, avg(w) w
                             FROM teams_enh
                            WHERE decadeid >= 1950
                            GROUP BY 1, 2, 3"
  )
  
})