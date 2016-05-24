context("showData")

pitching_info = dget("_pitchingInfo.dat")
batting_info = dget("_battingInfo.dat")

test_that("format 'boxplot' works", {
  p = ggplot_build(showData(tableName='pitching', tableInfo=pitching_info, format='boxplot'))
  
  notNullNumData = pitching_info[pitching_info$not_null_count > 0 & 
                                   pitching_info$COLUMN_NAME %in% getNumericColumns(pitching_info),]
  expect_equal(nrow(p$data[[1]]), length(getNumericColumns(notNullNumData)))
  expect_equal(setdiff(p$panel$ranges[[1]]$x.labels, getNumericColumns(notNullNumData)),
               character(0))
})

test_that("format 'boxplot' with facets works", {
  
  cols = c('bb','er','era','so','r')
  p = ggplot_build(showData(tableName='pitching', tableInfo=pitching_info, format='boxplot',
               include=cols, facet=TRUE))
  expect_equal(nrow(p$panel$layout), length(cols))
  expect_equal(nrow(p$data[[1]]), length(cols))
  
})

test_that("showData issues warnings", {
  
  # format='boxplot'
  expect_warning(showData(tableInfo=pitching_info, type='character', format='boxplot', test=TRUE),
                 "Automatically regressing to numerical types")
  
  # format='corr'
  expect_warning(showData(tableInfo=pitching_info, type='character', format='corr', tableName='pitching',
                          test=TRUE),
                 "Ignoring non-numeric types")
})


test_that("showData throws errors", {
  
  expect_error(showData(channel=NULL, tableName="fake", format='boxplot', test=TRUE),
               "Must provide tableInfo when test==TRUE")
})

test_that("showData format 'boxplot' throws errors", {
  
  # format='boxplot'
  expect_error(showData(tableInfo=batting_info, type='numeric', format='boxplot', include='nothing'),
               "Not all specified columns are in the table summary")

  expect_error(showData(tableInfo=batting_info, type='numeric', format='boxplot', 
                        include=c('ab','g','lgid','playerid'), except=c('ab','g')),
               "Nothing to show: check lists of columns")
  
  expect_error(showData(tableInfo=batting_info, type='numeric', format='boxplot', include=c('lgid','playerid')),
               "Nothing to show: check lists of columns")
})

test_that("showData format 'histogram' throws errors", {
  
  # format='histogram'
  expect_error(showData(tableInfo=batting_info, type='numeric', format='histogram'),
               "Must provide table name")
  expect_error(showData(tableInfo=batting_info, type='character', format='histogram', tableName='batting'),
               "Factor histograms are not supported")
  expect_error(showData(tableInfo=batting_info, type='temporal', format='histogram', tableName='batting'),
               "Datetime histograms are not supported")
})

test_that("showData format 'histogram' works", {
  
  expect_equal_normalized(
    showData(tableInfo=batting_info, tableName="batting", type='numeric', format='histogram', test=TRUE),
    "SELECT *
       FROM hist_reduce(
         ON hist_map(
         ON (SELECT cast(yearid as numeric) yearid FROM batting ) as data_input 
         PARTITION BY ANY  
         binsize('0.4')
         startvalue('2000')
         endvalue('2012')
         VALUE_COLUMN('yearid') )
       partition by  1 )")
})
  
  
test_that("showData format 'scatterplot' throws errors", {
  
  # format='scatterplot'
  expect_error(showData(tableInfo=batting_info, format='scatterplot'),
               "Must provide table name")
  expect_error(showData(tableInfo=batting_info, format='scatterplot', tableName='batting'),
               "define x and y coordiantes")
  expect_error(showData(tableInfo=batting_info, format='scatterplot', tableName='batting',
                        sampleSize=1, include=c('ba')), 
               "define x and y coordiantes")
  expect_error(showData(tableInfo=batting_info, format='scatterplot', tableName='batting',
                        include=c('lgid','playerid')),
               "Scatterplot format is valid for numerical data only.")
  expect_error(showData(tableInfo=batting_info, format='scatterplot', tableName='batting',
                        include=c('so','ba')), 
               "Sample fraction or sample size must be specified.")
  expect_error(showData(tableInfo=batting_info, format='scatterplot', tableName='batting',
                        sampleSize=1, include=c('lgid','playerid')), 
               "numerical data only")
  expect_error(showData(tableInfo=batting_info[batting_info$COLUMN_NAME %in% c('lgid')], format='scatterplot', 
                        tableName='batting', sampleSize=1, include=c('lgid','playerid')),
               "Not all specified columns are in the table summary")
})  

test_that("showData format 'scatterplot' works", {
  
  expect_equal_normalized(
    showData(tableInfo=batting_info, format='scatterplot', tableName="batting", 
             include=c('so','ba'), sampleSize=10, test=TRUE),
    "SELECT * 
       FROM sample(
         ON (SELECT so, ba FROM batting  )
           AS DATA PARTITION BY ANY
         ON (SELECT COUNT(*) as stratum_count FROM batting ) 
           AS SUMMARY DIMENSION
         ApproximateSampleSize('10'))")
  
  expect_equal_normalized(
    showData(tableInfo=batting_info, format='scatterplot', tableName="batting", 
             include=c('so','ba'), sampleFraction=0.10, test=TRUE),
    "SELECT *   
       FROM sample(
         ON (SELECT so, ba FROM batting  )
         SampleFraction('0.1'))")
  
  expect_equal_normalized(
    showData(tableInfo=pitching_info, format='scatterplot', tableName="pitching_enh",
             include=c('ktobb', 'fip', 'teamid','yearid'),
             sampleSize=100, facetName=c('teamid','yearid'), regressionLine=TRUE,
             where="yearid in (2010,2011,2012) and teamid in ('TEX','NYA')",
             test=TRUE),
    "SELECT * 
       FROM sample(
         ON (SELECT teamid, yearid, ktobb, fip FROM pitching_enh
              WHERE yearid in (2010,2011,2012) and teamid in ('TEX','NYA') 
            ) AS DATA PARTITION BY ANY
         ON (SELECT COUNT(*) as stratum_count FROM pitching_enh 
              WHERE yearid in (2010,2011,2012) and teamid in ('TEX','NYA') 
            ) AS SUMMARY DIMENSION
         ApproximateSampleSize('100'))"
    )
})
  

test_that("showData format 'corr' throws errors", {
  
  # format='corr'
  expect_error(showData(tableInfo=batting_info, format='corr'),
               "Must provide table name")
  expect_error(showData(tableInfo=batting_info, format='corr', tableName='batting',
                        include=c('lgid','playerid')),
               "Nothing to show: check lists of columns")
})

test_that("showData format 'corr' works", {
  expect_equal_normalized(
    showData(tableInfo=pitching_info, type='numeric', format='corr', tableName="pitching",
             include=c('era','ipouts','bb','so'), test=TRUE),
    "SELECT * FROM corr_reduce(
       ON corr_map(
         ON ( SELECT ipouts, bb, so, era FROM pitching  )
         columnpairs( 'bb:ipouts', 'era:ipouts', 'ipouts:so', 'bb:so', 'era:so', 'bb:era')
         key_name('key')
       )
       partition by key
     )")
})
  