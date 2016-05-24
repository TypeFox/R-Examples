context("computeSample")

test_that("computeSample throws errors", {
  
  expect_error(computeSample(channel=NULL),
               "Table name must be specified.")
  
  expect_error(computeSample(channel=NULL, tableName="faketable"),
               "Sample fraction or sample size must be specified.")
  
  expect_error(computeSample(channel=NULL, tableName="faketable", sampleFraction=0.1,
                             except = c("some", "columns", "defined", "without", "include"),
                             test=FALSE),
               "Can't test without include but with except.")
  
  expect_error(computeSample(channel=NULL, tableName="faketable", sampleFraction=1.000001,
                             test=FALSE),
               "sampleFraction <= 1 is not TRUE")
  
  expect_error(computeSample(channel=NULL, tableName="faketable", sampleFraction=-0.000001,
                             test=FALSE),
               "sampleFraction >= 0 is not TRUE")
  
  expect_error(computeSample(channel=NULL, tableName="faketable", sampleFraction=0.000001),
               "Connection is not valid RODBC object.")
})

test_that("computeSample SQL is correct", {
  
  expect_equal_normalized(computeSample(channel=NULL, tableName='public.batting_enh', 
                                        sampleFraction=0.1, test=TRUE),
                          "SELECT * FROM sample(
                                           ON (SELECT  *  FROM public.batting_enh )  
                                           SampleFraction('0.1'))")
  
  expect_equal_normalized(computeSample(channel=NULL, tableName='public.batting_enh', 
                                        sampleFraction=0.1, sampleSize=10, test=TRUE),
                          "SELECT * FROM sample(
                                           ON (SELECT  *  FROM public.batting_enh )  
                                           SampleFraction('0.1'))")
  
  expect_equal_normalized(computeSample(channel=NULL, tableName='public.batting_enh', 
                                        sampleFraction=0.1, where="lgid='AL'", test=TRUE),
                          "SELECT * FROM sample(
                                           ON (SELECT  *  FROM public.batting_enh WHERE lgid='AL' )  
                                           SampleFraction('0.1'))")
  
  expect_equal_normalized(computeSample(channel=NULL, tableName='public.batting_enh',
                                        sampleSize=1000, test=TRUE),
                          "SELECT * 
                             FROM sample(
                               ON (SELECT  *  FROM public.batting_enh  ) 
                                 AS DATA PARTITION BY ANY                           
                               ON (SELECT COUNT(*) as stratum_count FROM public.batting_enh )  
                                 AS SUMMARY DIMENSION
                               ApproximateSampleSize('1000'))")
  
  expect_equal_normalized(computeSample(channel=NULL, tableName='public.batting_enh',
                                        sampleSize=1000, where="lgid='NL'", test=TRUE),
                          "SELECT * 
                             FROM sample(
                               ON (SELECT  *  FROM public.batting_enh WHERE lgid='NL' )
                                 AS DATA PARTITION BY ANY                           
                               ON (SELECT COUNT(*) as stratum_count FROM public.batting_enh WHERE lgid='NL' )  
                                 AS SUMMARY DIMENSION
                               ApproximateSampleSize('1000'))")
})