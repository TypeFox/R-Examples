context("computePercentiles")

test_that("computePercentiles throws errros", {
  
  expect_error(computePercentiles(), 
               "Must provide connection.")
  
  expect_error(computePercentiles(channel=NULL),
               "Must provide table name.")
  
  expect_error(computePercentiles(channel=NULL, tableName=NULL),
               "Must provide table name.")
  
  expect_error(computePercentiles(channel=NULL, tableName="fake_table"),
               "Must provide at least one column name.")
  
  expect_error(computePercentiles(channel=NULL, tableName="fake_table", columns=NULL),
               "Must provide at least one column name.")
  
  expect_error(computePercentiles(channel=NULL, tableName="fake_table", columns=NULL),
               "Must provide at least one column name.")
  
  expect_error(computePercentiles(channel=NULL, tableName="fake_table", columns=character(0)),
               "Must provide at least one column name.")
  
  expect_error(computePercentiles(channel=NULL, tableName="fake_table", columns="column1"),
               "Connection is not valid RODBC object.")
})

test_that("computePercentiles numerical SQL is correct", {
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="pitching", columns="ipouts",
                       test=TRUE),
    "SELECT * FROM approxPercentileReduce(
       ON (
         SELECT * FROM approxPercentileMap(
           ON  ( SELECT * FROM pitching )                       
           TARGET_COLUMN( 'ipouts' )
           ERROR( 1 )                       
       ) )
       PARTITION BY  1                   
       PERCENTILE( 0,5,10,25,50,75,90,95,100 )
     )")
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="pitching", columns="ipouts",
                       test=TRUE),
    "SELECT * FROM approxPercentileReduce(
       ON (
         SELECT * FROM approxPercentileMap(
           ON  ( SELECT * FROM pitching )                       
           TARGET_COLUMN( 'ipouts' )
           ERROR( 1 )                       
       ) )
       PARTITION BY  1                   
       PERCENTILE( 0,5,10,25,50,75,90,95,100 )
     )")
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="pitching", columns="ipouts",
                       by="lgid", test=TRUE),
    "SELECT * FROM approxPercentileReduce(
       ON (
         SELECT * FROM approxPercentileMap(
           ON  ( SELECT * FROM pitching )                       
           TARGET_COLUMN( 'ipouts' )
           ERROR( 1 )  
           GROUP_COLUMNS( 'lgid' )
       ) )
       PARTITION BY lgid                  
       PERCENTILE( 0,5,10,25,50,75,90,95,100 ) 
       GROUP_COLUMNS( 'lgid' )                  
    )")
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="pitching", columns="ipouts",
                       percentiles = c(1,2,3,4),
                       test=TRUE),
    "SELECT * FROM approxPercentileReduce(
      ON (
        SELECT * FROM approxPercentileMap(
          ON  ( SELECT * FROM pitching )                       
          TARGET_COLUMN( 'ipouts' )
          ERROR( 1 )                       
        ) )
      PARTITION BY  1                   
      PERCENTILE( 0,1,2,3,4,25,50,75,100 )
    )")
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="pitching", columns="ipouts",
                       percentiles = c(10,90), where="lgid = 'NL'",
                       test=TRUE),
    "SELECT * FROM approxPercentileReduce(
      ON (
        SELECT * FROM approxPercentileMap(
          ON  ( SELECT * FROM pitching WHERE lgid = 'NL' )                       
          TARGET_COLUMN( 'ipouts' )
          ERROR( 1 )                       
        ) )
      PARTITION BY  1                   
      PERCENTILE( 0,10,25,50,75,90,100 )
    )")
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="pitching_enh", columns=c("ipouts","era"),
                       by=c("lgid","decadeid"), where="lgid IN ('NL','AL')",
                       test=TRUE),
    "SELECT * FROM approxPercentileReduce(
       ON (
         SELECT * FROM approxPercentileMap(
           ON  ( SELECT * FROM pitching_enh WHERE lgid IN ('NL','AL') )                       
           TARGET_COLUMN( 'ipouts' )
           ERROR( 1 )  
           GROUP_COLUMNS( 'lgid', 'decadeid' )
       ) )
       PARTITION BY lgid, decadeid                  
       PERCENTILE( 0,5,10,25,50,75,90,95,100 ) 
       GROUP_COLUMNS( 'lgid', 'decadeid' )                  
    )")
})

test_that("computePercentiles temporal SQL is correct", {
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="public.master_enh", columns="birthdate",
                       temporal=TRUE, test=TRUE),
    "SELECT percentile, MAX(birthdate)::varchar value, MAX(EXTRACT('EPOCH' FROM birthdate)) epoch FROM 
       (SELECT birthdate, ntile(100) OVER ( ORDER BY birthdate) percentile
          FROM public.master_enh WHERE birthdate IS NOT NULL) t
      WHERE percentile IN ( 5,10,25,50,75,90,95,100 )
      GROUP BY 1 
      ORDER BY 1")
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="public.master_enh", columns=c("birthdate", "deathdate"),
                       temporal=TRUE, by=c("bats", "throws"), test=TRUE),
    "SELECT bats, throws, percentile, MAX(birthdate)::varchar value, MAX(EXTRACT('EPOCH' FROM birthdate)) epoch FROM 
       (SELECT bats, throws, birthdate, ntile(100) OVER (PARTITION BY bats, throws ORDER BY birthdate) percentile
          FROM public.master_enh WHERE birthdate IS NOT NULL) t
      WHERE percentile IN ( 5,10,25,50,75,90,95,100 )
      GROUP BY 1, 2, 3
      ORDER BY 1, 2, 3")
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="public.master_enh", columns="birthdate",
                       temporal=TRUE, percentiles=c(0,5,10,25,50,75,90,95,100), test=TRUE),
    "SELECT percentile, MAX(birthdate)::varchar value, MAX(EXTRACT('EPOCH' FROM birthdate)) epoch FROM 
       (SELECT birthdate, ntile(100) OVER ( ORDER BY birthdate) percentile
          FROM public.master_enh WHERE birthdate IS NOT NULL) t
      WHERE percentile IN ( 0,5,10,25,50,75,90,95,100 )
      GROUP BY 1 
     UNION 
     SELECT 0, MIN(birthdate)::varchar, MIN(EXTRACT('EPOCH' FROM birthdate)) epoch
       FROM public.master_enh WHERE birthdate IS NOT NULL
      ORDER BY 1")
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="public.master_enh", columns=c("birthdate", "deathdate"),
                       temporal=TRUE, percentiles=c(0,5,10,25,50,75,90,95,100), 
                       by=c("bats", "throws"), test=TRUE),
    "SELECT bats, throws, percentile, MAX(birthdate)::varchar value, MAX(EXTRACT('EPOCH' FROM birthdate)) epoch FROM 
       (SELECT bats, throws, birthdate, ntile(100) OVER (PARTITION BY bats, throws ORDER BY birthdate) percentile
          FROM public.master_enh WHERE birthdate IS NOT NULL) t
      WHERE percentile IN ( 0,5,10,25,50,75,90,95,100 )
      GROUP BY 1, 2, 3
     UNION
     SELECT bats, throws, 0, MIN(birthdate)::varchar, MIN(EXTRACT('EPOCH' FROM birthdate)) epoch FROM public.master_enh WHERE birthdate IS NOT NULL
      GROUP BY 1, 2
      ORDER BY 1, 2, 3")
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="public.master_enh", columns="birthdate",
                       temporal=TRUE, percentiles=c(0,5,10,25,50,75,90,95,100), 
                       where="birthcountry = 'USA'", test=TRUE),
    "SELECT percentile, MAX(birthdate)::varchar value, MAX(EXTRACT('EPOCH' FROM birthdate)) epoch FROM 
       (SELECT birthdate, ntile(100) OVER ( ORDER BY birthdate) percentile
          FROM public.master_enh WHERE birthcountry = 'USA' AND birthdate IS NOT NULL) t
      WHERE percentile IN ( 0,5,10,25,50,75,90,95,100 )
      GROUP BY 1 
     UNION 
     SELECT 0, MIN(birthdate)::varchar, MIN(EXTRACT('EPOCH' FROM birthdate)) epoch
       FROM public.master_enh WHERE birthcountry = 'USA' AND birthdate IS NOT NULL
      ORDER BY 1")
  
  expect_equal_normalized(
    computePercentiles(channel=NULL, tableName="public.master_enh", columns=c("birthdate", "deathdate"),
                       temporal=TRUE, percentiles=c(0,5,10,25,50,75,90,95,100), 
                       by=c("bats", "throws"), where="birthcountry = 'USA'", test=TRUE),
    "SELECT bats, throws, percentile, MAX(birthdate)::varchar value, MAX(EXTRACT('EPOCH' FROM birthdate)) epoch FROM 
       (SELECT bats, throws, birthdate, ntile(100) OVER (PARTITION BY bats, throws ORDER BY birthdate) percentile
          FROM public.master_enh WHERE birthcountry = 'USA' AND birthdate IS NOT NULL) t
      WHERE percentile IN ( 0,5,10,25,50,75,90,95,100 )
      GROUP BY 1, 2, 3
     UNION
     SELECT bats, throws, 0, MIN(birthdate)::varchar, MIN(EXTRACT('EPOCH' FROM birthdate)) epoch 
       FROM public.master_enh WHERE birthcountry = 'USA' AND birthdate IS NOT NULL
      GROUP BY 1, 2
      ORDER BY 1, 2, 3")
  
})