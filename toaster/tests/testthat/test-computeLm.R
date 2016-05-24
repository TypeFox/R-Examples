context("computeLm")

batting_info = dget("_battingInfo.dat")
non_compliant_row = batting_info[1,]
non_compliant_row[1,c('COLUMN_NAME', 'TYPE_NAME')] = c('non_compliant', 'varchar')
batting_info2 = rbind(batting_info, non_compliant_row)
no_values_row = batting_info[1,]
no_values_row[1,c('COLUMN_NAME', 'TYPE_NAME')] = c('no_values', 'varchar')
batting_info2 = rbind(batting_info2, no_values_row)
teamid_row = batting_info[1,]
teamid_row[1,c('COLUMN_NAME', 'TYPE_NAME')] = c('teamid', 'varchar')
batting_info2 = rbind(batting_info2, teamid_row)
decadeid_row = batting_info[1,]
decadeid_row[1, c('COLUMN_NAME', 'TYPE_NAME')] = c('decadeid', 'integer')
batting_info2 = rbind(batting_info2, decadeid_row)

test_that("computeLm throws errors", {
  
  expect_error(computeLm(), 
               "Must provide connection.")
  
  expect_error(computeLm(channel=NULL), 
               "Must provide table and expression.")
  
  expect_error(computeLm(channel=NULL, tableName="fake_name"),
               "Must provide table and expression.")
  
  expect_error(computeLm(channel=NULL, tableName="fake_name", formula=y~1 , test=FALSE),
               "No predictors found in formula.")
  
  expect_error(computeLm(channel=NULL, tableName="fake_name", formula=~x, test=FALSE),
               "No predictors found in formula.")
  
  expect_error(computeLm(channel=NULL, tableName="fake_name", formula=y~x, test=TRUE),
               "Must provide tableInfo when test==TRUE.")
  
  expect_error(computeLm(channel=NULL, tableName="fake_name", formula=lgid~ba, tableInfo=batting_info, 
                         test=TRUE),
               "Response variable's column is not of numeric type.")
  
  expect_error(computeLm(channel=NULL, tableName="batting", formula=ba~rbi+stat2+stat3, tableInfo=batting_info, 
                         test=TRUE),
               "Columns stat2, stat3 are not found in table batting.")
  
  expect_error(computeLm(channel=NULL, tableName="fake_name", formula=ba~no_values+rbi, 
                         tableInfo=batting_info2, test=TRUE),
               "Categorical column 'no_values' has 1 or no values. Consider removing or replacing it.")
  
  expect_error(computeLm(channel=NULL, tableName="fake_name", formula=ba~non_compliant+rbi+so, 
                         tableInfo=batting_info2, test=TRUE),
               "Categorical column 'non_compliant' values are not valid strings. Consider replacing them with alpha-numeric versions.")
})

test_that("computeLm SQL is correct", {
  
  expect_equal_normalized(
    computeLm(channel=NULL, tableName="batting", formula = ba ~ rbi, tableInfo=batting_info, test=TRUE),
    "SELECT * 
       FROM linreg(
         ON linregmatrix(
           ON (SELECT rbi x1, ba y FROM batting ) 
         )
         PARTITION BY 1
       )")
  
  expect_equal_normalized(
    computeLm(channel=NULL, tableName="batting", formula = ba ~ rbi + bb + so, tableInfo=batting_info, 
              test=TRUE),
    "SELECT * 
       FROM linreg(
         ON linregmatrix(
           ON (SELECT rbi x1, bb x2, so x3, ba y FROM batting ) 
         )
         PARTITION BY 1
       )")
  
  expect_equal_normalized(
    computeLm(channel=NULL, tableName="batting", formula = ba ~ rbi + bb + so, tableInfo=batting_info, 
                  where = "lgid = 'AL'", test=TRUE),
    "SELECT * 
       FROM linreg(
         ON linregmatrix(
           ON (SELECT rbi x1, bb x2, so x3, ba y FROM batting WHERE lgid = 'AL' ) 
         )
         PARTITION BY 1
       )")
})

test_that("computeLm with categorical predictors SQL is correct", {
  
  expect_equal_normalized(
    computeLm(channel=NULL, tableName="batting", formula = ba ~ rbi + bb + lgid, tableInfo=batting_info,
              test=TRUE),
    "SELECT *
       FROM linreg(
         ON linregmatrix(
           ON (SELECT rbi x1, bb x2, CASE WHEN lgid = 'NL' THEN 1 ELSE 0 END x3, ba y FROM batting )
         )
       PARTITION BY 1
     )")
  
  expect_equal_normalized(
    computeLm(channel=NULL, tableName="batting", formula = ba ~ rbi + lgid + teamid, tableInfo=batting_info2,
              test=TRUE),
    "SELECT * 
       FROM linreg(
         ON linregmatrix(
           ON (SELECT rbi x1, CASE WHEN lgid = 'NL' THEN 1 ELSE 0 END x2, CASE WHEN teamid = 'NYY' THEN 1 ELSE 0 END x3, 
                      CASE WHEN teamid = 'TEX' THEN 1 ELSE 0 END x4, CASE WHEN teamid = 'TOR' THEN 1 ELSE 0 END x5, ba y  
                 FROM batting )
           )
         PARTITION BY 1
       )")
  
  expect_equal_normalized(
    computeLm(channel=NULL, tableName="batting", formula = ba ~ lgid + teamid, tableInfo=batting_info2,
              where="decadeid >= 1980", test=TRUE),
    "SELECT *
       FROM linreg(
         ON linregmatrix(
         ON (SELECT CASE WHEN lgid = 'NL' THEN 1 ELSE 0 END x1, CASE WHEN teamid = 'NYY' THEN 1 ELSE 0 END x2, 
                    CASE WHEN teamid = 'TEX' THEN 1 ELSE 0 END x3, CASE WHEN teamid = 'TOR' THEN 1 ELSE 0 END x4, ba y  
               FROM batting WHERE decadeid >= 1980  )
         )
         PARTITION BY 1
       )")
  
  expect_equal_normalized(
    computeLm(channel=NULL, tableName="batting", formula = ba ~ rbi + lgid + decadeid, 
              categories="decadeid", tableInfo=batting_info2,
              test=TRUE),
    "SELECT *
       FROM linreg(
         ON linregmatrix(
         ON (SELECT rbi x1, CASE WHEN lgid = 'NL' THEN 1 ELSE 0 END x2, CASE WHEN decadeid = '2000' THEN 1 ELSE 0 END x3, 
                      CASE WHEN decadeid = '2010' THEN 1 ELSE 0 END x4, ba y  
                 FROM batting )
           )
         PARTITION BY 1
       )")
  
})