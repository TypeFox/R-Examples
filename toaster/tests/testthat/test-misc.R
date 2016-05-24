context("misc-tests")

batting_info = dget("_battingInfo.dat")
pitching_info = dget("_pitchingInfo.dat")

test_that("type functions work", {
  expect_equal(getNumericTypes(), c('integer',
                                    'numeric',
                                    'bigint',
                                    'smallint',
                                    'real',
                                    'double precision',
                                    'serial',
                                    'bigserial',
                                    'float',
                                    'decimal'))
  
  expect_equal(getCharacterTypes(), c('varchar',
                                      'char',
                                      'character'))
  
  expect_equal(getTemporalTypes(), c('date', 
                                     'timestamp without time zone', 
                                     'timestamp with time zone',
                                     'time without time zone',
                                     'time with time zone'))
  
  expect_true(all(c('yearid','decadeid','w','g','gs','cg','sho','sv','ipouts','h','er',
                    'sh','sf','gidp','whip','ktobb','fip') %in% getNumericColumns(pitching_info)))
  
  expect_true(all(c('lgid','teamid','playerid') %in% getCharacterColumns(pitching_info)))
  
  expect_equal(getTemporalColumns(pitching_info), character(0))
})


test_that("temporary table name function", {
  
  expect_error(makeTempTableName("$"), "Prefix may contain alpha-numeric characters only")
  
  expect_error(makeTempTableName("tooobig", 60), "Too long prefix: 63 characters is Aster limit on table name length")
  
  expect_error(makeTempTableName("prefix", 10, "$%^^schema01"), "Schema may contain alpha-numeric characters only")
  
  expect_equal(nchar(makeTempTableName()), 28)
  
  expect_equal(nchar(makeTempTableName(NULL, 25)), 33)
  
  expect_equal(nchar(makeTempTableName('a', 25)), 35)
  
  expect_equal(charmatch("toa_temp_", makeTempTableName()), 1)
  
  expect_equal(charmatch("toa_temp_scaled_", makeTempTableName(prefix="scaled")), 1)
  
  expect_equal(charmatch("myschema.toa_temp_scaled_", makeTempTableName(prefix="scaled", schema="myschema")), 1)
  
})