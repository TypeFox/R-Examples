context("computeCorrelations")

pitching_info = dget("_pitchingInfo.dat")

test_that("computeCorrelations throws errors", {
  
  expect_error(computeCorrelations(test=TRUE),
               "Must provide tableInfo when test==TRUE")
  
  expect_error(computeCorrelations(NULL),
               "Connection is not valid RODBC object.")
  
  expect_error(computeCorrelations(channel=NULL, tableName="pitching_enh", tableInfo=pitching_info,
                                   include = c('abcdef', 'qwerty'), test=TRUE),
               "Must provide at least 2 numeric columns.")
  
})

test_that("computeCorrelations SQL is correct", {
  
  expect_equal_normalized(
    computeCorrelations(channel=conn, tableName="pitching_enh", tableInfo=pitching_info,
                        include = c('w', 'era'), test=TRUE),
    "SELECT * FROM corr_reduce(
                  ON corr_map(
                    ON ( SELECT w, era FROM pitching_enh  )
                    columnpairs( 'era:w')
                    key_name('key')
                  )
                  partition by key )")
  
  expect_equal_normalized(
    computeCorrelations(channel=conn, tableName="pitching_enh", tableInfo=pitching_info,
                        include = c('w','l','cg','sho','sv','ipouts','h','er','hr','bb','so','baopp',
                                    'era','whip','ktobb','fip'),
                        where = "decadeid = 2000", test=TRUE),
    "SELECT * FROM corr_reduce(
                  ON corr_map(
                    ON ( SELECT w, l, cg, sho, sv, ipouts, h, er, hr, bb, so, baopp, era, whip, ktobb, fip 
                           FROM pitching_enh WHERE decadeid = 2000 )
                    columnpairs( 'l:w', 'cg:w', 'sho:w', 'sv:w', 'ipouts:w', 'h:w', 'er:w', 'hr:w', 'bb:w', 'so:w', 'baopp:w', 'era:w', 'ktobb:w', 'fip:w', 
                                 'cg:l', 'ipouts:l', 'h:l', 'er:l', 'hr:l', 'bb:l', 'baopp:l', 'era:l', 'ktobb:l', 'fip:l', 'bb:cg', 'baopp:cg', 'l:sho', 
                                 'cg:sho', 'ipouts:sho', 'h:sho', 'er:sho', 'hr:sho', 'bb:sho', 'baopp:sho', 'era:sho', 'ktobb:sho', 'fip:sho', 'l:sv', 
                                 'cg:sv', 'sho:sv', 'ipouts:sv', 'h:sv', 'er:sv', 'hr:sv', 'bb:sv', 'so:sv', 'baopp:sv', 'era:sv', 'ktobb:sv', 'fip:sv', 
                                 'cg:ipouts', 'h:ipouts', 'er:ipouts', 'hr:ipouts', 'bb:ipouts', 'baopp:ipouts', 'era:ipouts', 'fip:ipouts', 'cg:h', 'er:h', 
                                 'bb:h', 'baopp:h', 'era:h', 'fip:h', 'cg:er', 'bb:er', 'baopp:er', 'cg:hr', 'h:hr', 'er:hr', 'bb:hr', 'baopp:hr', 'era:hr', 
                                 'fip:hr', 'baopp:bb', 'l:so', 'cg:so', 'sho:so', 'ipouts:so', 'h:so', 'er:so', 'hr:so', 'bb:so', 'baopp:so', 'era:so', 
                                 'ktobb:so', 'fip:so', 'cg:era', 'er:era', 'bb:era', 'baopp:era', 'w:whip', 'l:whip', 'cg:whip', 'sho:whip', 'sv:whip', 
                                 'ipouts:whip', 'h:whip', 'er:whip', 'hr:whip', 'bb:whip', 'so:whip', 'baopp:whip', 'era:whip', 'ktobb:whip', 'fip:whip', 
                                 'cg:ktobb', 'ipouts:ktobb', 'h:ktobb', 'er:ktobb', 'hr:ktobb', 'bb:ktobb', 'baopp:ktobb', 'era:ktobb', 'fip:ktobb', 'cg:fip', 
                                 'er:fip', 'bb:fip', 'baopp:fip', 'era:fip')
                    key_name('key')
                  )
                  partition by key )",
    "Error in pitching_enh correlations", "compute correlations in pitching_enh SQL")
  
  expect_equal_normalized(
    computeCorrelations(channel=conn, tableName="pitching_enh", tableInfo=pitching_info,
                        include=getNumericColumns(pitching_info), 
                        except=c("yearid", "decadeid", "sh", "sf", "gidp"), test=TRUE),
    "SELECT * FROM corr_reduce(
                  ON corr_map(
                    ON ( SELECT w, l, g, gs, cg, sho, sv, ipouts, h, er, hr, bb, so, baopp, era, ibb, wp, 
                                hbp, bk, bfp, gf, r, whip, ktobb, fip 
                           FROM pitching_enh  )
                    columnpairs( 'l:w', 'g:w', 'gs:w', 'cg:w', 'sho:w', 'sv:w', 'ipouts:w', 'h:w', 'er:w', 'hr:w', 'bb:w', 'so:w', 'baopp:w', 'era:w', 'ibb:w', 
                                 'hbp:w', 'bk:w', 'bfp:w', 'gf:w', 'r:w', 'ktobb:w', 'fip:w', 'g:l', 'gs:l', 'cg:l', 'ipouts:l', 'h:l', 'er:l', 'hr:l', 'bb:l', 
                                 'baopp:l', 'era:l', 'ibb:l', 'hbp:l', 'bk:l', 'bfp:l', 'gf:l', 'ktobb:l', 'fip:l', 'cg:g', 'er:g', 'bb:g', 'baopp:g', 'era:g', 
                                 'bk:g', 'bfp:g', 'fip:g', 'g:gs', 'cg:gs', 'er:gs', 'bb:gs', 'baopp:gs', 'era:gs', 'bk:gs', 'bfp:gs', 'gf:gs', 'fip:gs', 'bb:cg', 
                                 'baopp:cg', 'bk:cg', 'bfp:cg', 'l:sho', 'g:sho', 'gs:sho', 'cg:sho', 'ipouts:sho', 'h:sho', 'er:sho', 'hr:sho', 'bb:sho', 
                                 'baopp:sho', 'era:sho', 'ibb:sho', 'hbp:sho', 'bk:sho', 'bfp:sho', 'gf:sho', 'r:sho', 'ktobb:sho', 'fip:sho', 'l:sv', 'g:sv', 
                                 'gs:sv', 'cg:sv', 'sho:sv', 'ipouts:sv', 'h:sv', 'er:sv', 'hr:sv', 'bb:sv', 'so:sv', 'baopp:sv', 'era:sv', 'ibb:sv', 'hbp:sv', 
                                 'bk:sv', 'bfp:sv', 'gf:sv', 'r:sv', 'ktobb:sv', 'fip:sv', 'g:ipouts', 'gs:ipouts', 'cg:ipouts', 'h:ipouts', 'er:ipouts', 'hr:ipouts',
                                 'bb:ipouts', 'baopp:ipouts', 'era:ipouts', 'ibb:ipouts', 'hbp:ipouts', 'bk:ipouts', 'bfp:ipouts', 'gf:ipouts', 'fip:ipouts', 'g:h', 
                                 'gs:h', 'cg:h', 'er:h', 'bb:h', 'baopp:h', 'era:h', 'bk:h', 'bfp:h', 'gf:h', 'fip:h', 'cg:er', 'bb:er', 'baopp:er', 'bk:er', 
                                 'bfp:er', 'g:hr', 'gs:hr', 'cg:hr', 'h:hr', 'er:hr', 'bb:hr', 'baopp:hr', 'era:hr', 'hbp:hr', 'bk:hr', 'bfp:hr', 'gf:hr', 'fip:hr', 
                                 'baopp:bb', 'l:so', 'g:so', 'gs:so', 'cg:so', 'sho:so', 'ipouts:so', 'h:so', 'er:so', 'hr:so', 'bb:so', 'baopp:so', 'era:so', 
                                 'ibb:so', 'hbp:so', 'bk:so', 'bfp:so', 'gf:so', 'r:so', 'ktobb:so', 'fip:so', 'cg:era', 'er:era', 'bb:era', 'baopp:era', 'bk:era', 
                                 'bfp:era', 'g:ibb', 'gs:ibb', 'cg:ibb', 'h:ibb', 'er:ibb', 'hr:ibb', 'bb:ibb', 'baopp:ibb', 'era:ibb', 'hbp:ibb', 'bk:ibb', 'bfp:ibb',
                                 'gf:ibb', 'fip:ibb', 'w:wp', 'l:wp', 'g:wp', 'gs:wp', 'cg:wp', 'sho:wp', 'sv:wp', 'ipouts:wp', 'h:wp', 'er:wp', 'hr:wp', 'bb:wp', 
                                 'so:wp', 'baopp:wp', 'era:wp', 'ibb:wp', 'hbp:wp', 'bk:wp', 'bfp:wp', 'gf:wp', 'r:wp', 'whip:wp', 'ktobb:wp', 'fip:wp', 'g:hbp', 
                                 'gs:hbp', 'cg:hbp', 'h:hbp', 'er:hbp', 'bb:hbp', 'baopp:hbp', 'era:hbp', 'bk:hbp', 'bfp:hbp', 'gf:hbp', 'fip:hbp', 'bb:bk', 
                                 'baopp:bk', 'bfp:bk', 'bb:bfp', 'baopp:bfp', 'g:gf', 'cg:gf', 'er:gf', 'bb:gf', 'baopp:gf', 'era:gf', 'bk:gf', 'bfp:gf', 'fip:gf', 
                                 'l:r', 'g:r', 'gs:r', 'cg:r', 'ipouts:r', 'h:r', 'er:r', 'hr:r', 'bb:r', 'baopp:r', 'era:r', 'ibb:r', 'hbp:r', 'bk:r', 'bfp:r', 
                                 'gf:r', 'ktobb:r', 'fip:r', 'w:whip', 'l:whip', 'g:whip', 'gs:whip', 'cg:whip', 'sho:whip', 'sv:whip', 'ipouts:whip', 'h:whip', 
                                 'er:whip', 'hr:whip', 'bb:whip', 'so:whip', 'baopp:whip', 'era:whip', 'ibb:whip', 'hbp:whip', 'bk:whip', 'bfp:whip', 'gf:whip', 
                                 'r:whip', 'ktobb:whip', 'fip:whip', 'g:ktobb', 'gs:ktobb', 'cg:ktobb', 'ipouts:ktobb', 'h:ktobb', 'er:ktobb', 'hr:ktobb', 
                                 'bb:ktobb', 'baopp:ktobb', 'era:ktobb', 'ibb:ktobb', 'hbp:ktobb', 'bk:ktobb', 'bfp:ktobb', 'gf:ktobb', 'fip:ktobb', 'cg:fip', 
                                 'er:fip', 'bb:fip', 'baopp:fip', 'era:fip', 'bk:fip', 'bfp:fip')
                     key_name('key')
                   )
    partition by key )")
  
})