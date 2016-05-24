context("gzip utility functions")

test_that('digest(,algo="crc32") and gzip.crc agree',{
  # digest::digest("Hello", algo = 'crc32',serialize = F)
  expect_equal(gzip.crc('testdata/hello.gz'), "f7d18982")
})

test_that('is.gzip works',{
  expect_true(is.gzip(system.file('help/aliases.rds')))
  
  notgzipfile=tempfile()
  expect_true(is.na(is.gzip(notgzipfile)))
  writeLines('not a gzip', notgzipfile)
  expect_false(is.gzip(notgzipfile))
  
  expect_warning(rval<-gzip.crc(notgzipfile))
  expect_true(is.na(rval))
  # make an empty file
  writeBin(logical(), notgzipfile)
  expect_warning(gzip.crc(notgzipfile))
  
  con=gzfile(gzipfile<-tempfile(),open='wt')
  writeLines('This one is gzipped', con)
  close(con)
  expect_true(is.gzip(gzipfile))
  
  unlink(c(notgzipfile,gzipfile))
})
