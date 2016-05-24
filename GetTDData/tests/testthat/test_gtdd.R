library(testthat)
library(GetTDData)

dl.folder <- 'TD Data'

my.flag <- download.TD.data(asset.codes = 'LTN_2015', dl.folder = dl.folder)
test_that(desc = 'Test of download function',{
          expect_equal(my.flag , TRUE) } )


returned.rows <- nrow(read.TD.files(maturity = '010116', dl.folder = dl.folder ))

test_that(desc = 'Test of read function',{
  expect_equal(returned.rows>0, TRUE)
  } )

cat('\nDeleting test folder')
unlink(dl.folder, recursive = T)

