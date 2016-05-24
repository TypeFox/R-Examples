#  Tests for the neotoma package.  Mostly validating that changes to the functions
#  do not break the requirements for data formatting.

library("testthat")
library("neotoma")

context('get_contact work as expected')

block_one <- function(){
  marion <- get_site(sitename = 'Marion Lake%')
  louise <- get_site(sitename = 'Louise Pond%')
  western.sites <- rbind(marion, louise)
  western.data  <- get_dataset(western.sites)
  get_download(western.data)
  
}

test_that('some of the functions in the paper still work', 
{
  expect_true(nrow(get_site(sitename = 'Marion Lake%')) == 1)
  expect_true(nrow(get_site(sitename = 'Louise Pond%')) == 1)
  expect_is(block_one(), 'download_list')
  expect_is(get_site(loc = c(-140, 45, -110, 65)), 'site')
  expect_is(get_dataset(loc = c(-140, 45, -110, 65),
                        datasettype = 'pollen',
                        taxonname = 'Pinus%'), 'dataset_list')
  expect_is(get_dataset(loc = c(-120, 46, -110, 50),
                        datasettype = 'pollen',
                        taxonname = 'Pinus%'), 'dataset_list')
  #expect_is(compile_taxa(get_download(get_dataset(loc = c(-120, 46, -110, 50),
  #                                                datasettype = 'pollen',
  #                                                taxonname = 'Pinus%')), 'P25'), 'download_list')
})

