context("scatter_assemble")

df <- data.frame(class=c('Magnoliopsida','Magnoliopsida','Magnoliopsida',
                         'Magnoliopsida','Magnoliopsida','Magnoliopsida'),
         order=c('Asterales','Asterales','Fagales','Poales','Poales','Poales'),
         family=c('Asteraceae','Asteraceae','Fagaceae','Poaceae','Poaceae','Poaceae'),
         genus=c('Helianthus','Helianthus','Quercus','Poa','Festuca','Holodiscus'),
         stringsAsFactors = FALSE)
df2 <- taxon_df(df)

test_that("scatter works", {
  aa <- df2 %>% scatter()

  expect_is(df2, "taxondf")

  expect_is(aa, "taxa")
  expect_is(aa[[1]], "taxon")
  expect_is(aa[[2]], "taxon")

  expect_is(aa[[3]]$binomial, "binomial")
  expect_is(aa[[3]]$grouping, "grouping")
})

test_that("assemble works", {
  bb <- df2 %>% scatter() %>% assemble

  expect_is(bb, "taxondf")
  expect_equal(bb$clazz[1], "Magnoliopsida")
  expect_equal(bb$genus[1], "Helianthus")

  ## FIXME - ideally, these should be identical
  # cc <- bb
  # cc$species <- NULL
  # expect_identical(df2, cc)
})

test_that("scatter fails well", {
  expect_error(scatter(), "no applicable method")
  expect_error(scatter("ad"), "no applicable method")
})

test_that("assemble fails well", {
  expect_error(assemble(), "no applicable method")
  expect_error(assemble("ad"), "no applicable method")
})
