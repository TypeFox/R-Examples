context("span")

test_that("operating on taxon objects", {
  out <- make_taxon(genus="Poa", epithet="annua", authority="L.",
                    family='Poaceae', clazz='Poales', kingdom='Plantae', variety='annua')
  aa <- out %>% span(kingdom, genus)

  expect_is(aa, "taxon")
  expect_is(aa$binomial, "binomial")
  expect_is(aa$binomial$authority, "character")
  expect_is(aa$grouping, "grouping")
  expect_is(aa$grouping$kingdom, "taxonref")
  expect_is(aa$grouping$kingdom$rank, "character")
  expect_equal(length(aa), 2)

  expect_equal(length(out$grouping), 6)
  expect_equal(length(aa$grouping), 4)
})

df <- data.frame(class=c('Magnoliopsida','Magnoliopsida','Magnoliopsida',
                         'Magnoliopsida','Magnoliopsida','Magnoliopsida'),
                 order=c('Asterales','Asterales','Fagales','Poales','Poales','Poales'),
                 family=c('Asteraceae','Asteraceae','Fagaceae','Poaceae','Poaceae','Poaceae'),
                 genus=c('Helianthus','Helianthus','Quercus','Poa','Festuca','Holodiscus'),
                 stringsAsFactors = FALSE)
df2 <- taxon_df(df)

test_that("operating on taxonomic data.frames", {
  aa <- df2 %>% span(order, genus)
  bb <- df2 %>% span(family, genus)

  expect_is(aa, "taxondf")
  expect_named(aa, c('order', 'family', 'genus'))
  expect_more_than(NCOL(df2), NCOL(aa))

  expect_is(bb, "taxondf")
  expect_named(bb, c('family', 'genus'))
  expect_more_than(NCOL(df2), NCOL(bb))
  expect_more_than(NCOL(aa), NCOL(bb))
})

test_that("works on taxa object", {
  aa <- df2 %>% scatter %>% span(family, species)

  expect_is(aa, "list")
  expect_equal(length(aa), NROW(df2))

  expect_is(aa[[1]], "taxon")
  expect_is(aa[[1]]$binomial, "binomial")
  expect_is(aa[[1]]$binomial$authority, "character")
  expect_is(aa[[1]]$grouping, "grouping")
  expect_is(aa[[1]]$grouping$family, "taxonref")
  expect_is(aa[[1]]$grouping$family$rank, "character")
  expect_equal(length(aa[[1]]), 2)
})
